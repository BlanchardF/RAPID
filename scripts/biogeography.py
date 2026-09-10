#!/usr/bin/env python3
"""
RAPID — scripts/biogeography.py
Generate biogeography maps using GBIF tile layers.
"""

import sys
import time
import argparse
from pathlib import Path

try:
    import requests
except ImportError:
    sys.exit("ERROR: 'requests' not found. pip install requests")

try:
    import folium
except ImportError:
    sys.exit("ERROR: 'folium' not found. pip install folium")


GBIF_MATCH = "https://api.gbif.org/v1/species/match"
GBIF_TILE_TEMPLATE = ("https://api.gbif.org/v2/map/occurrence/density"
                      "/{z}/{x}/{y}@2x.png?taxonKey=__KEY__&bin=hex&hexPerTile=30&style=__STYLE__")


def gbif_tile_url(taxon_key, style):
    """
    Build a GBIF tile URL, keeping {z}/{x}/{y} as literal Leaflet placeholders
    (Folium/Leaflet substitutes those at render time in the browser).
    We must NOT use str.format() here since it would also try to fill {z}{x}{y}.
    """
    return (GBIF_TILE_TEMPLATE
            .replace("__KEY__", str(taxon_key))
            .replace("__STYLE__", style))

STYLES = [
    "classic.poly", "purpleYellow.poly", "fire.poly", "glacial.poly",
    "bluePurple.poly", "oranges.poly", "blues.poly", "greens.poly", "reds.poly",
]
LEGEND_COLORS = [
    "#1a9641", "#9e0142", "#d7191c", "#2c7bb6",
    "#7b2d8b", "#e6550d", "#3182bd", "#31a354", "#cb181d",
]


def get_taxon_key(species_name):
    try:
        r = requests.get(GBIF_MATCH,
                         params={"name": species_name, "strict": "false"},
                         timeout=15)
        r.raise_for_status()
        d = r.json()
    except Exception as e:
        return None, f"request error: {e}"
    key = d.get("usageKey") or d.get("speciesKey")
    if not key:
        return None, f"no match (matchType={d.get('matchType','?')})"
    return key, f"{d.get('matchType','?')} (confidence {d.get('confidence',0)}%)"


def _title(text):
    return folium.Element(
        f'<div style="position:fixed;top:10px;left:50%;transform:translateX(-50%);'
        f'background:white;padding:6px 16px;border-radius:6px;'
        f'box-shadow:2px 2px 8px rgba(0,0,0,.25);font-family:sans-serif;'
        f'font-size:13px;z-index:9999;">{text}</div>'
    )


def make_species_map(species, key, out_path):
    m = folium.Map(location=[20, 0], zoom_start=2, tiles="CartoDB positron")
    folium.TileLayer("OpenStreetMap", name="OpenStreetMap").add_to(m)
    folium.TileLayer(
        tiles=gbif_tile_url(key, "classic.poly"),
        attr='<a href="https://www.gbif.org">GBIF</a>',
        name=f"<i>{species}</i>", overlay=True, opacity=0.85,
    ).add_to(m)
    m.get_root().html.add_child(_title(f"<i>{species}</i> — GBIF occurrences"))
    folium.LayerControl().add_to(m)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    m.save(str(out_path))


def make_combined_map(pair_id, sp_keys, out_path):
    m = folium.Map(location=[20, 0], zoom_start=2, tiles="CartoDB positron")
    folium.TileLayer("OpenStreetMap", name="OpenStreetMap").add_to(m)
    legend_items = []
    for i, (species, key) in enumerate(sorted(sp_keys.items())):
        style = STYLES[i % len(STYLES)]
        color = LEGEND_COLORS[i % len(LEGEND_COLORS)]
        folium.TileLayer(
            tiles=gbif_tile_url(key, style),
            attr='<a href="https://www.gbif.org">GBIF</a>',
            name=f"<i>{species}</i>", overlay=True, opacity=0.80,
        ).add_to(m)
        legend_items.append(
            f'<div style="display:flex;align-items:center;margin:3px 0;">'
            f'<span style="background:{color};width:14px;height:14px;'
            f'border-radius:3px;display:inline-block;margin-right:6px;"></span>'
            f'<i style="font-size:12px;">{species}</i></div>'
        )
    legend = (
        '<div style="position:fixed;bottom:30px;left:15px;background:white;'
        'padding:10px 14px;border-radius:6px;box-shadow:2px 2px 8px rgba(0,0,0,.25);'
        'font-family:sans-serif;z-index:9999;max-height:300px;overflow-y:auto;">'
        f'<b style="font-size:13px;">{pair_id}</b><hr style="margin:4px 0;">'
        + "".join(legend_items) + '</div>'
    )
    m.get_root().html.add_child(folium.Element(legend))
    m.get_root().html.add_child(_title(f"Primer pair <b>{pair_id}</b> — {len(sp_keys)} species"))
    folium.LayerControl(collapsed=False).add_to(m)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    m.save(str(out_path))


def run(per_pair_dir, maps_out_dir):
    per_pair_dir = Path(per_pair_dir)
    maps_out_dir = Path(maps_out_dir)
    species_dir  = maps_out_dir / "species"
    pairs_dir    = maps_out_dir / "pairs"

    sp_files = sorted(per_pair_dir.glob("*_species.txt"))
    if not sp_files:
        print("[RAPID maps] No *_species.txt files found.", file=sys.stderr)
        return

    # Collect all unique species across all pairs
    all_species = set()
    pair_species = {}
    for sf in sp_files:
        pair_id = sf.stem.replace("_species", "")
        sps = [l.strip() for l in sf.read_text().splitlines() if l.strip()]
        pair_species[pair_id] = sps
        all_species.update(sps)

    print(f"[RAPID maps] {len(sp_files)} pair(s), {len(all_species)} unique species.")

    # Step 1: resolve taxon keys once per unique species
    print(f"[RAPID maps] Resolving taxon keys via GBIF …")
    taxon_keys = {}
    for species in sorted(all_species):
        key, info = get_taxon_key(species)
        taxon_keys[species] = key
        status = "OK" if key else "FAIL"
        print(f"  [{status}] {species:<35} {info}")
        time.sleep(0.15)

    # Step 2: per-species maps (shared, skip if already exists)
    print(f"\n[RAPID maps] Per-species maps → {species_dir}/")
    for species, key in taxon_keys.items():
        sp_safe = species.replace(" ", "_")
        out = species_dir / f"{sp_safe}.html"
        if out.exists():
            print(f"  (exists) {species}")
            continue
        if not key:
            print(f"  (skip)   {species} — no taxon key")
            continue
        make_species_map(species, key, out)
        print(f"  (new)    {species}")

    # Step 3: per-pair combined maps
    print(f"\n[RAPID maps] Combined maps → {pairs_dir}/")
    for pair_id, sps in pair_species.items():
        sp_keys = {sp: taxon_keys[sp] for sp in sps
                   if taxon_keys.get(sp) is not None}
        if not sp_keys:
            print(f"  (skip) {pair_id} — no valid taxon keys")
            continue
        out = pairs_dir / f"{pair_id}_combined.html"
        make_combined_map(pair_id, sp_keys, out)
        print(f"  (new)  {pair_id} ({len(sp_keys)} species)")

    print(f"\n[RAPID maps] Done.")


def main():
    p = argparse.ArgumentParser(description="GBIF tile biogeography maps for RAPID.")
    p.add_argument("per_pair_dir", help="Directory with *_species.txt files.")
    p.add_argument("maps_out_dir", help="Output directory for HTML maps.")
    args = p.parse_args()
    run(args.per_pair_dir, args.maps_out_dir)


if __name__ == "__main__":
    main()
