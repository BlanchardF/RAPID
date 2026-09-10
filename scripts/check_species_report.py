#!/usr/bin/env python3
"""
RAPID — scripts/check_species_report.py
========================================
Build a cross-pair species summary from ipcress results, plus a self-contained
interactive HTML page to explore which species each primer pair amplifies and
to filter species out ("I don't care about this species") to see which pairs
remain specific.

Species are parsed from ipcress "Target: ... TSA: Genus species ..." lines,
exactly as in split_ipcress_by_pair.py (only meaningful for TSA-style DBs).
"""

import sys
import json
from pathlib import Path
from collections import defaultdict


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

def species_from_target_line(stripped):
    """Return 'Genus species' from a 'Target: ... TSA: Genus species ...' line, else None."""
    if stripped.startswith("Target:") and "TSA:" in stripped:
        idx = stripped.index("TSA:") + len("TSA:")
        words = stripped[idx:].split()
        if len(words) >= 2:
            return f"{words[0]} {words[1]}"
    return None


def count_species_in_text(text):
    """Count species occurrences (amplicons) in a per-pair ipcress text block."""
    counts = defaultdict(int)
    for raw in text.splitlines():
        sp = species_from_target_line(raw.strip())
        if sp:
            counts[sp] += 1
    return counts


def parse_combined_file(text):
    """Parse a combined ipcress file: 'Experiment:' delimits pairs."""
    pairs = defaultdict(lambda: defaultdict(int))
    current = None
    for raw in text.splitlines():
        s = raw.strip()
        if s.startswith("Experiment:"):
            current = s.split("Experiment:", 1)[1].strip()
            pairs.setdefault(current, defaultdict(int))
        else:
            sp = species_from_target_line(s)
            if sp and current is not None:
                pairs[current][sp] += 1
    return {pid: dict(spd) for pid, spd in pairs.items()}


def load_pairs(input_path):
    """Return {pair_id: {species: count}} from a directory or a combined file."""
    p = Path(input_path)
    if p.is_dir():
        pairs = {}
        block_files = sorted(f for f in p.glob("*.txt")
                             if not f.name.endswith("_species.txt"))
        if block_files:
            for f in block_files:
                pairs[f.stem] = dict(count_species_in_text(f.read_text()))
        # Augment/fallback with *_species.txt (presence only, count treated as 1)
        for f in sorted(p.glob("*_species.txt")):
            pid = f.name[:-len("_species.txt")]
            entry = pairs.setdefault(pid, {})
            for line in f.read_text().splitlines():
                sp = line.strip()
                if sp and sp not in entry:
                    entry[sp] = 1
        if not pairs:
            print(f"[check_species_report] WARNING: no per-pair files found in {p}",
                  file=sys.stderr)
        return pairs
    if p.is_file():
        return parse_combined_file(p.read_text())
    print(f"[check_species_report] ERROR: input not found: {p}", file=sys.stderr)
    sys.exit(1)


# ---------------------------------------------------------------------------
# TSV matrix
# ---------------------------------------------------------------------------

def write_matrix_tsv(pairs, all_species, out_tsv):
    Path(out_tsv).parent.mkdir(parents=True, exist_ok=True)
    with open(out_tsv, "w") as fh:
        header = ["primer_pair"] + all_species + ["n_species"]
        fh.write("\t".join(header) + "\n")
        for pid in sorted(pairs):
            spd = pairs[pid]
            row = [pid]
            for sp in all_species:
                c = spd.get(sp, 0)
                row.append(str(c) if c else "")
            row.append(str(len([s for s in spd if spd[s] > 0])))
            fh.write("\t".join(row) + "\n")


# ---------------------------------------------------------------------------
# HTML explorer
# ---------------------------------------------------------------------------

HTML_TEMPLATE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>RAPID check — species / pair explorer</title>
<style>
  :root {
    --bg: #0f1115; --panel: #171a21; --line: #262b36; --fg: #e6e9ef;
    --muted: #98a2b3; --accent: #5b9dff; --good: #2ea043; --good-bg: #12331d;
    --warn: #d29922; --chip: #222834;
  }
  * { box-sizing: border-box; }
  body { margin: 0; font: 14px/1.5 -apple-system, Segoe UI, Roboto, sans-serif;
         background: var(--bg); color: var(--fg); }
  header { padding: 18px 22px; border-bottom: 1px solid var(--line); }
  h1 { margin: 0 0 4px; font-size: 18px; }
  .sub { color: var(--muted); font-size: 13px; }
  .layout { display: grid; grid-template-columns: 300px 1fr; gap: 0; height: calc(100vh - 74px); }
  .side { border-right: 1px solid var(--line); overflow-y: auto; padding: 14px; background: var(--panel); }
  .main { overflow-y: auto; padding: 16px 22px; }
  .toolbar { display: flex; flex-wrap: wrap; gap: 8px; align-items: center; margin-bottom: 12px; }
  input[type=search], select {
    background: #0d1017; color: var(--fg); border: 1px solid var(--line);
    border-radius: 8px; padding: 7px 10px; font-size: 13px;
  }
  button {
    background: #1f2530; color: var(--fg); border: 1px solid var(--line);
    border-radius: 8px; padding: 7px 12px; cursor: pointer; font-size: 13px;
  }
  button:hover { border-color: var(--accent); }
  .sp-item { display: flex; align-items: center; gap: 8px; padding: 4px 2px; }
  .sp-item label { flex: 1; cursor: pointer; font-style: italic; }
  .sp-item .cnt { color: var(--muted); font-size: 12px; }
  .legend { color: var(--muted); font-size: 12px; margin: 4px 0 10px; }
  .pair { border: 1px solid var(--line); border-radius: 10px; padding: 10px 12px;
          margin-bottom: 8px; background: var(--panel); }
  .pair.clean { border-color: var(--good); background: var(--good-bg); }
  .pair-head { display: flex; align-items: center; gap: 10px; }
  .pair-id { font-weight: 600; }
  .badge { margin-left: auto; font-size: 12px; padding: 2px 9px; border-radius: 999px;
           background: var(--chip); border: 1px solid var(--line); }
  .badge.zero { background: var(--good); color: #04220f; border-color: var(--good); font-weight: 700; }
  .chips { margin-top: 7px; display: flex; flex-wrap: wrap; gap: 6px; }
  .chip { font-size: 12px; padding: 2px 8px; border-radius: 999px; background: var(--chip);
          border: 1px solid var(--line); font-style: italic; }
  .chip.excluded { opacity: .35; text-decoration: line-through; }
  .count-tag { font-style: normal; color: var(--muted); }
  .muted { color: var(--muted); }
  .stat { display: inline-block; margin-right: 16px; }
  .stat b { color: var(--accent); }
  a { color: var(--accent); }
</style>
</head>
<body>
<header>
  <h1>RAPID check — species / pair explorer</h1>
  <div class="sub">Uncheck species on the left to ignore them ("don't count this one").
  A pair turns <b style="color:var(--good)">green</b> once it no longer amplifies any checked species.</div>
</header>
<div class="layout">
  <aside class="side">
    <div class="toolbar">
      <button onclick="setAll(true)">Check all</button>
      <button onclick="setAll(false)">Uncheck all</button>
    </div>
    <input type="search" id="spSearch" placeholder="Filter species..." oninput="renderSpecies()" style="width:100%; margin-bottom:10px;">
    <div class="legend" id="spCount"></div>
    <div id="speciesList"></div>
  </aside>
  <main class="main">
    <div class="toolbar">
      <label><input type="checkbox" id="onlyClean" onchange="renderPairs()"> Show only specific pairs (0 checked species)</label>
    </div>
    <div class="legend" id="stats"></div>
    <div id="pairsList"></div>
  </main>
</div>

<script>
const DATA = __DATA__;
// DATA = { pairs: { pairId: { "Genus species": count, ... } }, species: [ ... ] }
const excluded = new Set();   // species the user unchecked (ignored)

function speciesCounts() {
  // total amplicons per species across all pairs (for display)
  const t = {};
  for (const sp of DATA.species) t[sp] = 0;
  for (const pid in DATA.pairs)
    for (const sp in DATA.pairs[pid]) t[sp] += DATA.pairs[pid][sp];
  return t;
}
const SP_TOTAL = speciesCounts();

function setAll(checked) {
  excluded.clear();
  if (!checked) DATA.species.forEach(sp => excluded.add(sp));
  renderSpecies(); renderPairs();
}

function toggle(sp) {
  if (excluded.has(sp)) excluded.delete(sp); else excluded.add(sp);
  renderSpecies(); renderPairs();
}

function renderSpecies() {
  const q = (document.getElementById('spSearch').value || '').toLowerCase();
  const list = document.getElementById('speciesList');
  const shown = DATA.species.filter(sp => sp.toLowerCase().includes(q));
  list.innerHTML = shown.map(sp => {
    const checked = excluded.has(sp) ? '' : 'checked';
    return `<div class="sp-item">
      <input type="checkbox" ${checked} onchange="toggle('${sp.replace(/'/g,"\\'")}')" id="cb_${btoa(unescape(encodeURIComponent(sp)))}">
      <label for="cb_${btoa(unescape(encodeURIComponent(sp)))}">${sp}</label>
      <span class="cnt">${SP_TOTAL[sp]}</span>
    </div>`;
  }).join('');
  const active = DATA.species.length - excluded.size;
  document.getElementById('spCount').textContent =
    `${active} / ${DATA.species.length} species counted`;
}

function countedSpeciesForPair(pid) {
  // species amplified by this pair that are NOT excluded
  const res = [];
  for (const sp in DATA.pairs[pid])
    if (!excluded.has(sp)) res.push([sp, DATA.pairs[pid][sp]]);
  return res;
}

function renderPairs() {
  const onlyClean = document.getElementById('onlyClean').checked;
  const pids = Object.keys(DATA.pairs);
  const rows = pids.map(pid => {
    const counted = countedSpeciesForPair(pid);
    return { pid, counted, n: counted.length,
             all: Object.keys(DATA.pairs[pid]).length };
  }).sort((a, b) => a.n - b.n || a.pid.localeCompare(b.pid));

  const nClean = rows.filter(r => r.n === 0).length;
  document.getElementById('stats').innerHTML =
    `<span class="stat"><b>${pids.length}</b> pairs</span>
     <span class="stat"><b>${DATA.species.length}</b> species total</span>
     <span class="stat"><b style="color:var(--good)">${nClean}</b> pair(s) with no checked species</span>`;

  const visible = onlyClean ? rows.filter(r => r.n === 0) : rows;
  const html = visible.map(r => {
    const chips = Object.keys(DATA.pairs[r.pid]).sort().map(sp => {
      const ex = excluded.has(sp) ? ' excluded' : '';
      return `<span class="chip${ex}">${sp} <span class="count-tag">x${DATA.pairs[r.pid][sp]}</span></span>`;
    }).join('');
    const cls = r.n === 0 ? 'pair clean' : 'pair';
    const badge = r.n === 0 ? '<span class="badge zero">specific</span>'
                            : `<span class="badge">${r.n} checked species</span>`;
    const chipsHtml = chips ? `<div class="chips">${chips}</div>`
                            : `<div class="chips"><span class="muted">no species detected</span></div>`;
    return `<div class="${cls}">
      <div class="pair-head"><span class="pair-id">${r.pid}</span>${badge}</div>
      ${chipsHtml}
    </div>`;
  }).join('');
  document.getElementById('pairsList').innerHTML =
    html || '<p class="muted">No pair to display.</p>';
}

renderSpecies();
renderPairs();
</script>
</body>
</html>
"""


def write_html(pairs, all_species, out_html):
    data = {"pairs": pairs, "species": all_species}
    blob = json.dumps(data, ensure_ascii=False)
    blob = blob.replace("</", "<\\/")   # safe inside <script>
    html_text = HTML_TEMPLATE.replace("__DATA__", blob)
    Path(out_html).parent.mkdir(parents=True, exist_ok=True)
    Path(out_html).write_text(html_text, encoding="utf-8")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    if len(sys.argv) != 4:
        print(f"Usage: {sys.argv[0]} <input_dir_or_file> <out_tsv> <out_html>",
              file=sys.stderr)
        sys.exit(1)
    input_path, out_tsv, out_html = sys.argv[1:4]

    pairs = load_pairs(input_path)
    all_species = sorted({sp for spd in pairs.values() for sp in spd})

    write_matrix_tsv(pairs, all_species, out_tsv)
    write_html(pairs, all_species, out_html)

    print(f"[check_species_report] pairs   : {len(pairs)}", file=sys.stderr)
    print(f"[check_species_report] species : {len(all_species)}", file=sys.stderr)
    print(f"[check_species_report] matrix  -> {out_tsv}", file=sys.stderr)
    print(f"[check_species_report] explorer-> {out_html}", file=sys.stderr)


if __name__ == "__main__":
    main()
