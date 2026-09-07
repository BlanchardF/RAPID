#!/usr/bin/env python3
"""
RAPID — scripts/split_ipcress_by_pair.py
Split a combined ipcress result file into one file per primer pair,
optionally extracting TSA species names per pair.

Usage:
    python3 split_ipcress_by_pair.py <combined_raw> <per_pair_dir> <true|false>

The third argument controls whether species names are extracted from
'Target: ... TSA: Genus species ...' lines (only meaningful for TSA hits).
"""

import sys
from pathlib import Path


def split_by_pair(raw_text, per_pair_dir, extract_species):
    """
    Parse ipcress blocks and write one file per primer pair (Experiment ID).
    Also writes <id>_species.txt per pair when extract_species=True and
    at least one TSA species line was found.

    Returns a dict {pair_id: {"hits": N, "species": set(...)}}.
    """
    per_pair_dir = Path(per_pair_dir)
    per_pair_dir.mkdir(parents=True, exist_ok=True)

    pairs_data = {}
    current_id    = None
    current_lines = []

    def flush(exp_id, lines):
        if exp_id is None or not lines:
            return
        if exp_id not in pairs_data:
            pairs_data[exp_id] = {"lines": [], "species": set(), "hits": 0}
        pairs_data[exp_id]["lines"].extend(lines)
        pairs_data[exp_id]["lines"].append("")
        pairs_data[exp_id]["hits"] += 1

        if extract_species:
            for line in lines:
                stripped = line.strip()
                if stripped.startswith("Target:") and "TSA:" in stripped:
                    idx   = stripped.index("TSA:") + len("TSA:")
                    words = stripped[idx:].split()
                    if len(words) >= 2:
                        pairs_data[exp_id]["species"].add(f"{words[0]} {words[1]}")

    for raw_line in raw_text.splitlines():
        line = raw_line.rstrip()

        if line.strip() == "Ipcress result":
            flush(current_id, current_lines)
            current_lines = [line]
            current_id    = None
            continue

        stripped = line.strip()
        if stripped.startswith("Experiment:"):
            current_id = stripped.split("Experiment:", 1)[1].strip()

        current_lines.append(line)

    flush(current_id, current_lines)

    for pair_id, data in pairs_data.items():
        with open(per_pair_dir / f"{pair_id}.txt", "w") as fh:
            fh.write("\n".join(data["lines"]) + "\n")

        if extract_species and data["species"]:
            with open(per_pair_dir / f"{pair_id}_species.txt", "w") as fh:
                for sp in sorted(data["species"]):
                    fh.write(sp + "\n")

    return pairs_data


def main():
    if len(sys.argv) != 4:
        print(f"Usage: {sys.argv[0]} <combined_raw> <per_pair_dir> <true|false>",
              file=sys.stderr)
        sys.exit(1)

    combined_path, per_pair_dir, extract_flag = sys.argv[1:4]
    extract_species = extract_flag.strip().lower() == "true"

    raw_text = Path(combined_path).read_text()
    pairs_data = split_by_pair(raw_text, per_pair_dir, extract_species)

    if not pairs_data:
        print("[split_by_pair] No hits found in any primer pair.")
        return

    print(f"[split_by_pair] Results by primer pair → {per_pair_dir}\n")
    for pair_id, data in sorted(pairs_data.items()):
        species = sorted(data["species"])
        print(f"  {pair_id}")
        print(f"    hits    : {data['hits']}")
        if extract_species:
            print(f"    species : {len(species)}")
            for sp in species:
                print(f"              {sp}")
        print()


if __name__ == "__main__":
    main()
