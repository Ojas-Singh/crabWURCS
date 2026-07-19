#!/usr/bin/env python3
"""
Regenerate GlycoShape SNFG SVGs from all 4 notation formats.
Output: <id>_wurcs.svg, <id>_iupac.svg, <id>_iupac_ext.svg, <id>_glycam.svg
"""

import json
import subprocess
import sys
from pathlib import Path
from collections import Counter

CRAB = './target/release/crabwurcs'

CONFIG = [
    ('wurcs',     'wurcs',          'wurcs'),
    ('iupac',     'iupac',          'iupac-condensed'),
    ('iupac_ext', 'iupac_extended', 'iupac-extended'),
    ('glycam',    'glycam',         'glycam'),
]

def main():
    json_path = Path("GLYCOSHAPE.json")
    with open(json_path) as f:
        data = json.load(f)

    out_dir = Path("glycoshape")
    out_dir.mkdir(exist_ok=True)

    counts = Counter()

    for key, entry in sorted(data.items()):
        archetype = entry.get('archetype', {})

        for suffix, src_field, render_fmt in CONFIG:
            inp = archetype.get(src_field)
            if not inp:
                counts[f'{suffix}_no_field'] += 1
                continue

            res = subprocess.run(
                [CRAB, 'render', '--from', render_fmt],
                input=inp, capture_output=True, text=True
            )
            if res.returncode == 0:
                (out_dir / f"{key}_{suffix}.svg").write_text(res.stdout)
                counts[f'{suffix}_ok'] += 1
            else:
                counts[f'{suffix}_fail'] += 1

        if sum(1 for s, _, _ in CONFIG for c in ['_ok', '_no_field'] if counts[s+c]) == 0:
            continue
        total_done = sum(counts[f'{s}_ok'] + counts[f'{s}_no_field'] for s, _, _ in CONFIG)
        if total_done % 400 == 0:
            print(f"Progress: {key}")

    print(f"\n--- {len(data)} entries ---")
    for suffix, _, _ in CONFIG:
        ok = counts[f'{suffix}_ok']
        skip = counts[f'{suffix}_no_field']
        fail = counts[f'{suffix}_fail']
        print(f"  {suffix}: {ok} ok, {skip} no-field, {fail} fail")

if __name__ == '__main__':
    main()
