#!/usr/bin/env python3
"""
Regenerate all GlycoShape SNFG SVG files from WURCS data in GLYCOSHAPE.json
"""

import json
import subprocess
import sys
from pathlib import Path

def main():
    # Read GLYCOSHAPE.json
    glycoshape_path = Path("GLYCOSHAPE.json")
    if not glycoshape_path.exists():
        print(f"Error: {glycoshape_path} not found")
        sys.exit(1)

    with open(glycoshape_path, 'r') as f:
        data = json.load(f)

    print(f"Loaded {len(data)} entries from GLYCOSHAPE.json")

    success_count = 0
    failure_count = 0
    skipped_count = 0

    for key, entry in data.items():
        archetype = entry.get('archetype', {})
        wurcs = archetype.get('wurcs')

        if not wurcs:
            skipped_count += 1
            continue

        # Use crabwurcs CLI to render SNFG
        result = subprocess.run(
            ['./target/debug/crabwurcs', 'render', '--from', 'wurcs', '--input-file', 'false', wurcs],
            capture_output=True,
            text=True
        )

        if result.returncode == 0:
            svg_filename = Path(f"{key}.snfg.svg")
            try:
                with open(svg_filename, 'w') as f:
                    f.write(result.stdout)
                success_count += 1
                if success_count % 100 == 0:
                    print(f"Generated {success_count} SVG files so far...")
            except Exception as e:
                print(f"Error writing {svg_filename}: {e}")
                failure_count += 1
        else:
            print(f"Error rendering {key}: {result.stderr}")
            failure_count += 1

    print(f"\n--- Summary ---")
    print(f"Total entries: {len(data)}")
    print(f"Successfully generated: {success_count} SVG files")
    print(f"Failed: {failure_count}")
    print(f"Skipped (no WURCS): {skipped_count}")

if __name__ == '__main__':
    main()
