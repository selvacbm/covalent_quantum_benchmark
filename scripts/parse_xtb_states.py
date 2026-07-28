#!/usr/bin/env python3
"""Parse successful xTB state calculations into a provenance-preserving table."""
import argparse
import csv
import re
from pathlib import Path

ENERGY = re.compile(r"TOTAL ENERGY\s+([-0-9.]+)\s+Eh")
GAP = re.compile(r"HOMO-LUMO GAP\s+([-0-9.]+)\s+eV")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("states", type=Path)
    parser.add_argument("--results", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    output = []
    for row in csv.DictReader(args.states.open()):
        text_path = args.results / row["state_id"] / "xtb.out"
        text = text_path.read_text(errors="ignore") if text_path.exists() else ""
        energy = ENERGY.search(text)
        gap = GAP.search(text)
        success = "normal termination" in text and energy is not None and gap is not None
        output.append({
            **row,
            "total_energy_eh": energy.group(1) if energy else "",
            "homo_lumo_gap_ev": gap.group(1) if gap else "",
            "calculation_type": "GFN2-xTB_single_point_on_RDKit_MMFF_geometry",
            "xtb_status": "success" if success else "failed",
        })
    fields = list(output[0]) if output else ["state_id"]
    with args.out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(output)
    print(f"Wrote {len(output)} parsed state rows to {args.out}")


if __name__ == "__main__":
    main()
