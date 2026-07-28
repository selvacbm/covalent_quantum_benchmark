#!/usr/bin/env python3
"""Merge curated activities, RDKit descriptors, and xTB state statistics.

Raw xTB total energies are intentionally not used as machine-learning inputs:
they mostly scale with molecular size.  The quantum feature block consists of
state-aware HOMO-LUMO-gap summaries, which are comparable across compounds.
"""
import argparse
import csv
import statistics
from collections import defaultdict
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, Lipinski
from rdkit.Chem.Scaffolds import MurckoScaffold


def number(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def descriptor_row(molecule):
    return {
        "molecular_weight": Descriptors.MolWt(molecule),
        "tpsa": Descriptors.TPSA(molecule),
        "logp": Crippen.MolLogP(molecule),
        "hbond_donors": Lipinski.NumHDonors(molecule),
        "hbond_acceptors": Lipinski.NumHAcceptors(molecule),
        "rotatable_bonds": Lipinski.NumRotatableBonds(molecule),
        "aromatic_rings": Lipinski.NumAromaticRings(molecule),
        "fraction_csp3": Lipinski.FractionCSP3(molecule),
        "heavy_atom_count": Lipinski.HeavyAtomCount(molecule),
        "murcko_scaffold": MurckoScaffold.MurckoScaffoldSmiles(mol=molecule),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("benchmark", type=Path)
    parser.add_argument("quantum_states", type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--excluded-out", type=Path)
    args = parser.parse_args()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    excluded_path = args.excluded_out or args.out.parent / "quantum_excluded_records.csv"

    state_groups = defaultdict(list)
    all_state_counts = defaultdict(int)
    for state in csv.DictReader(args.quantum_states.open()):
        record_id = state.get("record_id", "")
        all_state_counts[record_id] += 1
        gap = number(state.get("homo_lumo_gap_ev"))
        if state.get("xtb_status") == "success" and gap is not None:
            state_groups[record_id].append(gap)

    output, excluded = [], []
    for row in csv.DictReader(args.benchmark.open()):
        record_id = row.get("record_id", "")
        gaps = state_groups.get(record_id, [])
        if not gaps:
            excluded.append({
                "record_id": record_id,
                "target_label": row.get("target_label", ""),
                "reason": "no_successful_quantum_state",
                "states_requested": all_state_counts.get(record_id, 0),
            })
            continue
        molecule = Chem.MolFromSmiles(row.get("canonical_smiles", ""))
        if molecule is None:
            excluded.append({
                "record_id": record_id,
                "target_label": row.get("target_label", ""),
                "reason": "invalid_smiles_after_benchmark",
                "states_requested": all_state_counts.get(record_id, 0),
            })
            continue
        descriptors = descriptor_row(molecule)
        gap_std = statistics.pstdev(gaps) if len(gaps) > 1 else 0.0
        item = dict(row)
        item.update({
            **{name: f"{value:.8f}" if isinstance(value, float) else value for name, value in descriptors.items()},
            "y": "1" if row.get("activity_class") == "active" else "0",
            "quantum_states_requested": str(all_state_counts.get(record_id, 0)),
            "quantum_states_success": str(len(gaps)),
            "mean_gap_ev": f"{statistics.mean(gaps):.8f}",
            "min_gap_ev": f"{min(gaps):.8f}",
            "max_gap_ev": f"{max(gaps):.8f}",
            "gap_spread_ev": f"{max(gaps) - min(gaps):.8f}",
            "gap_std_ev": f"{gap_std:.8f}",
            "quantum_method": "GFN2-xTB_single_point_on_RDKit_MMFF_geometry",
        })
        output.append(item)

    if not output:
        raise SystemExit("No records have successful quantum states; inspect xTB outputs before modelling.")
    fields = list(output[0].keys())
    with args.out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(output)
    with excluded_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["record_id", "target_label", "reason", "states_requested"])
        writer.writeheader()
        writer.writerows(excluded)
    print(f"Wrote {len(output)} feature rows to {args.out}")
    print(f"Excluded {len(excluded)} records without usable quantum descriptors: {excluded_path}")


if __name__ == "__main__":
    main()
