#!/usr/bin/env python3
"""Build a bounded tautomer-state ensemble and 3D xyz geometries."""
import argparse
import csv
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.MolStandardize import rdMolStandardize


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--xyz-dir", required=True, type=Path)
    parser.add_argument("--max-states", type=int, default=4)
    parser.add_argument("--skip-mmff", action="store_true")
    args = parser.parse_args()
    if not 1 <= args.max_states <= 8:
        raise SystemExit("--max-states must be between 1 and 8")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.xyz_dir.mkdir(parents=True, exist_ok=True)

    enumerator = rdMolStandardize.TautomerEnumerator()
    enumerator.SetMaxTautomers(12)
    enumerator.SetMaxTransforms(100)
    states_out = []
    for record_number, row in enumerate(csv.DictReader(args.input.open()), 1):
        molecule = Chem.MolFromSmiles(row.get("canonical_smiles", ""))
        if molecule is None:
            print(f"SKIP invalid SMILES: {row.get('record_id', record_number)}")
            continue
        candidates = [molecule] + list(enumerator.Enumerate(molecule))
        unique = []
        seen_smiles = set()
        for state in candidates:
            smiles = Chem.MolToSmiles(state, isomericSmiles=True)
            if smiles not in seen_smiles:
                seen_smiles.add(smiles)
                unique.append((smiles, state))
        for index, (smiles, state) in enumerate(unique[:args.max_states]):
            state_id = f"{row['record_id']}__state{index:02d}"
            mol3d = Chem.AddHs(Chem.MolFromSmiles(smiles))
            if AllChem.EmbedMolecule(mol3d, randomSeed=20260723 + index, useRandomCoords=True) != 0:
                print(f"SKIP embedding: {state_id}")
                continue
            if not args.skip_mmff:
                try:
                    AllChem.MMFFOptimizeMolecule(mol3d, maxIters=300)
                except Exception:
                    pass
            mol3d = Chem.RemoveHs(mol3d)
            xyz_path = args.xyz_dir / f"{state_id}.xyz"
            Chem.MolToXYZFile(mol3d, str(xyz_path))
            states_out.append({
                "state_id": state_id,
                "record_id": row["record_id"],
                "molecule_chembl_id": row["molecule_chembl_id"],
                "target_label": row["target_label"],
                "document_chembl_id": row["document_chembl_id"],
                "state_index": index,
                "state_smiles": smiles,
                "formal_charge": Chem.GetFormalCharge(mol3d),
                "xyz_path": str(xyz_path),
            })
        if record_number % 25 == 0:
            print(f"processed={record_number} states={len(states_out)}", flush=True)

    fields = [
        "state_id", "record_id", "molecule_chembl_id", "target_label",
        "document_chembl_id", "state_index", "state_smiles", "formal_charge", "xyz_path",
    ]
    with args.out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(states_out)
    print(f"Wrote {len(states_out)} quantum states to {args.out}")


if __name__ == "__main__":
    main()
