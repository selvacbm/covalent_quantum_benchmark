#!/usr/bin/env python3
"""Build a strict IC50-only, electrophile-enriched activity benchmark.

Warhead SMARTS identify plausible electrophilic groups; they are not evidence
that covalent bond formation occurred in the assay. Source provenance is kept
for later document-grouped validation and manual evidence checking.
"""
import argparse
import csv
import math
from collections import Counter, defaultdict
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import Descriptors

WARHEAD_SMARTS = {
    "acrylamide": "C=CC(=O)N",
    "propiolamide": "C#CC(=O)N",
    "chloroacetamide": "[Cl,Br]CC(=O)N",
    "vinyl_sulfonamide": "C=CS(=O)(=O)N",
    "sulfonyl_fluoride": "S(=O)(=O)F",
}
WARHEADS = {name: Chem.MolFromSmarts(pattern) for name, pattern in WARHEAD_SMARTS.items()}


def numeric_pchembl(row):
    try:
        return float(row.get("pchembl_value", ""))
    except ValueError:
        pass
    try:
        value = float(row.get("standard_value", ""))
        return 9.0 - math.log10(value)
    except (ValueError, TypeError):
        return None


def warhead_class(molecule):
    classes = [name for name, query in WARHEADS.items() if molecule.HasSubstructMatch(query)]
    return ";".join(classes)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--qc-out", required=True, type=Path)
    parser.add_argument("--min-per-class", type=int, default=25)
    parser.add_argument("--max-mw", type=float, default=900.0)
    parser.add_argument("--active-threshold", type=float, default=6.0)
    parser.add_argument("--inactive-threshold", type=float, default=5.0)
    args = parser.parse_args()

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.qc_out.parent.mkdir(parents=True, exist_ok=True)
    raw_rows = list(csv.DictReader(args.input.open()))
    counts = Counter()
    candidates = []

    for row in raw_rows:
        label = row.get("target_label", "")
        counts[(label, "raw")] += 1
        if row.get("standard_type", "").upper() != "IC50":
            continue
        if row.get("standard_relation", "") != "=":
            continue
        if row.get("standard_units", "").lower() != "nm":
            continue
        counts[(label, "strict_ic50")] += 1
        pchembl = numeric_pchembl(row)
        if pchembl is None:
            continue
        if pchembl >= args.active_threshold:
            activity_class = "active"
        elif pchembl < args.inactive_threshold:
            activity_class = "inactive"
        else:
            counts[(label, "borderline_excluded")] += 1
            continue
        molecule = Chem.MolFromSmiles(row.get("canonical_smiles", ""))
        if molecule is None:
            counts[(label, "invalid_structure")] += 1
            continue
        molecular_weight = Descriptors.MolWt(molecule)
        if molecular_weight > args.max_mw:
            counts[(label, "mw_excluded")] += 1
            continue
        warheads = warhead_class(molecule)
        if not warheads:
            counts[(label, "no_warhead")] += 1
            continue
        item = dict(row)
        item.update({
            "record_id": f"{label}__{row['molecule_chembl_id']}",
            "pchembl_value": f"{pchembl:.6f}",
            "activity_class": activity_class,
            "warhead_class": warheads,
            "molecular_weight": f"{molecular_weight:.6f}",
        })
        candidates.append(item)

    # One record per target-molecule prevents exact compound leakage.  Mixed
    # active/inactive calls are removed rather than resolved by selecting the
    # most favourable report.  For label-consistent records, retain the report
    # closest to the median potency and keep the number of supporting reports.
    # This is deliberately conservative: it avoids making an apparent model
    # success from a cherry-picked activity record.
    grouped = defaultdict(list)
    for row in candidates:
        grouped[(row["target_label"], row["molecule_chembl_id"])].append(row)

    selected = []
    duplicate_rows = []
    for (target, molecule_id), members in grouped.items():
        classes = {row["activity_class"] for row in members}
        if len(classes) > 1:
            counts[(target, "conflicting_duplicate_excluded")] += 1
            duplicate_rows.append({
                "target_label": target,
                "molecule_chembl_id": molecule_id,
                "resolution": "excluded_conflicting_labels",
                "n_records": len(members),
                "activity_classes": ";".join(sorted(classes)),
                "documents": ";".join(sorted({row["document_chembl_id"] for row in members if row["document_chembl_id"]})),
            })
            continue
        values = sorted(float(row["pchembl_value"]) for row in members)
        median = values[len(values) // 2]
        chosen = min(members, key=lambda row: abs(float(row["pchembl_value"]) - median))
        chosen = dict(chosen)
        chosen["source_record_count"] = str(len(members))
        chosen["source_document_count"] = str(len({row["document_chembl_id"] for row in members if row["document_chembl_id"]}))
        selected.append(chosen)
        duplicate_rows.append({
            "target_label": target,
            "molecule_chembl_id": molecule_id,
            "resolution": "median_potency_representative",
            "n_records": len(members),
            "activity_classes": next(iter(classes)),
            "documents": ";".join(sorted({row["document_chembl_id"] for row in members if row["document_chembl_id"]})),
        })

    by_target_class = Counter((row["target_label"], row["activity_class"]) for row in selected)
    eligible_targets = {
        target for target in {row["target_label"] for row in selected}
        if by_target_class[(target, "active")] >= args.min_per_class
        and by_target_class[(target, "inactive")] >= args.min_per_class
    }
    final_rows = [row for row in selected if row["target_label"] in eligible_targets]
    final_rows.sort(key=lambda row: (row["target_label"], row["molecule_chembl_id"]))

    fields = [
        "record_id", "target_label", "molecule_chembl_id", "target_chembl_id",
        "assay_chembl_id", "document_chembl_id", "standard_type", "standard_relation",
        "standard_value", "standard_units", "pchembl_value", "activity_class",
        "warhead_class", "molecular_weight", "canonical_smiles", "activity_comment",
        "source_record_count", "source_document_count",
    ]
    with args.out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows([{field: row.get(field, "") for field in fields} for row in final_rows])

    qc_rows = []
    targets = sorted({row.get("target_label", "") for row in raw_rows})
    for target in targets:
        qc_rows.append({
            "target_label": target,
            "raw_records": counts[(target, "raw")],
            "strict_ic50_records": counts[(target, "strict_ic50")],
            "borderline_excluded": counts[(target, "borderline_excluded")],
            "no_warhead": counts[(target, "no_warhead")],
            "mw_excluded": counts[(target, "mw_excluded")],
            "conflicting_duplicate_excluded": counts[(target, "conflicting_duplicate_excluded")],
            "unique_active": by_target_class[(target, "active")],
            "unique_inactive": by_target_class[(target, "inactive")],
            "unique_documents": len({row["document_chembl_id"] for row in selected if row["target_label"] == target}),
            "eligible": int(target in eligible_targets),
        })
    with args.qc_out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(qc_rows[0].keys()))
        writer.writeheader()
        writer.writerows(qc_rows)

    duplicate_path = args.qc_out.parent / "duplicate_resolution.csv"
    with duplicate_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=[
            "target_label", "molecule_chembl_id", "resolution", "n_records",
            "activity_classes", "documents",
        ])
        writer.writeheader()
        writer.writerows(sorted(duplicate_rows, key=lambda row: (row["target_label"], row["molecule_chembl_id"])))

    queue_path = args.qc_out.parent / "electrophile_evidence_queue.csv"
    queue_fields = ["target_label", "document_chembl_id", "warhead_class", "n_records", "note"]
    queue_groups = defaultdict(list)
    for row in final_rows:
        queue_groups[(row["target_label"], row["document_chembl_id"], row["warhead_class"])].append(row)
    with queue_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=queue_fields)
        writer.writeheader()
        for (target, document, warhead), members in sorted(queue_groups.items()):
            writer.writerow({
                "target_label": target,
                "document_chembl_id": document,
                "warhead_class": warhead,
                "n_records": len(members),
                "note": "Confirm covalent mechanism in the primary source before using covalent-inhibitor wording.",
            })

    print(f"Wrote {len(final_rows)} benchmark records to {args.out}")
    print(f"Eligible targets: {', '.join(sorted(eligible_targets)) or 'none'}")
    print(f"QC table: {args.qc_out}")
    print(f"Duplicate-resolution table: {duplicate_path}")


if __name__ == "__main__":
    main()
