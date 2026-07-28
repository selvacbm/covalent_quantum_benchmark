#!/usr/bin/env python3
"""Repeated source- or scaffold-grouped validation for CovalentScope-QM.

The primary analysis uses document groups.  A group is never shared by the
training and test sets in an individual fold.  The same valid folds are used
for every feature set, making paired model comparisons possible.
"""
import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score, f1_score, matthews_corrcoef, roc_auc_score
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

PHYSICOCHEMICAL = [
    "molecular_weight", "tpsa", "logp", "hbond_donors", "hbond_acceptors",
    "rotatable_bonds", "aromatic_rings", "fraction_csp3", "heavy_atom_count",
]
QUANTUM = ["mean_gap_ev", "min_gap_ev", "max_gap_ev", "gap_spread_ev", "gap_std_ev", "quantum_states_success"]


def as_float_matrix(rows, columns):
    return np.array([[float(row[column]) for column in columns] for row in rows], dtype=float)


def morgan_matrix(rows, radius=2, n_bits=1024):
    output = np.zeros((len(rows), n_bits), dtype=np.float32)
    for index, row in enumerate(rows):
        molecule = Chem.MolFromSmiles(row["canonical_smiles"])
        if molecule is None:
            raise ValueError(f"Invalid SMILES in modelling table: {row['record_id']}")
        fingerprint = AllChem.GetMorganFingerprintAsBitVect(molecule, radius, nBits=n_bits)
        DataStructs.ConvertToNumpyArray(fingerprint, output[index])
    return output


def target_one_hot(rows):
    labels = sorted({row["target_label"] for row in rows})
    positions = {label: index for index, label in enumerate(labels)}
    matrix = np.zeros((len(rows), len(labels)), dtype=float)
    for index, row in enumerate(rows):
        matrix[index, positions[row["target_label"]]] = 1.0
    return matrix


def bootstrap_interval(values, rng, draws=5000):
    values = np.asarray(values, dtype=float)
    if not len(values):
        return np.nan, np.nan
    samples = rng.choice(values, size=(draws, len(values)), replace=True).mean(axis=1)
    return float(np.percentile(samples, 2.5)), float(np.percentile(samples, 97.5))


def safe_auc(y_true, probability):
    return roc_auc_score(y_true, probability) if len(set(y_true)) == 2 else np.nan


def choose_groups(rows, group):
    if group == "document":
        return np.array([row.get("document_chembl_id") or f"missing_doc__{row['record_id']}" for row in rows])
    if group == "scaffold":
        return np.array([row.get("murcko_scaffold") or row["canonical_smiles"] for row in rows])
    raise ValueError("group must be document or scaffold")


def make_rows(table, scope):
    rows = []
    for row in csv.DictReader(table.open()):
        if row.get("activity_class") not in {"active", "inactive"}:
            continue
        if scope != "pooled" and row.get("target_label") != scope:
            continue
        try:
            int(row["y"])
            for value in PHYSICOCHEMICAL + QUANTUM:
                float(row[value])
        except (KeyError, ValueError):
            continue
        rows.append(row)
    return rows


def evaluate_scope(rows, scope, group_type, seeds, n_splits, output, oof_output, rng):
    y = np.array([int(row["y"]) for row in rows], dtype=int)
    groups = choose_groups(rows, group_type)
    if len(set(y)) != 2 or len(set(groups)) < 3:
        print(f"SKIP {scope}: requires both classes and at least three {group_type} groups")
        return
    local_splits = min(n_splits, len(set(groups)))
    matrices = {
        "RDKit_physicochemical": np.hstack([as_float_matrix(rows, PHYSICOCHEMICAL), target_one_hot(rows)]),
        "Morgan": np.hstack([morgan_matrix(rows), target_one_hot(rows)]),
        "Quantum": np.hstack([as_float_matrix(rows, QUANTUM), target_one_hot(rows)]),
        "Morgan_plus_quantum": np.hstack([morgan_matrix(rows), as_float_matrix(rows, QUANTUM), target_one_hot(rows)]),
    }
    for seed in range(seeds):
        cv = StratifiedGroupKFold(n_splits=local_splits, shuffle=True, random_state=seed + 1)
        for fold, (train, test) in enumerate(cv.split(matrices["Morgan"], y, groups), 1):
            # A fold without both labels has undefined AUC and can produce an
            # artificial balanced-accuracy result. It is omitted for all sets.
            if len(set(y[train])) < 2 or len(set(y[test])) < 2:
                continue
            for feature_set, matrix in matrices.items():
                model = make_pipeline(
                    StandardScaler(),
                    LogisticRegression(class_weight="balanced", max_iter=5000, random_state=42),
                )
                model.fit(matrix[train], y[train])
                prediction = model.predict(matrix[test])
                probability = model.predict_proba(matrix[test])[:, 1]
                result = {
                    "scope": scope,
                    "group_type": group_type,
                    "seed": seed + 1,
                    "fold": fold,
                    "feature_set": feature_set,
                    "n_train": len(train),
                    "n_test": len(test),
                    "n_train_groups": len(set(groups[train])),
                    "n_test_groups": len(set(groups[test])),
                    "balanced_accuracy": balanced_accuracy_score(y[test], prediction),
                    "macro_f1": f1_score(y[test], prediction, average="macro"),
                    "mcc": matthews_corrcoef(y[test], prediction),
                    "roc_auc": safe_auc(y[test], probability),
                }
                output.append(result)
                for local_index, row_index in enumerate(test):
                    row = rows[row_index]
                    oof_output.append({
                        "scope": scope,
                        "group_type": group_type,
                        "seed": seed + 1,
                        "fold": fold,
                        "feature_set": feature_set,
                        "record_id": row["record_id"],
                        "target_label": row["target_label"],
                        "document_chembl_id": row.get("document_chembl_id", ""),
                        "activity_class": row["activity_class"],
                        "y_true": int(y[row_index]),
                        "y_pred": int(prediction[local_index]),
                        "y_probability": float(probability[local_index]),
                        "gap_spread_ev": row["gap_spread_ev"],
                    })


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--group", choices=["document", "scaffold"], default="document")
    parser.add_argument("--seeds", type=int, default=30)
    parser.add_argument("--splits", type=int, default=3)
    parser.add_argument("--scope", choices=["pooled", "per_target", "both"], default="both")
    args = parser.parse_args()
    if args.seeds < 1 or args.splits < 2:
        raise SystemExit("--seeds must be >= 1 and --splits must be >= 2")
    args.outdir.mkdir(parents=True, exist_ok=True)
    all_rows = make_rows(args.input, "pooled")
    if not all_rows:
        raise SystemExit("No usable labelled rows in input table")
    scopes = []
    if args.scope in {"pooled", "both"}:
        scopes.append("pooled")
    if args.scope in {"per_target", "both"}:
        scopes.extend(sorted({row["target_label"] for row in all_rows}))

    results, oof = [], []
    rng = np.random.default_rng(20260723)
    for scope in scopes:
        rows = all_rows if scope == "pooled" else make_rows(args.input, scope)
        evaluate_scope(rows, scope, args.group, args.seeds, args.splits, results, oof, rng)
    if not results:
        raise SystemExit("No valid grouped folds were created. Inspect group and class counts.")

    result_fields = list(results[0])
    with (args.outdir / "repeated_folds.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=result_fields)
        writer.writeheader()
        writer.writerows(results)
    with (args.outdir / "oof_predictions.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(oof[0]))
        writer.writeheader()
        writer.writerows(oof)

    summary, paired = [], []
    by_scope_set = defaultdict(list)
    for row in results:
        by_scope_set[(row["scope"], row["feature_set"])].append(row)
    for (scope, feature_set), members in sorted(by_scope_set.items()):
        summary.append({
            "scope": scope,
            "feature_set": feature_set,
            "valid_folds": len(members),
            "balanced_accuracy_mean": float(np.mean([row["balanced_accuracy"] for row in members])),
            "balanced_accuracy_sd": float(np.std([row["balanced_accuracy"] for row in members], ddof=0)),
            "mcc_mean": float(np.mean([row["mcc"] for row in members])),
            "mcc_sd": float(np.std([row["mcc"] for row in members], ddof=0)),
            "roc_auc_mean": float(np.nanmean([row["roc_auc"] for row in members])),
            "roc_auc_sd": float(np.nanstd([row["roc_auc"] for row in members], ddof=0)),
        })
    with (args.outdir / "repeated_summary.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary[0]))
        writer.writeheader()
        writer.writerows(summary)

    # Paired increments are calculated at the repeat level, not by treating
    # all folds as independent observations.
    lookup = {(row["scope"], row["seed"], row["fold"], row["feature_set"]): row for row in results}
    keys = sorted({key[:3] for key in lookup})
    by_scope_metric = defaultdict(list)
    for scope, seed, fold in keys:
        baseline = lookup.get((scope, seed, fold, "Morgan"))
        combined = lookup.get((scope, seed, fold, "Morgan_plus_quantum"))
        if not baseline or not combined:
            continue
        paired.append({
            "scope": scope,
            "seed": seed,
            "fold": fold,
            "delta_balanced_accuracy": combined["balanced_accuracy"] - baseline["balanced_accuracy"],
            "delta_mcc": combined["mcc"] - baseline["mcc"],
            "delta_roc_auc": combined["roc_auc"] - baseline["roc_auc"],
        })
        by_scope_metric[scope].append(paired[-1])
    with (args.outdir / "paired_quantum_increment.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(paired[0]) if paired else ["scope"])
        writer.writeheader()
        writer.writerows(paired)

    increment_summary = []
    for scope, members in sorted(by_scope_metric.items()):
        # First average fold values within a repeat; then bootstrap repeats.
        by_seed = defaultdict(list)
        for row in members:
            by_seed[row["seed"]].append(row)
        for metric in ("delta_balanced_accuracy", "delta_mcc", "delta_roc_auc"):
            values = [np.mean([row[metric] for row in rows]) for rows in by_seed.values()]
            lo, hi = bootstrap_interval(values, rng)
            increment_summary.append({
                "scope": scope,
                "metric": metric,
                "n_repeats": len(values),
                "mean_delta": float(np.mean(values)),
                "bootstrap_ci_low": lo,
                "bootstrap_ci_high": hi,
                "quantum_increment_supported": int(lo > 0),
            })
    with (args.outdir / "quantum_increment_summary.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(increment_summary[0]) if increment_summary else ["scope"])
        writer.writeheader()
        writer.writerows(increment_summary)

    print(f"Wrote {len(results)} valid grouped-fold evaluations to {args.outdir}")
    for row in summary:
        if row["scope"] == "pooled":
            print(f"{row['feature_set']}: BA={row['balanced_accuracy_mean']:.3f}; MCC={row['mcc_mean']:.3f}")


if __name__ == "__main__":
    main()
