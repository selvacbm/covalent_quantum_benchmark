#!/usr/bin/env python3
"""Create manuscript-ready figures without overstating quantum performance."""
import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


PALETTE = {
    "RDKit_physicochemical": "#4C78A8",
    "Morgan": "#59A14F",
    "Quantum": "#F28E2B",
    "Morgan_plus_quantum": "#9C755F",
    "active": "#4C78A8",
    "inactive": "#F28E2B",
}


def read_csv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def grouped_values(rows, key, value, scope="pooled"):
    output = defaultdict(list)
    for row in rows:
        if row.get("scope") == scope:
            output[row[key]].append(float(row[value]))
    return output


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--table", type=Path)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    folds = read_csv(args.results / "repeated_folds.csv")
    increments = read_csv(args.results / "quantum_increment_summary.csv")
    if not folds:
        raise SystemExit("No repeated_folds.csv found or it contains no rows")

    ba = grouped_values(folds, "feature_set", "balanced_accuracy")
    mcc = grouped_values(folds, "feature_set", "mcc")
    labels = [label for label in ["RDKit_physicochemical", "Morgan", "Quantum", "Morgan_plus_quantum"] if label in ba]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), dpi=220)
    for axis, values, ylabel, baseline in [
        (axes[0], ba, "Balanced accuracy", 0.5),
        (axes[1], mcc, "Matthews correlation coefficient", 0.0),
    ]:
        axis.boxplot([values[label] for label in labels], tick_labels=labels, showfliers=False)
        axis.axhline(baseline, color="0.45", linestyle="--", linewidth=1)
        axis.set_ylabel(ylabel)
        axis.tick_params(axis="x", rotation=25)
    fig.suptitle("Repeated scaffold-held-out validation (pooled benchmark)")
    fig.tight_layout()
    fig.savefig(args.outdir / "primary_validation_boxplots.png", bbox_inches="tight")
    plt.close(fig)

    rows = [row for row in increments if row.get("scope") == "pooled"]
    metrics = ["delta_balanced_accuracy", "delta_mcc", "delta_roc_auc"]
    rows_by_metric = {row["metric"]: row for row in rows}
    fig, axis = plt.subplots(figsize=(7.5, 4.4), dpi=220)
    means = [float(rows_by_metric[m]["mean_delta"]) if m in rows_by_metric else np.nan for m in metrics]
    lows = [float(rows_by_metric[m]["bootstrap_ci_low"]) if m in rows_by_metric else np.nan for m in metrics]
    highs = [float(rows_by_metric[m]["bootstrap_ci_high"]) if m in rows_by_metric else np.nan for m in metrics]
    errors = np.array([np.array(means) - np.array(lows), np.array(highs) - np.array(means)])
    axis.errorbar(range(len(metrics)), means, yerr=errors, fmt="o", color="#9C755F", capsize=5)
    axis.axhline(0, color="0.45", linestyle="--", linewidth=1)
    axis.set_xticks(range(len(metrics)), ["Balanced\naccuracy", "MCC", "ROC-AUC"])
    axis.set_ylabel("Morgan + quantum minus Morgan")
    axis.set_title("Paired quantum increment across repeated scaffold-held-out validation")
    fig.tight_layout()
    fig.savefig(args.outdir / "paired_quantum_increment.png", bbox_inches="tight")
    plt.close(fig)

    if args.table:
        table = read_csv(args.table)
        composition = Counter((row["target_label"], row["activity_class"]) for row in table)
        targets = sorted({row["target_label"] for row in table})
        x = np.arange(len(targets))
        active = [composition[(target, "active")] for target in targets]
        inactive = [composition[(target, "inactive")] for target in targets]
        fig, axis = plt.subplots(figsize=(7.2, 4.4), dpi=220)
        axis.bar(x - 0.2, active, 0.4, label="active", color=PALETTE["active"])
        axis.bar(x + 0.2, inactive, 0.4, label="inactive", color=PALETTE["inactive"])
        axis.set_xticks(x, targets)
        axis.set_ylabel("Compounds after QC and quantum-state filtering")
        axis.set_title("Benchmark composition")
        axis.legend(frameon=False)
        fig.tight_layout()
        fig.savefig(args.outdir / "benchmark_composition.png", bbox_inches="tight")
        plt.close(fig)

        fig, axis = plt.subplots(figsize=(7.2, 4.4), dpi=220)
        for label in ["active", "inactive"]:
            members = [row for row in table if row["activity_class"] == label]
            axis.scatter(
                [float(row["quantum_states_success"]) for row in members],
                [float(row["gap_spread_ev"]) for row in members],
                alpha=0.7, s=30, label=label, color=PALETTE[label],
            )
        axis.set_xlabel("Successful quantum states")
        axis.set_ylabel("HOMO-LUMO gap spread across states (eV)")
        axis.set_title("State-aware electronic variability")
        axis.legend(frameon=False)
        fig.tight_layout()
        fig.savefig(args.outdir / "state_variability.png", bbox_inches="tight")
        plt.close(fig)
    print(f"Figures written to {args.outdir}")


if __name__ == "__main__":
    main()
