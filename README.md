# CovalentScope-QM: a source-held-out quantum-AI benchmark

## Scientific question

**Do tautomer-state-aware semiempirical quantum descriptors improve the
prediction of electrophile-enriched inhibitor activity beyond strong RDKit
molecular baselines when the test set is separated by experimental source?**

This is deliberately a benchmark, not a claim that a substructure alone proves
covalent binding.  Molecules are called *electrophile-enriched* when their
SMILES contains a predefined, plausible covalent-warhead pattern.  A paper may
call a molecule a covalent inhibitor only after the original source has been
checked for evidence of covalent mechanism.

## Why this is stronger than the KRAS pilot

* Six candidate kinase/oncology targets rather than a single target: EGFR,
  ERBB2, BTK, JAK3, FGFR4, and KRAS(G12C) pending assay-level confirmation.
* Strict, endpoint-harmonised labels: equality IC50 records in nM only.
* Borderline compounds are excluded: active pChEMBL >= 6; inactive pChEMBL < 5.
* Duplicate compound-target records are collapsed before validation.
* Validation is repeated, stratified, and **document-grouped**.  No source
  document is used in both training and testing in a fold.
* RDKit physicochemical, Morgan-fingerprint, quantum-only, and combined
  models are compared.  The conclusion is based on the comparison—not on a
  desired outcome.

## Project structure

```
config/targets.csv                  target definitions
scripts/fetch_multitarget_chembl.py provenance-preserving ChEMBL retrieval
scripts/build_electrophile_benchmark.py endpoint/QC/warhead filtering
scripts/fetch_source_metadata.py        literature-review sheet for source validation
scripts/fetch_assay_metadata.py         assay/variant review sheet for retained records
scripts/enumerate_states.py          bounded tautomer-state 3D generation
scripts/run_xtb_states.py            fault-tolerant parallel xTB execution
scripts/parse_xtb_states.py          descriptor extraction
scripts/aggregate_features.py        RDKit + state-aware quantum table
scripts/benchmark_repeated.py        repeated grouped validation
scripts/make_figures.py              manuscript-ready benchmark figures
```

## Installation

```bash
conda env create -f environment.yml
conda activate covalentscope-qm
```

Or, if you already have an environment containing `requests`, RDKit,
scikit-learn, matplotlib, pandas, and xTB, activate that environment instead.

## Run the study in stages

### 1. Retrieve activity records and create the quality-controlled benchmark

```bash
python scripts/fetch_multitarget_chembl.py \
  --targets config/targets.csv \
  --out data/raw/multitarget_activity.csv \
  --max-pages 30

python scripts/build_electrophile_benchmark.py \
  data/raw/multitarget_activity.csv \
  --out data/processed/electrophile_ic50_benchmark.csv \
  --qc-out results/qc/benchmark_qc.csv \
  --min-per-class 25
```

Review `results/qc/benchmark_qc.csv` first. A target with fewer than 25 active
and 25 inactive compounds is not taken forward. The build script also creates
`results/qc/electrophile_evidence_queue.csv`, which lists source documents and
warhead annotations for manual source checking.

Create the literature-review sheet and check the original papers before using
the term *covalent inhibitor* in a figure or manuscript:

```bash
python scripts/fetch_source_metadata.py \
  results/qc/electrophile_evidence_queue.csv \
  --out results/qc/electrophile_source_review.csv
```

### If fewer than two targets are eligible

Do not start xTB. First resolve additional human kinase targets from the
reviewable query list, then replace `config/targets.csv` with the retained
target IDs and repeat Stage 1:

```bash
python scripts/search_chembl_targets.py \
  --queries config/target_search_queries.csv \
  --out results/qc/target_candidates_for_review.csv
```

In `target_candidates_for_review.csv`, mark only the intended human
single-protein target candidates as `retain`. Keep the original target ID and
source data in the manuscript supplement.

For the retained benchmark, create assay-level metadata before calling a
mechanism or KRAS variant in the manuscript:

```bash
python scripts/fetch_assay_metadata.py \
  data/processed/electrophile_ic50_benchmark.csv \
  --out results/qc/assay_review.csv
```

### 2. Generate and calculate a bounded quantum-state ensemble

Start with four tautomer states per molecule. Do **not** run this stage until
the QC table is acceptable.

```bash
python scripts/enumerate_states.py \
  data/processed/electrophile_ic50_benchmark.csv \
  --out data/processed/tautomer_states.csv \
  --xyz-dir data/quantum/xyz \
  --max-states 4

python scripts/run_xtb_states.py \
  data/processed/tautomer_states.csv \
  --outdir data/quantum/results \
  --workers 3 \
  --timeout 420

python scripts/parse_xtb_states.py \
  data/processed/tautomer_states.csv \
  --results data/quantum/results \
  --out data/processed/quantum_state_descriptors.csv

python scripts/aggregate_features.py \
  data/processed/electrophile_ic50_benchmark.csv \
  data/processed/quantum_state_descriptors.csv \
  --out data/processed/covalentscope_ai_table.csv
```

### 3. Perform the primary validation and make figures

```bash
python scripts/benchmark_repeated.py \
  data/processed/covalentscope_ai_table.csv \
  --outdir results/primary_validation \
  --group document \
  --seeds 30

python scripts/benchmark_repeated.py \
  data/processed/covalentscope_ai_table.csv \
  --outdir results/scaffold_validation \
  --group scaffold \
  --seeds 30

python scripts/make_figures.py \
  --results results/primary_validation \
  --table data/processed/covalentscope_ai_table.csv \
  --outdir results/figures
```

## Pre-registered interpretation rules

The primary analysis is the repeated document-grouped validation. A quantum
increment is supported only if the **Morgan + quantum** model has a positive
paired improvement in balanced accuracy and MCC relative to Morgan alone,
with a 95% bootstrap interval excluding zero. If not, report the result as a
negative result: the chosen xTB state descriptors did not provide transferable
information beyond a strong molecular baseline.

## What is publishable

The credible paper is a rigorous benchmark of when quantum descriptors help,
or fail to help, across curated electrophile-enriched inhibitor data. Do not
claim that the workflow discovers covalent inhibitors without primary-source
evidence and prospective experimental testing.

Potential journal targets after successful validation: *Journal of
Cheminformatics*, *Digital Discovery*, *Journal of Chemical Information and
Modeling*, or *Journal of Computer-Aided Molecular Design*.
