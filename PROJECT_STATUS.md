# CovalentScope-QM project status

## Project started

The project is configured as a three-target, source-held-out benchmark for
electrophile-enriched inhibitor activity: KRAS(G12C), BTK, and EGFR.

## Completed framework

- Provenance-preserving ChEMBL retrieval.
- Strict equality-IC50/nM filtering and pChEMBL class thresholds.
- Predefined electrophile pattern annotation.
- Removal of contradictory activity labels for a target-molecule pair.
- Bounded tautomer-state generation and fault-tolerant xTB calculations.
- RDKit, Morgan, quantum-only, and Morgan-plus-quantum model comparisons.
- Repeated document-grouped primary validation and scaffold-grouped sensitivity analysis.
- Paired bootstrap test of the incremental value of quantum descriptors.
- Publication-quality, conclusion-safe figure generation.

## Next execution checkpoint

Run the two Stage 1 commands in the README. Inspect the QC table before
launching xTB calculations. The project should proceed only if at least two
targets have at least 25 active and 25 inactive eligible molecules, with a
reasonable number of independent source documents in each class.

## Claim boundary

This study tests prediction of activity in an electrophile-enriched benchmark.
It does not establish covalent binding, a new inhibitor, or clinical utility.
Those claims require manual primary-source evidence and prospective experiment.
