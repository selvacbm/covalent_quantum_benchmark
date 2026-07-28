# CovalentScope-QM study protocol

## Objective

Test whether state-aware GFN2-xTB electronic descriptors add transferable
predictive value beyond RDKit descriptors and Morgan fingerprints for
electrophile-enriched inhibitor activity across KRAS(G12C), BTK, and EGFR.

## Primary hypothesis

The Morgan-plus-quantum model will improve both balanced accuracy and Matthews
correlation coefficient relative to Morgan alone under repeated,
document-grouped validation.  The hypothesis is supported only when the 95%
bootstrap confidence interval for each paired increment excludes zero.

## Dataset rules

1. Keep only equality IC50 records reported in nM.
2. Define active as pChEMBL >= 6 and inactive as pChEMBL < 5; exclude the
   5 to <6 borderline interval.
3. Require a predefined plausible electrophile SMARTS pattern and molecular
   weight <=900 Da.
4. Exclude target-molecule pairs with contradictory active and inactive calls.
5. Keep one median-potency representative record per target-molecule pair.
6. Do not call a compound covalent solely from the SMARTS pattern. Confirm
   covalent mechanism in the cited primary source before making that claim.

## Quantum calculation

Generate up to four bounded tautomer states.  Use RDKit/MMFF conformers and
GFN2-xTB single-point calculations.  Retain only successful xTB outputs and
report the calculation mode exactly as recorded in the dataset.  The quantum
feature block is the state distribution of HOMO-LUMO gaps, not raw total
energy.

## Validation and interpretation

The primary validation is repeated, stratified, document-grouped cross
validation.  Scaffold-grouped validation is a secondary sensitivity analysis.
The models are RDKit physicochemical, Morgan, quantum-only, and
Morgan-plus-quantum.  Report a negative result if the paired quantum increment
does not robustly exceed zero.

## Minimum manuscript figures

1. Data curation and state-ensemble workflow.
2. Per-target active/inactive benchmark composition.
3. Repeated document-grouped performance distributions.
4. Paired Morgan-plus-quantum versus Morgan increments with 95% intervals.
5. State variability distribution as a mechanistic/uncertainty analysis.
