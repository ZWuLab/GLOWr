
# GLOWr `inst/extdata`

Genuine package assets referenced at runtime / used as worked examples and test fixtures.

## Files

* `FAVORdatabase_chrsplit.csv` — FAVOR database chromosome-split index; loaded via
  `system.file("extdata", ...)` in `R/favor_annotator.R`. Load-bearing.
* `ALS-known-SNPs-raw.xlsx` — raw known-ALS SNPs (GWAS Catalog / Open Target Genetics). The canonical
  **worked example** for the PI/B loaders: referenced in `prepare_PI_case_data()` `@examples` and via
  `system.file()` in the `get_PI_data_loader` tests.

## Bundled B-training datasets (provenance in the source studies)

The lazy-loaded package datasets `ALS_snvs_B_training` and `BMD_snvs_B_training` (in `data/`, via
`data(...)`) are documented in `R/data.R`. Their raw inputs + generator scripts live with the studies:

* ALS — raw input is `ALS-known-SNPs-raw.xlsx` (above); prepared by the study's
  B-training data-preparation scripts.
* BMD — prepared from BMD/osteoporosis literature summary statistics (raw inputs
  maintained with the source study).
