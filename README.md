# GLOWr: inteGrative anaLysis using Optimized Weights (GLOW) for Rare Variant Analysis

## Overview

**GLOWr** implements the GLOW (inteGrative anaLysis using Optimized Weights) methodology for
integrative variant-set association analysis with optimal weighting in
whole-genome sequencing (WGS) studies. It provides Burden, SKAT, Fisher, and
Omnibus variant-set tests with data-adaptive optimal weights, supporting both
continuous and binary phenotypes with covariate adjustment, and scales to
genome-wide analysis.

GLOWr is the methods-and-estimation layer of the GLOW software family: it
implements the single-region/single-variant tests, the optimal-weight training
(the effect-size model `B` and the variant-importance score `PI`), and the I/O
for GDS/aGDS inputs. It builds on the
[GFisher](https://github.com/ZWuLab/GFisher) package for accurate generalized
Fisher-type combination p-values.

## Key Features

- **GLOW variant-set tests**: Burden, SKAT, Fisher-combination, and the Omnibus
  test (`GLOW_Omni`) that combines them via the Cauchy Combination Test (CCT),
  for continuous and binary phenotypes with covariate adjustment.
- **Optimal data-adaptive weighting**: estimation of the allelic effect-size
  model `B` (`get_B`) and the variant-importance score `PI` (`get_PI`), combined
  into the per-variant weights used by the set-level tests
  (`Optimal_Weights_M`).
- **Single-variant score statistics**: per-variant score tests
  (`getZ_marg_score`) with saddlepoint (SPA) calibration for binary traits
  (`getZ_marg_score_binary_SPA`) and a genome-wide single-variant scan
  (`marginal_scan`).
- **Calibration and diagnostics**: p-value calibration
  (`calibrate_pvalues`), genomic-inflation estimation
  (`estimate_inflation_factor`, `compute_lambda_gc`), and LD-score regression
  utilities (`ldsc_regression`, `compute_ld_scores`).
- **Region and variant-set handling**: region definition
  (`define_regions_gene`, `define_regions_window`, `define_regions_custom`),
  variant-set extraction from GDS/aGDS (`extract_variant_set`), and FAVOR
  functional-annotation support (`annotate_favor`).
- **Built on GFisher**: uses the GFisher backend for the generalized Fisher-type
  combination tests.

## Installation

### From GitHub (recommended)

```r
# install.packages("remotes")
remotes::install_github("ZWuLab/GLOWr")
```

GLOWr depends on **GFisher** (`ZWuLab/GFisher`). The GitHub install resolves it
automatically via the `Remotes:` field in `DESCRIPTION`; no separate step is
needed. For a from-source install you may need to install GFisher first:

```r
remotes::install_github("ZWuLab/GFisher")
remotes::install_github("ZWuLab/GLOWr")
```

### System Requirements

- R >= 3.5.0
- A C++ compiler (the package compiles RcppArmadillo code)

### Optional Dependencies

Several capabilities rely on Suggested packages installed on demand:

- **GDS/aGDS I/O and annotation**: `gdsfmt`, `SeqArray`, `SNPRelate`
  (Bioconductor), and `STAAR` for the FAVOR/aGDS paths.
- **Reading the bundled ALS reference SNP list**: `readxl`.
- **Some simulations / utilities**: `MASS`, `mvtnorm`.

Install the Bioconductor dependencies with:

```r
# install.packages("BiocManager")
BiocManager::install(c("gdsfmt", "SeqArray", "SNPRelate"))
```

## Quick Start

The set-level tests take per-variant score statistics together with the
optimal-weight inputs `B` (allelic effect sizes) and `PI` (variant-importance
prior). In a real analysis `B` and `PI` are *estimated* from training data with
`get_B()` / `get_PI()` (see the vignettes); the toy values below stand in only to
show the call shape.

```r
library(GLOWr)

# Simulate a small variant set (n samples x m variants), covariates, and a
# binary phenotype.
set.seed(123)
n <- 500; m <- 20
G <- matrix(rbinom(n * m, 2, 0.1), n, m)   # genotypes (0/1/2)
X <- matrix(rnorm(n * 2), n, 2)            # covariates
Y <- rbinom(n, 1, 0.3)                     # binary phenotype

# Per-variant score statistics (fits the null model internally).
marg_stats <- getZ_marg_score(G, X, Y, trait = "binary")

# Optimal-weight inputs (normally estimated via get_B() / get_PI()).
B  <- rnorm(m, mean = 0, sd = 0.2)         # allelic effect sizes
PI <- runif(m, 0.1, 0.9)                   # variant-importance score

# GLOW omnibus test (combines Burden / SKAT / Fisher via CCT).
result <- GLOW_Omni(marg_stats, B, PI)

# The final omnibus p-value is the last row of the p-value matrix.
result$PVAL[nrow(result$PVAL), ]
```

For the full optimal-weighting workflow (estimating `B` and `PI`, then running
genome-wide), see the package vignettes (below) and the `GLOWpipeline`
orchestration package.

## Documentation

- **Vignettes**:
  - `vignette("estimating_effect_sizes", package = "GLOWr")` — estimating the
    allelic effect sizes `B` for optimal weighting.
  - `vignette("data-preparation-get-B", package = "GLOWr")` — preparing data for
    `B` estimation.
  - `vignette("B_estimation_diagnostics_example", package = "GLOWr")` —
    diagnostic tools for `B` estimation.
- **Function help**: `?GLOW_Omni`, `?get_B`, `?get_PI`, `?getZ_marg_score`, etc.
- **Package overview**: `?GLOWr`.

## Testing

The package includes a comprehensive `testthat` suite:

```r
devtools::test()
```

## Package Structure

```
GLOWr/
├── R/                  # R function implementations
├── src/                # C++ source (RcppArmadillo)
├── tests/testthat/     # Unit tests
├── vignettes/          # Tutorials and worked examples
├── data/               # Bundled training datasets
└── inst/extdata/       # Reference annotation / SNP-list resources
```

## License

GPL-3

## Citation

If you use GLOWr, please cite the GLOW methodology paper:

> Zhang, H., Liu, M., Landers, J. E., and Wu, Z. Integrated Weighted Association
> Test with Application to Genetic Association Studies. *Annals of Applied
> Statistics* (in revision).

See `inst/CITATION` (or run `citation("GLOWr")`).

## Related software / acknowledgments

GLOWr builds on and interoperates with:

- **FAVOR** — functional annotations used by the aGDS annotation paths
  (`annotate_favor`, `get_PI`). Zhou, H., Arapoglou, T., Li, X., et al. (2023).
  FAVOR: functional annotation of variants online resource and annotator for
  variation across the human genome. *Nucleic Acids Research*, 51(D1),
  D1300-D1311. doi:10.1093/nar/gkac966
- **STAAR / STAARpipeline** — used for the STAAR-O comparison/interop in the
  region tests and for the aGDS/genes_info conventions. Li, X., Li, Z., Zhou, H.,
  et al. (2020). Dynamic incorporation of multiple in silico functional
  annotations empowers rare variant association analysis of large whole-genome
  sequencing studies at scale. *Nature Genetics*, 52, 969-983.
  doi:10.1038/s41588-020-0676-4. Li, Z., Li, X., Zhou, H., et al. (2022). A
  framework for detecting noncoding rare-variant associations of large-scale
  whole-genome sequencing studies. *Nature Methods*, 19, 1599-1611.
  doi:10.1038/s41592-022-01640-x
- **SPA (via SPAtest)** — GLOWr uses the saddlepoint approximation (SPA), through
  the `SPAtest` package, for binary-trait single-variant score tests
  (`getZ_marg_score_binary_SPA`, `marginal_scan`), which calibrates the score
  statistic under case-control imbalance and for rare variants. Dey, R., Schmidt,
  E. M., Abecasis, G. R., and Lee, S. (2017). A fast and accurate algorithm to
  test for binary phenotypes and its application to PheWAS. American Journal of
  Human Genetics, 101(1), 37-49. doi:10.1016/j.ajhg.2017.05.014

## Contact

Maintainer: **Zheyang Wu** (zheyangwu@wpi.edu). Please file bugs and questions at
<https://github.com/ZWuLab/GLOWr/issues>.

## Authors

- **Hong Zhang** (consistencyzhang@gmail.com) — Author
- **Ming Liu** (mliu5@wpi.edu) — Author
- **Zheyang Wu** (zheyangwu@wpi.edu) — Author, Maintainer

> Developed with AI assistance; see [AI-USE.md](AI-USE.md).
