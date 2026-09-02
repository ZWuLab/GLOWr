# GLOWr 0.1.1

- Fixed a multi-chunk bug in `annotate_favor()`'s xsv fast path (exact
  matching): variants falling in the second and later FAVOR chunk files of a
  chromosome silently lost their annotations (the per-chunk left-joins made the
  downstream deduplication keep an empty row instead of the real match). The
  per-chunk joins are now inner joins, and every `xsv` invocation's exit status
  is checked, so a failed join aborts instead of writing a silently
  under-annotated aGDS.
- `annotate_favor()`'s default `features` now covers the full FAVOR Essential DB
  annotation content (30 columns, previously the 20 purely numeric ones), so an
  aGDS built with the defaults carries the categorical features (e.g. GENCODE
  categories) that the built-in coding masks filter on. String features are
  first-class on the R join path: flexible-match aggregation is type-aware
  (numeric features average across matches; character features take the first
  non-empty value), and `na_handling = "zero"` no longer writes `0` into
  character columns.
- New `annotation_predicate()` filter grammar for `variant_filter()`: a clause
  condition may now test `nonempty`/`empty`, `in`/`not_in`, or numeric
  `gt`/`ge`/`lt`/`le` alongside the existing set-membership conditions. Missing
  annotation values never satisfy a positive test.
- `calibrate_pvalues()` gains `inflation_only`: when `TRUE`, a calibration
  factor below 1 is clamped to 1, reproducing the field-standard one-sided
  genomic-control convention (correct inflation, leave deflation alone).
  Default `FALSE` keeps the previous two-sided behavior.
- The compiled code now links BLAS/LAPACK explicitly (`src/Makevars` and
  `src/Makevars.win`); the Rtools toolchain does not auto-link them, so this
  fixes the Windows `R CMD check`.
- Documentation: `PI` is described throughout as the variant-importance score,
  and the methodology reference now cites the GLOW methods paper (Zhang, Liu,
  Landers, and Wu, Annals of Applied Statistics, in revision).

# GLOWr 0.1.0

Initial public release.

- GLOW (inteGrative anaLysis using Optimized Weights) variant-set association tests:
  Burden, SKAT, Fisher-combination, and the Omnibus test that combines them,
  for both continuous and binary phenotypes with covariate adjustment.
- Optimal data-adaptive weighting: estimation of the effect-size model `B` and
  the variant-importance score `PI`, and their combination into per-variant
  weights used by the set-level tests.
- Single-variant (per-variant) score statistics with saddlepoint (SPA)
  calibration for binary traits.
- p-value calibration and genomic-inflation diagnostics, including LD-score
  regression utilities and Cauchy-combination (CCT) aggregation.
- Region definition and variant-set extraction from GDS/aGDS inputs, with FAVOR
  functional-annotation support.
- Built on the GFisher backend for the generalized Fisher-type combination
  p-values.
