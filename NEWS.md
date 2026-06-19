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
