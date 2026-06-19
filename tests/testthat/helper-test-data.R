# ==============================================================================
# Shared test helpers for GLOW testing pipeline tests
# ==============================================================================
#
# Used by: test-glow_input.R, test-glow_integration.R


#' Create a minimal glow_variant_set for unit testing
#'
#' Generates synthetic genotype data with random MAFs and annotations.
#'
#' @param n Integer. Number of samples.
#' @param p Integer. Number of variants.
#' @param seed Integer. Random seed for reproducibility.
#'
#' @return A \code{glow_variant_set} S3 object with synthetic data.
.make_test_variant_set <- function(n = 100, p = 20, seed = 42) {
  set.seed(seed)
  mafs <- runif(p, 0.01, 0.3)
  G <- matrix(rbinom(n * p, 2, rep(mafs, each = n)), n, p)
  vi <- data.frame(
    variant_id = seq_len(p),
    rsid = sprintf("rs%07d", seq_len(p)),
    chr = rep("22", p),
    pos = sort(sample(1e6:2e6, p)),
    ref = rep("A", p),
    alt = rep("G", p),
    MAF = colMeans(G) / 2,
    MAC = as.integer(colSums(G)),
    stringsAsFactors = FALSE
  )
  # Create annotation matrix using package-internal default PI features
  annot_names <- .default_PI_features()
  annot <- matrix(runif(p * length(annot_names)), p, length(annot_names),
                  dimnames = list(NULL, annot_names))

  structure(
    list(G = G, variant_info = vi, annotations = annot,
         region = list(chr = "22", start = 1e6, end = 2e6, label = "TEST_GENE"),
         filter_spec = NULL, n_samples = n, n_variants = p,
         n_total_in_region = p),
    class = "glow_variant_set"
  )
}
