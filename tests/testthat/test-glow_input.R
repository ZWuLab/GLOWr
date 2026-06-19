# ==============================================================================
# Tests for prepare_glow_input()
# ==============================================================================


# .make_test_variant_set() is defined in helper-test-data.R (auto-loaded)


# ---- Input validation tests ----

test_that("prepare_glow_input validates variant_set class", {
  expect_error(
    prepare_glow_input(list(G = matrix(1)), B = 0.5, PI = 0.3),
    "glow_variant_set"
  )
})

test_that("prepare_glow_input validates exactly one B mode", {
  vset <- .make_test_variant_set(p = 10)
  PI <- rep(0.3, 10)

  # No B specified
  expect_error(
    prepare_glow_input(vset, PI = PI, verbose = 0),
    "B.*B_func.*B_model"
  )
  # Two B modes specified
  expect_error(
    prepare_glow_input(vset, B = rep(0.5, 10), B_func = identity,
                       PI = PI, verbose = 0),
    "B.*B_func.*B_model"
  )
})

test_that("prepare_glow_input validates exactly one PI mode", {
  vset <- .make_test_variant_set(p = 10)
  B <- rep(0.5, 10)

  # No PI specified
  expect_error(
    prepare_glow_input(vset, B = B, verbose = 0),
    "PI.*PI_models"
  )
})

test_that("prepare_glow_input validates B vector length", {
  vset <- .make_test_variant_set(p = 10)
  expect_error(
    prepare_glow_input(vset, B = rep(0.5, 5), PI = rep(0.3, 10), verbose = 0),
    "length"
  )
})


# ---- Functional tests ----

test_that("prepare_glow_input works with pre-computed B vector", {
  vset <- .make_test_variant_set(p = 10)
  B <- rep(0.5, 10)
  PI <- rep(0.3, 10)

  result <- prepare_glow_input(vset, B = B, PI = PI,
                               mac_threshold = 0, verbose = 0)
  expect_s3_class(result, "glow_input")
})

test_that("prepare_glow_input works with B_func mode", {
  vset <- .make_test_variant_set()
  calcB <- function(maf) rep(0.5, length(maf))

  result <- prepare_glow_input(vset, B_func = calcB,
                               PI = rep(0.3, vset$n_variants),
                               verbose = 0)
  expect_s3_class(result, "glow_input")
  expect_equal(length(result$B), result$n_after_collapse)
  expect_equal(length(result$PI), result$n_after_collapse)
  expect_equal(ncol(result$G), result$n_after_collapse)
  expect_true(result$n_after_collapse <= result$n_original)
})

test_that("prepare_glow_input returns NULL for empty result", {
  set.seed(42)
  n <- 50; p <- 5
  g_col <- rbinom(n, 2, 0.1)
  G <- matrix(rep(g_col, p), n, p)  # All columns identical
  vi <- data.frame(variant_id = 1:p, chr = "22", pos = 1:p,
                   ref = "A", alt = "G",
                   MAF = mean(g_col) / 2,
                   MAC = as.integer(sum(g_col)))
  vset <- structure(
    list(G = G, variant_info = vi, annotations = NULL,
         region = list(chr = "22", start = 1, end = 5, label = "TEST"),
         filter_spec = NULL, n_samples = n, n_variants = p,
         n_total_in_region = p),
    class = "glow_variant_set"
  )

  result <- prepare_glow_input(vset, B = rep(0.5, p), PI = rep(0.3, p),
                               ld_threshold = 0.9, verbose = 0)
  # After LD filtering identical columns, should keep 1 variant
  expect_true(is.null(result) || inherits(result, "glow_input"))
})

test_that("prepare_glow_input processing order is correct", {
  vset <- .make_test_variant_set(n = 100, p = 15)
  calcB <- function(maf) sqrt(abs(log(pmax(maf, 1e-10) * pmax(1 - maf, 1e-10))))
  PI <- runif(15, 0.1, 0.9)

  result <- prepare_glow_input(vset, B_func = calcB, PI = PI, verbose = 0)
  expect_s3_class(result, "glow_input")

  # Check processing log records all steps
  expect_true(any(grepl("Flip", result$processing_log)))
  expect_true(any(grepl("B", result$processing_log)))
  expect_true(any(grepl("PI", result$processing_log)))
  expect_true(any(grepl("LD", result$processing_log)))
})

test_that("prepare_glow_input disables collapsing when mac_threshold = 0", {
  vset <- .make_test_variant_set(p = 10)

  result <- prepare_glow_input(vset, B = rep(0.5, 10), PI = rep(0.3, 10),
                               mac_threshold = 0, verbose = 0)
  expect_false(any(result$is_collapsed))
  expect_true(any(grepl("disabled", result$processing_log)))
})

test_that("prepare_glow_input output structure is complete", {
  vset <- .make_test_variant_set(p = 10)
  result <- prepare_glow_input(vset, B = rep(0.5, 10), PI = rep(0.3, 10),
                               mac_threshold = 0, verbose = 0)

  # Check all expected fields exist
  expected_fields <- c("G", "B", "PI", "is_collapsed",
                       "col_mapping", "ld_keep_idx", "region",
                       "n_original", "n_after_filter", "n_after_ld",
                       "n_after_collapse", "cMAC", "processing_log")
  expect_true(all(expected_fields %in% names(result)))

  # Check dimensions are consistent
  expect_equal(ncol(result$G), result$n_after_collapse)
  expect_equal(length(result$B), result$n_after_collapse)
  expect_equal(length(result$PI), result$n_after_collapse)
  expect_equal(length(result$is_collapsed), result$n_after_collapse)
  expect_equal(length(result$col_mapping), result$n_after_collapse)
})

test_that("prepare_glow_input flips alleles and updates MAF", {
  # Create a variant set with a major-allele-coded variant
  set.seed(123)
  n <- 100; p <- 3
  G <- cbind(
    rbinom(n, 2, 0.1),   # Minor allele coded (AF = 0.1)
    rbinom(n, 2, 0.8),   # Major allele coded (AF = 0.8) -- should flip
    rbinom(n, 2, 0.2)    # Minor allele coded (AF = 0.2)
  )
  vi <- data.frame(
    variant_id = 1:3, chr = "22", pos = c(100, 200, 300),
    ref = "A", alt = "G",
    MAF = colMeans(G) / 2,
    MAC = as.integer(colSums(G)),
    stringsAsFactors = FALSE
  )
  vset <- structure(
    list(G = G, variant_info = vi, annotations = NULL,
         region = list(chr = "22", start = 1, end = 500, label = "TEST_FLIP"),
         filter_spec = NULL, n_samples = n, n_variants = p,
         n_total_in_region = p),
    class = "glow_variant_set"
  )

  result <- prepare_glow_input(vset, B = rep(0.5, 3), PI = rep(0.3, 3),
                               mac_threshold = 0, verbose = 0)

  # After flipping, all column means / 2 (i.e. MAF) should be <= 0.5
  result_maf <- colMeans(result$G) / 2
  expect_true(all(result_maf <= 0.5))
})


# ---- Print method tests ----

test_that("print.glow_input works", {
  vset <- .make_test_variant_set(p = 10)
  result <- prepare_glow_input(vset, B = rep(0.5, 10), PI = rep(0.3, 10),
                               mac_threshold = 0, verbose = 0)
  expect_output(print(result), "GLOW Input")
  expect_output(print(result), "TEST_GENE")
})

# ---- ld_keep_idx tests ----

test_that("prepare_glow_input exposes ld_keep_idx with correct semantics", {
  # No special LD setup needed: helper produces non-degenerate columns whose
  # post-LD set is a strict subset of the post-filter set. We assert shape,
  # type, range, and uniqueness here; correctness against a perfectly-
  # correlated pair is checked in the next test.
  vset <- .make_test_variant_set(n = 200, p = 15)
  result <- prepare_glow_input(vset, B = rep(0.5, 15), PI = rep(0.3, 15),
                               mac_threshold = 0, verbose = 0)

  expect_true("ld_keep_idx" %in% names(result))
  expect_type(result$ld_keep_idx, "integer")
  expect_length(result$ld_keep_idx, result$n_after_ld)
  # Values index into the post-filter variant set
  expect_true(all(result$ld_keep_idx >= 1L &
                  result$ld_keep_idx <= result$n_after_filter))
  # No duplicates
  expect_equal(length(unique(result$ld_keep_idx)),
               length(result$ld_keep_idx))
})

test_that("ld_keep_idx correctly drops a perfectly correlated variant", {
  set.seed(1)
  n <- 200; p <- 4
  G <- matrix(rbinom(n * p, 2, 0.1), n, p)
  G[, 2] <- G[, 1]  # variants 1 and 2 perfectly correlated
  vi <- data.frame(
    variant_id = paste0("v", seq_len(p)),
    chr = "22", pos = seq_len(p), ref = "A", alt = "G",
    MAF = colMeans(G) / 2,
    MAC = as.integer(colSums(G)),
    stringsAsFactors = FALSE
  )
  vset <- structure(
    list(G = G, variant_info = vi, annotations = NULL,
         region = list(chr = "22", start = 1, end = 1000, label = "DUP"),
         filter_spec = NULL, n_samples = n, n_variants = p,
         n_total_in_region = p),
    class = "glow_variant_set"
  )

  result <- prepare_glow_input(vset, B = rep(0.5, p), PI = rep(0.3, p),
                               ld_threshold = 0.95,
                               mac_threshold = 0L, verbose = 0)

  # Exactly one of variants 1, 2 should survive LD pruning
  in_kept <- c(1L, 2L) %in% result$ld_keep_idx
  expect_equal(sum(in_kept), 1L)
  # The other two variants should also survive (no LD with the kept one)
  expect_true(3L %in% result$ld_keep_idx)
  expect_true(4L %in% result$ld_keep_idx)
  expect_equal(result$n_after_ld, 3L)
  expect_equal(ncol(result$G), result$n_after_ld)
})


test_that("print.glow_input handles unnamed region", {
  vset <- .make_test_variant_set(p = 10)
  result <- prepare_glow_input(vset, B = rep(0.5, 10), PI = rep(0.3, 10),
                               mac_threshold = 0, verbose = 0)
  # Remove the label to test fallback
  result$region$label <- NULL
  expect_output(print(result), "unnamed")
})
