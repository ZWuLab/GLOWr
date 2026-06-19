## Tests for marginal_plot.R
## Functions: plot_manhattan(), plot_qq(), plot.glow_marginal_scan()
##
## All plot tests use pdf(NULL) / dev.off() to suppress graphical output.

# Helper: open a null graphics device to capture plots without displaying
.open_null_dev <- function() {
  grDevices::pdf(file = NULL)
}
.close_null_dev <- function() {
  grDevices::dev.off()
}

# ==============================================================================
# Phase 1: plot_manhattan()
# ==============================================================================

test_that("plot_manhattan renders without error (full 22-chr dataset)", {
  sim_data <- data.frame(
    chr      = rep(1:22, each = 100),
    pos      = rep(seq(1e6, 1e8, length.out = 100), 22),
    pvalue   = runif(2200),
    variant_id = 1:2200,
    stringsAsFactors = FALSE
  )
  class(sim_data) <- c("glow_marginal_scan", "data.frame")

  .open_null_dev()
  on.exit(.close_null_dev(), add = TRUE)
  expect_silent(plot_manhattan(sim_data))
})

test_that("plot_manhattan handles single chromosome", {
  sim_data <- data.frame(
    chr    = rep(22L, 50),
    pos    = seq_len(50) * 1e6,
    pvalue = runif(50),
    stringsAsFactors = FALSE
  )

  .open_null_dev()
  on.exit(.close_null_dev(), add = TRUE)
  expect_silent(plot_manhattan(sim_data))
})

test_that("plot_manhattan handles NA p-values silently", {
  sim_data <- data.frame(
    chr    = rep(1L, 100),
    pos    = seq_len(100),
    pvalue = c(runif(90), rep(NA, 10)),
    stringsAsFactors = FALSE
  )

  .open_null_dev()
  on.exit(.close_null_dev(), add = TRUE)
  expect_silent(plot_manhattan(sim_data))
})

test_that("plot_manhattan highlights specified variants", {
  sim_data <- data.frame(
    chr        = rep(1L, 50),
    pos        = seq_len(50),
    pvalue     = runif(50),
    variant_id = seq_len(50),
    stringsAsFactors = FALSE
  )

  .open_null_dev()
  on.exit(.close_null_dev(), add = TRUE)
  expect_silent(plot_manhattan(sim_data, highlight = c(1L, 5L, 10L)))
})

test_that("plot_manhattan warns when highlight used without variant_id column", {
  sim_data <- data.frame(
    chr    = rep(1L, 30),
    pos    = seq_len(30),
    pvalue = runif(30),
    stringsAsFactors = FALSE
  )

  .open_null_dev()
  on.exit(.close_null_dev(), add = TRUE)
  expect_warning(
    plot_manhattan(sim_data, highlight = c(1L, 2L)),
    "highlight.*ignored"
  )
})

test_that("plot_manhattan errors when required column is missing", {
  bad_data <- data.frame(chromosome = 1:10, pos = 1:10, pvalue = runif(10))
  expect_error(plot_manhattan(bad_data), "not found in x")
})

test_that("plot_manhattan errors when all p-values are NA", {
  bad_data <- data.frame(chr = 1:5, pos = 1:5, pvalue = rep(NA_real_, 5))
  expect_error(plot_manhattan(bad_data), "No valid p-values found")
})

test_that("plot_manhattan errors on invalid significance line threshold", {
  sim_data <- data.frame(chr = rep(1L, 20), pos = 1:20, pvalue = runif(20))
  .open_null_dev()
  on.exit(.close_null_dev(), add = TRUE)
  expect_error(plot_manhattan(sim_data, suggestive_line = 5), "suggestive_line.*p-value")
  expect_error(plot_manhattan(sim_data, genome_wide_line = 10), "genome_wide_line.*p-value")
})

test_that("plot_manhattan handles single variant", {
  single <- data.frame(chr = 1L, pos = 1e6, pvalue = 0.001)
  .open_null_dev()
  on.exit(.close_null_dev(), add = TRUE)
  expect_silent(plot_manhattan(single))
})

test_that("plot_manhattan handles 'chr' prefix in chromosome column", {
  sim_data <- data.frame(
    chr    = paste0("chr", rep(1:5, each = 20)),
    pos    = rep(seq_len(20) * 1e6, 5),
    pvalue = runif(100),
    stringsAsFactors = FALSE
  )

  .open_null_dev()
  on.exit(.close_null_dev(), add = TRUE)
  expect_silent(plot_manhattan(sim_data))
})

test_that("plot_manhattan suppresses significance lines when NULL", {
  sim_data <- data.frame(
    chr    = rep(1L, 50),
    pos    = seq_len(50),
    pvalue = runif(50),
    stringsAsFactors = FALSE
  )

  .open_null_dev()
  on.exit(.close_null_dev(), add = TRUE)
  expect_silent(
    plot_manhattan(sim_data, suggestive_line = NULL, genome_wide_line = NULL)
  )
})


# ==============================================================================
# Phase 2: plot_qq()
# ==============================================================================

test_that("plot_qq returns a numeric lambda_gc", {
  set.seed(1)
  pvals  <- runif(1000)

  .open_null_dev()
  on.exit(.close_null_dev(), add = TRUE)
  lambda <- plot_qq(pvals)

  expect_true(is.numeric(lambda))
  expect_length(lambda, 1)
  expect_true(lambda > 0)
})

test_that("plot_qq lambda near 1.0 for uniform p-values", {
  set.seed(42)
  pvals  <- runif(10000)

  .open_null_dev()
  on.exit(.close_null_dev(), add = TRUE)
  lambda <- plot_qq(pvals, ci = FALSE)

  expect_true(abs(lambda - 1.0) < 0.1)
})

test_that("plot_qq handles data.frame input", {
  df <- data.frame(pvalue = runif(500))
  class(df) <- c("glow_marginal_scan", "data.frame")

  .open_null_dev()
  on.exit(.close_null_dev(), add = TRUE)
  lambda <- plot_qq(df)

  expect_true(is.numeric(lambda))
})

test_that("plot_qq handles alternative pvalue column name", {
  df <- data.frame(pvalue = runif(500), pvalue_SPA = runif(500))

  .open_null_dev()
  on.exit(.close_null_dev(), add = TRUE)
  lambda <- plot_qq(df, pvalue_col = "pvalue_SPA")

  expect_true(is.numeric(lambda))
})

test_that("plot_qq handles NA p-values without error", {
  pvals <- c(runif(90), rep(NA, 10))

  .open_null_dev()
  on.exit(.close_null_dev(), add = TRUE)
  expect_silent(plot_qq(pvals))
})

test_that("plot_qq detects inflation for inflated z-scores", {
  set.seed(42)
  z_inflated <- stats::rnorm(5000) * 1.2
  pvals      <- 2 * stats::pnorm(-abs(z_inflated))

  .open_null_dev()
  on.exit(.close_null_dev(), add = TRUE)
  lambda <- plot_qq(pvals)

  expect_true(lambda > 1.1)
})

test_that("plot_qq renders with ci = TRUE without error", {
  set.seed(7)
  pvals <- runif(300)

  .open_null_dev()
  on.exit(.close_null_dev(), add = TRUE)
  expect_silent(plot_qq(pvals, ci = TRUE))
})

test_that("plot_qq errors on non-data.frame, non-numeric input", {
  expect_error(plot_qq(list(a = 1:10)), "numeric vector or a data.frame")
})

test_that("plot_qq errors when pvalue_col missing from data.frame", {
  df <- data.frame(score = runif(100))
  expect_error(plot_qq(df, pvalue_col = "pvalue"), "not found in x")
})

test_that("plot_qq errors when all p-values are NA", {
  expect_error(plot_qq(rep(NA_real_, 10)), "No valid p-values after filtering")
})


# ==============================================================================
# Phase 3: plot.glow_marginal_scan() S3 dispatch
# ==============================================================================

test_that("plot.glow_marginal_scan dispatches to manhattan by default", {
  sim_data <- data.frame(
    chr    = rep(1:5, each = 100),
    pos    = rep(seq_len(100) * 1e6, 5),
    pvalue = runif(500),
    stringsAsFactors = FALSE
  )
  class(sim_data) <- c("glow_marginal_scan", "data.frame")

  .open_null_dev()
  on.exit(.close_null_dev(), add = TRUE)
  result <- plot(sim_data)
  expect_null(result)
})

test_that("plot.glow_marginal_scan dispatches to qq and returns lambda", {
  sim_data <- data.frame(
    chr    = rep(1:5, each = 100),
    pos    = rep(seq_len(100) * 1e6, 5),
    pvalue = runif(500),
    stringsAsFactors = FALSE
  )
  class(sim_data) <- c("glow_marginal_scan", "data.frame")

  .open_null_dev()
  on.exit(.close_null_dev(), add = TRUE)
  lambda <- plot(sim_data, type = "qq")

  expect_true(is.numeric(lambda))
  expect_true(lambda > 0)
})

test_that("plot.glow_marginal_scan type = 'both' returns lambda", {
  sim_data <- data.frame(
    chr    = rep(1:5, each = 100),
    pos    = rep(seq_len(100) * 1e6, 5),
    pvalue = runif(500),
    stringsAsFactors = FALSE
  )
  class(sim_data) <- c("glow_marginal_scan", "data.frame")

  .open_null_dev()
  on.exit(.close_null_dev(), add = TRUE)
  lambda <- plot(sim_data, type = "both")

  expect_true(is.numeric(lambda))
})

test_that("plot.glow_marginal_scan restores par after 'both'", {
  sim_data <- data.frame(
    chr    = rep(1:3, each = 50),
    pos    = rep(seq_len(50) * 1e6, 3),
    pvalue = runif(150),
    stringsAsFactors = FALSE
  )
  class(sim_data) <- c("glow_marginal_scan", "data.frame")

  .open_null_dev()
  on.exit(.close_null_dev(), add = TRUE)

  par_before <- graphics::par("mfrow")
  plot(sim_data, type = "both")
  par_after  <- graphics::par("mfrow")

  # par should be restored to what it was before the call
  expect_equal(par_before, par_after)
})

test_that("plot.glow_marginal_scan errors on invalid type", {
  sim_data <- data.frame(chr = 1, pos = 1, pvalue = 0.5)
  class(sim_data) <- c("glow_marginal_scan", "data.frame")
  expect_error(plot(sim_data, type = "invalid"), "should be one of")
})


# ==============================================================================
# Phase 4: Integration test with real GDS data
# ==============================================================================

test_that("integration: marginal_scan results can be plotted", {
  # Resolve test GDS path relative to testthat directory
  test_gds_path <- normalizePath(
    file.path(
      testthat::test_path(), "..", "..", "..", "..",
      "data", "local", "large-data", "test", "marginal_test_chr22.gds"
    ),
    mustWork = FALSE
  )
  skip_if(!file.exists(test_gds_path), "Test GDS file not available")

  # Run marginal scan (include 'sex' covariate so X is not NULL)
  pheno_covar <- extract_pheno_covar_gds(test_gds_path,
                                          pheno_name = "phenotype",
                                          covar_names = "sex",
                                          verbose = 0)
  nm <- fit_null_model(
    X         = pheno_covar$X,
    Y         = pheno_covar$Y,
    trait     = "binary",
    sample_id = pheno_covar$sample_id
  )
  results <- marginal_scan(test_gds_path, nm, mac_cutoff = 1, verbose = 0)

  # Class check
  expect_true(inherits(results, "glow_marginal_scan"))
  expect_true(is.data.frame(results))

  .open_null_dev()
  on.exit(.close_null_dev(), add = TRUE)

  # Manhattan
  expect_silent(plot_manhattan(results))

  # QQ
  lambda <- plot_qq(results)
  expect_true(is.numeric(lambda))
  expect_true(lambda > 0)

  # S3 dispatch "both"
  lambda2 <- plot(results, type = "both")
  expect_true(is.numeric(lambda2))
})
