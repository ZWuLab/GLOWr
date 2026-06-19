# Test file for core GLOW test calculation functions

# Load required packages
library(testthat)
library(GFisher)
library(MASS)  # For mvrnorm

# ============================================================================
# Helper Functions for Testing
# ============================================================================

# Generate test data with correlation
generate_test_data <- function(n = 10, rho = 0.3, seed = 123) {
  set.seed(seed)
  # Create correlation matrix with compound symmetry
  M <- matrix(rho, n, n) + diag(1 - rho, n, n)
  # Generate correlated Z-scores
  Z <- as.vector(mvrnorm(1, mu = rep(0, n), Sigma = M))
  list(Z = Z, M = M)
}

# ============================================================================
# Tests for burden_test()
# ============================================================================

test_that("burden_test uses normal distribution (NOT Davies method)", {
  # Generate test data
  data <- generate_test_data(n = 10, rho = 0.3)
  Z <- data$Z
  M <- data$M
  wts <- rep(1, 10)

  # Run burden test
  result <- burden_test(Z, M, wts, calc_p = TRUE)

  # Verify structure
  expect_type(result, "list")
  expect_named(result, c("S", "p"))

  # Verify S is a linear combination (NOT squared)
  S_expected <- sum(wts * Z) / sum(abs(wts))  # Weights are scaled
  expect_equal(result$S, S_expected, tolerance = 1e-10)

  # Verify p-value is from normal distribution
  # Manually compute p-value using normal distribution
  wts_scaled <- wts / sum(abs(wts))
  S_sd <- sqrt(t(wts_scaled) %*% M %*% wts_scaled)
  p_expected <- pnorm(abs(result$S), mean = 0, sd = S_sd, lower.tail = FALSE) * 2
  expect_equal(result$p, as.vector(p_expected), tolerance = 1e-10)

  # Verify p-value is in valid range
  expect_gte(result$p, 0)
  expect_lte(result$p, 1)
})

test_that("burden_test handles different weight schemes", {
  data <- generate_test_data(n = 5)
  Z <- data$Z
  M <- data$M

  # Equal weights
  result1 <- burden_test(Z, M, wts = rep(1, 5), calc_p = TRUE)

  # Non-equal weights
  result2 <- burden_test(Z, M, wts = c(1, 2, 3, 2, 1), calc_p = TRUE)

  # Verify both produce valid results
  expect_true(is.numeric(result1$S))
  expect_true(is.numeric(result2$S))
  expect_gte(result1$p, 0)
  expect_lte(result1$p, 1)
  expect_gte(result2$p, 0)
  expect_lte(result2$p, 1)

  # Different weights should give different statistics
  expect_false(isTRUE(all.equal(result1$S, result2$S)))
})

test_that("burden_test handles negative weights (forced to zero)", {
  data <- generate_test_data(n = 5)
  Z <- data$Z
  M <- data$M

  # Weights with negatives
  wts_neg <- c(1, -1, 2, -0.5, 3)

  # With is.posi.wts=TRUE (default), negatives become zero
  result <- burden_test(Z, M, wts = wts_neg, calc_p = TRUE, is.posi.wts = TRUE)

  # Verify it runs without error
  expect_true(is.numeric(result$S))
  expect_true(is.numeric(result$p))

  # The effective weights should be c(1, 0, 2, 0, 3)
  wts_effective <- pmax(wts_neg, 0)
  wts_effective <- wts_effective / sum(abs(wts_effective))
  S_expected <- sum(wts_effective * Z)
  expect_equal(result$S, S_expected, tolerance = 1e-10)
})

test_that("burden_test handles independence (identity correlation)", {
  set.seed(456)
  Z <- rnorm(10)
  M <- diag(10)
  wts <- rep(1, 10)

  result <- burden_test(Z, M, wts, calc_p = TRUE)

  # Under independence, variance is simply sum(wts^2) / sum(wts)^2
  wts_scaled <- wts / sum(abs(wts))
  S_var <- sum(wts_scaled^2)  # Since M = I
  S_sd <- sqrt(S_var)
  p_expected <- pnorm(abs(result$S), mean = 0, sd = S_sd, lower.tail = FALSE) * 2

  expect_equal(result$p, p_expected, tolerance = 1e-10)
})

# ============================================================================
# Tests for skat_test()
# ============================================================================

test_that("skat_test calls p.GFisher with df=1 (NOT Liu method)", {
  # Generate test data
  data <- generate_test_data(n = 8, rho = 0.4)
  Z <- data$Z
  M <- data$M
  wts <- rep(1, 8)

  # Run SKAT test
  result <- skat_test(Z, M_s = M, wts = wts, calc_p = TRUE)

  # Verify structure
  expect_type(result, "list")
  expect_named(result, c("S", "p"))

  # Verify S is a quadratic form
  wts_scaled <- wts / sum(abs(wts))
  S_expected <- sum(wts_scaled * Z^2)
  expect_equal(result$S, S_expected, tolerance = 1e-10)

  # Manually compute p-value using p.GFisher with df=1
  p_manual <- GFisher::p.GFisher(
    q = result$S,
    df = 1,
    w = wts_scaled,
    M = M,
    p.type = "two",
    method = "HYB"
  )
  expect_equal(result$p, p_manual, tolerance = 1e-10)

  # Verify p-value is in valid range
  expect_gte(result$p, 0)
  expect_lte(result$p, 1)
})

test_that("skat_test handles non-standardized scores", {
  # Generate covariance matrix (not correlation)
  set.seed(789)
  n <- 5
  SDs <- c(1, 2, 1.5, 3, 0.8)
  R <- matrix(0.3, n, n) + diag(0.7, n, n)  # Correlation
  M_cov <- diag(SDs) %*% R %*% diag(SDs)     # Covariance

  # Generate scores from this covariance
  scores <- as.vector(mvrnorm(1, mu = rep(0, n), Sigma = M_cov))
  wts <- rep(1, n)

  # SKAT should handle this correctly
  result <- skat_test(scores, M_s = M_cov, wts = wts, calc_p = TRUE)

  # Verify it runs and produces valid results
  expect_true(is.numeric(result$S))
  expect_true(is.numeric(result$p))
  expect_gte(result$p, 0)
  expect_lte(result$p, 1)

  # Verify internal standardization
  # After standardization: scores_std = scores / SDs
  # Adjusted weights: wts_adj = wts * SDs^2
  scores_std <- scores / SDs
  wts_adj <- wts * SDs^2
  wts_adj_scaled <- wts_adj / sum(abs(wts_adj))
  S_expected <- sum(wts_adj_scaled * scores_std^2)
  expect_equal(result$S, S_expected, tolerance = 1e-10)
})

test_that("skat_test handles different weight schemes", {
  data <- generate_test_data(n = 6, rho = 0.2)
  Z <- data$Z
  M <- data$M

  # Equal weights
  result1 <- skat_test(Z, M_s = M, wts = rep(1, 6), calc_p = TRUE)

  # Beta weights (common in rare variant analysis)
  beta_wts <- dbeta(seq(0.1, 0.6, length.out = 6), 1, 25)
  result2 <- skat_test(Z, M_s = M, wts = beta_wts, calc_p = TRUE)

  # Verify both produce valid results
  expect_true(is.numeric(result1$S))
  expect_true(is.numeric(result2$S))
  expect_gte(result1$p, 0)
  expect_lte(result1$p, 1)
  expect_gte(result2$p, 0)
  expect_lte(result2$p, 1)

  # Different weights should give different statistics
  expect_false(isTRUE(all.equal(result1$S, result2$S)))
})

# ============================================================================
# Tests for fisher_test_Z()
# ============================================================================

test_that("fisher_test_Z calls p.GFisher with df=2", {
  # Generate test data
  data <- generate_test_data(n = 10, rho = 0.25)
  Z <- data$Z
  M <- data$M
  wts <- rep(1, 10)

  # Run Fisher test
  result <- fisher_test_Z(Z, M, wts = wts, calc_p = TRUE, p.type = "two")

  # Verify structure
  expect_type(result, "list")
  expect_named(result, c("S", "p"))

  # Verify S uses g_GFisher transformation
  g <- function(x) g_GFisher(x, df = 2, p.type = "two")
  wts_scaled <- pmax(wts, 0) / sum(abs(pmax(wts, 0)))
  S_expected <- sum(wts_scaled * g(Z))
  expect_equal(result$S, S_expected, tolerance = 1e-10)

  # Manually compute p-value using p.GFisher with df=2
  p_manual <- GFisher::p.GFisher(
    q = result$S,
    df = 2,
    w = wts_scaled,
    M = M,
    p.type = "two",
    method = "HYB"
  )
  expect_equal(result$p, p_manual, tolerance = 1e-10)

  # Verify p-value is in valid range
  expect_gte(result$p, 0)
  expect_lte(result$p, 1)
})

test_that("fisher_test_Z handles one-sided vs two-sided", {
  data <- generate_test_data(n = 8, rho = 0.3)
  Z <- data$Z
  M <- data$M
  wts <- rep(1, 8)

  # Two-sided test (uses HYB method)
  result_two <- fisher_test_Z(Z, M, wts = wts, calc_p = TRUE, p.type = "two")

  # One-sided test (must use MR method since HYB only works for two-sided)
  result_one <- fisher_test_Z(Z, M, wts = wts, calc_p = TRUE, p.type = "one",
                              method = "MR", nsim = 1e4)

  # Both should be valid
  expect_true(is.numeric(result_two$S))
  expect_true(is.numeric(result_one$S))
  expect_gte(result_two$p, 0)
  expect_lte(result_two$p, 1)
  expect_gte(result_one$p, 0)
  expect_lte(result_one$p, 1)

  # Statistics should differ (different transformations)
  expect_false(isTRUE(all.equal(result_two$S, result_one$S)))
})

# ============================================================================
# Tests for calcu_SgZ_p() routing logic
# ============================================================================

test_that("calcu_SgZ_p correctly routes burden (df=Inf) to normal", {
  data <- generate_test_data(n = 5)
  Z <- data$Z
  M <- data$M
  wts <- rep(1, 5)

  # Test with df=Inf
  result_inf <- calcu_SgZ_p(
    g = function(x) x,
    Zscores = Z,
    wts = wts,
    calc_p = TRUE,
    M = M,
    df = Inf
  )

  # Should be identical to burden_test
  result_burden <- burden_test(Z, M, wts, calc_p = TRUE)

  expect_equal(result_inf$S, result_burden$S, tolerance = 1e-10)
  expect_equal(result_inf$p, result_burden$p, tolerance = 1e-10)
})

test_that("calcu_SgZ_p correctly routes SKAT (df=1) to p.GFisher", {
  data <- generate_test_data(n = 5)
  Z <- data$Z
  M <- data$M
  wts <- rep(1, 5)

  # Test with df=1
  result_df1 <- calcu_SgZ_p(
    g = function(x) x^2,
    Zscores = Z,
    wts = wts,
    calc_p = TRUE,
    M = M,
    df = 1
  )

  # Should be identical to skat_test
  result_skat <- skat_test(Z, M_s = M, wts, calc_p = TRUE)

  expect_equal(result_df1$S, result_skat$S, tolerance = 1e-10)
  expect_equal(result_df1$p, result_skat$p, tolerance = 1e-10)
})

test_that("calcu_SgZ_p correctly routes Fisher (df=2) to p.GFisher", {
  data <- generate_test_data(n = 5)
  Z <- data$Z
  M <- data$M
  wts <- rep(1, 5)

  # Test with df=2
  g <- function(x, df = 2, p.type = "two") g_GFisher(x, df, p.type)
  result_df2 <- calcu_SgZ_p(
    g = g,
    Zscores = Z,
    wts = wts,
    calc_p = TRUE,
    M = M,
    df = 2,
    p.type = "two"
  )

  # Should be identical to fisher_test_Z
  result_fisher <- fisher_test_Z(Z, M, wts, calc_p = TRUE, p.type = "two")

  expect_equal(result_df2$S, result_fisher$S, tolerance = 1e-10)
  expect_equal(result_df2$p, result_fisher$p, tolerance = 1e-10)
})

# ============================================================================
# Tests for multi_SgZ_test()
# ============================================================================

test_that("multi_SgZ_test computes multiple tests correctly", {
  data <- generate_test_data(n = 10, rho = 0.3)
  Z <- data$Z
  M <- data$M

  # Define three tests: Burden, SKAT, Fisher
  DF <- matrix(c(Inf, 1, 2), ncol = 1)
  W <- rbind(rep(1, 10), rep(1, 10), rep(1, 10))

  result <- multi_SgZ_test(Z, DF, W, M, p.type = "two", calcu_p = TRUE)

  # Verify structure
  expect_type(result, "list")
  expect_named(result, c("STAT", "PVAL"))
  expect_equal(nrow(result$STAT), 3)
  expect_equal(nrow(result$PVAL), 3)

  # Verify individual tests match
  burden_res <- burden_test(Z, M, wts = rep(1, 10), calc_p = TRUE)
  skat_res <- skat_test(Z, M_s = M, wts = rep(1, 10), calc_p = TRUE)
  fisher_res <- fisher_test_Z(Z, M, wts = rep(1, 10), calc_p = TRUE, p.type = "two")

  expect_equal(as.numeric(result$STAT[1, 1]), as.numeric(burden_res$S), tolerance = 1e-10)
  expect_equal(as.numeric(result$STAT[2, 1]), as.numeric(skat_res$S), tolerance = 1e-10)
  expect_equal(as.numeric(result$STAT[3, 1]), as.numeric(fisher_res$S), tolerance = 1e-10)

  expect_equal(as.numeric(result$PVAL[1, 1]), as.numeric(burden_res$p), tolerance = 1e-10)
  expect_equal(as.numeric(result$PVAL[2, 1]), as.numeric(skat_res$p), tolerance = 1e-10)
  expect_equal(as.numeric(result$PVAL[3, 1]), as.numeric(fisher_res$p), tolerance = 1e-10)
})

test_that("multi_SgZ_test handles multiple weight schemes", {
  data <- generate_test_data(n = 5, rho = 0.4)
  Z <- data$Z
  M <- data$M

  # Two SKAT tests with different weights
  DF <- matrix(c(1, 1), ncol = 1)
  W <- rbind(rep(1, 5), c(1, 2, 3, 2, 1))
  rownames(W) <- c("equal", "custom")

  result <- multi_SgZ_test(Z, DF, W, M, p.type = "two", calcu_p = TRUE)

  # Verify row names include weight names
  expect_match(rownames(result$STAT)[1], "equal")
  expect_match(rownames(result$STAT)[2], "custom")

  # Verify different weights give different results
  expect_false(isTRUE(all.equal(result$STAT[1, 1], result$STAT[2, 1])))
})

# ============================================================================
# Tests for omni_SgZ_test()
# ============================================================================

test_that("omni_SgZ_test combines tests using CCT", {
  data <- generate_test_data(n = 10, rho = 0.3)
  Z <- data$Z
  M <- data$M

  # Define three tests: Burden, SKAT, Fisher
  DF <- matrix(c(Inf, 1, 2), ncol = 1)
  W <- rbind(rep(1, 10), rep(1, 10), rep(1, 10))

  result <- omni_SgZ_test(Z, DF, W, M, p.type = "two", calcu_p = TRUE)

  # Verify structure
  expect_type(result, "list")
  expect_named(result, c("STAT", "PVAL", "cct", "pval_cct"))
  expect_equal(nrow(result$STAT), 4)  # 3 tests + 1 CCT
  expect_equal(nrow(result$PVAL), 4)

  # Verify CCT is computed correctly
  pvals <- result$PVAL[1:3, 1]
  cct_manual <- cct_test(pvals)

  expect_equal(result$cct, cct_manual$cct, tolerance = 1e-10)
  expect_equal(result$pval_cct, cct_manual$pval_cct, tolerance = 1e-10)

  # Verify CCT appended to matrices
  expect_equal(as.numeric(result$STAT[4, 1]), as.numeric(result$cct), tolerance = 1e-10)
  expect_equal(as.numeric(result$PVAL[4, 1]), as.numeric(result$pval_cct), tolerance = 1e-10)
})

# ============================================================================
# Edge Case Tests
# ============================================================================

test_that("All functions handle single variant case", {
  set.seed(999)
  Z <- rnorm(1)
  M <- matrix(1, 1, 1)
  wts <- 1

  # Burden test
  result_burden <- burden_test(Z, M, wts, calc_p = TRUE)
  expect_true(is.numeric(result_burden$S))
  expect_true(is.numeric(result_burden$p))

  # SKAT test
  result_skat <- skat_test(Z, M_s = M, wts, calc_p = TRUE)
  expect_true(is.numeric(result_skat$S))
  expect_true(is.numeric(result_skat$p))

  # Fisher test
  result_fisher <- fisher_test_Z(Z, M, wts, calc_p = TRUE, p.type = "two")
  expect_true(is.numeric(result_fisher$S))
  expect_true(is.numeric(result_fisher$p))
})

test_that("All functions handle perfect correlation", {
  # Perfect correlation (all variants identical)
  set.seed(888)
  n <- 5
  M <- matrix(1, n, n)  # Perfect correlation
  Z <- rep(rnorm(1), n)  # All identical Z-scores
  wts <- rep(1, n)

  # All tests should run without error
  result_burden <- burden_test(Z, M, wts, calc_p = TRUE)
  expect_true(is.numeric(result_burden$S))
  expect_true(is.numeric(result_burden$p))

  result_skat <- skat_test(Z, M_s = M, wts, calc_p = TRUE)
  expect_true(is.numeric(result_skat$S))
  expect_true(is.numeric(result_skat$p))

  result_fisher <- fisher_test_Z(Z, M, wts, calc_p = TRUE, p.type = "two")
  expect_true(is.numeric(result_fisher$S))
  expect_true(is.numeric(result_fisher$p))
})

test_that("All functions handle zero weights", {
  data <- generate_test_data(n = 5)
  Z <- data$Z
  M <- data$M
  wts <- c(0, 0, 0, 0, 0)  # All zero weights

  # With all zero weights, functions should return S=0 and p=1
  result_burden <- burden_test(Z, M, wts, calc_p = TRUE)
  expect_equal(result_burden$S, 0)
  expect_equal(result_burden$p, 1)

  result_skat <- skat_test(Z, M_s = M, wts, calc_p = TRUE)
  expect_equal(result_skat$S, 0)
  expect_equal(result_skat$p, 1)

  result_fisher <- fisher_test_Z(Z, M, wts, calc_p = TRUE, p.type = "two")
  expect_equal(result_fisher$S, 0)
  expect_equal(result_fisher$p, 1)
})

# ============================================================================
# Consistency Tests
# ============================================================================

test_that("calcu_p=FALSE returns only statistic", {
  data <- generate_test_data(n = 5)
  Z <- data$Z
  M <- data$M
  wts <- rep(1, 5)

  # Burden test without p-value
  result <- burden_test(Z, M, wts, calc_p = FALSE)
  expect_named(result, "S")
  expect_true(is.numeric(result$S))

  # SKAT test without p-value
  result <- skat_test(Z, M_s = M, wts, calc_p = FALSE)
  expect_named(result, "S")
  expect_true(is.numeric(result$S))

  # Fisher test without p-value
  result <- fisher_test_Z(Z, M, wts, calc_p = FALSE, p.type = "two")
  expect_named(result, "S")
  expect_true(is.numeric(result$S))
})

test_that("Results are reproducible with same seed", {
  # Test 1
  set.seed(777)
  data1 <- generate_test_data(n = 10, rho = 0.3, seed = 777)
  result1 <- burden_test(data1$Z, data1$M, rep(1, 10), calc_p = TRUE)

  # Test 2 with same seed
  set.seed(777)
  data2 <- generate_test_data(n = 10, rho = 0.3, seed = 777)
  result2 <- burden_test(data2$Z, data2$M, rep(1, 10), calc_p = TRUE)

  expect_equal(result1$S, result2$S)
  expect_equal(result1$p, result2$p)
})

# ============================================================================
# Print summary message
# ============================================================================

test_that("Test suite completion message", {
  message("\n================================")
  message("Core helpers test suite completed")
  message("All routing logic verified:")
  message("  - Burden -> Normal distribution")
  message("  - SKAT -> p.GFisher with df=1")
  message("  - Fisher -> p.GFisher with df=2")
  message("================================\n")
  expect_true(TRUE)
})
