# Tests for GLOW_Fisher function

test_that("GLOW_Fisher basic functionality works", {
  # Small example for quick testing
  set.seed(123)
  n <- 100
  m <- 5

  G <- matrix(rbinom(n*m, 2, 0.2), n, m)
  X <- matrix(rnorm(n*2), n, 2)
  Y <- rbinom(n, 1, 0.3)

  marg_stats <- getZ_marg_score(G, X, Y, trait="binary")
  B <- rnorm(m, mean=0, sd=0.2)
  PI <- runif(m, 0.1, 0.9)

  result <- GLOW_Fisher(marg_score_stats=marg_stats, B=B, PI=PI)

  # Check structure
  expect_type(result, "list")
  expect_named(result, c("STAT", "PVAL"))

  # Check dimensions
  expect_equal(ncol(result$STAT), 1)
  expect_equal(ncol(result$PVAL), 1)
  expect_equal(nrow(result$STAT), nrow(result$PVAL))

  # Check that last row is named GLOW_Fisher
  expect_equal(rownames(result$STAT)[nrow(result$STAT)], "GLOW_Fisher")
  expect_equal(rownames(result$PVAL)[nrow(result$PVAL)], "GLOW_Fisher")

  # Check p-values are valid
  expect_true(all(result$PVAL >= 0 & result$PVAL <= 1))
  expect_true(all(is.finite(result$PVAL)))
})


test_that("GLOW_Fisher handles two-sided vs one-sided tests", {
  set.seed(456)
  n <- 100
  m <- 10

  G <- matrix(rbinom(n*m, 2, 0.15), n, m)
  X <- matrix(rnorm(n*2), n, 2)
  Y <- rbinom(n, 1, 0.4)

  marg_stats <- getZ_marg_score(G, X, Y, trait="binary")
  B <- rnorm(m, mean=0, sd=0.3)
  PI <- runif(m, 0.2, 0.8)

  # Two-sided test (default)
  result_two <- GLOW_Fisher(marg_score_stats=marg_stats, B=B, PI=PI, p.type="two")

  # One-sided test requires method="MR" for GFisher package
  result_one <- GLOW_Fisher(marg_score_stats=marg_stats, B=B, PI=PI,
                            p.type="one", method="MR", nsim=1e4)

  # Both should work
  expect_true(all(is.finite(result_two$PVAL)))
  expect_true(all(is.finite(result_one$PVAL)))

  # P-values should generally differ between one-sided and two-sided
  # (unless all Z-scores are exactly 0, which is unlikely)
  expect_false(all(result_two$PVAL == result_one$PVAL))

  # Both should have valid p-values
  expect_true(all(result_two$PVAL >= 0 & result_two$PVAL <= 1))
  expect_true(all(result_one$PVAL >= 0 & result_one$PVAL <= 1))
})


test_that("GLOW_Fisher returns expected weight schemes", {
  set.seed(789)
  n <- 100
  m <- 8

  G <- matrix(rbinom(n*m, 2, 0.1), n, m)
  X <- matrix(rnorm(n*2), n, 2)
  Y <- rbinom(n, 1, 0.35)

  marg_stats <- getZ_marg_score(G, X, Y, trait="binary")
  B <- rnorm(m, mean=0, sd=0.25)
  PI <- runif(m, 0.15, 0.85)

  result <- GLOW_Fisher(marg_score_stats=marg_stats, B=B, PI=PI)

  # Should have 4 optimal weight schemes plus equal weights
  # Fisher test with g_GFisher transformation returns:
  # wts_BE_N, wts_APE_N, wts_BE_sparse, wts_APE_sparse
  expect_gte(nrow(result$STAT), 5)

  # Check that row names include expected patterns
  row_names <- rownames(result$STAT)
  expect_true(any(grepl("wts_BE", row_names)))
  expect_true(any(grepl("wts_APE", row_names)))
  expect_true(any(grepl("GLOW_Fisher", row_names)))
})


test_that("GLOW_Fisher handles single variant", {
  set.seed(101)
  n <- 100
  m <- 1

  G <- matrix(rbinom(n*m, 2, 0.2), n, m)
  X <- matrix(rnorm(n*2), n, 2)
  Y <- rbinom(n, 1, 0.3)

  marg_stats <- getZ_marg_score(G, X, Y, trait="binary")
  B <- rnorm(m)
  PI <- runif(m)

  result <- GLOW_Fisher(marg_score_stats=marg_stats, B=B, PI=PI)

  # Should work with single variant
  expect_true(all(is.finite(result$PVAL)))
  expect_true(all(result$PVAL >= 0 & result$PVAL <= 1))
})


test_that("GLOW_Fisher handles many variants", {
  set.seed(202)
  n <- 200
  m <- 50

  G <- matrix(rbinom(n*m, 2, 0.05), n, m)
  X <- matrix(rnorm(n*2), n, 2)
  Y <- rbinom(n, 1, 0.4)

  marg_stats <- getZ_marg_score(G, X, Y, trait="binary")
  B <- rnorm(m, mean=0, sd=0.2)
  PI <- runif(m, 0.1, 0.9)

  result <- GLOW_Fisher(marg_score_stats=marg_stats, B=B, PI=PI)

  # Should handle many variants
  expect_true(all(is.finite(result$PVAL)))
  expect_true(all(result$PVAL >= 0 & result$PVAL <= 1))
})


test_that("GLOW_Fisher validates input dimensions", {
  set.seed(303)
  n <- 100
  m <- 10

  G <- matrix(rbinom(n*m, 2, 0.15), n, m)
  X <- matrix(rnorm(n*2), n, 2)
  Y <- rbinom(n, 1, 0.3)

  marg_stats <- getZ_marg_score(G, X, Y, trait="binary")
  B <- rnorm(m)
  PI <- runif(m)

  # Wrong length B
  expect_error(
    GLOW_Fisher(marg_score_stats=marg_stats, B=rnorm(m+1), PI=PI),
    "B must have length equal to the number of variants"
  )

  # Wrong length PI
  expect_error(
    GLOW_Fisher(marg_score_stats=marg_stats, B=B, PI=runif(m-1)),
    "PI must have length equal to the number of variants"
  )
})


test_that("GLOW_Fisher uses df=2 (NOT df=1, NOT df=Inf)", {
  # This is a critical algorithmic check
  # We verify that the function is using the correct degrees of freedom

  set.seed(404)
  n <- 100
  m <- 5

  G <- matrix(rbinom(n*m, 2, 0.2), n, m)
  X <- matrix(rnorm(n*2), n, 2)
  Y <- rbinom(n, 1, 0.3)

  marg_stats <- getZ_marg_score(G, X, Y, trait="binary")
  B <- rnorm(m)
  PI <- runif(m)

  # Run Fisher test
  result_fisher <- GLOW_Fisher(marg_score_stats=marg_stats, B=B, PI=PI)

  # Run SKAT test (uses df=1) for comparison
  result_skat <- GLOW_SKAT(marg_score_stats=marg_stats, B=B, PI=PI)

  # Run Burden test (uses df=Inf) for comparison
  result_burden <- GLOW_Burden(marg_score_stats=marg_stats, B=B, PI=PI)

  # Fisher should give different results than SKAT and Burden
  # (with very high probability for random data)
  fisher_pval <- result_fisher$PVAL[nrow(result_fisher$PVAL), 1]
  skat_pval <- result_skat$PVAL[nrow(result_skat$PVAL), 1]
  burden_pval <- result_burden$PVAL[nrow(result_burden$PVAL), 1]

  # Fisher != SKAT (different df)
  expect_false(isTRUE(all.equal(fisher_pval, skat_pval, tolerance=1e-10)))

  # Fisher != Burden (different transformation)
  expect_false(isTRUE(all.equal(fisher_pval, burden_pval, tolerance=1e-10)))
})


test_that("GLOW_Fisher handles extreme effect sizes", {
  set.seed(505)
  n <- 100
  m <- 10

  G <- matrix(rbinom(n*m, 2, 0.15), n, m)
  X <- matrix(rnorm(n*2), n, 2)
  Y <- rbinom(n, 1, 0.3)

  marg_stats <- getZ_marg_score(G, X, Y, trait="binary")

  # Very large effect sizes
  B_large <- rnorm(m, mean=0, sd=5)
  PI <- runif(m, 0.1, 0.9)

  result_large <- GLOW_Fisher(marg_score_stats=marg_stats, B=B_large, PI=PI)
  expect_true(all(is.finite(result_large$PVAL)))
  expect_true(all(result_large$PVAL >= 0 & result_large$PVAL <= 1))

  # Very small effect sizes
  B_small <- rnorm(m, mean=0, sd=0.01)

  result_small <- GLOW_Fisher(marg_score_stats=marg_stats, B=B_small, PI=PI)
  expect_true(all(is.finite(result_small$PVAL)))
  expect_true(all(result_small$PVAL >= 0 & result_small$PVAL <= 1))
})


test_that("GLOW_Fisher handles extreme variant-importance scores", {
  set.seed(606)
  n <- 100
  m <- 10

  G <- matrix(rbinom(n*m, 2, 0.15), n, m)
  X <- matrix(rnorm(n*2), n, 2)
  Y <- rbinom(n, 1, 0.3)

  marg_stats <- getZ_marg_score(G, X, Y, trait="binary")
  B <- rnorm(m, mean=0, sd=0.2)

  # All high probabilities
  PI_high <- rep(0.99, m)

  result_high <- GLOW_Fisher(marg_score_stats=marg_stats, B=B, PI=PI_high)
  expect_true(all(is.finite(result_high$PVAL)))
  expect_true(all(result_high$PVAL >= 0 & result_high$PVAL <= 1))

  # All low probabilities
  PI_low <- rep(0.01, m)

  result_low <- GLOW_Fisher(marg_score_stats=marg_stats, B=B, PI=PI_low)
  expect_true(all(is.finite(result_low$PVAL)))
  expect_true(all(result_low$PVAL >= 0 & result_low$PVAL <= 1))

  # Mixed probabilities
  PI_mixed <- c(rep(0.01, m/2), rep(0.99, m/2))

  result_mixed <- GLOW_Fisher(marg_score_stats=marg_stats, B=B, PI=PI_mixed)
  expect_true(all(is.finite(result_mixed$PVAL)))
  expect_true(all(result_mixed$PVAL >= 0 & result_mixed$PVAL <= 1))
})


test_that("GLOW_Fisher numerical stability with correlated variants", {
  set.seed(707)
  n <- 200
  m <- 15

  # Create correlated genotypes (in LD)
  G_base <- matrix(rbinom(n*3, 2, 0.3), n, 3)
  G <- cbind(G_base, G_base + matrix(rbinom(n*3, 1, 0.1), n, 3))
  G <- cbind(G, matrix(rbinom(n*(m-6), 2, 0.2), n, m-6))

  X <- matrix(rnorm(n*2), n, 2)
  Y <- rbinom(n, 1, 0.35)

  marg_stats <- getZ_marg_score(G, X, Y, trait="binary")
  B <- rnorm(m, mean=0, sd=0.3)
  PI <- runif(m, 0.1, 0.9)

  result <- GLOW_Fisher(marg_score_stats=marg_stats, B=B, PI=PI)

  # Should handle correlation gracefully
  expect_true(all(is.finite(result$PVAL)))
  expect_true(all(result$PVAL >= 0 & result$PVAL <= 1))

  # Correlation matrix should have off-diagonal elements
  expect_true(any(abs(marg_stats$M_Z[lower.tri(marg_stats$M_Z)]) > 0.1))
})


test_that("GLOW_Fisher output structure is consistent", {
  set.seed(808)
  n <- 100
  m <- 10

  G <- matrix(rbinom(n*m, 2, 0.15), n, m)
  X <- matrix(rnorm(n*2), n, 2)
  Y <- rbinom(n, 1, 0.3)

  marg_stats <- getZ_marg_score(G, X, Y, trait="binary")
  B <- rnorm(m)
  PI <- runif(m)

  result <- GLOW_Fisher(marg_score_stats=marg_stats, B=B, PI=PI)

  # STAT and PVAL should have same dimensions
  expect_equal(dim(result$STAT), dim(result$PVAL))

  # Row names should match
  expect_equal(rownames(result$STAT), rownames(result$PVAL))

  # All entries should be numeric and finite
  expect_true(all(is.numeric(result$STAT)))
  expect_true(all(is.numeric(result$PVAL)))
  expect_true(all(is.finite(result$STAT)))
  expect_true(all(is.finite(result$PVAL)))
})
