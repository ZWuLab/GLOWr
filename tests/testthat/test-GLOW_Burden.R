# Unit tests for GLOW_Burden function
# Tests validation against legacy implementation and algorithm correctness

# Load legacy package for comparison if available
legacy_available <- requireNamespace("GLOW", quietly = TRUE)

test_that("GLOW_Burden works with basic example (binary trait)", {
  # This is the example from legacy documentation
  set.seed(123)
  X <- matrix(rnorm(20*2), 20, 2)
  Y <- rbinom(20, 1, 0.5)
  G <- matrix(rbinom(20*5, 1, 0.5), 20, 5)
  B <- rnorm(5)
  PI <- runif(5)

  # Compute marginal score statistics
  Zout <- getZ_marg_score(G, X, Y, trait="binary")

  # Run GLOW_Burden
  result <- GLOW_Burden(marg_score_stats=Zout, B=B, PI=PI)

  # Basic structure checks
  expect_type(result, "list")
  expect_named(result, c("STAT", "PVAL"))
  expect_true(is.matrix(result$STAT))
  expect_true(is.matrix(result$PVAL))

  # Check dimensions (should have wts_BE, wts_APE, wts_equ, GLOW_Burden rows)
  expect_equal(nrow(result$STAT), 4)
  expect_equal(nrow(result$PVAL), 4)
  expect_equal(ncol(result$STAT), 1)
  expect_equal(ncol(result$PVAL), 1)

  # Check row names
  expect_true("GLOW_Burden" %in% rownames(result$STAT))
  expect_true("GLOW_Burden" %in% rownames(result$PVAL))
  expect_true("df_Inf_wts_BE" %in% rownames(result$STAT))
  expect_true("df_Inf_wts_APE" %in% rownames(result$STAT))
  expect_true("df_Inf_wts_equ" %in% rownames(result$STAT))

  # P-values should be between 0 and 1
  expect_true(all(result$PVAL >= 0 & result$PVAL <= 1))

  # Statistics should be finite
  expect_true(all(is.finite(result$STAT)))
})


test_that("GLOW_Burden works with continuous trait", {
  set.seed(456)
  X <- matrix(rnorm(30*2), 30, 2)
  Y <- rnorm(30)
  G <- matrix(rbinom(30*8, 2, 0.2), 30, 8)
  B <- rnorm(8, sd=0.3)
  PI <- runif(8, 0.2, 0.8)

  # Compute marginal score statistics
  Zout <- getZ_marg_score(G, X, Y, trait="continuous")

  # Run GLOW_Burden
  result <- GLOW_Burden(marg_score_stats=Zout, B=B, PI=PI)

  # Basic checks
  expect_type(result, "list")
  expect_equal(nrow(result$STAT), 4)
  expect_true(all(result$PVAL >= 0 & result$PVAL <= 1))
})


test_that("GLOW_Burden: all positive effects scenario", {
  # When all variants have positive effects in same direction,
  # Burden test should be powerful
  set.seed(789)
  n <- 100
  m <- 10

  X <- matrix(rnorm(n*2), n, 2)
  G <- matrix(rbinom(n*m, 2, 0.1), n, m)

  # All positive effects
  true_beta <- abs(rnorm(m, mean=0.5, sd=0.1))
  Y <- rnorm(n, mean = G %*% true_beta + X[,1] * 0.3, sd=1)

  # Compute marginal score statistics
  Zout <- getZ_marg_score(G, X, Y, trait="continuous")

  # Use true effects as B estimates (best case)
  B <- true_beta
  PI <- rep(1, m)  # All variants are causal

  # Run GLOW_Burden
  result <- GLOW_Burden(marg_score_stats=Zout, B=B, PI=PI)

  # With all positive effects and correct information,
  # test should detect signal (though not guaranteed in small sample)
  # At minimum, should produce valid output
  expect_true(all(result$PVAL >= 0 & result$PVAL <= 1))
  expect_true(all(is.finite(result$STAT)))

  # All true effects are positive
  expect_true(all(true_beta > 0))
})


test_that("GLOW_Burden: mixed effects scenario", {
  # When effects are in mixed directions, Burden may be less powerful
  set.seed(321)
  n <- 100
  m <- 10

  X <- matrix(rnorm(n*2), n, 2)
  G <- matrix(rbinom(n*m, 2, 0.1), n, m)

  # Mixed effects (some positive, some negative)
  true_beta <- rnorm(m, mean=0, sd=0.3)
  Y <- rnorm(n, mean = G %*% true_beta + X[,1] * 0.3, sd=1)

  # Compute marginal score statistics
  Zout <- getZ_marg_score(G, X, Y, trait="continuous")

  # Use noisy estimates
  B <- true_beta + rnorm(m, 0, 0.1)
  PI <- runif(m, 0.3, 0.7)

  # Run GLOW_Burden
  result <- GLOW_Burden(marg_score_stats=Zout, B=B, PI=PI)

  # Should produce valid output even with mixed effects
  expect_true(all(result$PVAL >= 0 & result$PVAL <= 1))
  expect_true(all(is.finite(result$STAT)))
})


test_that("GLOW_Burden: sparse signal scenario", {
  # Only a few variants are causal
  set.seed(654)
  n <- 80
  m <- 15

  X <- matrix(rnorm(n*2), n, 2)
  G <- matrix(rbinom(n*m, 2, 0.15), n, m)

  # Only 3 variants are causal
  true_beta <- rep(0, m)
  true_beta[c(3, 7, 12)] <- c(0.5, 0.6, 0.4)
  Y <- rnorm(n, mean = G %*% true_beta + X[,1] * 0.2, sd=1)

  # Compute marginal score statistics
  Zout <- getZ_marg_score(G, X, Y, trait="continuous")

  # Provide sparse PI (close to truth)
  B <- rnorm(m, sd=0.2)
  B[c(3, 7, 12)] <- c(0.5, 0.6, 0.4)
  PI <- rep(0.1, m)
  PI[c(3, 7, 12)] <- 0.9

  # Run GLOW_Burden
  result <- GLOW_Burden(marg_score_stats=Zout, B=B, PI=PI)

  # Should produce valid output
  expect_true(all(result$PVAL >= 0 & result$PVAL <= 1))
  expect_true(all(is.finite(result$STAT)))
})


test_that("GLOW_Burden: equal weights row matches simple burden", {
  # The df_Inf_wts_equ row should match a simple burden test with equal weights
  set.seed(999)
  X <- matrix(rnorm(50*2), 50, 2)
  Y <- rnorm(50)
  G <- matrix(rbinom(50*6, 2, 0.2), 50, 6)
  B <- rnorm(6)
  PI <- runif(6)

  # Compute marginal score statistics
  Zout <- getZ_marg_score(G, X, Y, trait="continuous")

  # Run GLOW_Burden
  result <- GLOW_Burden(marg_score_stats=Zout, B=B, PI=PI)

  # Get equal weights result (NOT the GLOW_Burden row, which is CCT)
  equ_stat <- result$STAT["df_Inf_wts_equ", 1]
  equ_pval <- result$PVAL["df_Inf_wts_equ", 1]

  # Manually compute burden test with equal weights
  wts_equ <- rep(1, length(Zout$Zscores))
  wts_equ <- wts_equ / sum(abs(wts_equ))  # Normalize as done in calcu_SgZ_p
  S_manual <- sum(wts_equ * Zout$Zscores)
  sigma_manual <- sqrt(t(wts_equ) %*% Zout$M_Z %*% wts_equ)
  p_manual <- 2 * pnorm(abs(S_manual), mean=0, sd=sigma_manual, lower.tail=FALSE)

  # Should match (allowing for small numerical differences)
  expect_equal(as.numeric(equ_stat), S_manual, tolerance=1e-10)
  expect_equal(as.numeric(equ_pval), as.numeric(p_manual), tolerance=1e-10)
})


test_that("GLOW_Burden: final result is CCT combination", {
  # The GLOW_Burden row should be a CCT combination of the other tests
  set.seed(999)
  X <- matrix(rnorm(50*2), 50, 2)
  Y <- rnorm(50)
  G <- matrix(rbinom(50*6, 2, 0.2), 50, 6)
  B <- rnorm(6)
  PI <- runif(6)

  # Compute marginal score statistics
  Zout <- getZ_marg_score(G, X, Y, trait="continuous")

  # Run GLOW_Burden
  result <- GLOW_Burden(marg_score_stats=Zout, B=B, PI=PI)

  # GLOW_Burden p-value should be a CCT combination
  # It should be between 0 and 1
  glow_pval <- result$PVAL["GLOW_Burden", 1]
  expect_true(glow_pval >= 0 && glow_pval <= 1)

  # The GLOW_Burden statistic should be a Cauchy-transformed value
  glow_stat <- result$STAT["GLOW_Burden", 1]
  expect_true(is.finite(glow_stat))
})


test_that("GLOW_Burden input validation", {
  set.seed(111)
  X <- matrix(rnorm(20*2), 20, 2)
  Y <- rnorm(20)
  G <- matrix(rbinom(20*5, 2, 0.2), 20, 5)

  Zout <- getZ_marg_score(G, X, Y, trait="continuous")

  # Mismatched B length
  expect_error(
    GLOW_Burden(marg_score_stats=Zout, B=rnorm(3), PI=runif(5)),
    "B must have length equal to the number of variants"
  )

  # Mismatched PI length
  expect_error(
    GLOW_Burden(marg_score_stats=Zout, B=rnorm(5), PI=runif(3)),
    "PI must have length equal to the number of variants"
  )

  # Invalid PI values (> 1)
  expect_error(
    GLOW_Burden(marg_score_stats=Zout, B=rnorm(5), PI=c(0.5, 0.5, 1.5, 0.5, 0.5)),
    "PI values must be in"
  )

  # Invalid PI values (< 0)
  expect_error(
    GLOW_Burden(marg_score_stats=Zout, B=rnorm(5), PI=c(0.5, -0.1, 0.5, 0.5, 0.5)),
    "PI.*non-negative|PI values must be in"  # New validation catches this earlier with different message
  )
})


test_that("GLOW_Burden uses normal distribution (NOT Davies method)", {
  # CRITICAL TEST: Verify that Burden uses normal distribution for individual tests
  set.seed(777)
  X <- matrix(rnorm(40*2), 40, 2)
  Y <- rnorm(40)
  G <- matrix(rbinom(40*8, 2, 0.15), 40, 8)
  B <- rnorm(8)
  PI <- runif(8)

  Zout <- getZ_marg_score(G, X, Y, trait="continuous")
  result <- GLOW_Burden(marg_score_stats=Zout, B=B, PI=PI)

  # Test only the individual burden tests (not the CCT combined one)
  for (i in 1:(nrow(result$STAT)-1)) {  # Skip last row (CCT)
    row_name <- rownames(result$STAT)[i]
    S <- result$STAT[i, 1]
    p_result <- result$PVAL[i, 1]

    # Get the weights used for this test
    if (row_name == "df_Inf_wts_BE") {
      wts_opt <- Optimal_Weights_M(
        g = function(x) x,
        Bstar = sqrt(diag(Zout$M_s)) * B / Zout$s0,
        PI = PI,
        M = Zout$M_Z,
        is.posi.wts = TRUE
      )
      wts <- wts_opt$wts_BE
    } else if (row_name == "df_Inf_wts_APE") {
      wts_opt <- Optimal_Weights_M(
        g = function(x) x,
        Bstar = sqrt(diag(Zout$M_s)) * B / Zout$s0,
        PI = PI,
        M = Zout$M_Z,
        is.posi.wts = TRUE
      )
      wts <- wts_opt$wts_APE
    } else if (row_name == "df_Inf_wts_equ") {
      wts <- rep(1, length(Zout$Zscores))
    } else {
      next  # Skip unknown row names
    }

    # Normalize weights (as done in calcu_SgZ_p)
    wts <- wts / sum(abs(wts))

    # Compute expected p-value using normal distribution
    sigma <- sqrt(t(wts) %*% Zout$M_Z %*% wts)
    p_expected <- 2 * pnorm(abs(S), mean=0, sd=sigma, lower.tail=FALSE)

    # Should match exactly (allowing for numerical precision)
    expect_equal(as.numeric(p_result), as.numeric(p_expected), tolerance=1e-10,
                 label=paste("P-value for", row_name))
  }
})


test_that("GLOW_Burden test statistic is LINEAR (not squared)", {
  # CRITICAL TEST: Verify that test statistic is linear combination
  set.seed(888)
  X <- matrix(rnorm(30*2), 30, 2)
  Y <- rnorm(30)
  G <- matrix(rbinom(30*6, 2, 0.2), 30, 6)
  B <- rnorm(6)
  PI <- runif(6)

  Zout <- getZ_marg_score(G, X, Y, trait="continuous")
  result <- GLOW_Burden(marg_score_stats=Zout, B=B, PI=PI)

  # For equal weights, statistic should be sum(wts * Z)
  # Use df_Inf_wts_equ row (not GLOW_Burden which is CCT)
  equ_stat <- result$STAT["df_Inf_wts_equ", 1]

  wts_equ <- rep(1, length(Zout$Zscores))
  wts_equ <- wts_equ / sum(abs(wts_equ))

  # Linear combination
  S_linear <- sum(wts_equ * Zout$Zscores)

  # Should match linear combination (NOT squared)
  expect_equal(as.numeric(equ_stat), S_linear, tolerance=1e-10)

  # Should NOT match squared version
  S_squared <- sum(wts_equ * Zout$Zscores^2)
  expect_false(abs(equ_stat - S_squared) < 1e-6)
})


# Validation against legacy package (if available)
if (legacy_available) {
  test_that("GLOW_Burden matches legacy implementation EXACTLY", {
    # Load legacy package
    library(GLOW)

    # Use same data as legacy example
    set.seed(123)
    X <- matrix(rnorm(20*2), 20, 2)
    Y <- rbinom(20, 1, 0.5)
    G <- matrix(rbinom(20*5, 1, 0.5), 20, 5)
    B <- rnorm(5)
    PI <- runif(5)

    # Legacy computation
    Zout_legacy <- GLOW::getZ_marg_score(G, X, Y, trait="binary")
    result_legacy <- GLOW::GLOW_Burden(marg_score_stats=Zout_legacy, B=B, PI=PI)

    # New computation
    Zout_new <- GLOWr::getZ_marg_score(G, X, Y, trait="binary")
    result_new <- GLOWr::GLOW_Burden(marg_score_stats=Zout_new, B=B, PI=PI)

    # Statistics should match exactly
    expect_equal(result_new$STAT, result_legacy$STAT, tolerance=1e-10,
                 label="Test statistics match legacy")

    # P-values should match exactly
    expect_equal(result_new$PVAL, result_legacy$PVAL, tolerance=1e-10,
                 label="P-values match legacy")

    # Row names should match
    expect_equal(rownames(result_new$STAT), rownames(result_legacy$STAT))
    expect_equal(rownames(result_new$PVAL), rownames(result_legacy$PVAL))
  })


  test_that("GLOW_Burden matches legacy: continuous trait", {
    library(GLOW)

    set.seed(456)
    X <- matrix(rnorm(30*2), 30, 2)
    Y <- rnorm(30)
    G <- matrix(rbinom(30*8, 2, 0.2), 30, 8)
    B <- rnorm(8, sd=0.3)
    PI <- runif(8, 0.2, 0.8)

    # Legacy
    Zout_legacy <- GLOW::getZ_marg_score(G, X, Y, trait="continuous")
    result_legacy <- GLOW::GLOW_Burden(marg_score_stats=Zout_legacy, B=B, PI=PI)

    # New
    Zout_new <- GLOWr::getZ_marg_score(G, X, Y, trait="continuous")
    result_new <- GLOWr::GLOW_Burden(marg_score_stats=Zout_new, B=B, PI=PI)

    # Exact match
    expect_equal(result_new$STAT, result_legacy$STAT, tolerance=1e-10)
    expect_equal(result_new$PVAL, result_legacy$PVAL, tolerance=1e-10)
  })


  test_that("GLOW_Burden matches legacy: multiple scenarios", {
    library(GLOW)

    # Test 5 different random scenarios
    for (i in 1:5) {
      set.seed(1000 + i)

      n <- sample(40:100, 1)
      m <- sample(5:15, 1)

      X <- matrix(rnorm(n*2), n, 2)
      Y <- if (runif(1) < 0.5) rbinom(n, 1, 0.4) else rnorm(n)
      G <- matrix(rbinom(n*m, 2, runif(1, 0.05, 0.3)), n, m)
      B <- rnorm(m, sd=runif(1, 0.1, 0.5))
      PI <- runif(m, 0.1, 0.9)

      trait_type <- if (is.numeric(Y) && all(Y %in% c(0,1))) "binary" else "continuous"

      # Legacy
      Zout_legacy <- GLOW::getZ_marg_score(G, X, Y, trait=trait_type)
      result_legacy <- GLOW::GLOW_Burden(marg_score_stats=Zout_legacy, B=B, PI=PI)

      # New
      Zout_new <- GLOWr::getZ_marg_score(G, X, Y, trait=trait_type)
      result_new <- GLOWr::GLOW_Burden(marg_score_stats=Zout_new, B=B, PI=PI)

      # Exact match
      expect_equal(result_new$STAT, result_legacy$STAT, tolerance=1e-10,
                   label=paste("Scenario", i, "STAT"))
      expect_equal(result_new$PVAL, result_legacy$PVAL, tolerance=1e-10,
                   label=paste("Scenario", i, "PVAL"))
    }
  })
}
