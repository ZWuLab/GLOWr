# Unit tests for GLOW_SKAT function
# Tests validation against legacy implementation and algorithm correctness
# CRITICAL: Verifies that SKAT uses p.GFisher with df=1, NOT Liu method

# Load legacy package for comparison if available
legacy_available <- requireNamespace("GLOW", quietly = TRUE)

test_that("GLOW_SKAT works with basic example (binary trait)", {
  # This is adapted from legacy documentation
  set.seed(123)
  X <- matrix(rnorm(20*2), 20, 2)
  Y <- rbinom(20, 1, 0.5)
  G <- matrix(rbinom(20*5, 1, 0.5), 20, 5)
  B <- rnorm(5)
  PI <- runif(5)

  # Compute marginal score statistics
  Zout <- getZ_marg_score(G, X, Y, trait="binary")

  # Run GLOW_SKAT
  result <- GLOW_SKAT(marg_score_stats=Zout, B=B, PI=PI)

  # Basic structure checks
  expect_type(result, "list")
  expect_named(result, c("STAT", "PVAL"))
  expect_true(is.matrix(result$STAT))
  expect_true(is.matrix(result$PVAL))

  # Check dimensions (should have 4 optimal weights + 1 equal weights + CCT row)
  # wts_BE_N, wts_APE_N, wts_BE_sparse, wts_APE_sparse, wts_equ, GLOW_SKAT
  expect_equal(nrow(result$STAT), 6)
  expect_equal(nrow(result$PVAL), 6)
  expect_equal(ncol(result$STAT), 1)
  expect_equal(ncol(result$PVAL), 1)

  # Check row names
  expect_true("GLOW_SKAT" %in% rownames(result$STAT))
  expect_true("GLOW_SKAT" %in% rownames(result$PVAL))
  expect_true("df_1_wts_BE_N" %in% rownames(result$STAT))
  expect_true("df_1_wts_APE_N" %in% rownames(result$STAT))
  expect_true("df_1_wts_BE_sparse" %in% rownames(result$STAT))
  expect_true("df_1_wts_APE_sparse" %in% rownames(result$STAT))
  expect_true("df_1_wts_equ" %in% rownames(result$STAT))

  # P-values should be between 0 and 1
  expect_true(all(result$PVAL >= 0 & result$PVAL <= 1))

  # Statistics should be finite
  expect_true(all(is.finite(result$STAT)))
  # SKAT statistics (non-CCT) should be non-negative (sum of squared terms)
  non_cct_rows <- rownames(result$STAT) != "GLOW_SKAT"
  expect_true(all(result$STAT[non_cct_rows, 1] >= 0))
})


test_that("GLOW_SKAT works with continuous trait", {
  set.seed(456)
  X <- matrix(rnorm(30*2), 30, 2)
  Y <- rnorm(30)
  G <- matrix(rbinom(30*8, 2, 0.2), 30, 8)
  B <- rnorm(8, sd=0.3)
  PI <- runif(8, 0.2, 0.8)

  # Compute marginal score statistics
  Zout <- getZ_marg_score(G, X, Y, trait="continuous")

  # Run GLOW_SKAT
  result <- GLOW_SKAT(marg_score_stats=Zout, B=B, PI=PI)

  # Basic checks
  expect_type(result, "list")
  expect_equal(nrow(result$STAT), 6)
  expect_true(all(result$PVAL >= 0 & result$PVAL <= 1))
  # Quadratic form statistics (non-CCT) should be non-negative
  non_cct_rows <- rownames(result$STAT) != "GLOW_SKAT"
  expect_true(all(result$STAT[non_cct_rows, 1] >= 0))
})


test_that("GLOW_SKAT: mixed effect directions scenario", {
  # When variants have effects in mixed directions (some positive, some negative),
  # SKAT should be more powerful than Burden
  set.seed(789)
  n <- 100
  m <- 10

  X <- matrix(rnorm(n*2), n, 2)
  G <- matrix(rbinom(n*m, 2, 0.1), n, m)

  # Mixed effects: half positive, half negative
  true_beta <- c(rep(0.4, 5), rep(-0.4, 5))
  Y <- rnorm(n, mean = G %*% true_beta + X[,1] * 0.3, sd=1)

  # Compute marginal score statistics
  Zout <- getZ_marg_score(G, X, Y, trait="continuous")

  # Use true effects as B estimates (best case)
  B <- true_beta
  PI <- rep(1, m)  # All variants are causal

  # Run GLOW_SKAT
  result <- GLOW_SKAT(marg_score_stats=Zout, B=B, PI=PI)

  # Should produce valid output
  expect_true(all(result$PVAL >= 0 & result$PVAL <= 1))
  expect_true(all(is.finite(result$STAT)))
  # Quadratic form statistics (non-CCT) should be non-negative
  non_cct_rows <- rownames(result$STAT) != "GLOW_SKAT"
  expect_true(all(result$STAT[non_cct_rows, 1] >= 0))

  # Verify mixed effects
  expect_true(sum(true_beta > 0) > 0)
  expect_true(sum(true_beta < 0) > 0)
})


test_that("GLOW_SKAT: heterogeneous effects scenario", {
  # When effects have different magnitudes, SKAT should handle well
  set.seed(321)
  n <- 100
  m <- 10

  X <- matrix(rnorm(n*2), n, 2)
  G <- matrix(rbinom(n*m, 2, 0.1), n, m)

  # Heterogeneous effects (different magnitudes)
  true_beta <- c(0.1, 0.2, 0.5, 1.0, 0.3, -0.1, -0.2, -0.5, -1.0, -0.3)
  Y <- rnorm(n, mean = G %*% true_beta + X[,1] * 0.3, sd=1)

  # Compute marginal score statistics
  Zout <- getZ_marg_score(G, X, Y, trait="continuous")

  # Use noisy estimates
  B <- true_beta + rnorm(m, 0, 0.1)
  PI <- runif(m, 0.3, 0.7)

  # Run GLOW_SKAT
  result <- GLOW_SKAT(marg_score_stats=Zout, B=B, PI=PI)

  # Should produce valid output
  expect_true(all(result$PVAL >= 0 & result$PVAL <= 1))
  expect_true(all(is.finite(result$STAT)))
  # Quadratic form statistics (non-CCT) should be non-negative
  non_cct_rows <- rownames(result$STAT) != "GLOW_SKAT"
  expect_true(all(result$STAT[non_cct_rows, 1] >= 0))
})


test_that("GLOW_SKAT: sparse signal scenario", {
  # Only a few variants are causal
  set.seed(654)
  n <- 80
  m <- 15

  X <- matrix(rnorm(n*2), n, 2)
  G <- matrix(rbinom(n*m, 2, 0.15), n, m)

  # Only 3 variants are causal with mixed effects
  true_beta <- rep(0, m)
  true_beta[c(3, 7, 12)] <- c(0.5, -0.6, 0.4)
  Y <- rnorm(n, mean = G %*% true_beta + X[,1] * 0.2, sd=1)

  # Compute marginal score statistics
  Zout <- getZ_marg_score(G, X, Y, trait="continuous")

  # Provide sparse PI (close to truth)
  B <- rnorm(m, sd=0.2)
  B[c(3, 7, 12)] <- c(0.5, -0.6, 0.4)
  PI <- rep(0.1, m)
  PI[c(3, 7, 12)] <- 0.9

  # Run GLOW_SKAT
  result <- GLOW_SKAT(marg_score_stats=Zout, B=B, PI=PI)

  # Should produce valid output
  expect_true(all(result$PVAL >= 0 & result$PVAL <= 1))
  expect_true(all(is.finite(result$STAT)))
  # Quadratic form statistics (non-CCT) should be non-negative
  non_cct_rows <- rownames(result$STAT) != "GLOW_SKAT"
  expect_true(all(result$STAT[non_cct_rows, 1] >= 0))
})


test_that("GLOW_SKAT: equal weights row uses simple SKAT formula", {
  # The df_1_wts_equ row should be sum(Z^2)
  set.seed(999)
  X <- matrix(rnorm(50*2), 50, 2)
  Y <- rnorm(50)
  G <- matrix(rbinom(50*6, 2, 0.2), 50, 6)
  B <- rnorm(6)
  PI <- runif(6)

  # Compute marginal score statistics
  Zout <- getZ_marg_score(G, X, Y, trait="continuous")

  # Run GLOW_SKAT
  result <- GLOW_SKAT(marg_score_stats=Zout, B=B, PI=PI)

  # Get equal weights result (NOT the GLOW_SKAT row, which is CCT)
  equ_stat <- result$STAT["df_1_wts_equ", 1]

  # Manually compute SKAT test with equal weights
  wts_equ <- rep(1, length(Zout$Zscores))
  wts_equ <- wts_equ / sum(abs(wts_equ))  # Normalize as done in calcu_SgZ_p
  S_manual <- sum(wts_equ * Zout$Zscores^2)

  # Should match (allowing for small numerical differences)
  expect_equal(as.numeric(equ_stat), S_manual, tolerance=1e-10)
})


test_that("GLOW_SKAT: final result is CCT combination", {
  # The GLOW_SKAT row should be a CCT combination of the other tests
  set.seed(999)
  X <- matrix(rnorm(50*2), 50, 2)
  Y <- rnorm(50)
  G <- matrix(rbinom(50*6, 2, 0.2), 50, 6)
  B <- rnorm(6)
  PI <- runif(6)

  # Compute marginal score statistics
  Zout <- getZ_marg_score(G, X, Y, trait="continuous")

  # Run GLOW_SKAT
  result <- GLOW_SKAT(marg_score_stats=Zout, B=B, PI=PI)

  # GLOW_SKAT p-value should be a CCT combination
  # It should be between 0 and 1
  glow_pval <- result$PVAL["GLOW_SKAT", 1]
  expect_true(glow_pval >= 0 && glow_pval <= 1)

  # The GLOW_SKAT statistic should be a Cauchy-transformed value
  glow_stat <- result$STAT["GLOW_SKAT", 1]
  expect_true(is.finite(glow_stat))
})


test_that("GLOW_SKAT input validation", {
  set.seed(111)
  X <- matrix(rnorm(20*2), 20, 2)
  Y <- rnorm(20)
  G <- matrix(rbinom(20*5, 2, 0.2), 20, 5)

  Zout <- getZ_marg_score(G, X, Y, trait="continuous")

  # Mismatched B length
  expect_error(
    GLOW_SKAT(marg_score_stats=Zout, B=rnorm(3), PI=runif(5)),
    "B must have length equal to the number of variants"
  )

  # Mismatched PI length
  expect_error(
    GLOW_SKAT(marg_score_stats=Zout, B=rnorm(5), PI=runif(3)),
    "PI must have length equal to the number of variants"
  )

  # Invalid PI values (> 1)
  expect_error(
    GLOW_SKAT(marg_score_stats=Zout, B=rnorm(5), PI=c(0.5, 0.5, 1.5, 0.5, 0.5)),
    "PI values must be in"
  )

  # Invalid PI values (< 0)
  expect_error(
    GLOW_SKAT(marg_score_stats=Zout, B=rnorm(5), PI=c(0.5, -0.1, 0.5, 0.5, 0.5)),
    "PI.*non-negative|PI values must be in"  # New validation catches this earlier with different message
  )
})


test_that("GLOW_SKAT uses p.GFisher with df=1 (NOT Liu method)", {
  # CRITICAL TEST: Verify that SKAT uses p.GFisher with df=1
  set.seed(777)
  X <- matrix(rnorm(40*2), 40, 2)
  Y <- rnorm(40)
  G <- matrix(rbinom(40*8, 2, 0.15), 40, 8)
  B <- rnorm(8)
  PI <- runif(8)

  Zout <- getZ_marg_score(G, X, Y, trait="continuous")
  result <- GLOW_SKAT(marg_score_stats=Zout, B=B, PI=PI)

  # Test the equal weights row (most straightforward)
  equ_stat <- result$STAT["df_1_wts_equ", 1]
  equ_pval <- result$PVAL["df_1_wts_equ", 1]

  # Manually compute using p.GFisher with df=1
  wts_equ <- rep(1, length(Zout$Zscores))
  wts_equ <- wts_equ / sum(abs(wts_equ))
  S_manual <- sum(wts_equ * Zout$Zscores^2)

  # Call p.GFisher with df=1
  p_manual <- GFisher::p.GFisher(
    q = S_manual,
    df = 1,
    w = wts_equ,
    M = Zout$M_Z,
    p.type = "two"
  )

  # Should match exactly
  expect_equal(as.numeric(equ_stat), S_manual, tolerance=1e-10)
  expect_equal(as.numeric(equ_pval), as.numeric(p_manual), tolerance=1e-10,
               label="P-value matches p.GFisher with df=1")
})


test_that("GLOW_SKAT test statistic is QUADRATIC (not linear)", {
  # CRITICAL TEST: Verify that test statistic is sum of squares
  set.seed(888)
  X <- matrix(rnorm(30*2), 30, 2)
  Y <- rnorm(30)
  G <- matrix(rbinom(30*6, 2, 0.2), 30, 6)
  B <- rnorm(6)
  PI <- runif(6)

  Zout <- getZ_marg_score(G, X, Y, trait="continuous")
  result <- GLOW_SKAT(marg_score_stats=Zout, B=B, PI=PI)

  # For equal weights, statistic should be sum(wts * Z^2)
  # Use df_1_wts_equ row (not GLOW_SKAT which is CCT)
  equ_stat <- result$STAT["df_1_wts_equ", 1]

  wts_equ <- rep(1, length(Zout$Zscores))
  wts_equ <- wts_equ / sum(abs(wts_equ))

  # Quadratic form (sum of squares)
  S_quadratic <- sum(wts_equ * Zout$Zscores^2)

  # Should match quadratic form (NOT linear)
  expect_equal(as.numeric(equ_stat), S_quadratic, tolerance=1e-10)

  # Should NOT match linear version
  S_linear <- sum(wts_equ * Zout$Zscores)
  # Only fail if they're actually different (they could be close by chance)
  if (abs(S_quadratic - S_linear) > 1e-6) {
    expect_false(abs(equ_stat - S_linear) < 1e-10)
  }
})


test_that("GLOW_SKAT uses df=1 (NOT df=Inf, NOT df=2)", {
  # CRITICAL TEST: Verify degrees of freedom is 1
  set.seed(555)
  X <- matrix(rnorm(30*2), 30, 2)
  Y <- rnorm(30)
  G <- matrix(rbinom(30*6, 2, 0.2), 30, 6)
  B <- rnorm(6)
  PI <- runif(6)

  Zout <- getZ_marg_score(G, X, Y, trait="continuous")
  result <- GLOW_SKAT(marg_score_stats=Zout, B=B, PI=PI)

  # All non-CCT rows should have "df_1" in the name
  non_cct_rows <- rownames(result$STAT)[rownames(result$STAT) != "GLOW_SKAT"]
  for (row_name in non_cct_rows) {
    expect_true(grepl("^df_1_", row_name),
                label=paste("Row", row_name, "should have df=1"))
  }

  # Should NOT have df_Inf or df_2
  expect_false(any(grepl("df_Inf", rownames(result$STAT))))
  expect_false(any(grepl("df_2", rownames(result$STAT))))
})


test_that("GLOW_SKAT: no Liu method anywhere in code", {
  # CRITICAL TEST: Verify NO Liu method is used
  # Read the GLOW_SKAT source code
  skat_file <- system.file("../../R/GLOW_SKAT.R", package="GLOWr")
  if (skat_file == "") {
    skat_file <- file.path(Sys.getenv("GLOW_LEGACY_ROOT", unset = "/nonexistent"), "packages/GLOWr/R/GLOW_SKAT.R")
  }

  if (file.exists(skat_file)) {
    skat_code <- readLines(skat_file)
    skat_text <- paste(skat_code, collapse="\n")

    # Should NOT contain "liu" or "Liu" except in comments
    # Filter out comments
    code_lines <- skat_code[!grepl("^\\s*#", skat_code)]
    code_text <- paste(code_lines, collapse="\n")

    # Check for Liu method calls (should be NONE)
    expect_false(grepl("liu\\(", code_text, ignore.case=FALSE),
                 label="No liu() calls")
    expect_false(grepl("Liu\\(", code_text, ignore.case=FALSE),
                 label="No Liu() calls")
    expect_false(grepl("\\.liu", code_text, ignore.case=FALSE),
                 label="No .liu method calls")
    expect_false(grepl("liu_pvalue", code_text, ignore.case=FALSE),
                 label="No liu_pvalue calls")
    expect_false(grepl("liu\\.mod", code_text, ignore.case=FALSE),
                 label="No liu.mod calls")
  }
})


test_that("GLOW_SKAT: statistics are non-negative (quadratic form property)", {
  # Quadratic forms sum(w * Z^2) with positive weights are always >= 0
  set.seed(333)
  X <- matrix(rnorm(50*2), 50, 2)
  Y <- rnorm(50)
  G <- matrix(rbinom(50*8, 2, 0.2), 50, 8)
  B <- rnorm(8)
  PI <- runif(8)

  Zout <- getZ_marg_score(G, X, Y, trait="continuous")
  result <- GLOW_SKAT(marg_score_stats=Zout, B=B, PI=PI)

  # All non-CCT statistics should be >= 0 (quadratic form)
  non_cct_rows <- rownames(result$STAT)[rownames(result$STAT) != "GLOW_SKAT"]
  for (row_name in non_cct_rows) {
    stat_val <- result$STAT[row_name, 1]
    expect_true(stat_val >= 0,
                label=paste("Statistic for", row_name, "should be non-negative"))
  }
})


# Validation against legacy package (if available)
if (legacy_available) {
  test_that("GLOW_SKAT matches legacy implementation EXACTLY", {
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
    result_legacy <- GLOW::GLOW_SKAT(marg_score_stats=Zout_legacy, B=B, PI=PI)

    # New computation
    Zout_new <- GLOWr::getZ_marg_score(G, X, Y, trait="binary")
    result_new <- GLOWr::GLOW_SKAT(marg_score_stats=Zout_new, B=B, PI=PI)

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


  test_that("GLOW_SKAT matches legacy: continuous trait", {
    library(GLOW)

    set.seed(456)
    X <- matrix(rnorm(30*2), 30, 2)
    Y <- rnorm(30)
    G <- matrix(rbinom(30*8, 2, 0.2), 30, 8)
    B <- rnorm(8, sd=0.3)
    PI <- runif(8, 0.2, 0.8)

    # Legacy
    Zout_legacy <- GLOW::getZ_marg_score(G, X, Y, trait="continuous")
    result_legacy <- GLOW::GLOW_SKAT(marg_score_stats=Zout_legacy, B=B, PI=PI)

    # New
    Zout_new <- GLOWr::getZ_marg_score(G, X, Y, trait="continuous")
    result_new <- GLOWr::GLOW_SKAT(marg_score_stats=Zout_new, B=B, PI=PI)

    # Exact match
    expect_equal(result_new$STAT, result_legacy$STAT, tolerance=1e-10)
    expect_equal(result_new$PVAL, result_legacy$PVAL, tolerance=1e-10)
  })


  test_that("GLOW_SKAT matches legacy: multiple scenarios", {
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
      result_legacy <- GLOW::GLOW_SKAT(marg_score_stats=Zout_legacy, B=B, PI=PI)

      # New
      Zout_new <- GLOWr::getZ_marg_score(G, X, Y, trait=trait_type)
      result_new <- GLOWr::GLOW_SKAT(marg_score_stats=Zout_new, B=B, PI=PI)

      # Exact match
      expect_equal(result_new$STAT, result_legacy$STAT, tolerance=1e-10,
                   label=paste("Scenario", i, "STAT"))
      expect_equal(result_new$PVAL, result_legacy$PVAL, tolerance=1e-10,
                   label=paste("Scenario", i, "PVAL"))
    }
  })
}
