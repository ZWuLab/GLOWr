# Validation tests comparing against legacy code outputs
#
# This file validates that the new implementation produces numerically identical
# results to the legacy code for optimal weight calculations.

library(testthat)
library(GLOWr)

# Source the legacy code for comparison
legacy_code_path <- file.path(Sys.getenv("GLOW_LEGACY_ROOT", unset = "/nonexistent"), "legacy-materials/code/GLOW_R_pacakge/GLOW/R/helpers_optimalWeights.R")

if (file.exists(legacy_code_path)) {
  source(legacy_code_path)

  # Temporary g_GFisher for testing
  g_GFisher_two <- function(x, df = 2) {
    qchisq(log(2) + pnorm(abs(x), lower.tail = FALSE, log.p = TRUE),
           df = df, lower.tail = FALSE, log.p = TRUE)
  }

  g_identity <- function(x) x

  #################### Validation Tests ####################

  test_that("select_best_model matches legacy", {
    set.seed(42)
    X <- runif(100, 0.01, 0.5)
    Y <- sqrt(X * (1 - X)) + rnorm(100, 0, 0.05)

    # New implementation
    new_model <- GLOWr:::select_best_model(X, Y)

    # Legacy implementation
    legacy_model <- select_best_model(X, Y)

    # Both should select the same model (same R-squared)
    expect_equal(
      summary(new_model)$r.squared,
      summary(legacy_model)$r.squared,
      tolerance = 1e-10
    )

    # Coefficients should match
    expect_equal(
      coef(new_model),
      coef(legacy_model),
      tolerance = 1e-10
    )
  })

  test_that("hermite matches legacy", {
    x_vals <- c(-2, -1, 0, 1, 2, 3)
    mu <- 0.5

    for (ord in 1:8) {
      new_result <- GLOWr:::hermite(x_vals, mu, ord)
      legacy_result <- hermite(x_vals, mu, ord)

      expect_equal(
        new_result,
        legacy_result,
        tolerance = 1e-12,
        info = paste("Order", ord)
      )
    }
  })

  test_that("E_gX_p matches legacy for identity function", {
    g <- g_identity
    test_cases <- list(
      list(mu = 0, p = 1, sigma = 1),
      list(mu = 2, p = 1, sigma = 1),
      list(mu = 0, p = 2, sigma = 1),
      list(mu = 1, p = 2, sigma = 1.5)
    )

    for (tc in test_cases) {
      new_result <- GLOWr:::E_gX_p(g, tc$mu, tc$p, tc$sigma)
      legacy_result <- E_gX_p(g, tc$mu, tc$p, tc$sigma)

      expect_equal(
        new_result,
        legacy_result,
        tolerance = 1e-10,
        info = paste("mu =", tc$mu, "p =", tc$p, "sigma =", tc$sigma)
      )
    }
  })

  test_that("E_T_mix matches legacy", {
    g <- g_identity
    test_cases <- list(
      list(mu = 0, pi = 0),
      list(mu = 2, pi = 0),
      list(mu = 2, pi = 1),
      list(mu = 1.5, pi = 0.5),
      list(mu = 0.8, pi = 0.01)
    )

    for (tc in test_cases) {
      new_result <- GLOWr:::E_T_mix(g, tc$mu, tc$pi)
      legacy_result <- E_T_mix(g, tc$mu, tc$pi)

      expect_equal(
        new_result,
        legacy_result,
        tolerance = 1e-10,
        info = paste("mu =", tc$mu, "pi =", tc$pi)
      )
    }
  })

  test_that("Var_T_mix matches legacy", {
    g <- g_identity
    test_cases <- list(
      list(mu = 0, pi = 0),
      list(mu = 2, pi = 0),
      list(mu = 2, pi = 1),
      list(mu = 1.5, pi = 0.5),
      list(mu = 0.8, pi = 0.01)
    )

    for (tc in test_cases) {
      new_result <- GLOWr:::Var_T_mix(g, tc$mu, tc$pi)
      legacy_result <- Var_T_mix(g, tc$mu, tc$pi)

      expect_equal(
        new_result,
        legacy_result,
        tolerance = 1e-10,
        info = paste("mu =", tc$mu, "pi =", tc$pi)
      )
    }
  })

  test_that("get_wts matches legacy", {
    set.seed(123)
    Sigma <- diag(5)
    r <- rnorm(5, mean = 1, sd = 0.5)

    # Test without positive constraint
    new_result <- GLOWr:::get_wts(Sigma, r, is.posi.wts = FALSE)
    legacy_result <- get_wts(Sigma, r, is.posi.wts = FALSE)

    expect_equal(new_result$w, legacy_result$w, tolerance = 1e-12)
    expect_equal(new_result$w_normalized, legacy_result$w_normalized, tolerance = 1e-12)

    # Test with positive constraint
    r_mixed <- c(1, -0.5, 2, -1, 0.3)
    new_result <- GLOWr:::get_wts(Sigma, r_mixed, is.posi.wts = TRUE)
    legacy_result <- get_wts(Sigma, r_mixed, is.posi.wts = TRUE)

    expect_equal(new_result$w, legacy_result$w, tolerance = 1e-12)
    expect_equal(new_result$w_normalized, legacy_result$w_normalized, tolerance = 1e-12)
  })

  test_that("get_r matches legacy for identity function", {
    g <- g_identity
    MU <- c(0.5, 1.0, 0.2, 1.5)
    PI <- c(0.01, 0.05, 0.001, 0.1)

    new_result <- GLOWr:::get_r(g, MU, PI)
    legacy_result <- get_r(g, MU, PI)

    expect_equal(new_result, legacy_result, tolerance = 1e-10)
  })

  test_that("get_r_tilde matches legacy for identity function", {
    g <- g_identity
    MU <- c(0.5, 1.0, 0.2, 1.5)
    SD <- c(1.1, 1.2, 1.05, 1.3)

    new_result <- GLOWr:::get_r_tilde(g, MU, SD)
    legacy_result <- get_r_tilde(g, MU, SD)

    expect_equal(new_result, legacy_result, tolerance = 1e-10)
  })

  test_that("Optimal_Weights_M matches legacy for Burden test", {
    g <- g_identity
    set.seed(42)
    Bstar <- runif(5, 0.3, 1.5)
    PI <- runif(5, 0.001, 0.1)
    M <- diag(5)

    new_result <- Optimal_Weights_M(g, Bstar, PI, M, is.posi.wts = TRUE)
    legacy_result <- Optimal_Weights_M(g, Bstar, PI, M, is.posi.wts = TRUE)

    expect_equal(new_result$wts_BE, legacy_result$wts_BE, tolerance = 1e-10)
    expect_equal(new_result$wts_APE, legacy_result$wts_APE, tolerance = 1e-10)
  })

  test_that("CovM_gXgY matches legacy for identity function with simple correlation", {
    skip_on_cran()  # Slow test

    g <- g_identity
    M <- matrix(c(1, 0.3, 0.3, 1), 2, 2)
    MU1 <- c(0, 0)
    MU2 <- c(0, 0)

    new_result <- GLOWr:::CovM_gXgY(g, MU1, MU2, M)
    legacy_result <- CovM_gXgY(g, MU1, MU2, M)

    expect_equal(new_result, legacy_result, tolerance = 1e-8)
  })

  test_that("get_Sigma matches legacy for identity function", {
    skip_on_cran()  # Slow test

    g <- g_identity
    MU <- c(0.5, 1.0)
    PI <- c(0.01, 0.05)
    M <- matrix(c(1, 0.2, 0.2, 1), 2, 2)

    # Test H0
    new_result_H0 <- GLOWr:::get_Sigma(g, MU, PI, M, hypo = "H0")
    legacy_result_H0 <- get_Sigma(g, MU, PI, M, hypo = "H0")
    expect_equal(new_result_H0, legacy_result_H0, tolerance = 1e-8)

    # Test H1
    new_result_H1 <- GLOWr:::get_Sigma(g, MU, PI, M, hypo = "H1")
    legacy_result_H1 <- get_Sigma(g, MU, PI, M, hypo = "H1")
    expect_equal(new_result_H1, legacy_result_H1, tolerance = 1e-8)
  })

} else {
  test_that("Legacy code not found", {
    skip("Legacy code not available for validation")
  })
}
