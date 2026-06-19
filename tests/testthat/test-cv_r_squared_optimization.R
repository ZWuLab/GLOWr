# Test validation for PRESS-based cv_r_squared optimization
#
# This test verifies that the optimized PRESS-based implementation
# produces equivalent results to the original loop-based LOOCV implementation

test_that("PRESS-based cv_r_squared matches loop-based implementation", {

  # Reference implementation: Original loop-based LOOCV
  cv_r_squared_loop <- function(X, Y, model_formula, verbose = 0) {
    n <- length(Y)
    predictions <- numeric(n)

    # Create data frame for modeling
    if (is.matrix(X) || is.data.frame(X)) {
      model_data <- cbind(Y = Y, X)
    } else {
      model_data <- data.frame(Y = Y, X = X)
    }

    # Perform leave-one-out CV
    for (i in 1:n) {
      train_data <- model_data[-i, ]
      test_data <- model_data[i, , drop = FALSE]

      tryCatch({
        cv_model <- lm(model_formula, data = train_data)
        predictions[i] <- predict(cv_model, newdata = test_data)
      }, error = function(e) {
        predictions[i] <- mean(train_data$Y)
      })
    }

    # Calculate R-squared
    ss_res <- sum((Y - predictions)^2)
    ss_tot <- sum((Y - mean(Y))^2)
    r_squared <- 1 - (ss_res / ss_tot)

    # Handle edge cases
    if (is.na(r_squared) || is.infinite(r_squared)) {
      r_squared <- 0
    }

    return(r_squared)
  }

  # Test 1: Simple linear relationship
  set.seed(123)
  n <- 50
  X <- runif(n, 0.01, 0.5)
  Y <- 2 * X + rnorm(n, 0, 0.1)

  cv_old <- cv_r_squared_loop(X, Y, Y ~ X)
  cv_new <- cv_r_squared(X, Y, Y ~ X)

  expect_equal(cv_new, cv_old, tolerance = 1e-10,
               info = "Simple linear model: PRESS should match loop-based LOOCV")

  # Test 2: Quadratic relationship
  set.seed(456)
  X <- runif(n, 0.01, 0.5)
  Y <- X * (1 - X) + rnorm(n, 0, 0.05)

  cv_old <- cv_r_squared_loop(X, Y, Y ~ X + I(X^2))
  cv_new <- cv_r_squared(X, Y, Y ~ X + I(X^2))

  expect_equal(cv_new, cv_old, tolerance = 1e-10,
               info = "Quadratic model: PRESS should match loop-based LOOCV")

  # Test 3: Larger sample size
  set.seed(789)
  n <- 100
  X <- runif(n, 0.01, 0.5)
  Y <- X * (1 - X) + rnorm(n, 0, 0.1)

  cv_old <- cv_r_squared_loop(X, Y, Y ~ X + I(X^2))
  cv_new <- cv_r_squared(X, Y, Y ~ X + I(X^2))

  expect_equal(cv_new, cv_old, tolerance = 1e-10,
               info = "Larger sample (n=100): PRESS should match loop-based LOOCV")

  # Test 4: Very noisy data (low R²)
  set.seed(321)
  n <- 50
  X <- runif(n, 0.01, 0.5)
  Y <- X + rnorm(n, 0, 1)  # High noise

  cv_old <- cv_r_squared_loop(X, Y, Y ~ X)
  cv_new <- cv_r_squared(X, Y, Y ~ X)

  expect_equal(cv_new, cv_old, tolerance = 1e-10,
               info = "Noisy data: PRESS should match loop-based LOOCV")

  # Test 5: Near-perfect fit (high R²)
  set.seed(654)
  n <- 50
  X <- runif(n, 0.01, 0.5)
  Y <- 2 * X + rnorm(n, 0, 0.001)  # Very low noise

  cv_old <- cv_r_squared_loop(X, Y, Y ~ X)
  cv_new <- cv_r_squared(X, Y, Y ~ X)

  expect_equal(cv_new, cv_old, tolerance = 1e-10,
               info = "Near-perfect fit: PRESS should match loop-based LOOCV")
})

test_that("PRESS optimization provides expected speedup", {
  skip_on_cran()  # Performance tests can be skipped on CRAN
  skip_on_ci()    # absolute wall-clock bound is unreliable on loaded CI runners;
                  # this is a pure perf check (no correctness assertion). Local only.

  set.seed(999)
  n <- 100
  X <- runif(n, 0.01, 0.5)
  Y <- X * (1 - X) + rnorm(n, 0, 0.1)

  # Time the optimized version
  start_time <- Sys.time()
  for (i in 1:10) {
    cv_r_squared(X, Y, Y ~ X + I(X^2))
  }
  time_optimized <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

  # Verify it completes in reasonable time
  # 10 iterations with n=100 should take < 0.5 seconds with PRESS optimization
  expect_lt(time_optimized, 0.5,
            label = paste("Expected < 0.5s for 10 iterations, got",
                        round(time_optimized, 3), "s"))

  # Log the performance
  message("PRESS optimization: 10 iterations with n=100 took ",
          round(time_optimized, 3), " seconds")
})

test_that("cv_r_squared handles edge cases correctly", {

  # Edge case 1: Small sample size
  set.seed(111)
  X <- c(0.1, 0.2, 0.3)
  Y <- c(0.5, 0.6, 0.7)

  cv_result <- cv_r_squared(X, Y, Y ~ X)
  expect_true(is.numeric(cv_result))
  expect_true(cv_result >= 0 && cv_result <= 1)

  # Edge case 2: Perfect fit (but should still work)
  X <- 1:10
  Y <- 2 * X  # Perfect linear relationship

  cv_result <- cv_r_squared(X, Y, Y ~ X)
  expect_true(is.numeric(cv_result))
  expect_true(cv_result >= 0.9)  # Should have very high R²

  # Edge case 3: High leverage points
  set.seed(222)
  X <- c(runif(45, 0, 1), rep(10, 5))  # 5 extreme leverage points
  Y <- 2 * X + rnorm(50, 0, 0.5)

  # Should complete without error despite high leverage
  expect_no_error(cv_r_squared(X, Y, Y ~ X))
  cv_result <- cv_r_squared(X, Y, Y ~ X)
  expect_true(is.numeric(cv_result))
  expect_true(cv_result >= 0)
})
