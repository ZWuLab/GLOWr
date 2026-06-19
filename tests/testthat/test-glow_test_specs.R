# ==============================================================================
# Tests for glow_test_specs.R
# ==============================================================================
#
# Tests for default_test_specs(), validate_test_spec(), make_test_spec(),
# and .default_equal_weights().


# ==================== default_test_specs() ====================

test_that("default_test_specs returns a list of 3 specs", {
  specs <- default_test_specs()
  expect_type(specs, "list")
  expect_equal(length(specs), 3)
})

test_that("default_test_specs has correct family order: SKAT, Burden, Fisher", {
  specs <- default_test_specs()
  expect_equal(specs[[1]]$family, "SKAT")
  expect_equal(specs[[2]]$family, "Burden")
  expect_equal(specs[[3]]$family, "Fisher")
})

test_that("default_test_specs has correct df values", {
  specs <- default_test_specs()
  expect_equal(specs[[1]]$df, 1)      # SKAT
  expect_equal(specs[[2]]$df, Inf)    # Burden
  expect_equal(specs[[3]]$df, 2)      # Fisher
})

test_that("default_test_specs g functions match expected definitions", {
  specs <- default_test_specs()

  # SKAT: g(x) = x^2
  x_test <- c(-2, -1, 0, 1, 2)
  expect_equal(specs[[1]]$g(x_test), x_test^2)

  # Burden: g(x) = x (identity)
  expect_equal(specs[[2]]$g(x_test), x_test)

  # Fisher: g(x, df=2) = g_GFisher_two(x, df)
  # The Fisher g function must have a df formal with default 2
  expect_true("df" %in% names(formals(specs[[3]]$g)))
  expect_equal(formals(specs[[3]]$g)$df, 2)
  expect_equal(specs[[3]]$g(x_test), g_GFisher_two(x_test, df = 2))
})

test_that("default_test_specs includes weight_fns and p.type", {
  specs <- default_test_specs()
  for (i in seq_along(specs)) {
    expect_type(specs[[i]]$weight_fns, "list")
    expect_true(length(specs[[i]]$weight_fns) >= 2)
    expect_true("optimal" %in% names(specs[[i]]$weight_fns))
    expect_true("equal" %in% names(specs[[i]]$weight_fns))
    expect_equal(specs[[i]]$p.type, "two")
  }
})


# ==================== .default_equal_weights() ====================

test_that(".default_equal_weights returns named list with unit vector", {
  result <- GLOWr:::.default_equal_weights(5)
  expect_type(result, "list")
  expect_named(result, "equ")
  expect_equal(result$equ, rep(1, 5))
})

test_that(".default_equal_weights works for p=1", {
  result <- GLOWr:::.default_equal_weights(1)
  expect_equal(result$equ, 1)
})


# ==================== validate_test_spec() ====================

test_that("validate_test_spec applies defaults for missing weight_fns and p.type", {
  spec <- list(family = "Custom", g = function(x) x^3, df = 3)
  validated <- GLOWr:::validate_test_spec(spec)

  expect_type(validated$weight_fns, "list")
  expect_true("optimal" %in% names(validated$weight_fns))
  expect_true("equal" %in% names(validated$weight_fns))
  expect_equal(validated$p.type, "two")
})

test_that("validate_test_spec preserves user-supplied weight_fns and p.type", {
  custom_wfn <- function(p) list(custom = rep(2, p))
  spec <- list(
    family = "Custom",
    g = function(x) x^3,
    df = 3,
    weight_fns = list(mine = custom_wfn),
    p.type = "one"
  )
  validated <- GLOWr:::validate_test_spec(spec)

  expect_named(validated$weight_fns, "mine")
  expect_equal(validated$p.type, "one")
})

test_that("validate_test_spec rejects missing family", {
  expect_error(
    GLOWr:::validate_test_spec(list(g = function(x) x, df = 1)),
    "family.*non-empty character"
  )
})

test_that("validate_test_spec rejects empty family string", {
  expect_error(
    GLOWr:::validate_test_spec(list(family = "", g = function(x) x, df = 1)),
    "family.*non-empty character"
  )
})

test_that("validate_test_spec rejects non-function g", {
  expect_error(
    GLOWr:::validate_test_spec(list(family = "X", g = "not_a_function", df = 1)),
    "g.*must be a function"
  )
})

test_that("validate_test_spec rejects missing g", {
  expect_error(
    GLOWr:::validate_test_spec(list(family = "X", df = 1)),
    "g.*must be a function"
  )
})

test_that("validate_test_spec rejects non-numeric df", {
  expect_error(
    GLOWr:::validate_test_spec(list(family = "X", g = function(x) x, df = "two")),
    "df.*must be a numeric scalar"
  )
})

test_that("validate_test_spec rejects df <= 0", {
  expect_error(
    GLOWr:::validate_test_spec(list(family = "X", g = function(x) x, df = 0)),
    "df.*must be a numeric scalar > 0"
  )
  expect_error(
    GLOWr:::validate_test_spec(list(family = "X", g = function(x) x, df = -1)),
    "df.*must be a numeric scalar > 0"
  )
})

test_that("validate_test_spec accepts Inf df", {
  validated <- GLOWr:::validate_test_spec(
    list(family = "X", g = function(x) x, df = Inf)
  )
  expect_equal(validated$df, Inf)
})

test_that("validate_test_spec rejects unnamed weight_fns", {
  expect_error(
    GLOWr:::validate_test_spec(list(
      family = "X", g = function(x) x, df = 1,
      weight_fns = list(function(p) list(a = rep(1, p)))
    )),
    "weight_fns.*must be a named list"
  )
})

test_that("validate_test_spec rejects non-function weight_fns entry", {
  expect_error(
    GLOWr:::validate_test_spec(list(
      family = "X", g = function(x) x, df = 1,
      weight_fns = list(bad = "not_a_function")
    )),
    "must be a function"
  )
})


# ==================== make_test_spec() ====================

test_that("make_test_spec creates and validates a spec", {
  spec <- GLOWr:::make_test_spec(
    family = "MyTest",
    g = function(x) x^2,
    df = 1
  )
  expect_equal(spec$family, "MyTest")
  expect_equal(spec$df, 1)
  expect_equal(spec$p.type, "two")  # default
  expect_true("optimal" %in% names(spec$weight_fns))
})

test_that("make_test_spec passes through optional args", {
  custom_wfn <- function(p) list(custom = rep(3, p))
  spec <- GLOWr:::make_test_spec(
    family = "Custom",
    g = function(x) x,
    df = Inf,
    weight_fns = list(mine = custom_wfn),
    p.type = "one"
  )
  expect_named(spec$weight_fns, "mine")
  expect_equal(spec$p.type, "one")
})

test_that("make_test_spec rejects invalid input", {
  expect_error(
    GLOWr:::make_test_spec(family = "", g = function(x) x, df = 1),
    "family"
  )
})
