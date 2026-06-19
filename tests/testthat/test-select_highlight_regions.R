# ==============================================================================
# Tests for select_highlight_regions.R
# ==============================================================================
#
# Tests for select_highlight_regions() (threshold selection + min_n top-up).


make_df <- function(p) {
  data.frame(label = paste0("r", seq_along(p)),
             pval = p, stringsAsFactors = FALSE)
}

test_that("select_highlight_regions returns regions at or below threshold", {
  df <- make_df(c(1e-9, 1e-8, 0.5, 0.6, 0.7))
  res <- select_highlight_regions(df, "pval", threshold = 1e-6, min_n = 1L)
  # only r1, r2 pass; min_n = 1 does not force a top-up
  expect_equal(res, c("r1", "r2"))
})

test_that("select_highlight_regions tops up to min_n when too few pass", {
  df <- make_df(c(1e-9, 0.2, 0.3, 0.4, 0.5))
  res <- select_highlight_regions(df, "pval", threshold = 1e-6, min_n = 3L)
  # only r1 passes, top up to 3 by ascending p
  expect_equal(length(res), 3L)
  expect_equal(res, c("r1", "r2", "r3"))
})

test_that("select_highlight_regions returns ascending-p order", {
  df <- data.frame(label = c("a", "b", "c"),
                   pval = c(0.3, 0.1, 0.2), stringsAsFactors = FALSE)
  res <- select_highlight_regions(df, "pval", threshold = 1, min_n = 3L)
  expect_equal(res, c("b", "c", "a"))
})

test_that("select_highlight_regions caps at the number of non-NA p-values", {
  df <- make_df(c(0.01, NA, NA, NA, NA))
  res <- select_highlight_regions(df, "pval", threshold = 1e-6, min_n = 10L)
  # only one non-NA p-value, so capped at 1
  expect_equal(length(res), 1L)
  expect_equal(res, "r1")
})

test_that("select_highlight_regions returns empty when no usable p-values", {
  df <- make_df(rep(NA_real_, 4))
  res <- select_highlight_regions(df, "pval", threshold = 1e-6, min_n = 5L)
  expect_equal(res, character(0))
})

test_that("select_highlight_regions honors a custom label_col", {
  df <- data.frame(gene = c("X", "Y"), pval = c(0.001, 0.5),
                   stringsAsFactors = FALSE)
  res <- select_highlight_regions(df, "pval", threshold = 0.01,
                                  label_col = "gene", min_n = 1L)
  expect_equal(res, "X")
})
