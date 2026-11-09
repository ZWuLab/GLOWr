# ==============================================================================
# Tests for glow_preprocessing_defaults()

test_that("glow_preprocessing_defaults() is the one source of the function defaults", {
  d <- glow_preprocessing_defaults()
  expect_identical(d, list(ld_threshold = 0.9, mac_threshold = 10L, collapse_method = "mean"))
  expect_identical(eval(formals(prepare_glow_input)$ld_threshold), d$ld_threshold)
  expect_identical(eval(formals(prepare_glow_input)$mac_threshold), d$mac_threshold)
  expect_identical(eval(formals(prepare_glow_input)$collapse_method), d$collapse_method)
  expect_identical(eval(formals(glow_region_test)$ld_threshold), d$ld_threshold)
  expect_identical(eval(formals(glow_region_test)$mac_threshold), d$mac_threshold)
  expect_identical(eval(formals(filter_variants_ld)$ld_threshold), d$ld_threshold)
  expect_identical(eval(formals(collapse_rare_variants)$mac_threshold), d$mac_threshold)
  expect_identical(eval(formals(collapse_rare_variants)$agg_method), d$collapse_method)
})
