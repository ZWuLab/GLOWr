# ==============================================================================
# Tests for glow_test_families.R
# ==============================================================================
#
# Tests for glow_family_of() and group_glow_tests_by_family().


# ==================== glow_family_of() ====================

test_that("glow_family_of resolves the first token after GLOW_", {
  expect_equal(glow_family_of("GLOW_SKAT_BE_N"), "SKAT")
  expect_equal(glow_family_of("GLOW_Burden_equ"), "Burden")
  expect_equal(glow_family_of("GLOW_Fisher_APR_sparse"), "Fisher")
})

test_that("glow_family_of buckets the three omnibus columns into Omni", {
  expect_equal(
    glow_family_of(c("GLOW_BSF_Omni", "GLOW_SNV_CCT", "GLOW_Omni")),
    c("Omni", "Omni", "Omni"))
})

test_that("glow_family_of returns NA for non-GLOW names without overrides", {
  expect_true(is.na(glow_family_of("STAAR_O")))
})

test_that("glow_family_of is vectorised and order-preserving", {
  res <- glow_family_of(c("GLOW_SKAT_BE_N", "STAAR_O", "GLOW_Omni"))
  expect_equal(res, c("SKAT", NA, "Omni"))
})

test_that("glow_family_of applies family_overrides for non-GLOW names", {
  expect_equal(
    glow_family_of("STAAR_O", family_overrides = c(STAAR_O = "Omni")),
    "Omni")
})

test_that("glow_family_of rejects an unnamed family_overrides", {
  expect_error(
    glow_family_of("STAAR_O", family_overrides = "Omni"),
    "named character vector")
})


# ==================== group_glow_tests_by_family() ====================

test_that("group_glow_tests_by_family groups tests into a named family list", {
  tests <- c("GLOW_SKAT_BE_N", "GLOW_SKAT_equ",
             "GLOW_Burden_equ", "GLOW_Fisher_equ", "GLOW_Omni")
  fams <- group_glow_tests_by_family(tests)
  expect_type(fams, "list")
  expect_equal(fams$SKAT, c("GLOW_SKAT_BE_N", "GLOW_SKAT_equ"))
  expect_equal(fams$Burden, "GLOW_Burden_equ")
  expect_equal(fams$Fisher, "GLOW_Fisher_equ")
  expect_equal(fams$Omni, "GLOW_Omni")
})

test_that("group_glow_tests_by_family honors priority order and omni_last", {
  tests <- c("GLOW_Omni", "GLOW_Fisher_equ", "GLOW_SKAT_equ",
             "GLOW_Burden_equ")
  fams <- group_glow_tests_by_family(tests, omni_last = TRUE)
  # priority SKAT, Burden, Fisher first; Omni forced last
  expect_equal(names(fams), c("SKAT", "Burden", "Fisher", "Omni"))
})

test_that("group_glow_tests_by_family with omni_last=FALSE keeps Omni in order", {
  tests <- c("GLOW_SKAT_equ", "GLOW_Omni", "GLOW_Burden_equ")
  fams <- group_glow_tests_by_family(
    tests, priority = c("SKAT", "Burden", "Omni"), omni_last = FALSE)
  expect_equal(names(fams), c("SKAT", "Burden", "Omni"))
})

test_that("group_glow_tests_by_family errors on unrecognised names", {
  expect_error(
    group_glow_tests_by_family(c("GLOW_SKAT_equ", "STAAR_O")),
    "Unrecognised test column name")
})

test_that("group_glow_tests_by_family routes overrides into a family", {
  tests <- c("GLOW_SKAT_equ", "STAAR_O")
  fams <- group_glow_tests_by_family(
    tests, family_overrides = c(STAAR_O = "Omni"))
  expect_true("STAAR_O" %in% fams$Omni)
})
