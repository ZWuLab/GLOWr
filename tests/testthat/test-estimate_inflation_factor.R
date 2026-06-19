# ==============================================================================
# Tests for estimate_inflation_factor()
# ==============================================================================
#
# The "scan -> factor" wrapper: method dispatch (ldsc_intercept | lambda_gc),
# the MAF / MHC QC and scan<->ld_scores merge, consistency with the underlying
# ldsc_regression()/compute_lambda_gc(), the uniform return schema, and input
# validation.

.make_scan <- function(n = 6000L, seed = 1L, intercept = 1.1, slope = 5e-4) {
  set.seed(seed)
  chr <- as.character(rep(1:22, length.out = n))
  pos <- as.integer(seq_len(n) * 1000L)               # all < 25 Mb (no MHC by default)
  ld  <- stats::rgamma(n, shape = 2, scale = 30)
  chi2 <- (intercept + slope * ld) * stats::rchisq(n, df = 1)
  Z   <- sqrt(chi2) * sample(c(-1, 1), n, replace = TRUE)
  scan <- data.frame(chr = chr, pos = pos, ref = "A", alt = "G",
                     MAF = stats::runif(n, 0.001, 0.5), Z = Z,
                     pvalue = stats::pchisq(chi2, df = 1, lower.tail = FALSE),
                     stringsAsFactors = FALSE)
  ld_scores <- data.frame(chr = chr, pos = pos, ref = "A", alt = "G", ld = ld,
                          stringsAsFactors = FALSE)
  list(scan = scan, ld_scores = ld_scores)
}

.COLS <- c("method", "factor", "n_variants", "intercept", "intercept_se",
           "slope", "slope_se", "confounding_ratio", "ratio_se",
           "lambda_gc", "mean_chi2")


test_that("uniform 1-row schema for both methods", {
  d <- .make_scan()
  a <- estimate_inflation_factor(d$scan, "ldsc_intercept", ld_scores = d$ld_scores)
  b <- estimate_inflation_factor(d$scan, "lambda_gc")
  expect_identical(names(a), .COLS); expect_identical(names(b), .COLS)
  expect_equal(nrow(a), 1L); expect_equal(nrow(b), 1L)
  expect_equal(a$method, "ldsc_intercept"); expect_equal(b$method, "lambda_gc")
  expect_true(is.na(b$intercept) && is.na(b$confounding_ratio))   # ldsc fields NA
})

test_that("ldsc_intercept matches a manual QC + merge + ldsc_regression", {
  d <- .make_scan()
  est <- estimate_inflation_factor(d$scan, "ldsc_intercept", ld_scores = d$ld_scores)
  # manual replication of the wrapper's QC + merge
  s <- d$scan[d$scan$MAF >= 0.01, ]
  k <- function(x) paste(x$chr, x$pos, x$ref, x$alt, sep = ":")
  s$key <- k(s); L <- d$ld_scores; L$key <- k(L)
  m <- merge(data.frame(key = s$key, chr = s$chr, chi2 = s$Z^2),
             data.frame(key = L$key, ld = L$ld), by = "key")
  manual <- ldsc_regression(m$chi2, m$ld, m$chr)
  expect_equal(est$factor, manual$intercept)
  expect_equal(est$factor, est$intercept)
  expect_equal(est$n_variants, manual$n_variants)
  expect_equal(est$confounding_ratio, manual$confounding_ratio)
})

test_that("lambda_gc matches compute_lambda_gc on the MAF-QC'd p-values", {
  d <- .make_scan()
  est <- estimate_inflation_factor(d$scan, "lambda_gc")
  s   <- d$scan[d$scan$MAF >= 0.01, ]
  expect_equal(est$factor, compute_lambda_gc(s$pvalue))
  expect_equal(est$factor, est$lambda_gc)
})

test_that("MAF QC reduces the variant count; maf_min = 0 disables it", {
  d <- .make_scan()
  n_qc  <- estimate_inflation_factor(d$scan, "lambda_gc", maf_min = 0.01)$n_variants
  n_all <- estimate_inflation_factor(d$scan, "lambda_gc", maf_min = 0)$n_variants
  expect_lt(n_qc, n_all)
})

test_that("exclude_mhc drops MHC variants for ldsc_intercept", {
  d <- .make_scan()
  d$scan$chr[1:200] <- "6"; d$scan$pos[1:200] <- as.integer(30e6)   # into the MHC
  d$ld_scores$chr[1:200] <- "6"; d$ld_scores$pos[1:200] <- as.integer(30e6)
  with_mhc <- estimate_inflation_factor(d$scan, "ldsc_intercept",
                                        ld_scores = d$ld_scores, exclude_mhc = FALSE)
  no_mhc   <- estimate_inflation_factor(d$scan, "ldsc_intercept",
                                        ld_scores = d$ld_scores, exclude_mhc = TRUE)
  expect_gt(with_mhc$n_variants, no_mhc$n_variants)
})

test_that("input validation", {
  d <- .make_scan()
  expect_error(estimate_inflation_factor(d$scan, "ldsc_intercept"), "ld_scores")
  expect_error(estimate_inflation_factor(d$scan, "lambda_gc",
               pvalue_col = "nope"), "no 'nope' column")
  expect_error(estimate_inflation_factor(list(), "lambda_gc"), "data.frame")
})
