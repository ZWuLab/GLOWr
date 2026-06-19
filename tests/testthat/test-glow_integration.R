# ==============================================================================
# End-to-end Integration Tests for GLOW Testing Pipeline
# ==============================================================================
#
# Tests the full pipeline:
#   make_test_variant_set -> prepare_glow_input -> compute_score_stats
#   -> glow_test -> as.data.frame


# .make_test_variant_set() is defined in helper-test-data.R (auto-loaded)


# ---- Integration tests ----

test_that("full pipeline: prepare_glow_input -> compute_score_stats -> glow_test (binary)", {
  # 1. Create synthetic variant set
  vset <- .make_test_variant_set(n = 200, p = 20)

  # 2. Prepare phenotype and null model
  set.seed(123)
  Y <- rbinom(200, 1, 0.3)
  X <- cbind(1, rnorm(200))
  null_model <- fit_null_model(X, Y, trait = "binary")

  # 3. Prepare input
  calcB <- function(maf) rep(0.3, length(maf))
  input <- prepare_glow_input(vset, B_func = calcB,
                              PI = rep(0.5, 20), verbose = 0)
  expect_s3_class(input, "glow_input")

  # 4. Compute score stats
  stats <- compute_score_stats(input$G, null_model, verbose = 0)
  expect_true(is.list(stats))
  expect_equal(length(stats$Zscores), input$n_after_collapse)

  # 5. Run GLOW tests
  result <- glow_test(stats, input$B, input$PI,
                      region_info = input$region,
                      variant_summary = list(
                        n_original = input$n_original,
                        n_after_filter = input$n_after_filter,
                        n_after_ld = input$n_after_ld,
                        n_after_collapse = input$n_after_collapse,
                        cMAC = input$cMAC
                      ),
                      verbose = 0)
  expect_s3_class(result, "glow_test_result")
  expect_equal(length(result$pvalues), 16)
  expect_true(all(result$pvalues >= 0 & result$pvalues <= 1))

  # 6. Flatten to data.frame
  df <- as.data.frame(result)
  expect_equal(nrow(df), 1)
  expect_true("GLOW_Omni" %in% names(df))
  expect_true("label" %in% names(df))
})

test_that("full pipeline works with continuous trait", {
  vset <- .make_test_variant_set(n = 200, p = 15)

  set.seed(456)
  Y <- rnorm(200)
  X <- cbind(1, rnorm(200))
  null_model <- fit_null_model(X, Y, trait = "continuous")

  calcB <- function(maf) rep(0.3, length(maf))
  input <- prepare_glow_input(vset, B_func = calcB,
                              PI = rep(0.5, 15), verbose = 0)
  expect_s3_class(input, "glow_input")

  stats <- compute_score_stats(input$G, null_model, verbose = 0)
  result <- glow_test(stats, input$B, input$PI, verbose = 0)

  expect_s3_class(result, "glow_test_result")
  expect_equal(length(result$pvalues), 16)
  expect_true(all(result$pvalues >= 0 & result$pvalues <= 1))
})

test_that("multiple regions can be rbind-ed", {
  # Simulate analysis of two "genes" and verify results can be combined
  set.seed(789)
  n <- 200
  Y <- rbinom(n, 1, 0.3)
  X <- cbind(1, rnorm(n))
  null_model <- fit_null_model(X, Y, trait = "binary")
  calcB <- function(maf) rep(0.3, length(maf))

  results_list <- list()
  for (gene_i in 1:2) {
    vset <- .make_test_variant_set(n = n, p = 10 + gene_i * 2, seed = gene_i)
    input <- prepare_glow_input(vset, B_func = calcB,
                                PI = rep(0.5, vset$n_variants), verbose = 0)
    if (is.null(input)) next
    stats <- compute_score_stats(input$G, null_model, verbose = 0)
    result <- glow_test(stats, input$B, input$PI,
                        region_info = list(chr = "22",
                                           label = paste0("GENE", gene_i)),
                        variant_summary = list(
                          n_original = input$n_original,
                          n_after_collapse = input$n_after_collapse,
                          cMAC = input$cMAC
                        ),
                        verbose = 0)
    results_list[[gene_i]] <- as.data.frame(result)
  }

  combined <- do.call(rbind, results_list)
  expect_equal(nrow(combined), 2)
  expect_true("GLOW_Omni" %in% names(combined))
  expect_equal(combined$label[1], "GENE1")
  expect_equal(combined$label[2], "GENE2")
})

test_that("pipeline handles small variant set gracefully", {
  vset <- .make_test_variant_set(n = 100, p = 3, seed = 99)

  set.seed(100)
  Y <- rbinom(100, 1, 0.3)
  X <- cbind(1, rnorm(100))
  null_model <- fit_null_model(X, Y, trait = "binary")

  input <- prepare_glow_input(vset, B = rep(0.5, 3), PI = rep(0.3, 3),
                              mac_threshold = 0, verbose = 0)

  if (!is.null(input)) {
    stats <- compute_score_stats(input$G, null_model, verbose = 0)
    result <- glow_test(stats, input$B, input$PI, verbose = 0)
    expect_s3_class(result, "glow_test_result")
    expect_true(all(result$pvalues >= 0 & result$pvalues <= 1))

    # If single variant, verify all p-values are equal
    if (input$n_after_collapse == 1) {
      expect_true(length(unique(result$pvalues)) == 1)
    }
  }
})
