# ==============================================================================
# Tests for variant_processing.R
# ==============================================================================
# Tests for post-extraction processing functions:
#   - flip_alleles()
#   - filter_variants_ld()
#   - collapse_rare_variants()
#   - aggregate_B_PI()


# ---- Phase 1: flip_alleles() ------------------------------------------------

test_that("flip_alleles corrects MAF > 0.5", {
  # 10 samples, 3 variants: first two need flipping
  G <- matrix(c(
    2, 2, 2, 2, 2, 2, 2, 2, 1, 0,  # AF = 0.85, should flip
    2, 2, 2, 2, 2, 2, 2, 2, 2, 1,  # AF = 0.95, should flip
    0, 0, 0, 0, 0, 0, 0, 0, 1, 2   # AF = 0.15, no flip
  ), nrow = 10, ncol = 3)

  G_flipped <- flip_alleles(G)
  # After flipping, all columns should have AF <= 0.5
  AF_after <- colMeans(G_flipped) / 2
  expect_true(all(AF_after <= 0.5))

  # Third column unchanged
  expect_equal(G_flipped[, 3], G[, 3])
  # First column: 2 -> 0, 1 -> 1, 0 -> 2
  expect_equal(G_flipped[1, 1], 0)
  expect_equal(G_flipped[9, 1], 1)
  expect_equal(G_flipped[10, 1], 2)
})

test_that("flip_alleles handles no-flip case", {
  G <- matrix(c(0, 0, 0, 1, 0, 0, 0, 0, 0, 0), nrow = 10, ncol = 1)
  expect_equal(flip_alleles(G), G)
})

test_that("flip_alleles handles all-flip case", {
  G <- matrix(2, nrow = 10, ncol = 3)
  G_flipped <- flip_alleles(G)
  expect_true(all(G_flipped == 0))
})

test_that("flip_alleles handles AF exactly 0.5", {
  # AF = 0.5 should NOT flip (only > 0.5 triggers flip)
  G <- matrix(c(0, 0, 0, 0, 0, 2, 2, 2, 2, 2), nrow = 10, ncol = 1)
  expect_equal(flip_alleles(G), G)
})

test_that("flip_alleles validates input", {
  expect_error(flip_alleles(c(1, 2, 3)))  # not a matrix
  expect_error(flip_alleles(matrix("a", 2, 2)))  # not numeric
})


# ---- Phase 2: filter_variants_ld() ------------------------------------------

test_that("filter_variants_ld keeps uncorrelated variants", {
  set.seed(42)
  G <- matrix(rnorm(100 * 5), nrow = 100, ncol = 5)  # independent columns
  keep <- filter_variants_ld(G, ld_threshold = 0.9)
  expect_equal(length(keep), 5)  # all kept
})

test_that("filter_variants_ld removes highly correlated variants", {
  set.seed(42)
  x <- rnorm(100)
  G <- cbind(x, x + rnorm(100, sd = 0.01), rnorm(100))  # cols 1-2 near-identical
  keep <- filter_variants_ld(G, ld_threshold = 0.9)
  expect_equal(length(keep), 2)  # one of the pair removed
  expect_true(3 %in% keep)  # independent column always kept
})

test_that("filter_variants_ld prefer_keep = lower_maf keeps rarer variant", {
  n <- 200
  # Two near-identical variants with different MAF
  rare <- rbinom(n, 2, 0.02)   # MAF ~ 0.02
  common <- rbinom(n, 2, 0.10) # MAF ~ 0.10
  # Make them highly correlated by copying + small noise
  G <- cbind(rare, rare + rbinom(n, 1, 0.001), common)
  # Cols 1 and 2 are highly correlated; col 1 is rarer
  keep <- filter_variants_ld(G, ld_threshold = 0.5, prefer_keep = "lower_maf")
  # Should keep col 1 (rarer) and col 3 (independent)
  expect_true(1 %in% keep)
  expect_true(3 %in% keep)
})

test_that("filter_variants_ld removes linearly dependent columns", {
  set.seed(42)
  G <- matrix(rnorm(100 * 3), nrow = 100)
  # Add a column that is an exact linear combination of the first two
  G <- cbind(G, G[, 1] + G[, 2])
  # LD threshold is high (no LD removal), but linear dependence should be caught
  keep <- filter_variants_ld(G, ld_threshold = 1.0, remove_lindep = TRUE)
  expect_equal(length(keep), 3)  # 4th column removed
  # Verify remaining columns are full rank
  expect_equal(qr(G[, keep])$rank, 3)
})

test_that("filter_variants_ld remove_lindep = FALSE skips QR check", {
  set.seed(42)
  G <- matrix(rnorm(100 * 3), nrow = 100)
  G <- cbind(G, G[, 1] + G[, 2])
  keep <- filter_variants_ld(G, ld_threshold = 1.0, remove_lindep = FALSE)
  expect_equal(length(keep), 4)  # linearly dependent column kept
})

test_that("filter_variants_ld handles single column", {
  G <- matrix(rnorm(100), nrow = 100, ncol = 1)
  expect_equal(filter_variants_ld(G), 1L)
})

test_that("filter_variants_ld handles zero-variance columns", {
  G <- cbind(rep(0, 100), rnorm(100), rep(1, 100))
  keep <- filter_variants_ld(G, ld_threshold = 0.9)
  # Only the variable column should be kept
  expect_true(2 %in% keep)
  expect_equal(length(keep), 1)
})

test_that("filter_variants_ld validates inputs", {
  expect_error(filter_variants_ld(c(1, 2, 3)))  # not a matrix
  expect_error(filter_variants_ld(matrix(1, 5, 5), ld_threshold = 0))  # threshold must be > 0
  expect_error(filter_variants_ld(matrix(1, 5, 5), prefer_keep = "invalid"))
})

test_that("filter_variants_ld greedy pruning removes most-connected variant first", {
  # Construct a scenario: variant 1 is correlated with both 2 and 3
  # but 2 and 3 are NOT correlated with each other.
  # Use orthogonal components: x1 = a+b, x2 ~ a, x3 ~ b
  # cor(x1,x2) ~ 1/sqrt(2) ~ 0.71, cor(x1,x3) ~ 0.71, cor(x2,x3) ~ 0
  set.seed(123)
  n <- 1000
  a <- rnorm(n)
  b <- rnorm(n)
  x1 <- a + b                     # correlated with both x2 and x3
  x2 <- a + rnorm(n, sd = 0.05)   # correlated with x1 via 'a', not with x3
  x3 <- b + rnorm(n, sd = 0.05)   # correlated with x1 via 'b', not with x2
  x4 <- rnorm(n)                  # independent
  G <- cbind(x1, x2, x3, x4)

  keep <- filter_variants_ld(G, ld_threshold = 0.5, prefer_keep = "first",
                              remove_lindep = FALSE)
  # x1 has 2 partners (x2 and x3), x2 and x3 each have 1 partner
  # Greedy should remove x1 first, then x2 and x3 are uncorrelated so both kept
  expect_false(1 %in% keep)
  expect_true(all(c(2, 3, 4) %in% keep))
})

test_that("filter_variants_ld prefer_keep = higher_maf keeps common variant", {
  set.seed(42)
  n <- 500
  # Create two highly correlated columns with different frequencies
  base <- rbinom(n, 2, 0.3)
  G <- cbind(base, base, rbinom(n, 2, 0.05))
  keep <- filter_variants_ld(G, ld_threshold = 0.5, prefer_keep = "higher_maf",
                              remove_lindep = FALSE)
  # Should keep the one with higher MAF and the independent column
  expect_equal(length(keep), 2)
  expect_true(3 %in% keep)
})


# ---- Phase 3: collapse_rare_variants() + aggregate_B_PI() -------------------

test_that("collapse_rare_variants collapses correctly", {
  n <- 50
  # 5 variants: 2 common (MAC=20), 3 rare (MAC=2)
  G <- matrix(0, nrow = n, ncol = 5)
  G[1:10, 1] <- 2  # MAC = 20 (common)
  G[1:10, 3] <- 2  # MAC = 20 (common)
  G[1, 2] <- 2      # MAC = 2 (rare)
  G[2, 4] <- 2      # MAC = 2 (rare)
  G[3, 5] <- 2      # MAC = 2 (rare)

  result <- collapse_rare_variants(G, mac_threshold = 10)

  # 2 common + 1 collapsed group = 3 output columns
  expect_equal(ncol(result$G_collapsed), 3)
  expect_equal(sum(result$is_collapsed), 1)
  # Collapsed column = sum of 3 rare columns
  collapsed_col <- which(result$is_collapsed)
  expect_equal(result$G_collapsed[1, collapsed_col], 2)  # from col 2
  expect_equal(result$G_collapsed[2, collapsed_col], 2)  # from col 4
  expect_equal(result$G_collapsed[3, collapsed_col], 2)  # from col 5
})

test_that("collapse_rare_variants treats singleton rare group as passthrough", {
  # 1 rare + 3 common. The rare variant forms a "group" of size 1 in
  # .build_collapsed_result. No aggregation actually happens, so is_collapsed
  # must be FALSE for that column, and n_after_collapse must equal ncol(G).
  n <- 50
  G <- matrix(0, nrow = n, ncol = 4)
  G[1:15, 1] <- 2  # MAC = 30 (common)
  G[1, 2]    <- 2  # MAC = 2  (rare; will be the singleton group)
  G[1:15, 3] <- 2  # MAC = 30 (common)
  G[1:15, 4] <- 2  # MAC = 30 (common)

  result <- collapse_rare_variants(G, mac_threshold = 10)

  expect_equal(ncol(result$G_collapsed), 4)
  expect_equal(sum(result$is_collapsed), 0)
  # The output column for the rare variant must equal the original column
  # (no rowSums aggregation happened).
  rare_out_col <- which(vapply(result$col_mapping,
                               function(m) identical(as.integer(m), 2L),
                               logical(1)))
  expect_length(rare_out_col, 1)
  expect_equal(result$G_collapsed[, rare_out_col], G[, 2])
})

test_that("collapse_rare_variants with spatial_grouping", {
  n <- 50
  # cols 1(rare), 2(common), 3(rare), 4(rare), 5(common), 6(rare)
  G <- matrix(0, nrow = n, ncol = 6)
  G[1:15, 2] <- 2  # common
  G[1:15, 5] <- 2  # common
  G[1, 1] <- 2     # rare, before common[2]
  G[1, 3] <- 2     # rare, between common[2] and common[5]
  G[2, 4] <- 2     # rare, between common[2] and common[5]
  G[3, 6] <- 2     # rare, after common[5]

  result <- collapse_rare_variants(G, mac_threshold = 10, spatial_grouping = TRUE)

  # Should have: rare(col1, singleton), common(col2),
  #              collapsed_group(col3+4), common(col5), rare(col6, singleton)
  # Only the size>=2 group counts as a true merger, so sum(is_collapsed) == 1.
  expect_equal(ncol(result$G_collapsed), 5)
  expect_equal(sum(result$is_collapsed), 1)
})

test_that("collapse_rare_variants with max_group_size", {
  n <- 50
  G <- matrix(0, nrow = n, ncol = 6)
  # All rare
  for (j in 1:6) G[j, j] <- 2

  result <- collapse_rare_variants(G, mac_threshold = 10, max_group_size = 2)
  # 6 rare variants, groups of 2 -> 3 collapsed groups
  expect_equal(ncol(result$G_collapsed), 3)
})

test_that("collapse_rare_variants returns unchanged if no rare", {
  n <- 50
  G <- matrix(0, nrow = n, ncol = 3)
  G[1:15, ] <- 2  # all MAC = 30

  result <- collapse_rare_variants(G, mac_threshold = 10)
  expect_equal(ncol(result$G_collapsed), 3)
  expect_true(all(!result$is_collapsed))
})

test_that("collapse_rare_variants returns unchanged if all common", {
  n <- 50
  G <- matrix(0, nrow = n, ncol = 2)
  G[1:20, ] <- 2  # MAC = 40 > 10
  result <- collapse_rare_variants(G, mac_threshold = 10)
  expect_equal(result$G_collapsed, G)
  expect_equal(result$col_mapping, list(1L, 2L))
  expect_equal(result$is_collapsed, c(FALSE, FALSE))
})

test_that("collapse_rare_variants handles all rare case", {
  n <- 50
  G <- matrix(0, nrow = n, ncol = 4)
  G[1, 1] <- 2
  G[2, 2] <- 2
  G[3, 3] <- 2
  G[4, 4] <- 2
  result <- collapse_rare_variants(G, mac_threshold = 10)
  # All 4 rare variants -> 1 collapsed group
  expect_equal(ncol(result$G_collapsed), 1)
  expect_true(result$is_collapsed[1])
  # Collapsed column = rowSums of all 4
  expect_equal(result$G_collapsed[1, 1], 2)
  expect_equal(result$G_collapsed[2, 1], 2)
  expect_equal(result$G_collapsed[5, 1], 0)
})

test_that("collapse_rare_variants preserves positional order", {
  n <- 50
  # cols: rare(1), common(2), rare(3), common(4), rare(5)
  G <- matrix(0, nrow = n, ncol = 5)
  G[1:15, 2] <- 2  # common
  G[1:15, 4] <- 2  # common
  G[1, 1] <- 2     # rare
  G[1, 3] <- 2     # rare
  G[1, 5] <- 2     # rare

  result <- collapse_rare_variants(G, mac_threshold = 10)
  # Without spatial_grouping, all rare go into one group
  # Collapsed group position = min(1,3,5) = 1, then common 2, then common 4
  # So order: collapsed(1,3,5), common(2), common(4)
  expect_equal(ncol(result$G_collapsed), 3)
  # First column should be collapsed (position 1 from rare group)
  expect_true(result$is_collapsed[1])
  # col_mapping for the collapsed should contain 1, 3, 5
  expect_equal(sort(result$col_mapping[[1]]), c(1L, 3L, 5L))
})

test_that("collapse_rare_variants validates inputs", {
  expect_error(collapse_rare_variants(c(1, 2, 3)))  # not a matrix
  expect_error(collapse_rare_variants(matrix(1, 5, 5), agg_method = "invalid"))
})

test_that("aggregate_B_PI works correctly", {
  B <- c(1.0, 0.5, 0.8, 0.3, 0.9)
  PI <- c(0.1, 0.2, 0.3, 0.4, 0.5)

  collapse_result <- list(
    col_mapping = list(1, c(2, 3), 4, 5),  # col 2+3 collapsed
    is_collapsed = c(FALSE, TRUE, FALSE, FALSE),
    agg_method = "mean"
  )

  agg <- aggregate_B_PI(B, PI, collapse_result)
  expect_equal(agg$B_collapsed[1], 1.0)
  expect_equal(agg$B_collapsed[2], mean(c(0.5, 0.8)))
  expect_equal(agg$PI_collapsed[2], mean(c(0.2, 0.3)))
})

test_that("aggregate_B_PI respects agg_method = max", {
  B <- c(1.0, 0.5, 0.8)
  PI <- c(0.1, 0.2, 0.3)

  collapse_result <- list(
    col_mapping = list(1, c(2, 3)),
    is_collapsed = c(FALSE, TRUE),
    agg_method = "max"
  )

  agg <- aggregate_B_PI(B, PI, collapse_result)
  expect_equal(agg$B_collapsed[2], 0.8)
  expect_equal(agg$PI_collapsed[2], 0.3)
})

test_that("aggregate_B_PI respects agg_method = sum", {
  B <- c(1.0, 0.5, 0.8)
  PI <- c(0.1, 0.2, 0.3)

  collapse_result <- list(
    col_mapping = list(1, c(2, 3)),
    is_collapsed = c(FALSE, TRUE),
    agg_method = "sum"
  )

  agg <- aggregate_B_PI(B, PI, collapse_result)
  expect_equal(agg$B_collapsed[2], 0.5 + 0.8)
  expect_equal(agg$PI_collapsed[2], 0.2 + 0.3)
})

test_that("aggregate_B_PI with no collapsing returns original values", {
  B <- c(1.0, 2.0, 3.0)
  PI <- c(0.1, 0.2, 0.3)

  collapse_result <- list(
    col_mapping = list(1, 2, 3),
    is_collapsed = c(FALSE, FALSE, FALSE),
    agg_method = "mean"
  )

  agg <- aggregate_B_PI(B, PI, collapse_result)
  expect_equal(agg$B_collapsed, B)
  expect_equal(agg$PI_collapsed, PI)
})
