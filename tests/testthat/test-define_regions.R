# Tests for define_regions.R

test_that("define_regions_gene uses bundled data by default", {
  regions <- define_regions_gene(chr = 22)
  expect_true(nrow(regions) > 400)
  expect_equal(names(regions), c("region_id", "chr", "start", "end", "label"))
  expect_true(all(regions$chr == "22"))
})

test_that("define_regions_gene accepts STAARpipeline format", {
  custom_gi <- data.frame(
    hgnc_symbol = c("GENE_A", "GENE_B"),
    chromosome_name = c(1, 1),
    start_position = c(1000, 5000),
    end_position = c(2000, 6000)
  )
  regions <- define_regions_gene(custom_gi)
  expect_equal(nrow(regions), 2)
  expect_equal(regions$region_id, c("GENE_A", "GENE_B"))
})

test_that("define_regions_gene accepts simple 4-column format", {
  simple_format <- data.frame(
    V1 = c(1, 1), V2 = c(1000, 5000),
    V3 = c(2000, 6000), V4 = c("GENE_A", "GENE_B"),
    stringsAsFactors = FALSE
  )
  regions <- define_regions_gene(simple_format)
  expect_equal(nrow(regions), 2)
  expect_equal(regions$region_id, c("GENE_A", "GENE_B"))
})

test_that("define_regions_gene applies extension", {
  gi <- data.frame(
    hgnc_symbol = "TEST", chromosome_name = 1,
    start_position = 1000, end_position = 2000
  )
  regions <- define_regions_gene(gi, extend = 500)
  expect_equal(regions$start, 500)
  expect_equal(regions$end, 2500)
})

test_that("define_regions_gene extension clamps at 1", {
  gi <- data.frame(
    hgnc_symbol = "TEST", chromosome_name = 1,
    start_position = 100, end_position = 2000
  )
  regions <- define_regions_gene(gi, extend = 200)
  expect_equal(regions$start, 1)
})

test_that("define_regions_gene filters by chromosome", {
  gi <- data.frame(
    hgnc_symbol = c("A", "B", "C"),
    chromosome_name = c(1, 2, 1),
    start_position = c(100, 200, 300),
    end_position = c(200, 300, 400)
  )
  regions <- define_regions_gene(gi, chr = 1)
  expect_equal(nrow(regions), 2)
})

test_that("define_regions_gene warns for empty chromosome", {
  gi <- data.frame(
    hgnc_symbol = "A", chromosome_name = 1,
    start_position = 100, end_position = 200
  )
  expect_warning(define_regions_gene(gi, chr = 99), "No genes found")
})

test_that("define_regions_window generates correct windows", {
  windows <- define_regions_window(chr = 22, start = 1, end = 10000,
                                    window_size = 2000, step_size = 1000)
  expect_equal(names(windows), c("region_id", "chr", "start", "end", "label"))
  expect_true(all(windows$end - windows$start + 1 <= 2000))
  expect_equal(windows$start[2] - windows$start[1], 1000)
})

test_that("define_regions_window default step is 50% overlap", {
  windows <- define_regions_window(chr = 1, start = 1, end = 10000,
                                    window_size = 2000)
  expect_equal(windows$start[2] - windows$start[1], 1000)
})

test_that("define_regions_window handles region smaller than window", {
  windows <- define_regions_window(chr = 1, start = 1, end = 500,
                                    window_size = 2000)
  expect_equal(nrow(windows), 1)
  expect_equal(windows$end, 500)
})

test_that("define_regions_custom validates required columns", {
  expect_error(
    define_regions_custom(data.frame(chr = 1, start = 1, end = 2)),
    "Missing required columns"
  )
})

test_that("define_regions_custom works with label", {
  tbl <- data.frame(
    region_id = "r1", chr = "22", start = 100, end = 200,
    label = "my_region", stringsAsFactors = FALSE
  )
  result <- define_regions_custom(tbl)
  expect_equal(result$label, "my_region")
})

test_that("define_regions_custom defaults label to region_id", {
  tbl <- data.frame(
    region_id = "r1", chr = "22", start = 100, end = 200,
    stringsAsFactors = FALSE
  )
  result <- define_regions_custom(tbl)
  expect_equal(result$label, "r1")
})
