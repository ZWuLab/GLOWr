# Test file for GLOWr package

test_that("Package loads successfully", {
  expect_true(require(GLOWr, quietly = TRUE))
})

test_that("Package has proper structure", {
  # Check that package namespace exists
  expect_true("GLOWr" %in% loadedNamespaces() || "GLOWr" %in% .packages())
})

# Placeholder test - will be expanded as functions are added
test_that("Placeholder for future tests", {
  expect_true(TRUE)
})
