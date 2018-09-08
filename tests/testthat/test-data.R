context("data.R")

test_that("data sets are valid", {
  # IDs should always be unique
  expect_identical(nrow(antibiotics), length(unique(antibiotics$atc)))
  expect_identical(nrow(microorganisms), length(unique(microorganisms$mo)))
})
