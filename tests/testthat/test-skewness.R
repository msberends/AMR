context("skewness.R")

test_that("skewness works", {
  expect_equal(skewness(septic_patients$age),
               -0.8958019,
               tolerance = 0.00001)
  expect_equal(unname(skewness(data.frame(septic_patients$age))),
               -0.8958019,
               tolerance = 0.00001)
  expect_equal(skewness(matrix(septic_patients$age)),
               -0.8958019,
               tolerance = 0.00001)
})
