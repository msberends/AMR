context("kurtosis.R")

test_that("kurtosis works", {
  expect_equal(kurtosis(septic_patients$age),
               6.423118,
               tolerance = 0.00001)
  expect_equal(unname(kurtosis(data.frame(septic_patients$age))),
               6.423118,
               tolerance = 0.00001)
  expect_equal(kurtosis(matrix(septic_patients$age)),
               6.423118,
               tolerance = 0.00001)
})
