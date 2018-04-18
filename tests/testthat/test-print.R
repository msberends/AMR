context("print.R")


test_that("tibble printing works", {
  library(dplyr)
  expect_output(print(starwars))
  expect_output(print(septic_patients))
})
