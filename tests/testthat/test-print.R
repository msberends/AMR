context("print.R")


test_that("tibble printing works", {
  library(dplyr)
  library(data.table)
  expect_output(print(starwars))
  expect_output(print(starwars %>% group_by(homeworld, gender)))
  expect_output(print(starwars %>% as.data.table(), print.keys = TRUE))
  expect_output(print(septic_patients))
})
