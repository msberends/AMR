context("mdro.R")


test_that("MDRO works", {
  library(dplyr)

  outcome <- MDRO(septic_patients, "EUCAST", info = FALSE)
  # check class
  expect_equal(outcome %>% class(), c('ordered', 'factor'))

  outcome <- MDRO(septic_patients, "nl", info = FALSE)
  # check class
  expect_equal(outcome %>% class(), c('ordered', 'factor'))

  # septic_patients should have these finding using Dutch guidelines
  expect_equal(outcome %>% freq(toConsole = FALSE) %>% pull(Count), c(3, 21))

  expect_equal(BRMO(septic_patients, info = FALSE), MDRO(septic_patients, "nl", info = FALSE))

})
