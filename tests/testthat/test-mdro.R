context("mdro.R")


test_that("MDRO works", {
  library(dplyr)

  expect_error(suppressWarnings(MDRO(septic_patients, "invalid", col_bactid = "mo", info = TRUE)))
  expect_error(suppressWarnings(MDRO(septic_patients, "fr", col_bactid = "mo", info = TRUE)))
  expect_error(suppressWarnings(MDRO(septic_patients, country = c("de", "nl"), info = TRUE)))
  expect_error(suppressWarnings(MDRO(septic_patients, col_mo = "invalid", info = TRUE)))

  outcome <- suppressWarnings(MDRO(septic_patients))
  outcome <- suppressWarnings(EUCAST_exceptional_phenotypes(septic_patients, info = TRUE))
  # check class
  expect_equal(outcome %>% class(), c('ordered', 'factor'))

  outcome <- suppressWarnings(MDRO(septic_patients, "nl", info = TRUE))
  # check class
  expect_equal(outcome %>% class(), c('ordered', 'factor'))

  # septic_patients should have these finding using Dutch guidelines
  expect_equal(outcome %>% freq() %>% pull(count),
               c(2, 14)) # 2 unconfirmed, 14 positive

  expect_equal(BRMO(septic_patients, info = FALSE), MDRO(septic_patients, "nl", info = FALSE))

  # still working on German guidelines
  expect_error(suppressWarnings(MRGN(septic_patients, info = TRUE)))

})
