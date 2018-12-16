# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Analysis                              #
#                                                                      #
# AUTHORS                                                              #
# Berends MS (m.s.berends@umcg.nl), Luz CF (c.f.luz@umcg.nl)           #
#                                                                      #
# LICENCE                                                              #
# This package is free software; you can redistribute it and/or modify #
# it under the terms of the GNU General Public License version 2.0,    #
# as published by the Free Software Foundation.                        #
#                                                                      #
# This R package is distributed in the hope that it will be useful,    #
# but WITHOUT ANY WARRANTY; without even the implied warranty of       #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        #
# GNU General Public License version 2.0 for more details.             #
# ==================================================================== #

context("mdro.R")

test_that("mdro works", {
  library(dplyr)

  expect_error(suppressWarnings(mdro(septic_patients, "invalid", col_bactid = "mo", info = TRUE)))
  expect_error(suppressWarnings(mdro(septic_patients, "fr", col_bactid = "mo", info = TRUE)))
  expect_error(suppressWarnings(mdro(septic_patients, country = c("de", "nl"), info = TRUE)))
  expect_error(suppressWarnings(mdro(septic_patients, col_mo = "invalid", info = TRUE)))

  outcome <- suppressWarnings(mdro(septic_patients))
  outcome <- suppressWarnings(eucast_exceptional_phenotypes(septic_patients, info = TRUE))
  # check class
  expect_equal(outcome %>% class(), c('ordered', 'factor'))

  outcome <- suppressWarnings(mdro(septic_patients, "nl", info = TRUE))
  # check class
  expect_equal(outcome %>% class(), c('ordered', 'factor'))

  # septic_patients should have these finding using Dutch guidelines
  expect_equal(outcome %>% freq() %>% pull(count),
               c(1989, 9, 2)) # 1989 neg, 9 pos, 2 unconfirmed

  expect_equal(brmo(septic_patients, info = FALSE), mdro(septic_patients, "nl", info = FALSE))

  # still working on German guidelines
  expect_error(suppressWarnings(mrgn(septic_patients, info = TRUE)))

})
