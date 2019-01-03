# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Analysis                              #
#                                                                      #
# SOURCE                                                               #
# https://gitlab.com/msberends/AMR                                     #
#                                                                      #
# LICENCE                                                              #
# (c) 2019 Berends MS (m.s.berends@umcg.nl), Luz CF (c.f.luz@umcg.nl)  #
#                                                                      #
# This R package is free software; you can freely use and distribute   #
# it for both personal and commercial purposes under the terms of the  #
# GNU General Public License version 2.0 (GNU GPL-2), as published by  #
# the Free Software Foundation.                                        #
#                                                                      #
# This R package was created for academic research and was publicly    #
# released in the hope that it will be useful, but it comes WITHOUT    #
# ANY WARRANTY OR LIABILITY.                                           #
# Visit our website for more info: https://msberends.gitab.io/AMR.     #
# ==================================================================== #

context("mdro.R")

test_that("mdro works", {
  library(dplyr)

  expect_error(suppressWarnings(mdro(septic_patients, "invalid", col_mo = "mo", info = TRUE)))
  expect_error(suppressWarnings(mdro(septic_patients, "fr", info = TRUE)))
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

  expect_equal(
    suppressWarnings(
      brmo(septic_patients, info = FALSE)),
    suppressWarnings(
      mdro(septic_patients, "nl", info = FALSE)
    )
  )

  # still working on German guidelines
  expect_error(suppressWarnings(mrgn(septic_patients, info = TRUE)))

})
