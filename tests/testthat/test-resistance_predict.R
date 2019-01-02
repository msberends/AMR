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

context("portion.R")

test_that("prediction of rsi works", {
  amox_R <- septic_patients %>%
    filter(mo == "B_ESCHR_COL") %>%
    rsi_predict(col_ab = "amox",
                col_date = "date",
                minimum = 10,
                info = TRUE) %>%
    pull("value")
  # amox resistance will increase according to data set `septic_patients`
  expect_true(amox_R[3] < amox_R[20])

  library(dplyr)

  expect_output(rsi_predict(tbl = filter(septic_patients, mo == "B_ESCHR_COL"),
                            model = "binomial",
                            col_ab = "amox",
                            col_date = "date",
                            info = TRUE))
  expect_output(rsi_predict(tbl = filter(septic_patients, mo == "B_ESCHR_COL"),
                            model = "loglin",
                            col_ab = "amox",
                            col_date = "date",
                            info = TRUE))
  expect_output(rsi_predict(tbl = filter(septic_patients, mo == "B_ESCHR_COL"),
                            model = "lin",
                            col_ab = "amox",
                            col_date = "date",
                            info = TRUE))

  expect_error(rsi_predict(tbl = filter(septic_patients, mo == "B_ESCHR_COL"),
                           model = "INVALID MODEL",
                           col_ab = "amox",
                           col_date = "date",
                           info = TRUE))
  expect_error(rsi_predict(tbl = filter(septic_patients, mo == "B_ESCHR_COL"),
                           col_ab = "NOT EXISTING COLUMN",
                           col_date = "date",
                           info = TRUE))
  expect_error(rsi_predict(tbl = filter(septic_patients, mo == "B_ESCHR_COL"),
                           col_ab = "amox",
                           col_date = "NOT EXISTING COLUMN",
                           info = TRUE))
  # almost all E. coli are mero S in the Netherlands :)
  expect_error(resistance_predict(tbl = filter(septic_patients, mo == "B_ESCHR_COL"),
                                  col_ab = "mero",
                                  col_date = "date",
                                  info = TRUE))

  expect_error(portion_df(c("A", "B", "C")))
  expect_error(portion_df(septic_patients[,"date"]))
})
