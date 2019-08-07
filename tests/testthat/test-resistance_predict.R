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
# Visit our website for more info: https://msberends.gitlab.io/AMR.    #
# ==================================================================== #

context("portion.R")

test_that("prediction of rsi works", {
  AMX_R <- septic_patients %>%
    filter(mo == "B_ESCHR_COL") %>%
    rsi_predict(col_ab = "AMX",
                col_date = "date", 
                model = "binomial",
                minimum = 10,
                info = TRUE) %>%
    pull("value")
  # AMX resistance will increase according to data set `septic_patients`
  expect_true(AMX_R[3] < AMX_R[20])

  x <- resistance_predict(septic_patients, col_ab = "AMX", year_min = 2010, model = "binomial")
  plot(x)
  ggplot_rsi_predict(x)
  expect_error(ggplot_rsi_predict(septic_patients))

  library(dplyr)

  expect_output(rsi_predict(x = filter(septic_patients, mo == "B_ESCHR_COL"),
                            model = "binomial",
                            col_ab = "AMX",
                            col_date = "date",
                            info = TRUE))
  expect_output(rsi_predict(x = filter(septic_patients, mo == "B_ESCHR_COL"),
                            model = "loglin",
                            col_ab = "AMX",
                            col_date = "date",
                            info = TRUE))
  expect_output(rsi_predict(x = filter(septic_patients, mo == "B_ESCHR_COL"),
                            model = "lin",
                            col_ab = "AMX",
                            col_date = "date",
                            info = TRUE))

  expect_error(rsi_predict(x = filter(septic_patients, mo == "B_ESCHR_COL"),
                           model = "INVALID MODEL",
                           col_ab = "AMX",
                           col_date = "date",
                           info = TRUE))
  expect_error(rsi_predict(x = filter(septic_patients, mo == "B_ESCHR_COL"),
                           model = "binomial",
                           col_ab = "NOT EXISTING COLUMN",
                           col_date = "date",
                           info = TRUE))
  expect_error(rsi_predict(x = filter(septic_patients, mo == "B_ESCHR_COL"),
                           model = "binomial",
                           col_ab = "AMX",
                           col_date = "NOT EXISTING COLUMN",
                           info = TRUE))
  expect_error(rsi_predict(x = filter(septic_patients, mo == "B_ESCHR_COL"),
                           col_ab = "AMX",
                           col_date = "NOT EXISTING COLUMN",
                           info = TRUE))
  expect_error(rsi_predict(x = filter(septic_patients, mo == "B_ESCHR_COL"),
                           col_ab = "AMX",
                           col_date = "date",
                           info = TRUE))
  # almost all E. coli are MEM S in the Netherlands :)
  expect_error(resistance_predict(x = filter(septic_patients, mo == "B_ESCHR_COL"),
                                  model = "binomial",
                                  col_ab = "MEM",
                                  col_date = "date",
                                  info = TRUE))
})
