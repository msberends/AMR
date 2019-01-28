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

context("guess_ab_col.R")

test_that("guess_ab_col works", {

  expect_equal(guess_ab_col(septic_patients, "amox"),
               "amox")
  expect_equal(guess_ab_col(septic_patients, "amoxicillin"),
               "amox")
  expect_equal(guess_ab_col(septic_patients, "J01AA07"),
               "tetr")
  expect_equal(guess_ab_col(septic_patients, "tetracycline"),
               "tetr")
  expect_equal(guess_ab_col(septic_patients, "TETR"),
               "tetr")

  df <- data.frame(AMP_ND10 = "R",
                   AMC_ED20 = "S")
  expect_equal(guess_ab_col(df, "ampicillin"),
               "AMP_ND10")
  expect_equal(guess_ab_col(df, "J01CR02"),
               "AMC_ED20")
  expect_equal(guess_ab_col(df, as.atc("augmentin")),
               "AMC_ED20")
})
