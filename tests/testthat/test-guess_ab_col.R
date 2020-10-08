# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Analysis for R                        #
#                                                                      #
# SOURCE                                                               #
# https://github.com/msberends/AMR                                     #
#                                                                      #
# LICENCE                                                              #
# (c) 2018-2020 Berends MS, Luz CF et al.                              #
# Developed at the University of Groningen, the Netherlands, in        #
# collaboration with non-profit organisations Certe Medical            #
# Diagnostics & Advice, and University Medical Center Groningen.       # 
#                                                                      #
# This R package is free software; you can freely use and distribute   #
# it for both personal and commercial purposes under the terms of the  #
# GNU General Public License version 2.0 (GNU GPL-2), as published by  #
# the Free Software Foundation.                                        #
# We created this package for both routine data analysis and academic  #
# research and it was publicly released in the hope that it will be    #
# useful, but it comes WITHOUT ANY WARRANTY OR LIABILITY.              #
#                                                                      #
# Visit our website for the full manual and a complete tutorial about  #
# how to conduct AMR analysis: https://msberends.github.io/AMR/        #
# ==================================================================== #

context("guess_ab_col.R")

test_that("guess_ab_col works", {
  skip_on_cran()

  expect_equal(guess_ab_col(example_isolates, "amox"),
               "AMX")
  expect_equal(guess_ab_col(example_isolates, "amoxicillin"),
               "AMX")
  expect_equal(guess_ab_col(example_isolates, "J01AA07"),
               "TCY")
  expect_equal(guess_ab_col(example_isolates, "tetracycline"),
               "TCY")
  expect_equal(guess_ab_col(example_isolates, "TETR"),
               "TCY")

  df <- data.frame(AMP_ND10 = "R",
                   AMC_ED20 = "S")
  expect_equal(guess_ab_col(df, "ampicillin"),
               "AMP_ND10")
  expect_equal(guess_ab_col(df, "J01CR02"),
               "AMC_ED20")

})
