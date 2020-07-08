# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Analysis                              #
#                                                                      #
# SOURCE                                                               #
# https://github.com/msberends/AMR                                     #
#                                                                      #
# LICENCE                                                              #
# (c) 2018-2020 Berends MS, Luz CF et al.                              #
#                                                                      #
# This R package is free software; you can freely use and distribute   #
# it for both personal and commercial purposes under the terms of the  #
# GNU General Public License version 2.0 (GNU GPL-2), as published by  #
# the Free Software Foundation.                                        #
#                                                                      #
# We created this package for both routine data analysis and academic  #
# research and it was publicly released in the hope that it will be    #
# useful, but it comes WITHOUT ANY WARRANTY OR LIABILITY.              #
# Visit our website for more info: https://msberends.github.io/AMR.    #
# ==================================================================== #

context("bug_drug_combinations.R")

test_that("bug_drug_combinations works", {
  
  skip_on_cran()
  
  b <- suppressWarnings(bug_drug_combinations(example_isolates))
  expect_s3_class(b, "bug_drug_combinations")
  expect_output(print(b))
  expect_true(is.data.frame(format(b)))
  expect_true(is.data.frame(format(b, combine_IR = TRUE, add_ab_group = FALSE)))
})
