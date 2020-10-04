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

context("kurtosis.R")

test_that("kurtosis works", {
  skip_on_cran()
  expect_equal(kurtosis(example_isolates$age),
               3.549319,
               tolerance = 0.00001)
  
  expect_equal(unname(kurtosis(data.frame(example_isolates$age))),
               3.549319,
               tolerance = 0.00001)
  expect_equal(unname(kurtosis(data.frame(example_isolates$age), excess = TRUE)),
               0.549319,
               tolerance = 0.00001)
  
  expect_equal(kurtosis(matrix(example_isolates$age)),
               3.549319,
               tolerance = 0.00001)
  expect_equal(kurtosis(matrix(example_isolates$age), excess = TRUE),
               0.549319,
               tolerance = 0.00001)
})
