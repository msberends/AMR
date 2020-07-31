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

context("deprecated.R")

test_that("deprecated functions work", {
  skip_on_cran()
  expect_equal(suppressWarnings(portion_S(example_isolates$AMX)), proportion_S(example_isolates$AMX))
  expect_equal(suppressWarnings(portion_SI(example_isolates$AMX)), proportion_SI(example_isolates$AMX))
  expect_equal(suppressWarnings(portion_I(example_isolates$AMX)), proportion_I(example_isolates$AMX))
  expect_equal(suppressWarnings(portion_IR(example_isolates$AMX)), proportion_IR(example_isolates$AMX))
  expect_equal(suppressWarnings(portion_R(example_isolates$AMX)), proportion_R(example_isolates$AMX))
})
