# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Analysis                              #
#                                                                      #
# SOURCE                                                               #
# https://gitlab.com/msberends/AMR                                     #
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
# Visit our website for more info: https://msberends.gitlab.io/AMR.    #
# ==================================================================== #

context("key_antibiotics.R")

test_that("keyantibiotics work", {
  expect_equal(length(key_antibiotics(example_isolates, warnings = FALSE)), nrow(example_isolates))
  expect_false(all(is.na(key_antibiotics(example_isolates))))
  expect_true(key_antibiotics_equal("SSS", "SSS"))
  expect_false(key_antibiotics_equal("SSS", "SRS"))
  expect_true(key_antibiotics_equal("SSS", "SIS", ignore_I = TRUE))
  expect_false(key_antibiotics_equal("SSS", "SIS", ignore_I = FALSE))
  expect_true(key_antibiotics_equal(".SS", "SI.", ignore_I = TRUE))
  expect_false(key_antibiotics_equal(".SS", "SI.", ignore_I = FALSE))
  
  library(dplyr)
  expect_warning(key_antibiotics(example_isolates %>% slice(rep(1, 10))))
})
