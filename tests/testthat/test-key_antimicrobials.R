# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Data Analysis for R                   #
#                                                                      #
# SOURCE                                                               #
# https://github.com/msberends/AMR                                     #
#                                                                      #
# LICENCE                                                              #
# (c) 2018-2021 Berends MS, Luz CF et al.                              #
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
# how to conduct AMR data analysis: https://msberends.github.io/AMR/   #
# ==================================================================== #

context("key_antimcrobials.R")

test_that("key_antimcrobials work", {
  skip_on_cran()
  expect_equal(length(key_antimicrobials(example_isolates, antifungal = NULL)), nrow(example_isolates))
  expect_false(all(is.na(key_antimicrobials(example_isolates, antifungal = NULL))))
  expect_true(antimicrobials_equal("SSS", "SSS", type = "points"))
  expect_false(antimicrobials_equal("SSS", "SRS", type = "keyantimicrobials"))
  expect_true(antimicrobials_equal("SSS", "SRS", type = "points"))
  expect_true(antimicrobials_equal("SSS", "SIS", ignore_I = TRUE, type = "keyantimicrobials"))
  expect_false(antimicrobials_equal("SSS", "SIS", ignore_I = FALSE, type = "keyantimicrobials"))
  expect_true(antimicrobials_equal(".SS", "SI.", ignore_I = TRUE, type = "keyantimicrobials"))
  expect_false(antimicrobials_equal(".SS", "SI.", ignore_I = FALSE, type = "keyantimicrobials"))
  
  library(dplyr, warn.conflicts = FALSE)
  expect_warning(key_antimicrobials(example_isolates %>% slice(rep(1, 10))))
})
