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

context("isolate_identifier.R")

test_that("isolate_identifier works", {
  x <- suppressMessages(isolate_identifier(example_isolates))
  expect_s3_class(x, "isolate_identifier")
  expect_s3_class(x, "character")
  
  expect_equal(suppressMessages(
    unique(nchar(isolate_identifier(example_isolates, cols_ab = carbapenems(), col_mo = FALSE)))),
    2)
  
  expect_warning(isolate_identifier(example_isolates[, 1:3, drop = FALSE])) # without mo and without rsi
  expect_warning(isolate_identifier(example_isolates[, 1:9, drop = FALSE])) # only without rsi
  
  
  expect_output(print(x))
  expect_s3_class(unique(c(x, x)), "isolate_identifier")
  
})
