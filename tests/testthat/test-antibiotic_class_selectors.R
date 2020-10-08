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

context("ab_class_selectors.R")

test_that("Antibiotic class selectors work", {
  skip_on_cran()
  
  expect_lt(example_isolates %>% dplyr::select(aminoglycosides()) %>% ncol(), ncol(example_isolates))
  expect_lt(example_isolates %>% dplyr::select(carbapenems()) %>% ncol(), ncol(example_isolates))
  expect_lt(example_isolates %>% dplyr::select(cephalosporins()) %>% ncol(), ncol(example_isolates))
  expect_lt(example_isolates %>% dplyr::select(cephalosporins_1st()) %>% ncol(), ncol(example_isolates))
  expect_lt(example_isolates %>% dplyr::select(cephalosporins_2nd()) %>% ncol(), ncol(example_isolates))
  expect_lt(example_isolates %>% dplyr::select(cephalosporins_3rd()) %>% ncol(), ncol(example_isolates))
  expect_lt(example_isolates %>% dplyr::select(cephalosporins_4th()) %>% ncol(), ncol(example_isolates))
  expect_lt(example_isolates %>% dplyr::select(cephalosporins_5th()) %>% ncol(), ncol(example_isolates))
  expect_lt(example_isolates %>% dplyr::select(fluoroquinolones()) %>% ncol(), ncol(example_isolates))
  expect_lt(example_isolates %>% dplyr::select(glycopeptides()) %>% ncol(), ncol(example_isolates))
  expect_lt(example_isolates %>% dplyr::select(macrolides()) %>% ncol(), ncol(example_isolates))
  expect_lt(example_isolates %>% dplyr::select(penicillins()) %>% ncol(), ncol(example_isolates))
  expect_lt(example_isolates %>% dplyr::select(tetracyclines()) %>% ncol(), ncol(example_isolates))
  
})
