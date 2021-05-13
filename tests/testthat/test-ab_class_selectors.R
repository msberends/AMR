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

context("ab_class_selectors.R")

test_that("Antibiotic class selectors work", {
  skip_on_cran()
  
  if (suppressWarnings(require("dplyr"))) {
    expect_lt(example_isolates %>% select(aminoglycosides()) %>% ncol(), ncol(example_isolates))
    expect_lt(example_isolates %>% select(carbapenems()) %>% ncol(), ncol(example_isolates))
    expect_lt(example_isolates %>% select(cephalosporins()) %>% ncol(), ncol(example_isolates))
    expect_lt(example_isolates %>% select(cephalosporins_1st()) %>% ncol(), ncol(example_isolates))
    expect_lt(example_isolates %>% select(cephalosporins_2nd()) %>% ncol(), ncol(example_isolates))
    expect_lt(example_isolates %>% select(cephalosporins_3rd()) %>% ncol(), ncol(example_isolates))
    expect_lt(example_isolates %>% select(cephalosporins_4th()) %>% ncol(), ncol(example_isolates))
    expect_lt(example_isolates %>% select(cephalosporins_5th()) %>% ncol(), ncol(example_isolates))
    expect_lt(example_isolates %>% select(fluoroquinolones()) %>% ncol(), ncol(example_isolates))
    expect_lt(example_isolates %>% select(glycopeptides()) %>% ncol(), ncol(example_isolates))
    expect_lt(example_isolates %>% select(macrolides()) %>% ncol(), ncol(example_isolates))
    expect_lt(example_isolates %>% select(oxazolidinones()) %>% ncol(), ncol(example_isolates))
    expect_lt(example_isolates %>% select(penicillins()) %>% ncol(), ncol(example_isolates))
    expect_lt(example_isolates %>% select(tetracyclines()) %>% ncol(), ncol(example_isolates))
  }
})
