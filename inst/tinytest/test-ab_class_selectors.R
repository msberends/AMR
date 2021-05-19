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

if (!AMR:::current_R_older_than(3.2)) {
  # antibiotic class selectors require at least R-3.2
  expect_true(ncol(example_isolates[, aminoglycosides(), drop = FALSE]) < ncol(example_isolates))
  expect_true(ncol(example_isolates[, betalactams(), drop = FALSE]) < ncol(example_isolates))
  expect_true(ncol(example_isolates[, carbapenems(), drop = FALSE]) < ncol(example_isolates))
  expect_true(ncol(example_isolates[, cephalosporins(), drop = FALSE]) < ncol(example_isolates))
  expect_true(ncol(example_isolates[, cephalosporins_1st(), drop = FALSE]) < ncol(example_isolates))
  expect_true(ncol(example_isolates[, cephalosporins_2nd(), drop = FALSE]) < ncol(example_isolates))
  expect_true(ncol(example_isolates[, cephalosporins_3rd(), drop = FALSE]) < ncol(example_isolates))
  expect_true(ncol(example_isolates[, cephalosporins_4th(), drop = FALSE]) < ncol(example_isolates))
  expect_true(ncol(example_isolates[, cephalosporins_5th(), drop = FALSE]) < ncol(example_isolates))
  expect_true(ncol(example_isolates[, fluoroquinolones(), drop = FALSE]) < ncol(example_isolates))
  expect_true(ncol(example_isolates[, glycopeptides(), drop = FALSE]) < ncol(example_isolates))
  expect_true(ncol(example_isolates[, macrolides(), drop = FALSE]) < ncol(example_isolates))
  expect_true(ncol(example_isolates[, oxazolidinones(), drop = FALSE]) < ncol(example_isolates))
  expect_true(ncol(example_isolates[, penicillins(), drop = FALSE]) < ncol(example_isolates))
  expect_true(ncol(example_isolates[, tetracyclines(), drop = FALSE]) < ncol(example_isolates))
  
  # Examples:
  
  # select columns 'mo', 'AMK', 'GEN', 'KAN' and 'TOB'
  expect_equal(ncol(example_isolates[, c("mo", aminoglycosides())]), 5, tolerance = 0.5)
  
  # filter using any() or all()
  expect_equal(nrow(example_isolates[any(carbapenems() == "R"), ]), 55, tolerance = 0.5)
  expect_equal(nrow(subset(example_isolates, any(carbapenems() == "R"))), 55, tolerance = 0.5)
  
  # filter on any or all results in the carbapenem columns (i.e., IPM, MEM):
  expect_equal(nrow(example_isolates[any(carbapenems()), ]), 962, tolerance = 0.5)
  expect_equal(nrow(example_isolates[all(carbapenems()), ]), 756, tolerance = 0.5)
  
  # filter with multiple antibiotic selectors using c()
  expect_equal(nrow(example_isolates[all(c(carbapenems(), aminoglycosides()) == "R"), ]), 26, tolerance = 0.5)
  
  # filter + select in one go: get penicillins in carbapenems-resistant strains
  expect_equal(nrow(example_isolates[any(carbapenems() == "R"), penicillins()]), 55, tolerance = 0.5)
  expect_equal(ncol(example_isolates[any(carbapenems() == "R"), penicillins()]), 7, tolerance = 0.5)
}
