# ==================================================================== #
# TITLE:                                                               #
# AMR: An R Package for Working with Antimicrobial Resistance Data     #
#                                                                      #
# SOURCE CODE:                                                         #
# https://github.com/msberends/AMR                                     #
#                                                                      #
# PLEASE CITE THIS SOFTWARE AS:                                        #
# Berends MS, Luz CF, Friedrich AW, et al. (2022).                     #
# AMR: An R Package for Working with Antimicrobial Resistance Data.    #
# Journal of Statistical Software, 104(3), 1-31.                       #
# https://doi.org/10.18637/jss.v104.i03                                #
#                                                                      #
# Developed at the University of Groningen and the University Medical  #
# Center Groningen in The Netherlands, in collaboration with many      #
# colleagues from around the world, see our website.                   #
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
# how to conduct AMR data analysis: https://amr-for-r.org/             #
# ==================================================================== #

test_that("test-join_microorganisms.R", {
  skip_on_cran()

  unjoined <- example_isolates
  inner <- inner_join_microorganisms(example_isolates)
  left <- left_join_microorganisms(example_isolates)
  semi <- semi_join_microorganisms(example_isolates)
  anti <- anti_join_microorganisms(example_isolates)
  suppressWarnings(right <- right_join_microorganisms(example_isolates))
  suppressWarnings(full <- full_join_microorganisms(example_isolates))

  expect_true(ncol(unjoined) < ncol(inner))
  expect_true(nrow(unjoined) == nrow(inner))

  expect_true(ncol(unjoined) < ncol(left))
  expect_true(nrow(unjoined) == nrow(left))

  expect_true(ncol(semi) == ncol(semi))
  expect_true(nrow(semi) == nrow(semi))

  expect_true(nrow(anti) == 0)

  expect_true(nrow(unjoined) < nrow(right))
  expect_true(nrow(unjoined) < nrow(full))

  expect_equal(nrow(inner_join_microorganisms("B_ESCHR_COLI")), 1)
  expect_equal(nrow(inner_join_microorganisms("B_ESCHR_COLI", by = c("mo" = "mo"))), 1)

  expect_equal(nrow(left_join_microorganisms("B_ESCHR_COLI")), 1)

  expect_equal(nrow(semi_join_microorganisms("B_ESCHR_COLI")), 1)
  expect_equal(nrow(anti_join_microorganisms("B_ESCHR_COLI")), 0)

  # expect_warning(right_join_microorganisms("B_ESCHR_COLI"))
  # expect_warning(full_join_microorganisms("B_ESCHR_COLI"))
})
