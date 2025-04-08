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

test_that("test-top_n_microorganisms.R", {
  skip_on_cran()

  out1 <- top_n_microorganisms(example_isolates, n = 3)
  out2 <- top_n_microorganisms(example_isolates, n = 5, property = "genus")
  out3 <- top_n_microorganisms(example_isolates, n = 5, property = "genus", n_for_each = 3)

  expect_equal(NROW(out1), 1015, tolerance = 0.5)
  expect_equal(NROW(out2), 1742, tolerance = 0.5)
  expect_equal(NROW(out3), 1497, tolerance = 0.5)

  expect_equal(length(table(out1$mo)), 3, tolerance = 0.5)
  expect_equal(length(table(out2$mo)), 39, tolerance = 0.5)
  expect_equal(length(table(out3$mo)), 13, tolerance = 0.5)

  expect_equal(length(unique(mo_genus(out2$mo))), 5, tolerance = 0.5)
  expect_equal(length(unique(mo_genus(out3$mo))), 5, tolerance = 0.5)
})
