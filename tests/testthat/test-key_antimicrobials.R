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
# how to conduct AMR data analysis: https://msberends.github.io/AMR/   #
# ==================================================================== #

test_that("key_antimicrobials works", {
  expect_equal(length(key_antimicrobials(example_isolates, antifungal = NULL)), nrow(example_isolates))
  expect_false(all(is.na(key_antimicrobials(example_isolates, antifungal = NULL))))
  expect_true(antimicrobials_equal("SSS", "SSS", type = "points"))
  expect_false(antimicrobials_equal("SSS", "SRS", type = "keyantimicrobials"))
  expect_true(antimicrobials_equal("SSS", "SRS", type = "points"))
  expect_true(antimicrobials_equal("SSS", "SIS", ignore_I = TRUE, type = "keyantimicrobials"))
  expect_false(antimicrobials_equal("SSS", "SIS", ignore_I = FALSE, type = "keyantimicrobials"))
  expect_true(antimicrobials_equal(".SS", "SI.", ignore_I = TRUE, type = "keyantimicrobials"))
  expect_false(antimicrobials_equal(".SS", "SI.", ignore_I = FALSE, type = "keyantimicrobials"))

  # expect_warning(key_antimicrobials(example_isolates[rep(1, 10), , drop = FALSE]))
})
