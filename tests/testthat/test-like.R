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
# how to conduct AMR data analysis: https://amr-for-r.org              #
# ==================================================================== #

test_that("test-like.R", {
  skip_on_cran()

  expect_true(sum("test" %like% c("^t", "^s")) == 1)

  expect_true("test" %like% "test")
  expect_false("test" %like_case% "TEST")
  expect_true(factor("test") %like% factor("t"))
  expect_true(factor("test") %like% "t")
  expect_true("test" %like% factor("t"))

  expect_true(as.factor("test") %like% "TEST")
  expect_identical(
    factor(c("Test case", "Something different", "Yet another thing")) %like% c("case", "diff", "yet"),
    c(TRUE, TRUE, TRUE)
  )
  expect_identical(
    "test" %like% c("t", "e", "s", "t"),
    c(TRUE, TRUE, TRUE, TRUE)
  )
  expect_identical(
    factor("test") %like% factor(c("t", "e", "s", "t")),
    c(TRUE, TRUE, TRUE, TRUE)
  )
})
