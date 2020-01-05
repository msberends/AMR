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

context("like.R")

test_that("`like` works", {
  expect_true(sum("test" %like% c("^t", "^s")) == 1)
  expect_true("test" %like% "test")
  expect_true("test" %like% "TEST")
  expect_true(as.factor("test") %like% "TEST")
  expect_identical(factor(c("Test case", "Something different", "Yet another thing")) %like% c("case", "diff", "yet"),
                   c(TRUE, TRUE, TRUE))
})
