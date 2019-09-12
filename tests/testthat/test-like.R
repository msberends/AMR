# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Analysis                              #
#                                                                      #
# SOURCE                                                               #
# https://gitlab.com/msberends/AMR                                     #
#                                                                      #
# LICENCE                                                              #
# (c) 2019 Berends MS (m.s.berends@umcg.nl), Luz CF (c.f.luz@umcg.nl)  #
#                                                                      #
# This R package is free software; you can freely use and distribute   #
# it for both personal and commercial purposes under the terms of the  #
# GNU General Public License version 2.0 (GNU GPL-2), as published by  #
# the Free Software Foundation.                                        #
#                                                                      #
# This R package was created for academic research and was publicly    #
# released in the hope that it will be useful, but it comes WITHOUT    #
# ANY WARRANTY OR LIABILITY.                                           #
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
