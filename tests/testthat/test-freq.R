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

context("freq.R")

test_that("frequency table works", {
  library(clean)
  # mo
  expect_true(is.freq(freq(example_isolates$mo)))
  # for this to work, the output of mo_gramstain() is to be expected as follows:
  expect_equal(mo_gramstain("B_ESCHR_COL", language = NULL), "Gram-negative")
  expect_equal(mo_gramstain("B_STPHY_AUR", language = NULL), "Gram-positive")
  
  # rsi
  expect_true(is.freq(freq(example_isolates$AMX)))
  library(dplyr)
  expect_true(is.freq(example_isolates %>% freq(AMX)))
})

