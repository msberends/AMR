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
# Visit our website for more info: https://msberends.gitab.io/AMR.     #
# ==================================================================== #

context("mo_history.R")

test_that("mo_history works", {
  clean_mo_history(force = TRUE)
  expect_equal(read_mo_history(force = TRUE),
               NULL)

  expect_equal(as.character(suppressWarnings(as.mo("testsubject"))), "UNKNOWN")

  set_mo_history("testsubject", "B_ESCHR_COL",
                 uncertainty_level = translate_allow_uncertain(TRUE),
                 force = TRUE)

  expect_equal(get_mo_history("testsubject",
                              uncertainty_level = translate_allow_uncertain(TRUE),
                              force = TRUE),
               "B_ESCHR_COL")

  expect_equal(as.character(suppressWarnings(as.mo("testsubject"))), "B_ESCHR_COL")

  expect_equal(colnames(read_mo_history(force = TRUE)),
               c("x", "mo", "uncertainty_level", "package_version"))
})
