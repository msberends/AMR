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

context("deprecated.R")

test_that("deprecated functions work", {

  # first 5 chars of official name
  expect_equal(suppressWarnings(as.character(as.atc(c("nitro", "cipro")))),
               c("J01XE01", "J01MA02"))

  # EARS-Net
  expect_equal(suppressWarnings(as.character(as.atc("AMX"))),
               "J01CA04")

  expect_equal(suppressWarnings(guess_ab_col(data.frame(AMP_ND10 = "R",
                                                        AMC_ED20 = "S"),
                                             as.atc("augmentin"))),
               "AMC_ED20")
})
