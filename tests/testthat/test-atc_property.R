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

context("atc_property.R")

test_that("atc_property works", {
  expect_equal(atc_certe("amox"), "amox")
  expect_equal(atc_name("amox", language = "en"), "Amoxicillin")
  expect_equal(atc_name("amox", language = "nl"), "Amoxicilline")
  expect_equal(atc_official("amox", language = "en"), "Amoxicillin")
  expect_equal(atc_trivial_nl("amox"), "Amoxicilline")
  expect_equal(atc_umcg("amox"), "AMOX")
  expect_equal(class(atc_tradenames("amox")), "character")
  expect_equal(class(atc_tradenames(c("amox", "amox"))), "list")

  expect_error(atc_property("amox", "invalid property"))
  expect_error(atc_name("amox", language = "INVALID"))
  expect_output(print(atc_name("amox", language = NULL)))
})
