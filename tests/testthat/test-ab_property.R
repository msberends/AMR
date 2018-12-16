# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Analysis                              #
#                                                                      #
# AUTHORS                                                              #
# Berends MS (m.s.berends@umcg.nl), Luz CF (c.f.luz@umcg.nl)           #
#                                                                      #
# LICENCE                                                              #
# This package is free software; you can redistribute it and/or modify #
# it under the terms of the GNU General Public License version 2.0,    #
# as published by the Free Software Foundation.                        #
#                                                                      #
# This R package is distributed in the hope that it will be useful,    #
# but WITHOUT ANY WARRANTY; without even the implied warranty of       #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        #
# GNU General Public License version 2.0 for more details.             #
# ==================================================================== #

context("ab_property.R")

test_that("ab_property works", {
  expect_equal(ab_certe("amox"), "amox")
  expect_equal(ab_name("amox", language = "en"), "Amoxicillin")
  expect_equal(ab_name("amox", language = "nl"), "Amoxicilline")
  expect_equal(ab_official("amox", language = "en"), "Amoxicillin")
  expect_equal(ab_trivial_nl("amox"), "Amoxicilline")
  expect_equal(ab_umcg("amox"), "AMOX")
  expect_equal(class(ab_tradenames("amox")), "character")
  expect_equal(class(ab_tradenames(c("amox", "amox"))), "list")
  expect_equal(ab_atc("amox"), as.character(as.atc("amox")))

  expect_error(ab_property("amox", "invalid property"))
  expect_error(ab_name("amox", language = "INVALID"))
  expect_output(print(ab_name("amox", language = NULL)))
})
