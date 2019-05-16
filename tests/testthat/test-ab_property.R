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

context("ab_property.R")

test_that("ab_property works", {

  expect_identical(ab_name("AMX"), "Amoxicillin")
  expect_identical(as.character(ab_atc("AMX")), "J01CA04")
  expect_identical(ab_cid("AMX"), as.integer(33613))

  expect_equal(class(ab_tradenames("AMX")), "character")
  expect_equal(class(ab_tradenames(c("AMX", "AMX"))), "list")

  expect_identical(ab_group("AMX"), "Beta-lactams/penicillins")
  expect_identical(ab_atc_group1("AMX"), "Beta-lactam antibacterials, penicillins")
  expect_identical(ab_atc_group2("AMX"), "Penicillins with extended spectrum")

  expect_identical(ab_name("Fluclox"), "Flucloxacillin")
  expect_identical(ab_name("fluklox"), "Flucloxacillin")
  expect_identical(ab_name("floxapen"), "Flucloxacillin")
  expect_identical(ab_name(21319) , "Flucloxacillin")
  expect_identical(ab_name("J01CF05"), "Flucloxacillin")

  expect_identical(ab_ddd("AMX", "oral"), 1)
  expect_identical(ab_ddd("AMX", "oral", units = TRUE) , "g")
  expect_identical(ab_ddd("AMX", "iv"), 1)
  expect_identical(ab_ddd("AMX", "iv", units = TRUE) , "g")

  expect_identical(ab_name(x = c("AMC", "PLB")), c("Amoxicillin/clavulanic acid", "Polymyxin B"))
  expect_identical(ab_name(x = c("AMC", "PLB"), tolower = TRUE),
                   c("amoxicillin/clavulanic acid", "polymyxin B"))

  expect_equal(class(ab_info("AMX")), "list")

  expect_error(ab_property("amox", "invalid property"))
  expect_error(ab_name("amox", language = "INVALID"))
  expect_output(print(ab_name("amox", language = NULL)))
})
