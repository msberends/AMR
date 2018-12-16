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

context("abname.R")

test_that("abname works", {
  expect_equal(abname("AMOX"), "Amoxicillin")
  expect_equal(abname(c("AMOX", "GENT")), c("Amoxicillin", "Gentamicin"))
  expect_equal(abname(c("AMOX+GENT")), "Amoxicillin + gentamicin")
  expect_equal(abname("AMOX", from = 'umcg'), "Amoxicillin")
  expect_equal(abname("amox", from = 'certe', tolower = TRUE), "amoxicillin")
  expect_equal(abname("J01CA04", from = 'atc'), "Amoxicillin")
  expect_equal(abname(c("amox", "J01CA04", "Trimox", "dispermox", "Amoxil")),
                      rep("Amoxicillin", 5))
  expect_equal(abname("AMOX", to = 'atc'), "J01CA04")

  expect_error(abname("AMOX", to = c(1:3)))
  expect_error(abname("AMOX", to = "test"))
  expect_warning(abname("NOTEXISTING
       "))
  expect_warning(abname("AMOX or GENT"))

  # this one is being found with as.atc internally
  expect_equal(abname("flu_clox123"), "Flucloxacillin")
})
