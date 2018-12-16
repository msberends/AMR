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

context("key_antibiotics.R")

test_that("keyantibiotics work", {
  expect_equal(length(key_antibiotics(septic_patients, warnings = FALSE)), nrow(septic_patients))
  expect_false(all(is.na(key_antibiotics(septic_patients))))
  expect_true(key_antibiotics_equal("SSS", "SSS"))
  expect_false(key_antibiotics_equal("SSS", "SRS"))
  expect_true(key_antibiotics_equal("SSS", "SIS", ignore_I = TRUE))
  expect_false(key_antibiotics_equal("SSS", "SIS", ignore_I = FALSE))
  expect_true(key_antibiotics_equal(".SS", "SI.", ignore_I = TRUE))
  expect_false(key_antibiotics_equal(".SS", "SI.", ignore_I = FALSE))
})
