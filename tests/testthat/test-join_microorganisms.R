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

context("join_microorganisms.R")

test_that("joins work", {
  unjoined <- septic_patients
  inner <- septic_patients %>% inner_join_microorganisms()
  left <- septic_patients %>% left_join_microorganisms()
  semi <- septic_patients %>% semi_join_microorganisms()
  anti <- septic_patients %>% anti_join_microorganisms()
  suppressWarnings(right <- septic_patients %>% right_join_microorganisms())
  suppressWarnings(full <- septic_patients %>% full_join_microorganisms())

  expect_true(ncol(unjoined) < ncol(inner))
  expect_true(nrow(unjoined) == nrow(inner))

  expect_true(ncol(unjoined) < ncol(left))
  expect_true(nrow(unjoined) == nrow(left))

  expect_true(ncol(semi) == ncol(semi))
  expect_true(nrow(semi) == nrow(semi))

  expect_true(nrow(anti) == 0)

  expect_true(nrow(unjoined) < nrow(right))
  expect_true(nrow(unjoined) < nrow(full))


  expect_equal(nrow(inner_join_microorganisms("B_ESCHR_COL")), 1)
  expect_equal(nrow(inner_join_microorganisms("B_ESCHR_COL", by = c("mo" = "mo"))), 1)
  expect_warning(inner_join_microorganisms("Escherichia", by = c("mo" = "genus")))

  expect_equal(nrow(left_join_microorganisms("B_ESCHR_COL")), 1)
  expect_warning(left_join_microorganisms("Escherichia", by = c("mo" = "genus")))

  expect_equal(nrow(semi_join_microorganisms("B_ESCHR_COL")), 1)
  expect_equal(nrow(anti_join_microorganisms("B_ESCHR_COL")), 0)

  expect_warning(right_join_microorganisms("B_ESCHR_COL"))
  expect_warning(full_join_microorganisms("B_ESCHR_COL"))

})
