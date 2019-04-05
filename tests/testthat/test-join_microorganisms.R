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
