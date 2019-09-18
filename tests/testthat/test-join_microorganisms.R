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
  unjoined <- example_isolates
  inner <- example_isolates %>% inner_join_microorganisms()
  left <- example_isolates %>% left_join_microorganisms()
  semi <- example_isolates %>% semi_join_microorganisms()
  anti <- example_isolates %>% anti_join_microorganisms()
  suppressWarnings(right <- example_isolates %>% right_join_microorganisms())
  suppressWarnings(full <- example_isolates %>% full_join_microorganisms())

  expect_true(ncol(unjoined) < ncol(inner))
  expect_true(nrow(unjoined) == nrow(inner))

  expect_true(ncol(unjoined) < ncol(left))
  expect_true(nrow(unjoined) == nrow(left))

  expect_true(ncol(semi) == ncol(semi))
  expect_true(nrow(semi) == nrow(semi))

  expect_true(nrow(anti) == 0)

  expect_true(nrow(unjoined) < nrow(right))
  expect_true(nrow(unjoined) < nrow(full))


  expect_equal(nrow(inner_join_microorganisms("B_ESCHR_COLI")), 1)
  expect_equal(nrow(inner_join_microorganisms("B_ESCHR_COLI", by = c("mo" = "mo"))), 1)
  expect_warning(inner_join_microorganisms("Escherichia", by = c("mo" = "genus")))

  expect_equal(nrow(left_join_microorganisms("B_ESCHR_COLI")), 1)
  expect_warning(left_join_microorganisms("Escherichia", by = c("mo" = "genus")))

  expect_equal(nrow(semi_join_microorganisms("B_ESCHR_COLI")), 1)
  expect_equal(nrow(anti_join_microorganisms("B_ESCHR_COLI")), 0)

  expect_warning(right_join_microorganisms("B_ESCHR_COLI"))
  expect_warning(full_join_microorganisms("B_ESCHR_COLI"))

})
