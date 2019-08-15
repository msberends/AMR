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

context("misc.R")

test_that("percentages works", {
  expect_equal(percent(0.25), "25%")
  expect_equal(percent(0.5), "50%")
  expect_equal(percent(0.500, force_zero = TRUE), "50.0%")
  expect_equal(percent(0.1234), "12.3%")
  # round up 0.5
  expect_equal(percent(0.0054), "0.5%")
  expect_equal(percent(0.0055), "0.6%")
})

test_that("functions missing in older R versions work", {
  expect_equal(strrep("A", 5), "AAAAA")
  expect_equal(strrep(c("A", "B"), c(5, 2)), c("AAAAA", "BB"))
  expect_equal(trimws(" test "), "test")
  expect_equal(trimws(" test ", "l"), "test ")
  expect_equal(trimws(" test ", "r"), " test")
})

test_that("looking up ab columns works", {
  expect_warning(generate_warning_abs_missing(c("AMP", "AMX")))
  expect_warning(generate_warning_abs_missing(c("AMP", "AMX"), any = TRUE))
  expect_warning(get_column_abx(septic_patients, hard_dependencies = "FUS"))
  expect_message(get_column_abx(septic_patients, soft_dependencies = "FUS"))
  expect_message(get_column_abx(dplyr::rename(septic_patients, thisone = AMX), amox = "thisone", tmp = "thisone", verbose = TRUE))
  expect_warning(get_column_abx(dplyr::rename(septic_patients, thisone = AMX), amox = "thisone", tmp = "thisone", verbose = FALSE))
})
