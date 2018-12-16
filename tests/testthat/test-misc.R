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

test_that("size format works", {
  expect_equal(size_humanreadable(123456), "121 kB")
})

test_that("functions missing in older R versions work", {
  expect_equal(strrep("A", 5), "AAAAA")
  expect_equal(strrep(c("A", "B"), c(5, 2)), c("AAAAA", "BB"))
  expect_equal(trimws(" test "), "test")
  expect_equal(trimws(" test ", "l"), "test ")
  expect_equal(trimws(" test ", "r"), " test")
})
