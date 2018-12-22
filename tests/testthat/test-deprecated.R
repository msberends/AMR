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

context("deprecated.R")

test_that("deprecated functions work", {

  expect_error(suppressWarnings(ratio("A")))
  expect_error(suppressWarnings(ratio(1, ratio = "abc")))
  expect_error(suppressWarnings(ratio(c(1, 2), ratio = c(1, 2, 3))))
  expect_warning(ratio(c(772, 1611, 737), ratio = "1:2:1"))
  expect_identical(suppressWarnings(ratio(c(772, 1611, 737), ratio = "1:2:1")), c(780, 1560,  780))
  expect_identical(suppressWarnings(ratio(c(1752, 1895), ratio = c(1, 1))), c(1823.5, 1823.5))

})
