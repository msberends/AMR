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

context("kurtosis.R")

test_that("kurtosis works", {
  expect_equal(kurtosis(septic_patients$age),
               3.549319,
               tolerance = 0.00001)
  expect_equal(unname(kurtosis(data.frame(septic_patients$age))),
               3.549319,
               tolerance = 0.00001)
  expect_equal(kurtosis(matrix(septic_patients$age)),
               3.549319,
               tolerance = 0.00001)
})
