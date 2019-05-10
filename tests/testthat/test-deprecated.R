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

context("deprecated.R")

test_that("deprecated functions work", {

  expect_error(suppressWarnings(ratio("A")))
  expect_error(suppressWarnings(ratio(1, ratio = "abc")))
  expect_error(suppressWarnings(ratio(c(1, 2), ratio = c(1, 2, 3))))
  expect_warning(ratio(c(772, 1611, 737), ratio = "1:2:1"))
  expect_identical(suppressWarnings(ratio(c(772, 1611, 737), ratio = "1:2:1")), c(780, 1560,  780))
  expect_identical(suppressWarnings(ratio(c(1752, 1895), ratio = c(1, 1))), c(1823.5, 1823.5))

  expect_warning(atc_property("amox"))
  expect_warning(atc_official("amox"))
  expect_warning(ab_official("amox"))
  expect_warning(atc_name("amox"))
  expect_warning(atc_trivial_nl("amox"))
  expect_warning(atc_tradenames("amox"))

})
