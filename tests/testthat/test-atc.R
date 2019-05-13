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

context("ab.R")

test_that("as.atc works", {
  expect_identical(class(as.atc("amox")), "atc")
  expect_true(is.atc(as.atc("amox")))
  expect_output(print(as.atc("amox")))
  expect_output(print(data.frame(a = as.atc("amox"))))

  expect_identical(class(pull(antibiotics, atc)), "atc")

  expect_warning(as.atc("Z00ZZ00")) # not yet availatcle in data set
  expect_warning(as.atc("UNKNOWN"))

  expect_output(print(as.atc("amox")))



})
