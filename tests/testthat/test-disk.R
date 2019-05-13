# ==================================================================== #
# TITLE                                                                #
# Antidiskrobial Resistance (AMR) Analysis                              #
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

context("disk.R")

test_that("disk works", {
  expect_true(as.disk(8) == as.disk("8"))
  expect_true(is.disk(as.disk(8)))

  expect_equal(suppressWarnings(as.logical(as.disk("INVALID VALUE"))), NA)

  # all levels should be valid disks
  expect_silent(as.disk(levels(as.disk(15))))

  expect_warning(as.disk("INVALID VALUE"))

})
