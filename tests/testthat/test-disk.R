# ==================================================================== #
# TITLE                                                                #
# Antidiskrobial Resistance (AMR) Analysis                              #
#                                                                      #
# SOURCE                                                               #
# https://github.com/msberends/AMR                                     #
#                                                                      #
# LICENCE                                                              #
# (c) 2018-2020 Berends MS, Luz CF et al.                              #
#                                                                      #
# This R package is free software; you can freely use and distribute   #
# it for both personal and commercial purposes under the terms of the  #
# GNU General Public License version 2.0 (GNU GPL-2), as published by  #
# the Free Software Foundation.                                        #
#                                                                      #
# We created this package for both routine data analysis and academic  #
# research and it was publicly released in the hope that it will be    #
# useful, but it comes WITHOUT ANY WARRANTY OR LIABILITY.              #
# Visit our website for more info: https://msberends.github.io/AMR.    #
# ==================================================================== #

context("disk.R")

test_that("disk works", {
  skip_on_cran()
  expect_true(as.disk(8) == as.disk("8"))
  expect_true(is.disk(as.disk(8)))

  expect_equal(suppressWarnings(as.logical(as.disk("INVALID VALUE"))), NA)

  # all levels should be valid disks
  expect_silent(as.disk(levels(as.disk(15))))

  expect_warning(as.disk("INVALID VALUE"))
  
  expect_output(print(as.disk(12)))
  library(dplyr)
  expect_output(print(tibble(d = as.disk(12))))

})
