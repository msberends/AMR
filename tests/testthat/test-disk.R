# ==================================================================== #
# TITLE                                                                #
# Antidiskrobial Resistance (AMR) Analysis                              #
#                                                                      #
# SOURCE                                                               #
# https://github.com/msberends/AMR                                     #
#                                                                      #
# LICENCE                                                              #
# (c) 2018-2021 Berends MS, Luz CF et al.                              #
# Developed at the University of Groningen, the Netherlands, in        #
# collaboration with non-profit organisations Certe Medical            #
# Diagnostics & Advice, and University Medical Center Groningen.       # 
#                                                                      #
# This R package is free software; you can freely use and distribute   #
# it for both personal and commercial purposes under the terms of the  #
# GNU General Public License version 2.0 (GNU GPL-2), as published by  #
# the Free Software Foundation.                                        #
# We created this package for both routine data analysis and academic  #
# research and it was publicly released in the hope that it will be    #
# useful, but it comes WITHOUT ANY WARRANTY OR LIABILITY.              #
#                                                                      #
# Visit our website for the full manual and a complete tutorial about  #
# how to conduct AMR data analysis: https://msberends.github.io/AMR/   #
# ==================================================================== #

context("disk.R")

test_that("disk works", {
  skip_on_cran()
  expect_true(as.disk(8) == as.disk("8"))
  expect_true(is.disk(as.disk(8)))

  expect_equal(suppressWarnings(as.logical(as.disk("INVALID VALUE"))), NA)

  # all levels should be valid disks
  x <- as.disk(c(20, 40))
  expect_s3_class(x[1], "disk")
  expect_s3_class(x[[1]], "disk")
  expect_s3_class(c(x[1], x[9]), "disk")
  expect_s3_class(unique(x[1], x[9]), "disk")
  expect_warning(as.disk("INVALID VALUE"))
  x[2] <- 32
  expect_s3_class(x, "disk")
  
  pdf(NULL) # prevent Rplots.pdf being created
  expect_silent(plot(as.disk(c(10, 20, 40))))

  expect_output(print(as.disk(12)))
  library(dplyr, warn.conflicts = FALSE)
  expect_output(print(tibble(d = as.disk(12))))

})
