# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Analysis                              #
#                                                                      #
# SOURCE                                                               #
# https://gitlab.com/msberends/AMR                                     #
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
# Visit our website for more info: https://msberends.gitlab.io/AMR.    #
# ==================================================================== #

context("read.4d.R")

test_that("read 4D works", {

  library(dplyr)
  test1 <- data.frame(Patientnr = "ABC",
                      MV = "M",
                      Monsternr = "0123",
                      Afnamedat = "10-11-12",
                      Bepaling = "bk",
                      Afd. = "ABC",
                      Spec = "ABC",
                      Matbijz. = "ABC",
                      Mat = "ABC",
                      Mocode = "esccol",
                      PENI = "R",
                      stringsAsFactors = FALSE)
  tf <- tempfile()
  write.table(test1, file = tf, quote = F, sep = "\t")

  x <- read.4D(tf, skip = 0, info = TRUE)
  unlink(tf)

  expect_equal(ncol(x), 11)
  expect_equal(class(x$date_received), "Date")
  expect_equal(class(x$mo), "mo")
  expect_equal(as.character(x$mo), "B_ESCHR_COLI")
  expect_equal(is.rsi(x$peni), TRUE)

})
