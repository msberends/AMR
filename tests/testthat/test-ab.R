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

context("ab.R")

test_that("as.ab works", {
  expect_equal(as.character(as.ab(c("J01FA01",
                                    "J 01 FA 01",
                                    "Erythromycin",
                                    "eryt",
                                    "   eryt 123",
                                    "ERYT",
                                    "ERY",
                                    "erytromicine",
                                    "Erythrocin",
                                    "Romycin"))),
               rep("ERY", 10))

  expect_identical(class(as.ab("amox")), "ab")
  expect_identical(class(pull(antibiotics, ab)), "ab")
  expect_true(is.ab(as.ab("amox")))
  expect_output(print(as.ab("amox")))
  expect_output(print(data.frame(a = as.ab("amox"))))

  expect_warning(as.ab("Z00ZZ00")) # not yet available in data set
  expect_warning(as.ab("UNKNOWN"))
  expect_warning(as.ab(""))

  expect_output(print(as.ab("amox")))

  expect_identical(class(pull(antibiotics, ab)), "ab")

  expect_equal(as.character(as.ab("Phloxapen")),
               "FLC")

  expect_equal(suppressWarnings(as.character(as.ab(c("Bacteria", "Bacterial")))),
               c(NA, "TMP"))
  
  expect_equal(as.character(as.ab("Amoxy + clavulaanzuur")),
               "AMC")

  # assigning and subsetting
  x <- antibiotics$ab
  expect_s3_class(x[1], "ab")
  expect_s3_class(x[[1]], "ab")
  expect_s3_class(c(x[1], x[9]), "ab")
  expect_warning(x[1] <- "invalid code")
  expect_warning(x[[1]] <- "invalid code")
  expect_warning(c(x[1], "test"))
})
