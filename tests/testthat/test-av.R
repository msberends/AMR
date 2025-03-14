# ==================================================================== #
# TITLE:                                                               #
# AMR: An R Package for Working with Antimicrobial Resistance Data     #
#                                                                      #
# SOURCE CODE:                                                         #
# https://github.com/msberends/AMR                                     #
#                                                                      #
# PLEASE CITE THIS SOFTWARE AS:                                        #
# Berends MS, Luz CF, Friedrich AW, et al. (2022).                     #
# AMR: An R Package for Working with Antimicrobial Resistance Data.    #
# Journal of Statistical Software, 104(3), 1-31.                       #
# https://doi.org/10.18637/jss.v104.i03                                #
#                                                                      #
# Developed at the University of Groningen and the University Medical  #
# Center Groningen in The Netherlands, in collaboration with many      #
# colleagues from around the world, see our website.                   #
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

test_that("test-av.R", {
  expect_equal(
    as.character(as.av(c(
      "J05AB01",
      "J 05 AB 01",
      "Aciclovir",
      "aciclo",
      "   aciclo 123",
      "ACICL",
      "ACI",
      "Virorax",
      "Zovirax"
    ))),
    rep("ACI", 9)
  )

  expect_identical(class(as.av("acic")), c("av", "character"))
  expect_identical(class(antivirals$av), c("av", "character"))
  expect_true(is.av(as.av("acic")))
  expect_output(print(as.av("acic")))
  expect_output(print(data.frame(a = as.av("acic"))))

  # expect_warning(as.av("J00AA00")) # ATC not yet available in data set
  # expect_warning(as.av("UNKNOWN"))

  expect_output(print(as.av("acic")))

  expect_equal(
    as.character(as.av("zovirax")),
    "ACI"
  )

  expect_equal(
    as.character(as.av(c("Abacaivr", "Celvudine"))),
    c("ABA", "CLE")
  )

  # expect_warning(as.av("Abacavir Clevudine"))

  # based on Levenshtein distance
  expect_identical(av_name("adevofir dypifo", language = NULL), "Adefovir dipivoxil")

  # assigning and subsetting
  x <- antivirals$av
  expect_inherits(x[1], "av")
  expect_inherits(x[[1]], "av")
  expect_inherits(c(x[1], x[9]), "av")
  expect_inherits(unique(x[1], x[9]), "av")
  expect_inherits(rep(x[1], 2), "av")
  # expect_warning(x[1] <- "invalid code")
  # expect_warning(x[[1]] <- "invalid code")
  # expect_warning(c(x[1], "test"))
})
