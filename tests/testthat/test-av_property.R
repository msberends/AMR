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

test_that("test-ab_property.R", {
  skip_on_cran()

  expect_identical(av_name("ACI", language = NULL), "Aciclovir")
  expect_identical(av_atc("ACI"), "J05AB01")
  expect_identical(av_cid("ACI"), as.integer(135398513))

  expect_inherits(av_tradenames("ACI"), "character")
  expect_inherits(av_tradenames(c("ACI", "ACI")), "list")

  expect_identical(av_group("ACI", language = NULL), "Nucleosides and nucleotides excl. reverse transcriptase inhibitors")

  expect_identical(av_name(135398513, language = NULL), "Aciclovir")
  expect_identical(av_name("J05AB01", language = NULL), "Aciclovir")

  expect_identical(av_ddd("ACI", "oral"), 4)
  expect_identical(av_ddd_units("ACI", "iv"), "g")
  expect_identical(av_ddd("ACI", "iv"), 4)

  expect_identical(
    av_name(x = c("ACI", "VALA"), tolower = TRUE, language = NULL),
    c("aciclovir", "valaciclovir")
  )

  expect_inherits(av_info("ACI"), "list")

  expect_error(av_property("acic", "invalid property"))
  expect_error(av_name("acic", language = "INVALID"))
  expect_output(print(av_name("acic", language = NULL)))

  expect_equal(av_name("29113-8", language = NULL), "Abacavir")
  expect_equal(
    av_loinc("Abacavir"),
    c("29113-8", "30273-7", "30287-7", "30303-2", "78772-1", "78773-9", "79134-3", "80118-3")
  )

  expect_true(av_url("ACI") %like% "fhi[.]no")
})
