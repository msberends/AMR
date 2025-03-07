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


test_that("custom mo works", {
  expect_identical(
    as.mo("Enterobacter asburiae/cloacae"),
    as.mo("Enterobacter asburiae")
  )

  suppressMessages(
    add_custom_microorganisms(
      data.frame(
        mo = "ENT_ASB_CLO",
        genus = "Enterobacter",
        species = "asburiae/cloacae"
      )
    )
  )

  expect_identical(as.character(as.mo("ENT_ASB_CLO")), "ENT_ASB_CLO")
  expect_identical(mo_name("ENT_ASB_CLO"), "Enterobacter asburiae/cloacae")
  expect_identical(mo_gramstain("ENT_ASB_CLO", language = NULL), "Gram-negative")

  if (getRversion() >= "3.3.0") {
    # until R 3.2, abbreviate() used a completely different algorithm, making these tests unreproducible
    expect_identical(
      paste("B", AMR:::abbreviate_mo("Klebsiella"), AMR:::abbreviate_mo("pneumoniae", 4), sep = "_"),
      as.character(as.mo("Klebsiella pneumoniae"))
    )
    expect_identical(
      paste("B", AMR:::abbreviate_mo("Aerococcus"), AMR:::abbreviate_mo("urinae", 4), sep = "_"),
      as.character(as.mo("Aerococcus urinae"))
    )
  }
})
