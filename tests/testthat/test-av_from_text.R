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
# how to conduct AMR data analysis: https://amr-for-r.org              #
# ==================================================================== #

test_that("test-av_from_text.R", {
  skip_on_cran()

  expect_identical(
    av_from_text("28/03/2020 regular aciclovir 500mg po tds")[[1]],
    as.av("Aciclovir")
  )
  expect_identical(
    av_from_text("28/03/2020 regular aciclovir 500mg po tds", thorough_search = TRUE)[[1]],
    as.av("Aciclovir")
  )
  expect_identical(
    av_from_text("28/03/2020 regular aciclovir 500mg po tds", thorough_search = FALSE)[[1]],
    as.av("Aciclovir")
  )
  expect_identical(
    av_from_text("28/03/2020 regular aciclovir 500mg po tds", translate_av = TRUE)[[1]],
    "Aciclovir"
  )
  expect_identical(
    av_from_text("administered aciclo and valaciclo", collapse = ", ")[[1]],
    "ACI, VALA"
  )

  expect_identical(
    av_from_text("28/03/2020 regular aciclo 500mg po tds", type = "dose")[[1]],
    500
  )
  expect_identical(
    av_from_text("28/03/2020 regular aciclo 500mg po tds", type = "admin")[[1]],
    "oral"
  )
})
