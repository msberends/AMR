# ==================================================================== #
# TITLE                                                                #
# AMR: An R Package for Working with Antimicrobial Resistance Data     #
#                                                                      #
# SOURCE                                                               #
# https://github.com/msberends/AMR                                     #
#                                                                      #
# CITE AS                                                              #
# Berends MS, Luz CF, Friedrich AW, Sinha BNM, Albers CJ, Glasner C    #
# (2022). AMR: An R Package for Working with Antimicrobial Resistance  #
# Data. Journal of Statistical Software, 104(3), 1-31.                 #
# doi:10.18637/jss.v104.i03                                            #
#                                                                      #
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

# extra tests for {vctrs} pkg support
if (pkg_is_available("dplyr", also_load = FALSE)) {
  test <- dplyr::tibble(ab = as.ab("CIP"),
                        mo = as.mo("Escherichia coli"),
                        mic = as.mic(2),
                        disk = as.disk(20),
                        rsi = as.rsi("S"))
  check1 <- lapply(test, class)
  test[1, "ab"] <- "GEN"
  test[1, "mo"] <- "B_KLBSL_PNMN"
  test[1, "mic"] <- ">=32"
  test[1, "mic"] <- 32
  test[1, "disk"] <- "35"
  test[1, "disk"] <- 25
  test[1, "disk"] <- 26L
  test[1, "rsi"] <- "R"
  check2 <- lapply(test, class)
  expect_identical(check1, check2)
  
  test <- dplyr::tibble(cipro = as.rsi("S"),
                        variable = "test")
  expect_equal(nrow(test[quinolones() == "S", ]), 1)
  expect_equal(nrow(test[quinolones() == "R", ]), 0)
}
