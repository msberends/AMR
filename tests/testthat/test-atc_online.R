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

test_that("test-atc_online.R", {
  skip_on_cran()

  if (AMR:::pkg_is_available("curl") &&
    AMR:::pkg_is_available("rvest") &&
    AMR:::pkg_is_available("xml2") &&
    tryCatch(curl::has_internet(), error = function(e) FALSE)) {
    expect_true(length(atc_online_groups(ab_atc("AMX"))) >= 1)
    expect_equal(atc_online_ddd(ab_atc("AMX"), administration = "O"), 1.5)
    expect_equal(atc_online_ddd(ab_atc("AMX"), administration = "P"), 3)
    expect_equal(atc_online_ddd_units("AMX", administration = "P"), "g")
  }
})
