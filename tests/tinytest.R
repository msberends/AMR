# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Data Analysis for R                   #
#                                                                      #
# SOURCE                                                               #
# https://github.com/msberends/AMR                                     #
#                                                                      #
# LICENCE                                                              #
# (c) 2018-2022 Berends MS, Luz CF et al.                              #
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

# we use {tinytest} instead of {testthat} because it does not rely on recent R versions - we want to test on R >= 3.0.

# Run them in RStudio using:
# rstudioapi::jobRunScript("tests/tinytest.R", name = "Tinytest Unit Tests", workingDir = getwd(), exportEnv = "tinytest_results")

# test only on GitHub Actions and at home - not on CRAN as tests are lengthy
if (identical(Sys.getenv("R_RUN_TINYTEST"), "true")) {
  # env var 'R_LIBS_USER' got overwritten during 'R CMD check' in GitHub Actions, so:
  .libPaths(c(Sys.getenv("R_LIBS_USER_GH_ACTIONS"), .libPaths()))
  if (AMR:::pkg_is_available("tinytest", also_load = TRUE)) {
    library(AMR)
    # set language
    set_AMR_locale("English")
    # get trimws() and strrep() if on old R
    if (getRversion() < "3.3.0") {
      trimws <- AMR:::trimws
      strrep <- AMR:::strrep
    }
    # start the unit tests
    out <- test_package("AMR",
      testdir = ifelse(dir.exists("inst/tinytest"),
        "inst/tinytest",
        "tinytest"
      ),
      verbose = 99,
      color = FALSE
    )
    cat("\n\nSUMMARY:\n")
    print(summary(out))
  }
}

