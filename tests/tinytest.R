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

# we use {tinytest} instead of {testthat} because it does not rely on recent R versions - we want to test on R >= 3.0.

# Run them in RStudio using:
# rstudioapi::jobRunScript("tests/tinytest.R", name = "AMR Unit Tests", workingDir = getwd(), exportEnv = "tinytest_results")

# test only on GitHub Actions and at using RStudio jobs - not on CRAN as tests are lengthy
if (tryCatch(isTRUE(AMR:::import_fn("isJob", "rstudioapi")()), error = function(e) FALSE) ||
  identical(Sys.getenv("R_RUN_TINYTEST"), "true")) {
  # env var 'R_LIBS_USER' got overwritten during 'R CMD check' in GitHub Actions, so:
  .libPaths(c(Sys.getenv("R_LIBS_USER_GH_ACTIONS"), .libPaths()))
  if (AMR:::pkg_is_available("tinytest", also_load = TRUE)) {
    library(AMR)
    # set language
    set_AMR_locale("English")
    # set some functions if on old R
    if (getRversion() < "3.2.0") {
      anyNA <- AMR:::anyNA
      dir.exists <- AMR:::dir.exists
      file.size <- AMR:::file.size
      file.mtime <- AMR:::file.mtime
      isNamespaceLoaded <- AMR:::isNamespaceLoaded
      lengths <- AMR:::lengths
    }
    if (getRversion() < "3.3.0") {
      strrep <- AMR:::strrep
    }
    if (getRversion() < "3.5.0") {
      isFALSE <- AMR:::isFALSE
    }
    if (getRversion() < "3.6.0") {
      str2lang <- AMR:::str2lang
      # trims() was introduced in 3.3.0, but its argument `whitespace` only in 3.6.0
      trimws <- AMR:::trimws
    }
    if (getRversion() < "4.0.0") {
      deparse1 <- AMR:::deparse1
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
