# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Data Analysis for R                   #
#                                                                      #
# SOURCE                                                               #
# https://github.com/msberends/AMR                                     #
#                                                                      #
# LICENCE                                                              #
# (c) 2018-2021 Berends MS, Luz CF et al.                              #
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

# test only on GitHub Actions and at home - not on CRAN as tests are lengthy
if (identical(Sys.getenv("R_RUN_TINYTEST"), "true")) {
  # env var 'R_LIBS_USER' got overwritten during 'R CMD check' in GitHub Actions, so:
  .libPaths(c(Sys.getenv("R_LIBS_USER_GH_ACTIONS"), .libPaths()))
  # helper function
  pkg_is_available <- function(pkg, also_load = TRUE) {
    if (also_load == TRUE) {
      out <- suppressWarnings(require(pkg, character.only = TRUE, warn.conflicts = FALSE, quietly = TRUE))
    } else {
      out <- requireNamespace(pkg, quietly = TRUE)
    }
    isTRUE(out)
  }
  
  if (pkg_is_available("tinytest")) {
    library(AMR)
    out <- test_package("AMR")
  }
}
