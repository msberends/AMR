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
# Center Groningen in the Netherlands, in collaboration with many      #
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

# Runs all unit tests from source, see run_tests.sh. Quits with status 1 on any failure or error, and with status 2 if the tests could not be run at all.
# Every top-level error must end in quit(), since an interactive session would otherwise continue reading "1" from stdin.

status <- tryCatch(
  {
    if (identical(Sys.getenv("AMR_TEST_INTERACTIVE"), "true") && !interactive()) {
      stop("AMR_TEST_INTERACTIVE is set, but this R session is not interactive")
    }
    cat("Interactive:", interactive(), "| LANG:", Sys.getenv("LANG"), "| AMR options:", paste(names(options())[grepl("^AMR_", names(options()))], collapse = ", "), "\n")
    # report all failures, not only the first 10
    options(warn = 1, testthat.summary.max_reports = Inf, testthat.progress.max_fails = Inf)
    res <- as.data.frame(testthat::test_local(Sys.getenv("AMR_TEST_PKG_DIR", "."), reporter = "summary", stop_on_failure = FALSE))
    cat("\nFiles:", length(unique(res$file)), "| failed:", sum(res$failed), "| errors:", sum(res$error), "| skipped:", sum(res$skipped), "| passed:", sum(res$passed), "\n")
    as.integer(sum(res$failed) > 0 || any(res$error))
  },
  error = function(e) {
    cat("Tests could not be run:", conditionMessage(e), "\n")
    2L
  }
)
quit(save = "no", status = status)
