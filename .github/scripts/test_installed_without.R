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

# Run all unit tests and all examples against the INSTALLED AMR package, in a library where some suggested packages are not available.
# Used by .github/workflows/check-without-dplyr.yaml; run from the root of the repository:
#   Rscript .github/scripts/test_installed_without.R dplyr tidyselect
# The package must have been installed BEFORE the listed packages were made unavailable, as binary packages from CRAN are built on machines that have all suggested packages: code that depends on what was available at installation (instead of at runtime) then fails, which R CMD check cannot detect as it reinstalls the package.

unavailable <- commandArgs(trailingOnly = TRUE)
if (length(unavailable) == 0) {
  stop("provide the package(s) that must be unavailable, e.g.: Rscript test_installed_without.R dplyr tidyselect")
}
still_available <- unavailable[vapply(unavailable, requireNamespace, logical(1), quietly = TRUE)]
if (length(still_available) > 0) {
  stop("these packages must not be available, but are: ", paste(still_available, collapse = ", "))
}
cat("Not available:", paste(unavailable, collapse = ", "), "\n")
cat("AMR", format(utils::packageVersion("AMR")), "installed in", find.package("AMR"), "\n\n")

Sys.setenv(NOT_CRAN = "true")
library(AMR)

# 1. all unit tests, against the installed package
cat("== Unit tests ==\n")
results <- as.data.frame(testthat::test_dir(
  "tests/testthat",
  package = "AMR",
  load_package = "installed",
  reporter = "summary",
  stop_on_failure = FALSE
))
n_failed <- sum(results$failed) + sum(results$error)
cat("\nTests: ", sum(results$nb), " expectations, ", sum(results$failed), " failed, ", sum(results$error), " errored, ", sum(results$skipped), " skipped\n", sep = "")

# 2. all examples, including \donttest{} (not \dontrun{}), against the installed package
cat("\n== Examples ==\n")
rd_db <- tools::Rd_db("AMR")
failed_examples <- character(0)
n_examples <- 0
for (rd in sort(names(rd_db))) {
  ex_file <- tempfile(fileext = ".R")
  tools::Rd2ex(rd_db[[rd]], ex_file, commentDonttest = FALSE)
  if (!file.exists(ex_file)) next
  n_examples <- n_examples + 1
  # print the name first, so that a crash of R itself shows which help page caused it
  cat(rd, "\n")
  ok <- tryCatch(
    {
      utils::capture.output(suppressMessages(suppressWarnings(source(ex_file, local = new.env(parent = globalenv()), echo = FALSE))))
      TRUE
    },
    error = function(e) {
      cat("FAILED: ", rd, ": ", conditionMessage(e), "\n", sep = "")
      FALSE
    }
  )
  if (!ok) failed_examples <- c(failed_examples, rd)
  grDevices::graphics.off()
}
cat("Examples: ", length(failed_examples), " of ", n_examples, " help pages failed\n", sep = "")

if (n_failed > 0 || length(failed_examples) > 0) {
  stop(n_failed, " test(s) failed and ", length(failed_examples), " example(s) failed without ", paste(unavailable, collapse = " and "), call. = FALSE)
}
cat("\nAll tests and examples passed without ", paste(unavailable, collapse = " and "), ".\n", sep = "")
