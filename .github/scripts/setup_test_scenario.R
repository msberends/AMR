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

# Sets up the user environment of a test scenario in the home directory, before run_tests.sh runs the tests. Run it through run_tests.sh with AMR_TEST_INTERACTIVE=true, since set_mo_source() only works in interactive sessions.
# AMR_TEST_SCENARIO can contain one or more of these words, separated by a comma:
# - "user": an existing ~/mo_source.rds, and the options AMR_custom_ab and AMR_custom_mo set to RDS files
# - "clsi": the option AMR_guideline set to "CLSI"
# The options are written to ~/.Rprofile, so that they are set at start-up as with users.

status <- tryCatch(
  {
    stopifnot(interactive())
    scenarios <- trimws(strsplit(Sys.getenv("AMR_TEST_SCENARIO"), ",")[[1]])
    stopifnot(length(scenarios) > 0, all(scenarios %in% c("user", "clsi")))
    pkgload::load_all(Sys.getenv("AMR_TEST_PKG_DIR", "."), quiet = TRUE)
    opts <- character(0)
    if ("user" %in% scenarios) {
      dir.create("~/amr_fixtures", showWarnings = FALSE)
      write.csv(data.frame(lab_code = c("AMRTEST_ECOLI", "AMRTEST_SAUR"), mo = c("B_ESCHR_COLI", "B_STPHY_AURS")), "~/amr_fixtures/mo_source.csv", row.names = FALSE)
      set_mo_source("~/amr_fixtures/mo_source.csv")
      stopifnot(file.exists("~/mo_source.rds"))
      saveRDS(data.frame(ab = "AMRTST", name = "Amrtestacillin", group = "Test antimicrobials"), "~/amr_fixtures/custom_ab.rds")
      saveRDS(data.frame(genus = "Amrtestbacterium", species = "scenarii"), "~/amr_fixtures/custom_mo.rds")
      opts <- c(opts, paste0("AMR_custom_ab = '", path.expand("~/amr_fixtures/custom_ab.rds"), "'"), paste0("AMR_custom_mo = '", path.expand("~/amr_fixtures/custom_mo.rds"), "'"))
    }
    if ("clsi" %in% scenarios) {
      opts <- c(opts, "AMR_guideline = 'CLSI'")
    }
    if (length(opts) > 0) {
      write(paste0("options(", paste(opts, collapse = ", "), ")"), "~/.Rprofile", append = TRUE)
    }
    cat("Scenario set up:", paste(scenarios, collapse = ", "), "\n")
    if (file.exists("~/.Rprofile")) cat(readLines("~/.Rprofile"), sep = "\n")
    0L
  },
  error = function(e) {
    cat("Scenario could not be set up:", conditionMessage(e), "\n")
    2L
  }
)
quit(save = "no", status = status)
