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

# Recreates the data sets that are computed from the microorganisms data set, in the right order. Run this after every
# rebuild of the microorganisms data set (run_microorganisms_build.R), from the root of the repository:
#   Rscript data-raw/_reproduction_scripts/run_after_microorganisms_build.R
#
# 1. microorganisms.groups: its members are taken from the current names in the new data set
# 2. clinical_breakpoints: the human EUCAST breakpoints are parsed again from the EUCAST workbooks, as their organisms
#    are resolved against the new data set (the WHONET-based rows are kept, with their rank_index updated)
# 3. intrinsic_resistant: computed from all codes of the new data set and the interpretive rules, which use the
#    species groups of step 1 (so this must be last); it is never an input of the microorganisms build itself

if (!file.exists("data-raw/_reproduction_scripts/run_after_microorganisms_build.R")) {
  stop("Run this from the root of the AMR repository", call. = FALSE)
}
step <- function(name, file) {
  message("\n=== ", name, " ===")
  started <- Sys.time()
  exit <- system2("Rscript", file)
  if (exit != 0) {
    stop(name, " failed (exit code ", exit, ")", call. = FALSE)
  }
  message(name, " done in ", round(difftime(Sys.time(), started, units = "mins"), 1), " minutes")
}

# 1. species groups
step("microorganisms.groups", "data-raw/_reproduction_scripts/reproduction_of_microorganisms.groups.R")

# 2. EUCAST breakpoints, merged into clinical_breakpoints (as at the end of reproduction_of_clinical_breakpoints.R)
merge_script <- tempfile(fileext = ".R")
writeLines(c(
  'devtools::load_all()',
  'library(dplyr, warn.conflicts = FALSE)',
  'source("data-raw/_reproduction_scripts/reproduction_of_clinical_breakpoints_eucast.R")',
  '# the rank of an organism can have changed with the new data set',
  'breakpoints_now <- clinical_breakpoints %>%',
  '  mutate(rank_index = case_when(',
  '    mo_rank(mo, keep_synonyms = TRUE) %like% "(infra|sub)" ~ 1,',
  '    mo_rank(mo, keep_synonyms = TRUE) == "species" ~ 2,',
  '    mo_rank(mo, keep_synonyms = TRUE) == "species group" ~ 2.5,',
  '    mo_rank(mo, keep_synonyms = TRUE) == "genus" ~ 3,',
  '    mo_rank(mo, keep_synonyms = TRUE) == "family" ~ 4,',
  '    mo_rank(mo, keep_synonyms = TRUE) == "order" ~ 5,',
  '    mo != "UNKNOWN" ~ 6,',
  '    TRUE ~ 7',
  '  ) + if_else(is.na(breakpoint_S) & is.na(breakpoint_R) & guideline %like% "EUCAST", 0.1, 0))',
  'clinical_breakpoints <- merge_eucast_breakpoints(breakpoints_now, breakpoints_eucast, eucast_coverage)',
  'usethis::use_data(clinical_breakpoints, overwrite = TRUE, compress = "xz", version = 2)'
), merge_script)
step("clinical_breakpoints (EUCAST)", merge_script)

# 3. intrinsic resistance
step("intrinsic_resistant", "data-raw/_reproduction_scripts/reproduction_of_intrinsic_resistant.R")

message("\nAll data sets that depend on the microorganisms data set are recreated. Run the unit tests next.")
