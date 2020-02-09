# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Analysis                              #
#                                                                      #
# SOURCE                                                               #
# https://gitlab.com/msberends/AMR                                     #
#                                                                      #
# LICENCE                                                              #
# (c) 2018-2020 Berends MS, Luz CF et al.                              #
#                                                                      #
# This R package is free software; you can freely use and distribute   #
# it for both personal and commercial purposes under the terms of the  #
# GNU General Public License version 2.0 (GNU GPL-2), as published by  #
# the Free Software Foundation.                                        #
#                                                                      #
# We created this package for both routine data analysis and academic  #
# research and it was publicly released in the hope that it will be    #
# useful, but it comes WITHOUT ANY WARRANTY OR LIABILITY.              #
# Visit our website for more info: https://msberends.gitlab.io/AMR.    #
# ==================================================================== #

# Run this file to update the package using: -------------------------------
# source("data-raw/internals.R")
# --------------------------------------------------------------------------

# See 'data-raw/eucast_rules.tsv' for the EUCAST reference file
eucast_rules_file <- utils::read.delim(file = "data-raw/eucast_rules.tsv",
                            skip = 10,
                            sep = "\t",
                            stringsAsFactors = FALSE,
                            header = TRUE,
                            strip.white = TRUE,
                            na = c(NA, "", NULL))
# take the order of the reference.rule_group column in the orginal data file
eucast_rules_file$reference.rule_group <- factor(eucast_rules_file$reference.rule_group,
                                                 levels = unique(eucast_rules_file$reference.rule_group),
                                                 ordered = TRUE)
eucast_rules_file <- dplyr::arrange(eucast_rules_file,
  reference.rule_group,
  reference.rule)
eucast_rules_file$reference.rule_group <- as.character(eucast_rules_file$reference.rule_group)

# Translations ----
translations_file <- utils::read.delim(file = "data-raw/translations.tsv",
                                       sep = "\t",
                                       stringsAsFactors = FALSE,
                                       header = TRUE,
                                       blank.lines.skip = TRUE,
                                       fill = TRUE,
                                       strip.white = TRUE,
                                       encoding = "UTF-8",
                                       fileEncoding = "UTF-8",
                                       na.strings = c(NA, "", NULL),
                                       allowEscapes = TRUE, # else "\\1" will be imported as "\\\\1"
                                       quote = "")

# Old microorganism codes -------------------------------------------------

microorganisms.translation <- readRDS("data-raw/microorganisms.translation.rds")

# Export to package as internal data ----
usethis::use_data(eucast_rules_file, translations_file, microorganisms.translation,
                  internal = TRUE,
                  overwrite = TRUE,
                  version = 2)

# Remove from global environment ----
rm(eucast_rules_file)
rm(translations_file)
rm(microorganisms.translation)

# Save to raw data to repository ----
library(dplyr, warn.conflicts = FALSE, quietly = TRUE)
usethis::ui_done(paste0("Saving raw data to {usethis::ui_value('/data-raw/')}"))
devtools::load_all(quiet = TRUE)
write.table(AMR::rsi_translation,
            "data-raw/rsi_translation.txt", sep = "\t", na = "", row.names = FALSE)
write.table(microorganisms %>% mutate_if(~!is.numeric(.), as.character),
            "data-raw/microorganisms.txt", sep = "\t", na = "", row.names = FALSE)
write.table(antibiotics %>% mutate_if(~!is.numeric(.), as.character),
            "data-raw/antibiotics.txt", sep = "\t", na = "", row.names = FALSE)
write.table(antivirals %>% mutate_all(as.character),
            "data-raw/antivirals.txt", sep = "\t", na = "", row.names = FALSE)
