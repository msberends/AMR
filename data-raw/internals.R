# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Analysis                              #
#                                                                      #
# SOURCE                                                               #
# https://github.com/msberends/AMR                                     #
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
# Visit our website for more info: https://msberends.github.io/AMR.    #
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
# take the order of the reference.rule_group column in the original data file
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
usethis::ui_done(paste0("Saving raw data to {usethis::ui_value('/data-raw/')}"))
devtools::load_all(quiet = TRUE)
# give official names to ABs and MOs
rsi <- dplyr::mutate(rsi_translation, ab = ab_name(ab), mo = mo_name(mo))
saveRDS(rsi, "data-raw/rsi_translation.rds", version = 2)
write.table(rsi, "data-raw/rsi_translation.txt", sep = "\t", na = "", row.names = FALSE)
haven::write_sas(rsi, "data-raw/rsi_translation.sas")
haven::write_sav(rsi, "data-raw/rsi_translation.sav")
haven::write_dta(rsi, "data-raw/rsi_translation.dta")
openxlsx::write.xlsx(rsi, "data-raw/rsi_translation.xlsx")

mo <- dplyr::mutate_if(microorganisms, ~!is.numeric(.), as.character)
saveRDS(mo, "data-raw/microorganisms.rds", version = 2)
write.table(mo, "data-raw/microorganisms.txt", sep = "\t", na = "", row.names = FALSE)
haven::write_sas(mo, "data-raw/microorganisms.sas")
haven::write_sav(mo, "data-raw/microorganisms.sav")
haven::write_dta(mo, "data-raw/microorganisms.dta")
openxlsx::write.xlsx(mo, "data-raw/microorganisms.xlsx")

saveRDS(microorganisms.old, "data-raw/microorganisms.old.rds", version = 2)
write.table(microorganisms.old, "data-raw/microorganisms.old.txt", sep = "\t", na = "", row.names = FALSE)
haven::write_sas(microorganisms.old, "data-raw/microorganisms.old.sas")
haven::write_sav(microorganisms.old, "data-raw/microorganisms.old.sav")
haven::write_dta(microorganisms.old, "data-raw/microorganisms.old.dta")
openxlsx::write.xlsx(microorganisms.old, "data-raw/microorganisms.old.xlsx")

ab <- dplyr::mutate_if(antibiotics, ~!is.numeric(.), as.character)
saveRDS(ab, "data-raw/antibiotics.rds", version = 2)
write.table(ab, "data-raw/antibiotics.txt", sep = "\t", na = "", row.names = FALSE)
haven::write_sas(ab, "data-raw/antibiotics.sas")
haven::write_sav(ab, "data-raw/antibiotics.sav")
haven::write_dta(ab, "data-raw/antibiotics.dta")
openxlsx::write.xlsx(ab, "data-raw/antibiotics.xlsx")

av <- dplyr::mutate_if(antivirals, ~!is.numeric(.), as.character)
saveRDS(av, "data-raw/antivirals.rds", version = 2)
write.table(av, "data-raw/antivirals.txt", sep = "\t", na = "", row.names = FALSE)
haven::write_sas(av, "data-raw/antivirals.sas")
haven::write_sav(av, "data-raw/antivirals.sav")
haven::write_dta(av, "data-raw/antivirals.dta")
openxlsx::write.xlsx(av, "data-raw/antivirals.xlsx")

saveRDS(intrinsic_resistant, "data-raw/intrinsic_resistant.rds", version = 2)
write.table(intrinsic_resistant, "data-raw/intrinsic_resistant.txt", sep = "\t", na = "", row.names = FALSE)
haven::write_sas(intrinsic_resistant, "data-raw/intrinsic_resistant.sas")
haven::write_sav(intrinsic_resistant, "data-raw/intrinsic_resistant.sav")
haven::write_dta(intrinsic_resistant, "data-raw/intrinsic_resistant.dta")
openxlsx::write.xlsx(intrinsic_resistant, "data-raw/intrinsic_resistant.xlsx")
