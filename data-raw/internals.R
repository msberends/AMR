# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Analysis for R                        #
#                                                                      #
# SOURCE                                                               #
# https://github.com/msberends/AMR                                     #
#                                                                      #
# LICENCE                                                              #
# (c) 2018-2020 Berends MS, Luz CF et al.                              #
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
# how to conduct AMR analysis: https://msberends.github.io/AMR/        #
# ==================================================================== #

# Run this file to update the package using: -------------------------------
# source("data-raw/internals.R")
# --------------------------------------------------------------------------

# See 'data-raw/eucast_rules.tsv' for the EUCAST reference file
library(dplyr, warn.conflicts = FALSE)
eucast_rules_file <- utils::read.delim(file = "data-raw/eucast_rules.tsv",
                                       skip = 10,
                                       sep = "\t",
                                       stringsAsFactors = FALSE,
                                       header = TRUE,
                                       strip.white = TRUE,
                                       na = c(NA, "", NULL)) %>% 
  # take the order of the reference.rule_group column in the original data file
  mutate(reference.rule_group = factor(reference.rule_group,
                                       levels = unique(reference.rule_group),
                                       ordered = TRUE),
         sorting_rule = ifelse(grepl("^Table", reference.rule, ignore.case = TRUE), 1, 2)) %>% 
  arrange(reference.rule_group,
          reference.version,
          sorting_rule,
          reference.rule) %>% 
  mutate(reference.rule_group = as.character(reference.rule_group)) %>% 
  select(-sorting_rule)

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
                  version = 2,
                  compress = "xz")

# Remove from global environment ----
rm(eucast_rules_file)
rm(translations_file)
rm(microorganisms.translation)

# Save to raw data to repository ----
write_md5 <- function(object) {
  conn <- file(paste0("data-raw/", deparse(substitute(object)), ".md5"))
  writeLines(digest::digest(object, "md5"), conn)
  close(conn)
}
changed_md5 <- function(object) {
  tryCatch({
    conn <- file(paste0("data-raw/", deparse(substitute(object)), ".md5"))
    compared <- digest::digest(object, "md5") != readLines(con = conn)
    close(conn)
    compared
  }, error = function(e) TRUE)
}
usethis::ui_done(paste0("Saving raw data to {usethis::ui_value('/data-raw/')}"))
devtools::load_all(quiet = TRUE)
# give official names to ABs and MOs
rsi <- dplyr::mutate(rsi_translation, ab = ab_name(ab), mo = mo_name(mo))
if (changed_md5(rsi)) {
  write_md5(rsi)
  try(saveRDS(rsi, "data-raw/rsi_translation.rds", version = 2, compress = "xz"), silent = TRUE)
  try(write.table(rsi, "data-raw/rsi_translation.txt", sep = "\t", na = "", row.names = FALSE), silent = TRUE)
  try(haven::write_sas(rsi, "data-raw/rsi_translation.sas"), silent = TRUE)
  try(haven::write_sav(rsi, "data-raw/rsi_translation.sav"), silent = TRUE)
  try(haven::write_dta(rsi, "data-raw/rsi_translation.dta"), silent = TRUE)
  try(openxlsx::write.xlsx(rsi, "data-raw/rsi_translation.xlsx"), silent = TRUE)
}

mo <- dplyr::mutate_if(microorganisms, ~!is.numeric(.), as.character)
if (changed_md5(mo)) {
  write_md5(mo)
  try(saveRDS(mo, "data-raw/microorganisms.rds", version = 2, compress = "xz"), silent = TRUE)
  try(write.table(mo, "data-raw/microorganisms.txt", sep = "\t", na = "", row.names = FALSE), silent = TRUE)
  try(haven::write_sas(mo, "data-raw/microorganisms.sas"), silent = TRUE)
  try(haven::write_sav(mo, "data-raw/microorganisms.sav"), silent = TRUE)
  try(haven::write_dta(mo, "data-raw/microorganisms.dta"), silent = TRUE)
  try(openxlsx::write.xlsx(mo, "data-raw/microorganisms.xlsx"), silent = TRUE)
}

if (changed_md5(microorganisms.old)) {
  write_md5(microorganisms.old)
  try(saveRDS(microorganisms.old, "data-raw/microorganisms.old.rds", version = 2, compress = "xz"), silent = TRUE)
  try(write.table(microorganisms.old, "data-raw/microorganisms.old.txt", sep = "\t", na = "", row.names = FALSE), silent = TRUE)
  try(haven::write_sas(microorganisms.old, "data-raw/microorganisms.old.sas"), silent = TRUE)
  try(haven::write_sav(microorganisms.old, "data-raw/microorganisms.old.sav"), silent = TRUE)
  try(haven::write_dta(microorganisms.old, "data-raw/microorganisms.old.dta"), silent = TRUE)
  try(openxlsx::write.xlsx(microorganisms.old, "data-raw/microorganisms.old.xlsx"), silent = TRUE)
}

ab <- dplyr::mutate_if(antibiotics, ~!is.numeric(.), as.character)
if (changed_md5(ab)) {
  write_md5(ab)
  try(saveRDS(ab, "data-raw/antibiotics.rds", version = 2, compress = "xz"), silent = TRUE)
  try(write.table(ab, "data-raw/antibiotics.txt", sep = "\t", na = "", row.names = FALSE), silent = TRUE)
  try(haven::write_sas(ab, "data-raw/antibiotics.sas"), silent = TRUE)
  try(haven::write_sav(ab, "data-raw/antibiotics.sav"), silent = TRUE)
  try(haven::write_dta(ab, "data-raw/antibiotics.dta"), silent = TRUE)
  try(openxlsx::write.xlsx(ab, "data-raw/antibiotics.xlsx"), silent = TRUE)
}

av <- dplyr::mutate_if(antivirals, ~!is.numeric(.), as.character)
if (changed_md5(av)) {
  write_md5(av)
  try(saveRDS(av, "data-raw/antivirals.rds", version = 2, compress = "xz"), silent = TRUE)
  try(write.table(av, "data-raw/antivirals.txt", sep = "\t", na = "", row.names = FALSE), silent = TRUE)
  try(haven::write_sas(av, "data-raw/antivirals.sas"), silent = TRUE)
  try(haven::write_sav(av, "data-raw/antivirals.sav"), silent = TRUE)
  try(haven::write_dta(av, "data-raw/antivirals.dta"), silent = TRUE)
  try(openxlsx::write.xlsx(av, "data-raw/antivirals.xlsx"), silent = TRUE)
}

if (changed_md5(intrinsic_resistant)) {
  write_md5(intrinsic_resistant)
  try(saveRDS(intrinsic_resistant, "data-raw/intrinsic_resistant.rds", version = 2, compress = "xz"), silent = TRUE)
  try(write.table(intrinsic_resistant, "data-raw/intrinsic_resistant.txt", sep = "\t", na = "", row.names = FALSE), silent = TRUE)
  try(haven::write_sas(intrinsic_resistant, "data-raw/intrinsic_resistant.sas"), silent = TRUE)
  try(haven::write_sav(intrinsic_resistant, "data-raw/intrinsic_resistant.sav"), silent = TRUE)
  try(haven::write_dta(intrinsic_resistant, "data-raw/intrinsic_resistant.dta"), silent = TRUE)
  try(openxlsx::write.xlsx(intrinsic_resistant, "data-raw/intrinsic_resistant.xlsx"), silent = TRUE)
}

rm(write_md5)
rm(changed_md5)
