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

# Run this file to update the package using:
# source("data-raw/_internals.R")

library(dplyr, warn.conflicts = FALSE)
devtools::load_all(quiet = TRUE)

old_globalenv <- ls(envir = globalenv())

# Helper functions --------------------------------------------------------

create_species_cons_cops <- function(type = c("CoNS", "CoPS")) {
  # Determination of which staphylococcal species are CoNS/CoPS according to:
  # - Becker et al. 2014, PMID 25278577
  # - Becker et al. 2019, PMID 30872103
  # - Becker et al. 2020, PMID 32056452
  # this function returns class <mo>
  MO_staph <- AMR::microorganisms
  MO_staph <- MO_staph[which(MO_staph$genus == "Staphylococcus"), , drop = FALSE]
  if (type == "CoNS") {
    MO_staph[which(MO_staph$species %in% c("coagulase-negative", "argensis", "arlettae",
                                           "auricularis", "borealis", "caeli", "capitis", "caprae", 
                                           "carnosus", "casei", "chromogenes", "cohnii", "condimenti",
                                           "croceilyticus",
                                           "debuckii", "devriesei", "edaphicus", "epidermidis",
                                           "equorum", "felis", "fleurettii", "gallinarum",
                                           "haemolyticus", "hominis", "jettensis", "kloosii",
                                           "lentus", "lugdunensis", "massiliensis", "microti",
                                           "muscae", "nepalensis", "pasteuri", "petrasii",
                                           "pettenkoferi", "piscifermentans", "pragensis", "pseudoxylosus",
                                           "pulvereri", "rostri", "saccharolyticus", "saprophyticus",
                                           "sciuri", "simulans", "stepanovicii", "succinus",
                                           "ureilyticus",
                                           "vitulinus", "vitulus", "warneri", "xylosus")
                   | (MO_staph$species == "schleiferi" & MO_staph$subspecies %in% c("schleiferi", ""))),
             "mo", drop = TRUE]
  } else if (type == "CoPS") {
    MO_staph[which(MO_staph$species %in% c("coagulase-positive", "coagulans",
                                           "agnetis", "argenteus",
                                           "cornubiensis",
                                           "delphini", "lutrae",
                                           "hyicus", "intermedius",
                                           "pseudintermedius", "pseudointermedius",
                                           "schweitzeri", "simiae")
                   | (MO_staph$species == "schleiferi" & MO_staph$subspecies == "coagulans")),
             "mo", drop = TRUE]
  }
}

create_AB_lookup <- function() {
  AB_lookup <- AMR::antibiotics
  AB_lookup$generalised_name <- generalise_antibiotic_name(AB_lookup$name)
  AB_lookup$generalised_synonyms <- lapply(AB_lookup$synonyms, generalise_antibiotic_name)
  AB_lookup$generalised_abbreviations <- lapply(AB_lookup$abbreviations, generalise_antibiotic_name)
  AB_lookup$generalised_loinc <- lapply(AB_lookup$loinc, generalise_antibiotic_name)
  AB_lookup$generalised_all <- unname(lapply(as.list(as.data.frame(t(AB_lookup[, 
                                                                               c("ab", "atc", "cid", "name",
                                                                                 colnames(AB_lookup)[colnames(AB_lookup) %like% "generalised"]),
                                                                               drop = FALSE]),
                                                                   stringsAsFactors = FALSE)),
                                             function(x) {
                                               x <- generalise_antibiotic_name(unname(unlist(x)))
                                               x[x != ""]
                                             }))
  AB_lookup
}

create_MO_lookup <- function() {
  MO_lookup <- AMR::microorganisms
  
  MO_lookup$kingdom_index <- NA_real_
  MO_lookup[which(MO_lookup$kingdom == "Bacteria" | MO_lookup$mo == "UNKNOWN"), "kingdom_index"] <- 1
  MO_lookup[which(MO_lookup$kingdom == "Fungi"), "kingdom_index"] <- 2
  MO_lookup[which(MO_lookup$kingdom == "Protozoa"), "kingdom_index"] <- 3
  MO_lookup[which(MO_lookup$kingdom == "Archaea"), "kingdom_index"] <- 4
  # all the rest
  MO_lookup[which(is.na(MO_lookup$kingdom_index)), "kingdom_index"] <- 5
  
  # use this paste instead of `fullname` to work with Viridans Group Streptococci, etc.
  MO_lookup$fullname_lower <- tolower(trimws(paste(MO_lookup$genus, 
                                                   MO_lookup$species,
                                                   MO_lookup$subspecies)))
  ind <- MO_lookup$genus == "" | grepl("^[(]unknown ", MO_lookup$fullname)
  MO_lookup[ind, "fullname_lower"] <- tolower(MO_lookup[ind, "fullname"])
  MO_lookup$fullname_lower <- trimws(gsub("[^.a-z0-9/ \\-]+", "", MO_lookup$fullname_lower, perl = TRUE))
  
  # add a column with only "e coli" like combinations
  MO_lookup$g_species <- gsub("^([a-z])[a-z]+ ([a-z]+) ?.*", "\\1 \\2", MO_lookup$fullname_lower, perl = TRUE)
  
  # so arrange data on prevalence first, then kingdom, then full name
  MO_lookup[order(MO_lookup$prevalence, MO_lookup$kingdom_index, MO_lookup$fullname_lower), ]
}

create_MO.old_lookup <- function() {
  MO.old_lookup <- AMR::microorganisms.old
  MO.old_lookup$fullname_lower <- trimws(gsub("[^.a-z0-9/ \\-]+", "", tolower(trimws(MO.old_lookup$fullname))))
  
  # add a column with only "e coli"-like combinations
  MO.old_lookup$g_species <- trimws(gsub("^([a-z])[a-z]+ ([a-z]+) ?.*", "\\1 \\2", MO.old_lookup$fullname_lower))
  
  # so arrange data on prevalence first, then full name
  MO.old_lookup[order(MO.old_lookup$prevalence, MO.old_lookup$fullname_lower), ]
}

create_intr_resistance <- function() {
  # for mo_is_intrinsic_resistant() - saves a lot of time when executed on this vector
  paste(AMR::microorganisms[match(AMR::intrinsic_resistant$microorganism, AMR::microorganisms$fullname), "mo", drop = TRUE],
        AMR::antibiotics[match(AMR::intrinsic_resistant$antibiotic, AMR::antibiotics$name), "ab", drop = TRUE])
}



# Save internal data sets to R/sysdata.rda --------------------------------

# See 'data-raw/eucast_rules.tsv' for the EUCAST reference file
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

# Translations
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

# Old microorganism codes
microorganisms.translation <- readRDS("data-raw/microorganisms.translation.rds")

# for mo_is_intrinsic_resistant() - saves a lot of time when executed on this vector
INTRINSIC_R <- create_intr_resistance()

# for checking input in `language` argument in e.g. mo_*() and ab_*() functions
LANGUAGES_SUPPORTED <- sort(c("en", colnames(translations_file)[nchar(colnames(translations_file)) == 2]))

# vectors of CoNS and CoPS, improves speed in as.mo()
MO_CONS <- create_species_cons_cops("CoNS")
MO_COPS <- create_species_cons_cops("CoPS")

# reference data - they have additional columns compared to `antibiotics` and `microorganisms` to improve speed
AB_lookup <- create_AB_lookup()
MO_lookup <- create_MO_lookup()
MO.old_lookup <- create_MO.old_lookup()

# Export to package as internal data ----
usethis::use_data(eucast_rules_file, 
                  translations_file,
                  microorganisms.translation,
                  INTRINSIC_R,
                  LANGUAGES_SUPPORTED,
                  MO_CONS,
                  MO_COPS,
                  AB_lookup,
                  MO_lookup,
                  MO.old_lookup,
                  internal = TRUE,
                  overwrite = TRUE,
                  version = 2,
                  compress = "xz")

# Export data sets to the repository in different formats -----------------

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

# give official names to ABs and MOs
rsi <- dplyr::mutate(rsi_translation, ab = ab_name(ab), mo = mo_name(mo))
if (changed_md5(rsi)) {
  usethis::ui_info(paste0("Saving {usethis::ui_value('rsi_translation')} to {usethis::ui_value('/data-raw/')}"))
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
  usethis::ui_info(paste0("Saving {usethis::ui_value('microorganisms')} to {usethis::ui_value('/data-raw/')}"))
  write_md5(mo)
  try(saveRDS(mo, "data-raw/microorganisms.rds", version = 2, compress = "xz"), silent = TRUE)
  try(write.table(mo, "data-raw/microorganisms.txt", sep = "\t", na = "", row.names = FALSE), silent = TRUE)
  try(haven::write_sas(dplyr::select(mo, -snomed), "data-raw/microorganisms.sas"), silent = TRUE)
  try(haven::write_sav(dplyr::select(mo, -snomed), "data-raw/microorganisms.sav"), silent = TRUE)
  try(haven::write_dta(dplyr::select(mo, -snomed), "data-raw/microorganisms.dta"), silent = TRUE)
  try(openxlsx::write.xlsx(mo, "data-raw/microorganisms.xlsx"), silent = TRUE)
}

if (changed_md5(microorganisms.old)) {
  usethis::ui_info(paste0("Saving {usethis::ui_value('microorganisms.old')} to {usethis::ui_value('/data-raw/')}"))
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
  usethis::ui_info(paste0("Saving {usethis::ui_value('antibiotics')} to {usethis::ui_value('/data-raw/')}"))
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
  usethis::ui_info(paste0("Saving {usethis::ui_value('antivirals')} to {usethis::ui_value('/data-raw/')}"))
  write_md5(av)
  try(saveRDS(av, "data-raw/antivirals.rds", version = 2, compress = "xz"), silent = TRUE)
  try(write.table(av, "data-raw/antivirals.txt", sep = "\t", na = "", row.names = FALSE), silent = TRUE)
  try(haven::write_sas(av, "data-raw/antivirals.sas"), silent = TRUE)
  try(haven::write_sav(av, "data-raw/antivirals.sav"), silent = TRUE)
  try(haven::write_dta(av, "data-raw/antivirals.dta"), silent = TRUE)
  try(openxlsx::write.xlsx(av, "data-raw/antivirals.xlsx"), silent = TRUE)
}

if (changed_md5(intrinsic_resistant)) {
  usethis::ui_info(paste0("Saving {usethis::ui_value('intrinsic_resistant')} to {usethis::ui_value('/data-raw/')}"))
  write_md5(intrinsic_resistant)
  try(saveRDS(intrinsic_resistant, "data-raw/intrinsic_resistant.rds", version = 2, compress = "xz"), silent = TRUE)
  try(write.table(intrinsic_resistant, "data-raw/intrinsic_resistant.txt", sep = "\t", na = "", row.names = FALSE), silent = TRUE)
  try(haven::write_sas(intrinsic_resistant, "data-raw/intrinsic_resistant.sas"), silent = TRUE)
  try(haven::write_sav(intrinsic_resistant, "data-raw/intrinsic_resistant.sav"), silent = TRUE)
  try(haven::write_dta(intrinsic_resistant, "data-raw/intrinsic_resistant.dta"), silent = TRUE)
  try(openxlsx::write.xlsx(intrinsic_resistant, "data-raw/intrinsic_resistant.xlsx"), silent = TRUE)
}

if (changed_md5(dosage)) {
  usethis::ui_info(paste0("Saving {usethis::ui_value('dosage')} to {usethis::ui_value('/data-raw/')}"))
  write_md5(dosage)
  try(saveRDS(dosage, "data-raw/dosage.rds", version = 2, compress = "xz"), silent = TRUE)
  try(write.table(dosage, "data-raw/dosage.txt", sep = "\t", na = "", row.names = FALSE), silent = TRUE)
  try(haven::write_sas(dosage, "data-raw/dosage.sas"), silent = TRUE)
  try(haven::write_sav(dosage, "data-raw/dosage.sav"), silent = TRUE)
  try(haven::write_dta(dosage, "data-raw/dosage.dta"), silent = TRUE)
  try(openxlsx::write.xlsx(dosage, "data-raw/dosage.xlsx"), silent = TRUE)
}

# remove leftovers from global env
current_globalenv <- ls(envir = globalenv())
rm(list = current_globalenv[!current_globalenv %in% old_globalenv])
rm(current_globalenv)
