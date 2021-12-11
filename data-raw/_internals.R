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

# Save internal data to R/sysdata.rda -------------------------------------

# See 'data-raw/eucast_rules.tsv' for the EUCAST reference file
EUCAST_RULES_DF <- utils::read.delim(file = "data-raw/eucast_rules.tsv",
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
TRANSLATIONS <- utils::read.delim(file = "data-raw/translations.tsv",
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

# for checking input in `language` argument in e.g. mo_*() and ab_*() functions
LANGUAGES_SUPPORTED <- sort(c("en", colnames(TRANSLATIONS)[nchar(colnames(TRANSLATIONS)) == 2]))

# EXAMPLE_ISOLATES <- readRDS("data-raw/example_isolates.rds")

# vectors of CoNS and CoPS, improves speed in as.mo()
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
                                           "vitulinus", "vitulus", "warneri", "xylosus",
                                           "caledonicus", "canis",
                                           "durrellii", "lloydii")
                   | (MO_staph$species == "schleiferi" & MO_staph$subspecies %in% c("schleiferi", ""))),
             "mo", drop = TRUE]
  } else if (type == "CoPS") {
    MO_staph[which(MO_staph$species %in% c("coagulase-positive", "coagulans",
                                           "agnetis", "argenteus",
                                           "cornubiensis",
                                           "delphini", "lutrae",
                                           "hyicus", "intermedius",
                                           "pseudintermedius", "pseudointermedius",
                                           "schweitzeri", "simiae",
                                           "roterodami")
                   | (MO_staph$species == "schleiferi" & MO_staph$subspecies == "coagulans")),
             "mo", drop = TRUE]
  }
}
MO_CONS <- create_species_cons_cops("CoNS")
MO_COPS <- create_species_cons_cops("CoPS")
MO_STREP_ABCG <- as.mo(MO_lookup[which(MO_lookup$genus == "Streptococcus"), "mo", drop = TRUE], Lancefield = TRUE) %in% c("B_STRPT_GRPA", "B_STRPT_GRPB", "B_STRPT_GRPC", "B_STRPT_GRPG")

# antibiotic groups
# (these will also be used for eucast_rules() and understanding data-raw/eucast_rules.tsv)
globalenv_before_ab <- c(ls(envir = globalenv()), "globalenv_before_ab")
AB_AMINOGLYCOSIDES <- antibiotics %>% filter(group %like% "aminoglycoside") %>% pull(ab)
AB_AMINOPENICILLINS <- as.ab(c("AMP", "AMX"))
AB_ANTIFUNGALS <- AB_lookup %>% filter(group %like% "antifungal") %>% pull(ab)
AB_ANTIMYCOBACTERIALS <- AB_lookup %>% filter(group %like% "antimycobacterial") %>% pull(ab)
AB_CARBAPENEMS <- antibiotics %>% filter(group %like% "carbapenem") %>% pull(ab)
AB_CEPHALOSPORINS <- antibiotics %>% filter(group %like% "cephalosporin") %>% pull(ab)
AB_CEPHALOSPORINS_1ST <- antibiotics %>% filter(group %like% "cephalosporin.*1") %>% pull(ab)
AB_CEPHALOSPORINS_2ND <- antibiotics %>% filter(group %like% "cephalosporin.*2") %>% pull(ab)
AB_CEPHALOSPORINS_3RD <- antibiotics %>% filter(group %like% "cephalosporin.*3") %>% pull(ab)
AB_CEPHALOSPORINS_4TH <- antibiotics %>% filter(group %like% "cephalosporin.*4") %>% pull(ab)
AB_CEPHALOSPORINS_5TH <- antibiotics %>% filter(group %like% "cephalosporin.*5") %>% pull(ab)
AB_CEPHALOSPORINS_EXCEPT_CAZ <- AB_CEPHALOSPORINS[AB_CEPHALOSPORINS != "CAZ"]
AB_FLUOROQUINOLONES <- antibiotics %>% filter(atc_group2 %like% "fluoroquinolone" | (group %like% "quinolone" & is.na(atc_group2))) %>% pull(ab)
AB_GLYCOPEPTIDES <- antibiotics %>% filter(group %like% "glycopeptide") %>% pull(ab)
AB_LIPOGLYCOPEPTIDES <- as.ab(c("DAL", "ORI", "TLV")) # dalba/orita/tela
AB_GLYCOPEPTIDES_EXCEPT_LIPO <- AB_GLYCOPEPTIDES[!AB_GLYCOPEPTIDES %in% AB_LIPOGLYCOPEPTIDES]
AB_LINCOSAMIDES <- antibiotics %>% filter(atc_group2 %like% "lincosamide" | (group %like% "lincosamide" & is.na(atc_group2))) %>% pull(ab)
AB_MACROLIDES <- antibiotics %>% filter(atc_group2 %like% "macrolide" | (group %like% "macrolide" & is.na(atc_group2))) %>% pull(ab)
AB_OXAZOLIDINONES <- antibiotics %>% filter(group %like% "oxazolidinone") %>% pull(ab)
AB_PENICILLINS <- antibiotics %>% filter(group %like% "penicillin") %>% pull(ab)
AB_POLYMYXINS <- antibiotics %>% filter(group %like% "polymyxin") %>% pull(ab)
AB_QUINOLONES <- antibiotics %>% filter(group %like% "quinolone") %>% pull(ab)
AB_STREPTOGRAMINS <- antibiotics %>% filter(atc_group2 %like% "streptogramin") %>% pull(ab)
AB_TETRACYCLINES <- antibiotics %>% filter(group %like% "tetracycline") %>% pull(ab)
AB_TETRACYCLINES_EXCEPT_TGC <- AB_TETRACYCLINES[AB_TETRACYCLINES != "TGC"]
AB_TRIMETHOPRIMS <- antibiotics %>% filter(group %like% "trimethoprim") %>% pull(ab)
AB_UREIDOPENICILLINS <- as.ab(c("PIP", "TZP", "AZL", "MEZ"))
AB_BETALACTAMS <- c(AB_PENICILLINS, AB_CEPHALOSPORINS, AB_CARBAPENEMS)
# this will be used for documentation:
DEFINED_AB_GROUPS <- ls(envir = globalenv())
DEFINED_AB_GROUPS <- DEFINED_AB_GROUPS[!DEFINED_AB_GROUPS %in% globalenv_before_ab]

# Export to package as internal data ----
usethis::use_data(EUCAST_RULES_DF, 
                  TRANSLATIONS,
                  LANGUAGES_SUPPORTED,
                  # EXAMPLE_ISOLATES,
                  MO_CONS,
                  MO_COPS,
                  MO_STREP_ABCG,
                  AB_AMINOGLYCOSIDES,
                  AB_AMINOPENICILLINS,
                  AB_ANTIFUNGALS,
                  AB_ANTIMYCOBACTERIALS,
                  AB_CARBAPENEMS,
                  AB_CEPHALOSPORINS,
                  AB_CEPHALOSPORINS_1ST,
                  AB_CEPHALOSPORINS_2ND,
                  AB_CEPHALOSPORINS_3RD,
                  AB_CEPHALOSPORINS_4TH,
                  AB_CEPHALOSPORINS_5TH,
                  AB_CEPHALOSPORINS_EXCEPT_CAZ,
                  AB_FLUOROQUINOLONES,
                  AB_LIPOGLYCOPEPTIDES,
                  AB_GLYCOPEPTIDES,
                  AB_GLYCOPEPTIDES_EXCEPT_LIPO,
                  AB_LINCOSAMIDES,
                  AB_MACROLIDES,
                  AB_OXAZOLIDINONES,
                  AB_PENICILLINS,
                  AB_POLYMYXINS,
                  AB_QUINOLONES,
                  AB_STREPTOGRAMINS,
                  AB_TETRACYCLINES,
                  AB_TETRACYCLINES_EXCEPT_TGC,
                  AB_TRIMETHOPRIMS,
                  AB_UREIDOPENICILLINS,
                  AB_BETALACTAMS,
                  DEFINED_AB_GROUPS,
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
