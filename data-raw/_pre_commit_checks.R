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
# Center Groningen in The Netherlands, in collaboration with many      #
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
# how to conduct AMR data analysis: https://msberends.github.io/AMR/   #
# ==================================================================== #

# Run this file to update the package using:
# source("data-raw/_pre_commit_checks.R")

library(dplyr, warn.conflicts = FALSE)
try(detach("package:data.table", unload = TRUE), silent = TRUE) # to prevent like() to precede over AMR::like
devtools::load_all(quiet = TRUE)

suppressMessages(set_AMR_locale("English"))

pre_commit_lst <- list()

# Save internal data to R/sysdata.rda -------------------------------------

usethis::ui_info(paste0("Updating internal package data"))

# See 'data-raw/eucast_rules.tsv' for the EUCAST reference file
pre_commit_lst$EUCAST_RULES_DF <- utils::read.delim(
  file = "data-raw/eucast_rules.tsv",
  skip = 9,
  sep = "\t",
  stringsAsFactors = FALSE,
  header = TRUE,
  strip.white = TRUE,
  na = c(NA, "", NULL)
) %>%
  # take the order of the reference.rule_group column in the original data file
  mutate(
    reference.rule_group = factor(reference.rule_group,
      levels = unique(reference.rule_group),
      ordered = TRUE
    ),
    sorting_rule = ifelse(grepl("^Table", reference.rule, ignore.case = TRUE), 1, 2)
  ) %>%
  arrange(
    reference.rule_group,
    reference.version,
    sorting_rule,
    reference.rule
  ) %>%
  mutate(reference.rule_group = as.character(reference.rule_group)) %>%
  select(-sorting_rule)

pre_commit_lst$TRANSLATIONS <- utils::read.delim(
  file = "data-raw/translations.tsv",
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
  quote = ""
)

pre_commit_lst$LANGUAGES_SUPPORTED_NAMES <- c(
  list(en = list(exonym = "English", endonym = "English")),
  lapply(
    TRANSLATIONS[, which(nchar(colnames(pre_commit_lst$TRANSLATIONS)) == 2), drop = FALSE],
    function(x) list(exonym = x[1], endonym = x[2])
  )
)

pre_commit_lst$LANGUAGES_SUPPORTED <- names(pre_commit_lst$LANGUAGES_SUPPORTED_NAMES)

# vectors of CoNS and CoPS, improves speed in as.mo()
create_species_cons_cops <- function(type = c("CoNS", "CoPS")) {
  # Determination of which staphylococcal species are CoNS/CoPS according to:
  # - Becker et al. 2014, PMID 25278577
  # - Becker et al. 2019, PMID 30872103
  # - Becker et al. 2020, PMID 32056452
  # this function returns class <mo>
  MO_staph <- microorganisms
  MO_staph <- MO_staph[which(MO_staph$genus == "Staphylococcus"), , drop = FALSE]
  if (type == "CoNS") {
    MO_staph[
      which(MO_staph$species %in% c(
        "coagulase-negative", "argensis", "arlettae",
        "auricularis", "borealis", "caeli", "capitis", "caprae",
        "carnosus", "casei", "caseolyticus", "chromogenes", "cohnii", "condimenti",
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
        "durrellii", "lloydii",
        "ratti", "taiwanensis", "veratri", "urealyticus",
        "americanisciuri", "marylandisciuri", "shinii", "brunensis"
      ) |
        # old, now renamed to S. schleiferi (but still as synonym in our data of course):
        (MO_staph$species == "schleiferi" & MO_staph$subspecies %in% c("schleiferi", ""))),
      "mo",
      drop = TRUE
    ]
  } else if (type == "CoPS") {
    MO_staph[
      which(MO_staph$species %in% c(
        "coagulase-positive", "coagulans",
        "agnetis", "argenteus",
        "cornubiensis",
        "delphini", "lutrae",
        "hyicus", "intermedius",
        "pseudintermedius", "pseudointermedius",
        "schweitzeri", "simiae",
        "roterodami",
        "singaporensis"
      ) |
        # old, now renamed to S. coagulans (but still as synonym in our data of course):
        (MO_staph$species == "schleiferi" & MO_staph$subspecies == "coagulans")),
      "mo",
      drop = TRUE
    ]
  }
}
pre_commit_lst$MO_CONS <- create_species_cons_cops("CoNS")
pre_commit_lst$MO_COPS <- create_species_cons_cops("CoPS")
pre_commit_lst$MO_STREP_ABCG <- microorganisms$mo[which(microorganisms$genus == "Streptococcus" &
  tolower(microorganisms$species) %in% c(
    "pyogenes", "agalactiae", "dysgalactiae", "equi", "canis",
    "group a", "group b", "group c", "group g"
  ))]
pre_commit_lst$MO_LANCEFIELD <- microorganisms$mo[which(microorganisms$mo %like% "^(B_STRPT_PYGN(_|$)|B_STRPT_AGLC(_|$)|B_STRPT_(DYSG|EQUI)(_|$)|B_STRPT_ANGN(_|$)|B_STRPT_(DYSG|CANS)(_|$)|B_STRPT_SNGN(_|$)|B_STRPT_SLVR(_|$))")]
pre_commit_lst$MO_WHO_PRIORITY_GENERA <- c(
  # World Health Organization's (WHO) Priority Pathogen List (some are from the group Enterobacteriaceae)
  "Acinetobacter",
  "Aspergillus",
  "Blastomyces",
  "Campylobacter",
  "Candida",
  "Citrobacter",
  "Clostridioides",
  "Coccidioides",
  "Cryptococcus",
  "Edwardsiella",
  "Enterobacter",
  "Enterococcus",
  "Escherichia",
  "Fusarium",
  "Haemophilus",
  "Helicobacter",
  "Histoplasma",
  "Klebsiella",
  "Morganella",
  "Mycobacterium",
  "Neisseria",
  "Paracoccidioides",
  "Pneumocystis",
  "Proteus",
  "Providencia",
  "Pseudomonas",
  "Salmonella",
  "Serratia",
  "Shigella",
  "Staphylococcus",
  "Streptococcus",
  "Yersinia"
)
pre_commit_lst$MO_RELEVANT_GENERA <- c(
  "Absidia",
  "Acanthamoeba",
  "Acremonium",
  "Actinomucor",
  "Aedes",
  "Alternaria",
  "Amoeba",
  "Ancylostoma",
  "Angiostrongylus",
  "Anisakis",
  "Anopheles",
  "Apophysomyces",
  "Arthroderma",
  "Aspergillus",
  "Aureobasidium",
  "Basidiobolus",
  "Beauveria",
  "Bipolaris",
  "Blastobotrys",
  "Blastocystis",
  "Blastomyces",
  "Candida",
  "Capillaria",
  "Chaetomium",
  "Chilomastix",
  "Chrysonilia",
  "Chrysosporium",
  "Cladophialophora",
  "Cladosporium",
  "Clavispora",
  "Coccidioides",
  "Cokeromyces",
  "Conidiobolus",
  "Coniochaeta",
  "Contracaecum",
  "Cordylobia",
  "Cryptococcus",
  "Cryptosporidium",
  "Cunninghamella",
  "Curvularia",
  "Cyberlindnera",
  "Debaryozyma",
  "Demodex",
  "Dermatobia",
  "Dientamoeba",
  "Diphyllobothrium",
  "Dirofilaria",
  "Echinostoma",
  "Entamoeba",
  "Enterobius",
  "Epidermophyton",
  "Exidia",
  "Exophiala",
  "Exserohilum",
  "Fasciola",
  "Fonsecaea",
  "Fusarium",
  "Geotrichum",
  "Giardia",
  "Graphium",
  "Haloarcula",
  "Halobacterium",
  "Halococcus",
  "Hansenula",
  "Hendersonula",
  "Heterophyes",
  "Histomonas",
  "Histoplasma",
  "Hortaea",
  "Hymenolepis",
  "Hypomyces",
  "Hysterothylacium",
  "Kloeckera",
  "Kluyveromyces",
  "Kodamaea",
  "Lacazia",
  "Leishmania",
  "Lichtheimia",
  "Lodderomyces",
  "Lomentospora",
  "Madurella",
  "Malassezia",
  "Malbranchea",
  "Metagonimus",
  "Meyerozyma",
  "Microsporidium",
  "Microsporum",
  "Millerozyma",
  "Mortierella",
  "Mucor",
  "Mycocentrospora",
  "Nannizzia",
  "Necator",
  "Nectria",
  "Ochroconis",
  "Oesophagostomum",
  "Oidiodendron",
  "Opisthorchis",
  "Paecilomyces",
  "Paracoccidioides",
  "Pediculus",
  "Penicillium",
  "Phaeoacremonium",
  "Phaeomoniella",
  "Phialophora",
  "Phlebotomus",
  "Phoma",
  "Pichia",
  "Piedraia",
  "Pithomyces",
  "Pityrosporum",
  "Pneumocystis",
  "Pseudallescheria",
  "Pseudoscopulariopsis",
  "Pseudoterranova",
  "Pulex",
  "Purpureocillium",
  "Quambalaria",
  "Rhinocladiella",
  "Rhizomucor",
  "Rhizopus",
  "Rhodotorula",
  "Saccharomyces",
  "Saksenaea",
  "Saprochaete",
  "Sarcoptes",
  "Scedosporium",
  "Schistosoma",
  "Schizosaccharomyces",
  "Scolecobasidium",
  "Scopulariopsis",
  "Scytalidium",
  "Spirometra",
  "Sporobolomyces",
  "Sporopachydermia",
  "Sporothrix",
  "Sporotrichum",
  "Stachybotrys",
  "Strongyloides",
  "Syncephalastrum",
  "Syngamus",
  "Taenia",
  "Talaromyces",
  "Teleomorph",
  "Toxocara",
  "Trichinella",
  "Trichobilharzia",
  "Trichoderma",
  "Trichomonas",
  "Trichophyton",
  "Trichosporon",
  "Trichostrongylus",
  "Trichuris",
  "Tritirachium",
  "Trombicula",
  "Trypanosoma",
  "Tunga",
  "Ulocladium",
  "Ustilago",
  "Verticillium",
  "Wallemia",
  "Wangiella",
  "Wickerhamomyces",
  "Wuchereria",
  "Yarrowia",
  "Zygosaccharomyces"
)

# antibiotic groups
# (these will also be used for eucast_rules() and understanding data-raw/eucast_rules.tsv)
pre_commit_lst$AB_AMINOGLYCOSIDES <- antibiotics %>%
  filter(group %like% "aminoglycoside") %>%
  pull(ab)
pre_commit_lst$AB_AMINOPENICILLINS <- as.ab(c("AMP", "AMX"))
pre_commit_lst$AB_ANTIFUNGALS <- antibiotics %>%
  filter(group %like% "antifungal") %>%
  pull(ab)
pre_commit_lst$AB_ANTIMYCOBACTERIALS <- antibiotics %>%
  filter(group %like% "antimycobacterial") %>%
  pull(ab)
pre_commit_lst$AB_CARBAPENEMS <- antibiotics %>%
  filter(group %like% "carbapenem") %>%
  pull(ab)
pre_commit_lst$AB_CEPHALOSPORINS <- antibiotics %>%
  filter(group %like% "cephalosporin") %>%
  pull(ab)
pre_commit_lst$AB_CEPHALOSPORINS_1ST <- antibiotics %>%
  filter(group %like% "cephalosporin.*1") %>%
  pull(ab)
pre_commit_lst$AB_CEPHALOSPORINS_2ND <- antibiotics %>%
  filter(group %like% "cephalosporin.*2") %>%
  pull(ab)
pre_commit_lst$AB_CEPHALOSPORINS_3RD <- antibiotics %>%
  filter(group %like% "cephalosporin.*3") %>%
  pull(ab)
pre_commit_lst$AB_CEPHALOSPORINS_4TH <- antibiotics %>%
  filter(group %like% "cephalosporin.*4") %>%
  pull(ab)
pre_commit_lst$AB_CEPHALOSPORINS_5TH <- antibiotics %>%
  filter(group %like% "cephalosporin.*5") %>%
  pull(ab)
pre_commit_lst$AB_CEPHALOSPORINS_EXCEPT_CAZ <- pre_commit_lst$AB_CEPHALOSPORINS[pre_commit_lst$AB_CEPHALOSPORINS != "CAZ"]
pre_commit_lst$AB_FLUOROQUINOLONES <- antibiotics %>%
  filter(atc_group2 %like% "fluoroquinolone" | (group %like% "quinolone" & is.na(atc_group2))) %>%
  pull(ab)
pre_commit_lst$AB_GLYCOPEPTIDES <- antibiotics %>%
  filter(group %like% "glycopeptide") %>%
  pull(ab)
pre_commit_lst$AB_ISOXAZOLYLPENICILLINS <- antibiotics %>%
  filter(name %like% "oxacillin|cloxacillin|dicloxacillin|flucloxacillin|meth?icillin") %>%
  pull(ab)
pre_commit_lst$AB_LIPOGLYCOPEPTIDES <- as.ab(c("DAL", "ORI", "TLV")) # dalba/orita/tela
pre_commit_lst$AB_GLYCOPEPTIDES_EXCEPT_LIPO <- pre_commit_lst$AB_GLYCOPEPTIDES[!pre_commit_lst$AB_GLYCOPEPTIDES %in% pre_commit_lst$AB_LIPOGLYCOPEPTIDES]
pre_commit_lst$AB_LINCOSAMIDES <- antibiotics %>%
  filter(atc_group2 %like% "lincosamide" | (group %like% "lincosamide" & is.na(atc_group2) & name %like% "^(pirlimycin)" & name %unlike% "screening|inducible")) %>%
  pull(ab)
pre_commit_lst$AB_MACROLIDES <- antibiotics %>%
  filter(atc_group2 %like% "macrolide" | (group %like% "macrolide" & is.na(atc_group2) & name %like% "^(acetylmidecamycin|acetylspiramycin|gamith?romycin|kitasamycin|meleumycin|nafith?romycin|solith?romycin|tildipirosin|tilmicosin|tulath?romycin|tylosin|tylvalosin)" & name %unlike% "screening|inducible")) %>%
  pull(ab)
pre_commit_lst$AB_MONOBACTAMS <- antibiotics %>%
  filter(group %like% "monobactam") %>%
  pull(ab)
pre_commit_lst$AB_NITROFURANS <- antibiotics %>%
  filter(name %like% "^furaz|nitrofura" | atc_group2 %like% "nitrofuran") %>%
  pull(ab)
pre_commit_lst$AB_OXAZOLIDINONES <- antibiotics %>%
  filter(group %like% "oxazolidinone") %>%
  pull(ab)
pre_commit_lst$AB_PENICILLINS <- antibiotics %>%
  filter(group %like% "penicillin" & !(name %unlike% "/" & name %like% ".*bactam$")) %>%
  pull(ab)
pre_commit_lst$AB_PHENICOLS <- antibiotics %>%
  filter(group %like% "phenicol" | atc_group1 %like% "phenicol" | atc_group2 %like% "phenicol") %>%
  pull(ab)
pre_commit_lst$AB_POLYMYXINS <- antibiotics %>%
  filter(group %like% "polymyxin") %>%
  pull(ab)
pre_commit_lst$AB_QUINOLONES <- antibiotics %>%
  filter(group %like% "quinolone") %>%
  pull(ab)
pre_commit_lst$AB_RIFAMYCINS <- antibiotics %>%
  filter(name %like% "Rifampi|Rifabutin|Rifapentine|rifamy") %>%
  pull(ab)
pre_commit_lst$AB_STREPTOGRAMINS <- antibiotics %>%
  filter(atc_group2 %like% "streptogramin") %>%
  pull(ab)
pre_commit_lst$AB_TETRACYCLINES <- antibiotics %>%
  filter(group %like% "tetracycline") %>%
  pull(ab)
pre_commit_lst$AB_TETRACYCLINES_EXCEPT_TGC <- pre_commit_lst$AB_TETRACYCLINES[pre_commit_lst$AB_TETRACYCLINES != "TGC"]
pre_commit_lst$AB_TRIMETHOPRIMS <- antibiotics %>%
  filter(group %like% "trimethoprim") %>%
  pull(ab)
pre_commit_lst$AB_UREIDOPENICILLINS <- as.ab(c("PIP", "TZP", "AZL", "MEZ"))
pre_commit_lst$AB_BETALACTAMS <- sort(c(pre_commit_lst$AB_PENICILLINS, pre_commit_lst$AB_CEPHALOSPORINS, pre_commit_lst$AB_CARBAPENEMS, pre_commit_lst$AB_MONOBACTAMS))
pre_commit_lst$AB_BETALACTAMS_WITH_INHIBITOR <- antibiotics %>%
  filter(name %like% "/" & name %unlike% "EDTA" & ab %in% pre_commit_lst$AB_BETALACTAMS) %>%
  pull(ab)
# this will be used for documentation:
pre_commit_lst$DEFINED_AB_GROUPS <- sort(names(pre_commit_lst)[names(pre_commit_lst) %like% "^AB_" & names(pre_commit_lst) != "AB_LOOKUP"])

pre_commit_lst$AB_LOOKUP <- create_AB_AV_lookup(antibiotics)
pre_commit_lst$AV_LOOKUP <- create_AB_AV_lookup(antivirals)

# Export to package as internal data ----
# usethis::use_data() must receive unquoted object names, which is not flexible at all.
# we'll use good old base::save() instead
save(list = names(pre_commit_lst),
     file = "R/sysdata.rda",
     envir = as.environment(pre_commit_lst),
     compress = "xz",
     version = 2,
     ascii = FALSE)
usethis::ui_done("Saved to {usethis::ui_value('R/sysdata.rda')}")




# Export data sets to the repository in different formats -----------------

for (pkg in c("haven", "openxlsx2", "arrow")) {
  if (!pkg %in% rownames(utils::installed.packages())) {
    message("NOTE: package '", pkg, "' not installed! Ignoring export where this package is required.")
  }
}
if ("digest" %in% rownames(utils::installed.packages())) {
  md5 <- function(object) digest::digest(object, "md5")
} else {
  # will write all files anyway, since MD5 hash cannot be determined
  md5 <- function(object) "unknown-md5-hash"
}

write_md5 <- function(object) {
  conn <- file(paste0("data-raw/", deparse(substitute(object)), ".md5"))
  writeLines(md5(object), conn)
  close(conn)
}
changed_md5 <- function(object) {
  tryCatch(
    {
      conn <- file(paste0("data-raw/", deparse(substitute(object)), ".md5"))
      compared <- md5(object) != readLines(con = conn)
      close(conn)
      compared
    },
    error = function(e) TRUE
  )
}

# give official names to ABs and MOs
clin_break <- clinical_breakpoints %>%
  mutate(mo_name = microorganisms$fullname[match(mo, microorganisms$mo)], .after = mo) %>%
  mutate(ab_name = antibiotics$name[match(ab, antibiotics$ab)], .after = ab)
if (changed_md5(clin_break)) {
  usethis::ui_info(paste0("Saving {usethis::ui_value('clinical_breakpoints')} to {usethis::ui_value('data-raw/')}"))
  write_md5(clin_break)
  try(saveRDS(clin_break, "data-raw/clinical_breakpoints.rds", version = 2, compress = "xz"), silent = TRUE)
  try(write.table(clinical_breakpoints, "data-raw/clinical_breakpoints.txt", sep = "\t", na = "", row.names = FALSE), silent = TRUE)
  try(haven::write_sav(clin_break, "data-raw/clinical_breakpoints.sav"), silent = TRUE)
  try(haven::write_dta(clin_break, "data-raw/clinical_breakpoints.dta"), silent = TRUE)
  try(openxlsx2::write_xlsx(clin_break, "data-raw/clinical_breakpoints.xlsx"), silent = TRUE)
  try(arrow::write_feather(clin_break, "data-raw/clinical_breakpoints.feather"), silent = TRUE)
  try(arrow::write_parquet(clin_break, "data-raw/clinical_breakpoints.parquet"), silent = TRUE)
}

if (changed_md5(microorganisms)) {
  usethis::ui_info(paste0("Saving {usethis::ui_value('microorganisms')} to {usethis::ui_value('data-raw/')}"))
  write_md5(microorganisms)
  try(saveRDS(microorganisms, "data-raw/microorganisms.rds", version = 2, compress = "xz"), silent = TRUE)
  max_50_snomed <- sapply(microorganisms$snomed, function(x) paste(x[seq_len(min(50, length(x), na.rm = TRUE))], collapse = " "))
  mo <- microorganisms
  mo$snomed <- max_50_snomed
  mo <- dplyr::mutate_if(mo, ~ !is.numeric(.), as.character)
  try(haven::write_sav(mo, "data-raw/microorganisms.sav"), silent = TRUE)
  try(haven::write_dta(mo, "data-raw/microorganisms.dta"), silent = TRUE)
  mo_all_snomed <- microorganisms %>% mutate_if(is.list, function(x) sapply(x, paste, collapse = ","))
  try(write.table(mo_all_snomed, "data-raw/microorganisms.txt", sep = "\t", na = "", row.names = FALSE), silent = TRUE)
  try(openxlsx2::write_xlsx(mo_all_snomed, "data-raw/microorganisms.xlsx"), silent = TRUE)
  try(arrow::write_feather(microorganisms, "data-raw/microorganisms.feather"), silent = TRUE)
  try(arrow::write_parquet(microorganisms, "data-raw/microorganisms.parquet"), silent = TRUE)
}

if (changed_md5(microorganisms.codes)) {
  usethis::ui_info(paste0("Saving {usethis::ui_value('microorganisms.codes')} to {usethis::ui_value('data-raw/')}"))
  write_md5(microorganisms.codes)
  try(saveRDS(microorganisms.codes, "data-raw/microorganisms.codes.rds", version = 2, compress = "xz"), silent = TRUE)
  try(write.table(microorganisms.codes, "data-raw/microorganisms.codes.txt", sep = "\t", na = "", row.names = FALSE), silent = TRUE)
  try(haven::write_sav(microorganisms.codes, "data-raw/microorganisms.codes.sav"), silent = TRUE)
  try(haven::write_dta(microorganisms.codes, "data-raw/microorganisms.codes.dta"), silent = TRUE)
  try(openxlsx2::write_xlsx(microorganisms.codes, "data-raw/microorganisms.codes.xlsx"), silent = TRUE)
  try(arrow::write_feather(microorganisms.codes, "data-raw/microorganisms.codes.feather"), silent = TRUE)
  try(arrow::write_parquet(microorganisms.codes, "data-raw/microorganisms.codes.parquet"), silent = TRUE)
}

if (changed_md5(microorganisms.groups)) {
  usethis::ui_info(paste0("Saving {usethis::ui_value('microorganisms.groups')} to {usethis::ui_value('data-raw/')}"))
  write_md5(microorganisms.groups)
  try(saveRDS(microorganisms.groups, "data-raw/microorganisms.groups.rds", version = 2, compress = "xz"), silent = TRUE)
  try(write.table(microorganisms.groups, "data-raw/microorganisms.groups.txt", sep = "\t", na = "", row.names = FALSE), silent = TRUE)
  try(haven::write_sav(microorganisms.groups, "data-raw/microorganisms.groups.sav"), silent = TRUE)
  try(haven::write_dta(microorganisms.groups, "data-raw/microorganisms.groups.dta"), silent = TRUE)
  try(openxlsx2::write_xlsx(microorganisms.groups, "data-raw/microorganisms.groups.xlsx"), silent = TRUE)
  try(arrow::write_feather(microorganisms.groups, "data-raw/microorganisms.groups.feather"), silent = TRUE)
  try(arrow::write_parquet(microorganisms.groups, "data-raw/microorganisms.groups.parquet"), silent = TRUE)
}

ab <- dplyr::mutate_if(antibiotics, ~ !is.numeric(.), as.character)
if (changed_md5(ab)) {
  usethis::ui_info(paste0("Saving {usethis::ui_value('antibiotics')} to {usethis::ui_value('data-raw/')}"))
  write_md5(ab)
  try(saveRDS(antibiotics, "data-raw/antibiotics.rds", version = 2, compress = "xz"), silent = TRUE)
  try(haven::write_sav(ab, "data-raw/antibiotics.sav"), silent = TRUE)
  try(haven::write_dta(ab, "data-raw/antibiotics.dta"), silent = TRUE)
  ab_lists <- antibiotics %>% mutate_if(is.list, function(x) sapply(x, paste, collapse = ","))
  try(write.table(ab_lists, "data-raw/antibiotics.txt", sep = "\t", na = "", row.names = FALSE), silent = TRUE)
  try(openxlsx2::write_xlsx(ab_lists, "data-raw/antibiotics.xlsx"), silent = TRUE)
  try(arrow::write_feather(antibiotics, "data-raw/antibiotics.feather"), silent = TRUE)
  try(arrow::write_parquet(antibiotics, "data-raw/antibiotics.parquet"), silent = TRUE)
}

av <- dplyr::mutate_if(antivirals, ~ !is.numeric(.), as.character)
if (changed_md5(av)) {
  usethis::ui_info(paste0("Saving {usethis::ui_value('antivirals')} to {usethis::ui_value('data-raw/')}"))
  write_md5(av)
  try(saveRDS(antivirals, "data-raw/antivirals.rds", version = 2, compress = "xz"), silent = TRUE)
  try(haven::write_sav(av, "data-raw/antivirals.sav"), silent = TRUE)
  try(haven::write_dta(av, "data-raw/antivirals.dta"), silent = TRUE)
  av_lists <- antivirals %>% mutate_if(is.list, function(x) sapply(x, paste, collapse = ","))
  try(write.table(av_lists, "data-raw/antivirals.txt", sep = "\t", na = "", row.names = FALSE), silent = TRUE)
  try(openxlsx2::write_xlsx(av_lists, "data-raw/antivirals.xlsx"), silent = TRUE)
  try(arrow::write_feather(antivirals, "data-raw/antivirals.feather"), silent = TRUE)
  try(arrow::write_parquet(antivirals, "data-raw/antivirals.parquet"), silent = TRUE)
}

# give official names to ABs and MOs
intrinsicR <- data.frame(
  microorganism = mo_name(intrinsic_resistant$mo, language = NULL, keep_synonyms = TRUE, info = FALSE),
  antibiotic = ab_name(intrinsic_resistant$ab, language = NULL),
  stringsAsFactors = FALSE
)
if (changed_md5(intrinsicR)) {
  usethis::ui_info(paste0("Saving {usethis::ui_value('intrinsic_resistant')} to {usethis::ui_value('data-raw/')}"))
  write_md5(intrinsicR)
  try(saveRDS(intrinsicR, "data-raw/intrinsic_resistant.rds", version = 2, compress = "xz"), silent = TRUE)
  try(write.table(intrinsicR, "data-raw/intrinsic_resistant.txt", sep = "\t", na = "", row.names = FALSE), silent = TRUE)
  try(haven::write_sav(intrinsicR, "data-raw/intrinsic_resistant.sav"), silent = TRUE)
  try(haven::write_dta(intrinsicR, "data-raw/intrinsic_resistant.dta"), silent = TRUE)
  try(openxlsx2::write_xlsx(intrinsicR, "data-raw/intrinsic_resistant.xlsx"), silent = TRUE)
  try(arrow::write_feather(intrinsicR, "data-raw/intrinsic_resistant.feather"), silent = TRUE)
  try(arrow::write_parquet(intrinsicR, "data-raw/intrinsic_resistant.parquet"), silent = TRUE)
}

if (changed_md5(dosage)) {
  usethis::ui_info(paste0("Saving {usethis::ui_value('dosage')} to {usethis::ui_value('data-raw/')}"))
  write_md5(dosage)
  try(saveRDS(dosage, "data-raw/dosage.rds", version = 2, compress = "xz"), silent = TRUE)
  try(write.table(dosage, "data-raw/dosage.txt", sep = "\t", na = "", row.names = FALSE), silent = TRUE)
  try(haven::write_sav(dosage, "data-raw/dosage.sav"), silent = TRUE)
  try(haven::write_dta(dosage, "data-raw/dosage.dta"), silent = TRUE)
  try(openxlsx2::write_xlsx(dosage, "data-raw/dosage.xlsx"), silent = TRUE)
  try(arrow::write_feather(dosage, "data-raw/dosage.feather"), silent = TRUE)
  try(arrow::write_parquet(dosage, "data-raw/dosage.parquet"), silent = TRUE)
}

suppressMessages(reset_AMR_locale())

devtools::load_all(quiet = TRUE)
suppressMessages(set_AMR_locale("English"))

# Update URLs -------------------------------------------------------------
usethis::ui_info("Checking URLs for redirects")
invisible(urlchecker::url_update("."))

# Style pkg ---------------------------------------------------------------
usethis::ui_info("Styling package")
styler::style_pkg(include_roxygen_examples = FALSE,
                  exclude_dirs = list.dirs(full.names = FALSE, recursive = FALSE)[!list.dirs(full.names = FALSE, recursive = FALSE) %in% c("R", "tests")])

# Document pkg ------------------------------------------------------------
usethis::ui_info("Documenting package")
suppressMessages(devtools::document(quiet = TRUE))

# Finished ----------------------------------------------------------------
usethis::ui_done("All done")
suppressMessages(reset_AMR_locale())
