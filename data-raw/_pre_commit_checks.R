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
# how to conduct AMR data analysis: https://amr-for-r.org              #
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
pre_commit_lst$TRANSLATIONS <- pre_commit_lst$TRANSLATIONS[, which(colnames(pre_commit_lst$TRANSLATIONS) != "en"), drop = FALSE]

pre_commit_lst$LANGUAGES_SUPPORTED_NAMES <- c(
  list(en = list(exonym = "English", endonym = "English")),
  lapply(
    pre_commit_lst$TRANSLATIONS[, which(nchar(colnames(pre_commit_lst$TRANSLATIONS)) == 2), drop = FALSE],
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
pre_commit_lst$AB_AMINOGLYCOSIDES <- antimicrobials %>%
  filter(group %like% "aminoglycoside") %>%
  pull(ab)
pre_commit_lst$AB_AMINOPENICILLINS <- as.ab(c("AMP", "AMX", "AMC"))
pre_commit_lst$AB_ANTIFUNGALS <- antimicrobials %>%
  filter(group %like% "antifungal") %>%
  pull(ab)
pre_commit_lst$AB_ANTIMYCOBACTERIALS <- antimicrobials %>%
  filter(group %like% "antimycobacterial") %>%
  pull(ab)
pre_commit_lst$AB_CARBAPENEMS <- antimicrobials %>%
  filter(group %like% "carbapenem") %>%
  pull(ab)
pre_commit_lst$AB_CEPHALOSPORINS <- antimicrobials %>%
  filter(group %like% "cephalosporin") %>%
  pull(ab)
pre_commit_lst$AB_CEPHALOSPORINS_1ST <- antimicrobials %>%
  filter(group %like% "cephalosporin.*1") %>%
  pull(ab)
pre_commit_lst$AB_CEPHALOSPORINS_2ND <- antimicrobials %>%
  filter(group %like% "cephalosporin.*2") %>%
  pull(ab)
pre_commit_lst$AB_CEPHALOSPORINS_3RD <- antimicrobials %>%
  filter(group %like% "cephalosporin.*3") %>%
  pull(ab)
pre_commit_lst$AB_CEPHALOSPORINS_4TH <- antimicrobials %>%
  filter(group %like% "cephalosporin.*4") %>%
  pull(ab)
pre_commit_lst$AB_CEPHALOSPORINS_5TH <- antimicrobials %>%
  filter(group %like% "cephalosporin.*5") %>%
  pull(ab)
pre_commit_lst$AB_CEPHALOSPORINS_EXCEPT_CAZ <- pre_commit_lst$AB_CEPHALOSPORINS[pre_commit_lst$AB_CEPHALOSPORINS != "CAZ"]
pre_commit_lst$AB_GLYCOPEPTIDES <- antimicrobials %>%
  filter(group %like% "glycopeptide") %>%
  pull(ab)
pre_commit_lst$AB_ISOXAZOLYLPENICILLINS <- antimicrobials %>%
  filter(name %like% "oxacillin|cloxacillin|dicloxacillin|flucloxacillin|meth?icillin") %>%
  pull(ab)
pre_commit_lst$AB_LIPOGLYCOPEPTIDES <- as.ab(c("DAL", "ORI", "TLV")) # dalba/orita/tela
pre_commit_lst$AB_GLYCOPEPTIDES_EXCEPT_LIPO <- pre_commit_lst$AB_GLYCOPEPTIDES[!pre_commit_lst$AB_GLYCOPEPTIDES %in% pre_commit_lst$AB_LIPOGLYCOPEPTIDES]
pre_commit_lst$AB_LINCOSAMIDES <- antimicrobials %>%
  filter(atc_group2 %like% "lincosamide" | (group %like% "lincosamide" & is.na(atc_group2) & name %like% "^(pirlimycin|clinda)")) %>%
  pull(ab)
pre_commit_lst$AB_MACROLIDES <- antimicrobials %>%
  filter(atc_group2 %like% "macrolide" | (group %like% "macrolide" & is.na(atc_group2)) | name %like% "^(acetylmidecamycin|acetylspiramycin|gamith?romycin|kitasamycin|meleumycin|nafith?romycin|primycin|solith?romycin|tildipirosin|tilmicosin|tulath?romycin|tylosin|tylvalosin)") %>%
  pull(ab)
pre_commit_lst$AB_MONOBACTAMS <- antimicrobials %>%
  filter(group %like% "monobactam") %>%
  pull(ab)
pre_commit_lst$AB_NITROFURANS <- antimicrobials %>%
  filter(name %like% "^furaz|nitrofura" | atc_group2 %like% "nitrofuran") %>%
  pull(ab)
pre_commit_lst$AB_OXAZOLIDINONES <- antimicrobials %>%
  filter(group %like% "oxazolidinone") %>%
  pull(ab)
pre_commit_lst$AB_PENICILLINS <- antimicrobials %>%
  filter(group %like% "penicillin" & !(name %unlike% "/" & name %like% ".*bactam$")) %>%
  pull(ab)
pre_commit_lst$AB_PHENICOLS <- antimicrobials %>%
  filter(group %like% "phenicol" | atc_group1 %like% "phenicol" | atc_group2 %like% "phenicol") %>%
  pull(ab)
pre_commit_lst$AB_PHOSPHONICS <- antimicrobials %>%
  filter(group %like% "phosphonic" | name %like% "fosfo") %>%
  pull(ab)
pre_commit_lst$AB_POLYMYXINS <- antimicrobials %>%
  filter(group %like% "polymyxin") %>%
  pull(ab)
pre_commit_lst$AB_QUINOLONES <- antimicrobials %>%
  filter(group %like% "quinolone" | atc_group1 %like% "quinolone" | atc_group2 %like% "quinolone" | name %like% "ozenoxacin") %>%
  pull(ab)
pre_commit_lst$AB_FLUOROQUINOLONES <- antimicrobials %>%
  # see DOI 10.23937/2378-3656/1410369, more specifically this table: https://www.clinmedjournals.org/articles/cmrcr/cmrcr-8-369-table1.html
  filter(ab %in% pre_commit_lst$AB_QUINOLONES & name %unlike% " acid|nalidixic|cinoxacin|flumequine|oxolinic|ozenoxacin|piromidic|pipemidic|rosoxacin") %>%
  pull(ab)
pre_commit_lst$AB_RIFAMYCINS <- antimicrobials %>%
  filter(name %like% "Rifampi|Rifabutin|Rifapentine|rifamy") %>%
  pull(ab)
pre_commit_lst$AB_SPIROPYRIMIDINETRIONES <- antimicrobials %>%
  filter(name %like% "zoliflodacin") %>%
  pull(ab)
pre_commit_lst$AB_STREPTOGRAMINS <- antimicrobials %>%
  filter(atc_group2 %like% "streptogramin") %>%
  pull(ab)
pre_commit_lst$AB_TETRACYCLINES <- antimicrobials %>%
  filter(atc_group1 %like% "tetracycline" | atc_group2 %like% "tetracycline" | name %like% "chlortetracycline|cetocycline|demeclocycline|doxycycline|eravacycline|lymecycline|meclocycline|meth?acycline|minocycline|omadacycline|oxytetracycline|rolitetracycline|sarecycline|tetracycline|tigecycline") %>%
  pull(ab)
pre_commit_lst$AB_TETRACYCLINES_EXCEPT_TGC <- pre_commit_lst$AB_TETRACYCLINES[pre_commit_lst$AB_TETRACYCLINES != "TGC"]
pre_commit_lst$AB_TRIMETHOPRIMS <- antimicrobials %>%
  filter(atc_group1 %like% "trimethoprim" | atc_group2 %like% "trimethoprim" | name %like% "trimethoprim|ormetroprim|iclaprim") %>%
  pull(ab)
pre_commit_lst$AB_SULFONAMIDES <- antimicrobials %>%
  filter(name %like% "(^|/)sulf[oai]") %>%
  pull(ab)
pre_commit_lst$AB_UREIDOPENICILLINS <- as.ab(c("PIP", "TZP", "AZL", "MEZ"))
pre_commit_lst$AB_BETALACTAMS <- sort(c(
  pre_commit_lst$AB_PENICILLINS,
  pre_commit_lst$AB_CEPHALOSPORINS,
  pre_commit_lst$AB_CARBAPENEMS,
  pre_commit_lst$AB_MONOBACTAMS))
pre_commit_lst$AB_BETALACTAMASE_INHIBITORS <- antimicrobials %>%
  filter(atc_group2 %like% "Beta-lactamase inhibitors" | name %like% "bactam") %>%
  pull(ab)
# for EUCAST:
pre_commit_lst$AB_BETALACTAMS_WITH_INHIBITOR <- antimicrobials %>%
  filter(ab %in% pre_commit_lst$AB_BETALACTAMS & name %like% "/" & name %unlike% "EDTA") %>%
  pull(ab)
# this will be used for documentation:
pre_commit_lst$DEFINED_AB_GROUPS <- sort(names(pre_commit_lst)[names(pre_commit_lst) %like% "^AB_" & names(pre_commit_lst) != "AB_LOOKUP"])

# Update the antimicrobials$group column
usethis::ui_info("Updating 'group' column in antimicrobials data set from AB_* vectors")
prettify_group_name <- function(name) {
  raw <- gsub("^AB_", "", name)
  pretty <- tools::toTitleCase(gsub("_", " ", tolower(raw)))
  pretty[pretty %like% " (except|with) "] <- ""
  pretty <- gsub(" (1st|2nd|3rd|4th|5th|6th)", " (\\1 gen.)", pretty)
  pretty <- gsub("([Bb])eta[-]?", "\\1eta-", pretty)
  pretty <- gsub(" Inhibitor", " inhibitor", pretty)
  pretty <- pretty[pretty != ""]
  return(pretty)
}
group_map <- vector("list", length = nrow(antimicrobials))
names(group_map) <- antimicrobials$ab
for (group_name in pre_commit_lst$DEFINED_AB_GROUPS) {
  ab_vector <- pre_commit_lst[[group_name]]
  pretty_name <- prettify_group_name(group_name)
  for (ab in ab_vector) {
    ab_chr <- as.character(ab)
    group_map[[ab_chr]] <- sort(unique(c(group_map[[ab_chr]], pretty_name)))
  }
}
for (i in seq_along(group_map)) {
  if (is.null(group_map[[i]])) {
    group_map[[i]] <- "Other"
    if (antimicrobials$group[i] %unlike% "other") {
      usethis::ui_warn("AB had a group but not anymore: ", antimicrobials$name[i], " (", antimicrobials$ab[i], "), was ", toString(antimicrobials$group[i]))
    }
  }
  group_map[[i]] <- group_map[[i]][order(nchar(group_map[[i]]))]
}

# create priority list for ab_group()
pre_commit_lst$ABX_PRIORITY_LIST <- c("Aminopenicillins",
                                      "Isoxazolylpenicillins",
                                      "Ureidopenicillins",
                                      "Oxazolidinones",
                                      "Carbapenems",
                                      "Cephalosporins (1st gen.)",
                                      "Cephalosporins (2nd gen.)",
                                      "Cephalosporins (3rd gen.)",
                                      "Cephalosporins (4th gen.)",
                                      "Cephalosporins (5th gen.)",
                                      "Cephalosporins",
                                      "Penicillins",
                                      "Monobactams",
                                      "Aminoglycosides",
                                      "Lipoglycopeptides",
                                      "Glycopeptides",
                                      "Lincosamides",
                                      "Streptogramins",
                                      "Macrolides",
                                      "Nitrofurans",
                                      "Phenicols",
                                      "Phosphonics",
                                      "Polymyxins",
                                      "Fluoroquinolones",
                                      "Quinolones",
                                      "Rifamycins",
                                      "Spiropyrimidinetriones",
                                      "Trimethoprims",
                                      "Sulfonamides",
                                      "Tetracyclines",
                                      "Antifungals",
                                      "Antimycobacterials",
                                      "Beta-lactams",
                                      "Beta-lactamase inhibitors",
                                      "Other")
if (!all(unlist(antimicrobials$group) %in% pre_commit_lst$ABX_PRIORITY_LIST)) {
  stop("Missing group(s) in priority list: ", paste(setdiff(unlist(antimicrobials$group), pre_commit_lst$ABX_PRIORITY_LIST), collapse = ", "))
}
for (i in seq_along(group_map)) {
  group_map[[i]] <- intersect(pre_commit_lst$ABX_PRIORITY_LIST, group_map[[i]])
}
antimicrobials$group <- unname(group_map)
usethis::use_data(antimicrobials, overwrite = TRUE, version = 2, compress = "xz")

pre_commit_lst$AB_LOOKUP <- create_AB_AV_lookup(antimicrobials)
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
  path <- paste0("data-raw/", deparse(substitute(object)), ".md5")
  if (!file.exists(path)) return(TRUE)
  tryCatch(
    {
      conn <- file(path)
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
  mutate(ab_name = antimicrobials$name[match(ab, antimicrobials$ab)], .after = ab)
if (changed_md5(clin_break)) {
  usethis::ui_info(paste0("Saving {usethis::ui_value('clinical_breakpoints')} to {usethis::ui_value('data-raw/datasets/')}"))
  write_md5(clin_break)
  try(saveRDS(clin_break, "data-raw/datasets/clinical_breakpoints.rds", version = 2, compress = "xz"), silent = TRUE)
  try(write.table(clinical_breakpoints, "data-raw/datasets/clinical_breakpoints.txt", sep = "\t", na = "", row.names = FALSE), silent = TRUE)
  try(haven::write_sav(clin_break, "data-raw/datasets/clinical_breakpoints.sav"), silent = TRUE)
  try(haven::write_dta(clin_break, "data-raw/datasets/clinical_breakpoints.dta"), silent = TRUE)
  try(openxlsx2::write_xlsx(clin_break, "data-raw/datasets/clinical_breakpoints.xlsx"), silent = TRUE)
  try(arrow::write_feather(clin_break, "data-raw/datasets/clinical_breakpoints.feather"), silent = TRUE)
  try(arrow::write_parquet(clin_break, "data-raw/datasets/clinical_breakpoints.parquet"), silent = TRUE)
}

if (changed_md5(microorganisms)) {
  usethis::ui_info(paste0("Saving {usethis::ui_value('microorganisms')} to {usethis::ui_value('data-raw/datasets/')}"))
  write_md5(microorganisms)
  try(saveRDS(microorganisms, "data-raw/datasets/microorganisms.rds", version = 2, compress = "xz"), silent = TRUE)
  max_50_snomed <- sapply(microorganisms$snomed, function(x) paste(x[seq_len(min(50, length(x), na.rm = TRUE))], collapse = " "))
  mo <- microorganisms
  mo$snomed <- max_50_snomed
  mo <- dplyr::mutate_if(mo, ~ !is.numeric(.), as.character)
  try(haven::write_sav(mo, "data-raw/datasets/microorganisms.sav"), silent = TRUE)
  try(haven::write_dta(mo, "data-raw/datasets/microorganisms.dta"), silent = TRUE)
  mo_all_snomed <- microorganisms %>% mutate_if(is.list, function(x) sapply(x, paste, collapse = ","))
  try(write.table(mo_all_snomed, "data-raw/datasets/microorganisms.txt", sep = "\t", na = "", row.names = FALSE), silent = TRUE)
  try(openxlsx2::write_xlsx(mo_all_snomed, "data-raw/datasets/microorganisms.xlsx"), silent = TRUE)
  try(arrow::write_feather(microorganisms, "data-raw/datasets/microorganisms.feather"), silent = TRUE)
  try(arrow::write_parquet(microorganisms, "data-raw/datasets/microorganisms.parquet"), silent = TRUE)
}

if (changed_md5(microorganisms.codes)) {
  usethis::ui_info(paste0("Saving {usethis::ui_value('microorganisms.codes')} to {usethis::ui_value('data-raw/datasets/')}"))
  write_md5(microorganisms.codes)
  try(saveRDS(microorganisms.codes, "data-raw/datasets/microorganisms.codes.rds", version = 2, compress = "xz"), silent = TRUE)
  try(write.table(microorganisms.codes, "data-raw/datasets/microorganisms.codes.txt", sep = "\t", na = "", row.names = FALSE), silent = TRUE)
  try(haven::write_sav(microorganisms.codes, "data-raw/datasets/microorganisms.codes.sav"), silent = TRUE)
  try(haven::write_dta(microorganisms.codes, "data-raw/datasets/microorganisms.codes.dta"), silent = TRUE)
  try(openxlsx2::write_xlsx(microorganisms.codes, "data-raw/datasets/microorganisms.codes.xlsx"), silent = TRUE)
  try(arrow::write_feather(microorganisms.codes, "data-raw/datasets/microorganisms.codes.feather"), silent = TRUE)
  try(arrow::write_parquet(microorganisms.codes, "data-raw/datasets/microorganisms.codes.parquet"), silent = TRUE)
}

if (changed_md5(microorganisms.groups)) {
  usethis::ui_info(paste0("Saving {usethis::ui_value('microorganisms.groups')} to {usethis::ui_value('data-raw/datasets/')}"))
  write_md5(microorganisms.groups)
  try(saveRDS(microorganisms.groups, "data-raw/datasets/microorganisms.groups.rds", version = 2, compress = "xz"), silent = TRUE)
  try(write.table(microorganisms.groups, "data-raw/datasets/microorganisms.groups.txt", sep = "\t", na = "", row.names = FALSE), silent = TRUE)
  try(haven::write_sav(microorganisms.groups, "data-raw/datasets/microorganisms.groups.sav"), silent = TRUE)
  try(haven::write_dta(microorganisms.groups, "data-raw/datasets/microorganisms.groups.dta"), silent = TRUE)
  try(openxlsx2::write_xlsx(microorganisms.groups, "data-raw/datasets/microorganisms.groups.xlsx"), silent = TRUE)
  try(arrow::write_feather(microorganisms.groups, "data-raw/datasets/microorganisms.groups.feather"), silent = TRUE)
  try(arrow::write_parquet(microorganisms.groups, "data-raw/datasets/microorganisms.groups.parquet"), silent = TRUE)
}

ab <- dplyr::mutate_if(antimicrobials, ~ !is.numeric(.), as.character)
if (changed_md5(ab)) {
  usethis::ui_info(paste0("Saving {usethis::ui_value('antimicrobials')} to {usethis::ui_value('data-raw/datasets/')}"))
  write_md5(ab)
  try(saveRDS(antimicrobials, "data-raw/datasets/antimicrobials.rds", version = 2, compress = "xz"), silent = TRUE)
  try(haven::write_sav(ab, "data-raw/datasets/antimicrobials.sav"), silent = TRUE)
  try(haven::write_dta(ab, "data-raw/datasets/antimicrobials.dta"), silent = TRUE)
  ab_lists <- antimicrobials %>% mutate_if(is.list, function(x) sapply(x, paste, collapse = ","))
  try(write.table(ab_lists, "data-raw/datasets/antimicrobials.txt", sep = "\t", na = "", row.names = FALSE), silent = TRUE)
  try(openxlsx2::write_xlsx(ab_lists, "data-raw/datasets/antimicrobials.xlsx"), silent = TRUE)
  try(arrow::write_feather(antimicrobials, "data-raw/datasets/antimicrobials.feather"), silent = TRUE)
  try(arrow::write_parquet(antimicrobials, "data-raw/datasets/antimicrobials.parquet"), silent = TRUE)
}

av <- dplyr::mutate_if(antivirals, ~ !is.numeric(.), as.character)
if (changed_md5(av)) {
  usethis::ui_info(paste0("Saving {usethis::ui_value('antivirals')} to {usethis::ui_value('data-raw/datasets/')}"))
  write_md5(av)
  try(saveRDS(antivirals, "data-raw/datasets/antivirals.rds", version = 2, compress = "xz"), silent = TRUE)
  try(haven::write_sav(av, "data-raw/datasets/antivirals.sav"), silent = TRUE)
  try(haven::write_dta(av, "data-raw/datasets/antivirals.dta"), silent = TRUE)
  av_lists <- antivirals %>% mutate_if(is.list, function(x) sapply(x, paste, collapse = ","))
  try(write.table(av_lists, "data-raw/datasets/antivirals.txt", sep = "\t", na = "", row.names = FALSE), silent = TRUE)
  try(openxlsx2::write_xlsx(av_lists, "data-raw/datasets/antivirals.xlsx"), silent = TRUE)
  try(arrow::write_feather(antivirals, "data-raw/datasets/antivirals.feather"), silent = TRUE)
  try(arrow::write_parquet(antivirals, "data-raw/datasets/antivirals.parquet"), silent = TRUE)
}

# give official names to ABs and MOs
intrinsicR <- data.frame(
  microorganism = mo_name(intrinsic_resistant$mo, language = NULL, keep_synonyms = TRUE, info = FALSE),
  antibiotic = ab_name(intrinsic_resistant$ab, language = NULL),
  stringsAsFactors = FALSE
)
if (changed_md5(intrinsicR)) {
  usethis::ui_info(paste0("Saving {usethis::ui_value('intrinsic_resistant')} to {usethis::ui_value('data-raw/datasets/')}"))
  write_md5(intrinsicR)
  try(saveRDS(intrinsicR, "data-raw/datasets/intrinsic_resistant.rds", version = 2, compress = "xz"), silent = TRUE)
  try(write.table(intrinsicR, "data-raw/datasets/intrinsic_resistant.txt", sep = "\t", na = "", row.names = FALSE), silent = TRUE)
  try(haven::write_sav(intrinsicR, "data-raw/datasets/intrinsic_resistant.sav"), silent = TRUE)
  try(haven::write_dta(intrinsicR, "data-raw/datasets/intrinsic_resistant.dta"), silent = TRUE)
  try(openxlsx2::write_xlsx(intrinsicR, "data-raw/datasets/intrinsic_resistant.xlsx"), silent = TRUE)
  try(arrow::write_feather(intrinsicR, "data-raw/datasets/intrinsic_resistant.feather"), silent = TRUE)
  try(arrow::write_parquet(intrinsicR, "data-raw/datasets/intrinsic_resistant.parquet"), silent = TRUE)
}

if (changed_md5(dosage)) {
  usethis::ui_info(paste0("Saving {usethis::ui_value('dosage')} to {usethis::ui_value('data-raw/datasets/')}"))
  write_md5(dosage)
  try(saveRDS(dosage, "data-raw/datasets/dosage.rds", version = 2, compress = "xz"), silent = TRUE)
  try(write.table(dosage, "data-raw/datasets/dosage.txt", sep = "\t", na = "", row.names = FALSE), silent = TRUE)
  try(haven::write_sav(dosage, "data-raw/datasets/dosage.sav"), silent = TRUE)
  try(haven::write_dta(dosage, "data-raw/datasets/dosage.dta"), silent = TRUE)
  try(openxlsx2::write_xlsx(dosage, "data-raw/datasets/dosage.xlsx"), silent = TRUE)
  try(arrow::write_feather(dosage, "data-raw/datasets/dosage.feather"), silent = TRUE)
  try(arrow::write_parquet(dosage, "data-raw/datasets/dosage.parquet"), silent = TRUE)
}

# Set `antibiotics` as a deprecated data set
antibiotics <- structure(antimicrobials, class = c("deprecated_amr_dataset", class(antimicrobials)))
usethis::use_data(antibiotics, internal = FALSE, overwrite = TRUE, compress = "xz", version = 2)
rm(antibiotics)

suppressMessages(reset_AMR_locale())

devtools::load_all(quiet = TRUE)
suppressMessages(set_AMR_locale("English"))

files_changed <- function(paths = "^(R|data)/") {
  tryCatch({
    changed_files <- system("git status", intern = TRUE)
    changed_files <- unlist(strsplit(changed_files, " "))
    any(changed_files %like% paths[paths != "R/sysdata.rda"])
  }, error = function(e) TRUE)
}

# Update URLs -------------------------------------------------------------
if (files_changed()) {
  usethis::ui_info("Checking URLs for redirects")
  # Step 1: Get sources from tools (excluding man/)
  sources <- tools:::url_db_from_package_sources(".")
  sources <- sources[!grepl("^man/", sources$Parent), ]
  # Step 2: Get URLs from .R files in R/
  r_files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
  # Function to extract URLs from a file
  extract_urls_from_file <- function(file_path) {
    lines <- readLines(file_path, warn = FALSE)
    urls <- stringr::str_extract_all(lines, "https?://[^\\s)\"'>]+")
    urls <- unlist(urls)
    if (length(urls) == 0) {
      return(NULL)
    }
    # Remove trailing punctuation (e.g., .,), etc.)
    urls <- stringr::str_replace(urls, "[\\.,;)]+$", "")
    data.frame(
      URL = urls,
      Parent = gsub("^\\./", "", file_path),
      stringsAsFactors = FALSE
    )
  }
  r_file_urls <- do.call(rbind, lapply(r_files, extract_urls_from_file))
  # Step 3: Combine the two sources
  total <- rbind(sources, r_file_urls)
  # Step 4: Check URLs and update
  results <- urlchecker::url_check(db = total)
  invisible(urlchecker::url_update(results = results))
}

# Style pkg ---------------------------------------------------------------
if (files_changed(paths = "^(R|tests)/")) {
  usethis::ui_info("Styling package")
  styler::style_pkg(include_roxygen_examples = FALSE,
                    exclude_dirs = list.dirs(full.names = FALSE, recursive = FALSE)[!list.dirs(full.names = FALSE, recursive = FALSE) %in% c("R", "tests")])
}

# Document pkg ------------------------------------------------------------
if (files_changed()) {
  usethis::ui_info("Documenting package")
  suppressMessages(devtools::document(quiet = TRUE))
}

# Update index.md and README.md -------------------------------------------
if (files_changed("README.Rmd") ||
    files_changed("index.Rmd") ||
    files_changed("man/microorganisms.Rd") ||
    files_changed("man/antimicrobials.Rd") ||
    files_changed("man/clinical_breakpoints.Rd") ||
    files_changed("man/antibiogram.Rd") ||
    files_changed("R/antibiogram.R") ||
    files_changed("data-raw/translations.tsv")) {
  usethis::ui_info("Rendering {usethis::ui_field('index.md')} and {usethis::ui_field('README.md')}")
  suppressWarnings(rmarkdown::render("index.Rmd", quiet = TRUE))
  suppressWarnings(rmarkdown::render("README.Rmd", quiet = TRUE))
  unlink("index.html") # remove previews from folder
  unlink("README.html")
}

# Finished ----------------------------------------------------------------
rm(antimicrobials)
usethis::ui_done("All done")
suppressMessages(reset_AMR_locale())
