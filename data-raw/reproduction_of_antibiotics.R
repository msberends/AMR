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

library(dplyr)

# got EARS-Net codes (= ECDC/WHO codes) from here:

# Installed WHONET 2019 software on Windows (http://www.whonet.org/software.html),
#    opened C:\WHONET\Codes\WHONETCodes.mdb in MS Access
#    and exported table 'DRGLST' to MS Excel
library(readxl)
DRGLST <- read_excel("DRGLST.xlsx")
abx <- DRGLST %>%
  select(ab = WHON5_CODE,
         name = ANTIBIOTIC) %>%
  # remove the ones without WHONET code
  filter(!is.na(ab)) %>%
  distinct(name, .keep_all = TRUE) %>%
  # add the ones without WHONET code
  bind_rows(
    DRGLST %>%
      select(ab = WHON5_CODE,
             name = ANTIBIOTIC) %>%
      filter(is.na(ab)) %>%
      distinct(name, .keep_all = TRUE)
      # add new ab code later
  ) %>%
  arrange(name)

# add old ATC codes
ab_old <- AMR::antibiotics %>%
  mutate(official = gsub("( and |, )", "/", official),
         abbr = tolower(paste(ifelse(is.na(abbr), "", abbr),
                      ifelse(is.na(certe), "", certe),
                      ifelse(is.na(umcg), "", umcg),
                      sep = "|")))
for (i in 1:nrow(ab_old)) {
  abbr <- ab_old[i, "abbr"]
  abbr <- strsplit(abbr, "|", fixed = TRUE) %>% unlist() %>% unique()
  abbr <- abbr[abbr != ""]
  #print(abbr)
  if (length(abbr) == 0) {
    ab_old[i, "abbr"] <- NA_character_
  } else {
    ab_old[i, "abbr"] <- paste(abbr, collapse = "|")
  }
}

# create reference data set: to be able to map ab to atc
abx_atc1 <- abx %>%
  mutate(name_lower = tolower(name)) %>%
  left_join(ab_old %>%
              select(ears_net, atc), by = c(ab = "ears_net")) %>%
  rename(atc1 = atc) %>%
  left_join(ab_old %>%
              mutate(official = gsub(", combinations", "", official, fixed = TRUE)) %>%
              transmute(official = tolower(official), atc), by = c(name_lower = "official")) %>%
  rename(atc2 = atc) %>%
  left_join(ab_old %>%
              mutate(official = gsub(", combinations", "", official, fixed = TRUE)) %>%
              mutate(official = gsub("f", "ph", official)) %>%
              transmute(official = tolower(official), atc), by = c(name_lower = "official")) %>%
  rename(atc3 = atc) %>%
  left_join(ab_old %>%
              mutate(official = gsub(", combinations", "", official, fixed = TRUE)) %>%
              mutate(official = gsub("t", "th", official)) %>%
              transmute(official = tolower(official), atc), by = c(name_lower = "official")) %>%
  rename(atc4 = atc) %>%
  left_join(ab_old %>%
              mutate(official = gsub(", combinations", "", official, fixed = TRUE)) %>%
              mutate(official = gsub("f", "ph", official)) %>%
              mutate(official = gsub("t", "th", official)) %>%
              transmute(official = tolower(official), atc), by = c(name_lower = "official")) %>%
  rename(atc5 = atc) %>%
  left_join(ab_old %>%
              mutate(official = gsub(", combinations", "", official, fixed = TRUE)) %>%
              mutate(official = gsub("f", "ph", official)) %>%
              mutate(official = gsub("t", "th", official)) %>%
              mutate(official = gsub("ine$", "in", official)) %>%
              transmute(official = tolower(official), atc), by = c(name_lower = "official")) %>%
  rename(atc6 = atc) %>%
  mutate(atc = case_when(!is.na(atc1) ~ atc1,
                         !is.na(atc2) ~ atc2,
                         !is.na(atc3) ~ atc3,
                         !is.na(atc4) ~ atc4,
                         !is.na(atc4) ~ atc5,
                         TRUE ~ atc6)) %>%
  distinct(ab, name, .keep_all = TRUE) %>%
  select(ab, atc, name)

abx_atc2 <- ab_old %>%
  filter(!atc %in% abx_atc1$atc,
         is.na(ears_net),
         !is.na(atc_group1),
         !atc_group1 %like% ("virus|vaccin|viral|immun"),
         !official %like% "(combinations| with )") %>%
  mutate(ab = NA_character_) %>%
  as.data.frame(stringsAsFactors = FALSE) %>%
  select(ab, atc, name = official)

abx2 <- bind_rows(abx_atc1, abx_atc2)

rm(abx_atc1)
rm(abx_atc2)

abx2$ab[is.na(abx2$ab)] <- toupper(abbreviate(gsub("[/0-9-]",
                                                   " ",
                                                   abx2$name[is.na(abx2$ab)]),
                                              minlength = 3,
                                              method = "left.kept",
                                              strict = TRUE))

n_distinct(abx2$ab)

abx2 <- abx2 %>% arrange(ab)
seqnr <- 0
# add follow up nrs
for (i in 2:nrow(abx2)) {
  if (abx2[i, "ab"] == abx2[i - 1, "ab"]) {
    seqnr <- seqnr + 1
    abx2[i, "seqnr"] <- seqnr
  } else {
    seqnr <- 0
  }
}
for (i in 2:nrow(abx2)) {
  if (!is.na(abx2[i, "seqnr"])) {
    abx2[i, "ab"] <- paste0(abx2[i, "ab"], abx2[i, "seqnr"])
  }
}
abx2 <- abx2 %>% select(-seqnr) %>% arrange(name)

# everything unique??
nrow(abx2) == n_distinct(abx2$ab)

# get ATC properties
abx2 <- abx2 %>%
  left_join(ab_old %>%
              select(atc, abbr, atc_group1, atc_group2,
                     oral_ddd, oral_units, iv_ddd, iv_units))

abx2$abbr <- lapply(as.list(abx2$abbr), function(x) unlist(strsplit(x, "|", fixed = TRUE)))

# vector with official names, returns vector with CIDs
get_CID <- function(ab) {
  CID <- rep(NA_integer_, length(ab))
  p <- progress_estimated(n = length(ab), min_time = 0)
  for (i in 1:length(ab)) {
    p$tick()$print()

    CID[i] <- tryCatch(
      data.table::fread(paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/",
                               URLencode(ab[i], reserved = TRUE),
                               "/cids/TXT?name_type=complete"),
                        showProgress = FALSE)[[1]][1],
      error = function(e) NA_integer_)
    if (is.na(CID[i])) {
      # try with removing the text in brackets
      CID[i] <- tryCatch(
        data.table::fread(paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/",
                                 URLencode(trimws(gsub("[(].*[)]", "", ab[i])), reserved = TRUE),
                                 "/cids/TXT?name_type=complete"),
                          showProgress = FALSE)[[1]][1],
        error = function(e) NA_integer_)
    }
    if (is.na(CID[i])) {
      # try match on word and take the lowest CID value (sorted)
      ab[i] <- gsub("[^a-z0-9]+", " ", ab[i], ignore.case = TRUE)
      CID[i] <- tryCatch(
        data.table::fread(paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/",
                                 URLencode(ab[i], reserved = TRUE),
                                 "/cids/TXT?name_type=word"),
                          showProgress = FALSE)[[1]][1],
        error = function(e) NA_integer_)
    }
    Sys.sleep(0.1)
  }
  CID
}

# get CIDs (2-3 min)
CIDs <- get_CID(abx2$name)
# These could not be found:
abx2[is.na(CIDs),] %>% View()

# returns list with synonyms (brand names), with CIDs as names
get_synonyms <- function(CID, clean = TRUE) {
  synonyms <- rep(NA_character_, length(CID))
  p <- progress_estimated(n = length(CID), min_time = 0)

  for (i in 1:length(CID)) {
    p$tick()$print()

    synonyms_txt <- ""

    if (is.na(CID[i])) {
      next
    }

    synonyms_txt <- tryCatch(
      data.table::fread(paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastidentity/cid/",
                               CID[i],
                               "/synonyms/TXT"),
                        sep = "\n",
                        showProgress = FALSE)[[1]],
      error = function(e) NA_character_)

    Sys.sleep(0.1)

    if (clean == TRUE) {
      # remove text between brackets
      synonyms_txt <- trimws(gsub("[(].*[)]", "",
                                  gsub("[[].*[]]", "",
                                       gsub("[(].*[]]", "",
                                            gsub("[[].*[)]", "", synonyms_txt)))))
      synonyms_txt <- gsub("Co-", "Co", synonyms_txt, fixed = TRUE)
      # only length 6 to 20 and no txt with reading marks or numbers and must start with capital letter (= brand)
      synonyms_txt <- synonyms_txt[nchar(synonyms_txt) %in% c(6:20)
                                   & !grepl("[-&{},_0-9/]", synonyms_txt)
                                   & grepl("^[A-Z]", synonyms_txt, ignore.case = FALSE)]
      synonyms_txt <- unlist(strsplit(synonyms_txt,  ";", fixed = TRUE))
    }
    synonyms_txt <- unique(trimws(synonyms_txt[tolower(synonyms_txt) %in% unique(tolower(synonyms_txt))]))
    synonyms[i] <- list(sort(synonyms_txt))
  }
  names(synonyms) <- CID
  synonyms
}

# get brand names from PubChem (2-3 min)
synonyms <- get_synonyms(CIDs)
synonyms <- lapply(synonyms,
                   function(x) {
                     if (length(x) == 0 | all(is.na(x))) {
                       ""
                     } else {
                       x
                     }})

# add them to data set
antibiotics <- abx2 %>%
  left_join(DRGLST %>%
              select(ab = WHON5_CODE, CLASS, SUBCLASS) %>%
              distinct(ab, .keep_all = TRUE), by = "ab") %>%
  transmute(ab,
            atc,
            cid = CIDs,
            # no capital after a slash: Ampicillin/Sulbactam -> Ampicillin/sulbactam
            name = name %>%
              gsub("([/-])([A-Z])", "\\1\\L\\2", ., perl = TRUE) %>%
              gsub("edta", "EDTA", ., ignore.case = TRUE),
            group = case_when(
              paste(atc_group1, atc_group2, CLASS, SUBCLASS) %like% "am(ph|f)enicol" ~ "Amphenicols",
              paste(atc_group1, atc_group2, CLASS, SUBCLASS) %like% "aminoglycoside" ~ "Aminoglycosides",
              paste(atc_group1, atc_group2, CLASS, SUBCLASS) %like% "carbapenem" | name %like% "(imipenem|meropenem)" ~ "Carbapenems",
              paste(atc_group1, atc_group2, CLASS, SUBCLASS) %like% "First-generation cephalosporin" ~ "Cephalosporins (1st gen.)",
              paste(atc_group1, atc_group2, CLASS, SUBCLASS) %like% "Second-generation cephalosporin" ~ "Cephalosporins (2nd gen.)",
              paste(atc_group1, atc_group2, CLASS, SUBCLASS) %like% "Third-generation cephalosporin" ~ "Cephalosporins (3rd gen.)",
              paste(atc_group1, atc_group2, CLASS, SUBCLASS) %like% "Fourth-generation cephalosporin" ~ "Cephalosporins (4th gen.)",
              paste(atc_group1, atc_group2, CLASS, SUBCLASS) %like% "(tuberculosis|mycobacter)" ~ "Antimycobacterials",
              paste(atc_group1, atc_group2, CLASS, SUBCLASS) %like% "cephalosporin" ~ "Cephalosporins",
              name %like% "^Ce" & is.na(atc_group1) & paste(atc_group1, atc_group2, CLASS, SUBCLASS) %like% "beta-?lactam" ~ "Cephalosporins",
              paste(atc_group1, atc_group2, CLASS, SUBCLASS) %like% "(beta-?lactam|penicillin)" ~ "Beta-lactams/penicillins",
              paste(atc_group1, atc_group2, CLASS, SUBCLASS) %like% "quinolone" ~ "Quinolones",
              paste(atc_group1, atc_group2, CLASS, SUBCLASS) %like% "glycopeptide" ~ "Glycopeptides",
              paste(atc_group1, atc_group2, CLASS, SUBCLASS) %like% "macrolide" ~ "Macrolides/lincosamides",
              paste(atc_group1, atc_group2, CLASS, SUBCLASS) %like% "tetracycline" ~ "Tetracyclines",
              paste(atc_group1, atc_group2, CLASS, SUBCLASS) %like% "trimethoprim" ~ "Trimethoprims",
              paste(atc_group1, atc_group2, CLASS, SUBCLASS) %like% "polymyxin" ~ "Polymyxins",
              paste(atc_group1, atc_group2, CLASS, SUBCLASS) %like% "(fungal|mycot)" ~ "Antifungals/antimycotics",
              TRUE ~ "Other antibacterials"
            ),
            atc_group1, atc_group2,
            abbreviations = unname(abbr),
            synonyms = unname(synonyms),
            oral_ddd, oral_units,
            iv_ddd, iv_units) %>%
  as.data.frame(stringsAsFactors = FALSE)

# some exceptions
antibiotics[which(antibiotics$ab == "DOX"), "abbreviations"][[1]] <- list(c("dox", "doxy"))
antibiotics[which(antibiotics$ab == "FLC"), "abbreviations"][[1]] <- list(c("clox"))
antibiotics[which(antibiotics$ab == "CEC"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "CEC"), "abbreviations"][[1]], "CFC")) # cefaclor old WHONET4 code
antibiotics[which(antibiotics$ab == "AMX"), "synonyms"][[1]] <- list(sort(c(antibiotics[which(antibiotics$ab == "AMX"), "synonyms"][[1]], "Amoxy")))
# 'Polymixin B' (POL) and 'Polymyxin B' (PLB) both exist, so:
antibiotics[which(antibiotics$ab == "PLB"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "PLB"), "abbreviations"][[1]], "POL", "Polymixin", "Polymixin B", "Poly B"))
antibiotics <- filter(antibiotics, ab != "POL")
# 'Latamoxef' (LTM) and 'Moxalactam (Latamoxef)' (MOX) both exist, so:
antibiotics[which(antibiotics$ab == "LTM"), "abbreviations"][[1]] <- list(c("MOX", "moxa"))
antibiotics <- filter(antibiotics, ab != "MOX")
# RFP and RFP1 (the J0 one) both mean 'rifapentine', although 'rifp' is not recognised, so:
antibiotics <- filter(antibiotics, ab != "RFP")
antibiotics[which(antibiotics$ab == "RFP1"), "ab"] <- "RFP"
antibiotics[which(antibiotics$ab == "RFP"), "abbreviations"][[1]] <- list(c("rifp"))
# Rifampicin is better known as a drug than Rifampin (Rifampin is still listed as a brand name), so:
antibiotics[which(antibiotics$ab == "RIF"), "name"] <- "Rifampicin"
# PME and PVM1 (the J0 one) both mean 'Pivmecillinam', so:
antibiotics <- filter(antibiotics, ab != "PME")
antibiotics[which(antibiotics$ab == "PVM1"), "ab"] <- "PME"
# Remove Sinecatechins
antibiotics <- filter(antibiotics, ab != "SNC")
# GLIMS codes
antibiotics[which(antibiotics$ab == as.ab("cefuroxim")), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == as.ab("cefuroxim")), "abbreviations"][[1]], "cfrx"))
antibiotics[which(antibiotics$ab == as.ab("cefotaxim")), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == as.ab("cefotaxim")), "abbreviations"][[1]], "cftx"))
antibiotics[which(antibiotics$ab == as.ab("ceftazidime")), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == as.ab("ceftazidime")), "abbreviations"][[1]], "cftz"))
antibiotics[which(antibiotics$ab == as.ab("cefepime")), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == as.ab("cefepime")), "abbreviations"][[1]], "cfpi"))
antibiotics[which(antibiotics$ab == as.ab("cefoxitin")), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == as.ab("cefoxitin")), "abbreviations"][[1]], "cfxt", "cfsc"))
# More GLIMS codes
antibiotics[which(antibiotics$ab == "CAZ"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "CAZ"), "abbreviations"][[1]], "cftz"))
antibiotics[which(antibiotics$ab == "CRO"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "CRO"), "abbreviations"][[1]], "cftr"))
antibiotics[which(antibiotics$ab == "CTX"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "CTX"), "abbreviations"][[1]], "cftx"))
antibiotics[which(antibiotics$ab == "CXM"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "CXM"), "abbreviations"][[1]], "cfrx"))
antibiotics[which(antibiotics$ab == "CZO"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "CZO"), "abbreviations"][[1]], "cfzl"))
antibiotics[which(antibiotics$ab == "FOX"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "FOX"), "abbreviations"][[1]], "cfxt"))
antibiotics[which(antibiotics$ab == "PIP"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "PIP"), "abbreviations"][[1]], "pipc"))
antibiotics[which(antibiotics$ab == "PIP"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "PIP"), "abbreviations"][[1]], "PIPC"))
antibiotics[which(antibiotics$ab == "SXT"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "SXT"), "abbreviations"][[1]], "COTRIM"))
antibiotics[which(antibiotics$ab == "SXT"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "SXT"), "abbreviations"][[1]], "trsx"))
antibiotics[which(antibiotics$ab == "TZP"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "TZP"), "abbreviations"][[1]], "PIPTAZ"))
antibiotics[which(antibiotics$ab == "TZP"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "TZP"), "abbreviations"][[1]], "pita"))
# official RIVM codes (Dutch National Health Institute)
# https://www.rivm.nl/sites/default/files/2019-09/Bijlage_4_Lijst_antibiotica%202020%201.0.pdf
antibiotics[which(antibiotics$ab == "FCT"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "FCT"), "abbreviations"][[1]], "5flc"))
antibiotics[which(antibiotics$ab == "AMC"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "AMC"), "abbreviations"][[1]], "amcl"))
antibiotics[which(antibiotics$ab == "AMB"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "AMB"), "abbreviations"][[1]], "amfb"))
antibiotics[which(antibiotics$ab == "AMH"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "AMH"), "abbreviations"][[1]], "amhl"))
antibiotics[which(antibiotics$ab == "AMK"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "AMK"), "abbreviations"][[1]], "amik"))
antibiotics[which(antibiotics$ab == "AMX"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "AMX"), "abbreviations"][[1]], "amox"))
antibiotics[which(antibiotics$ab == "AMP"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "AMP"), "abbreviations"][[1]], "ampi"))
antibiotics[which(antibiotics$ab == "SAM"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "SAM"), "abbreviations"][[1]], "amsu"))
antibiotics[which(antibiotics$ab == "ANI"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "ANI"), "abbreviations"][[1]], "anid"))
antibiotics[which(antibiotics$ab == "SAM"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "SAM"), "abbreviations"][[1]], "apsu"))
antibiotics[which(antibiotics$ab == "AZM"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "AZM"), "abbreviations"][[1]], "azit"))
antibiotics[which(antibiotics$ab == "AZL"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "AZL"), "abbreviations"][[1]], "azlo"))
antibiotics[which(antibiotics$ab == "ATM"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "ATM"), "abbreviations"][[1]], "aztr"))
antibiotics[which(antibiotics$ab == "PNV"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "PNV"), "abbreviations"][[1]], "bepe"))
antibiotics[which(antibiotics$ab == "CAP"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "CAP"), "abbreviations"][[1]], "capr"))
antibiotics[which(antibiotics$ab == "CRB"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "CRB"), "abbreviations"][[1]], "carb"))
antibiotics[which(antibiotics$ab == "CAS"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "CAS"), "abbreviations"][[1]], "casp"))
antibiotics[which(antibiotics$ab == "CDC"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "CDC"), "abbreviations"][[1]], "cecl"))
antibiotics[which(antibiotics$ab == "CXA"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "CXA"), "abbreviations"][[1]], "cfax"))
antibiotics[which(antibiotics$ab == "CTB"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "CTB"), "abbreviations"][[1]], "cfbu"))
antibiotics[which(antibiotics$ab == "CEC"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "CEC"), "abbreviations"][[1]], "cfcl"))
antibiotics[which(antibiotics$ab == "CFR"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "CFR"), "abbreviations"][[1]], "cfdx"))
antibiotics[which(antibiotics$ab == "CEP"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "CEP"), "abbreviations"][[1]], "cflt"))
antibiotics[which(antibiotics$ab == "LEX"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "LEX"), "abbreviations"][[1]], "cflx"))
antibiotics[which(antibiotics$ab == "MAN"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "MAN"), "abbreviations"][[1]], "cfmn"))
antibiotics[which(antibiotics$ab == "CPD"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "CPD"), "abbreviations"][[1]], "cfpd"))
antibiotics[which(antibiotics$ab == "FEP"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "FEP"), "abbreviations"][[1]], "cfpi"))
antibiotics[which(antibiotics$ab == "CPO"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "CPO"), "abbreviations"][[1]], "cfpr"))
antibiotics[which(antibiotics$ab == "CFP"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "CFP"), "abbreviations"][[1]], "cfpz"))
antibiotics[which(antibiotics$ab == "CED"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "CED"), "abbreviations"][[1]], "cfrd"))
antibiotics[which(antibiotics$ab == "CPT"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "CPT"), "abbreviations"][[1]], "cfro"))
antibiotics[which(antibiotics$ab == "CXM"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "CXM"), "abbreviations"][[1]], "cfrx"))
antibiotics[which(antibiotics$ab == "CFS"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "CFS"), "abbreviations"][[1]], "cfsl"))
antibiotics[which(antibiotics$ab == "CRO"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "CRO"), "abbreviations"][[1]], "cftr"))
antibiotics[which(antibiotics$ab == "CTT"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "CTT"), "abbreviations"][[1]], "cftt"))
antibiotics[which(antibiotics$ab == "CTX"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "CTX"), "abbreviations"][[1]], "cftx"))
antibiotics[which(antibiotics$ab == "CAZ"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "CAZ"), "abbreviations"][[1]], "cftz"))
antibiotics[which(antibiotics$ab == "CFM"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "CFM"), "abbreviations"][[1]], "cfxm"))
antibiotics[which(antibiotics$ab == "FOX"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "FOX"), "abbreviations"][[1]], "cfxt"))
antibiotics[which(antibiotics$ab == "CZA"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "CZA"), "abbreviations"][[1]], "cfav"))
antibiotics[which(antibiotics$ab == "CZO"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "CZO"), "abbreviations"][[1]], "cfzl"))
antibiotics[which(antibiotics$ab == "CZX"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "CZX"), "abbreviations"][[1]], "cfzx"))
antibiotics[which(antibiotics$ab == "CHL"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "CHL"), "abbreviations"][[1]], "chlo"))
antibiotics[which(antibiotics$ab == "CPC"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "CPC"), "abbreviations"][[1]], "cicl"))
antibiotics[which(antibiotics$ab == "CIN"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "CIN"), "abbreviations"][[1]], "cino"))
antibiotics[which(antibiotics$ab == "CIP"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "CIP"), "abbreviations"][[1]], "cipr"))
antibiotics[which(antibiotics$ab == "CIX"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "CIX"), "abbreviations"][[1]], "cipx"))
antibiotics[which(antibiotics$ab == "CLR"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "CLR"), "abbreviations"][[1]], "clar"))
antibiotics[which(antibiotics$ab == "CLI"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "CLI"), "abbreviations"][[1]], "clin"))
antibiotics[which(antibiotics$ab == "CTR"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "CTR"), "abbreviations"][[1]], "clot"))
antibiotics[which(antibiotics$ab == "CLO"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "CLO"), "abbreviations"][[1]], "clox"))
antibiotics[which(antibiotics$ab == "COL"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "COL"), "abbreviations"][[1]], "coli"))
antibiotics[which(antibiotics$ab == "CTC"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "CTC"), "abbreviations"][[1]], "cxcl"))
antibiotics[which(antibiotics$ab == "CYC"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "CYC"), "abbreviations"][[1]], "cycl"))
antibiotics[which(antibiotics$ab == "CCV"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "CCV"), "abbreviations"][[1]], "czcl"))
antibiotics[which(antibiotics$ab == "DAP"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "DAP"), "abbreviations"][[1]], "dapt"))
antibiotics[which(antibiotics$ab == "DIC"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "DIC"), "abbreviations"][[1]], "dicl"))
antibiotics[which(antibiotics$ab == "DOR"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "DOR"), "abbreviations"][[1]], "dori"))
antibiotics[which(antibiotics$ab == "DOX"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "DOX"), "abbreviations"][[1]], "doxy"))
antibiotics[which(antibiotics$ab == "ENX"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "ENX"), "abbreviations"][[1]], "enox"))
antibiotics[which(antibiotics$ab == "ETP"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "ETP"), "abbreviations"][[1]], "erta"))
antibiotics[which(antibiotics$ab == "ERY"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "ERY"), "abbreviations"][[1]], "eryt"))
antibiotics[which(antibiotics$ab == "PHE"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "PHE"), "abbreviations"][[1]], "fene"))
antibiotics[which(antibiotics$ab == "PHN"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "PHN"), "abbreviations"][[1]], "fepe"))
antibiotics[which(antibiotics$ab == "FLE"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "FLE"), "abbreviations"][[1]], "fler"))
antibiotics[which(antibiotics$ab == "FLU"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "FLU"), "abbreviations"][[1]], "fluc"))
antibiotics[which(antibiotics$ab == "FLC"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "FLC"), "abbreviations"][[1]], "flux"))
antibiotics[which(antibiotics$ab == "FOS"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "FOS"), "abbreviations"][[1]], "fosf"))
antibiotics[which(antibiotics$ab == "FRM"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "FRM"), "abbreviations"][[1]], "fram"))
antibiotics[which(antibiotics$ab == "FUS"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "FUS"), "abbreviations"][[1]], "fusi"))
antibiotics[which(antibiotics$ab == "GAT"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "GAT"), "abbreviations"][[1]], "gati"))
antibiotics[which(antibiotics$ab == "GEH"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "GEH"), "abbreviations"][[1]], "gehl"))
antibiotics[which(antibiotics$ab == "GEN"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "GEN"), "abbreviations"][[1]], "gent"))
antibiotics[which(antibiotics$ab == "GRX"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "GRX"), "abbreviations"][[1]], "grep"))
antibiotics[which(antibiotics$ab == "IPM"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "IPM"), "abbreviations"][[1]], "imci"))
antibiotics[which(antibiotics$ab == "IPM"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "IPM"), "abbreviations"][[1]], "imip"))
antibiotics[which(antibiotics$ab == "ISV"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "ISV"), "abbreviations"][[1]], "isav"))
antibiotics[which(antibiotics$ab == "ITR"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "ITR"), "abbreviations"][[1]], "itra"))
antibiotics[which(antibiotics$ab == "KAH"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "KAH"), "abbreviations"][[1]], "kahl"))
antibiotics[which(antibiotics$ab == "KAN"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "KAN"), "abbreviations"][[1]], "kana"))
antibiotics[which(antibiotics$ab == "KET"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "KET"), "abbreviations"][[1]], "keto"))
antibiotics[which(antibiotics$ab == "LVX"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "LVX"), "abbreviations"][[1]], "levo"))
antibiotics[which(antibiotics$ab == "LIN"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "LIN"), "abbreviations"][[1]], "linc"))
antibiotics[which(antibiotics$ab == "LNZ"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "LNZ"), "abbreviations"][[1]], "line"))
antibiotics[which(antibiotics$ab == "LOR"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "LOR"), "abbreviations"][[1]], "lora"))
antibiotics[which(antibiotics$ab == "MEM"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "MEM"), "abbreviations"][[1]], "mero"))
antibiotics[which(antibiotics$ab == "MET"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "MET"), "abbreviations"][[1]], "meti"))
antibiotics[which(antibiotics$ab == "MTR"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "MTR"), "abbreviations"][[1]], "metr"))
antibiotics[which(antibiotics$ab == "MEZ"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "MEZ"), "abbreviations"][[1]], "mezl"))
antibiotics[which(antibiotics$ab == "MIF"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "MIF"), "abbreviations"][[1]], "mica"))
antibiotics[which(antibiotics$ab == "MCZ"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "MCZ"), "abbreviations"][[1]], "mico"))
antibiotics[which(antibiotics$ab == "MNO"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "MNO"), "abbreviations"][[1]], "mino"))
antibiotics[which(antibiotics$ab == "LTM"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "LTM"), "abbreviations"][[1]], "moxa"))
antibiotics[which(antibiotics$ab == "MFX"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "MFX"), "abbreviations"][[1]], "moxi"))
antibiotics[which(antibiotics$ab == "NAL"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "NAL"), "abbreviations"][[1]], "nali"))
antibiotics[which(antibiotics$ab == "NEO"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "NEO"), "abbreviations"][[1]], "neom"))
antibiotics[which(antibiotics$ab == "NET"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "NET"), "abbreviations"][[1]], "neti"))
antibiotics[which(antibiotics$ab == "NIT"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "NIT"), "abbreviations"][[1]], "nitr"))
antibiotics[which(antibiotics$ab == "NOR"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "NOR"), "abbreviations"][[1]], "norf"))
antibiotics[which(antibiotics$ab == "NYS"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "NYS"), "abbreviations"][[1]], "nyst"))
antibiotics[which(antibiotics$ab == "OFX"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "OFX"), "abbreviations"][[1]], "oflo"))
antibiotics[which(antibiotics$ab == "OXA"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "OXA"), "abbreviations"][[1]], "oxal"))
antibiotics[which(antibiotics$ab == "PEF"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "PEF"), "abbreviations"][[1]], "pefl"))
antibiotics[which(antibiotics$ab == "PEN"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "PEN"), "abbreviations"][[1]], "peni"))
antibiotics[which(antibiotics$ab == "PIP"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "PIP"), "abbreviations"][[1]], "pipc"))
antibiotics[which(antibiotics$ab == "PPA"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "PPA"), "abbreviations"][[1]], "pipz"))
antibiotics[which(antibiotics$ab == "TZP"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "TZP"), "abbreviations"][[1]], "pita"))
antibiotics[which(antibiotics$ab == "PLB"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "PLB"), "abbreviations"][[1]], "polb"))
antibiotics[which(antibiotics$ab == "POS"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "POS"), "abbreviations"][[1]], "posa"))
antibiotics[which(antibiotics$ab == "PRI"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "PRI"), "abbreviations"][[1]], "pris"))
antibiotics[which(antibiotics$ab == "QDA"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "QDA"), "abbreviations"][[1]], "quda"))
antibiotics[which(antibiotics$ab == "RIF"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "RIF"), "abbreviations"][[1]], "rifa"))
antibiotics[which(antibiotics$ab == "RXT"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "RXT"), "abbreviations"][[1]], "roxi"))
antibiotics[which(antibiotics$ab == "SMX"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "SMX"), "abbreviations"][[1]], "sfmx"))
antibiotics[which(antibiotics$ab == "SLF4"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "SLF4"), "abbreviations"][[1]], "sfmz"))
antibiotics[which(antibiotics$ab == "SSS"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "SSS"), "abbreviations"][[1]], "sfna"))
antibiotics[which(antibiotics$ab == "SLF"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "SLF"), "abbreviations"][[1]], "sfsz"))
antibiotics[which(antibiotics$ab == "SPX"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "SPX"), "abbreviations"][[1]], "spar"))
antibiotics[which(antibiotics$ab == "SPT"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "SPT"), "abbreviations"][[1]], "spec"))
antibiotics[which(antibiotics$ab == "SPI"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "SPI"), "abbreviations"][[1]], "spir"))
antibiotics[which(antibiotics$ab == "STH"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "STH"), "abbreviations"][[1]], "sthl"))
antibiotics[which(antibiotics$ab == "STR1"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "STR1"), "abbreviations"][[1]], "stre"))
antibiotics[which(antibiotics$ab == "TAZ"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "TAZ"), "abbreviations"][[1]], "tazo"))
antibiotics[which(antibiotics$ab == "TEC"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "TEC"), "abbreviations"][[1]], "teic"))
antibiotics[which(antibiotics$ab == "TLT"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "TLT"), "abbreviations"][[1]], "teli"))
antibiotics[which(antibiotics$ab == "TMX"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "TMX"), "abbreviations"][[1]], "tema"))
antibiotics[which(antibiotics$ab == "TEM"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "TEM"), "abbreviations"][[1]], "temo"))
antibiotics[which(antibiotics$ab == "TRB"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "TRB"), "abbreviations"][[1]], "terb"))
antibiotics[which(antibiotics$ab == "TCY"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "TCY"), "abbreviations"][[1]], "tetr"))
antibiotics[which(antibiotics$ab == "TIC"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "TIC"), "abbreviations"][[1]], "tica"))
antibiotics[which(antibiotics$ab == "TCC"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "TCC"), "abbreviations"][[1]], "ticl"))
antibiotics[which(antibiotics$ab == "TGC"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "TGC"), "abbreviations"][[1]], "tige"))
antibiotics[which(antibiotics$ab == "TIN"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "TIN"), "abbreviations"][[1]], "tini"))
antibiotics[which(antibiotics$ab == "TOB"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "TOB"), "abbreviations"][[1]], "tobr"))
antibiotics[which(antibiotics$ab == "TOH"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "TOH"), "abbreviations"][[1]], "tohl"))
antibiotics[which(antibiotics$ab == "TMP"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "TMP"), "abbreviations"][[1]], "trim"))
antibiotics[which(antibiotics$ab == "TVA"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "TVA"), "abbreviations"][[1]], "trov"))
antibiotics[which(antibiotics$ab == "SLT4"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "SLT4"), "abbreviations"][[1]], "trsm"))
antibiotics[which(antibiotics$ab == "SXT"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "SXT"), "abbreviations"][[1]], "trsx"))
antibiotics[which(antibiotics$ab == "VAN"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "VAN"), "abbreviations"][[1]], "vanc"))
antibiotics[which(antibiotics$ab == "VOR"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "VOR"), "abbreviations"][[1]], "vori"))
antibiotics <- antibiotics %>% 
  mutate(ab = as.character(ab)) %>% 
  rbind(antibiotics %>% 
          filter(ab == "GEH") %>% 
          mutate(ab = "AMH",
                 name = "Amphotericin B-high",
                 abbreviations = list(c("amhl", "amfo b high", "ampho b high", "amphotericin high")))) %>% 
  rbind(antibiotics %>% 
          filter(ab == "GEH") %>% 
          mutate(ab = "TOH",
                 name = "Tobramycin-high",
                 abbreviations = list(c("tohl", "tobra high", "tobramycin high")))) %>% 
  rbind(antibiotics %>% 
          filter(ab == "BUT") %>% 
          mutate(ab = "CIX",
                 atc = "D01AE14",
                 name = "Ciclopirox",
                 group = "Antifungals/antimycotics",
                 atc_group1 = "Antifungals for topical use",
                 atc_group2 = "Other antifungals for topical use",
                 abbreviations = list(c("cipx"))))
antibiotics[which(antibiotics$ab == "SSS"), "name"] <- "Sulfonamide"
# ESBL E-test codes:
antibiotics[which(antibiotics$ab == "CCV"), "abbreviations"][[1]] <- list(c("xtzl"))
antibiotics[which(antibiotics$ab == "CAZ"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "CAZ"), "abbreviations"][[1]], "xtz", "cefta"))
antibiotics[which(antibiotics$ab == "CPC"), "abbreviations"][[1]] <- list(c("xpml"))
antibiotics[which(antibiotics$ab == "FEP"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "FEP"), "abbreviations"][[1]], "xpm"))
antibiotics[which(antibiotics$ab == "CTC"), "abbreviations"][[1]] <- list(c("xctl"))
antibiotics[which(antibiotics$ab == "CTX"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "CTX"), "abbreviations"][[1]], "xct"))
# High level Gentamcin and Streptomycin
antibiotics[which(antibiotics$ab == "GEH"), "abbreviations"][[1]] <- list(c("gehl", "gentamicin high", "genta high"))
antibiotics[which(antibiotics$ab == "STH"), "abbreviations"][[1]] <- list(c("sthl", "streptomycin high", "strepto high"))
# add imi and "imipenem/cilastatine" to imipenem
antibiotics[which(antibiotics$ab == "IPM"), "abbreviations"][[1]] <- list(c("imip", "imi", "imp"))
antibiotics[which(antibiotics$ab == "IPM"), "synonyms"][[1]] <- list(sort(c(antibiotics[which(antibiotics$ab == "IPM"), "synonyms"][[1]], "imipenem/cilastatin")))
# add synonyms of ones not found
antibiotics[which(antibiotics$ab == "TZP"), "synonyms"][[1]] <- list(sort(c(antibiotics[which(antibiotics$ab == "TZP"), "synonyms"][[1]], "Tazocel", "tazocillin", "Tazocin", "Zosyn")))
antibiotics[which(antibiotics$ab == "COL"), "synonyms"][[1]] <- list(sort(unique(c(antibiotics[which(antibiotics$ab == "COL"), "synonyms"][[1]], "Colisticin", "Polymyxin E", "Colimycin", "Coly-Mycin", "Totazina", "Colistimethate", "Promixin", "Colistimethate Sodium"))))
# remove incorrect synonyms from rifampicin (RIF) and add them to the combination rifampicin/isoniazid (RFI)
old_sym <- antibiotics[which(antibiotics$ab == "RIF"), "synonyms"][[1]]
old_sym <- old_sym[!old_sym %in% c("Rifinah", "Rimactazid")]
antibiotics[which(antibiotics$ab == "RIF"), "synonyms"][[1]] <- list(old_sym)
antibiotics[which(antibiotics$ab == "RFI"), "synonyms"][[1]] <- list(sort(c("Rifinah", "Rimactazid")))
# remove incorrect synonyms from sulfamethoxazole (SMX) and add them to the combination trimethoprim/sulfamethoxazole (SXT)
old_sym <- antibiotics[which(antibiotics$ab == "SMX"), "synonyms"][[1]]
old_sym <- old_sym[!old_sym %in% c("Cotrimoxazole", "Bactrimel")]
antibiotics[which(antibiotics$ab == "SMX"), "synonyms"][[1]] <- list(old_sym)
antibiotics[which(antibiotics$ab == "SXT"), "synonyms"][[1]] <- list(sort(unique(c(antibiotics[which(antibiotics$ab == "COL"), "synonyms"][[1]], "Cotrimoxazole", "Bactrimel", "Septra", "Bactrim", "Cotrimazole"))))

## new ATC codes
# ceftaroline
antibiotics[which(antibiotics$ab == "CPT"), "atc"] <- "J01DI02"
# faropenem
antibiotics[which(antibiotics$ab == "FAR"), "atc"] <- "J01DI03"
# ceftobiprole
antibiotics[which(antibiotics$ab == "BPR"), "atc"] <- "J01DI01"

# typo
antibiotics[which(antibiotics$ab == "RXT"), "name"] <- "Roxithromycin"
antibiotics[which(antibiotics$ab == "PEN"), "atc"] <- "J01CE01"

# set cephalosporins groups for the ones that could not be determined automatically:
antibiotics <- antibiotics %>% 
  mutate(group = case_when(
    name == "Cefcapene" ~ "Cephalosporins (3rd gen.)",
    name == "Cefcapene pivoxil" ~ "Cephalosporins (3rd gen.)",
    name == "Cefditoren pivoxil" ~ "Cephalosporins (3rd gen.)",
    name == "Cefepime/clavulanic acid" ~ "Cephalosporins (4th gen.)",
    name == "Cefepime/tazobactam" ~ "Cephalosporins (4th gen.)",
    name == "Cefetamet pivoxil" ~ "Cephalosporins (3rd gen.)",
    name == "Cefetecol (Cefcatacol)" ~ "Cephalosporins (4th gen.)",
    name == "Cefetrizole" ~ "Cephalosporins (unclassified gen.)",
    name == "Cefoselis" ~ "Cephalosporins (4th gen.)",
    name == "Cefotaxime/clavulanic acid" ~ "Cephalosporins (3rd gen.)",
    name == "Cefotaxime/sulbactam" ~ "Cephalosporins (3rd gen.)",
    name == "Cefotiam hexetil" ~ "Cephalosporins (3rd gen.)",
    name == "Cefovecin" ~ "Cephalosporins (3rd gen.)",
    name == "Cefozopran" ~ "Cephalosporins (4th gen.)",
    name == "Cefpimizole" ~ "Cephalosporins (3rd gen.)",
    name == "Cefpodoxime proxetil" ~ "Cephalosporins (3rd gen.)",
    name == "Cefpodoxime/clavulanic acid" ~ "Cephalosporins (3rd gen.)",
    name == "Cefquinome" ~ "Cephalosporins (4th gen.)",
    name == "Cefsumide" ~ "Cephalosporins (unclassified gen.)",
    name == "Ceftaroline" ~ "Cephalosporins (5th gen.)",
    name == "Ceftaroline/avibactam" ~ "Cephalosporins (5th gen.)",
    name == "Ceftazidime/avibactam" ~ "Cephalosporins (3rd gen.)",
    name == "Cefteram" ~ "Cephalosporins (3rd gen.)",
    name == "Cefteram pivoxil" ~ "Cephalosporins (3rd gen.)",
    name == "Ceftiofur" ~ "Cephalosporins (3rd gen.)",
    name == "Ceftizoxime alapivoxil" ~ "Cephalosporins (3rd gen.)",
    name == "Ceftobiprole" ~ "Cephalosporins (5th gen.)",
    name == "Ceftobiprole medocaril" ~ "Cephalosporins (5th gen.)",
    name == "Ceftolozane/enzyme inhibitor" ~ "Cephalosporins (5th gen.)",
    name == "Ceftolozane/tazobactam" ~ "Cephalosporins (5th gen.)",
    name == "Cefuroxime axetil" ~ "Cephalosporins (2nd gen.)",
    TRUE ~ group))

# set as data.frame again
antibiotics <- as.data.frame(antibiotics, stringsAsFactors = FALSE)
class(antibiotics$ab) <- c("ab", "character")
antibiotics <- antibiotics %>% arrange(name)

# make all abbreviations and synonyms lower case, unique and alphabetically sorted ----
for (i in 1:nrow(antibiotics)) {
  abb <- sort(unique(tolower(antibiotics[i, "abbreviations"][[1]])))
  syn <- sort(unique(tolower(antibiotics[i, "synonyms"][[1]])))
  antibiotics[i, "abbreviations"][[1]] <- ifelse(length(abb[!abb == ""]) == 0, list(""), list(abb))
  antibiotics[i, "synonyms"][[1]] <- ifelse(length(syn[!syn == ""]) == 0, list(""), list(syn))
}

# REFER TO data-raw/loinc.R FOR ADDING LOINC CODES

usethis::use_data(antibiotics, overwrite = TRUE)
rm(antibiotics)
