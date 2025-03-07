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

library(dplyr)

# got EARS-Net codes (= ECDC/WHO codes) from here:

# Installed WHONET 2019 software on Windows (http://www.whonet.org/software.html),
#    opened C:\WHONET\Codes\WHONETCodes.mdb in MS Access
#    and exported table 'DRGLST' to MS Excel
library(readxl)
DRGLST <- read_excel("DRGLST.xlsx")
abx <- DRGLST %>%
  select(
    ab = WHON5_CODE,
    name = ANTIBIOTIC
  ) %>%
  # remove the ones without WHONET code
  filter(!is.na(ab)) %>%
  distinct(name, .keep_all = TRUE) %>%
  # add the ones without WHONET code
  bind_rows(
    DRGLST %>%
      select(
        ab = WHON5_CODE,
        name = ANTIBIOTIC
      ) %>%
      filter(is.na(ab)) %>%
      distinct(name, .keep_all = TRUE)
    # add new ab code later
  ) %>%
  arrange(name)

# add old ATC codes
ab_old <- AMR::antimicrobials %>%
  mutate(
    official = gsub("( and |, )", "/", official),
    abbr = tolower(paste(ifelse(is.na(abbr), "", abbr),
      ifelse(is.na(certe), "", certe),
      ifelse(is.na(umcg), "", umcg),
      sep = "|"
    ))
  )
for (i in 1:nrow(ab_old)) {
  abbr <- ab_old[i, "abbr"]
  abbr <- strsplit(abbr, "|", fixed = TRUE) %>%
    unlist() %>%
    unique()
  abbr <- abbr[abbr != ""]
  # print(abbr)
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
  mutate(atc = case_when(
    !is.na(atc1) ~ atc1,
    !is.na(atc2) ~ atc2,
    !is.na(atc3) ~ atc3,
    !is.na(atc4) ~ atc4,
    !is.na(atc4) ~ atc5,
    TRUE ~ atc6
  )) %>%
  distinct(ab, name, .keep_all = TRUE) %>%
  select(ab, atc, name)

abx_atc2 <- ab_old %>%
  filter(
    !atc %in% abx_atc1$atc,
    is.na(ears_net),
    !is.na(atc_group1),
    atc_group1 %unlike% ("virus|vaccin|viral|immun"),
    official %unlike% "(combinations| with )"
  ) %>%
  mutate(ab = NA_character_) %>%
  as.data.frame(stringsAsFactors = FALSE) %>%
  select(ab, atc, name = official)

abx2 <- bind_rows(abx_atc1, abx_atc2)

rm(abx_atc1)
rm(abx_atc2)

abx2$ab[is.na(abx2$ab)] <- toupper(abbreviate(
  gsub(
    "[/0-9-]",
    " ",
    abx2$name[is.na(abx2$ab)]
  ),
  minlength = 3,
  method = "left.kept",
  strict = TRUE
))

n_distinct(abx2$ab)

abx2 <- abx2 %>% arrange(ab)
seqnr <- 0
# add follow up nrs
for (i in 2:nrow(abx2)) {
  if (abx2[i, "ab", drop = TRUE] == abx2[i - 1, "ab", drop = TRUE]) {
    seqnr <- seqnr + 1
    abx2[i, "seqnr"] <- seqnr
  } else {
    seqnr <- 0
  }
}
for (i in 2:nrow(abx2)) {
  if (!is.na(abx2[i, "seqnr"])) {
    abx2[i, "ab"] <- paste0(abx2[i, "ab", drop = TRUE], abx2[i, "seqnr", drop = TRUE])
  }
}
abx2 <- abx2 %>%
  select(-seqnr) %>%
  arrange(name)

# everything unique??
nrow(abx2) == n_distinct(abx2$ab)

# get ATC properties
abx2 <- abx2 %>%
  left_join(ab_old %>%
    select(
      atc, abbr, atc_group1, atc_group2,
      oral_ddd, oral_units, iv_ddd, iv_units
    ))

abx2$abbr <- lapply(as.list(abx2$abbr), function(x) unlist(strsplit(x, "|", fixed = TRUE)))

# Update Compound IDs and Synonyms ----

# vector with official names, returns vector with CIDs
get_CID <- function(ab) {
  CID <- rep(NA_integer_, length(ab))
  p <- AMR:::progress_ticker(n = length(ab), min_time = 0)
  for (i in 1:length(ab)) {
    p$tick()

    CID[i] <- tryCatch(
      data.table::fread(
        paste0(
          "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/",
          URLencode(ab[i], reserved = TRUE),
          "/cids/TXT?name_type=complete"
        ),
        showProgress = FALSE
      )[[1]][1],
      error = function(e) NA_integer_
    )
    if (is.na(CID[i])) {
      # try with removing the text in brackets
      CID[i] <- tryCatch(
        data.table::fread(
          paste0(
            "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/",
            URLencode(trimws(gsub("[(].*[)]", "", ab[i])), reserved = TRUE),
            "/cids/TXT?name_type=complete"
          ),
          showProgress = FALSE
        )[[1]][1],
        error = function(e) NA_integer_
      )
    }
    if (is.na(CID[i])) {
      # try match on word and take the lowest CID value (sorted)
      ab[i] <- gsub("[^a-z0-9]+", " ", ab[i], ignore.case = TRUE)
      CID[i] <- tryCatch(
        data.table::fread(
          paste0(
            "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/",
            URLencode(ab[i], reserved = TRUE),
            "/cids/TXT?name_type=word"
          ),
          showProgress = FALSE
        )[[1]][1],
        error = function(e) NA_integer_
      )
    }
    Sys.sleep(0.1)
  }
  CID
}

# get CIDs (4-5 min)
CIDs <- get_CID(antimicrobials$name)
# take missing from previously found CIDs
CIDs[is.na(CIDs) & !is.na(antimicrobials$cid)] <- antimicrobials$cid[is.na(CIDs) & !is.na(antimicrobials$cid)]
# These could not be found:
antimicrobials[is.na(CIDs), ] %>% View()

# returns list with synonyms (brand names), with CIDs as names
get_synonyms <- function(CID, clean = TRUE) {
  synonyms <- rep(NA_character_, length(CID))
  p <- AMR:::progress_ticker(n = length(CID), min_time = 0)

  for (i in 1:length(CID)) {
    p$tick()

    synonyms_txt <- ""

    if (is.na(CID[i])) {
      next
    }
    
    # we will now get the closest compounds with a 96% threshold
    similar_cids <- tryCatch(
      data.table::fread(
        paste0(
          "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastsimilarity_2d/cid/",
          CID[i],
          "/cids/TXT?Threshold=96&MaxRecords=5"
        ),
        sep = "\n",
        showProgress = FALSE
      )[[1]],
      error = function(e) NA_character_
    )
    # include the current CID of course
    all_cids <- unique(c(CID[i], similar_cids))
    # but leave out all CIDs that we have in our antimicrobials dataset to prevent duplication
    all_cids <- all_cids[!all_cids %in% antimicrobials$cid[!is.na(antimicrobials$cid)]]
    # for each one, we are getting the synonyms
    current_syns <- character(0)
    for (j in seq_len(length(all_cids))) {
      synonyms_txt <- tryCatch(
        data.table::fread(
          paste0(
            "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastidentity/cid/",
            all_cids[j],
            "/synonyms/TXT"
          ),
          sep = "\n",
          showProgress = FALSE
        )[[1]],
        error = function(e) NA_character_
      )
      
      Sys.sleep(0.05)
      
      if (clean == TRUE) {
        # remove text between brackets
        synonyms_txt <- trimws(gsub(
          "[(].*[)]", "",
          gsub(
            "[[].*[]]", "",
            gsub(
              "[(].*[]]", "",
              gsub("[[].*[)]", "", synonyms_txt)
            )
          )
        ))
        synonyms_txt <- gsub("Co-", "Co", synonyms_txt, fixed = TRUE)
        synonyms_txt <- gsub(" ?(mono)?sodium ?", "", ignore.case = TRUE, synonyms_txt)
        synonyms_txt <- gsub(" ?(injection|pediatric) ?", "", ignore.case = TRUE, synonyms_txt)
        # only length 6 to 20 and no txt with reading marks or numbers and must start with capital letter (= brand)
        synonyms_txt <- synonyms_txt[nchar(synonyms_txt) %in% c(5:20) &
                                       !grepl("[-&{},_0-9/:]", synonyms_txt) &
                                       grepl("^[A-Z]", synonyms_txt, ignore.case = FALSE)]
        synonyms_txt <- unlist(strsplit(synonyms_txt, ";", fixed = TRUE))
      }
      
      current_syns <- c(current_syns, synonyms_txt)
    }
    
    current_syns <- unique(trimws(current_syns[tolower(current_syns) %in% unique(tolower(current_syns))]))
    synonyms[i] <- list(sort(current_syns))
  }
  names(synonyms) <- CID
  synonyms
}

# get brand names from PubChem (3-4 min)
synonyms <- get_synonyms(CIDs)
synonyms.bak <- synonyms
synonyms <- synonyms.bak

# add existing ones (will be cleaned later)
for (i in seq_len(length(synonyms))) {
  old <- unname(unlist(AMR::antimicrobials[i, "synonyms", drop = TRUE]))
  synonyms[[i]] <- c(unname(synonyms[[i]]), old)
}

antimicrobials$synonyms <- synonyms

stop("remember to remove co-trimoxazole as synonyms from SMX (Sulfamethoxazole), so it only exists in SXT!")
sulfa <- antimicrobials[which(antimicrobials$ab == "SMX"), "synonyms", drop = TRUE][[1]]
cotrim <- antimicrobials[which(antimicrobials$ab == "SXT"), "synonyms", drop = TRUE][[1]]
# 2024-10-06 not the case anymore, no overlapping names: sulfa[sulfa %in% cotrim]
sulfa <- sulfa[!sulfa %in% cotrim]
antimicrobials[which(antimicrobials$ab == "SMX"), "synonyms"][[1]][[1]] <- sulfa


# now go to end of this file


# -----

# add them to data set
antimicrobials <- abx2 %>%
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
    iv_ddd, iv_units
  ) %>%
  as.data.frame(stringsAsFactors = FALSE)

# some exceptions
antimicrobials[which(antimicrobials$ab == "DOX"), "abbreviations"][[1]] <- list(c("dox", "doxy"))
antimicrobials[which(antimicrobials$ab == "FLC"), "abbreviations"][[1]] <- list(c("clox"))
antimicrobials[which(antimicrobials$ab == "CEC"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "CEC"), "abbreviations"][[1]], "CFC")) # cefaclor old WHONET4 code
antimicrobials[which(antimicrobials$ab == "AMX"), "synonyms"][[1]] <- list(sort(c(antimicrobials[which(antimicrobials$ab == "AMX"), "synonyms"][[1]], "Amoxy")))
# 'Polymixin B' (POL) and 'Polymyxin B' (PLB) both exist, so:
antimicrobials[which(antimicrobials$ab == "PLB"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "PLB"), "abbreviations"][[1]], "POL", "Polymixin", "Polymixin B", "Poly B"))
antimicrobials <- filter(antimicrobials, ab != "POL")
# 'Latamoxef' (LTM) and 'Moxalactam (Latamoxef)' (MOX) both exist, so:
antimicrobials[which(antimicrobials$ab == "LTM"), "abbreviations"][[1]] <- list(c("MOX", "moxa"))
antimicrobials <- filter(antimicrobials, ab != "MOX")
# RFP and RFP1 (the J0 one) both mean 'rifapentine', although 'rifp' is not recognised, so:
antimicrobials <- filter(antimicrobials, ab != "RFP")
antimicrobials[which(antimicrobials$ab == "RFP1"), "ab"] <- "RFP"
antimicrobials[which(antimicrobials$ab == "RFP"), "abbreviations"][[1]] <- list(c("rifp"))
# Rifampicin is better known as a drug than Rifampin (Rifampin is still listed as a brand name), so:
antimicrobials[which(antimicrobials$ab == "RIF"), "name"] <- "Rifampicin"
# PME and PVM1 (the J0 one) both mean 'Pivmecillinam', so:
antimicrobials <- filter(antimicrobials, ab != "PME")
antimicrobials[which(antimicrobials$ab == "PVM1"), "ab"] <- "PME"
# Remove Sinecatechins
antimicrobials <- filter(antimicrobials, ab != "SNC")
# GLIMS codes
antimicrobials[which(antimicrobials$ab == as.ab("cefuroxim")), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == as.ab("cefuroxim")), "abbreviations"][[1]], "cfrx"))
antimicrobials[which(antimicrobials$ab == as.ab("cefotaxim")), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == as.ab("cefotaxim")), "abbreviations"][[1]], "cftx"))
antimicrobials[which(antimicrobials$ab == as.ab("ceftazidime")), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == as.ab("ceftazidime")), "abbreviations"][[1]], "cftz"))
antimicrobials[which(antimicrobials$ab == as.ab("cefepime")), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == as.ab("cefepime")), "abbreviations"][[1]], "cfpi"))
antimicrobials[which(antimicrobials$ab == as.ab("cefoxitin")), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == as.ab("cefoxitin")), "abbreviations"][[1]], "cfxt"))
# Add cefoxitin screening
class(antimicrobials$ab) <- "character"
antimicrobials <- rbind(antimicrobials, data.frame(
  ab = "FOX1", atc = NA, cid = NA,
  name = "Cefoxitin screening",
  group = "Cephalosporins (2nd gen.)", atc_group1 = NA, atc_group2 = NA,
  abbreviations = "cfsc", synonyms = NA,
  oral_ddd = NA, oral_units = NA, iv_ddd = NA, iv_units = NA,
  loinc = NA,
  stringsAsFactors = FALSE
))
# More GLIMS codes
antimicrobials[which(antimicrobials$ab == "AMB"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "AMB"), "abbreviations"][[1]], "amf"))
antimicrobials[which(antimicrobials$ab == "CAZ"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "CAZ"), "abbreviations"][[1]], "cftz"))
antimicrobials[which(antimicrobials$ab == "COL"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "COL"), "abbreviations"][[1]], "cst"))
antimicrobials[which(antimicrobials$ab == "CRO"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "CRO"), "abbreviations"][[1]], "cftr"))
antimicrobials[which(antimicrobials$ab == "CTX"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "CTX"), "abbreviations"][[1]], "cftx"))
antimicrobials[which(antimicrobials$ab == "CXM"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "CXM"), "abbreviations"][[1]], "cfrx"))
antimicrobials[which(antimicrobials$ab == "CZO"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "CZO"), "abbreviations"][[1]], "cfzl"))
antimicrobials[which(antimicrobials$ab == "FCT"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "FCT"), "abbreviations"][[1]], "fcu"))
antimicrobials[which(antimicrobials$ab == "FCT"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "FCT"), "abbreviations"][[1]], "fluy"))
antimicrobials[which(antimicrobials$ab == "FLU"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "FLU"), "abbreviations"][[1]], "flz"))
antimicrobials[which(antimicrobials$ab == "FOS"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "FOS"), "abbreviations"][[1]], "fof"))
antimicrobials[which(antimicrobials$ab == "FOX"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "FOX"), "abbreviations"][[1]], "cfxt"))
antimicrobials[which(antimicrobials$ab == "FUS"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "FUS"), "abbreviations"][[1]], "fa"))
antimicrobials[which(antimicrobials$ab == "GEH"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "GEH"), "abbreviations"][[1]], "g_h"))
antimicrobials[which(antimicrobials$ab == "KAH"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "KAH"), "abbreviations"][[1]], "k_h"))
antimicrobials[which(antimicrobials$ab == "KET"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "KET"), "abbreviations"][[1]], "ktc"))
antimicrobials[which(antimicrobials$ab == "PIP"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "PIP"), "abbreviations"][[1]], "pipc"))
antimicrobials[which(antimicrobials$ab == "PIP"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "PIP"), "abbreviations"][[1]], "PIPC"))
antimicrobials[which(antimicrobials$ab == "SPX"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "SPX"), "abbreviations"][[1]], "spa"))
antimicrobials[which(antimicrobials$ab == "STH"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "STH"), "abbreviations"][[1]], "s_h"))
antimicrobials[which(antimicrobials$ab == "STR1"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "STR1"), "abbreviations"][[1]], "stm"))
antimicrobials[which(antimicrobials$ab == "SXT"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "SXT"), "abbreviations"][[1]], "COTRIM"))
antimicrobials[which(antimicrobials$ab == "SXT"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "SXT"), "abbreviations"][[1]], "trsx"))
antimicrobials[which(antimicrobials$ab == "TGC"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "TGC"), "abbreviations"][[1]], "tig"))
antimicrobials[which(antimicrobials$ab == "TMP"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "TMP"), "abbreviations"][[1]], "tri"))
antimicrobials[which(antimicrobials$ab == "TZP"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "TZP"), "abbreviations"][[1]], "PIPTAZ"))
antimicrobials[which(antimicrobials$ab == "TZP"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "TZP"), "abbreviations"][[1]], "pit"))
antimicrobials[which(antimicrobials$ab == "TZP"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "TZP"), "abbreviations"][[1]], "pita"))
antimicrobials[which(antimicrobials$ab == "VOR"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "VOR"), "abbreviations"][[1]], "vrc"))

# official RIVM codes (Dutch National Health Institute)
# https://www.rivm.nl/sites/default/files/2019-09/Bijlage_4_Lijst_antibiotica%202020%201.0.pdf
antimicrobials[which(antimicrobials$ab == "FCT"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "FCT"), "abbreviations"][[1]], "5flc"))
antimicrobials[which(antimicrobials$ab == "AMC"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "AMC"), "abbreviations"][[1]], "amcl"))
antimicrobials[which(antimicrobials$ab == "AMB"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "AMB"), "abbreviations"][[1]], "amfb"))
antimicrobials[which(antimicrobials$ab == "AMH"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "AMH"), "abbreviations"][[1]], "amhl"))
antimicrobials[which(antimicrobials$ab == "AMK"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "AMK"), "abbreviations"][[1]], "amik"))
antimicrobials[which(antimicrobials$ab == "AMX"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "AMX"), "abbreviations"][[1]], "amox"))
antimicrobials[which(antimicrobials$ab == "AMP"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "AMP"), "abbreviations"][[1]], "ampi"))
antimicrobials[which(antimicrobials$ab == "SAM"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "SAM"), "abbreviations"][[1]], "amsu"))
antimicrobials[which(antimicrobials$ab == "ANI"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "ANI"), "abbreviations"][[1]], "anid"))
antimicrobials[which(antimicrobials$ab == "SAM"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "SAM"), "abbreviations"][[1]], "apsu"))
antimicrobials[which(antimicrobials$ab == "AZM"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "AZM"), "abbreviations"][[1]], "azit"))
antimicrobials[which(antimicrobials$ab == "AZL"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "AZL"), "abbreviations"][[1]], "azlo"))
antimicrobials[which(antimicrobials$ab == "ATM"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "ATM"), "abbreviations"][[1]], "aztr"))
antimicrobials[which(antimicrobials$ab == "PNV"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "PNV"), "abbreviations"][[1]], "bepe"))
antimicrobials[which(antimicrobials$ab == "CAP"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "CAP"), "abbreviations"][[1]], "capr"))
antimicrobials[which(antimicrobials$ab == "CRB"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "CRB"), "abbreviations"][[1]], "carb"))
antimicrobials[which(antimicrobials$ab == "CAS"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "CAS"), "abbreviations"][[1]], "casp"))
antimicrobials[which(antimicrobials$ab == "CDC"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "CDC"), "abbreviations"][[1]], "cecl"))
antimicrobials[which(antimicrobials$ab == "CXA"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "CXA"), "abbreviations"][[1]], "cfax"))
antimicrobials[which(antimicrobials$ab == "CTB"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "CTB"), "abbreviations"][[1]], "cfbu"))
antimicrobials[which(antimicrobials$ab == "CEC"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "CEC"), "abbreviations"][[1]], "cfcl"))
antimicrobials[which(antimicrobials$ab == "CFR"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "CFR"), "abbreviations"][[1]], "cfdx"))
antimicrobials[which(antimicrobials$ab == "CEP"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "CEP"), "abbreviations"][[1]], "cflt"))
antimicrobials[which(antimicrobials$ab == "LEX"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "LEX"), "abbreviations"][[1]], "cflx"))
antimicrobials[which(antimicrobials$ab == "MAN"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "MAN"), "abbreviations"][[1]], "cfmn"))
antimicrobials[which(antimicrobials$ab == "CPD"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "CPD"), "abbreviations"][[1]], "cfpd"))
antimicrobials[which(antimicrobials$ab == "FEP"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "FEP"), "abbreviations"][[1]], "cfpi"))
antimicrobials[which(antimicrobials$ab == "CPO"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "CPO"), "abbreviations"][[1]], "cfpr"))
antimicrobials[which(antimicrobials$ab == "CFP"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "CFP"), "abbreviations"][[1]], "cfpz"))
antimicrobials[which(antimicrobials$ab == "CED"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "CED"), "abbreviations"][[1]], "cfrd"))
antimicrobials[which(antimicrobials$ab == "CPT"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "CPT"), "abbreviations"][[1]], "cfro"))
antimicrobials[which(antimicrobials$ab == "CXM"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "CXM"), "abbreviations"][[1]], "cfrx"))
antimicrobials[which(antimicrobials$ab == "CFS"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "CFS"), "abbreviations"][[1]], "cfsl"))
antimicrobials[which(antimicrobials$ab == "CRO"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "CRO"), "abbreviations"][[1]], "cftr"))
antimicrobials[which(antimicrobials$ab == "CTT"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "CTT"), "abbreviations"][[1]], "cftt"))
antimicrobials[which(antimicrobials$ab == "CTX"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "CTX"), "abbreviations"][[1]], "cftx"))
antimicrobials[which(antimicrobials$ab == "CAZ"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "CAZ"), "abbreviations"][[1]], "cftz"))
antimicrobials[which(antimicrobials$ab == "CFM"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "CFM"), "abbreviations"][[1]], "cfxm"))
antimicrobials[which(antimicrobials$ab == "FOX"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "FOX"), "abbreviations"][[1]], "cfxt"))
antimicrobials[which(antimicrobials$ab == "CZA"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "CZA"), "abbreviations"][[1]], "cfav"))
antimicrobials[which(antimicrobials$ab == "CZO"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "CZO"), "abbreviations"][[1]], "cfzl"))
antimicrobials[which(antimicrobials$ab == "CZX"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "CZX"), "abbreviations"][[1]], "cfzx"))
antimicrobials[which(antimicrobials$ab == "CHL"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "CHL"), "abbreviations"][[1]], "chlo"))
antimicrobials[which(antimicrobials$ab == "CPC"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "CPC"), "abbreviations"][[1]], "cicl"))
antimicrobials[which(antimicrobials$ab == "CIN"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "CIN"), "abbreviations"][[1]], "cino"))
antimicrobials[which(antimicrobials$ab == "CIP"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "CIP"), "abbreviations"][[1]], "cipr"))
antimicrobials[which(antimicrobials$ab == "CIX"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "CIX"), "abbreviations"][[1]], "cipx"))
antimicrobials[which(antimicrobials$ab == "CLR"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "CLR"), "abbreviations"][[1]], "clar"))
antimicrobials[which(antimicrobials$ab == "CLI"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "CLI"), "abbreviations"][[1]], "clin"))
antimicrobials[which(antimicrobials$ab == "CTR"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "CTR"), "abbreviations"][[1]], "clot"))
antimicrobials[which(antimicrobials$ab == "CLO"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "CLO"), "abbreviations"][[1]], "clox"))
antimicrobials[which(antimicrobials$ab == "COL"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "COL"), "abbreviations"][[1]], "coli"))
antimicrobials[which(antimicrobials$ab == "CTC"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "CTC"), "abbreviations"][[1]], "cxcl"))
antimicrobials[which(antimicrobials$ab == "CYC"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "CYC"), "abbreviations"][[1]], "cycl"))
antimicrobials[which(antimicrobials$ab == "CCV"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "CCV"), "abbreviations"][[1]], "czcl"))
antimicrobials[which(antimicrobials$ab == "DAP"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "DAP"), "abbreviations"][[1]], "dapt"))
antimicrobials[which(antimicrobials$ab == "DIC"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "DIC"), "abbreviations"][[1]], "dicl"))
antimicrobials[which(antimicrobials$ab == "DOR"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "DOR"), "abbreviations"][[1]], "dori"))
antimicrobials[which(antimicrobials$ab == "DOX"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "DOX"), "abbreviations"][[1]], "doxy"))
antimicrobials[which(antimicrobials$ab == "ENX"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "ENX"), "abbreviations"][[1]], "enox"))
antimicrobials[which(antimicrobials$ab == "ETP"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "ETP"), "abbreviations"][[1]], "erta"))
antimicrobials[which(antimicrobials$ab == "ERY"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "ERY"), "abbreviations"][[1]], "eryt"))
antimicrobials[which(antimicrobials$ab == "PHE"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "PHE"), "abbreviations"][[1]], "fene"))
antimicrobials[which(antimicrobials$ab == "PHN"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "PHN"), "abbreviations"][[1]], "fepe"))
antimicrobials[which(antimicrobials$ab == "FLE"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "FLE"), "abbreviations"][[1]], "fler"))
antimicrobials[which(antimicrobials$ab == "FLU"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "FLU"), "abbreviations"][[1]], "fluc"))
antimicrobials[which(antimicrobials$ab == "FLC"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "FLC"), "abbreviations"][[1]], "flux"))
antimicrobials[which(antimicrobials$ab == "FOS"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "FOS"), "abbreviations"][[1]], "fosf"))
antimicrobials[which(antimicrobials$ab == "FRM"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "FRM"), "abbreviations"][[1]], "fram"))
antimicrobials[which(antimicrobials$ab == "FUS"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "FUS"), "abbreviations"][[1]], "fusi"))
antimicrobials[which(antimicrobials$ab == "GAT"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "GAT"), "abbreviations"][[1]], "gati"))
antimicrobials[which(antimicrobials$ab == "GEH"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "GEH"), "abbreviations"][[1]], "gehl"))
antimicrobials[which(antimicrobials$ab == "GEN"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "GEN"), "abbreviations"][[1]], "gent"))
antimicrobials[which(antimicrobials$ab == "GRX"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "GRX"), "abbreviations"][[1]], "grep"))
antimicrobials[which(antimicrobials$ab == "IPM"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "IPM"), "abbreviations"][[1]], "imci"))
antimicrobials[which(antimicrobials$ab == "IPM"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "IPM"), "abbreviations"][[1]], "imip"))
antimicrobials[which(antimicrobials$ab == "ISV"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "ISV"), "abbreviations"][[1]], "isav"))
antimicrobials[which(antimicrobials$ab == "ITR"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "ITR"), "abbreviations"][[1]], "itra"))
antimicrobials[which(antimicrobials$ab == "KAH"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "KAH"), "abbreviations"][[1]], "kahl"))
antimicrobials[which(antimicrobials$ab == "KAN"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "KAN"), "abbreviations"][[1]], "kana"))
antimicrobials[which(antimicrobials$ab == "KET"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "KET"), "abbreviations"][[1]], "keto"))
antimicrobials[which(antimicrobials$ab == "LVX"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "LVX"), "abbreviations"][[1]], "levo"))
antimicrobials[which(antimicrobials$ab == "LIN"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "LIN"), "abbreviations"][[1]], "linc"))
antimicrobials[which(antimicrobials$ab == "LNZ"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "LNZ"), "abbreviations"][[1]], "line"))
antimicrobials[which(antimicrobials$ab == "LOR"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "LOR"), "abbreviations"][[1]], "lora"))
antimicrobials[which(antimicrobials$ab == "MEM"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "MEM"), "abbreviations"][[1]], "mero"))
antimicrobials[which(antimicrobials$ab == "MET"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "MET"), "abbreviations"][[1]], "meti"))
antimicrobials[which(antimicrobials$ab == "MTR"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "MTR"), "abbreviations"][[1]], "metr"))
antimicrobials[which(antimicrobials$ab == "MEZ"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "MEZ"), "abbreviations"][[1]], "mezl"))
antimicrobials[which(antimicrobials$ab == "MIF"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "MIF"), "abbreviations"][[1]], "mica"))
antimicrobials[which(antimicrobials$ab == "MCZ"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "MCZ"), "abbreviations"][[1]], "mico"))
antimicrobials[which(antimicrobials$ab == "MNO"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "MNO"), "abbreviations"][[1]], "mino"))
antimicrobials[which(antimicrobials$ab == "LTM"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "LTM"), "abbreviations"][[1]], "moxa", "moxalactam"))
antimicrobials[which(antimicrobials$ab == "MFX"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "MFX"), "abbreviations"][[1]], "moxi"))
antimicrobials[which(antimicrobials$ab == "NAL"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "NAL"), "abbreviations"][[1]], "nali"))
antimicrobials[which(antimicrobials$ab == "NEO"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "NEO"), "abbreviations"][[1]], "neom"))
antimicrobials[which(antimicrobials$ab == "NET"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "NET"), "abbreviations"][[1]], "neti"))
antimicrobials[which(antimicrobials$ab == "NIT"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "NIT"), "abbreviations"][[1]], "nitr"))
antimicrobials[which(antimicrobials$ab == "NOR"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "NOR"), "abbreviations"][[1]], "norf"))
antimicrobials[which(antimicrobials$ab == "NYS"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "NYS"), "abbreviations"][[1]], "nyst"))
antimicrobials[which(antimicrobials$ab == "OFX"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "OFX"), "abbreviations"][[1]], "oflo"))
antimicrobials[which(antimicrobials$ab == "OXA"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "OXA"), "abbreviations"][[1]], "oxal"))
antimicrobials[which(antimicrobials$ab == "PEF"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "PEF"), "abbreviations"][[1]], "pefl"))
antimicrobials[which(antimicrobials$ab == "PEN"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "PEN"), "abbreviations"][[1]], "peni"))
antimicrobials[which(antimicrobials$ab == "PIP"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "PIP"), "abbreviations"][[1]], "pipc"))
antimicrobials[which(antimicrobials$ab == "PPA"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "PPA"), "abbreviations"][[1]], "pipz"))
antimicrobials[which(antimicrobials$ab == "TZP"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "TZP"), "abbreviations"][[1]], "pita"))
antimicrobials[which(antimicrobials$ab == "PLB"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "PLB"), "abbreviations"][[1]], "polb"))
antimicrobials[which(antimicrobials$ab == "POS"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "POS"), "abbreviations"][[1]], "posa"))
antimicrobials[which(antimicrobials$ab == "PRI"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "PRI"), "abbreviations"][[1]], "pris"))
antimicrobials[which(antimicrobials$ab == "QDA"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "QDA"), "abbreviations"][[1]], "quda"))
antimicrobials[which(antimicrobials$ab == "RIF"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "RIF"), "abbreviations"][[1]], "rifa"))
antimicrobials[which(antimicrobials$ab == "RXT"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "RXT"), "abbreviations"][[1]], "roxi"))
antimicrobials[which(antimicrobials$ab == "SMX"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "SMX"), "abbreviations"][[1]], "sfmx"))
antimicrobials[which(antimicrobials$ab == "SLF4"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "SLF4"), "abbreviations"][[1]], "sfmz"))
antimicrobials[which(antimicrobials$ab == "SSS"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "SSS"), "abbreviations"][[1]], "sfna"))
antimicrobials[which(antimicrobials$ab == "SLF"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "SLF"), "abbreviations"][[1]], "sfsz"))
antimicrobials[which(antimicrobials$ab == "SPX"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "SPX"), "abbreviations"][[1]], "spar"))
antimicrobials[which(antimicrobials$ab == "SPT"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "SPT"), "abbreviations"][[1]], "spec"))
antimicrobials[which(antimicrobials$ab == "SPI"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "SPI"), "abbreviations"][[1]], "spir"))
antimicrobials[which(antimicrobials$ab == "STH"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "STH"), "abbreviations"][[1]], "sthl"))
antimicrobials[which(antimicrobials$ab == "STR1"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "STR1"), "abbreviations"][[1]], "stre"))
antimicrobials[which(antimicrobials$ab == "TAZ"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "TAZ"), "abbreviations"][[1]], "tazo"))
antimicrobials[which(antimicrobials$ab == "TEC"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "TEC"), "abbreviations"][[1]], "teic"))
antimicrobials[which(antimicrobials$ab == "TLT"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "TLT"), "abbreviations"][[1]], "teli"))
antimicrobials[which(antimicrobials$ab == "TMX"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "TMX"), "abbreviations"][[1]], "tema"))
antimicrobials[which(antimicrobials$ab == "TEM"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "TEM"), "abbreviations"][[1]], "temo"))
antimicrobials[which(antimicrobials$ab == "TRB"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "TRB"), "abbreviations"][[1]], "terb"))
antimicrobials[which(antimicrobials$ab == "TCY"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "TCY"), "abbreviations"][[1]], "tetr"))
antimicrobials[which(antimicrobials$ab == "TIC"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "TIC"), "abbreviations"][[1]], "tica"))
antimicrobials[which(antimicrobials$ab == "TCC"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "TCC"), "abbreviations"][[1]], "ticl"))
antimicrobials[which(antimicrobials$ab == "TGC"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "TGC"), "abbreviations"][[1]], "tige"))
antimicrobials[which(antimicrobials$ab == "TIN"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "TIN"), "abbreviations"][[1]], "tini"))
antimicrobials[which(antimicrobials$ab == "TOB"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "TOB"), "abbreviations"][[1]], "tobr"))
antimicrobials[which(antimicrobials$ab == "TOH"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "TOH"), "abbreviations"][[1]], "tohl"))
antimicrobials[which(antimicrobials$ab == "TMP"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "TMP"), "abbreviations"][[1]], "trim"))
antimicrobials[which(antimicrobials$ab == "TVA"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "TVA"), "abbreviations"][[1]], "trov"))
antimicrobials[which(antimicrobials$ab == "SLT4"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "SLT4"), "abbreviations"][[1]], "trsm"))
antimicrobials[which(antimicrobials$ab == "SXT"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "SXT"), "abbreviations"][[1]], "trsx"))
antimicrobials[which(antimicrobials$ab == "VAN"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "VAN"), "abbreviations"][[1]], "vanc"))
antimicrobials[which(antimicrobials$ab == "VOR"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "VOR"), "abbreviations"][[1]], "vori"))

antimicrobials[which(antimicrobials$ab == "FOS"), "synonyms"][[1]] <- list(sort(c(antimicrobials[which(antimicrobials$ab == "FOS"), "synonyms"][[1]], "Monuril")))
antimicrobials[which(antimicrobials$ab == "FOS"), "synonyms"][[1]] <- list(sort(c(antimicrobials[which(antimicrobials$ab == "FOS"), "synonyms"][[1]], "Monurol")))

antimicrobials[which(antimicrobials$ab == "TZP"), "abbreviations"][[1]] <- list(sort(c(antimicrobials[which(antimicrobials$ab == "TZP"), "abbreviations"][[1]], "piptazo")))

antimicrobials[which(antimicrobials$ab == "RFP"), "abbreviations"][[1]] <- list(sort(c(antimicrobials[which(antimicrobials$ab == "RFP"), "abbreviations"][[1]], "RPT")))
antimicrobials[which(antimicrobials$ab == "RTP"), "abbreviations"][[1]] <- list(sort(c(antimicrobials[which(antimicrobials$ab == "RTP"), "abbreviations"][[1]], "RET")))
antimicrobials[which(antimicrobials$ab == "TYL1"), "abbreviations"][[1]] <- list(sort(c(antimicrobials[which(antimicrobials$ab == "TYL1"), "abbreviations"][[1]], "TVN")))

antimicrobials <- antimicrobials %>%
  mutate(ab = as.character(ab)) %>%
  rbind(antimicrobials %>%
    filter(ab == "GEH") %>%
    mutate(
      ab = "AMH",
      name = "Amphotericin B-high",
      abbreviations = list(c("amhl", "amfo b high", "ampho b high", "amphotericin high"))
    )) %>%
  rbind(antimicrobials %>%
    filter(ab == "GEH") %>%
    mutate(
      ab = "TOH",
      name = "Tobramycin-high",
      abbreviations = list(c("tohl", "tobra high", "tobramycin high"))
    )) %>%
  rbind(antimicrobials %>%
    filter(ab == "BUT") %>%
    mutate(
      ab = "CIX",
      atc = "D01AE14",
      name = "Ciclopirox",
      group = "Antifungals/antimycotics",
      atc_group1 = "Antifungals for topical use",
      atc_group2 = "Other antifungals for topical use",
      abbreviations = list(c("cipx"))
    ))
antimicrobials[which(antimicrobials$ab == "SSS"), "name"] <- "Sulfonamide"
# ESBL E-test codes:
antimicrobials[which(antimicrobials$ab == "CCV"), "abbreviations"][[1]] <- list(c("xtzl"))
antimicrobials[which(antimicrobials$ab == "CAZ"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "CAZ"), "abbreviations"][[1]], "xtz", "cefta"))
antimicrobials[which(antimicrobials$ab == "CPC"), "abbreviations"][[1]] <- list(c("xpml"))
antimicrobials[which(antimicrobials$ab == "FEP"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "FEP"), "abbreviations"][[1]], "xpm"))
antimicrobials[which(antimicrobials$ab == "CTC"), "abbreviations"][[1]] <- list(c("xctl"))
antimicrobials[which(antimicrobials$ab == "CTX"), "abbreviations"][[1]] <- list(c(antimicrobials[which(antimicrobials$ab == "CTX"), "abbreviations"][[1]], "xct"))
# High level Gentamcin and Streptomycin
antimicrobials[which(antimicrobials$ab == "GEH"), "abbreviations"][[1]] <- list(c("gehl", "gentamicin high", "genta high", "gehi"))
antimicrobials[which(antimicrobials$ab == "STH"), "abbreviations"][[1]] <- list(c("sthl", "streptomycin high", "strepto high", "sthi"))
# add imi and "imipenem/cilastatine" to imipenem
antimicrobials[which(antimicrobials$ab == "IPM"), "abbreviations"][[1]] <- list(c("imip", "imi", "imp"))
antimicrobials[which(antimicrobials$ab == "IPM"), "synonyms"][[1]] <- list(sort(c(antimicrobials[which(antimicrobials$ab == "IPM"), "synonyms"][[1]], "imipenem/cilastatin")))
# add synonyms of ones not found
antimicrobials[which(antimicrobials$ab == "TZP"), "synonyms"][[1]] <- list(sort(c(antimicrobials[which(antimicrobials$ab == "TZP"), "synonyms"][[1]], "Tazocel", "tazocillin", "Tazocin", "Zosyn")))
antimicrobials[which(antimicrobials$ab == "COL"), "synonyms"][[1]] <- list(sort(unique(c(antimicrobials[which(antimicrobials$ab == "COL"), "synonyms"][[1]], "Colisticin", "Polymyxin E", "Colimycin", "Coly-Mycin", "Totazina", "Colistimethate", "Promixin", "Colistimethate Sodium"))))
# remove incorrect synonyms from rifampicin (RIF) and add them to the combination rifampicin/isoniazid (RFI)
old_sym <- antimicrobials[which(antimicrobials$ab == "RIF"), "synonyms"][[1]]
old_sym <- old_sym[!old_sym %in% c("Rifinah", "Rimactazid")]
antimicrobials[which(antimicrobials$ab == "RIF"), "synonyms"][[1]] <- list(old_sym)
antimicrobials[which(antimicrobials$ab == "RFI"), "synonyms"][[1]] <- list(sort(c("Rifinah", "Rimactazid")))
# remove incorrect synonyms from sulfamethoxazole (SMX) and add them to the combination trimethoprim/sulfamethoxazole (SXT)
old_sym <- antimicrobials[which(antimicrobials$ab == "SMX"), "synonyms"][[1]]
old_sym <- old_sym[!old_sym %in% c("Cotrimoxazole", "Bactrimel")]
antimicrobials[which(antimicrobials$ab == "SMX"), "synonyms"][[1]] <- list(old_sym)
antimicrobials[which(antimicrobials$ab == "SXT"), "synonyms"][[1]] <- list(sort(unique(c(antimicrobials[which(antimicrobials$ab == "COL"), "synonyms"][[1]], "Cotrimoxazole", "Bactrimel", "Septra", "Bactrim", "Cotrimazole"))))

# Fix penicillins
antimicrobials[which(antimicrobials$ab == "PEN"), "abbreviations"][[1]] <- list(c("bepe", "pg", "pen", "peni", "peni g", "penicillin", "penicillin g"))
antimicrobials[which(antimicrobials$ab == "PEN"), "name"] <- "Benzylpenicillin"
antimicrobials[which(antimicrobials$ab == "PHN"), "abbreviations"][[1]] <- list(c("fepe", "peni v", "pv", "penicillin v", "PNV"))
antimicrobials <- subset(antimicrobials, antimicrobials$ab != "PNV")

# New DDDs
antimicrobials[which(antimicrobials$ab == "PEN"), "iv_ddd"] <- 3.6
antimicrobials[which(antimicrobials$ab == "PEN"), "iv_units"] <- "g"

## new ATC codes
# ceftaroline
antimicrobials[which(antimicrobials$ab == "CPT"), "atc"] <- "J01DI02"
# faropenem
antimicrobials[which(antimicrobials$ab == "FAR"), "atc"] <- "J01DI03"
# ceftobiprole
antimicrobials[which(antimicrobials$ab == "BPR"), "atc"] <- "J01DI01"
# ceftazidime / avibactam
antimicrobials[which(antimicrobials$ab == "CZA"), "atc"] <- "J01DD52"
antimicrobials[which(antimicrobials$ab == "CZA"), "cid"] <- 90643431
antimicrobials[which(antimicrobials$ab == "CZA"), "atc_group1"] <- "Other beta-lactam antibacterials"
antimicrobials[which(antimicrobials$ab == "CZA"), "atc_group2"] <- "Third-generation cephalosporins"
antimicrobials[which(antimicrobials$ab == "CZA"), "iv_ddd"] <- 6
antimicrobials[which(antimicrobials$ab == "CZA"), "iv_units"] <- "g"
antimicrobials[which(antimicrobials$ab == "CZA"), "synonyms"] <- list(c("Avycaz", "Zavicefta"))

# typo
antimicrobials[which(antimicrobials$ab == "RXT"), "name"] <- "Roxithromycin"
antimicrobials[which(antimicrobials$ab == "PEN"), "atc"] <- "J01CE01"

# WHONET cleanup
antimicrobials[which(antimicrobials$ab == "BCZ"), "name"] <- "Bicyclomycin"
antimicrobials[which(antimicrobials$ab == "CCL"), "name"] <- "Cefetecol"
antimicrobials[which(antimicrobials$ab == "ENV"), "name"] <- "Enviomycin"
antimicrobials[which(antimicrobials$ab == "KIT"), "name"] <- "Kitasamycin"
antimicrobials[which(antimicrobials$ab == "LSP"), "name"] <- "Linco-spectin"
antimicrobials[which(antimicrobials$ab == "MEC"), "name"] <- "Mecillinam"
antimicrobials[which(antimicrobials$ab == "PMR"), "name"] <- "Pimaricin"
antimicrobials[which(antimicrobials$ab == "BCZ"), "abbreviations"][[1]] <- list(sort(unique(c(antimicrobials[which(antimicrobials$ab == "BCZ"), "abbreviations"][[1]], "Bicozamycin"))))
antimicrobials[which(antimicrobials$ab == "CCL"), "abbreviations"][[1]] <- list(sort(unique(c(antimicrobials[which(antimicrobials$ab == "CCL"), "abbreviations"][[1]], "Cefcatacol"))))
antimicrobials[which(antimicrobials$ab == "ENV"), "abbreviations"][[1]] <- list(sort(unique(c(antimicrobials[which(antimicrobials$ab == "ENV"), "abbreviations"][[1]], "Tuberactinomycin"))))
antimicrobials[which(antimicrobials$ab == "KIT"), "abbreviations"][[1]] <- list(sort(unique(c(antimicrobials[which(antimicrobials$ab == "KIT"), "abbreviations"][[1]], "Leucomycin"))))
antimicrobials[which(antimicrobials$ab == "LSP"), "abbreviations"][[1]] <- list(sort(unique(c(antimicrobials[which(antimicrobials$ab == "LSP"), "abbreviations"][[1]], "lincomycin/spectinomycin"))))
antimicrobials[which(antimicrobials$ab == "MEC"), "abbreviations"][[1]] <- list(sort(unique(c(antimicrobials[which(antimicrobials$ab == "MEC"), "abbreviations"][[1]], "Amdinocillin"))))
antimicrobials[which(antimicrobials$ab == "PMR"), "abbreviations"][[1]] <- list(sort(unique(c(antimicrobials[which(antimicrobials$ab == "PMR"), "abbreviations"][[1]], "Natamycin"))))


# set cephalosporins groups for the ones that could not be determined automatically:
antimicrobials <- antimicrobials %>%
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
    TRUE ~ group
  ))
antimicrobials[which(antimicrobials$ab %in% c("CYC", "LNZ", "THA", "TZD")), "group"] <- "Oxazolidinones"

# add efflux
effl <- antimicrobials |>
  filter(ab == "ACM") |>
  mutate(ab = as.character("EFF"),
         cid = NA_real_,
         name = "Efflux",
         group = "Other")
antimicrobials <- antimicrobials |>
  mutate(ab = as.character(ab)) |>
  bind_rows(effl)
class(antimicrobials$ab) <- c("ab", "character")
antimicrobials[which(antimicrobials$ab == "EFF"), "abbreviations"][[1]] <- list(c("effflux pump"))


# add clindamycin inducible screening
clin <- antimicrobials |>
  filter(ab == "FOX1") |>
  mutate(ab = as.character("CLI1"),
         name = "Clindamycin inducible screening",
         group = "Macrolides/lincosamides")
antimicrobials <- antimicrobials |>
  mutate(ab = as.character(ab)) |>
  bind_rows(clin)
class(antimicrobials$ab) <- c("ab", "character")
antimicrobials[which(antimicrobials$ab == "CLI1"), "abbreviations"][[1]] <- list(c("clindamycin inducible", "clinda inducible", "clin inducible"))

# add pretomanid
antimicrobials <- antimicrobials %>%
  mutate(ab = as.character(ab)) %>%
  bind_rows(antimicrobials %>%
    mutate(ab = as.character(ab)) %>%
    filter(ab == "SMF") %>%
    mutate(
      ab = "PMD",
      atc = "J04AK08",
      cid = 456199,
      name = "Pretomanid",
      abbreviations = list(""),
      oral_ddd = NA_real_
    ))



# update ATC codes from WHOCC website -------------------------------------

# last time checked: 2024-02-22

library(rvest)
updated_atc <- as.list(antimicrobials$atc)

get_atcs <- function(ab_name, type = "human") {
  if (type == "human") {
    url <- "https://atcddd.fhi.no/atc_ddd_index/"
  } else if (type == "veterinary") {
    url <- "https://atcddd.fhi.no/atcvet/atcvet_index/"  
  } else {
    stop("invalid type")
  }
  
  ab_name <- gsub("/", " and ", tolower(ab_name), fixed = TRUE)

  # we will do a search on their website, which means:

  # go to the url
  atc_tbl <- read_html(url) %>%
    # get all forms
    html_form() %>%
    # get the second form (the first form is a global website form)
    .[[2]] %>%
    # set the name input box to our search parameter
    html_form_set(name = ab_name) %>%
    # hit Submit
    html_form_submit() %>%
    # read the resulting page
    read_html() %>%
    # retrieve the table on it
    html_node("table") %>%
    # transform it to an R data set
    html_table(header = FALSE)
  
  # and get the ATCs (first column) of only exact hits
  unique(as.character(atc_tbl[which(tolower(atc_tbl[, 2, drop = TRUE]) == ab_name), 1, drop = TRUE]))
}

# this takes around 4 minutes (some are skipped and go faster)
for (i in seq_len(nrow(antimicrobials))) {
  message(percentage(i / nrow(antimicrobials), digits = 1),
    " - Downloading ", antimicrobials$name[i],
    appendLF = FALSE
  )
  atcs <- get_atcs(antimicrobials$name[i], type = "human")
  if (all(is.na(atcs))) {
    atcs <- get_atcs(antimicrobials$name[i], type = "veterinary")
  }
  if (length(atcs) > 0) {
    updated_atc[[i]] <- atcs
    message(" (", length(atcs), " results)")
    # let the WHO server rest for a second - they might have a limitation on the queries per second
    Sys.sleep(1)
  } else {
    message(" (skipping)")
  }
}

antimicrobials$atc <- updated_atc

# update DDDs from WHOCC website ------------------------------------------

# last time checked: 2024-02-22
ddd_oral <- rep(NA_real_, nrow(antimicrobials))
ddd_oral_units <- rep(NA_character_, nrow(antimicrobials))
ddd_iv <- rep(NA_real_, nrow(antimicrobials))
ddd_iv_units <- rep(NA_character_, nrow(antimicrobials))
progress <- progress_ticker(nrow(antimicrobials))
for (i in seq_len(nrow(antimicrobials))) {
  on.exit(close(progress))
  progress$tick()
  atcs <- antimicrobials$atc[[i]]
  if (!all(is.na(atcs))) {
    for (j in seq_len(length(atcs))) {
      # oral
      if (is.na(ddd_oral[i])) {
        ddd_oral[i] <- atc_online_ddd(atcs[j], administration = "O")
        if (!is.na(ddd_oral[i])) {
          ddd_oral_units[i] <- atc_online_ddd_units(atcs[j], administration = "O")
        }
      }
      # parenteral
      if (is.na(ddd_iv[i])) {
        ddd_iv[i] <- atc_online_ddd(atcs[j], administration = "P")
        if (!is.na(ddd_iv[i])) {
          ddd_iv_units[i] <- atc_online_ddd_units(atcs[j], administration = "P")
        }
      }
    }
  }
  if (!is.na(ddd_oral[i]) | !is.na(ddd_iv[i])) {
    # let the WHO server rest for 0.25 second - they might have a limitation on the queries per second
    Sys.sleep(0.25)
  }
}

antimicrobials$oral_ddd <- ddd_oral
antimicrobials$oral_units <- ddd_oral_units
antimicrobials$iv_ddd <- ddd_iv
antimicrobials$iv_units <- ddd_iv_units

# Wrap up -----------------------------------------------------------------

# set as data.frame again
antimicrobials <- dataset_UTF8_to_ASCII(as.data.frame(antimicrobials, stringsAsFactors = FALSE))
class(antimicrobials$ab) <- c("ab", "character")
antimicrobials <- dplyr::arrange(antimicrobials, name)

# REFER TO data-raw/loinc.R FOR ADDING LOINC CODES

# make all abbreviations and synonyms lower case, unique and alphabetically sorted ----
for (i in 1:nrow(antimicrobials)) {
  abb <- as.character(sort(unique(tolower(antimicrobials[i, "abbreviations", drop = TRUE][[1]]))))
  abb <- abb[abb != "" & abb %unlike% ":"]
  syn <- as.character(sort(unique(tolower(unname(unlist(antimicrobials[i, "synonyms", drop = TRUE]))))))
  syn <- gsub("[^a-z]", "", syn)
  syn <- gsub(" +", " ", syn)
  pharm_terms <- "(pa?ediatric|injection|oral|inhale|otic|sulfate|sulphate|sodium|base|anhydrous|anhydrate|stearate|syrup|natrium|hydrate|x?hcl|gsalt|vet[.]?)"
  syn <- gsub(paste0(" ", pharm_terms, "$"), "", syn)
  syn <- gsub(paste0("^", pharm_terms, " "), "", syn)
  syn <- trimws(syn)
  syn <- gsub(" [a-z]{1,3}$", "", syn, perl = TRUE)
  syn <- trimws(syn)
  syn <- syn[syn != "" & syn %unlike% ":" & !syn %in% tolower(antimicrobials$name)]
  syn <- unique(syn)
  # special cases
  if (antimicrobials$ab[i] == "VAN") syn <- syn[syn %unlike% "^tei?ch?o"]
  if (antimicrobials$ab[i] == "CLR") syn <- syn[syn %unlike% "^ery"]
  antimicrobials[i, "abbreviations"][[1]] <- ifelse(length(abb) == 0, list(""), list(abb))
  antimicrobials[i, "synonyms"][[1]] <- ifelse(length(syn) == 0, list(""), list(syn))
  if ("loinc" %in% colnames(antimicrobials)) {
    loinc <- as.character(sort(unique(tolower(antimicrobials[i, "loinc", drop = TRUE][[1]]))))
    loinc <- loinc[loinc != ""]
    antimicrobials[i, "loinc"][[1]] <- ifelse(length(loinc) == 0, list(""), list(loinc))
  }
}


usethis::use_data(antimicrobials, overwrite = TRUE, version = 2, compress = "xz")
rm(antimicrobials)
