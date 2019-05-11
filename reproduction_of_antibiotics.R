
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

# get brand names (2-3 min)
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
            name = gsub("edta", "EDTA", gsub("/([A-Z])", "/\\L\\1", name, perl = TRUE), ignore.case = TRUE),
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
# 'Polymixin B' (POL) and 'Polymyxin B' (PLB) both exist, so:
antibiotics[which(antibiotics$ab == "PLB"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "PLB"), "abbreviations"][[1]], "POL", "Polymixin", "Polymixin B"))
antibiotics <- filter(antibiotics, ab != "POL")
# 'Latamoxef' (LTM) and 'Moxalactam (Latamoxef)' (MOX) both exist, so:
antibiotics[which(antibiotics$ab == "LTM"), "abbreviations"][[1]] <- list(c("MOX", "moxa"))
antibiotics <- filter(antibiotics, ab != "MOX")
# ESBL E-test codes:
antibiotics[which(antibiotics$ab == "CCV"), "abbreviations"][[1]] <- list(c("xtzl"))
antibiotics[which(antibiotics$ab == "CAZ"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "CAZ"), "abbreviations"][[1]], "xtz", "cefta"))
antibiotics[which(antibiotics$ab == "CPC"), "abbreviations"][[1]] <- list(c("xpml"))
antibiotics[which(antibiotics$ab == "FEP"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "FEP"), "abbreviations"][[1]], "xpm"))
antibiotics[which(antibiotics$ab == "CTC"), "abbreviations"][[1]] <- list(c("xctl"))
antibiotics[which(antibiotics$ab == "CTX"), "abbreviations"][[1]] <- list(c(antibiotics[which(antibiotics$ab == "CTX"), "abbreviations"][[1]], "xct"))

class(antibiotics$ab) <- "ab"
class(antibiotics$atc) <- "atc"

usethis::use_data(antibiotics, overwrite = TRUE)
rm(antibiotics)
