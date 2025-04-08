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
# how to conduct AMR data analysis: https://amr-for-r.org/             #
# ==================================================================== #

library(dplyr)
library(tidyr)
library(rvest)

# get all data from the WHOCC website
get_atc_table <- function(atc_group) {
  # give as input J0XXX, like atc_group = "J05AB"
  downloaded <- read_html(paste0("https://atcddd.fhi.no/atc_ddd_index/?code=", atc_group, "&showdescription=no"))
  table_title <- downloaded %>%
    html_nodes(paste0('a[href^="./?code=', atc_group, '&"]')) %>%
    html_text()
  table_title <- table_title[tolower(table_title) != "show text from guidelines"][1]
  table_content <- downloaded %>%
    html_nodes("table") %>%
    html_table(header = TRUE) %>%
    # returns list, so make data.frame out of it
    as.data.frame(stringsAsFactors = FALSE) %>%
    # select right columns
    select(atc = ATC.code, name = Name, ddd = DDD, unit = U, ddd_type = Adm.R) %>%
    # fill empty rows
    mutate(atc = ifelse(atc == "", lag(atc), atc), name = ifelse(name == "", lag(name), name)) %>%
    pivot_wider(names_from = ddd_type, values_from = c(ddd, unit)) %>%
    mutate(atc_group = table_title)
  if (!"ddd_O" %in% colnames(table_content)) {
    table_content <- table_content %>% mutate(ddd_O = NA_real_, unit_O = NA_character_)
  }
  if (!"ddd_P" %in% colnames(table_content)) {
    table_content <- table_content %>% mutate(ddd_P = NA_real_, unit_P = NA_character_)
  }
  table_content %>% select(atc, name, atc_group,
    oral_ddd = ddd_O, oral_units = unit_O,
    iv_ddd = ddd_P, iv_units = unit_P
  )
}

# these are the relevant groups for input: https://atcddd.fhi.no/atc_ddd_index/?code=J05A (J05 only contains J05A)
atc_groups <- c("J05AA", "J05AB", "J05AC", "J05AD", "J05AE", "J05AF", "J05AG", "J05AH", "J05AJ", "J05AP", "J05AR", "J05AX")

# get the first
antivirals <- get_atc_table(atc_groups[1])
# bind all others to it
for (i in 2:length(atc_groups)) {
  message(atc_groups[i], "...")
  antivirals <- rbind(antivirals, get_atc_table(atc_groups[i]))
}

# arrange on name, untibble it
antivirals <- antivirals %>%
  arrange(name) %>%
  as.data.frame(stringsAsFactors = FALSE)

# add PubChem Compound ID (cid) and their trade names
# see `data-raw/_reproduction_scripts/reproduction_of_antimicrobials.R` for get_CID() and get_synonyms()
CIDs <- get_CID(antivirals$name)
# these could not be found:
antivirals[is.na(CIDs), ] %>% View()
# get brand names from PubChem
synonyms <- get_synonyms(CIDs)
synonyms <- lapply(
  synonyms,
  function(x) {
    if (length(x) == 0 | all(is.na(x))) {
      ""
    } else {
      x
    }
  }
)

antivirals <- antivirals %>%
  transmute(atc,
    cid = as.double(CIDs),
    name,
    atc_group,
    synonyms = unname(synonyms),
    oral_ddd,
    oral_units,
    iv_ddd,
    iv_units
  ) %>%
  AMR:::dataset_UTF8_to_ASCII()

av_codes <- tibble(name = antivirals$name %>%
  strsplit("(, | and )") %>%
  unlist() %>%
  unique() %>%
  sort()) %>%
  mutate(av_1st = toupper(abbreviate(name, minlength = 3, use.classes = FALSE))) %>%
  filter(!name %in% c("acid", "dipivoxil", "disoproxil", "marboxil", "alafenamide"))

replace_with_av_code <- function(name) {
  unname(av_codes$av_1st[match(name, av_codes$name)])
}

names_codes <- antivirals %>%
  separate(name,
    into = paste0("name", c(1:7)),
    sep = "(, | and )",
    remove = FALSE,
    fill = "right"
  ) %>%
  # remove empty columns
  select(!where(function(x) all(is.na(x)))) %>%
  mutate_at(vars(matches("name[1-9]")), replace_with_av_code) %>%
  unite(av, matches("name[1-9]"), sep = "+", na.rm = TRUE) %>%
  mutate(name = gsub("(, | and )", "/", name))
substr(names_codes$name, 1, 1) <- toupper(substr(names_codes$name, 1, 1))

antivirals <- bind_cols(
  names_codes %>% select(av, name),
  antivirals %>% select(-name)
)
class(antivirals$av) <- c("av", "character")
antivirals <- antivirals %>% AMR:::dataset_UTF8_to_ASCII()

# ! add loinc, run 'data-raw/loinc.R' !

# de-duplicate synonyms
for (i in 1:nrow(antivirals)) {
  syn <- as.character(sort(unique(tolower(antivirals[i, "synonyms", drop = TRUE][[1]]))))
  syn <- syn[!syn %in% tolower(antivirals[i, "name", drop = TRUE])]
  antivirals[i, "synonyms"][[1]] <- ifelse(length(syn[!syn == ""]) == 0, list(""), list(syn))
}

antivirals <- antivirals %>% AMR:::dataset_UTF8_to_ASCII()

# check it
antivirals

# save it
usethis::use_data(antivirals, overwrite = TRUE, internal = FALSE, compress = "xz", version = 2)
