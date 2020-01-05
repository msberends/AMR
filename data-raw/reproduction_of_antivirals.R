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

# get all data from the WHOCC website

get_atc_table <- function(atc_group) {
  # give as input J0XXX, like atc_group = "J05AB"
  downloaded <- read_html(paste0("https://www.whocc.no/atc_ddd_index/?code=", atc_group, "&showdescription=no"))
  table_title <- downloaded %>% html_nodes(paste0('a[href="./?code=', atc_group, '"]')) %>% html_text()
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
                           iv_ddd = ddd_P, iv_units = unit_P)
}

# these are the relevant groups for input: https://www.whocc.no/atc_ddd_index/?code=J05A (J05 only contains J05A)
atc_groups <- c("J05AA", "J05AB", "J05AC", "J05AD", "J05AE", "J05AF", "J05AG", "J05AH", "J05AP", "J05AR", "J05AX")

# get the first
antivirals <- get_atc_table(atc_groups[1])
# bind all others to it
for (i in 2:length(atc_groups)) {
  antivirals <- rbind(antivirals, get_atc_table(atc_groups[i]))
}

# arrange on name, untibble it
antivirals <- antivirals %>% arrange(name) %>% as.data.frame(stringsAsFactors = FALSE)

# add PubChem Compound ID (cid) and their trade names - functions are in file to create `antibiotics` data set
CIDs <- get_CID(antivirals$name)
# these could not be found:
antivirals[is.na(CIDs),] %>% View()
# get brand names from PubChem 
synonyms <- get_synonyms(CIDs)
synonyms <- lapply(synonyms,
                   function(x) {
                     if (length(x) == 0 | all(is.na(x))) {
                       ""
                     } else {
                       x
                     }})

antivirals <- antivirals %>%
  transmute(atc,
            cid = CIDs,
            name,
            atc_group,
            synonyms = unname(synonyms),
            oral_ddd,
            oral_units,
            iv_ddd,
            iv_units)

# save it
usethis::use_data(antivirals, overwrite = TRUE)
