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

# last updated: 30 October 2022 - Loinc_2.73

# Steps to reproduce:
# 1. Create a fake account at https://loinc.org (sad you have to create one...)
# 2. Download the CSV from https://loinc.org/download/loinc-complete/ 
# 3. Read file LoincTable/Loinc.csv
loinc_df <- read.csv("data-raw/Loinc.csv",
  row.names = NULL,
  stringsAsFactors = FALSE
)

# 4. Clean and add
library(dplyr)
library(cleaner)
library(AMR)
# to find the drugs:
loinc_df %>%
  filter(COMPONENT %like% "ampicillin|fluconazol|meropenem") %>%
  count(CLASS, sort = TRUE)
loinc_df <- loinc_df %>%
  filter(CLASS %in% c("DRUG/TOX", "ABXBACT")) %>% 
  mutate(name = generalise_antibiotic_name(COMPONENT), .before = 1)

# antibiotics
antibiotics$loinc <- as.list(rep(NA_character_, nrow(antibiotics)))
for (i in seq_len(nrow(antibiotics))) {
  message(i)
  loinc_ab <- loinc_df %>%
    filter(name %like% paste0("^", generalise_antibiotic_name(antibiotics$name[i]))) %>%
    pull(LOINC_NUM)
  if (length(loinc_ab) > 0) {
    antibiotics$loinc[i] <- list(loinc_ab)
  }
}

# antivirals
antivirals$loinc <- as.list(rep(NA_character_, nrow(antivirals)))
for (i in seq_len(nrow(antivirals))) {
  message(i)
  loinc_ab <- loinc_df %>%
    filter(name %like% paste0("^", generalise_antibiotic_name(antivirals$name[i]))) %>%
    pull(LOINC_NUM)
  if (length(loinc_ab) > 0) {
    antivirals$loinc[i] <- list(loinc_ab)
  }
}

# sort and fix for empty values
for (i in 1:nrow(antibiotics)) {
  loinc <- as.character(sort(unique(tolower(antibiotics[i, "loinc", drop = TRUE][[1]]))))
  loinc <- loinc[loinc != ""]
  antibiotics[i, "loinc"][[1]] <- ifelse(length(loinc) == 0, list(""), list(loinc))
}
for (i in 1:nrow(antivirals)) {
  loinc <- as.character(sort(unique(tolower(antivirals[i, "loinc", drop = TRUE][[1]]))))
  loinc <- loinc[loinc != ""]
  antivirals[i, "loinc"][[1]] <- ifelse(length(loinc) == 0, list(""), list(loinc))
}

antibiotics <- dataset_UTF8_to_ASCII(as.data.frame(antibiotics, stringsAsFactors = FALSE))
antibiotics <- dplyr::arrange(antibiotics, name)

antivirals <- dataset_UTF8_to_ASCII(as.data.frame(antivirals, stringsAsFactors = FALSE))
antivirals <- dplyr::arrange(antivirals, name)

# remember to update R/aa_globals.R for the documentation

dim(antibiotics) # for R/data.R
usethis::use_data(antibiotics, internal = FALSE, overwrite = TRUE, compress = "xz", version = 2)
rm(antibiotics)

dim(antivirals) # for R/data.R
usethis::use_data(antivirals, internal = FALSE, overwrite = TRUE, compress = "xz", version = 2)
rm(antivirals)
