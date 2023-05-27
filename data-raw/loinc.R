# ==================================================================== #
# TITLE                                                                #
# AMR: An R Package for Working with Antimicrobial Resistance Data     #
#                                                                      #
# SOURCE                                                               #
# https://github.com/msberends/AMR                                     #
#                                                                      #
# CITE AS                                                              #
# Berends MS, Luz CF, Friedrich AW, Sinha BNM, Albers CJ, Glasner C    #
# (2022). AMR: An R Package for Working with Antimicrobial Resistance  #
# Data. Journal of Statistical Software, 104(3), 1-31.                 #
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
# 2. Download the CSV from https://loinc.org/download/loinc-complete/ (Loinc_2.67_Text_2.67.zip)
# 3. Read Loinc.csv that's in zip folder LoincTable
loinc_df <- read.csv("data-raw/Loinc.csv",
  row.names = NULL,
  stringsAsFactors = FALSE
)

# 4. Clean and add
library(dplyr)
library(cleaner)
library(AMR)
loinc_df %>% freq(CLASS) # to find the drugs
loinc_df <- loinc_df %>% filter(CLASS == "DRUG/TOX")
ab_names <- antibiotics %>%
  pull(name) %>%
  paste0(collapse = "|") %>%
  paste0("(", ., ")")

antibiotics$loinc <- as.list(rep(NA_character_, nrow(antibiotics)))
for (i in seq_len(nrow(antibiotics))) {
  message(i)
  loinc_ab <- loinc_df %>%
    filter(COMPONENT %like% paste0("^", antibiotics$name[i])) %>%
    pull(LOINC_NUM)
  if (length(loinc_ab) > 0) {
    antibiotics$loinc[i] <- list(loinc_ab)
  }
}
# sort and fix for empty values
for (i in 1:nrow(antibiotics)) {
  loinc <- as.character(sort(unique(tolower(antibiotics[i, "loinc"][[1]]))))
  antibiotics[i, "loinc"][[1]] <- ifelse(length(loinc[!loinc == ""]) == 0, list(""), list(loinc))
}

# remember to update R/aa_globals.R for the documentation

dim(antibiotics) # for R/data.R
usethis::use_data(antibiotics, overwrite = TRUE)
rm(antibiotics)
