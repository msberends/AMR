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

# last updated: 20 January 2020 - Loinc_2.67

# Steps to reproduce:
# 1. Create a fake account at https://loinc.org (sad you have to create one...)
# 2. Download the CSV from https://loinc.org/download/loinc-table-file-csv/ (Loinc_2.67_Text_2.67.zip)
# 3. Read Loinc.csv that's in this zip file
loinc_df <- read.csv("data-raw/Loinc.csv",
                     row.names = NULL,
                     stringsAsFactors = FALSE)

# 4. Clean and add
library(dplyr)
library(cleaner)
library(AMR)
loinc_df %>% freq(CLASS) # to find the drugs
loinc_df <- loinc_df %>% filter(CLASS == "DRUG/TOX")
ab_names <- antibiotics %>% pull(name) %>% paste0(collapse = "|") %>% paste0("(", ., ")")

antibiotics$loinc <- as.list(rep(NA_character_, nrow(antibiotics)))
for (i in seq_len(nrow(antibiotics))) {
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
  antibiotics[i, "loinc"][[1]] <- ifelse(length(syn[!syn == ""]) == 0, list(""), list(loinc))
}

dim(antibiotics) # for R/data.R
usethis::use_data(antibiotics, overwrite = TRUE)
rm(antibiotics)
