# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Data Analysis for R                   #
#                                                                      #
# SOURCE                                                               #
# https://github.com/msberends/AMR                                     #
#                                                                      #
# LICENCE                                                              #
# (c) 2018-2022 Berends MS, Luz CF et al.                              #
# Developed at the University of Groningen, the Netherlands, in        #
# collaboration with non-profit organisations Certe Medical            #
# Diagnostics & Advice, and University Medical Center Groningen.       # 
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

library(AMR)
library(tidyverse)

# we will use Public Health Information Network Vocabulary Access and Distribution System (PHIN VADS)
# as a source, which copies directly from the latest US SNOMED CT version
# - go to https://phinvads.cdc.gov/vads/ViewValueSet.action?oid=2.16.840.1.114222.4.11.1009
# - check that current online version is higher than SNOMED_VERSION$current_version
# - if so, click on 'Download Value Set', choose 'TXT'
snomed <- read_tsv("data-raw/SNOMED_PHVS_Microorganism_CDC_V12.txt", skip = 3) %>% 
  select(1:2) %>% 
  set_names(c("snomed", "mo"))

# save all valid genera, species and subspecies
vctr <- unique(unlist(strsplit(c(microorganisms$fullname, microorganisms.old$fullname), " ")))
vctr <- tolower(vctr[vctr %like% "^[a-z]+$"])

# remove all parts of the name that are no valid values in genera, species or subspecies
# this takes ~20 seconds
snomed <- snomed %>% 
  mutate(fullname = vapply(FUN.VALUE = character(1),
                           # split on space and/or comma
                           strsplit(tolower(mo), "[ ,]"),
                           function(x) trimws(paste0(x[x %in% vctr], collapse = " "))),
         # remove " group"
         fullname = gsub(" group", "", fullname, fixed = TRUE))

snomed_keep <- snomed %>% 
  filter(fullname %in% tolower(c(microorganisms$fullname, microorganisms.old$fullname))) %>% 
  group_by(fullname_lower = fullname) %>% 
  summarise(snomed = list(snomed))

message(nrow(snomed_keep), " MO's will get a SNOMED code.")

# save to microorganisms data set
microorganisms <- microorganisms %>%
  # remove old snomed
  select(-snomed) %>%
  # create dummy var for joining
  mutate(fullname_lower = tolower(fullname)) %>%
  # join new snomed
  left_join(snomed_keep) %>%
  # remove dummy var
  select(-fullname_lower) %>% 
  AMR:::dataset_UTF8_to_ASCII()

# don't forget to update the version number in SNOMED_VERSION in ./R/globals.R!

# usethis::use_data(microorganisms, overwrite = TRUE, version = 2, compress = "xz")

