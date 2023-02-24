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
# doi:10.18637/jss.v104.i03                                            #
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

# This script runs in under a minute and renews all guidelines of CLSI and EUCAST!
# Run it with source("data-raw/reproduction_of_clinical_breakpoints.R")

library(dplyr)
library(readr)
library(tidyr)
library(AMR)

# Install the WHONET 2022 software on Windows (http://www.whonet.org/software.html),
# and copy the folder C:\WHONET\Resources to the data-raw/WHONET/ folder
# (for ASIARS-Net update, also copy C:\WHONET\Codes to the data-raw/WHONET/ folder)

# Load source data ----
whonet_organisms <- read_tsv("data-raw/WHONET/Resources/Organisms.txt", na = c("", "NA", "-"), show_col_types = FALSE) %>%
  # remove old taxonomic names
  filter(TAXONOMIC_STATUS == "C") %>%
  transmute(ORGANISM_CODE = tolower(WHONET_ORG_CODE), ORGANISM) %>%
  # what's wrong here? 'sau' is both S. areus and S. aureus sp. aureus
  mutate(
    ORGANISM = if_else(ORGANISM_CODE == "sau", "Staphylococcus aureus", ORGANISM),
    ORGANISM = if_else(ORGANISM_CODE == "pam", "Pasteurella multocida", ORGANISM)
  )
whonet_breakpoints <- read_tsv("data-raw/WHONET/Resources/Breakpoints.txt", na = c("", "NA", "-"), show_col_types = FALSE) %>%
  filter(BREAKPOINT_TYPE == "Human", GUIDELINES %in% c("CLSI", "EUCAST"))
whonet_antibiotics <- read_tsv("data-raw/WHONET/Resources/Antibiotics.txt", na = c("", "NA", "-"), show_col_types = FALSE) %>%
  arrange(WHONET_ABX_CODE) %>%
  distinct(WHONET_ABX_CODE, .keep_all = TRUE)


# Transform data ----

whonet_organisms <- whonet_organisms %>%
  bind_rows(data.frame(
    ORGANISM_CODE = c("ebc", "cof"),
    ORGANISM = c("Enterobacterales", "Campylobacter")
  ))

breakpoints <- whonet_breakpoints %>%
  mutate(ORGANISM_CODE = tolower(ORGANISM_CODE)) %>%
  left_join(whonet_organisms) %>%
  filter(ORGANISM %unlike% "(^cdc |Gram.*variable|virus)")
# this ones lack a MO name, they will become "UNKNOWN":
breakpoints %>%
  filter(is.na(ORGANISM)) %>%
  pull(ORGANISM_CODE) %>%
  unique()


# Generate new lookup table for microorganisms ----

new_mo_codes <- breakpoints %>%
  distinct(ORGANISM_CODE, ORGANISM) %>%
  mutate(ORGANISM = ORGANISM %>%
    gsub("Issatchenkia orientalis", "Candida krusei", .) %>%
    gsub(", nutritionally variant", "", .) %>%
    gsub(", toxin-.*producing", "", .)) %>%
  mutate(
    mo = as.mo(ORGANISM, language = NULL, keep_synonyms = FALSE),
    mo_name = mo_name(mo, language = NULL)
  )


# Update microorganisms.codes with the latest WHONET codes ----

# these will be changed :
new_mo_codes %>%
  mutate(code = toupper(ORGANISM_CODE)) %>%
  rename(mo_new = mo) %>%
  left_join(microorganisms.codes %>% rename(mo_old = mo)) %>%
  filter(mo_old != mo_new)

microorganisms.codes <- microorganisms.codes %>%
  filter(!code %in% toupper(new_mo_codes$ORGANISM_CODE)) %>%
  bind_rows(new_mo_codes %>% transmute(code = toupper(ORGANISM_CODE), mo = mo) %>% filter(!is.na(mo))) %>%
  arrange(code) %>%
  as_tibble()
usethis::use_data(microorganisms.codes, overwrite = TRUE, compress = "xz", version = 2)
rm(microorganisms.codes)
devtools::load_all()

# update ASIARS-Net?
asiarsnet <- read_tsv("data-raw/WHONET/Codes/ASIARS_Net_Organisms_ForwardLookup.txt")
asiarsnet <- asiarsnet %>%
  mutate(WHONET_Code = toupper(WHONET_Code)) %>%
  left_join(whonet_organisms %>% mutate(WHONET_Code = toupper(ORGANISM_CODE))) %>%
  mutate(
    mo1 = as.mo(ORGANISM_CODE),
    mo2 = as.mo(ORGANISM)
  ) %>%
  mutate(mo = if_else(mo2 == "UNKNOWN" | is.na(mo2), mo1, mo2)) %>%
  filter(!is.na(mo))
insert1 <- asiarsnet %>% transmute(code = WHONET_Code, mo)
insert2 <- asiarsnet %>% transmute(code = as.character(ASIARS_Net_Code), mo)
# these will be updated
bind_rows(insert1, insert2) %>%
  rename(mo_new = mo) %>%
  left_join(microorganisms.codes) %>%
  filter(mo != mo_new)
microorganisms.codes <- microorganisms.codes %>%
  filter(!code %in% c(insert1$code, insert2$code)) %>%
  bind_rows(insert1, insert2) %>%
  arrange(code)


# Create new breakpoint table ----

breakpoints_new <- breakpoints %>%
  # only last 10 years
  filter(YEAR > as.double(format(Sys.Date(), "%Y")) - 10) %>%
  # "all" and "gen" (general) must become UNKNOWNs:
  mutate(ORGANISM_CODE = if_else(ORGANISM_CODE %in% c("all", "gen"), "UNKNOWN", ORGANISM_CODE)) %>%
  transmute(
    guideline = paste(GUIDELINES, YEAR),
    method = TEST_METHOD,
    site = gsub("Urinary tract infection", "UTI", SITE_OF_INFECTION),
    mo = as.mo(ORGANISM_CODE, keep_synonyms = FALSE),
    rank_index = case_when(
      mo_rank(mo) %like% "(infra|sub)" ~ 1,
      mo_rank(mo) == "species" ~ 2,
      mo_rank(mo) == "genus" ~ 3,
      mo_rank(mo) == "family" ~ 4,
      mo_rank(mo) == "order" ~ 5,
      TRUE ~ 6
    ),
    ab = as.ab(WHONET_ABX_CODE),
    ref_tbl = REFERENCE_TABLE,
    disk_dose = POTENCY,
    breakpoint_S = S,
    breakpoint_R = R,
    uti = SITE_OF_INFECTION %like% "(UTI|urinary|urine)"
  ) %>%
  # Greek symbols and EM dash symbols are not allowed by CRAN, so replace them with ASCII:
  mutate(disk_dose = disk_dose %>%
    gsub("μ", "u", ., fixed = TRUE) %>%
    gsub("µ", "u", ., fixed = TRUE) %>% # this is another micro sign, although we cannot see it
    gsub("–", "-", ., fixed = TRUE)) %>%
  arrange(desc(guideline), ab, mo, method) %>%
  filter(!(is.na(breakpoint_S) & is.na(breakpoint_R)) & !is.na(mo) & !is.na(ab)) %>%
  distinct(guideline, ab, mo, method, site, breakpoint_S, .keep_all = TRUE)

# clean disk zones and MICs
breakpoints_new[which(breakpoints_new$method == "DISK"), "breakpoint_S"] <- as.double(as.disk(breakpoints_new[which(breakpoints_new$method == "DISK"), "breakpoint_S", drop = TRUE]))
breakpoints_new[which(breakpoints_new$method == "DISK"), "breakpoint_R"] <- as.double(as.disk(breakpoints_new[which(breakpoints_new$method == "DISK"), "breakpoint_R", drop = TRUE]))

# WHONET has no >1024 but instead uses 1025, 513, etc, so as.mic() cannot be used to clean.
# instead, clean based on MIC factor levels
m <- unique(as.double(as.mic(levels(as.mic(1)))))
breakpoints_new[which(breakpoints_new$method == "MIC" &
  is.na(breakpoints_new$breakpoint_S)), "breakpoint_S"] <- min(m)
breakpoints_new[which(breakpoints_new$method == "MIC" &
  is.na(breakpoints_new$breakpoint_R)), "breakpoint_R"] <- max(m)
# raise these one higher valid MIC factor level:
breakpoints_new[which(breakpoints_new$breakpoint_R == 129), "breakpoint_R"] <- m[which(m == 128) + 1]
breakpoints_new[which(breakpoints_new$breakpoint_R == 257), "breakpoint_R"] <- m[which(m == 256) + 1]
breakpoints_new[which(breakpoints_new$breakpoint_R == 513), "breakpoint_R"] <- m[which(m == 512) + 1]
breakpoints_new[which(breakpoints_new$breakpoint_R == 1025), "breakpoint_R"] <- m[which(m == 1024) + 1]

# WHONET adds one log2 level to the R breakpoint for their software, e.g. in AMC in Enterobacterales:
# EUCAST 2021 guideline: S <= 8 and R > 8
#           WHONET file: S <= 8 and R >= 16
breakpoints_new %>% filter(guideline == "EUCAST 2022", ab == "AMC", mo == "B_[ORD]_ENTRBCTR", method == "MIC")
# this will make an MIC of 12 I, which should be R, so:
breakpoints_new <- breakpoints_new %>%
  mutate(breakpoint_R = ifelse(guideline %like% "EUCAST" & method == "MIC" & log2(breakpoint_R) - log2(breakpoint_S) != 0,
    pmax(breakpoint_S, breakpoint_R / 2),
    breakpoint_R
  ))
# fix disks as well
breakpoints_new <- breakpoints_new %>%
  mutate(breakpoint_R = ifelse(guideline %like% "EUCAST" & method == "DISK" & breakpoint_S - breakpoint_R != 0,
    breakpoint_R + 1,
    breakpoint_R
  ))
# fix missing R breakpoint where there is an S breakpoint
breakpoints_new[which(is.na(breakpoints_new$breakpoint_R)), "breakpoint_R"] <- breakpoints_new[which(is.na(breakpoints_new$breakpoint_R)), "breakpoint_S"]

# check again
breakpoints_new %>% filter(guideline == "EUCAST 2022", ab == "AMC", mo == "B_[ORD]_ENTRBCTR", method == "MIC")
# compare with current version
clinical_breakpoints %>% filter(guideline == "EUCAST 2022", ab == "AMC", mo == "B_[ORD]_ENTRBCTR", method == "MIC")

# Save to package ----

clinical_breakpoints <- breakpoints_new
usethis::use_data(clinical_breakpoints, overwrite = TRUE, compress = "xz", version = 2)
rm(clinical_breakpoints)
devtools::load_all(".")
