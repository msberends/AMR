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

# This script runs in under a minute and renews all guidelines of CLSI and EUCAST!

library(dplyr)
library(readr)
library(tidyr)
library(AMR)

# Install the WHONET software on Windows (http://www.whonet.org/software.html),
# and copy the folder C:\WHONET\Codes to data-raw/WHONET/Codes
DRGLST <- readr::read_tsv("data-raw/WHONET/Codes/DRGLST.txt", na = c("", "NA", "-"))
DRGLST1 <- readr::read_tsv("data-raw/WHONET/Codes/DRGLST1.txt", na = c("", "NA", "-"))
ORGLIST <- readr::read_tsv("data-raw/WHONET/Codes/ORGLIST.txt", na = c("", "NA", "-"))

# create data set for generic rules (i.e., AB-specific but not MO-specific)
rsi_generic <- DRGLST %>%
  filter(CLSI == "X" | EUCST == "X") %>%
  select(ab = ANTIBIOTIC, disk_dose = POTENCY, matches("^(CLSI|EUCST)[0-9]")) %>% 
  mutate(ab = as.ab(ab),
         across(matches("(CLSI|EUCST)"), as.double)) %>%
  pivot_longer(-c(ab, disk_dose), names_to = "method") %>% 
  separate(method, into = c("guideline", "method"), sep = "_") %>% 
  mutate(method = ifelse(method %like% "D",
                         gsub("D", "DISK_", method, fixed = TRUE),
                         gsub("M", "MIC_", method, fixed = TRUE))) %>% 
  separate(method, into = c("method", "rsi"), sep = "_") %>% 
  # I is in the middle, so we only need R and S (saves data)
  filter(rsi %in% c("R", "S")) %>% 
  pivot_wider(names_from = rsi, values_from = value) %>%
  transmute(guideline = gsub("([0-9]+)$", " 20\\1", gsub("EUCST", "EUCAST", guideline)),
            method,
            site = NA_character_,
            mo = as.mo("UNKNOWN"),
            ab,
            ref_tbl = "Generic rules",
            disk_dose,
            breakpoint_S = S,
            breakpoint_R = R,
            uti = FALSE) %>% 
  filter(!(is.na(breakpoint_S) & is.na(breakpoint_R)), !is.na(mo), !is.na(ab))
rsi_generic

# create data set for AB-specific and MO-specific rules
rsi_specific <- DRGLST1 %>% 
  # only support guidelines for humans (for now)
  filter(HOST == "Human" & SITE_INF %unlike% "canine|feline",
         # only CLSI and EUCAST
         GUIDELINES %like% "(CLSI|EUCST)") %>% 
  # get microorganism names from another WHONET table
  mutate(ORG_CODE = tolower(ORG_CODE)) %>% 
  left_join(ORGLIST %>%
              transmute(ORG_CODE = tolower(ORG),
                        SCT_TEXT = case_when(is.na(SCT_TEXT) & is.na(ORGANISM) ~ ORG_CODE,
                                             is.na(SCT_TEXT) ~ ORGANISM,
                                             TRUE ~ SCT_TEXT)) %>% 
              # WHO for 'Generic'
              bind_rows(tibble(ORG_CODE = "gen", SCT_TEXT = "Unknown")) %>% 
              # WHO for 'Enterobacterales'
              bind_rows(tibble(ORG_CODE = "ebc", SCT_TEXT = "Enterobacterales"))
  ) %>% 
  # still some manual cleaning required
  filter(!SCT_TEXT %in% c("Anaerobic Actinomycetes")) %>% 
  transmute(guideline = gsub("([0-9]+)$", " 20\\1", gsub("EUCST", "EUCAST", GUIDELINES)),
            method = toupper(TESTMETHOD),
            site = SITE_INF,
            mo = as.mo(SCT_TEXT),
            ab = as.ab(WHON5_CODE),
            ref_tbl = REF_TABLE,
            disk_dose = POTENCY,
            breakpoint_S = as.double(ifelse(method == "DISK", DISK_S, MIC_S)),
            breakpoint_R = as.double(ifelse(method == "DISK", DISK_R, MIC_R)),
            uti = site %like% "(UTI|urinary|urine)") %>% 
  filter(!(is.na(breakpoint_S) & is.na(breakpoint_R)), !is.na(mo), !is.na(ab))
rsi_specific

rsi_translation <- rsi_generic %>% 
  bind_rows(rsi_specific) %>% 
  # add the taxonomic rank index, used for sorting (so subspecies match first, order matches last)
  mutate(rank_index = case_when(mo_rank(mo) %like% "(infra|sub)" ~ 1,
                                mo_rank(mo) == "species" ~ 2,
                                mo_rank(mo) == "genus" ~ 3,
                                mo_rank(mo) == "family" ~ 4,
                                mo_rank(mo) == "order" ~ 5,
                                TRUE ~ 6),
         .after = mo) %>% 
  arrange(desc(guideline), ab, mo, method) %>% 
  distinct(guideline, ab, mo, method, site, .keep_all = TRUE) %>% 
  as.data.frame(stringsAsFactors = FALSE)

# disks MUST be 6-50 mm, so correct where that is wrong:
rsi_translation[which(rsi_translation$method == "DISK" &
                        (is.na(rsi_translation$breakpoint_S) | rsi_translation$breakpoint_S > 50)), "breakpoint_S"] <- 50
rsi_translation[which(rsi_translation$method == "DISK" &
                        (is.na(rsi_translation$breakpoint_R) | rsi_translation$breakpoint_R < 6)), "breakpoint_R"] <- 6
m <- unique(as.double(as.mic(levels(as.mic(1)))))
rsi_translation[which(rsi_translation$method == "MIC" &
                        is.na(rsi_translation$breakpoint_S)), "breakpoint_S"] <- min(m)
rsi_translation[which(rsi_translation$method == "MIC" &
                        is.na(rsi_translation$breakpoint_R)), "breakpoint_R"] <- max(m)

# WHONET has no >1024 but instead uses 1025, 513, etc, so raise these one higher valid MIC factor level:
rsi_translation[which(rsi_translation$breakpoint_R == 129), "breakpoint_R"] <- m[which(m == 128) + 1]
rsi_translation[which(rsi_translation$breakpoint_R == 257), "breakpoint_R"] <- m[which(m == 256) + 1]
rsi_translation[which(rsi_translation$breakpoint_R == 513), "breakpoint_R"] <- m[which(m == 512) + 1]
rsi_translation[which(rsi_translation$breakpoint_R == 1025), "breakpoint_R"] <- m[which(m == 1024) + 1]

# Greek symbols and EM dash symbols are not allowed by CRAN, so replace them with ASCII:
rsi_translation$disk_dose <- gsub("μ", "u", rsi_translation$disk_dose, fixed = TRUE)
rsi_translation$disk_dose <- gsub("–", "-", rsi_translation$disk_dose, fixed = TRUE)

# save to package
usethis::use_data(rsi_translation, overwrite = TRUE)
rm(rsi_translation)
devtools::load_all(".")
