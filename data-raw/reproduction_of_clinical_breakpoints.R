# ==================================================================== #
# TITLE:                                                               #
# AMR: An R Package for Working with Antimicrobial Resistance Data     #
#                                                                      #
# SOURCE CODE:                                                         #
# https://github.com/msberends/AMR                                     #
#                                                                      #
# PLEASE CITE THIS SOFTWARE AS:                                        #
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

# This script runs in 20-30 minutes and renews all guidelines of CLSI and EUCAST!
# Run it with source("data-raw/reproduction_of_clinical_breakpoints.R")

library(dplyr)
library(readr)
library(tidyr)
devtools::load_all()

# Install the WHONET 2022 software on Windows (http://www.whonet.org/software.html),
# and copy the folder C:\WHONET\Resources to the data-raw/WHONET/ folder
# (for ASIARS-Net update, also copy C:\WHONET\Codes to the data-raw/WHONET/ folder)

# READ DATA ----

whonet_organisms <- read_tsv("data-raw/WHONET/Resources/Organisms.txt", na = c("", "NA", "-"), show_col_types = FALSE) %>%
  # remove old taxonomic names
  filter(TAXONOMIC_STATUS == "C") %>%
  mutate(ORGANISM_CODE = toupper(WHONET_ORG_CODE))

whonet_breakpoints <- read_tsv("data-raw/WHONET/Resources/Breakpoints.txt", na = c("", "NA", "-"),
                               show_col_types = FALSE, guess_max = Inf) %>%
  filter(GUIDELINES %in% c("CLSI", "EUCAST"))

whonet_antibiotics <- read_tsv("data-raw/WHONET/Resources/Antibiotics.txt", na = c("", "NA", "-"), show_col_types = FALSE) %>%
  arrange(WHONET_ABX_CODE) %>%
  distinct(WHONET_ABX_CODE, .keep_all = TRUE)

# MICROORGANISMS WHONET CODES ----

whonet_organisms <- whonet_organisms %>%
  select(ORGANISM_CODE, ORGANISM, SPECIES_GROUP, GBIF_TAXON_ID) %>%
  mutate(
    # this one was called Issatchenkia orientalis, but it should be:
    ORGANISM = if_else(ORGANISM_CODE == "ckr", "Candida krusei", ORGANISM)
  ) %>% 
  # try to match on GBIF identifier
  left_join(microorganisms %>% distinct(mo, gbif, status) %>% filter(!is.na(gbif)), by = c("GBIF_TAXON_ID" = "gbif")) %>% 
  # remove duplicates
  arrange(ORGANISM_CODE, GBIF_TAXON_ID, status) %>%
  distinct(ORGANISM_CODE, .keep_all = TRUE) %>% 
  # add Enterobacterales, which is a subkingdom code in their data
  bind_rows(data.frame(ORGANISM_CODE = "ebc", ORGANISM = "Enterobacterales", mo = as.mo("Enterobacterales"))) %>% 
  arrange(ORGANISM)

# check non-existing species groups in the microorganisms table
complexes <- whonet_organisms %>%
  filter(ORGANISM %like% " (group|complex)$" & toupper(SPECIES_GROUP) == toupper(ORGANISM_CODE)) %>%
  mutate(mo = as.mo(ORGANISM, minimum_matching_score = 0.6, keep_synonyms = TRUE)) %>% 
  mutate(mo_new = paste0(mo, "-C"))
complexes[which(!complexes$mo_new %in% AMR::microorganisms$mo), ]


## Add new WHO codes to microorganisms.codes ----

matched <- whonet_organisms %>% filter(!is.na(mo))
unmatched <- whonet_organisms %>% filter(is.na(mo))

# generate the mo codes and add their names
message("Getting MO codes for WHONET input...")
unmatched <- unmatched %>% 
  mutate(mo = as.mo(gsub("(sero[a-z]*| complex| group| nontypable| non[-][a-zA-Z]+|var[.]| not .*|sp[.],.*|, .*variant.*|, .*toxin.*|, microaer.*| beta-haem[.])", "", ORGANISM),
                    minimum_matching_score = 0.6,
                    keep_synonyms = TRUE,
                    language = "en"),
         mo = case_when(ORGANISM %like% "Anaerobic" & ORGANISM %like% "negative" ~ as.mo("B_ANAER-NEG"),
                        ORGANISM %like% "Anaerobic" & ORGANISM %like% "positive" ~ as.mo("B_ANAER-POS"),
                        ORGANISM %like% "Anaerobic" ~ as.mo("B_ANAER"),
                        TRUE ~ mo),
         mo_name = mo_name(mo,
                           keep_synonyms = TRUE,
                           language = "en"))
# check if coercion at least resembles the first part (genus)
unmatched <- unmatched %>% 
  mutate(
    first_part = sapply(ORGANISM, function(x) strsplit(gsub("[^a-zA-Z _-]+", "", x), " ")[[1]][1], USE.NAMES = FALSE),
    keep = mo_name %like_case% first_part | ORGANISM %like% "Gram " | ORGANISM == "Other" | ORGANISM %like% "anaerobic") %>% 
  filter(keep == TRUE)

organisms <- matched %>% transmute(code = toupper(ORGANISM_CODE), group = SPECIES_GROUP, mo) %>% 
  bind_rows(unmatched %>% transmute(code = toupper(ORGANISM_CODE), group = SPECIES_GROUP, mo)) %>% 
  mutate(name = mo_name(mo, keep_synonyms = TRUE)) %>% 
  # remove the species groups themselves, we'll look the species within these groups later
  filter(is.na(group) | code != group) %>% 
  arrange(code)

# some subspecies exist, while their upper species do not, add them as the species level:
subspp <- organisms %>%
  filter(mo_species(mo, keep_synonyms = TRUE) == mo_subspecies(mo, keep_synonyms = TRUE) &
           mo_species(mo, keep_synonyms = TRUE) != "" &
           mo_genus(mo, keep_synonyms = TRUE) != "Salmonella") %>% 
  mutate(mo = as.mo(paste(mo_genus(mo, keep_synonyms = TRUE),
                          mo_species(mo, keep_synonyms = TRUE)),
                    keep_synonyms = TRUE),
         name = mo_name(mo, keep_synonyms = TRUE))
organisms <- organisms %>%
  filter(!code %in% subspp$code) %>%
  bind_rows(subspp) %>%
  arrange(code)

organism_groups <- organisms %>%
  filter(!is.na(group)) %>%
  arrange(group, name)
saveRDS(organism_groups, "data-raw/organism_groups.rds", version = 2)

#---
# AT THIS POINT, `organisms` is clean and all entries have an mo code
#---

# update microorganisms.codes with the latest WHONET codes
microorganisms.codes2 <- microorganisms.codes %>% 
  # remove all old WHONET codes, whether we (in the end) keep them or not
  filter(!toupper(code) %in% toupper(organisms$code) | toupper(code) %in% complexes$SPECIES_GROUP) %>% 
  # and add the new ones
  bind_rows(organisms %>% filter(code != group) %>% select(code, mo)) %>% 
  arrange(code)
# new codes:
microorganisms.codes2$code[which(!microorganisms.codes2$code %in% microorganisms.codes$code)]
mo_name(microorganisms.codes2$mo[which(!microorganisms.codes2$code %in% microorganisms.codes$code)], keep_synonyms = TRUE)
microorganisms.codes <- microorganisms.codes2

# Run this part to update ASIARS-Net:
# # start
# asiarsnet <- read_tsv("data-raw/WHONET/Codes/ASIARS_Net_Organisms_ForwardLookup.txt")
# asiarsnet <- asiarsnet %>%
#   mutate(WHONET_Code = toupper(WHONET_Code)) %>%
#   left_join(whonet_organisms %>% mutate(WHONET_Code = toupper(ORGANISM_CODE))) %>%
#   mutate(
#     mo1 = as.mo(ORGANISM_CODE),
#     mo2 = as.mo(ORGANISM)
#   ) %>%
#   mutate(mo = if_else(mo2 == "UNKNOWN" | is.na(mo2), mo1, mo2)) %>%
#   filter(!is.na(mo))
# insert1 <- asiarsnet %>% transmute(code = WHONET_Code, mo)
# insert2 <- asiarsnet %>% transmute(code = as.character(ASIARS_Net_Code), mo)
# # these will be updated
# bind_rows(insert1, insert2) %>%
#   rename(mo_new = mo) %>%
#   left_join(microorganisms.codes) %>%
#   filter(mo != mo_new)
# microorganisms.codes <- microorganisms.codes %>%
#   filter(!code %in% c(insert1$code, insert2$code)) %>%
#   bind_rows(insert1, insert2) %>%
#   arrange(code)
# # end

## Save to package ----
usethis::use_data(microorganisms.codes, overwrite = TRUE, compress = "xz", version = 2)
rm(microorganisms.codes)
devtools::load_all()


# BREAKPOINTS ----

# now that we have the right MO codes, get the breakpoints and convert them

whonet_breakpoints %>% 
  count(GUIDELINES, BREAKPOINT_TYPE) %>% 
  pivot_wider(names_from = BREAKPOINT_TYPE, values_from = n) %>% 
  janitor::adorn_totals(where = c("row", "col"))

breakpoints <- whonet_breakpoints %>%
  mutate(code = toupper(ORGANISM_CODE)) %>%
  left_join(bind_rows(microorganisms.codes %>% filter(!code %in% c("ALL", "GEN")),
                      # GEN (Generic) and ALL (All) are PK/PD codes
                      data.frame(code = c("ALL", "GEN"),
                                 mo = rep(as.mo("UNKNOWN"), 2))))
# these ones lack an MO name, they cannot be used:
unknown <- breakpoints %>%
  filter(is.na(mo) & !ORGANISM_CODE %in% organism_groups$group) %>%
  pull(code) %>%
  unique()
breakpoints %>% 
  filter(code %in% unknown) %>% 
  count(GUIDELINES, YEAR, ORGANISM_CODE, BREAKPOINT_TYPE, sort = TRUE)
breakpoints <- breakpoints %>% 
  filter(!is.na(mo) | ORGANISM_CODE %in% organism_groups$group)

# and these ones have unknown antibiotics according to WHONET itself:
breakpoints %>% 
  filter(!WHONET_ABX_CODE %in% whonet_antibiotics$WHONET_ABX_CODE) %>% 
  count(YEAR, GUIDELINES, WHONET_ABX_CODE) %>% 
  arrange(desc(YEAR))
breakpoints %>% 
  filter(!WHONET_ABX_CODE %in% whonet_antibiotics$WHONET_ABX_CODE) %>%
  pull(WHONET_ABX_CODE) %>%
  unique()
# they are at the moment all old codes that have right replacements in `antibiotics`, so we can use as.ab()


## Extend species groups ----
message("Extending breakpoints table on species groups...")
# get the species groups, they must be converted to rules per-species and added
breakpoints.bak <- breakpoints
spp_groups <- breakpoints %>% filter(ORGANISM_CODE_TYPE == "SPECIES_GROUP")
p <- progress_ticker(nrow(spp_groups))
for (i in seq_len(nrow(spp_groups))) {
  p$tick()
  mos <- organism_groups %>% filter(group == spp_groups[i, "ORGANISM_CODE", drop = TRUE]) %>% pull(mo)
  for (m in mos) {
    breakpoints <- breakpoints %>% 
      bind_rows(spp_groups %>% 
                  slice(i) %>% 
                  mutate(mo = m))
  }
}
close(p)

# extend all group A and B streptococci
breakpoints_new <- breakpoints_new %>% 
  bind_rows(breakpoints_new %>%
              filter(mo == as.mo("Streptococcus Group A")) %>% 
              mutate(mo = as.mo("Streptococcus pyogenes"))) %>% 
  bind_rows(breakpoints_new %>%
              filter(mo == as.mo("Streptococcus Group B")) %>% 
              mutate(mo = as.mo("Streptococcus agalactiae")))
# remove duplicates again for CoNS/CoPS/GBS and arrange
breakpoints_new <- breakpoints_new %>% 
  mutate(mo = as.mo(mo, keep_synonyms = TRUE)) %>% 
  distinct(guideline, ab, mo, method, site, .keep_all = TRUE) %>% 
  arrange(desc(guideline), ab, mo, method)

## Build new breakpoints table ----

breakpoints_new <- breakpoints %>%
  filter(!is.na(WHONET_ABX_CODE)) %>% 
  transmute(
    guideline = paste(GUIDELINES, YEAR),
    type = ifelse(BREAKPOINT_TYPE == "ECOFF", "ECOFF", tolower(BREAKPOINT_TYPE)),
    method = TEST_METHOD,
    site = SITE_OF_INFECTION,
    mo,
    rank_index = case_when(
      is.na(mo_rank(mo, keep_synonyms = TRUE)) ~ 6, # for UNKNOWN, B_GRAMN, B_ANAER, B_ANAER-NEG, etc.
      mo_rank(mo, keep_synonyms = TRUE) %like% "(infra|sub)" ~ 1,
      mo_rank(mo, keep_synonyms = TRUE) == "species" ~ 2,
      mo_rank(mo, keep_synonyms = TRUE) == "genus" ~ 3,
      mo_rank(mo, keep_synonyms = TRUE) == "family" ~ 4,
      mo_rank(mo, keep_synonyms = TRUE) == "order" ~ 5,
      TRUE ~ 6
    ),
    ab = as.ab(WHONET_ABX_CODE),
    ref_tbl = ifelse(type == "ECOFF" & is.na(REFERENCE_TABLE), "ECOFF", REFERENCE_TABLE),
    disk_dose = POTENCY,
    breakpoint_S = ifelse(type == "ECOFF" & is.na(S) & !is.na(ECV_ECOFF), ECV_ECOFF, S),
    breakpoint_R = ifelse(type == "ECOFF" & is.na(R) & !is.na(ECV_ECOFF), ECV_ECOFF, R),
    uti = ifelse(is.na(site), FALSE, gsub(".*(UTI|urinary|urine).*", "UTI", site) == "UTI")
  ) %>%
  # Greek symbols and EM dash symbols are not allowed by CRAN, so replace them with ASCII:
  mutate(disk_dose = disk_dose %>%
    gsub("μ", "u", ., fixed = TRUE) %>% # this is 'mu', \u03bc
    gsub("µ", "u", ., fixed = TRUE) %>% # this is 'micro', u00b5 (yes, they look the same)
    gsub("–", "-", ., fixed = TRUE)) %>%
  arrange(desc(guideline), ab, mo, method) %>%
  filter(!(is.na(breakpoint_S) & is.na(breakpoint_R)) & !is.na(mo) & !is.na(ab)) %>%
  distinct(guideline, ab, mo, method, site, breakpoint_S, .keep_all = TRUE)

# check the strange duplicates
breakpoints_new %>% 
  mutate(id = paste(guideline, type, ab, mo, method, site)) %>% 
  filter(id %in% .$id[which(duplicated(id))])
# remove duplicates
breakpoints_new <- breakpoints_new %>% 
  distinct(guideline, ab, mo, method, site, .keep_all = TRUE)

# fix reference table names
breakpoints_new %>% filter(guideline %like% "EUCAST", is.na(ref_tbl))
breakpoints_new <- breakpoints_new %>% 
  mutate(ref_tbl = case_when(is.na(ref_tbl) & guideline %like% "EUCAST 202" ~ lead(ref_tbl),
                             is.na(ref_tbl) ~ "Unknown",
                             TRUE ~ ref_tbl))

# clean disk zones
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
# EUCAST 2022 guideline: S <= 8 and R > 8
#           WHONET file: S <= 8 and R >= 16
breakpoints_new %>% filter(guideline == "EUCAST 2023", ab == "AMC", mo == "B_[ORD]_ENTRBCTR", method == "MIC")
# this will make an MIC of 12 I, which should be R, so:
breakpoints_new <- breakpoints_new %>%
  mutate(breakpoint_R = ifelse(guideline %like% "EUCAST" & method == "MIC" & log2(breakpoint_R) - log2(breakpoint_S) != 0,
    pmax(breakpoint_S, breakpoint_R / 2),
    breakpoint_R
  ))
# fix disks as well
breakpoints_new %>% filter(guideline == "EUCAST 2023", ab == "AMC", mo == "B_[ORD]_ENTRBCTR", method == "DISK")
breakpoints_new <- breakpoints_new %>%
  mutate(breakpoint_R = ifelse(guideline %like% "EUCAST" & method == "DISK" & breakpoint_S - breakpoint_R != 0,
    breakpoint_R + 1,
    breakpoint_R
  ))
# fix missing R breakpoint where there is an S breakpoint
breakpoints_new[which(is.na(breakpoints_new$breakpoint_R)), "breakpoint_R"] <- breakpoints_new[which(is.na(breakpoints_new$breakpoint_R)), "breakpoint_S"]

# check again
breakpoints_new %>% filter(guideline == "EUCAST 2023", ab == "AMC", mo == "B_[ORD]_ENTRBCTR", method == "MIC")
# compare with current version
clinical_breakpoints %>% filter(guideline == "EUCAST 2022", ab == "AMC", mo == "B_[ORD]_ENTRBCTR", method == "MIC")

# check dimensions
dim(breakpoints_new)
dim(clinical_breakpoints)

# SAVE TO PACKAGE ----

clinical_breakpoints <- breakpoints_new
usethis::use_data(clinical_breakpoints, overwrite = TRUE, compress = "xz", version = 2)
rm(clinical_breakpoints)
devtools::load_all(".")
