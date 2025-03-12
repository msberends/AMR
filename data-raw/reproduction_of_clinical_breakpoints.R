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

# This script runs in 20-30 minutes and renews all guidelines of CLSI and EUCAST!
# Run it with source("data-raw/reproduction_of_clinical_breakpoints.R")

library(dplyr)
library(readr)
library(tidyr)
devtools::load_all()

# Install the WHONET software on Windows (http://www.whonet.org/software.html),
# and copy the folder C:\WHONET\Resources to the data-raw/WHONET/ folder
# (for ASIARS-Net update, also copy C:\WHONET\Codes to the data-raw/WHONET/ folder)

# BE SURE TO RUN data-raw/reproduction_of_microorganisms.groups.R FIRST TO GET THE GROUPS!

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


## Add new WHO codes to microorganisms.codes ----

matched <- whonet_organisms %>% filter(!is.na(mo))
unmatched <- whonet_organisms %>% filter(is.na(mo))

# generate the mo codes and add their names
message("Getting MO codes for WHONET input...")
unmatched <- unmatched %>% 
  mutate(mo = as.mo(gsub("(sero[a-z]*| nontypable| non[-][a-zA-Z]+|var[.]| not .*|sp[.],.*|, .*variant.*|, .*toxin.*|, microaer.*| beta-haem[.])", "", ORGANISM),
                    minimum_matching_score = 0.55,
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
  arrange(keep)
unmatched %>% 
  View()
unmatched <- unmatched %>% 
  filter(keep == TRUE)

organisms <- matched %>% transmute(code = toupper(ORGANISM_CODE), group = SPECIES_GROUP, mo) %>% 
  bind_rows(unmatched %>% transmute(code = toupper(ORGANISM_CODE), group = SPECIES_GROUP, mo)) %>% 
  mutate(name = mo_name(mo, keep_synonyms = TRUE)) %>% 
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

# add the groups
organisms <- organisms %>% 
  bind_rows(tibble(code = organisms %>% filter(!is.na(group)) %>% pull(group) %>% unique(),
                   group = NA,
                   mo = organisms %>% filter(!is.na(group)) %>% pull(group) %>% unique() %>% as.mo(keep_synonyms = TRUE),
                   name = mo_name(mo, keep_synonyms = TRUE))) %>% 
  arrange(code, group) %>% 
  select(-group) %>% 
  distinct()

# 2023-07-08 SGM is also Strep gamma in WHONET, must only be Slowly-growing Mycobacterium
# 2024-06-14 still the case
organisms <- organisms %>% 
  filter(!(code == "SGM" & name %like% "Streptococcus"))
# this must be empty:
organisms$code[organisms$code %>% duplicated()]

saveRDS(organisms, "data-raw/organisms.rds", version = 2)

#---
# AT THIS POINT, `organisms` is clean and all entries have an mo code
#---

# update microorganisms.codes with the latest WHONET codes
microorganisms.codes2 <- microorganisms.codes %>% 
  # remove all old WHONET codes, whether we (in the end) keep them or not
  filter(!toupper(code) %in% toupper(organisms$code)) %>% 
  # and add the new ones
  bind_rows(organisms %>% select(code, mo)) %>% 
  arrange(code) %>% 
  distinct(code, .keep_all = TRUE)
# new codes:
microorganisms.codes2$code[which(!microorganisms.codes2$code %in% microorganisms.codes$code)]
mo_name(microorganisms.codes2$mo[which(!microorganisms.codes2$code %in% microorganisms.codes$code)], keep_synonyms = TRUE)
microorganisms.codes <- microorganisms.codes2

# Run this part to update ASIARS-Net:
# 2024-06-14: file not available anymore
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
class(microorganisms.codes$mo) <- c("mo", "character")
usethis::use_data(microorganisms.codes, overwrite = TRUE, compress = "xz", version = 2)
rm(microorganisms.codes)
devtools::load_all()


# BREAKPOINTS ----

# now that we have the correct MO codes, get the breakpoints and convert them

whonet_breakpoints %>% 
  count(GUIDELINES, BREAKPOINT_TYPE) %>% 
  pivot_wider(names_from = BREAKPOINT_TYPE, values_from = n) %>% 
  janitor::adorn_totals(where = c("row", "col"))
# compared to current
AMR::clinical_breakpoints %>%
  count(GUIDELINES = gsub("[^a-zA-Z]", "", guideline), type) %>%
  arrange(tolower(type)) %>%
  pivot_wider(names_from = type, values_from = n) %>% 
  as.data.frame() %>%
  janitor::adorn_totals(where = c("row", "col"))

breakpoints <- whonet_breakpoints %>%
  mutate(code = toupper(ORGANISM_CODE)) %>%
  left_join(bind_rows(microorganisms.codes %>% filter(!code %in% c("ALL", "GEN")),
                      # GEN (Generic) and ALL (All) are PK/PD codes
                      data.frame(code = c("ALL", "GEN"),
                                 mo = rep(as.mo("UNKNOWN"), 2))))
# these ones lack an MO name, they cannot be used:
unknown <- breakpoints %>%
  filter(is.na(mo)) %>%
  pull(code) %>%
  unique()
breakpoints %>% 
  filter(code %in% unknown) %>% 
  count(GUIDELINES, YEAR, ORGANISM_CODE, BREAKPOINT_TYPE, sort = TRUE)
# 2025-03-11: these codes are currently: clu, kma, fso, tyi. No clue (are not in MO list of WHONET), and they are only ECOFFs, so remove them:
breakpoints <- breakpoints %>% 
  filter(!is.na(mo))

# and these ones have unknown antibiotics according to WHONET itself:
breakpoints %>% 
  filter(!WHONET_ABX_CODE %in% whonet_antibiotics$WHONET_ABX_CODE) %>% 
  count(YEAR, GUIDELINES, WHONET_ABX_CODE) %>% 
  arrange(desc(YEAR))
breakpoints %>% 
  filter(!WHONET_ABX_CODE %in% whonet_antibiotics$WHONET_ABX_CODE) %>%
  pull(WHONET_ABX_CODE) %>%
  unique()
# they are at the moment all old codes ("CFC", "ROX", "FIX") that have the right replacements in `antimicrobials`, so we can use as.ab()


## Build new breakpoints table ----

breakpoints_new <- breakpoints %>%
  filter(!is.na(WHONET_ABX_CODE)) %>% 
  transmute(
    guideline = paste(GUIDELINES, YEAR),
    type = ifelse(BREAKPOINT_TYPE == "ECOFF", "ECOFF", tolower(BREAKPOINT_TYPE)),
    host = ifelse(BREAKPOINT_TYPE == "ECOFF", "ECOFF", tolower(HOST)),
    method = TEST_METHOD,
    site = SITE_OF_INFECTION,
    mo,
    rank_index = case_when(
      is.na(mo_rank(mo, keep_synonyms = TRUE)) ~ 6, # for UNKNOWN, B_GRAMN, B_ANAER, B_ANAER-NEG, etc.
      mo_rank(mo, keep_synonyms = TRUE) %like% "(infra|sub)" ~ 1,
      mo_rank(mo, keep_synonyms = TRUE) == "species" ~ 2,
      mo_rank(mo, keep_synonyms = TRUE) == "species group" ~ 2.5,
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
    uti = ifelse(is.na(site), FALSE, gsub(".*(UTI|urinary|urine).*", "UTI", site) == "UTI"),
    is_SDD = !is.na(SDD)
  ) %>%
  # Greek symbols and EM dash symbols are not allowed by CRAN, so replace them with ASCII:
  mutate(disk_dose = disk_dose %>%
    gsub("μ", "mc", ., fixed = TRUE) %>% # this is 'mu', \u03bc
    gsub("µ", "mc", ., fixed = TRUE) %>% # this is 'micro', \u00b5 (yes, they look the same)
    gsub("–", "-", ., fixed = TRUE)) %>%
  arrange(desc(guideline), mo, ab, type, method) %>%
  filter(!(is.na(breakpoint_S) & is.na(breakpoint_R)) & !is.na(mo) & !is.na(ab)) %>%
  distinct(guideline, type, host, ab, mo, method, site, breakpoint_S, .keep_all = TRUE)

# fix reference table names
breakpoints_new %>% filter(guideline %like% "EUCAST", is.na(ref_tbl)) %>% View()
breakpoints_new <- breakpoints_new %>% 
  mutate(ref_tbl = case_when(is.na(ref_tbl) & guideline %like% "EUCAST 202" ~ lead(ref_tbl),
                             is.na(ref_tbl) ~ "Unknown",
                             TRUE ~ ref_tbl))

# clean disk zones
breakpoints_new[which(breakpoints_new$method == "DISK"), "breakpoint_S"] <- as.double(as.disk(breakpoints_new[which(breakpoints_new$method == "DISK"), "breakpoint_S", drop = TRUE]))
breakpoints_new[which(breakpoints_new$method == "DISK"), "breakpoint_R"] <- as.double(as.disk(breakpoints_new[which(breakpoints_new$method == "DISK"), "breakpoint_R", drop = TRUE]))

# regarding animal breakpoints, CLSI has adults and foals for horses, but only for amikacin - only keep adult horses
breakpoints_new %>% 
  filter(host %like% "foal") %>%
  count(guideline, host)
breakpoints_new <- breakpoints_new %>% 
  filter(host %unlike% "foal") %>% 
  mutate(host = ifelse(host %like% "horse", "horse", host))

# FIXES FOR WHONET ERRORS ----
m <- unique(as.double(as.mic(levels(as.mic(1)))))

# WHONET has no >1024 but instead uses 1025, 513, etc, so as.mic() cannot be used to clean.
# instead, clean based on MIC factor levels
breakpoints_new[which(breakpoints_new$breakpoint_R == 129), "breakpoint_R"] <- m[which(m == 128) + 1]
breakpoints_new[which(breakpoints_new$breakpoint_R == 257), "breakpoint_R"] <- m[which(m == 256) + 1]
breakpoints_new[which(breakpoints_new$breakpoint_R == 513), "breakpoint_R"] <- m[which(m == 512) + 1]
breakpoints_new[which(breakpoints_new$breakpoint_R == 1025), "breakpoint_R"] <- m[which(m == 1024) + 1]

# a lot of R breakpoints are missing, though none of the S breakpoints are missing:
anyNA(breakpoints_new$breakpoint_S)

breakpoints_new %>% 
  filter(is.na(breakpoint_R)) %>% 
  count(guideline, host) |>
  pivot_wider(names_from = host,
              values_from = n,
              values_fill = list(n = 0)) |>
  View()

breakpoints_new[which(breakpoints_new$method == "MIC" &
                        is.na(breakpoints_new$breakpoint_R)), "breakpoint_R"] <- max(m)
# raise these one higher valid MIC factor level:


# fix streptococci in WHONET table of EUCAST: Strep A, B, C and G must only include these groups and not all streptococci:
breakpoints_new$mo[breakpoints_new$mo == "B_STRPT" & breakpoints_new$ref_tbl %like% "^strep.* a.* b.*c.*g"] <- as.mo("B_STRPT_ABCG")
# Haemophilus same error (must only be H. influenzae)
breakpoints_new$mo[breakpoints_new$mo == "B_HMPHL" & breakpoints_new$ref_tbl %like% "^h.* influenzae"] <- as.mo("B_HMPHL_INFL")
# EUCAST says that for H. parainfluenzae the H. influenza rules can be used, so add them
breakpoints_new <- breakpoints_new %>% 
  bind_rows(
    breakpoints_new %>%
      filter(guideline %like% "EUCAST", mo == "B_HMPHL_INFL") %>% 
      mutate(mo = as.mo("B_HMPHL_PRNF"))
  ) %>% 
  arrange(desc(guideline), mo, ab, type, host, method)
# Achromobacter denitrificans is in WHONET included in their A. xylosoxidans table, must be removed
breakpoints_new <- breakpoints_new %>% filter(mo != as.mo("Achromobacter denitrificans"))
# WHONET contains gentamicin breakpoints for viridans streptocci, which are intrinsic R - they meant genta-high, which is ALSO in their table, so we just remove gentamicin in viridans streptococci
breakpoints_new <- breakpoints_new %>% filter(!(mo == as.mo("Streptococcus viridans") & ab == "GEN"))
# Nitrofurantoin in Staph (EUCAST) only applies to S. saprophyticus, while WHONET has the DISK correct but the MIC on genus level
breakpoints_new$mo[breakpoints_new$mo == "B_STPHY" & breakpoints_new$ab == "NIT" & breakpoints_new$guideline %like% "EUCAST"] <- as.mo("B_STPHY_SPRP")
# WHONET sets the 2023 breakpoints for SAM to MIC of 16/32 for Enterobacterales, should be MIC 8/32 like AMC (see issue #123 on github.com/msberends/AMR)
# UPDATE 2024-02-22: fixed now

# determine rank again now that some changes were made on taxonomic level (genus -> species)
breakpoints_new <- breakpoints_new %>% 
  mutate(rank_index = case_when(
    is.na(mo_rank(mo, keep_synonyms = TRUE)) ~ 6, # for UNKNOWN, B_GRAMN, B_ANAER, B_ANAER-NEG, etc.
    mo_rank(mo, keep_synonyms = TRUE) %like% "(infra|sub)" ~ 1,
    mo_rank(mo, keep_synonyms = TRUE) == "species" ~ 2,
    mo_rank(mo, keep_synonyms = TRUE) == "species group" ~ 2.5,
    mo_rank(mo, keep_synonyms = TRUE) == "genus" ~ 3,
    mo_rank(mo, keep_synonyms = TRUE) == "family" ~ 4,
    mo_rank(mo, keep_synonyms = TRUE) == "order" ~ 5,
    TRUE ~ 6
  ))

# WHONET adds one log2 level to the R breakpoint for their software, e.g. in AMC in Enterobacterales:
# EUCAST 2023 guideline: S <= 8 and R > 8
#           WHONET file: S <= 8 and R >= 16
breakpoints_new %>% filter(guideline == "EUCAST 2023", ab == "AMC", mo == "B_[ORD]_ENTRBCTR", method == "MIC")
# but this will make an MIC of 12 I, which should be R according to EUCAST, so:
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
# fill missing R breakpoint where there is an S breakpoint
breakpoints_new[which(is.na(breakpoints_new$breakpoint_R)), "breakpoint_R"] <- breakpoints_new[which(is.na(breakpoints_new$breakpoint_R)), "breakpoint_S"]


# check the strange duplicates
breakpoints_new %>% 
  mutate(id = paste(guideline, type, host, method, site, mo, ab, uti)) %>% 
  filter(id %in% .$id[which(duplicated(id))]) %>% 
  arrange(desc(guideline))
# 2024-06-19 mostly ECOFFs, but there's no explanation in the whonet_breakpoints file, we have to remove duplicates
# remove duplicates
breakpoints_new <- breakpoints_new %>% 
  distinct(guideline, type, host, method, site, mo, ab, uti, .keep_all = TRUE)


# CHECKS AND SAVE TO PACKAGE ----

# check again
breakpoints_new %>% filter(guideline == "EUCAST 2024", ab == "AMC", mo == "B_[ORD]_ENTRBCTR", method == "MIC")
# compare with current version
clinical_breakpoints %>% filter(guideline == "EUCAST 2023", ab == "AMC", mo == "B_[ORD]_ENTRBCTR", method == "MIC")

# must have "human" and "ECOFF"
breakpoints_new %>% filter(mo == "B_STRPT_PNMN", ab == "AMP", guideline == "EUCAST 2020", method == "MIC")

# check dimensions
dim(breakpoints_new)
dim(clinical_breakpoints)

clinical_breakpoints <- breakpoints_new
clinical_breakpoints <- clinical_breakpoints %>% dataset_UTF8_to_ASCII()
usethis::use_data(clinical_breakpoints, overwrite = TRUE, compress = "xz", version = 2)
rm(clinical_breakpoints)
devtools::load_all(".")
