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

# This script runs in under a minute and renews all guidelines of CLSI and EUCAST!
# Run it with source("data-raw/reproduction_of_clinical_breakpoints.R")

library(dplyr)
library(readr)
library(tidyr)
devtools::load_all()

# Install the WHONET 2022 software on Windows (http://www.whonet.org/software.html),
# and copy the folder C:\WHONET\Resources to the data-raw/WHONET/ folder
# (for ASIARS-Net update, also copy C:\WHONET\Codes to the data-raw/WHONET/ folder)


# MICROORGANISMS WHONET CODES ----

whonet_organisms <- read_tsv("data-raw/WHONET/Resources/Organisms.txt", na = c("", "NA", "-"), show_col_types = FALSE) %>%
  # remove old taxonomic names
  filter(TAXONOMIC_STATUS == "C") %>%
  transmute(ORGANISM_CODE = tolower(WHONET_ORG_CODE), ORGANISM) %>%
  mutate(
    # what's wrong here? all these are only in the table on subspecies level (where species == subspecies), not on species level
    ORGANISM = if_else(ORGANISM_CODE == "sau", "Staphylococcus aureus", ORGANISM),
    ORGANISM = if_else(ORGANISM_CODE == "pam", "Pasteurella multocida", ORGANISM),
    ORGANISM = if_else(ORGANISM_CODE == "kpn", "Klebsiella pneumoniae", ORGANISM),
    ORGANISM = if_else(ORGANISM_CODE == "caj", "Campylobacter jejuni", ORGANISM),
    ORGANISM = if_else(ORGANISM_CODE == "mmo", "Morganella morganii", ORGANISM),
    ORGANISM = if_else(ORGANISM_CODE == "sap", "Staphylococcus saprophyticus", ORGANISM),
    ORGANISM = if_else(ORGANISM_CODE == "fne", "Fusobacterium necrophorum", ORGANISM),
    ORGANISM = if_else(ORGANISM_CODE == "fnu", "Fusobacterium nucleatum", ORGANISM),
    ORGANISM = if_else(ORGANISM_CODE == "sdy", "Streptococcus dysgalactiae", ORGANISM),
    ORGANISM = if_else(ORGANISM_CODE == "axy", "Achromobacter xylosoxidans", ORGANISM),
    # and this one was called Issatchenkia orientalis, but it should be:
    ORGANISM = if_else(ORGANISM_CODE == "ckr", "Candida krusei", ORGANISM)
  )

# add some general codes
whonet_organisms <- whonet_organisms %>%
  bind_rows(data.frame(
    ORGANISM_CODE = c("ebc", "cof"),
    ORGANISM = c("Enterobacterales", "Campylobacter")
  ))

whonet_organisms.bak <- whonet_organisms
# generate the mo codes and add their names
whonet_organisms <- whonet_organisms.bak %>% 
  mutate(mo = as.mo(gsub("(sero[a-z]*| complex| nontypable| non[-][a-zA-Z]+|var[.]| not .*|sp[.],.*|, .*variant.*|, .*toxin.*|, microaer.*| beta-haem[.])", "", ORGANISM),
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
new_mo_codes <- whonet_organisms %>% 
  mutate(
    first_part = sapply(ORGANISM, function(x) strsplit(gsub("[^a-zA-Z _-]+", "", x), " ")[[1]][1], USE.NAMES = FALSE),
    keep = mo_name %like_case% first_part | ORGANISM %like% "Gram " | ORGANISM == "Other" | ORGANISM %like% "anaerobic") %>% 
  filter(keep == TRUE) %>% 
  transmute(code = toupper(ORGANISM_CODE),
            mo = mo)
# update microorganisms.codes with the latest WHONET codes
microorganisms.codes2 <- microorganisms.codes %>% 
  # remove all old WHONET codes, whether we (in the end) keep them or not
  filter(!toupper(code) %in% toupper(new_mo_codes$code)) %>% 
  # and add the new ones
  bind_rows(new_mo_codes) %>% 
  arrange(code)
# new codes:
microorganisms.codes2$code[which(!microorganisms.codes2$code %in% microorganisms.codes$code)]
mo_name(microorganisms.codes2$mo[which(!microorganisms.codes2$code %in% microorganisms.codes$code)], keep_synonyms = TRUE)
microorganisms.codes <- microorganisms.codes2

# Run this part to update ASIARS-Net:
# start
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
# end

# save to package
usethis::use_data(microorganisms.codes, overwrite = TRUE, compress = "xz", version = 2)
rm(microorganisms.codes)
devtools::load_all()


# BREAKPOINTS ----

# now that we have the right MO codes, get the breakpoints and convert them
whonet_breakpoints <- read_tsv("data-raw/WHONET/Resources/Breakpoints.txt", na = c("", "NA", "-"), show_col_types = FALSE) %>%
  filter(BREAKPOINT_TYPE == "Human", GUIDELINES %in% c("CLSI", "EUCAST"))
whonet_antibiotics <- read_tsv("data-raw/WHONET/Resources/Antibiotics.txt", na = c("", "NA", "-"), show_col_types = FALSE) %>%
  arrange(WHONET_ABX_CODE) %>%
  distinct(WHONET_ABX_CODE, .keep_all = TRUE)

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
  filter(code %in% unknown)
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
# we cannot use them
# breakpoints <- breakpoints %>% 
#   filter(WHONET_ABX_CODE %in% whonet_antibiotics$WHONET_ABX_CODE)
# now check with our own antibiotics
breakpoints %>% 
  filter(!toupper(WHONET_ABX_CODE) %in% antibiotics$ab) %>% 
  pull(WHONET_ABX_CODE) %>% 
  unique()
# they are at the moment all old codes that have right replacements in `antibiotics`, so we can use as.ab()
  
breakpoints_new <- breakpoints %>%
  # only last available 10 years
  # filter(YEAR > max(YEAR) - 10) %>%
  transmute(
    guideline = paste(GUIDELINES, YEAR),
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
    ref_tbl = REFERENCE_TABLE,
    disk_dose = POTENCY,
    breakpoint_S = S,
    breakpoint_R = R,
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
  mutate(id = paste(guideline, ab, mo, method, site)) %>% 
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

# ECOFFs ----

# ECOFF = Epidemiological Cut-Off
whonet_ecoff <- read_tsv("data-raw/WHONET/Resources/Breakpoints.txt", na = c("", "NA", "-"), show_col_types = FALSE) %>%
  filter(BREAKPOINT_TYPE == "ECOFF", GUIDELINES %in% c("CLSI", "EUCAST"))

ecoff <- whonet_ecoff %>% 
  filter(!ORGANISM_CODE %in% c("clu", "BFX", "PFX", "kma", "cdh")) %>% 
  transmute(guideline = paste(GUIDELINES, YEAR),
            mo = as.mo(ORGANISM_CODE, keep_synonyms = TRUE),
            ab = as.ab(WHONET_ABX_CODE),
            method = TEST_METHOD,
            ecoff = as.double(ECV_ECOFF)) %>% 
  filter(!is.na(ecoff)) %>% 
  distinct()

# join to breakpoints
breakpoints_new <- breakpoints_new %>%
  bind_rows(breakpoints_new %>%
              right_join(ecoff, by = c("guideline", "mo", "ab", "method"))) %>% 
  mutate(ref_tbl = ifelse(is.na(ref_tbl), "ECOFF", ref_tbl)) %>% 
  distinct(guideline, ab, mo, method, site, .keep_all = TRUE) %>% 
  arrange(desc(guideline), ab, mo, method) %>% 
  mutate(rank_index = case_when(
    is.na(mo_rank(mo, keep_synonyms = TRUE)) ~ 6, # for UNKNOWN, B_GRAMN, B_ANAER, B_ANAER-NEG, etc.
    mo_rank(mo, keep_synonyms = TRUE) %like% "(infra|sub)" ~ 1,
    mo_rank(mo, keep_synonyms = TRUE) == "species" ~ 2,
    mo_rank(mo, keep_synonyms = TRUE) == "genus" ~ 3,
    mo_rank(mo, keep_synonyms = TRUE) == "family" ~ 4,
    mo_rank(mo, keep_synonyms = TRUE) == "order" ~ 5,
    TRUE ~ 6
  )) %>% 
  mutate(uti = ifelse(is.na(uti), FALSE, uti)) %>% 
  relocate(ecoff, .after = breakpoint_R)

breakpoints_new.bak <- mutate(uti = ifelse(is.na(uti), FALSE, uti), .after = ecoff)

# EXTEND CoNS/CoPS/GAS/GBS ----

# extend all coagulase-postive/-negative staphylococci
CoNS <- breakpoints_new %>% filter(mo == as.mo("CoNS"))
for (m in MO_CONS[mo_subspecies(MO_CONS, keep_synonyms = TRUE) == ""]) {
  breakpoints_new <- breakpoints_new %>% 
    bind_rows(CoNS %>% 
                mutate(mo = m))
}
CoPS <- breakpoints_new %>% filter(mo == as.mo("CoPS"))
for (m in MO_COPS[mo_subspecies(MO_COPS, keep_synonyms = TRUE) == ""]) {
  breakpoints_new <- breakpoints_new %>% 
    bind_rows(CoPS %>% 
                mutate(mo = m))
}
# do the same for group A and B streptococci
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


# Save to package ----

clinical_breakpoints <- breakpoints_new
usethis::use_data(clinical_breakpoints, overwrite = TRUE, compress = "xz", version = 2)
rm(clinical_breakpoints)
devtools::load_all(".")
