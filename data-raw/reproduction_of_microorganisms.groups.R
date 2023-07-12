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

# this data set is being used in the clinical_breakpoints data set, and thus by as.sir().
# it prevents the breakpoints table from being extremely long for species that are part of a species group.

library(dplyr)
library(readr)
library(tidyr)
devtools::load_all()

# Install the WHONET software on Windows (http://www.whonet.org/software.html),
# and copy the folder C:\WHONET\Resources to the data-raw/WHONET/ folder

# READ DATA ----

whonet_organisms <- read_tsv("data-raw/WHONET/Resources/Organisms.txt", na = c("", "NA", "-"), show_col_types = FALSE) %>%
  # remove old taxonomic names
  filter(TAXONOMIC_STATUS == "C") %>%
  mutate(ORGANISM_CODE = toupper(WHONET_ORG_CODE))

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
whonet_organisms %>%
  filter(!is.na(SPECIES_GROUP)) %>%
  group_by(SPECIES_GROUP) %>%
  summarise(complex = ORGANISM[ORGANISM %like% " (group|complex)"][1],
            organisms = paste0(n(), ": ", paste(sort(unique(ORGANISM)), collapse = ", "))) %>% 
  filter(!SPECIES_GROUP %in% AMR::microorganisms.codes$code)

# create the species group data set ----
microorganisms.groups <- whonet_organisms %>%
  # these will not be translated well
  filter(!ORGANISM %in% c("Trueperella pyogenes-like bacteria",
                          "Mycobacterium suricattae",
                          "Mycobacterium canetti")) %>% 
  filter(!is.na(SPECIES_GROUP), SPECIES_GROUP != ORGANISM_CODE) %>%
  transmute(mo_group = as.mo(SPECIES_GROUP),
            mo = ifelse(is.na(mo),
                        as.character(as.mo(ORGANISM, keep_synonyms = TRUE, minimum_matching_score = 0)),
                        mo)) %>% 
  # add our own CoNS and CoPS, WHONET does not strictly follow Becker et al (2014, 2019, 2020)
  filter(mo_group != as.mo("CoNS")) %>% 
  bind_rows(tibble(mo_group = as.mo("CoNS"), mo = MO_CONS)) %>% 
  filter(mo_group != as.mo("CoPS")) %>% 
  bind_rows(tibble(mo_group = as.mo("CoPS"), mo = MO_COPS)) %>% 
  # at least all our Lancefield-grouped streptococci must be in the beta-haemolytic group:
  bind_rows(tibble(mo_group = as.mo("Beta-haemolytic streptococcus"), 
                   mo = c(MO_LANCEFIELD,
                          microorganisms %>% filter(fullname %like% "^Streptococcus Group") %>% pull(mo)))) %>% 
  # and per Streptococcus group as well:
  # group A - S. pyogenes
  bind_rows(tibble(mo_group = as.mo("Streptococcus Group A"),
                   mo = microorganisms$mo[which(microorganisms$mo %like% "^B_STRPT_PYGN(_|$)")])) %>% 
  # group B - S. agalactiae
  bind_rows(tibble(mo_group = as.mo("Streptococcus Group B"),
                   mo = microorganisms$mo[which(microorganisms$mo %like% "^B_STRPT_AGLC(_|$)")])) %>% 
  # group C - all subspecies within S. dysgalactiae and S. equi (such as S. equi zooepidemicus)
  bind_rows(tibble(mo_group = as.mo("Streptococcus Group C"),
                   mo = microorganisms$mo[which(microorganisms$mo %like% "^B_STRPT_(DYSG|EQUI)(_|$)")])) %>% 
  # group F - S. anginosus, incl. S. anginosus anginosus and S. anginosus whileyi
  bind_rows(tibble(mo_group = as.mo("Streptococcus Group F"),
                   mo = microorganisms$mo[which(microorganisms$mo %like% "^B_STRPT_ANGN(_|$)")])) %>% 
  # group G - S. dysgalactiae and S. canis (though dysgalactiae is also group C and will be matched there)
  bind_rows(tibble(mo_group = as.mo("Streptococcus Group G"),
                   mo = microorganisms$mo[which(microorganisms$mo %like% "^B_STRPT_(DYSG|CANS)(_|$)")])) %>% 
  # group H - S. sanguinis
  bind_rows(tibble(mo_group = as.mo("Streptococcus Group H"),
                   mo = microorganisms$mo[which(microorganisms$mo %like% "^B_STRPT_SNGN(_|$)")])) %>% 
  # group K - S. salivarius, incl. S. salivarius salivariuss and S. salivarius thermophilus
  bind_rows(tibble(mo_group = as.mo("Streptococcus Group K"),
                   mo = microorganisms$mo[which(microorganisms$mo %like% "^B_STRPT_SLVR(_|$)")])) %>%
  # and for EUCAST: Strep group A, B, C, G
  bind_rows(tibble(mo_group = as.mo("Streptococcus Group A, B, C, G"),
                   mo = microorganisms$mo[which(microorganisms$mo %like% "^B_STRPT_(PYGN|AGLC|DYSG|EQUI|CANS|GRPA|GRPB|GRPC|GRPG)(_|$)")])) %>%
  # HACEK is:
  # - Haemophilus species
  # - Aggregatibacter species
  # - Cardiobacterium hominis
  # - Eikenella corrodens
  # - Kingella species
  # - and previously Actinobacillus actinomycetemcomitans
  # see https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3656887/
  filter(mo_group != as.mo("HACEK")) %>% 
  bind_rows(tibble(mo_group = as.mo("HACEK"), mo = microorganisms %>% filter(genus == "Haemophilus") %>% pull(mo))) %>% 
  bind_rows(tibble(mo_group = as.mo("HACEK"), mo = microorganisms %>% filter(genus == "Aggregatibacter") %>% pull(mo))) %>% 
  bind_rows(tibble(mo_group = as.mo("HACEK"), mo = as.mo("Cardiobacterium hominis", keep_synonyms = TRUE))) %>%
  bind_rows(tibble(mo_group = as.mo("HACEK"), mo = as.mo("Eikenella corrodens", keep_synonyms = TRUE))) %>%
  bind_rows(tibble(mo_group = as.mo("HACEK"), mo = microorganisms %>% filter(genus == "Kingella") %>% pull(mo))) %>%
  bind_rows(tibble(mo_group = as.mo("HACEK"), mo = as.mo("Actinobacillus actinomycetemcomitans", keep_synonyms = TRUE))) %>%
  # RGM are Rapidly-grwoing Mycobacteria, see https://pubmed.ncbi.nlm.nih.gov/28084211/
  filter(mo_group != "B_MYCBC_RGM") %>% 
  bind_rows(tibble(mo_group = as.mo("B_MYCBC_RGM"),
                   mo = paste("Mycobacterium", c( "abscessus abscessus", "abscessus bolletii", "abscessus massiliense", "agri", "aichiense", "algericum", "alvei", "anyangense", "arabiense", "aromaticivorans", "aubagnense", "aubagnense", "aurum", "austroafricanum", "bacteremicum", "boenickei", "bourgelatii", "brisbanense", "brumae", "canariasense", "celeriflavum", "chelonae", "chitae", "chlorophenolicum", "chubuense", "confluentis", "cosmeticum", "crocinum", "diernhoferi", "duvalii", "elephantis", "fallax", "flavescens", "fluoranthenivorans", "fortuitum", "franklinii", "frederiksbergense", "gadium", "gilvum", "goodii", "hassiacum", "hippocampi", "hodleri", "holsaticum", "houstonense", "immunogenum", "insubricum", "iranicum", "komossense", "litorale", "llatzerense", "madagascariense", "mageritense", "monacense", "moriokaense", "mucogenicum", "mucogenicum", "murale", "neoaurum", "neworleansense", "novocastrense", "obuense", "pallens", "parafortuitum", "peregrinum", "phlei", "phocaicum", "phocaicum", "porcinum", "poriferae", "psychrotolerans", "pyrenivorans", "rhodesiae", "rufum", "rutilum", "salmoniphilum", "sediminis", "senegalense", "septicum", "setense", "smegmatis", "sphagni", "thermoresistibile", "tokaiense", "vaccae", "vanbaalenii", "wolinskyi")) %>% as.mo(keep_synonyms = TRUE))) %>%
  # add full names
  mutate(mo_group_name = mo_name(mo_group, keep_synonyms = TRUE, language = NULL),
         mo_name = mo_name(mo, keep_synonyms = TRUE, language = NULL)) %>% 
  arrange(mo_group_name, mo_name) %>% 
  filter(mo_group != mo) %>% 
  distinct() %>% 
  dataset_UTF8_to_ASCII()
mo_uncertainties()

class(microorganisms.groups$mo_group) <- c("mo", "character")
class(microorganisms.groups$mo) <- c("mo", "character")
usethis::use_data(microorganisms.groups, internal = FALSE, overwrite = TRUE, compress = "xz", version = 2)
rm(microorganisms.groups)
devtools::load_all()
