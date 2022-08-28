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

# Go to https://lpsn.dsmz.de/downloads (register first) and download the latest CSV file.

file_location <- "data-raw/taxonomy.csv"

library(tidyverse)
library(AMR)

# these should still work after this update
test_fullname <- microorganisms$fullname
test_mo <- microorganisms$mo


# Helper functions --------------------------------------------------------

get_author_year <- function(ref) {
  # Only keep first author, e.g. transform 'Smith, Jones, 2011' to 'Smith et al., 2011'

  authors2 <- iconv(ref, from = "UTF-8", to = "ASCII//TRANSLIT")
  authors2 <- gsub(" ?\\(Approved Lists [0-9]+\\) ?", " () ", authors2)
  authors2 <- gsub(" [)(]+ $", "", authors2)
  # remove leading and trailing brackets
  authors2 <- trimws(gsub("^[(](.*)[)]$", "\\1", authors2))
  # only take part after brackets if there's a name
  authors2 <- ifelse(grepl(".*[)] [a-zA-Z]+.*", authors2),
    gsub(".*[)] (.*)", "\\1", authors2),
    authors2
  )
  # get year from last 4 digits
  lastyear <- as.integer(gsub(".*([0-9]{4})$", "\\1", authors2))
  # can never be later than now
  lastyear <- ifelse(lastyear > as.integer(format(Sys.Date(), "%Y")),
    NA,
    lastyear
  )
  # get authors without last year
  authors <- gsub("(.*)[0-9]{4}$", "\\1", authors2)
  # remove nonsense characters from names
  authors <- gsub("[^a-zA-Z,'& -]", "", authors)
  # remove trailing and leading spaces
  authors <- trimws(authors)
  # only keep first author and replace all others by 'et al'
  authors <- gsub("(,| and| et| &| ex| emend\\.?) .*", " et al.", authors)
  # et al. always with ending dot
  authors <- gsub(" et al\\.?", " et al.", authors)
  authors <- gsub(" ?,$", "", authors)
  # don't start with 'sensu' or 'ehrenb'
  authors <- gsub("^(sensu|Ehrenb.?) ", "", authors, ignore.case = TRUE)
  # no initials, only surname
  authors <- gsub("^([A-Z]+ )+", "", authors, ignore.case = FALSE)
  # combine author and year if year is available
  ref <- ifelse(!is.na(lastyear),
    paste0(authors, ", ", lastyear),
    authors
  )
  # fix beginning and ending
  ref <- gsub(", $", "", ref)
  ref <- gsub("^, ", "", ref)
  ref <- gsub("^(emend|et al.,?)", "", ref)
  ref <- trimws(ref)
  ref <- gsub("'", "", ref)

  # a lot start with a lowercase character - fix that
  ref[!grepl("^d[A-Z]", ref)] <- gsub("^([a-z])", "\\U\\1", ref[!grepl("^d[A-Z]", ref)], perl = TRUE)
  # specific one for the French that are named dOrbigny
  ref[grepl("^d[A-Z]", ref)] <- gsub("^d", "d'", ref[grepl("^d[A-Z]", ref)])
  ref <- gsub(" +", " ", ref)
  ref
}

df_remove_nonASCII <- function(df) {
  # Remove non-ASCII characters (these are not allowed by CRAN)
  df %>%
    mutate_if(is.character, iconv, from = "UTF-8", to = "ASCII//TRANSLIT") %>%
    # also remove invalid characters
    mutate_if(is.character, ~ gsub("[\"'`]+", "", .)) %>%
    AMR:::dataset_UTF8_to_ASCII()
}

abbreviate_mo <- function(x, minlength = 5, prefix = "", ...) {
  # keep a starting Latin ae
  suppressWarnings(
    gsub("^ae", "\u00E6\u00E6", x, ignore.case = TRUE) %>%
      abbreviate(
        minlength = minlength,
        use.classes = TRUE,
        method = "both.sides", ...
      ) %>%
      paste0(prefix, .) %>%
      toupper() %>%
      gsub("(\u00C6|\u00E6)+", "AE", .)
  )
}

# Read data ---------------------------------------------------------------

taxonomy <- read_csv(file_location)

# Create synonyms ---------------------------------------------------------

new_synonyms <- taxonomy %>%
  left_join(taxonomy,
    by = c("record_lnk" = "record_no"),
    suffix = c("", ".new")
  ) %>%
  filter(!is.na(record_lnk)) %>%
  mutate_all(~ ifelse(is.na(.), "", .)) %>%
  transmute(
    fullname = trimws(paste(genus_name, sp_epithet, subsp_epithet)),
    fullname_new = trimws(paste(genus_name.new, sp_epithet.new, subsp_epithet.new)),
    ref = get_author_year(authors),
    prevalence = 0
  ) %>%
  distinct(fullname, .keep_all = TRUE) %>%
  filter(fullname != fullname_new) %>%
  # this part joins this table to itself to correct for entries that had >1 renames,
  # such as:
  # Bacteroides tectum -> Bacteroides tectus -> Bacteroides pyogenes
  left_join(., .,
    by = c("fullname_new" = "fullname"),
    suffix = c("", ".2")
  ) %>%
  mutate(
    fullname_new = ifelse(!is.na(fullname_new.2), fullname_new.2, fullname_new),
    ref = ifelse(!is.na(ref.2), ref.2, ref)
  ) %>%
  select(-ends_with(".2"))

mo_became_synonym <- microorganisms %>%
  filter(fullname %in% new_synonyms$fullname)

updated_microorganisms <- taxonomy %>%
  filter(is.na(record_lnk)) %>%
  mutate_all(~ ifelse(is.na(.), "", .)) %>%
  transmute(
    mo = "",
    fullname = trimws(paste(genus_name, sp_epithet, subsp_epithet)),
    kingdom = "Bacteria",
    phylum = "",
    class = "",
    order = "",
    family = "",
    genus = trimws(genus_name),
    species = trimws(replace_na(sp_epithet, "")),
    subspecies = trimws(replace_na(subsp_epithet, "")),
    rank = case_when(
      subspecies == "" & species == "" ~ "genus",
      subspecies == "" ~ "species",
      TRUE ~ "subsp."
    ),
    ref = get_author_year(authors),
    species_id = as.character(record_no),
    source = "LPSN",
    prevalence = 0,
    snomed = NA
  )

new_microorganisms <- updated_microorganisms %>%
  filter(!fullname %in% microorganisms$fullname)

genera_with_mo_code <- updated_microorganisms %>%
  filter(genus %in% (microorganisms %>% filter(kingdom == "Bacteria", rank == "genus") %>% pull(genus))) %>%
  distinct(genus) %>%
  left_join(microorganisms %>% filter(kingdom == "Bacteria", rank == "genus") %>% select(mo, genus),
    by = "genus"
  )

genera_without_mo_code <- updated_microorganisms %>%
  filter(!genus %in% genera_with_mo_code$genus) %>%
  pull(genus) %>%
  unique()

genera_without_mo_code_abbr <- genera_without_mo_code %>% abbreviate_mo(5, prefix = "B_")
genera_without_mo_code_abbr[genera_without_mo_code_abbr %in% microorganisms$mo] <- abbreviate_mo(genera_without_mo_code[genera_without_mo_code_abbr %in% microorganisms$mo], 6, prefix = "B_")
genera_without_mo_code_abbr[genera_without_mo_code_abbr %in% microorganisms$mo] <- abbreviate_mo(genera_without_mo_code[genera_without_mo_code_abbr %in% microorganisms$mo], 7, prefix = "B_")
# all unique??
sum(genera_without_mo_code_abbr %in% microorganisms$mo) == 0

genus_abb <- tibble(
  genus = genera_without_mo_code,
  abbr = genera_without_mo_code_abbr
) %>%
  bind_rows(microorganisms %>%
    filter(kingdom == "Bacteria", rank == "genus", !genus %in% genera_without_mo_code) %>%
    transmute(genus, abbr = as.character(mo))) %>%
  arrange(genus)


# Update taxonomy ---------------------------------------------------------

# fill in the taxonomy of new genera
updated_taxonomy <- tibble(
  phylum = character(0),
  class = character(0),
  order = character(0),
  family = character(0),
  genus = character(0)
)
for (page in LETTERS) {
  message("Downloading page ", page, "... ", appendLF = FALSE)
  url <- paste0("https://lpsn.dsmz.de/genus?page=", page)

  x <- xml2::read_html(url) %>%
    rvest::html_node(".main-list") %>%
    # evety list element with a set <id> attribute
    rvest::html_nodes("li[id]")
  for (i in seq_len(length(x))) {
    txt <- x %>%
      magrittr::extract2(i) %>%
      rvest::html_text() %>%
      gsub("\\[[A-Za-z]+, no [a-z]+\\]", "NA", .) %>%
      gsub("Candidatus ", "", ., fixed = TRUE) %>%
      gsub("[ \t\r\n\"]+", "|", .) %>%
      gsub("\\|ShowHide.*", "", .) %>%
      gsub("[\\[\\]]", "", ., fixed = TRUE) %>%
      gsub("^\\|", "", .) %>%
      strsplit("|", fixed = TRUE) %>%
      unlist()
    txt[txt == "NA"] <- ""
    txt <- gsub("[^A-Za-z]+", "", txt)
    updated_taxonomy <- updated_taxonomy %>%
      bind_rows(tibble(
        phylum = txt[2],
        class = txt[3],
        order = txt[4],
        family = txt[5],
        genus = txt[6]
      ))
  }
  message(length(x), " entries (total ", nrow(updated_taxonomy), ")")
}

# Create new microorganisms -----------------------------------------------

new_microorganisms <- new_microorganisms %>%
  left_join(genus_abb, by = "genus") %>%
  group_by(genus) %>%
  mutate(species_abb = abbreviate_mo(species, 4)) %>%
  group_by(genus, species) %>%
  mutate(subspecies_abb = abbreviate_mo(subspecies, 4)) %>%
  ungroup() %>%
  mutate(
    mo = paste(abbr, species_abb, subspecies_abb, sep = "_"),
    mo = gsub("_+$", "", mo)
  ) %>%
  select(-matches("abb"))

# add taxonomy new microorganisms
MOs <- microorganisms %>%
  mutate(mo = as.character(mo)) %>%
  bind_rows(new_microorganisms) %>%
  arrange(fullname)

# unique MO codes
MOs$mo[which(duplicated(MOs$mo))] <- paste0(MOs$mo[which(duplicated(MOs$mo))], 1)
# all unique?
!any(duplicated(MOs$mo))

MOs <- MOs %>%
  # remove entries that are now a synonym
  filter(!fullname %in% new_synonyms$fullname) %>%
  # update the taxonomy
  left_join(updated_taxonomy, by = "genus", suffix = c("", ".new")) %>%
  mutate(
    phylum = ifelse(!is.na(phylum.new), phylum.new, phylum),
    class = ifelse(!is.na(class.new), class.new, class),
    order = ifelse(!is.na(order.new), order.new, order),
    family = ifelse(!is.na(family.new), family.new, family)
  ) %>%
  select(-ends_with(".new")) %>%
  # update prevalence based on taxonomy (Berends et al., 2021)
  mutate(prevalence = case_when(
    class == "Gammaproteobacteria" |
      genus %in% c("Enterococcus", "Staphylococcus", "Streptococcus")
    ~ 1,
    kingdom %in% c("Archaea", "Bacteria", "Chromista", "Fungi") &
      (phylum %in% c(
        "Proteobacteria",
        "Firmicutes",
        "Actinobacteria",
        "Sarcomastigophora"
      ) |
        genus %in% MO_PREVALENT_GENERA |
        rank %in% c("kingdom", "phylum", "class", "order", "family"))
    ~ 2,
    TRUE ~ 3
  ))

# add all mssing genera, families and orders
MOs <- MOs %>%
  bind_rows(
    MOs %>%
      arrange(genus, species) %>%
      distinct(genus, .keep_all = TRUE) %>%
      filter(rank == "species", source != "manually added") %>%
      mutate(
        mo = gsub("^([A-Z]_[A-Z]+)_.*", "\\1", mo),
        fullname = genus,
        species = "",
        subspecies = "",
        rank = "genus",
        species_id = "",
        snomed = NA,
        ref = NA_character_
      ),
    MOs %>%
      group_by(family) %>%
      filter(!any(rank == "family") & n() > 1) %>%
      ungroup() %>%
      arrange(family) %>%
      distinct(family, .keep_all = TRUE) %>%
      filter(!family %in% c("", NA), source != "manually added") %>%
      mutate(
        mo = paste0(
          substr(kingdom, 1, 1), "_[FAM]_",
          abbreviate(family,
            minlength = 8,
            use.classes = TRUE,
            method = "both.sides",
            strict = FALSE
          )
        ),
        mo = toupper(mo),
        fullname = family,
        genus = "",
        species = "",
        subspecies = "",
        rank = "family",
        species_id = "",
        snomed = NA,
        ref = NA_character_
      ),
    MOs %>%
      group_by(order) %>%
      filter(!any(rank == "order") & n() > 1) %>%
      ungroup() %>%
      arrange(order) %>%
      distinct(order, .keep_all = TRUE) %>%
      filter(!order %in% c("", NA), source != "manually added") %>%
      mutate(
        mo = paste0(
          substr(kingdom, 1, 1), "_[ORD]_",
          abbreviate(order,
            minlength = 8,
            use.classes = TRUE,
            method = "both.sides",
            strict = FALSE
          )
        ),
        mo = toupper(mo),
        fullname = order,
        family = "",
        genus = "",
        species = "",
        subspecies = "",
        rank = "order",
        species_id = "",
        snomed = NA,
        ref = NA_character_
      )
  ) %>%
  arrange(fullname)

# clean up
MOs <- MOs %>%
  df_remove_nonASCII()

# Add LPSN record IDs -----------------------------------------------------

records_ids <- taxonomy %>%
  mutate(across(1:3, function(x) {
    x[is.na(x)] <- ""
    x
  }),
  fullname = trimws(paste(genus_name, sp_epithet, subsp_epithet))
  ) %>%
  transmute(fullname, species_id = as.numeric(record_no)) %>%
  arrange(fullname, species_id) %>%
  distinct(fullname, .keep_all = TRUE)
message("Adding ", sum(records_ids$fullname %in% microorganisms$fullname), " LPSN record IDs")
MOs <- MOs %>%
  select(-species_id) %>%
  left_join(records_ids, by = "fullname") %>%
  relocate(species_id, .after = ref) %>%
  mutate(source = case_when(
    !is.na(species_id) ~ "LPSN",
    source %unlike% "manual" ~ "CoL",
    TRUE ~ source
  ))

# Merge synonyms ----------------------------------------------------------

# remove synonyms that are now valid names
MOs.old <- microorganisms.old %>%
  # add new synonyms
  bind_rows(new_synonyms) %>%
  filter(!fullname %in% MOs$fullname) %>%
  arrange(fullname) %>%
  distinct(fullname, fullname_new, .keep_all = TRUE) %>%
  # add prevalence to old taxonomic names
  select(-prevalence) %>%
  left_join(MOs %>% select(fullname, prevalence), by = c("fullname_new" = "fullname")) %>%
  # clean up
  df_remove_nonASCII()

message("microorganisms new:     ", sum(!MOs$fullname %in% c(microorganisms$fullname, MOs.old$fullname)))
message("microorganisms renamed: ", sum(!MOs.old$fullname %in% microorganisms.old$fullname))


# Save --------------------------------------------------------------------

# class <mo>
class(MOs$mo) <- c("mo", "character")

microorganisms <- MOs
microorganisms.old <- MOs.old

# --- Moraxella catarrhalis was named Branhamella catarrhalis (Catlin, 1970), but this is unaccepted in clinical microbiology
# we keep them both
microorganisms <- microorganisms %>%
  bind_rows(microorganisms %>%
    filter(fullname == "Branhamella catarrhalis") %>%
    mutate(
      mo = "B_MRXLL_CTRR",
      fullname = "Moraxella catarrhalis",
      genus = "Moraxella",
      ref = "Henriksen et al., 1968",
      species_id = "a374f6f0868e05f9c0f5077b60ee0a6c",
      snomed = as.list(24226003)
    )) %>%
  arrange(fullname) %>%
  df_remove_nonASCII()
microorganisms.old <- microorganisms.old %>%
  filter(fullname != "Moraxella catarrhalis")
# ---

# (this would be a great moment to run data-raw/snomed.R as well)

# on the server, do:
usethis::use_data(microorganisms, overwrite = TRUE, version = 2, compress = "xz")
usethis::use_data(microorganisms.old, overwrite = TRUE, version = 2, compress = "xz")
rm(microorganisms)
rm(microorganisms.old)

# DON'T FORGET TO UPDATE R/globals.R!

# load new data sets
devtools::load_all(".")

# reset previously changed mo codes
rsi_translation$mo <- as.mo(rsi_translation$mo, language = NULL)
usethis::use_data(rsi_translation, overwrite = TRUE, version = 2, compress = "xz")
rm(rsi_translation)

microorganisms.codes$mo <- as.mo(microorganisms.codes$mo, language = NULL)
# new NAs introduced?
any(is.na(microorganisms.codes$mo))
usethis::use_data(microorganisms.codes, overwrite = TRUE, version = 2, compress = "xz")
rm(microorganisms.codes)

example_isolates$mo <- as.mo(example_isolates$mo, language = NULL)
usethis::use_data(example_isolates, overwrite = TRUE, version = 2)
rm(example_isolates)

intrinsic_resistant$microorganism <- suppressMessages(mo_name(intrinsic_resistant$microorganism))
usethis::use_data(intrinsic_resistant, overwrite = TRUE, version = 2)
rm(intrinsic_resistant)

# load new data sets again
devtools::load_all(".")
source("data-raw/_pre_commit_hook.R")
devtools::load_all(".")


# Test updates ------------------------------------------------------------

# and check: these codes should not be missing (will otherwise throw a unit test error):
AMR::microorganisms.codes %>% filter(!mo %in% MOs$mo)
AMR::rsi_translation %>% filter(!mo %in% MOs$mo)
AMR:::microorganisms.translation %>% filter(!mo_new %in% MOs$mo)
AMR::example_isolates %>% filter(!mo %in% MOs$mo)

# Don't forget to add SNOMED codes! (data-raw/snomed.R)

# run the unit tests
Sys.setenv(NOT_CRAN = "true")
testthat::test_file("tests/testthat/test-data.R")
testthat::test_file("tests/testthat/test-mo.R")
testthat::test_file("tests/testthat/test-mo_property.R")
