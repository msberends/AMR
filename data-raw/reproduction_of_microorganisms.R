# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Analysis                              #
#                                                                      #
# SOURCE                                                               #
# https://github.com/msberends/AMR                                     #
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
# Visit our website for more info: https://msberends.github.io/AMR.    #
# ==================================================================== #

# Reproduction of the `microorganisms` data set

# Data retrieved from the Catalogue of Life (CoL) through the Encyclopaedia of Life:
# https://opendata.eol.org/dataset/catalogue-of-life/
# Data retrieved from the Global Biodiversity Information Facility (GBIF): 
# https://doi.org/10.15468/rffz4x
#
# And from the Leibniz Institute: German Collection of Microorganisms and Cell Cultures (DSMZ)
# (register first at https://bacdive.dsmz.de/api/pnu/registration/register/ and use API as done below)

library(dplyr)
library(AMR)
# also needed: data.table, httr, jsonlite, cleaner, stringr

# unzip and extract taxa.txt (both around 1.5 GB, 3.7-3.9M rows) from Col and GBIF, then:
data_col_raw <- data.table::fread("data-raw/taxon.tab", quote = "")
data_gbif <- data.table::fread("data-raw/taxa.txt", quote = "")

# merge the two
data_col <- data_gbif %>% 
  rename(referenceID = identifier) %>% 
  bind_rows(data_col_raw) %>% 
  distinct(scientificName, kingdom, genus, specificEpithet, infraspecificEpithet, .keep_all = TRUE)
rm(data_col_raw)
rm(data_gbif)


# read the data from the DSMZ API (around 19000 rows)
dsmz_username <- ""
dsmz_password <- ""
GET_df <- function(url) {
  result <- httr::GET(url, httr::authenticate(dsmz_username, dsmz_password))
  httr::stop_for_status(result)
  result %>%
    httr::content(type = "text", encoding = "UTF-8") %>%
    jsonlite::fromJSON(flatten = TRUE)
}
dsmz_first <- GET_df("https://bacdive.dsmz.de/api/pnu/species?page=1&format=json")
data_dsmz <- dsmz_first$results
# this next process will take appr. `dsmz_first$count / 100 * 5 / 60` minutes
for (i in 2:round((dsmz_first$count / 100) + 0.5)) {
  data_dsmz <<- rbind(data_dsmz,
                      GET_df(paste0("https://bacdive.dsmz.de/api/pnu/species/?page=", i, "&format=json"))$results)
  cat(i, "-", AMR:::percentage(i / round((dsmz_first$count / 100) + 0.5)), "\n")
}
rm(dsmz_first)

# the CoL data is over 3.7M rows:
data_col %>% cleaner::freq(kingdom)
#      Item             Count   Percent   Cum. Count   Cum. Percent
# ---  ----------  ----------  --------  -----------  -------------
# 1    Animalia     2,494,992    55.43%    2,494,992         55.43%
# 2    Plantae      1,379,674    30.65%    3,874,666         86.08%
# 3    Fungi          547,619    12.17%    4,422,285         98.24%
# 4    Chromista       51,475     1.14%    4,473,760         99.39%
# 5    Bacteria        14,442     0.32%    4,488,202         99.71%
# 6    Protozoa         8,750     0.19%    4,496,952         99.90%
# 7    Viruses          3,805     0.08%    4,500,757         99.99%
# 8    Archaea            609     0.01%    4,501,366        100.00%

# clean data_col
data_col.bak <- data_col
data_col_old <- data_col %>%
  # filter: has new accepted name
  filter(!is.na(acceptedNameUsageID)) %>% 
  as_tibble() %>%
  transmute(fullname = trimws(stringr::str_replace(scientificName, 
                                                   pattern = stringr::fixed(scientificNameAuthorship), 
                                                   replacement = "")),
            fullname_new = trimws(paste(ifelse(is.na(genus), "", genus), 
                                        ifelse(is.na(specificEpithet), "", specificEpithet), 
                                        ifelse(is.na(infraspecificEpithet), "", infraspecificEpithet))),
            ref = scientificNameAuthorship,
            prevalence = NA_integer_)
data_col <- data_col %>%
  # filter: has no new accepted name
  filter(is.na(acceptedNameUsageID)) %>% 
  as_tibble() %>%
  transmute(fullname = "",
            kingdom,
            phylum,
            class,
            order,
            family,
            genus,
            species = specificEpithet,
            subspecies = infraspecificEpithet,
            rank = taxonRank,
            ref = scientificNameAuthorship,
            species_id = referenceID,
            source = "CoL")

# clean data_dsmz
data_dsmz.bak <- data_dsmz
data_dsmz_old <- data_dsmz %>%
  # filter: correct name is not NULL
  filter(!sapply(correct_name, is.null)) %>% 
  as_tibble() %>%
  transmute(fullname = trimws(paste(ifelse(is.na(genus), "", genus), 
                                    ifelse(is.na(species_epithet), "", species_epithet), 
                                    ifelse(is.na(subspecies_epithet), "", subspecies_epithet))),
            fullname_new = sapply(correct_name, function(x) x[2L]),
            ref = authors,
            prevalence = NA_integer_)

data_dsmz <- data_dsmz %>%
  # filter: correct name is NULL
  filter(sapply(correct_name, is.null)) %>% 
  as_tibble() %>%
  transmute(fullname = "",
            kingdom = regio,
            phylum,
            class = classis,
            # order = "", # does not contain order, will add later based on CoL
            family = familia,
            genus = ifelse(is.na(genus), "", genus),
            species = ifelse(is.na(species_epithet), "", species_epithet),
            subspecies = ifelse(is.na(subspecies_epithet), "", subspecies_epithet),
            rank = ifelse(species == "", "genus", "species"),
            ref = authors,
            species_id = as.character(pnu_no),
            source = "DSMZ")

# DSMZ only contains genus/(sub)species, try to find taxonomic properties based on genus and data_col
ref_taxonomy <- data_col %>%
  filter(family %in% data_dsmz$family & family != "") %>% 
  arrange(kingdom) %>% 
  distinct(family, .keep_all = TRUE) %>% 
  select(family, order)

data_dsmz <- data_dsmz %>%
  left_join(ref_taxonomy, by = "family") # NAs will later become "(unknown ...)"

# combine everything
data_total <- data_col %>%
  bind_rows(data_dsmz)

rm(data_col)
rm(data_dsmz)
rm(ref_taxonomy)
rm(data_col.bak)
rm(data_dsmz.bak)

mo_found_in_NL <- c("Absidia", "Acremonium", "Actinotignum", "Aedes", "Alternaria", "Anaerosalibacter", "Ancylostoma", 
                    "Angiostrongylus", "Anisakis", "Anopheles", "Apophysomyces", "Arachnia", "Ascaris", "Aspergillus", 
                    "Aureobacterium", "Aureobasidium", "Bacteroides", "Balantidum", "Basidiobolus", "Beauveria", 
                    "Bilophilia", "Blastocystis", "Branhamella", "Brochontrix", "Brugia", "Calymmatobacterium", "Candida", "Capillaria",
                    "Capnocytophaga", "Catabacter", "Cdc", "Chaetomium", "Chilomastix", "Chryseobacterium",
                    "Chryseomonas", "Chrysonilia", "Cladophialophora", "Cladosporium", "Clonorchis", "Conidiobolus",
                    "Contracaecum", "Cordylobia", "Cryptococcus", "Curvularia", "Demodex", "Dermatobia", "Dicrocoelium", 
                    "Dioctophyma", "Diphyllobothrium", "Dipylidium", "Dirofilaria", "Dracunculus", "Echinococcus", 
                    "Echinostoma", "Elisabethkingia", "Enterobius", "Enteromonas", "Euascomycetes", "Exophiala",
                    "Exserohilum", "Fasciola", "Fasciolopsis", "Flavobacterium", "Fonsecaea", "Fusarium", "Fusobacterium",
                    "Giardia", "Gnathostoma", "Hendersonula", "Heterophyes", "Hymenolepis", "Hypomyces", 
                    "Hysterothylacium", "Kloeckera", "Koserella", "Larva", "Lecythophora", "Leishmania", "Lelliottia",
                    "Leptomyxida", "Leptosphaeria", "Leptotrichia", "Loa", "Lucilia", "Lumbricus", "Malassezia", 
                    "Malbranchea", "Mansonella", "Mesocestoides", "Metagonimus", "Metarrhizium", "Molonomonas", 
                    "Mortierella", "Mucor", "Multiceps", "Mycocentrospora", "Mycoplasma", "Nanophetus", "Nattrassia",
                    "Necator", "Nectria", "Novospingobium", "Ochroconis", "Oesophagostomum", "Oidiodendron", "Onchocerca",
                    "Opisthorchis", "Opistorchis", "Paragonimus", "Paramyxovirus", "Pediculus", "Phlebotomus",
                    "Phocanema", "Phoma", "Phthirus", "Piedraia", "Pithomyces", "Pityrosporum", "Prevotella", 
                    "Pseudallescheria", "Pseudoterranova", "Pulex", "Retortamonas", "Rhizomucor", "Rhizopus",
                    "Rhodotorula", "Salinococcus", "Sanguibacteroides", "Sarcophagidae", "Sarcoptes", "Schistosoma", 
                    "Scolecobasidium", "Scopulariopsis", "Scytalidium", "Spirometra", "Sporobolomyces", "Stachybotrys",
                    "Stenotrophomononas", "Stomatococcus", "Strongyloides", "Syncephalastraceae", "Syngamus", "Taenia",
                    "Ternidens", "Torulopsis", "Toxocara", "Toxoplasma", "Treponema", "Trichinella", "Trichobilharzia", "Trichoderma",
                    "Trichomonas", "Trichophyton", "Trichosporon", "Trichostrongylus", "Trichuris", "Tritirachium",
                    "Trombicula", "Trypanosoma", "Tunga", "Ureaplasma", "Wuchereria")

MOs <- data_total %>%
  filter(
    (
      # we only want all MICROorganisms and no viruses
      !kingdom %in% c("Animalia", "Plantae", "Viruses")
      # and not all fungi: Aspergillus, Candida, Trichphyton and Pneumocystis are the most important,
      # so only keep these orders from the fungi:
      & !(kingdom == "Fungi"
          & !order %in% c("Eurotiales", "Microascales", "Mucorales", "Saccharomycetales", "Schizosaccharomycetales", "Tremellales", "Onygenales", "Pneumocystales"))
    )
    # or the genus has to be one of the genera we found in our hospitals last decades (Northern Netherlands, 2002-2018)
    | genus %in% mo_found_in_NL
  ) %>%
  # really no Plantae (e.g. Dracunculus exist both as worm and as plant)
  filter(kingdom != "Plantae") %>% 
  filter(!rank %in% c("kingdom", "phylum", "class", "order", "family", "genus"))

# include all ranks other than species for the included species
MOs <- MOs %>% bind_rows(data_total %>% 
                           filter((kingdom %in% MOs$kingdom & rank == "kingdom")
                                  | (phylum %in% MOs$phylum & rank == "phylum")
                                  | (class %in% MOs$class & rank == "class")
                                  | (order %in% MOs$order & rank == "order")
                                  | (family %in% MOs$family & rank == "family")
                                  | (genus %in% MOs$genus & rank == "genus")))

get_author_year <- function(ref) {
  # Only keep first author, e.g. transform 'Smith, Jones, 2011' to 'Smith et al., 2011'
  
  authors2 <- iconv(ref, from = "UTF-8", to = "ASCII//TRANSLIT")
  # remove leading and trailing brackets
  authors2 <- gsub("^[(](.*)[)]$", "\\1", authors2)
  # only take part after brackets if there's a name
  authors2 <- ifelse(grepl(".*[)] [a-zA-Z]+.*", authors2),
                     gsub(".*[)] (.*)", "\\1", authors2),
                     authors2)
  # get year from last 4 digits
  lastyear = as.integer(gsub(".*([0-9]{4})$", "\\1", authors2))
  # can never be later than now
  lastyear = ifelse(lastyear > as.integer(format(Sys.Date(), "%Y")),
                    NA,
                    lastyear)
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
                authors)
  # fix beginning and ending
  ref <- gsub(", $", "", ref)
  ref <- gsub("^, ", "", ref)
  ref <- gsub("^(emend|et al.,?)", "", ref)
  ref <- trimws(ref)
  
  # a lot start with a lowercase character - fix that
  ref[!grepl("^d[A-Z]", ref)] <- gsub("^([a-z])", "\\U\\1", ref[!grepl("^d[A-Z]", ref)], perl = TRUE)
  # specific one for the French that are named dOrbigny 
  ref[grepl("^d[A-Z]", ref)] <- gsub("^d", "d'", ref[grepl("^d[A-Z]", ref)])
  ref <- gsub(" +", " ", ref)
  ref
}

MOs <- MOs %>% mutate(ref = get_author_year(ref))

# Remove non-ASCII characters (these are not allowed by CRAN)
MOs <- MOs %>%
  lapply(iconv, from = "UTF-8", to = "ASCII//TRANSLIT") %>% 
  as_tibble(stringsAsFactors = FALSE) %>% 
  # remove invalid characters
  mutate_all(~gsub("[\"'`]+", "", .))

# set new fullnames
MOs <- MOs %>%
  mutate(fullname = trimws(case_when(rank == "family" ~ family,
                                     rank == "order" ~ order,
                                     rank == "class" ~ class,
                                     rank == "phylum" ~ phylum,
                                     rank == "kingdom" ~ kingdom,
                                     TRUE ~ paste(genus, species, subspecies))),
         fullname = gsub(" (var|f|subsp)[.]", "", fullname)) %>% 
  # remove text if it contains 'Not assigned', etc.
  mutate_all(function(x) ifelse(x %like% "(not assigned|homonym|mistake)", NA, x)) %>% 
  # clean taxonomy
  mutate(kingdom = ifelse(is.na(kingdom) | trimws(kingdom) == "", "(unknown kingdom)", trimws(kingdom)),
         phylum = ifelse(is.na(phylum) | trimws(phylum) == "", "(unknown phylum)", trimws(phylum)),
         class = ifelse(is.na(class) | trimws(class) == "", "(unknown class)", trimws(class)),
         order = ifelse(is.na(order) | trimws(order) == "", "(unknown order)", trimws(order)),
         family = ifelse(is.na(family) | trimws(family) == "", "(unknown family)", trimws(family)))

# Split old taxonomic names
MOs.old <- data_col_old %>% 
  filter(!gsub(" (var|f|subsp)[.]", "", fullname_new) %in% data_dsmz_old$fullname) %>% 
  bind_rows(data_dsmz_old) %>%
  mutate(fullname_new = gsub(" (var|f|subsp)[.]", "", fullname_new),
         fullname = gsub(" (var|f|subsp)[.]", "", fullname)) %>% 
  # for cases like Chlamydia pneumoniae -> Chlamydophila pneumoniae -> Chlamydia pneumoniae:
  filter(!fullname %in% fullname_new &
           fullname_new %in% MOs$fullname &
           !is.na(fullname) & 
           fullname != fullname_new) %>%
  distinct(fullname, .keep_all = TRUE) %>%
  arrange(fullname) %>% 
  mutate(ref = get_author_year(ref))

MOs <- MOs %>%
  # remove entries that are old and in MOs.old
  filter(!fullname %in% MOs.old$fullname) %>% 
  # mark up
  transmute(fullname,
            kingdom,
            phylum,
            class,
            order,
            family,
            genus,
            species,
            subspecies,
            rank,
            ref,
            species_id = gsub("[^a-zA-Z0-9].*", "", species_id),
            source) %>%
  # prefer known taxonomy over unknown taxonomy, then DSMZ over CoL (= desc)
  arrange(desc(kingdom, genus, species, source)) %>%
  distinct(kingdom, fullname, .keep_all = TRUE)

# remove all genera that have no species - they are irrelevant for microbiology and almost all from the kingdom of Animalia
to_remove <- MOs %>%
  filter(!kingdom %in% c("Bacteria", "Protozoa")) %>% 
  group_by(kingdom, genus) %>%
  count() %>%
  filter(n == 1) %>%
  ungroup() %>%
  mutate(kingdom_genus = paste(kingdom, genus)) %>% 
  pull(kingdom_genus)
MOs <- MOs %>% filter(!(paste(kingdom, genus) %in% to_remove))
rm(to_remove)

# add all mssing genera, families and orders
MOs <- MOs %>% 
  bind_rows(MOs %>% 
              arrange(genus, species) %>%
              distinct(genus, .keep_all = TRUE) %>%
              filter(rank == "species") %>%
              mutate(fullname = genus, 
                     species = "", 
                     rank = "genus", 
                     species_id = "",
                     ref = NA_character_)) %>% 
  bind_rows(MOs %>% 
              arrange(family, genus) %>%
              distinct(family, .keep_all = TRUE) %>%
              filter(rank == "genus") %>%
              mutate(fullname = family, 
                     genus = "",
                     rank = "family", 
                     species_id = "",
                     ref = NA_character_)) %>% 
  bind_rows(MOs %>% 
              arrange(order, family) %>%
              distinct(family, .keep_all = TRUE) %>%
              filter(rank == "family") %>%
              mutate(fullname = order, 
                     family = "",
                     rank = "order", 
                     species_id = "",
                     ref = NA_character_))

# remove the empty ones
MOs <- MOs %>%
  mutate(fullname = gsub(",.*", "", fullname)) %>% 
  distinct(kingdom, fullname, .keep_all = TRUE) %>% 
  filter(fullname != "")

# what characters are in the fullnames?
table(sort(unlist(strsplit(x = paste(MOs$fullname, collapse = ""), split = ""))))
MOs %>% filter(!fullname %like% "^[a-z ]+$") %>% arrange(fullname) %>% View()

table(MOs$kingdom, MOs$rank)
table(AMR::microorganisms$kingdom, AMR::microorganisms$rank)

# set prevalence per species
MOs <- MOs %>%
  mutate(prevalence = case_when(
    class == "Gammaproteobacteria"
    | genus %in% c("Enterococcus", "Staphylococcus", "Streptococcus")
    ~ 1,
    kingdom %in% c("Archaea", "Bacteria", "Chromista", "Fungi")
    & (phylum %in% c("Proteobacteria",
                     "Firmicutes",
                     "Actinobacteria",
                     "Sarcomastigophora")
       | genus %in% mo_found_in_NL
       | rank %in% c("kingdom", "phylum", "class", "order", "family"))
    ~ 2,
    TRUE ~ 3
  ))

# Add abbreviations so we can easily know which ones are which ones.
# These will become valid and unique microbial IDs for the AMR package.
MOs <- MOs %>%
  arrange(prevalence, genus, species, subspecies) %>% 
  group_by(kingdom) %>%
  mutate(abbr_other = case_when(
    rank == "family" ~ paste0("[FAM]_",
                              abbreviate(family,
                                         minlength = 8,
                                         use.classes = TRUE,
                                         method = "both.sides",
                                         strict = FALSE)),
    rank == "order" ~ paste0("[ORD]_",
                             abbreviate(order,
                                        minlength = 8,
                                        use.classes = TRUE,
                                        method = "both.sides",
                                        strict = FALSE)),
    rank == "class" ~ paste0("[CLS]_",
                             abbreviate(class,
                                        minlength = 8,
                                        use.classes = TRUE,
                                        method = "both.sides",
                                        strict = FALSE)),
    rank == "phylum" ~ paste0("[PHL]_",
                              abbreviate(phylum,
                                         minlength = 8,
                                         use.classes = TRUE,
                                         method = "both.sides",
                                         strict = FALSE)),
    rank == "kingdom" ~ paste0("[KNG]_", kingdom),
    TRUE ~ NA_character_
  )) %>%
  # abbreviations may be same for genera between kingdoms,
  # because each abbreviation starts with the the first character(s) of the kingdom
  mutate(abbr_genus = abbreviate(gsub("^ae", "\u00E6\u00E6", genus, ignore.case = TRUE), # keep a starting Latin ae
                                 minlength = 5,
                                 use.classes = TRUE,
                                 method = "both.sides")) %>%
  ungroup() %>%
  group_by(genus) %>%
  # species abbreviations may be the same between genera
  # because the genus abbreviation is part of the abbreviation
  mutate(abbr_species = abbreviate(gsub("^ae", "\u00E6\u00E6", species),
                                   minlength = 4,
                                   use.classes = TRUE,
                                   method = "both.sides")) %>%
  ungroup() %>%
  group_by(genus, species) %>%
  mutate(abbr_subspecies = abbreviate(gsub("^ae", "\u00E6\u00E6", subspecies),
                                      minlength = 4,
                                      use.classes = TRUE,
                                      method = "both.sides")) %>%
  ungroup() %>%
  # remove trailing underscores
  mutate(mo = gsub("_+$", "",
                   toupper(paste(ifelse(kingdom %in% c("Animalia", "Plantae"),
                                        substr(kingdom, 1, 2),
                                        substr(kingdom, 1, 1)),
                                 ifelse(is.na(abbr_other),
                                        paste(abbr_genus,
                                              abbr_species,
                                              abbr_subspecies,
                                              sep = "_"),
                                        abbr_other),
                                 sep = "_"))),
         mo = gsub("(\u00C6|\u00E6)+", "AE", mo)) %>%
  mutate(mo = ifelse(duplicated(.$mo),
                     # these one or two must be unique too
                     paste0(mo, "1"),
                     mo),
         fullname = ifelse(fullname == "",
                           trimws(paste(genus, species, subspecies)),
                           fullname)) %>%
  # put `mo` in front, followed by the rest
  select(mo, everything(), -abbr_other, -abbr_genus, -abbr_species, -abbr_subspecies)

# add non-taxonomic entries
MOs <- MOs %>%
  bind_rows(
    # Unknowns
    data.frame(mo = "UNKNOWN",
               fullname = "(unknown name)",
               kingdom = "(unknown kingdom)",
               phylum = "(unknown phylum)",
               class = "(unknown class)",
               order = "(unknown order)",
               family = "(unknown family)",
               genus = "(unknown genus)",
               species = "(unknown species)",
               subspecies = "(unknown subspecies)",
               rank = "(unknown rank)",
               ref = NA_character_,
               species_id = "",
               source = "manually added",
               prevalence = 1,
               stringsAsFactors = FALSE),
    data.frame(mo = "B_GRAMN",
               fullname = "(unknown Gram-negatives)",
               kingdom = "Bacteria",
               phylum = "(unknown phylum)",
               class = "(unknown class)",
               order = "(unknown order)",
               family = "(unknown family)",
               genus = "(unknown Gram-negatives)",
               species = "(unknown species)",
               subspecies = "(unknown subspecies)",
               rank = "species",
               ref = NA_character_,
               species_id = "",
               source = "manually added",
               prevalence = 1,
               stringsAsFactors = FALSE),
    data.frame(mo = "B_GRAMP",
               fullname = "(unknown Gram-positives)",
               kingdom = "Bacteria",
               phylum = "(unknown phylum)",
               class = "(unknown class)",
               order = "(unknown order)",
               family = "(unknown family)",
               genus = "(unknown Gram-positives)",
               species = "(unknown species)",
               subspecies = "(unknown subspecies)",
               rank = "species",
               ref = NA_character_,
               species_id = "",
               source = "manually added",
               prevalence = 1,
               stringsAsFactors = FALSE),
    data.frame(mo = "F_YEAST",
               fullname = "(unknown yeast)",
               kingdom = "Fungi",
               phylum = "(unknown phylum)",
               class = "(unknown class)",
               order = "(unknown order)",
               family = "(unknown family)",
               genus = "(unknown genus)",
               species = "(unknown species)",
               subspecies = "(unknown subspecies)",
               rank = "species",
               ref = NA_character_,
               species_id = "",
               source = "manually added",
               prevalence = 2,
               stringsAsFactors = FALSE),
    data.frame(mo = "F_FUNGUS",
               fullname = "(unknown fungus)",
               kingdom = "Fungi",
               phylum = "(unknown phylum)",
               class = "(unknown class)",
               order = "(unknown order)",
               family = "(unknown family)",
               genus = "(unknown genus)",
               species = "(unknown species)",
               subspecies = "(unknown subspecies)",
               rank = "species",
               ref = NA_character_,
               species_id = "",
               source = "manually added",
               prevalence = 2,
               stringsAsFactors = FALSE),
    # CoNS
    MOs %>%
      filter(genus == "Staphylococcus", species == "epidermidis") %>% .[1,] %>%
      mutate(mo = gsub("(.*)_(.*)_.*", "\\1_\\2_CONS", mo),
             species = "coagulase-negative",
             fullname = "Coagulase-negative Staphylococcus (CoNS)",
             ref = NA_character_,
             species_id = "",
             source = "manually added"),
    # CoPS
    MOs %>%
      filter(genus == "Staphylococcus", species == "epidermidis") %>% .[1,] %>%
      mutate(mo = gsub("(.*)_(.*)_.*", "\\1_\\2_COPS", mo),
             species = "coagulase-positive",
             fullname = "Coagulase-positive Staphylococcus (CoPS)",
             ref = NA_character_,
             species_id = "",
             source = "manually added"),
    # Streptococci groups A, B, C, F, H, K
    MOs %>%
      filter(genus == "Streptococcus", species == "pyogenes") %>% .[1,] %>%
      # we can keep all other details, since S. pyogenes is the only member of group A
      mutate(mo = gsub("(.*)_(.*)_.*", "\\1_\\2_GRPA", mo),
             species = "group A" ,
             fullname = "Streptococcus group A",
             source = "manually added"),
    MOs %>%
      filter(genus == "Streptococcus", species == "agalactiae") %>% .[1,] %>%
      # we can keep all other details, since S. agalactiae is the only member of group B
      mutate(mo = gsub("(.*)_(.*)_.*", "\\1_\\2_GRPB", mo),
             species = "group B" ,
             fullname = "Streptococcus group B",
             source = "manually added"),
    MOs %>%
      filter(genus == "Streptococcus", species == "dysgalactiae") %>% .[1,] %>%
      mutate(mo = gsub("(.*)_(.*)_.*", "\\1_\\2_GRPC", mo),
             species = "group C" ,
             fullname = "Streptococcus group C",
             ref = NA_character_,
             species_id = "",
             source = "manually added"),
    MOs %>%
      filter(genus == "Streptococcus", species == "agalactiae") %>% .[1,] %>%
      mutate(mo = gsub("(.*)_(.*)_.*", "\\1_\\2_GRPD", mo),
             species = "group D" ,
             fullname = "Streptococcus group D",
             ref = NA_character_,
             species_id = "",
             source = "manually added"),
    MOs %>%
      filter(genus == "Streptococcus", species == "agalactiae") %>% .[1,] %>%
      mutate(mo = gsub("(.*)_(.*)_.*", "\\1_\\2_GRPF", mo),
             species = "group F" ,
             fullname = "Streptococcus group F",
             ref = NA_character_,
             species_id = "",
             source = "manually added"),
    MOs %>%
      filter(genus == "Streptococcus", species == "agalactiae") %>% .[1,] %>%
      mutate(mo = gsub("(.*)_(.*)_.*", "\\1_\\2_GRPG", mo),
             species = "group G" ,
             fullname = "Streptococcus group G",
             ref = NA_character_,
             species_id = "",
             source = "manually added"),
    MOs %>%
      filter(genus == "Streptococcus", species == "agalactiae") %>% .[1,] %>%
      mutate(mo = gsub("(.*)_(.*)_.*", "\\1_\\2_GRPH", mo),
             species = "group H" ,
             fullname = "Streptococcus group H",
             ref = NA_character_,
             species_id = "",
             source = "manually added"),
    MOs %>%
      filter(genus == "Streptococcus", species == "agalactiae") %>% .[1,] %>%
      mutate(mo = gsub("(.*)_(.*)_.*", "\\1_\\2_GRPK", mo),
             species = "group K" ,
             fullname = "Streptococcus group K",
             ref = NA_character_,
             species_id = "",
             source = "manually added"),
    # Beta haemolytic Streptococci
    MOs %>%
      filter(genus == "Streptococcus", species == "agalactiae") %>% .[1,] %>%
      mutate(mo = gsub("(.*)_(.*)_.*", "\\1_\\2_HAEM", mo),
             species = "beta-haemolytic" ,
             fullname = "Beta-haemolytic Streptococcus",
             ref = NA_character_,
             species_id = "",
             source = "manually added"),
    # Viridans Streptococci
    MOs %>%
      filter(genus == "Streptococcus", species == "agalactiae") %>% .[1,] %>%
      mutate(mo = gsub("(.*)_(.*)_.*", "\\1_\\2_VIRI", mo),
             species = "viridans" ,
             fullname = "Viridans Group Streptococcus (VGS)",
             ref = NA_character_,
             species_id = "",
             source = "manually added"),
    # Milleri Streptococci
    MOs %>%
      filter(genus == "Streptococcus", species == "agalactiae") %>% .[1,] %>%
      mutate(mo = gsub("(.*)_(.*)_.*", "\\1_\\2_MILL", mo),
             species = "milleri" ,
             fullname = "Milleri Group Streptococcus (MGS)",
             ref = NA_character_,
             species_id = "",
             source = "manually added"),
    # Candida krusei
    MOs %>%
      filter(genus == "Candida", species == "glabrata") %>% .[1,] %>%
      mutate(mo = gsub("(.*)_(.*)_.*", "\\1_\\2_KRUS", mo),
             species = "krusei" ,
             fullname = "Candida krusei",
             ref = NA_character_,
             species_id = "",
             source = "manually added"),
    # Blastocystis hominis does not exist (it means 'got a Blastocystis from humans', PMID 15634993)
    # but let's be nice to the clinical people in microbiology
    MOs %>%
      filter(fullname == "Blastocystis") %>%
      mutate(mo = paste0(mo, "_HMNS"),
             fullname = paste(fullname, "hominis"),
             species = "hominis",
             source = "manually added",
             ref = NA_character_,
             species_id = ""),
    # Trichomonas vaginalis is missing, same order as Dientamoeba
    MOs %>%
      filter(fullname == "Dientamoeba") %>%
      mutate(mo = gsub("(.*?)_.*", "\\1_THMNS", mo),
             fullname = "Trichomonas",
             family = "Trichomonadidae",
             genus = "Trichomonas",
             source = "manually added",
             ref = "Donne, 1836",
             species_id = ""),
    MOs %>%
      filter(fullname == "Dientamoeba fragilis") %>%
      mutate(mo = gsub("(.*?)_.*", "\\1_THMNS_VAG", mo),
             fullname = "Trichomonas vaginalis",
             family = "Trichomonadidae",
             genus = "Trichomonas",
             species = "vaginalis",
             source = "manually added",
             ref = "Donne, 1836",
             species_id = ""),
    MOs %>% # add family as such too
      filter(fullname == "Monocercomonadidae") %>%
      mutate(mo = gsub("(.*)_(.*)_.*", "\\1_\\2_TRCHMNDD", mo),
             fullname = "Trichomonadidae",
             family = "Trichomonadidae",
             rank = "family",
             genus = "",
             species = "",
             source = "manually added",
             ref = "",
             species_id = ""),
  )

# Incorporate new microbial order for Gammaproteobacteria - Adeolu et al. (2016), PMID 27620848
MOs[which(MOs$family == "Enterobacteriaceae"), "family"] <- ""
MOs[which(MOs$genus %in% c("Escherichia",
                           "Atlantibacter",
                           "Biostraticola",
                           "Buttiauxella",
                           "Cedecea",
                           "Citrobacter",
                           "Cronobacter",
                           "Enterobacillus",
                           "Enterobacter",
                           "Franconibacter",
                           "Gibbsiella",
                           "Izhakiella",
                           "Klebsiella",
                           "Kluyvera",
                           "Kosakonia",
                           "Leclercia",
                           "Lelliottia",
                           "Mangrovibacter",
                           "Pluralibacter",
                           "Pseudocitrobacter",
                           "Raoultella",
                           "Rosenbergiella",
                           "Saccharobacter",
                           "Salmonella",
                           "Shigella",
                           "Shimwellia",
                           "Siccibacter",
                           "Trabulsiella",
                           "Yokenella")), "family"] <- "Enterobacteriaceae"
MOs[which(MOs$genus %in% c("Erwinia",
                           "Buchnera",
                           "Pantoea",
                           "Phaseolibacter",
                           "Tatumella",
                           "Wigglesworthia")), "family"] <- "Erwiniaceae"
MOs[which(MOs$genus %in% c("Pectobacterium",
                           "Brenneria",
                           "Dickeya",
                           "Lonsdalea",
                           "Sodalis")), "family"] <- "Pectobacteriaceae"
MOs[which(MOs$genus %in% c("Yersinia",
                           "Chania",
                           "Ewingella",
                           "Rahnella",
                           "Rouxiella",
                           "Samsonia",
                           "Serratia")), "family"] <- "Yersiniaceae"
MOs[which(MOs$genus %in% c("Hafnia",
                           "Edwardsiella",
                           "Obesumbacterium")), "family"] <- "Hafniaceae"
MOs[which(MOs$genus %in% c("Morganella",
                           "Arsenophonus",
                           "Cosenzaea",
                           "Moellerella",
                           "Photorhabdus",
                           "Proteus",
                           "Providencia",
                           "Xenorhabdus")), "family"] <- "Morganellaceae"
MOs[which(MOs$genus %in% c("Budvicia",
                           "Leminorella",
                           "Pragia")), "family"] <- "Budviciaceae"
MOs[which(MOs$family %in% c("Enterobacteriaceae",
                            "Erwiniaceae",
                            "Pectobacteriaceae",
                            "Yersiniaceae",
                            "Hafniaceae",
                            "Morganellaceae",
                            "Budviciaceae")), "order"] <- "Enterobacterales"
new_families <- MOs %>%
  filter(order == "Enterobacterales") %>%
  pull(family) %>%
  unique()

MOs <- MOs %>% 
  filter(!(rank == "family" & fullname %in% new_families)) %>% 
  bind_rows(tibble(mo = paste0("B_[FAM]_",
                               toupper(abbreviate(new_families,
                                                  minlength = 8,
                                                  use.classes = TRUE,
                                                  method = "both.sides",
                                                  strict = FALSE))),
                   fullname = new_families,
                   kingdom = "Bacteria",
                   phylum = "Proteobacteria",
                   class = "Gammaproteobacteria",
                   order = "Enterobacterales",
                   family = new_families,
                   genus = "",
                   species = "",
                   subspecies = "",
                   rank = "family",
                   ref = "Adeolu et al., 2016",
                   species_id = NA_character_,
                   source = "manually added",
                   prevalence = 1))

MOs[which(MOs$order == "Enterobacteriales"), "order"] <- "Enterobacterales"
MOs[which(MOs$fullname == "Enterobacteriales"), "fullname"] <- "Enterobacterales"

# add prevalence to old taxonomic names
MOs.old <- MOs.old %>% 
  select(-prevalence) %>% 
  left_join(MOs %>% select(fullname, prevalence), by = c("fullname_new" = "fullname"))

# everything distinct?
sum(duplicated(MOs$mo))
sum(duplicated(MOs$fullname))
colnames(MOs)

# add the ones we would delete now, that have unexisting codes and names (also in the old names)
MOs <- MOs %>% 
  mutate(mo = as.character(mo)) %>% 
  bind_rows(
    AMR::microorganisms %>%
      mutate(mo = as.character(mo)) %>% 
      filter(genus %in% gen & !fullname %in% AMR::microorganisms$fullname & 
               !fullname %in% AMR::microorganisms.old$fullname &
               !mo  %in% microorganisms$mo) %>%
      select(all_of(colnames(AMR::microorganisms)))
  )

# here we welcome the new ones:
MOs %>% arrange(fullname) %>% filter(!fullname %in% AMR::microorganisms$fullname) %>% View()
MOs.old %>% arrange(fullname) %>% filter(!fullname %in% AMR::microorganisms.old$fullname) %>% View()
# and the ones we lost:
# AMR::microorganisms %>% filter(!fullname %in% MOs$fullname) %>% View() # based on fullname
AMR::microorganisms %>% filter(!fullname %in% c(MOs$fullname, MOs.old$fullname)) %>% View() # excluding renamed ones
# AMR::microorganisms %>% filter(!mo %in% MOs$mo) %>% View()             # based on mo
# AMR::microorganisms %>% filter(!mo %in% MOs$mo & !fullname %in% MOs$fullname) %>% View() 
# and these IDs have changed:
old_new <- MOs %>%
  mutate(kingdom_fullname = paste(kingdom, fullname)) %>% 
  filter(kingdom_fullname %in% (AMR::microorganisms %>% 
                                  mutate(kingdom_fullname = paste(kingdom, fullname)) %>%
                                  pull(kingdom_fullname))) %>%
  left_join(AMR::microorganisms %>% 
              mutate(kingdom_fullname = paste(kingdom, fullname)) %>%
              select(mo, kingdom_fullname), by = "kingdom_fullname", suffix = c("_new", "_old")) %>% 
  filter(mo_new != mo_old) %>% 
  select(mo_old, mo_new, everything())
View(old_new)

# set new MO codes as names to existing data sets
rsi_translation$mo <- mo_name(rsi_translation$mo, language = NULL)
microorganisms.codes$mo <- mo_name(microorganisms.codes$mo, language = NULL)
microorganisms.translation <- AMR:::microorganisms.translation %>%
  bind_rows(tibble(mo_old = AMR:::microorganisms.translation$mo_new, mo_new = mo_old)) %>%
  filter(!mo_old %in% MOs$mo) %>% 
  mutate(mo_new = mo_name(mo_new, language = NULL)) %>% 
  bind_rows(old_new %>% select(mo_old, mo_new)) %>% 
  distinct(mo_old, .keep_all = TRUE)

# arrange the data sets to save
MOs <- MOs %>% arrange(fullname)
MOs.old <- MOs.old %>% arrange(fullname)

# transform
MOs <- as.data.frame(MOs, stringsAsFactors = FALSE)
MOs.old <- as.data.frame(MOs.old, stringsAsFactors = FALSE)
microorganisms.codes <- as.data.frame(microorganisms.codes, stringsAsFactors = FALSE)
class(MOs$mo) <- c("mo", "character")

# SAVE
### for same server
microorganisms <- dataset_UTF8_to_ASCII(MOs)
microorganisms.old <- dataset_UTF8_to_ASCII(MOs.old)
### for other server
saveRDS(MOs, "microorganisms.rds")
saveRDS(MOs.old, "microorganisms.old.rds")
saveRDS(microorganisms.codes, "microorganisms.codes.rds")

# on the server, do:
usethis::use_data(microorganisms, overwrite = TRUE, version = 2)
usethis::use_data(microorganisms.old, overwrite = TRUE, version = 2)
rm(microorganisms)
rm(microorganisms.old)

# load new data sets
devtools::load_all(".")

# reset previously changed mo codes
rsi_translation$mo <- as.mo(rsi_translation$mo)
microorganisms.codes$mo <- as.mo(microorganisms.codes$mo)
class(microorganisms.codes$mo) <- c("mo", "character")
microorganisms.translation <- microorganisms.translation %>%
  # (to do: add last package version to column pkg_version)
  left_join(microorganisms.old[, c("fullname", "fullname_new")], # microorganisms.old is now new and loaded
            by = c("mo_new" = "fullname")) %>%
  mutate(name = ifelse(!is.na(fullname_new), fullname_new, mo_new)) %>% 
  left_join(microorganisms[, c("fullname", "mo")],               # as is microorganisms
            by = c("name" = "fullname")) %>% 
  select(mo_old, mo_new = mo) %>% 
  filter(!is.na(mo_old), !is.na(mo_new))
class(microorganisms.translation$mo_old) <- "character" # no class <mo> since those aren't valid MO codes
class(microorganisms.translation$mo_new) <- c("mo", "character")
# save those to the package
usethis::use_data(rsi_translation, overwrite = TRUE, version = 2)
usethis::use_data(microorganisms.codes, overwrite = TRUE, version = 2)
saveRDS(microorganisms.translation, file = "data-raw/microorganisms.translation.rds", version = 2)
# to save microorganisms.translation internally to the package
source("data-raw/internals.R")

# load new data sets again
devtools::load_all(".")

# and check: these codes should not be missing (will otherwise throw a unit test error):
AMR::microorganisms.codes %>% filter(!mo %in% MOs$mo)
AMR::rsi_translation %>% filter(!mo %in% MOs$mo)
AMR:::microorganisms.translation %>% filter(!mo_new %in% MOs$mo)

# update the example_isolates data set
example_isolates$mo <- as.mo(example_isolates$mo)
usethis::use_data(example_isolates, overwrite = TRUE)

# Don't forget to add SNOMED codes! (data-raw/snomed.R)

# run the unit tests
testthat::test_file("tests/testthat/test-data.R")
testthat::test_file("tests/testthat/test-mo.R")
testthat::test_file("tests/testthat/test-mo_property.R")

# edit 2020-05-28
# Not sure why it now says M. tuberculosis was renamed to M. africanum (B_MYCBC_AFRC), but that's not true
microorganisms <- microorganisms %>% 
  bind_rows(microorganisms %>% 
              filter(mo == "B_MYCBC_AFRC") %>%
              mutate(mo = "B_MYCBC_TBRC", snomed = list(c("113861009", "113858008")),
                     ref = "Lehmann et al., 2018",species_id = "778540",
                     source = "DSMZ", species = "tuberculosis",
                     fullname = "Mycobacterium tuberculosis")) %>% 
  arrange(fullname)
class(microorganisms$mo) <- c("mo", "character")
microorganisms.old <- microorganisms.old %>% filter(fullname != "Mycobacterium tuberculosis")

usethis::use_data(microorganisms, overwrite = TRUE, version = 2, compress = "xz")
usethis::use_data(microorganisms.old, overwrite = TRUE, version = 2)


# OLD CODE ----------------------------------------------------------------

# to keep all the old IDs:
# MOs <- MOs %>% filter(!mo %in% old_new$mo_new) %>%
#   rbind(microorganisms %>%
#           filter(mo %in% old_new$mo_old) %>%
#           select(mo, fullname) %>%
#           left_join(MOs %>%
#                       select(-mo), by = "fullname"))
# this is how to fix it
# microorganisms.codes <- AMR::microorganisms.codes %>% 
#   left_join(MOs %>%
#               mutate(kingdom_fullname = paste(kingdom, fullname)) %>% 
#               left_join(AMR::microorganisms %>%
#                           transmute(mo, kingdom_fullname = paste(kingdom, fullname)),
#                         by = "kingdom_fullname", suffix = c("_new", "_old")) %>%
#               select(mo_old, mo_new),
#             by = c("mo" = "mo_old")) %>% 
#   select(code, mo = mo_new) %>% 
#   filter(!is.na(mo))
# microorganisms.codes %>% filter(!mo %in% MOs$mo)
# # and for microorganisms.translation:
# microorganisms.translation <- AMR:::microorganisms.translation %>% 
#   select(mo = mo_new) %>% 
#   left_join(AMR::microorganisms %>%
#               transmute(mo, kingdom_fullname = paste(kingdom, fullname)),
#             by = "kingdom_fullname", suffix = c("_new", "_old")) %>%
#   select(mo_old, mo_new)
#   left_join(MOs %>%
#               mutate(kingdom_fullname = paste(kingdom, fullname)) %>% 
#               left_join(AMR::microorganisms %>%
#                           transmute(mo, kingdom_fullname = paste(kingdom, fullname)),
#                         by = "kingdom_fullname", suffix = c("_new", "_old")) %>%
#               select(mo_old, mo_new),
#             by = c("mo" = "mo_old")) %>% 
#   select(code, mo = mo_new) %>% 
#   filter(!is.na(mo))
# microorganisms.codes %>% filter(!mo %in% MOs$mo)
