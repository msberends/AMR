# ---------------------------------------------------------------------------------
# Reproduction of the `microorganisms` data set
# ---------------------------------------------------------------------------------
# Data retrieved from:
#
# [1] Catalogue of Life (CoL) through the Encyclopaedia of Life
#     https://opendata.eol.org/dataset/catalogue-of-life/
#      * Download the resource file with a name like "Catalogue of Life yyyy-mm-dd"
#      * Extract "taxon.tab"
#
# [2] Global Biodiversity Information Facility (GBIF)
#     https://doi.org/10.15468/39omei
#     * Extract "Taxon.tsv"
#
# [3] Deutsche Sammlung von Mikroorganismen und Zellkulturen (DSMZ)
#     https://www.dsmz.de/support/bacterial-nomenclature-up-to-date-downloads.html
#     * Download the latest "Complete List" as xlsx file (DSMZ_bactnames.xlsx)
# ---------------------------------------------------------------------------------

library(dplyr)
library(AMR)

data_col <- data.table::fread("Documents/taxon.tab")
data_gbif <- data.table::fread("Documents/Taxon.tsv")

# read the xlsx file from DSMZ (only around 2.5 MB):
data_dsmz <- readxl::read_xlsx("Downloads/DSMZ_bactnames.xlsx")

# the CoL data is over 3.7M rows:
data_col %>% freq(kingdom)
#      Item             Count   Percent   Cum. Count   Cum. Percent
# ---  ----------  ----------  --------  -----------  -------------
# 1    Animalia     2,225,627     59.1%    2,225,627          59.1%
# 2    Plantae      1,177,412     31.3%    3,403,039          90.4%
# 3    Fungi          290,145      7.7%    3,693,184          98.1%
# 4    Chromista       47,126      1.3%    3,740,310          99.3%
# 5    Bacteria        14,478      0.4%    3,754,788          99.7%
# 6    Protozoa         6,060      0.2%    3,760,848          99.9%
# 7    Viruses          3,827      0.1%    3,764,675         100.0%
# 8    Archaea            610      0.0%    3,765,285         100.0%

# the GBIF data is over 5.8M rows:
data_gbif %>% freq(kingdom)
#      Item                  Count   Percent   Cum. Count   Cum. Percent
# ---  ---------------  ----------  --------  -----------  -------------
# 1    Animalia          3,264,138     55.7%    3,264,138          55.7%
# 2    Plantae           1,814,962     31.0%    5,079,100          86.7%
# 3    Fungi               538,086      9.2%    5,617,186          95.9%
# 4    Chromista           181,374      3.1%    5,798,560          99.0%
# 5    Bacteria             24,048      0.4%    5,822,608          99.4%
# 6    Protozoa             15,138      0.3%    5,837,746          99.7%
# 7    incertae sedis        9,995      0.2%    5,847,741          99.8%
# 8    Viruses               9,630      0.2%    5,857,371         100.0%
# 9    Archaea                 771      0.0%    5,858,142         100.0%


# Clean up helper function ------------------------------------------------
clean_new <- function(new) {
  new %>%
    # only the ones that have no new ID to refer to a newer name
    filter(is.na(col_id_new)) %>%
    filter(
      (
        # we only want all MICROorganisms and no viruses
        !kingdom %in% c("Animalia", "Chromista", "Plantae", "Viruses")
        # and not all fungi: Aspergillus, Candida, Trichphyton and Pneumocystis are the most important,
        # so only keep these orders from the fungi:
        & !(kingdom == "Fungi"
            & !order %in% c("Eurotiales", "Saccharomycetales", "Schizosaccharomycetales", "Tremellales", "Onygenales", "Pneumocystales"))
      )
      # or the family has to contain a genus we found in our hospitals last decades (Northern Netherlands, 2002-2018)
      | genus %in% c("Absidia", "Acremonium", "Actinotignum", "Alternaria", "Anaerosalibacter", "Ancylostoma", "Anisakis", "Apophysomyces",
                     "Arachnia", "Ascaris", "Aureobacterium", "Aureobasidium", "Balantidum", "Bilophilia", "Branhamella", "Brochontrix",
                     "Brugia", "Calymmatobacterium", "Catabacter", "Chilomastix", "Chryseomonas", "Cladophialophora", "Cladosporium",
                     "Clonorchis", "Cordylobia", "Curvularia", "Demodex", "Dermatobia", "Diphyllobothrium", "Dracunculus", "Echinococcus",
                     "Enterobius", "Euascomycetes", "Exophiala", "Fasciola", "Fusarium", "Hendersonula", "Hymenolepis", "Kloeckera",
                     "Koserella", "Larva", "Leishmania", "Lelliottia", "Loa", "Lumbricus", "Malassezia", "Metagonimus", "Molonomonas",
                     "Mucor", "Nattrassia", "Necator", "Novospingobium", "Onchocerca", "Opistorchis", "Paragonimus", "Paramyxovirus",
                     "Pediculus", "Phoma", "Phthirus", "Pityrosporum", "Pseudallescheria", "Pulex", "Rhizomucor", "Rhizopus", "Rhodotorula",
                     "Salinococcus", "Sanguibacteroides", "Schistosoma", "Scopulariopsis", "Scytalidium", "Sporobolomyces", "Stomatococcus",
                     "Strongyloides", "Syncephalastraceae", "Taenia", "Torulopsis", "Trichinella", "Trichobilharzia", "Trichomonas",
                     "Trichosporon", "Trichuris", "Trypanosoma", "Wuchereria")) %>%
    mutate(
      authors2 = iconv(ref, from = "UTF-8", to = "ASCII//TRANSLIT"),
      # remove leading and trailing brackets
      authors2 = gsub("^[(](.*)[)]$", "\\1", authors2),
      # only take part after brackets if there's a name
      authors2 = ifelse(grepl(".*[)] [a-zA-Z]+.*", authors2),
                        gsub(".*[)] (.*)", "\\1", authors2),
                        authors2),
      # get year from last 4 digits
      lastyear = as.integer(gsub(".*([0-9]{4})$", "\\1", authors2)),
      # can never be later than now
      lastyear = ifelse(lastyear > as.integer(format(Sys.Date(), "%Y")),
                        NA,
                        lastyear),
      # get authors without last year
      authors = gsub("(.*)[0-9]{4}$", "\\1", authors2),
      # remove nonsense characters from names
      authors = gsub("[^a-zA-Z,'& -]", "", authors),
      # remove trailing and leading spaces
      authors = trimws(authors),
      # only keep first author and replace all others by 'et al'
      authors = gsub("(,| and| et| &| ex| emend\\.?) .*", " et al.", authors),
      # et al. always with ending dot
      authors = gsub(" et al\\.?", " et al.", authors),
      authors = gsub(" ?,$", "", authors),
      # don't start with 'sensu' or 'ehrenb'
      authors = gsub("^(sensu|Ehrenb.?) ", "", authors, ignore.case = TRUE),
      # no initials, only surname
      authors = gsub("^([A-Z]+ )+", "", authors, ignore.case = FALSE),
      # combine author and year if year is available
      ref = ifelse(!is.na(lastyear),
                   paste0(authors, ", ", lastyear),
                   authors),
      # fix beginning and ending
      ref = gsub(", $", "", ref),
      ref = gsub("^, ", "", ref)) %>%
    # remove text if it contains 'Not assigned' like phylum in viruses
    mutate_all(~gsub("Not assigned", "", .)) %>%
    # Remove non-ASCII characters (these are not allowed by CRAN)
    lapply(iconv, from = "UTF-8", to = "ASCII//TRANSLIT") %>%
    as_tibble(stringsAsFactors = FALSE) %>%
    mutate(fullname = trimws(case_when(rank == "family" ~ family,
                                       rank == "order" ~ order,
                                       rank == "class" ~ class,
                                       rank == "phylum" ~ phylum,
                                       rank == "kingdom" ~ kingdom,
                                       TRUE ~ paste(genus, species, subspecies))))
}
clean_old <- function(old, new) {
  old %>%
    # only the ones that exist in the new data set
    filter(col_id_new %in% new$col_id) %>%
    mutate(
      authors2 = iconv(ref, from = "UTF-8", to = "ASCII//TRANSLIT"),
      # remove leading and trailing brackets
      authors2 = gsub("^[(](.*)[)]$", "\\1", authors2),
      # only take part after brackets if there's a name
      authors2 = ifelse(grepl(".*[)] [a-zA-Z]+.*", authors2),
                        gsub(".*[)] (.*)", "\\1", authors2),
                        authors2),
      # get year from last 4 digits
      lastyear = as.integer(gsub(".*([0-9]{4})$", "\\1", authors2)),
      # can never be later than now
      lastyear = ifelse(lastyear > as.integer(format(Sys.Date(), "%Y")),
                        NA,
                        lastyear),
      # get authors without last year
      authors = gsub("(.*)[0-9]{4}$", "\\1", authors2),
      # remove nonsense characters from names
      authors = gsub("[^a-zA-Z,'& -]", "", authors),
      # remove trailing and leading spaces
      authors = trimws(authors),
      # only keep first author and replace all others by 'et al'
      authors = gsub("(,| and| et| &| ex| emend\\.?) .*", " et al.", authors),
      # et al. always with ending dot
      authors = gsub(" et al\\.?", " et al.", authors),
      authors = gsub(" ?,$", "", authors),
      # don't start with 'sensu' or 'ehrenb'
      authors = gsub("^(sensu|Ehrenb.?) ", "", authors, ignore.case = TRUE),
      # no initials, only surname
      authors = gsub("^([A-Z]+ )+", "", authors, ignore.case = FALSE),
      # combine author and year if year is available
      ref = ifelse(!is.na(lastyear),
                   paste0(authors, ", ", lastyear),
                   authors),
      # fix beginning and ending
      ref = gsub(", $", "", ref),
      ref = gsub("^, ", "", ref)) %>%
    # remove text if it contains 'Not assigned' like phylum in viruses
    mutate_all(~gsub("Not assigned", "", .)) %>%
    # Remove non-ASCII characters (these are not allowed by CRAN)
    lapply(iconv, from = "UTF-8", to = "ASCII//TRANSLIT") %>%
    as_tibble(stringsAsFactors = FALSE) %>%
    select(col_id_new, fullname, ref, authors2) %>%
    left_join(new %>% select(col_id, fullname_new = fullname), by = c(col_id_new = "col_id")) %>%
    mutate(fullname = trimws(
      gsub("(.*)[(].*", "\\1",
           stringr::str_replace(
             string = fullname,
             pattern = stringr::fixed(authors2),
             replacement = "")) %>%
        gsub(" (var|f|subsp)[.]", "", .))) %>%
    select(-c("col_id_new", "authors2")) %>%
    filter(!is.na(fullname), !is.na(fullname_new)) %>%
    filter(fullname != fullname_new, !fullname %like% "^[?]")
}

# clean CoL and GBIF ----
# clean data_col
data_col <- data_col %>%
  as_tibble() %>%
  select(col_id = taxonID,
         col_id_new = acceptedNameUsageID,
         fullname = scientificName,
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
         species_id = furtherInformationURL) %>%
  mutate(source = "CoL")
# split into old and new
data_col.new <- data_col %>% clean_new()
data_col.old <- data_col %>% clean_old(new = data_col.new)
rm(data_col)

# clean data_gbif
data_gbif <- data_gbif %>%
  as_tibble() %>%
  filter(
    # no uncertain taxonomic placements
    taxonRemarks != "doubtful",
    kingdom != "incertae sedis",
    taxonRank != "unranked") %>%
  transmute(col_id = taxonID,
            col_id_new = acceptedNameUsageID,
            fullname = scientificName,
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
            species_id = as.character(parentNameUsageID)) %>%
  mutate(source = "GBIF")
# split into old and new
data_gbif.new <- data_gbif %>% clean_new()
data_gbif.old <- data_gbif %>% clean_old(new = data_gbif.new)
rm(data_gbif)

# put CoL and GBIF together ----
MOs.new <- bind_rows(data_col.new,
                     data_gbif.new) %>%
  mutate(taxonomic_tree_length = nchar(trimws(paste(kingdom, phylum, class, order, family, genus, species, subspecies)))) %>%
  arrange(desc(taxonomic_tree_length)) %>%
  distinct(fullname, .keep_all = TRUE) %>%
  select(-c("col_id_new", "authors2", "authors", "lastyear", "taxonomic_tree_length")) %>%
  arrange(fullname)
MOs.old <- bind_rows(data_col.old,
                     data_gbif.old) %>%
  distinct(fullname, .keep_all = TRUE) %>%
  arrange(fullname)

# clean up DSMZ ---
data_dsmz <- data_dsmz %>%
  as_tibble() %>%
  transmute(col_id = NA_integer_,
            col_id_new = NA_integer_,
            fullname = "",
            # kingdom = "",
            # phylum = "",
            # class = "",
            # order = "",
            # family = "",
            genus = ifelse(is.na(GENUS), "", GENUS),
            species = ifelse(is.na(SPECIES), "", SPECIES),
            subspecies = ifelse(is.na(SUBSPECIES), "", SUBSPECIES),
            rank = ifelse(species == "", "genus", "species"),
            ref = AUTHORS,
            species_id = as.character(RECORD_NO),
            source = "DSMZ")

# DSMZ only contains genus/(sub)species, try to find taxonomic properties based on genus and data_col
ref_taxonomy <- MOs.new %>%
  distinct(genus, .keep_all = TRUE) %>%
  filter(family != "") %>%
  filter(genus %in% data_dsmz$genus) %>%
  distinct(genus, .keep_all = TRUE) %>%
  select(kingdom, phylum, class, order, family, genus)

data_dsmz <- data_dsmz %>%
  left_join(ref_taxonomy, by = "genus") %>%
  mutate(kingdom = "Bacteria")

data_dsmz.new <- data_dsmz %>%
  clean_new() %>%
  distinct(fullname, .keep_all = TRUE) %>%
  select(colnames(MOs.new)) %>%
  arrange(fullname)

# combine everything ----
MOs <- bind_rows(MOs.new,
                 data_dsmz.new) %>%
  distinct(fullname, .keep_all = TRUE) %>%
  # not the ones that are old
  filter(!fullname %in% MOs.old$fullname) %>%
  arrange(fullname) %>%
  mutate(col_id = ifelse(source != "CoL", NA_integer_, col_id)) %>%
  filter(fullname != "")

rm(data_col.new)
rm(data_col.old)
rm(data_gbif.new)
rm(data_gbif.old)
rm(data_dsmz)
rm(data_dsmz.new)
rm(ref_taxonomy)
rm(MOs.new)

MOs.bak <- MOs

# Trichomonas trick ----
# for species in Trypanosoma and Trichomonas we observe al lot of taxonomic info missing
MOs %>% filter(genus %in% c("Trypanosoma", "Trichomonas")) %>% View()
MOs[which(MOs$genus == "Trypanosoma"), "kingdom"] <- MOs[which(MOs$fullname == "Trypanosoma"),]$kingdom
MOs[which(MOs$genus == "Trypanosoma"), "phylum"] <- MOs[which(MOs$fullname == "Trypanosoma"),]$phylum
MOs[which(MOs$genus == "Trypanosoma"), "class"] <- MOs[which(MOs$fullname == "Trypanosoma"),]$class
MOs[which(MOs$genus == "Trypanosoma"), "order"] <- MOs[which(MOs$fullname == "Trypanosoma"),]$order
MOs[which(MOs$genus == "Trypanosoma"), "family"] <- MOs[which(MOs$fullname == "Trypanosoma"),]$family
MOs[which(MOs$genus == "Trichomonas"), "kingdom"] <- MOs[which(MOs$fullname == "Trichomonas"),]$kingdom
MOs[which(MOs$genus == "Trichomonas"), "phylum"] <- MOs[which(MOs$fullname == "Trichomonas"),]$phylum
MOs[which(MOs$genus == "Trichomonas"), "class"] <- MOs[which(MOs$fullname == "Trichomonas"),]$class
MOs[which(MOs$genus == "Trichomonas"), "order"] <- MOs[which(MOs$fullname == "Trichomonas"),]$order
MOs[which(MOs$genus == "Trichomonas"), "family"] <- MOs[which(MOs$fullname == "Trichomonas"),]$family

# fill taxonomic properties that are missing
MOs <- MOs %>%
  mutate(phylum = ifelse(phylum %in% c(NA, ""), "(unknown phylum)", phylum),
         class = ifelse(class %in% c(NA, ""), "(unknown class)", class),
         order = ifelse(order %in% c(NA, ""), "(unknown order)", order),
         family = ifelse(family %in% c(NA, ""), "(unknown family)", family))

# Abbreviations ----
# Add abbreviations so we can easily know which ones are which ones.
# These will become valid and unique microbial IDs for the AMR package.
MOs <- MOs %>%
  arrange(kingdom, fullname) %>%
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
  # abbreviations determined per kingdom and family
  # becuase they are part of the abbreviation
  mutate(abbr_genus = abbreviate(genus,
                                 minlength = 7,
                                 use.classes = TRUE,
                                 method = "both.sides",
                                 strict = FALSE)) %>%
  ungroup() %>%
  group_by(genus) %>%
  # species abbreviations may be the same between genera
  # because the genus abbreviation is part of the abbreviation
  mutate(abbr_species = abbreviate(stringr::str_to_title(species),
                                   minlength = 3,
                                   use.classes = FALSE,
                                   method = "both.sides")) %>%
  ungroup() %>%
  group_by(genus, species) %>%
  mutate(abbr_subspecies = abbreviate(stringr::str_to_title(subspecies),
                                      minlength = 3,
                                      use.classes = FALSE,
                                      method = "both.sides")) %>%
  ungroup() %>%
  # remove trailing underscores
  mutate(mo = gsub("_+$", "",
                   toupper(paste(
                     # first character: kingdom
                     ifelse(kingdom %in% c("Animalia", "Plantae"),
                            substr(kingdom, 1, 2),
                            substr(kingdom, 1, 1)),
                     # next: genus, species, subspecies
                     ifelse(is.na(abbr_other),
                            paste(abbr_genus,
                                  abbr_species,
                                  abbr_subspecies,
                                  sep = "_"),
                            abbr_other),
                     sep = "_")))) %>%
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
               col_id = NA_integer_,
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
               stringsAsFactors = FALSE),
    data.frame(mo = "B_GRAMN",
               col_id = NA_integer_,
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
               stringsAsFactors = FALSE),
    data.frame(mo = "B_GRAMP",
               col_id = NA_integer_,
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
               stringsAsFactors = FALSE),
    # CoNS
    MOs %>%
      filter(genus == "Staphylococcus", species == "") %>% .[1,] %>%
      mutate(mo = paste(mo, "CNS", sep = "_"),
             rank = "species",
             col_id = NA_integer_,
             species = "coagulase-negative",
             fullname = "Coagulase-negative Staphylococcus (CoNS)",
             ref = NA_character_,
             species_id = "",
             source = "manually added"),
    # CoPS
    MOs %>%
      filter(genus == "Staphylococcus", species == "") %>% .[1,] %>%
      mutate(mo = paste(mo, "CPS", sep = "_"),
             rank = "species",
             col_id = NA_integer_,
             species = "coagulase-positive",
             fullname = "Coagulase-positive Staphylococcus (CoPS)",
             ref = NA_character_,
             species_id = "",
             source = "manually added"),
    # Streptococci groups A, B, C, F, H, K
    MOs %>%
      filter(genus == "Streptococcus", species == "pyogenes") %>% .[1,] %>%
      # we can keep all other details, since S. pyogenes is the only member of group A
      mutate(mo = paste(MOs[MOs$fullname == "Streptococcus",]$mo, "GRA", sep = "_"),
             species = "group A" ,
             fullname = "Streptococcus group A"),
    MOs %>%
      filter(genus == "Streptococcus", species == "agalactiae") %>% .[1,] %>%
      # we can keep all other details, since S. agalactiae is the only member of group B
      mutate(mo = paste(MOs[MOs$fullname == "Streptococcus",]$mo, "GRB", sep = "_"),
             species = "group B" ,
             fullname = "Streptococcus group B"),
    MOs %>%
      filter(genus == "Streptococcus", species == "dysgalactiae") %>% .[1,] %>%
      mutate(mo = paste(MOs[MOs$fullname == "Streptococcus",]$mo, "GRC", sep = "_"),
             col_id = NA_integer_,
             species = "group C" ,
             fullname = "Streptococcus group C",
             ref = NA_character_,
             species_id = "",
             source = "manually added"),
    MOs %>%
      filter(genus == "Streptococcus", species == "agalactiae") %>% .[1,] %>%
      mutate(mo = paste(MOs[MOs$fullname == "Streptococcus",]$mo, "GRD", sep = "_"),
             col_id = NA_integer_,
             species = "group D" ,
             fullname = "Streptococcus group D",
             ref = NA_character_,
             species_id = "",
             source = "manually added"),
    MOs %>%
      filter(genus == "Streptococcus", species == "agalactiae") %>% .[1,] %>%
      mutate(mo = paste(MOs[MOs$fullname == "Streptococcus",]$mo, "GRF", sep = "_"),
             col_id = NA_integer_,
             species = "group F" ,
             fullname = "Streptococcus group F",
             ref = NA_character_,
             species_id = "",
             source = "manually added"),
    MOs %>%
      filter(genus == "Streptococcus", species == "agalactiae") %>% .[1,] %>%
      mutate(mo = paste(MOs[MOs$fullname == "Streptococcus",]$mo, "GRG", sep = "_"),
             col_id = NA_integer_,
             species = "group G" ,
             fullname = "Streptococcus group G",
             ref = NA_character_,
             species_id = "",
             source = "manually added"),
    MOs %>%
      filter(genus == "Streptococcus", species == "agalactiae") %>% .[1,] %>%
      mutate(mo = paste(MOs[MOs$fullname == "Streptococcus",]$mo, "GRH", sep = "_"),
             col_id = NA_integer_,
             species = "group H" ,
             fullname = "Streptococcus group H",
             ref = NA_character_,
             species_id = "",
             source = "manually added"),
    MOs %>%
      filter(genus == "Streptococcus", species == "agalactiae") %>% .[1,] %>%
      mutate(mo = paste(MOs[MOs$fullname == "Streptococcus",]$mo, "GRK", sep = "_"),
             col_id = NA_integer_,
             species = "group K" ,
             fullname = "Streptococcus group K",
             ref = NA_character_,
             species_id = "",
             source = "manually added"),
    # Beta-haemolytic Streptococci
    MOs %>%
      filter(genus == "Streptococcus", species == "agalactiae") %>% .[1,] %>%
      mutate(mo = paste(MOs[MOs$fullname == "Streptococcus",]$mo, "HAE", sep = "_"),
             col_id = NA_integer_,
             species = "beta-haemolytic" ,
             fullname = "Beta-haemolytic Streptococcus",
             ref = NA_character_,
             species_id = "",
             source = "manually added")
  )


# everything distinct?
sum(duplicated(MOs$mo))
colnames(MOs)

# set prevalence per species
MOs <- MOs %>%
  mutate(prevalence = case_when(
    class == "Gammaproteobacteria"
    | genus %in% c("Enterococcus", "Staphylococcus", "Streptococcus")
    | mo %in% c("UNKNOWN", "B_GRAMN", "B_GRAMP")
    ~ 1,
    phylum %in% c("Proteobacteria",
                  "Firmicutes",
                  "Actinobacteria",
                  "Sarcomastigophora")
    | genus %in% c("Aspergillus",
                   "Bacteroides",
                   "Candida",
                   "Capnocytophaga",
                   "Chryseobacterium",
                   "Cryptococcus",
                   "Elisabethkingia",
                   "Flavobacterium",
                   "Fusobacterium",
                   "Giardia",
                   "Leptotrichia",
                   "Mycoplasma",
                   "Prevotella",
                   "Rhodotorula",
                   "Treponema",
                   "Trichophyton",
                   "Trichomonas",
                   "Ureaplasma")
    | rank %in% c("kingdom", "phylum", "class", "order", "family")
    ~ 2,
    TRUE ~ 3
  ))

# arrange
MOs <- MOs %>% arrange(fullname)

# transform
MOs <- as.data.frame(MOs, stringsAsFactors = FALSE)
MOs.old <- as.data.frame(MOs.old, stringsAsFactors = FALSE)
class(MOs$mo) <- "mo"
MOs$col_id <- as.integer(MOs$col_id)

# get differences in MO codes between this data and the package version
MO_diff <- AMR::microorganisms %>%
  mutate(pastedtext = paste(mo, fullname)) %>%
  filter(!pastedtext %in% (MOs %>% mutate(pastedtext = paste(mo, fullname)) %>% pull(pastedtext))) %>%
  select(mo_old = mo, fullname, pastedtext) %>%
  left_join(MOs %>%
              transmute(mo_new = mo, fullname_new = fullname, pastedtext = paste(mo, fullname)), "pastedtext") %>%
  select(mo_old, mo_new, fullname_new)

mo_diff2 <- AMR::microorganisms %>%
  select(mo, fullname) %>%
  left_join(MOs %>%
              select(mo, fullname),
            by = "fullname",
            suffix = c("_old", "_new")) %>%
  filter(mo_old != mo_new,
         #!mo_new %in% mo_old,
         !mo_old %like% "\\[")

mo_diff3 <- tibble(previous_old = names(AMR:::make_trans_tbl()),
                   previous_new = AMR:::make_trans_tbl()) %>%
  left_join(AMR::microorganisms %>% select(mo, fullname), by = c(previous_new = "mo")) %>%
  left_join(MOs %>% select(mo_new = mo, fullname), by = "fullname")

# what did we win most?
MOs %>% filter(!fullname %in% AMR::microorganisms$fullname) %>% freq(genus)
# what did we lose most?
AMR::microorganisms %>%
  filter(kingdom !=  "Chromista" & !fullname %in% MOs$fullname & !fullname %in% MOs.old$fullname) %>%
  freq(genus)


# save
saveRDS(MOs, "microorganisms.rds")
saveRDS(MOs.old, "microorganisms.old.rds")

# on the server, do:
usethis::use_data(microorganisms, overwrite = TRUE, version = 2)
usethis::use_data(microorganisms.old, overwrite = TRUE, version = 2)
rm(microorganisms)
rm(microorganisms.old)

# TO DO AFTER THIS
# * Update the year and dim()s in R/data.R
# * Rerun data-raw/reproduction_of_rsi_translation.R
# * Run unit tests
