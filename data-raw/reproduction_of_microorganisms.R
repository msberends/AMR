# Reproduction of the `microorganisms` data set

# Data retrieved from the Catalogue of Life (CoL) through the Encyclopaedia of Life:
# https://opendata.eol.org/dataset/catalogue-of-life/
# (download the resource file with a name like "Catalogue of Life yyyy-mm-dd")
# and from the Leibniz Institute DSMZ-German Collection of Microorganisms and Cell Cultures
# https://www.dsmz.de/support/bacterial-nomenclature-up-to-date-downloads.html
# (download the latest "Complete List" as xlsx file)

library(dplyr)
library(AMR)

# unzip and extract taxon.tab (around 1.5 GB) from the CoL archive, then:
data_col <- data.table::fread("Downloads/taxon.tab")

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
         species_id = furtherInformationURL)
data_col$source <- "CoL"

# clean data_dsmz
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
ref_taxonomy <- data_col %>%
  filter(genus %in% data_dsmz$genus,
         family != "") %>%
  distinct(genus, .keep_all = TRUE) %>%
  select(kingdom, phylum, class, order, family, genus)

data_dsmz <- data_dsmz %>%
  left_join(ref_taxonomy, by = "genus") %>%
  mutate(kingdom = "Bacteria",
         phylum = ifelse(is.na(phylum), "(unknown phylum)", phylum),
         class = ifelse(is.na(class), "(unknown class)", class),
         order = ifelse(is.na(order), "(unknown order)", order),
         family = ifelse(is.na(family), "(unknown family)", family),
  )

# combine everything
data_total <- data_col %>%
  bind_rows(data_dsmz)

rm(data_col)
rm(data_dsmz)
rm(ref_taxonomy)

MOs <- data_total %>%
  filter(
    (
      # we only want all MICROorganisms and no viruses
      !kingdom %in% c("Animalia", "Plantae", "Viruses")
      # and not all fungi: Aspergillus, Candida, Trichphyton and Pneumocystis are the most important,
      # so only keep these orders from the fungi:
      & !(kingdom == "Fungi"
          & !order %in% c("Eurotiales", "Saccharomycetales", "Schizosaccharomycetales", "Tremellales", "Onygenales", "Pneumocystales"))
    )
    # or the genus has to be one of the genera we found in our hospitals last decades (Northern Netherlands, 2002-2018)
    | genus %in% c("Absidia", "Acremonium", "Actinotignum", "Alternaria", "Anaerosalibacter", "Ancylostoma", "Anisakis", "Apophysomyces",
                   "Arachnia", "Ascaris", "Aureobacterium", "Aureobasidium", "Balantidum", "Bilophilia", "Branhamella", "Brochontrix",
                   "Brugia", "Calymmatobacterium", "Catabacter", "Cdc", "Chilomastix", "Chryseomonas", "Cladophialophora", "Cladosporium",
                   "Clonorchis", "Cordylobia", "Curvularia", "Demodex", "Dermatobia", "Diphyllobothrium", "Dracunculus", "Echinococcus",
                   "Enterobius", "Euascomycetes", "Exophiala", "Fasciola", "Fusarium", "Hendersonula", "Hymenolepis", "Kloeckera",
                   "Koserella", "Larva", "Leishmania", "Lelliottia", "Loa", "Lumbricus", "Malassezia", "Metagonimus", "Molonomonas",
                   "Mucor", "Nattrassia", "Necator", "Novospingobium", "Onchocerca", "Opistorchis", "Paragonimus", "Paramyxovirus",
                   "Pediculus", "Phoma", "Phthirus", "Pityrosporum", "Pseudallescheria", "Pulex", "Rhizomucor", "Rhizopus", "Rhodotorula",
                   "Salinococcus", "Sanguibacteroides", "Schistosoma", "Scopulariopsis", "Scytalidium", "Sporobolomyces", "Stomatococcus",
                   "Strongyloides", "Syncephalastraceae", "Taenia", "Torulopsis", "Trichinella", "Trichobilharzia", "Trichomonas",
                   "Trichosporon", "Trichuris", "Trypanosoma", "Wuchereria")
    # or the taxonomic entry is old - the species was renamed
    | !is.na(col_id_new)
  )

# filter old taxonomic names so only the ones with an existing reference will be kept
MOs <- MOs %>%
  filter(is.na(col_id_new) | (!is.na(col_id_new) & col_id_new %in% MOs$col_id))

MOs <- MOs %>%
  # remove text if it contains 'Not assigned' like phylum in viruses
  mutate_all(~gsub("Not assigned", "", .))

MOs <- MOs %>%
  # Only keep first author, e.g. transform 'Smith, Jones, 2011' to 'Smith et al., 2011':
  mutate(authors2 = iconv(ref, from = "UTF-8", to = "ASCII//TRANSLIT"),
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
         ref = gsub("^, ", "", ref)
  )

# Remove non-ASCII characters (these are not allowed by CRAN)
MOs <- MOs %>%
  lapply(iconv, from = "UTF-8", to = "ASCII//TRANSLIT") %>%
  as_tibble(stringsAsFactors = FALSE)

# Split old taxonomic names - they refer in the original data to a new `taxonID` with `acceptedNameUsageID`
MOs.old <- MOs %>%
  filter(!is.na(col_id_new),
         ref != "",
         source != "DSMZ") %>%
  transmute(col_id,
            col_id_new,
            fullname =
              trimws(
                gsub("(.*)[(].*", "\\1",
                     stringr::str_replace(
                       string = fullname,
                       pattern = stringr::fixed(authors2),
                       replacement = "")) %>%
                  gsub(" (var|f|subsp)[.]", "", .)),
            ref) %>%
  filter(!is.na(fullname)) %>%
  distinct(fullname, .keep_all = TRUE) %>%
  arrange(col_id)

MOs <- MOs %>%
  filter(is.na(col_id_new) | source == "DSMZ") %>%
  transmute(col_id,
            fullname = trimws(case_when(rank == "family" ~ family,
                                        rank == "order" ~ order,
                                        rank == "class" ~ class,
                                        rank == "phylum" ~ phylum,
                                        rank == "kingdom" ~ kingdom,
                                        TRUE ~ paste(genus, species, subspecies))),
            kingdom,
            phylum,
            class,
            order,
            family,
            genus = gsub(":", "", genus),
            species,
            subspecies,
            rank,
            ref,
            species_id = gsub(".*/([a-f0-9]+)", "\\1", species_id),
            source) %>%
  #distinct(fullname, .keep_all = TRUE) %>%
  filter(!grepl("unassigned", fullname, ignore.case = TRUE))

# Filter out the DSMZ records that were renamed and are now in MOs.old
MOs <- MOs %>%
  filter(!(source == "DSMZ" & fullname %in% MOs.old$fullname),
         !(source == "DSMZ" & fullname %in% (MOs %>% filter(source == "CoL") %>% pull(fullname)))) %>%
  distinct(fullname, .keep_all = TRUE)

# Add abbreviations so we can easily know which ones are which ones.
# These will become valid and unique microbial IDs for the AMR package.
MOs <- MOs %>%
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
  mutate(abbr_genus = abbreviate(genus,
                                 minlength = 5,
                                 use.classes = TRUE,
                                 method = "both.sides",
                                 strict = FALSE)) %>%
  ungroup() %>%
  group_by(genus) %>%
  # species abbreviations may be the same between genera
  # because the genus abbreviation is part of the abbreviation
  mutate(abbr_species = abbreviate(species,
                                   minlength = 3,
                                   use.classes = FALSE,
                                   method = "both.sides")) %>%
  ungroup() %>%
  group_by(genus, species) %>%
  mutate(abbr_subspecies = abbreviate(subspecies,
                                      minlength = 3,
                                      use.classes = FALSE,
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
      filter(genus == "Staphylococcus", species == "epidermidis") %>% .[1,] %>%
      mutate(mo = gsub("EPI", "CNS", mo),
             col_id = NA_integer_,
             species = "coagulase-negative",
             fullname = "Coagulase-negative Staphylococcus (CoNS)",
             ref = NA_character_,
             species_id = "",
             source = "manually added"),
    # CoPS
    MOs %>%
      filter(genus == "Staphylococcus", species == "epidermidis") %>% .[1,] %>%
      mutate(mo = gsub("EPI", "CPS", mo),
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
      mutate(mo = gsub("PYO", "GRA", mo),
             species = "group A" ,
             fullname = "Streptococcus group A"),
    MOs %>%
      filter(genus == "Streptococcus", species == "agalactiae") %>% .[1,] %>%
      # we can keep all other details, since S. agalactiae is the only member of group B
      mutate(mo = gsub("AGA", "GRB", mo),
             species = "group B" ,
             fullname = "Streptococcus group B"),
    MOs %>%
      filter(genus == "Streptococcus", species == "dysgalactiae") %>% .[1,] %>%
      mutate(mo = gsub("DYS", "GRC", mo),
             col_id = NA_integer_,
             species = "group C" ,
             fullname = "Streptococcus group C",
             ref = NA_character_,
             species_id = "",
             source = "manually added"),
    MOs %>%
      filter(genus == "Streptococcus", species == "agalactiae") %>% .[1,] %>%
      mutate(mo = gsub("AGA", "GRD", mo),
             col_id = NA_integer_,
             species = "group D" ,
             fullname = "Streptococcus group D",
             ref = NA_character_,
             species_id = "",
             source = "manually added"),
    MOs %>%
      filter(genus == "Streptococcus", species == "agalactiae") %>% .[1,] %>%
      mutate(mo = gsub("AGA", "GRF", mo),
             col_id = NA_integer_,
             species = "group F" ,
             fullname = "Streptococcus group F",
             ref = NA_character_,
             species_id = "",
             source = "manually added"),
    MOs %>%
      filter(genus == "Streptococcus", species == "agalactiae") %>% .[1,] %>%
      mutate(mo = gsub("AGA", "GRG", mo),
             col_id = NA_integer_,
             species = "group G" ,
             fullname = "Streptococcus group G",
             ref = NA_character_,
             species_id = "",
             source = "manually added"),
    MOs %>%
      filter(genus == "Streptococcus", species == "agalactiae") %>% .[1,] %>%
      mutate(mo = gsub("AGA", "GRH", mo),
             col_id = NA_integer_,
             species = "group H" ,
             fullname = "Streptococcus group H",
             ref = NA_character_,
             species_id = "",
             source = "manually added"),
    MOs %>%
      filter(genus == "Streptococcus", species == "agalactiae") %>% .[1,] %>%
      mutate(mo = gsub("AGA", "GRK", mo),
             col_id = NA_integer_,
             species = "group K" ,
             fullname = "Streptococcus group K",
             ref = NA_character_,
             species_id = "",
             source = "manually added"),
    # Beta haemolytic Streptococci
    MOs %>%
      filter(genus == "Streptococcus", species == "agalactiae") %>% .[1,] %>%
      mutate(mo = gsub("AGA", "HAE", mo),
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
                   "Ureaplasma")
    | rank %in% c("kingdom", "phylum", "class", "order", "family")
    ~ 2,
    TRUE ~ 3
  ))

# save it
MOs <- as.data.frame(MOs %>% arrange(fullname), stringsAsFactors = FALSE)
MOs.old <- as.data.frame(MOs.old, stringsAsFactors = FALSE)
class(MOs$mo) <- "mo"

saveRDS(MOs, "microorganisms.rds")
saveRDS(MOs.old, "microorganisms.old.rds")

# on the server:
usethis::use_data(microorganisms, overwrite = TRUE, version = 2)
usethis::use_data(microorganisms.old, overwrite = TRUE, version = 2)
rm(microorganisms)
rm(microorganisms.old)
# and update the year in R/data.R
