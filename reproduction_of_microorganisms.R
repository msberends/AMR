# Catalogue of Life
# Data retrieved from Encyclopaedia of Life:
# https://opendata.eol.org/dataset/catalogue-of-life/

# unzip and extract taxon.tab, then:
taxon <- data.table::fread("taxon.tab")
# result is over 3.7M rows
taxon %>% freq(kingdom)
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

MOs <- taxon %>%
  # tibble for future transformations
  as_tibble() %>%
  filter(
    # we only want all microorganisms and viruses
    !kingdom %in% c("Animalia", "Plantae"),
    # and no entries above genus - they all already have a taxonomic tree
    !taxonRank %in% c("kingdom", "phylum", "superfamily", "class", "order", "family"),
    # not all fungi: Aspergillus, Candida, Trichphyton and Pneumocystis are the most important,
    # so only keep these orders from the fungi:
    !(kingdom == "Fungi" & !order %in% c("Eurotiales", "Saccharomycetales", "Schizosaccharomycetales", "Onygenales", "Pneumocystales"))) %>%
  # remove text if it contains 'Not assigned' like phylum in viruses
  mutate_all(funs(gsub("Not assigned", "", .))) %>%
  # only latest ref, not original authors
  mutate(scientificNameAuthorship = trimws(gsub(".*[)] ", "", scientificNameAuthorship)),
         scientificNameAuthorship = ifelse(grepl(" emend[. ]", scientificNameAuthorship, ignore.case = TRUE),
                                           gsub("(.*)emend[. ]+(.*)", "\\2", scientificNameAuthorship, ignore.case = TRUE),
                                           scientificNameAuthorship),
         scientificNameAuthorship = gsub(".", "", scientificNameAuthorship, fixed = TRUE),
         scientificNameAuthorship = gsub(",? et al", " et al.", scientificNameAuthorship, fixed = FALSE, ignore.case = TRUE),
         scientificNameAuthorship = gsub("[()]", "", scientificNameAuthorship),
         # year always preceded by comma
         scientificNameAuthorship = gsub(" ([0-9]{4})$", ", \\1", scientificNameAuthorship),
         scientificNameAuthorship = gsub(",,", ",", scientificNameAuthorship, fixed = TRUE),
         # only first author with *et al.*
         scientificNameAuthorship = gsub(",.*,", " et al.,", scientificNameAuthorship),
         scientificNameAuthorship = gsub(" (and|&) .*,", " et al.,", scientificNameAuthorship),
         scientificNameAuthorship = gsub(", [^0-9]+", ", ", scientificNameAuthorship),
         scientificNameAuthorship = gsub(", $", "", scientificNameAuthorship)
  )

# remove non-ASCII characters (not allowed by CRAN)
MOs <- MOs %>%
  lapply(iconv, from = "UTF-8", to = "ASCII//TRANSLIT") %>%
  as_tibble(stringsAsFactors = FALSE)

# split old taxonomic names - they refer to a new `taxonID` with `acceptedNameUsageID`
MOs.old <- MOs %>%
  filter(!is.na(acceptedNameUsageID),
         scientificNameAuthorship != "") %>%
  transmute(col_id = taxonID,
            col_id_new = acceptedNameUsageID,
            fullname =
              trimws(
                gsub("(.*)[(].*", "\\1",
                     stringr::str_replace(
                       string = scientificName,
                       pattern = stringr::fixed(scientificNameAuthorship),
                       replacement = ""))),
            ref = scientificNameAuthorship) %>%
  filter(!is.na(fullname)) %>%
  distinct(fullname, .keep_all = TRUE) %>%
  arrange(col_id)

MOs <- MOs %>%
  filter(is.na(acceptedNameUsageID)) %>%
  transmute(col_id = taxonID,
            fullname = trimws(ifelse(kingdom == "Viruses",
                                     paste(specificEpithet, infraspecificEpithet),
                                     paste(genus, specificEpithet, infraspecificEpithet))),
            kingdom,
            phylum,
            class,
            order,
            family,
            genus = gsub(":", "", genus),
            species = specificEpithet,
            subspecies = infraspecificEpithet,
            rank = taxonRank,
            ref = scientificNameAuthorship,
            species_id = gsub(".*/([a-f0-9]+)", "\\1", furtherInformationURL)) %>%
  distinct(fullname, .keep_all = TRUE) %>%
  filter(!grepl("unassigned", fullname, ignore.case = TRUE))

# only old names of species that are in MOs:
MOs.old <- MOs.old %>% filter(col_id_new %in% MOs$col_id)

MOs <- MOs %>%
  group_by(kingdom) %>%
  # abbreviations may be same for genera between kingdoms,
  # because each abbreviation starts with the the first character of the kingdom
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
                   toupper(paste(substr(kingdom, 1, 1),
                                 abbr_genus,
                                 abbr_species,
                                 abbr_subspecies,
                                 sep = "_")))) %>%
  mutate(mo = ifelse(duplicated(.$mo), paste0(mo, "1"), mo)) %>%
  select(mo, everything(), -abbr_genus, -abbr_species, -abbr_subspecies)

# everything distinct?
sum(duplicated(MOs$mo))

# add non-taxonomic entries
MOs <- MOs %>%
  bind_rows(
    # CoNS
    MOs %>%
      filter(genus == "Staphylococcus", species == "epidermidis") %>% .[1,] %>%
      mutate(mo = gsub("EPI", "CNS", mo),
             col_id = NA_integer_,
             species = "coagulase negative",
             fullname = "Coagulase Negative Staphylococcus (CoNS)",
             ref = NA_character_),
    # CoPS
    MOs %>%
      filter(genus == "Staphylococcus", species == "epidermidis") %>% .[1,] %>%
      mutate(mo = gsub("EPI", "CPS", mo),
             col_id = NA_integer_,
             species = "coagulase positive",
             fullname = "Coagulase Positive Staphylococcus (CoPS)",
             ref = NA_character_),
    # Streptococci groups A, B, C, F, H, K
    MOs %>%
      filter(genus == "Streptococcus", species == "agalactiae") %>% .[1,] %>%
      mutate(mo = gsub("AGA", "GRA", mo),
             col_id = NA_integer_,
             species = "group A" ,
             fullname = "Streptococcus group A"),
    MOs %>%
      filter(genus == "Streptococcus", species == "dysgalactiae") %>% .[1,] %>%
      mutate(mo = gsub("DYS", "GRB", mo),
             col_id = NA_integer_,
             species = "group B" ,
             fullname = "Streptococcus group B"),
    MOs %>%
      filter(genus == "Streptococcus", species == "agalactiae") %>% .[1,] %>%
      mutate(mo = gsub("AGA", "GRC", mo),
             col_id = NA_integer_,
             species = "group C" ,
             fullname = "Streptococcus group C",
             ref = NA_character_),
    MOs %>%
      filter(genus == "Streptococcus", species == "agalactiae") %>% .[1,] %>%
      mutate(mo = gsub("AGA", "GRD", mo),
             col_id = NA_integer_,
             species = "group D" ,
             fullname = "Streptococcus group D",
             ref = NA_character_),
    MOs %>%
      filter(genus == "Streptococcus", species == "agalactiae") %>% .[1,] %>%
      mutate(mo = gsub("AGA", "GRF", mo),
             col_id = NA_integer_,
             species = "group F" ,
             fullname = "Streptococcus group F",
             ref = NA_character_),
    MOs %>%
      filter(genus == "Streptococcus", species == "agalactiae") %>% .[1,] %>%
      mutate(mo = gsub("AGA", "GRG", mo),
             col_id = NA_integer_,
             species = "group F" ,
             fullname = "Streptococcus group G",
             ref = NA_character_),
    MOs %>%
      filter(genus == "Streptococcus", species == "agalactiae") %>% .[1,] %>%
      mutate(mo = gsub("AGA", "GRH", mo),
             col_id = NA_integer_,
             species = "group H" ,
             fullname = "Streptococcus group H",
             ref = NA_character_),
    MOs %>%
      filter(genus == "Streptococcus", species == "agalactiae") %>% .[1,] %>%
      mutate(mo = gsub("AGA", "GRK", mo),
             col_id = NA_integer_,
             species = "group K" ,
             fullname = "Streptococcus group K",
             ref = NA_character_),
    # Beta haemolytic Streptococci
    MOs %>%
      filter(genus == "Streptococcus", species == "agalactiae") %>% .[1,] %>%
      mutate(mo = gsub("AGA", "HAE", mo),
             col_id = NA_integer_,
             species = "beta-haemolytic" ,
             fullname = "Beta-haemolytic Streptococcus",
             ref = NA_character_),
    # unknowns
    data.frame(mo = "B_GRAMN",
               col_id = NA_integer_,
               fullname = "(unknown Gram negatives)",
               kingdom = "Bacteria",
               phylum = NA_character_,
               class = NA_character_,
               order = NA_character_,
               family = NA_character_,
               genus = "(unknown Gram negatives)",
               species = NA_character_,
               subspecies = NA_character_,
               rank = "species",
               ref = NA_character_,
               stringsAsFactors = FALSE),
    data.frame(mo = "B_GRAMP",
               col_id = NA_integer_,
               fullname = "(unknown Gram positives)",
               kingdom = "Bacteria",
               phylum = NA_character_,
               class = NA_character_,
               order = NA_character_,
               family = NA_character_,
               genus = "(unknown Gram positives)",
               species = NA_character_,
               subspecies = NA_character_,
               rank = "species",
               ref = NA_character_,
               stringsAsFactors = FALSE)
  )

# save it
MOs <- as.data.frame(MOs %>% arrange(mo), stringsAsFactors = FALSE)
MOs.old <- as.data.frame(MOs.old, stringsAsFactors = FALSE)
class(MOs$mo) <- "mo"

saveRDS(MOs, "microorganisms.rds")
saveRDS(MOs.old, "microorganisms.old.rds")
# on the server:
# usethis::use_data(microorganisms, overwrite = TRUE)
# usethis::use_data(microorganisms.old, overwrite = TRUE)
