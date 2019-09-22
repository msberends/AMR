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
data_col <- data.table::fread("data-raw/taxon.tab")

# read the xlsx file from DSMZ (only around 2.5 MB):
data_dsmz <- readxl::read_xlsx("data-raw/DSMZ_bactnames.xlsx")

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
    # or the taxonomic entry is old - the species was renamed
    | !is.na(col_id_new)
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

# filter old taxonomic names so only the ones with an existing reference will be kept
MOs <- MOs %>%
  filter(is.na(col_id_new) | (!is.na(col_id_new) & col_id_new %in% MOs$col_id))

MOs <- MOs %>%
  # remove text if it contains 'Not assigned' like phylum in viruses
  mutate_all(~gsub("(Not assigned|\\[homonym\\]|\\[mistake\\])", "", ., ignore.case = TRUE))

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
  as_tibble(stringsAsFactors = FALSE) %>% 
  # remove invalid characters
  mutate_all(~gsub("[\"'`]+", "", .))

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

MO.bak <- MOs

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
  filter(!grepl("unassigned", fullname, ignore.case = TRUE)) %>% 
  # prefer DSMZ over CoL, since that's more recent
  arrange(desc(source)) %>% 
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

# add CoL's col_id, source and ref from MOs.bak, for the cases where DSMZ took preference
MOs <- MOs %>% 
  mutate(kingdom_fullname = paste(kingdom, fullname)) %>% 
  left_join(MO.bak %>%
              filter(is.na(col_id_new), !is.na(col_id)) %>%
              transmute(col_id, species_id, source, ref, kingdom_fullname = trimws(paste(kingdom, genus, species, subspecies))), 
            by = "kingdom_fullname",
            suffix = c("_dsmz", "_col")) %>% 
  mutate(col_id = col_id_col, 
         species_id = ifelse(!is.na(species_id_col), gsub(".*/(.*)$", "\\1", species_id_col), species_id_dsmz),
         source = ifelse(!is.na(species_id_col), source_col, source_dsmz), 
         ref = ifelse(!is.na(species_id_col) & ref_col != "", ref_col, ref_dsmz)) %>% 
  select(-matches("(_col|_dsmz|kingdom_fullname)"))


MOs.old <- MOs.old %>%
  # remove the ones that are in the MOs data set
  filter(col_id_new %in% MOs$col_id) %>% 
  # and remove the ones that have the exact same fullname in the MOs data set, like Moraxella catarrhalis
  left_join(MOs, by = "fullname") %>%
  filter(col_id_new != col_id.y | is.na(col_id.y)) %>% 
  select(col_id = col_id.x, col_id_new, fullname, ref = ref.x)

# remove the records that are in MOs.old
sum(MOs.old$fullname %in% MOs$fullname)
MOs <- MOs %>% filter(!fullname %in% MOs.old$fullname)
sum(MOs.old$fullname %in% MOs$fullname)

# what characters are in the fullnames?
table(sort(unlist(strsplit(x = paste(MOs$fullname, collapse = ""), split = ""))))

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
               prevalence = 1,
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
               prevalence = 1,
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
               prevalence = 1,
               stringsAsFactors = FALSE),
    data.frame(mo = "F_YEAST",
               col_id = NA_integer_,
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
               col_id = NA_integer_,
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
             col_id = NA_integer_,
             species = "coagulase-negative",
             fullname = "Coagulase-negative Staphylococcus (CoNS)",
             ref = NA_character_,
             species_id = "",
             source = "manually added"),
    # CoPS
    MOs %>%
      filter(genus == "Staphylococcus", species == "epidermidis") %>% .[1,] %>%
      mutate(mo = gsub("(.*)_(.*)_.*", "\\1_\\2_COPS", mo),
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
             col_id = NA_integer_,
             species = "group C" ,
             fullname = "Streptococcus group C",
             ref = NA_character_,
             species_id = "",
             source = "manually added"),
    MOs %>%
      filter(genus == "Streptococcus", species == "agalactiae") %>% .[1,] %>%
      mutate(mo = gsub("(.*)_(.*)_.*", "\\1_\\2_GRPD", mo),
             col_id = NA_integer_,
             species = "group D" ,
             fullname = "Streptococcus group D",
             ref = NA_character_,
             species_id = "",
             source = "manually added"),
    MOs %>%
      filter(genus == "Streptococcus", species == "agalactiae") %>% .[1,] %>%
      mutate(mo = gsub("(.*)_(.*)_.*", "\\1_\\2_GRPF", mo),
             col_id = NA_integer_,
             species = "group F" ,
             fullname = "Streptococcus group F",
             ref = NA_character_,
             species_id = "",
             source = "manually added"),
    MOs %>%
      filter(genus == "Streptococcus", species == "agalactiae") %>% .[1,] %>%
      mutate(mo = gsub("(.*)_(.*)_.*", "\\1_\\2_GRPG", mo),
             col_id = NA_integer_,
             species = "group G" ,
             fullname = "Streptococcus group G",
             ref = NA_character_,
             species_id = "",
             source = "manually added"),
    MOs %>%
      filter(genus == "Streptococcus", species == "agalactiae") %>% .[1,] %>%
      mutate(mo = gsub("(.*)_(.*)_.*", "\\1_\\2_GRPH", mo),
             col_id = NA_integer_,
             species = "group H" ,
             fullname = "Streptococcus group H",
             ref = NA_character_,
             species_id = "",
             source = "manually added"),
    MOs %>%
      filter(genus == "Streptococcus", species == "agalactiae") %>% .[1,] %>%
      mutate(mo = gsub("(.*)_(.*)_.*", "\\1_\\2_GRPK", mo),
             col_id = NA_integer_,
             species = "group K" ,
             fullname = "Streptococcus group K",
             ref = NA_character_,
             species_id = "",
             source = "manually added"),
    # Beta haemolytic Streptococci
    MOs %>%
      filter(genus == "Streptococcus", species == "agalactiae") %>% .[1,] %>%
      mutate(mo = gsub("(.*)_(.*)_.*", "\\1_\\2_HAEM", mo),
             col_id = NA_integer_,
             species = "beta-haemolytic" ,
             fullname = "Beta-haemolytic Streptococcus",
             ref = NA_character_,
             species_id = "",
             source = "manually added"),
    # Viridans Streptococci
    MOs %>%
      filter(genus == "Streptococcus", species == "agalactiae") %>% .[1,] %>%
      mutate(mo = gsub("(.*)_(.*)_.*", "\\1_\\2_VIRI", mo),
             col_id = NA_integer_,
             species = "viridans" ,
             fullname = "Viridans Group Streptococcus (VGS)",
             ref = NA_character_,
             species_id = "",
             source = "manually added"),
    # Milleri Streptococci
    MOs %>%
      filter(genus == "Streptococcus", species == "agalactiae") %>% .[1,] %>%
      mutate(mo = gsub("(.*)_(.*)_.*", "\\1_\\2_MILL", mo),
             col_id = NA_integer_,
             species = "milleri" ,
             fullname = "Milleri Group Streptococcus (MGS)",
             ref = NA_character_,
             species_id = "",
             source = "manually added"),
    # Blastocystis hominis does not exist (it means 'got a Bastocystis from humans', PMID 15634993)
    # but let's be nice to the clinical people in microbiology
    MOs %>%
      filter(fullname == "Blastocystis") %>%
      mutate(mo = paste0(mo, "_HMNS"),
             fullname = paste(fullname, "hominis"),
             species = "hominis",
             col_id = NA,
             source = "manually added",
             ref = NA_character_,
             species_id = ""),
    # Trichomonas vaginalis is missing, same order as Dientamoeba
    MOs %>%
      filter(fullname == "Dientamoeba") %>%
      mutate(mo = gsub("(.*?)_.*", "\\1_THMNS", mo),
             col_id = NA,
             fullname = "Trichomonas",
             family = "Trichomonadidae",
             genus = "Trichomonas",
             source = "manually added",
             ref = "Donne, 1836",
             species_id = ""),
    MOs %>%
      filter(fullname == "Dientamoeba fragilis") %>%
      mutate(mo = gsub("(.*?)_.*", "\\1_THMNS_VAG", mo),
             col_id = NA,
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
             col_id = NA,
             fullname = "Trichomonadidae",
             family = "Trichomonadidae",
             rank = "family",
             genus = "",
             species = "",
             source = "manually added",
             ref = "",
             species_id = ""),
  )

MOs <- MOs %>%
  group_by(kingdom) %>%
  distinct(fullname, .keep_all = TRUE) %>% 
  ungroup() %>% 
  filter(fullname != "")

# add prevalence to old taxonomic names
MOs.old <- MOs.old %>% 
  left_join(MOs %>% select(col_id, prevalence), by = c("col_id_new" = "col_id"))

# everything distinct?
sum(duplicated(MOs$mo))
sum(duplicated(MOs$fullname))
colnames(MOs)

# here we welcome the new ones:
MOs %>% arrange(fullname) %>% filter(!fullname %in% AMR::microorganisms$fullname) %>% View()
MOs.old %>% arrange(fullname) %>% filter(!fullname %in% AMR::microorganisms.old$fullname) %>% View()
# and the ones we lost:
AMR::microorganisms %>% filter(!fullname %in% MOs$fullname) %>% View()
# and these IDs have changed:
old_new <- MOs %>%
  mutate(kingdom_fullname = paste(kingdom, fullname)) %>% 
  filter(kingdom_fullname %in% (AMR::microorganisms %>% mutate(kingdom_fullname = paste(kingdom, fullname)) %>% pull(kingdom_fullname))) %>%
  left_join(AMR::microorganisms %>% mutate(kingdom_fullname = paste(kingdom, fullname)) %>% select(mo, kingdom_fullname), by = "kingdom_fullname", suffix = c("_new", "_old")) %>% 
  filter(mo_new != mo_old) %>% 
  select(mo_old, mo_new, everything())
old_new %>% 
  View()
# and these codes are now missing (which will throw a unit test error):
AMR::microorganisms.codes %>% filter(!mo %in% MOs$mo)
AMR::rsi_translation %>% filter(!mo %in% MOs$mo)
AMR::microorganisms.translation %>% filter(!mo_new %in% MOs$mo)
# this is how to fix it
microorganisms.codes <- AMR::microorganisms.codes %>% 
  left_join(MOs %>%
              mutate(kingdom_fullname = paste(kingdom, fullname)) %>% 
              left_join(AMR::microorganisms  %>%
                          mutate(kingdom_fullname = paste(kingdom, fullname)) %>% 
                          select(mo, kingdom_fullname), by = "kingdom_fullname", suffix = c("_new", "_old")) %>%
              select(mo_old, mo_new),
            by = c("mo" = "mo_old")) %>% 
  select(code, mo = mo_new) %>% 
  filter(!is.na(mo))
microorganisms.codes %>% filter(!mo %in% MOs$mo)

# arrange
MOs <- MOs %>% arrange(fullname)
MOs.old <- MOs.old %>% arrange(fullname)
microorganisms.codes <- microorganisms.codes %>% arrange(code)

# transform
MOs <- as.data.frame(MOs, stringsAsFactors = FALSE)
MOs.old <- as.data.frame(MOs.old, stringsAsFactors = FALSE)
microorganisms.codes <- as.data.frame(microorganisms.codes, stringsAsFactors = FALSE)
class(MOs$mo) <- "mo"
class(microorganisms.codes$mo) <- "mo"
MOs$col_id <- as.integer(MOs$col_id)
MOs.old$col_id <- as.integer(MOs.old$col_id)
MOs.old$col_id_new <- as.integer(MOs.old$col_id_new)

# SAVE
### for other server
saveRDS(MOs, "microorganisms.rds")
saveRDS(MOs.old, "microorganisms.old.rds")
saveRDS(microorganisms.codes, "microorganisms.codes.rds")
### for same server
microorganisms <- MOs
microorganisms.old <- MOs.old
microorganisms.translation <- old_new %>% select(mo_old, mo_new)
class(microorganisms.translation$mo_old) <- "mo"
class(microorganisms.translation$mo_new) <- "mo"

# on the server, do:
usethis::use_data(microorganisms, overwrite = TRUE, version = 2)
usethis::use_data(microorganisms.old, overwrite = TRUE, version = 2)
usethis::use_data(microorganisms.codes, overwrite = TRUE, version = 2)
saveRDS(microorganisms.translation, file = "data-raw/microorganisms.translation.rds", version = 2) # this one will be covered in data-raw/internals.R
rm(microorganisms)
rm(microorganisms.old)
rm(microorganisms.codes)
rm(microorganisms.translation)
devtools::load_all(".")

# TO DO AFTER THIS
# * Update the year and dim()s in R/data.R
# * Rerun data-raw/reproduction_of_rsi_translation.R
# * Run unit tests
