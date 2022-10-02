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
# doi:10.18637/jss.v104.i03                                            #
#                                                                      #
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

# THIS SCRIPT REQUIRES AT LEAST 16 GB RAM
# (at least 10 GB will be used by the R session for the size of the files)

# 1. Go to https://doi.org/10.15468/39omei and find the download link for the
#    latest GBIF backbone taxonony ZIP file - unpack Taxon.tsv from it (~2.2 GB)
#    ALSO BE SURE to get the date of release and update R/aa_globals.R later!
# 2. Go to https://lpsn.dsmz.de/downloads (register first) and download the latest
#    CSV file. It should be (re)named "taxonomy.csv". Their API unfortunately does
#    not include the full taxonomy and is currently (2022) pretty worthless.
# 3. Set this folder_location to the path where these two files are:
folder_location <- "~/Downloads/backbone/"
file_gbif <- paste0(folder_location, "Taxon.tsv")
file_lpsn <- paste0(folder_location, "taxonomy.csv")

# 4. Run the rest of this script line by line and check everything :)

if (!file.exists(file_gbif)) stop("GBIF file not found")
if (!file.exists(file_lpsn)) stop("LPSN file not found")

library(dplyr)
library(vroom)
library(AMR)
# also requires 'rvest' and 'progress' packages

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
  # replace parentheses with emend. to get the latest authors
  authors2 <- gsub("(", " emend. ", authors2, fixed = TRUE)
  authors2 <- gsub(")", "", authors2, fixed = TRUE)
  authors2 <- gsub(" +", " ", authors2)
  authors2 <- trimws(authors2)

  # get year from last 4 digits
  lastyear <- as.integer(gsub(".*([0-9]{4})$", "\\1", authors2))
  # can never be later than now
  lastyear <- ifelse(lastyear > as.integer(format(Sys.Date(), "%Y")),
    NA,
    lastyear
  )
  # get authors without last year
  authors <- gsub("(.*)[0-9]{4}$", "\\1", authors2)
  # not sure what this is
  authors <- gsub("(Saito)", "", authors, fixed = TRUE)
  authors <- gsub("(Oudem.)", "", authors, fixed = TRUE)
  # remove nonsense characters from names
  authors <- gsub("[^a-zA-Z,'&. -]", "", authors)
  # no initials, only surname
  authors <- gsub("[A-Z][.]", "", authors, ignore.case = FALSE)
  # remove trailing and leading spaces
  authors <- trimws(authors)
  # keep only the part after last 'emend.' to get the latest authors
  authors <- gsub(".*emend[.] ?", "", authors)
  # only keep first author and replace all others by 'et al'
  authors <- gsub("(,| and| et| &| ex| emend\\.?) .*", " et al.", authors)
  # et al. always with ending dot
  authors <- gsub(" et al\\.?", " et al.", authors)
  authors <- gsub(" ?,$", "", authors)
  # don't start with 'sensu' or 'ehrenb'
  authors <- gsub("^(sensu|Ehrenb.?|corrig.?) ", "", authors, ignore.case = TRUE)
  # no initials, only surname
  authors <- trimws(authors)
  authors <- gsub("^([A-Z][.])+( & ?)?", "", authors, ignore.case = FALSE)
  authors <- gsub("^([A-Z]+ )+", "", authors, ignore.case = FALSE)
  # remove dots
  authors <- gsub(".", "", authors, fixed = TRUE)
  authors <- gsub("et al", "et al.", authors, fixed = TRUE)
  authors[nchar(authors) <= 3] <- ""
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
  ref[ref == ""] <- NA_character_
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

abbreviate_mo <- function(x, minlength = 5, prefix = "", hyphen_as_space = FALSE, ...) {
  if (hyphen_as_space == TRUE) {
    x <- gsub("-", " ", x, fixed = TRUE)
  }
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

# to retrieve LPSN and authors from LPSN website
get_lpsn_and_author <- function(rank, name) {
  url <- paste0("https://lpsn.dsmz.de/", tolower(rank), "/", tolower(name))
  page_txt <- tryCatch(rvest::read_html(url), error = function(e) NULL)
  if (is.null(page_txt)) {
    warning("No LPSN found for ", tolower(rank), " '", name, "'")
    lpsn <- NA_character_
    ref <- NA_character_
  } else {
    page_txt <- page_txt %>%
      rvest::html_element("#detail-page") %>%
      rvest::html_text()
    lpsn <- gsub(".*Record number:[\r\n\t ]*([0-9]+).*", "\\1", page_txt, perl = FALSE)
    ref <- page_txt %>%
      gsub(".*?Name: (.*[0-9]{4}?).*", "\\1", ., perl = FALSE) %>%
      gsub(name, "", ., fixed = TRUE) %>%
      trimws()
  }
  c("lpsn" = lpsn, "ref" = ref)
}

# MB/ August 2022: useless, does not contain full taxonomy, e.g. LPSN::request(cred, category = "family") is empty.
# get_from_lpsn <- function (user, pw) {
#   if (!"LPSN" %in% rownames(utils::installed.packages())) {
#     stop("Install the official LPSN package for R using: install.packages('LPSN', repos = 'https://r-forge.r-project.org')")
#   }
#   cred <- LPSN::open_lpsn(user, pw)
#
#   lpsn_genus <- LPSN::request(cred, category = "genus")
#   message("Downloading genus data (n = ", lpsn_genus$count, ") from LPSN API...")
#   lpsn_genus <- as.data.frame(LPSN::retrieve(cred, category = "genus"))
#
#   lpsn_species <- LPSN::request(cred, category = "species")
#   message("Downloading species data (n = ", lpsn_species$count, ") from LPSN API...")
#   lpsn_species <- as.data.frame(LPSN::retrieve(cred, category = "species"))
#
#   lpsn_subspecies <- LPSN::request(cred, category = "subspecies")
#   message("Downloading subspecies data (n = ", lpsn_subspecies$count, ") from LPSN API...")
#   lpsn_subspecies <- as.data.frame(LPSN::retrieve(cred, category = "subspecies"))
#
#   message("Binding rows...")
#   lpsn_total <- bind_rows(lpsn_genus, lpsn_species, lpsn_subspecies)
#   message("Done.")
#   lpsn_total
# }

# Read GBIF data ----------------------------------------------------------

taxonomy_gbif.bak <- vroom(file_gbif)
taxonomy_gbif <- taxonomy_gbif.bak %>%
  # immediately filter rows we really never want
  filter(
    # never doubtful status, only accepted and all synonyms, and only ranked items
    taxonomicStatus != "doubtful",
    taxonRank != "unranked",
    #  include these kingdoms (no Chromista)
    kingdom %in% c("Archaea", "Bacteria", "Protozoa") |
      # include all of these fungal orders
      order %in% c(
        "Eurotiales", "Microascales", "Mucorales", "Saccharomycetales",
        "Schizosaccharomycetales", "Tremellales", "Onygenales", "Pneumocystales"
      ) |
      # and all of these important genera (see "data-raw/_pre_commit_hook.R")
      # (they also contain bacteria and protozoa, but these will get prevalence = 2 later on)
      genus %in% AMR:::MO_PREVALENT_GENERA
  ) %>%
  select(
    kingdom,
    phylum,
    class,
    order,
    family,
    genus,
    species = specificEpithet,
    subspecies = infraspecificEpithet,
    rank = taxonRank,
    status = taxonomicStatus,
    ref = scientificNameAuthorship,
    gbif = taxonID,
    gbif_parent = parentNameUsageID,
    gbif_renamed_to = acceptedNameUsageID
  ) %>%
  mutate(
    # do this mutate after the original selection/filtering, as it decreases computing time tremendously
    status = ifelse(status == "accepted", "accepted", "synonym"),
    # checked taxonRank - the "form" and "variety" always have a subspecies, so:
    rank = ifelse(rank %in% c("form", "variety"), "subspecies", rank),
    source = "GBIF"
  ) %>%
  filter(
    # their data is messy - keep only these:
    rank == "kingdom" & !is.na(kingdom) |
      rank == "phylum" & !is.na(phylum) |
      rank == "class" & !is.na(class) |
      rank == "order" & !is.na(order) |
      rank == "family" & !is.na(family) |
      rank == "genus" & !is.na(genus) |
      rank == "species" & !is.na(species) |
      rank == "subspecies" & !is.na(subspecies)
  ) %>%
  # some items end with _A or _B... why??
  mutate_all(~ gsub("_[A-Z]$", "", .x, perl = TRUE)) %>%
  # now we have duplicates, remove these, but prioritise "accepted" status and highest taxon ID
  arrange(status, gbif) %>%
  distinct(kingdom, phylum, class, order, family, genus, species, subspecies, .keep_all = TRUE) %>%
  filter(
    kingdom %unlike% "[0-9]",
    phylum %unlike% "[0-9]",
    class %unlike% "[0-9]",
    order %unlike% "[0-9]",
    family %unlike% "[0-9]",
    genus %unlike% "[0-9]"
  )

# integrity tests
sort(table(taxonomy_gbif$rank))
sort(table(taxonomy_gbif$status))


# Read LPSN data ----------------------------------------------------------

taxonomy_lpsn.bak <- vroom(file_lpsn)
taxonomy_lpsn <- taxonomy_lpsn.bak %>%
  transmute(
    genus = genus_name,
    species = sp_epithet,
    subspecies = subsp_epithet,
    rank = case_when(
      !is.na(subsp_epithet) ~ "subspecies",
      !is.na(sp_epithet) ~ "species",
      TRUE ~ "genus"
    ),
    status = ifelse(is.na(record_lnk), "accepted", "synonym"),
    ref = authors,
    lpsn = as.character(record_no),
    lpsn_parent = NA_character_,
    lpsn_renamed_to = as.character(record_lnk)
  ) %>%
  mutate(source = "LPSN")
taxonomy_lpsn

# download additional taxonomy to the domain/kingdom level (their API is not sufficient...)
taxonomy_lpsn_missing <- tibble(
  kingdom = character(0),
  phylum = character(0),
  class = character(0),
  order = character(0),
  family = character(0),
  genus = character(0)
)
for (page in LETTERS) {
  message("Downloading page ", page, appendLF = FALSE)
  url <- paste0("https://lpsn.dsmz.de/genus?page=", page)
  x <- xml2::read_html(url) %>%
    # class "main-list" is the main table
    rvest::html_element(".main-list") %>%
    # get every list element with a set <id> attribute
    rvest::html_elements("li[id]")
  for (i in seq_len(length(x))) {
    if (i %% 25 == 0) {
      message(".", appendLF = FALSE)
    }
    elements <- x[[i]] %>% rvest::html_elements("a")
    hrefs <- elements %>% rvest::html_attr("href")
    ranks <- hrefs %>% gsub(".*/(.*?)/.*", "\\1", .)
    names <- elements %>%
      rvest::html_text() %>%
      gsub('"', "", ., fixed = TRUE)
    # no species, this must be until genus level
    hrefs <- hrefs[ranks != "species"]
    names <- names[ranks != "species"]
    ranks <- ranks[ranks != "species"]
    ranks[ranks == "domain"] <- "kingdom"

    df <- names %>%
      tibble() %>%
      t() %>%
      as_tibble() %>%
      setNames(ranks) %>%
      # no candidates please
      filter(genus %unlike% "^(Candidatus|\\[)")

    taxonomy_lpsn_missing <- taxonomy_lpsn_missing %>%
      bind_rows(df)
  }
  message(length(x), " entries incl. candidates (cleaned total: ", nrow(taxonomy_lpsn_missing), ")")
}

taxonomy_lpsn <- taxonomy_lpsn %>%
  left_join(taxonomy_lpsn_missing, by = "genus") %>%
  select(kingdom:family, everything()) %>%
  # remove entries like "[Bacteria, no family]" and "[Bacteria, no class]"
  mutate_all(function(x) ifelse(x %like_case% " no ", NA_character_, x))

taxonomy_lpsn.bak2 <- taxonomy_lpsn
# add family
pb <- progress::progress_bar$new(total = length(unique(taxonomy_lpsn$family)))
for (f in unique(taxonomy_lpsn$family)) {
  pb$tick()
  if (is.na(f)) next
  tax_info <- get_lpsn_and_author("Family", f)
  taxonomy_lpsn <- taxonomy_lpsn %>%
    bind_rows(tibble(
      kingdom = taxonomy_lpsn$kingdom[which(taxonomy_lpsn$family == f)[1]],
      phylum = taxonomy_lpsn$phylum[which(taxonomy_lpsn$family == f)[1]],
      class = taxonomy_lpsn$class[which(taxonomy_lpsn$family == f)[1]],
      order = taxonomy_lpsn$order[which(taxonomy_lpsn$family == f)[1]],
      family = f,
      rank = "family",
      status = "accepted",
      source = "LPSN",
      lpsn = unname(tax_info["lpsn"]),
      ref = unname(tax_info["ref"])
    ))
}
# add order
pb <- progress::progress_bar$new(total = length(unique(taxonomy_lpsn$order)))
for (o in unique(taxonomy_lpsn$order)) {
  pb$tick()
  if (is.na(o)) next
  tax_info <- get_lpsn_and_author("Order", o)
  taxonomy_lpsn <- taxonomy_lpsn %>%
    bind_rows(tibble(
      kingdom = taxonomy_lpsn$kingdom[which(taxonomy_lpsn$order == o)[1]],
      phylum = taxonomy_lpsn$phylum[which(taxonomy_lpsn$order == o)[1]],
      class = taxonomy_lpsn$class[which(taxonomy_lpsn$order == o)[1]],
      order = o,
      rank = "order",
      status = "accepted",
      source = "LPSN",
      lpsn = unname(tax_info["lpsn"]),
      ref = unname(tax_info["ref"])
    ))
}
# add class
pb <- progress::progress_bar$new(total = length(unique(taxonomy_lpsn$class)))
for (cc in unique(taxonomy_lpsn$class)) {
  pb$tick()
  if (is.na(cc)) next
  tax_info <- get_lpsn_and_author("Class", cc)
  taxonomy_lpsn <- taxonomy_lpsn %>%
    bind_rows(tibble(
      kingdom = taxonomy_lpsn$kingdom[which(taxonomy_lpsn$class == cc)[1]],
      phylum = taxonomy_lpsn$phylum[which(taxonomy_lpsn$class == cc)[1]],
      class = cc,
      rank = "class",
      status = "accepted",
      source = "LPSN",
      lpsn = unname(tax_info["lpsn"]),
      ref = unname(tax_info["ref"])
    ))
}
# add phylum
pb <- progress::progress_bar$new(total = length(unique(taxonomy_lpsn$phylum)))
for (p in unique(taxonomy_lpsn$phylum)) {
  pb$tick()
  if (is.na(p)) next
  tax_info <- get_lpsn_and_author("Phylum", p)
  taxonomy_lpsn <- taxonomy_lpsn %>%
    bind_rows(tibble(
      kingdom = taxonomy_lpsn$kingdom[which(taxonomy_lpsn$phylum == p)[1]],
      phylum = p,
      rank = "phylum",
      status = "accepted",
      source = "LPSN",
      lpsn = unname(tax_info["lpsn"]),
      ref = unname(tax_info["ref"])
    ))
}
# add kingdom
pb <- progress::progress_bar$new(total = length(unique(taxonomy_lpsn$kingdom)))
for (k in unique(taxonomy_lpsn$kingdom)) {
  pb$tick()
  if (is.na(k)) next
  tax_info <- get_lpsn_and_author("Domain", k)
  taxonomy_lpsn <- taxonomy_lpsn %>%
    bind_rows(tibble(
      kingdom = k,
      rank = "kingdom",
      status = "accepted",
      source = "LPSN",
      lpsn = unname(tax_info["lpsn"]),
      ref = unname(tax_info["ref"])
    ))
}

# integrity tests
sort(table(taxonomy_lpsn$rank))
sort(table(taxonomy_lpsn$status))


# Save intermediate results -----------------------------------------------

saveRDS(taxonomy_gbif, "data-raw/taxonomy_gbif.rds", version = 2)
saveRDS(taxonomy_lpsn, "data-raw/taxonomy_lpsn.rds", version = 2)


# Add full names ----------------------------------------------------------

taxonomy_gbif <- taxonomy_gbif %>%
  # clean NAs and add fullname
  mutate(across(kingdom:subspecies, function(x) ifelse(is.na(x), "", x)),
    fullname = trimws(case_when(
      rank == "family" ~ family,
      rank == "order" ~ order,
      rank == "class" ~ class,
      rank == "phylum" ~ phylum,
      rank == "kingdom" ~ kingdom,
      TRUE ~ paste(genus, species, subspecies)
    )), .before = 1
  ) %>%
  # keep only one GBIF taxon ID per full name
  arrange(fullname, gbif) %>%
  distinct(kingdom, fullname, .keep_all = TRUE)

taxonomy_lpsn <- taxonomy_lpsn %>%
  # clean NAs and add fullname
  mutate(across(kingdom:subspecies, function(x) ifelse(is.na(x), "", x)),
    fullname = trimws(case_when(
      rank == "family" ~ family,
      rank == "order" ~ order,
      rank == "class" ~ class,
      rank == "phylum" ~ phylum,
      rank == "kingdom" ~ kingdom,
      TRUE ~ paste(genus, species, subspecies)
    )), .before = 1
  ) %>%
  # keep only one LPSN record ID per full name
  arrange(fullname, lpsn) %>%
  distinct(kingdom, fullname, .keep_all = TRUE)

# set parent LPSN IDs, requires full name
taxonomy_lpsn$lpsn_parent[taxonomy_lpsn$rank == "phylum"] <- taxonomy_lpsn$lpsn[match(taxonomy_lpsn$kingdom[taxonomy_lpsn$rank == "phylum"], taxonomy_lpsn$fullname)]
taxonomy_lpsn$lpsn_parent[taxonomy_lpsn$rank == "class"] <- taxonomy_lpsn$lpsn[match(taxonomy_lpsn$phylum[taxonomy_lpsn$rank == "class"], taxonomy_lpsn$fullname)]
taxonomy_lpsn$lpsn_parent[taxonomy_lpsn$rank == "order"] <- taxonomy_lpsn$lpsn[match(taxonomy_lpsn$class[taxonomy_lpsn$rank == "order"], taxonomy_lpsn$fullname)]
taxonomy_lpsn$lpsn_parent[taxonomy_lpsn$rank == "family"] <- taxonomy_lpsn$lpsn[match(taxonomy_lpsn$order[taxonomy_lpsn$rank == "family"], taxonomy_lpsn$fullname)]
taxonomy_lpsn$lpsn_parent[taxonomy_lpsn$rank == "genus"] <- taxonomy_lpsn$lpsn[match(taxonomy_lpsn$family[taxonomy_lpsn$rank == "genus"], taxonomy_lpsn$fullname)]
taxonomy_lpsn$lpsn_parent[taxonomy_lpsn$rank == "species"] <- taxonomy_lpsn$lpsn[match(taxonomy_lpsn$genus[taxonomy_lpsn$rank == "species"], taxonomy_lpsn$fullname)]
taxonomy_lpsn$lpsn_parent[taxonomy_lpsn$rank == "subspecies"] <- taxonomy_lpsn$lpsn[match(paste(taxonomy_lpsn$genus[taxonomy_lpsn$rank == "subspecies"], taxonomy_lpsn$species[taxonomy_lpsn$rank == "subspecies"]), taxonomy_lpsn$fullname)]


# Combine the datasets ----------------------------------------------------

# basis must be LPSN as it's most recent
taxonomy <- taxonomy_lpsn %>%
  # join GBIF identifiers to them
  left_join(taxonomy_gbif %>% select(kingdom, fullname, starts_with("gbif")),
    by = c("kingdom", "fullname")
  )

# for everything else, add the GBIF data
taxonomy <- taxonomy %>%
  bind_rows(taxonomy_gbif %>%
    filter(!paste(kingdom, fullname) %in% paste(taxonomy$kingdom, taxonomy$fullname))) %>%
  arrange(fullname) %>%
  filter(fullname != "")

# fix rank
taxonomy <- taxonomy %>%
  mutate(rank = case_when(
    subspecies != "" ~ "subspecies",
    species != "" ~ "species",
    genus != "" ~ "genus",
    family != "" ~ "family",
    order != "" ~ "order",
    class != "" ~ "class",
    phylum != "" ~ "phylum",
    kingdom != "" ~ "kingdom",
    TRUE ~ NA_character_
  ))

table(taxonomy$rank, useNA = "always")

# get the latest upper taxonomy from LPSN to update the GBIF data
# (e.g., phylum above class "Bacilli" was still "Firmicutes", should be "Bacillota" in 2022)
for (k in unique(taxonomy$kingdom[taxonomy$kingdom != ""])) {
  message("Fixing GBIF taxonomy for kingdom ", k, "...")
  for (g in unique(taxonomy$genus[taxonomy$genus != "" & taxonomy$kingdom == k & taxonomy$source == "LPSN"])) {
    taxonomy$family[which(taxonomy$genus == g & taxonomy$kingdom == k)] <- taxonomy$family[which(taxonomy$genus == g & taxonomy$kingdom == k & taxonomy$source == "LPSN")][1]
  }
  for (f in unique(taxonomy$family[taxonomy$family != "" & taxonomy$kingdom == k & taxonomy$source == "LPSN"])) {
    taxonomy$order[which(taxonomy$family == f & taxonomy$kingdom == k)] <- taxonomy$order[which(taxonomy$family == f & taxonomy$kingdom == k & taxonomy$source == "LPSN")][1]
  }
  for (o in unique(taxonomy$order[taxonomy$order != "" & taxonomy$kingdom == k & taxonomy$source == "LPSN"])) {
    taxonomy$class[which(taxonomy$order == o & taxonomy$kingdom == k)] <- taxonomy$class[which(taxonomy$order == o & taxonomy$kingdom == k & taxonomy$source == "LPSN")][1]
  }
  for (cc in unique(taxonomy$class[taxonomy$class != "" & taxonomy$kingdom == k & taxonomy$source == "LPSN"])) {
    taxonomy$phylum[which(taxonomy$class == cc & taxonomy$kingdom == k)] <- taxonomy$phylum[which(taxonomy$class == cc & taxonomy$kingdom == k & taxonomy$source == "LPSN")][1]
  }
}


# Add missing taxonomic entries -------------------------------------------

# this part will make sure that the whole taxonomy of every included species exists, so no missing genera, classes, etc.

current_gbif <- taxonomy_gbif.bak %>%
  filter(is.na(acceptedNameUsageID)) %>%
  mutate(
    taxonID = as.character(taxonID),
    parentNameUsageID = as.character(parentNameUsageID)
  )

# add missing kingdoms
taxonomy <- taxonomy %>%
  bind_rows(
    taxonomy %>%
      filter(kingdom != "") %>%
      distinct(kingdom) %>%
      mutate(
        fullname = kingdom,
        rank = "kingdom",
        status = "accepted",
        source = "manually added"
      ) %>%
      filter(!paste(kingdom, rank) %in% paste(taxonomy$kingdom, taxonomy$rank)) %>%
      left_join(current_gbif %>%
        select(kingdom, rank = taxonRank, ref = scientificNameAuthorship, gbif = taxonID, gbif_parent = parentNameUsageID),
      by = c("kingdom", "rank")
      ) %>%
      mutate(source = ifelse(!is.na(gbif), "GBIF", source))
  )

# 2 = phylum ... 6 = genus
for (i in 2:6) {
  i_name <- colnames(taxonomy)[i + 1]
  message("Adding missing: ", i_name, "... ", appendLF = FALSE)
  to_add <- taxonomy %>%
    filter(.[[i + 1]] != "") %>%
    distinct(kingdom, .[[i + 1]], .keep_all = TRUE) %>%
    select(kingdom:(i + 1)) %>%
    mutate(
      fullname = .[[ncol(.)]],
      rank = i_name,
      status = "accepted",
      source = "manually added"
    ) %>%
    filter(!paste(kingdom, .[[ncol(.) - 4]], rank) %in% paste(taxonomy$kingdom, taxonomy[[i + 1]], taxonomy$rank)) %>%
    # get GBIF identifier where available
    left_join(current_gbif %>%
      select(kingdom, all_of(i_name), rank = taxonRank, ref = scientificNameAuthorship, gbif = taxonID, gbif_parent = parentNameUsageID),
    by = c("kingdom", "rank", i_name)
    ) %>%
    mutate(source = ifelse(!is.na(gbif), "GBIF", source))
  message("n = ", nrow(to_add))
  taxonomy <- taxonomy %>%
    bind_rows(to_add)
}

# species (requires combination with genus)
taxonomy <- taxonomy %>%
  bind_rows(taxonomy %>%
    filter(species != "") %>%
    distinct(kingdom, genus, species, .keep_all = TRUE) %>%
    select(kingdom:species) %>%
    mutate(
      fullname = paste(genus, species),
      rank = "species",
      status = "accepted",
      source = "manually added"
    ) %>%
    filter(!paste(kingdom, genus, species, rank) %in% paste(taxonomy$kingdom, taxonomy$genus, taxonomy$species, taxonomy$rank)) %>%
    # get GBIF identifier where available
    left_join(current_gbif %>%
      select(kingdom, genus, species = specificEpithet, rank = taxonRank, ref = scientificNameAuthorship, gbif = taxonID, gbif_parent = parentNameUsageID),
    by = c("kingdom", "rank", "genus", "species")
    ) %>%
    mutate(source = ifelse(!is.na(gbif), "GBIF", source)))


# remove NAs from taxonomy again, and keep unique full names
taxonomy <- taxonomy %>%
  mutate(across(kingdom:subspecies, function(x) ifelse(is.na(x), "", x))) %>%
  distinct(kingdom, fullname, .keep_all = TRUE) %>%
  filter(kingdom != "")

# Save intermediate results -----------------------------------------------

saveRDS(taxonomy, "data-raw/taxonomy.rds")


# Get previously manually added entries -----------------------------------

manually_added <- AMR::microorganisms %>%
  filter(source == "manually added", !fullname %in% taxonomy$fullname) %>%
  select(fullname:subspecies, ref, source, rank)

# get latest taxonomy for those entries
for (g in unique(manually_added$genus[manually_added$genus != "" & manually_added$genus %in% taxonomy$genus])) {
  manually_added$family[which(manually_added$genus == g)] <- taxonomy$family[which(taxonomy$genus == g & is.na(taxonomy$lpsn))][1]
}
for (f in unique(manually_added$family[manually_added$family != "" & manually_added$family %in% taxonomy$family])) {
  manually_added$order[which(manually_added$family == f)] <- taxonomy$order[which(taxonomy$family == f & is.na(taxonomy$lpsn))][1]
}
for (o in unique(manually_added$order[manually_added$order != "" & manually_added$order %in% taxonomy$order])) {
  manually_added$class[which(manually_added$order == o)] <- taxonomy$class[which(taxonomy$order == o & is.na(taxonomy$lpsn))][1]
}
for (cc in unique(manually_added$class[manually_added$class != "" & manually_added$class %in% taxonomy$class])) {
  manually_added$phylum[which(manually_added$class == cc)] <- taxonomy$phylum[which(taxonomy$class == cc & is.na(taxonomy$lpsn))][1]
}

# previously required:
taxonomy$ref[which(taxonomy$genus == "Streptococcus" & taxonomy$species %like% "group")] <- "Lancefield, 1933"

taxonomy <- taxonomy %>%
  bind_rows(manually_added %>%
    mutate(
      status = "accepted",
      rank = ifelse(fullname %like% "unknown", "(unknown rank)", rank)
    )) %>%
  arrange(fullname)

table(taxonomy$rank, useNA = "always")


# Clean scientific reference ----------------------------------------------

taxonomy <- taxonomy %>%
  mutate(ref = get_author_year(ref))


# Add prevalence ----------------------------------------------------------

# update prevalence based on taxonomy (our own JSS paper: Berends et al., 2022)
taxonomy <- taxonomy %>%
  mutate(prevalence = case_when(
    class == "Gammaproteobacteria" |
      genus %in% c("Enterococcus", "Staphylococcus", "Streptococcus")
    ~ 1,
    kingdom %in% c("Archaea", "Bacteria", "Chromista", "Fungi") &
      (phylum %in% c(
        "Proteobacteria",
        "Firmicutes",
        "Bacillota", # this one is new - this was renamed from Firmicutes by Gibbons et al., 2021
        "Actinobacteria",
        "Sarcomastigophora"
      ) |
        genus %in% MO_PREVALENT_GENERA)
    ~ 2,
    TRUE ~ 3
  ))
table(taxonomy$prevalence, useNA = "always")


# Add microbial IDs -------------------------------------------------------

# MO codes in the AMR package have the form KINGDOM_GENUS_SPECIES_SUBSPECIES where all are abbreviated.

# Kingdom is abbreviated with 1 character, with exceptions for Animalia and Plantae
mo_kingdom <- taxonomy %>%
  filter(rank == "kingdom") %>%
  select(kingdom) %>%
  mutate(mo_kingdom = case_when(
    kingdom == "Animalia" ~ "AN",
    kingdom == "Archaea" ~ "A",
    kingdom == "Bacteria" ~ "B",
    kingdom == "Chromista" ~ "C",
    kingdom == "Fungi" ~ "F",
    kingdom == "Plantae" ~ "PL",
    kingdom == "Protozoa" ~ "P",
    TRUE ~ ""
  ))
# phylum until family are abbreviated with 8 characters and prefixed with their rank

# Phylum - keep old and fill up for new ones
mo_phylum <- taxonomy %>%
  filter(rank == "phylum") %>%
  distinct(kingdom, phylum) %>%
  left_join(AMR::microorganisms %>%
    filter(rank == "phylum") %>%
    transmute(kingdom,
      phylum = fullname,
      mo_old = gsub("[A-Z]{1,2}_", "", as.character(mo))
    ),
  by = c("kingdom", "phylum")
  ) %>%
  group_by(kingdom) %>%
  mutate(
    mo_phylum8 = abbreviate_mo(phylum, minlength = 8, prefix = "[PHL]_"),
    mo_phylum9 = abbreviate_mo(phylum, minlength = 9, prefix = "[PHL]_"),
    mo_phylum = ifelse(!is.na(mo_old), mo_old, mo_phylum8),
    mo_duplicated = duplicated(mo_phylum),
    mo_phylum = ifelse(mo_duplicated, mo_phylum9, mo_phylum),
    mo_duplicated = duplicated(mo_phylum)
  ) %>%
  ungroup()
if (any(mo_phylum$mo_duplicated, na.rm = TRUE)) stop("Duplicate MO codes for phylum!")
mo_phylum <- mo_phylum %>%
  select(kingdom, phylum, mo_phylum)

# Class - keep old and fill up for new ones
mo_class <- taxonomy %>%
  filter(rank == "class") %>%
  distinct(kingdom, class) %>%
  left_join(AMR::microorganisms %>%
    filter(rank == "class") %>%
    transmute(kingdom,
      class = fullname,
      mo_old = gsub("[A-Z]{1,2}_", "", as.character(mo))
    ),
  by = c("kingdom", "class")
  ) %>%
  group_by(kingdom) %>%
  mutate(
    mo_class8 = abbreviate_mo(class, minlength = 8, prefix = "[CLS]_"),
    mo_class9 = abbreviate_mo(class, minlength = 9, prefix = "[CLS]_"),
    mo_class = ifelse(!is.na(mo_old), mo_old, mo_class8),
    mo_duplicated = duplicated(mo_class),
    mo_class = ifelse(mo_duplicated, mo_class9, mo_class),
    mo_duplicated = duplicated(mo_class)
  ) %>%
  ungroup()
if (any(mo_class$mo_duplicated, na.rm = TRUE)) stop("Duplicate MO codes for class!")
mo_class <- mo_class %>%
  select(kingdom, class, mo_class)

# Order - keep old and fill up for new ones
mo_order <- taxonomy %>%
  filter(rank == "order") %>%
  distinct(kingdom, order) %>%
  left_join(AMR::microorganisms %>%
    filter(rank == "order") %>%
    transmute(kingdom,
      order = fullname,
      mo_old = gsub("[A-Z]{1,2}_", "", as.character(mo))
    ),
  by = c("kingdom", "order")
  ) %>%
  group_by(kingdom) %>%
  mutate(
    mo_order8 = abbreviate_mo(order, minlength = 8, prefix = "[ORD]_"),
    mo_order9 = abbreviate_mo(order, minlength = 9, prefix = "[ORD]_"),
    mo_order = ifelse(!is.na(mo_old), mo_old, mo_order8),
    mo_duplicated = duplicated(mo_order),
    mo_order = ifelse(mo_duplicated, mo_order9, mo_order),
    mo_duplicated = duplicated(mo_order)
  ) %>%
  ungroup()
if (any(mo_order$mo_duplicated, na.rm = TRUE)) stop("Duplicate MO codes for order!")
mo_order <- mo_order %>%
  select(kingdom, order, mo_order)

# Family - keep old and fill up for new ones
mo_family <- taxonomy %>%
  filter(rank == "family") %>%
  distinct(kingdom, family) %>%
  left_join(AMR::microorganisms %>%
    filter(rank == "family") %>%
    transmute(kingdom,
      family = fullname,
      mo_old = gsub("[A-Z]{1,2}_", "", as.character(mo))
    ),
  by = c("kingdom", "family")
  ) %>%
  group_by(kingdom) %>%
  mutate(
    mo_family8 = abbreviate_mo(family, minlength = 8, prefix = "[FAM]_"),
    mo_family9 = abbreviate_mo(family, minlength = 9, prefix = "[FAM]_"),
    mo_family = ifelse(!is.na(mo_old), mo_old, mo_family8),
    mo_duplicated = duplicated(mo_family),
    mo_family = ifelse(mo_duplicated, mo_family9, mo_family),
    mo_duplicated = duplicated(mo_family)
  ) %>%
  ungroup()
if (any(mo_family$mo_duplicated, na.rm = TRUE)) stop("Duplicate MO codes for family!")
mo_family <- mo_family %>%
  select(kingdom, family, mo_family)

# construct code part for genus - keep old code where available and generate new ones where needed
mo_genus <- taxonomy %>%
  filter(rank == "genus") %>%
  distinct(kingdom, genus) %>%
  # get available old MO codes
  left_join(AMR::microorganisms %>%
    filter(rank == "genus") %>%
    transmute(mo_genus_old = gsub("^[A-Z]+_", "", as.character(mo)), kingdom, genus) %>%
    distinct(kingdom, genus, .keep_all = TRUE),
  by = c("kingdom", "genus")
  ) %>%
  distinct(kingdom, genus, .keep_all = TRUE) %>%
  # since kingdom is part of the code, genus abbreviations may be duplicated between kingdoms
  group_by(kingdom) %>%
  # generate new MO codes for genus and set the right one
  mutate(
    mo_genus_new5 = abbreviate_mo(genus, 5),
    mo_genus_new5b = paste0(abbreviate_mo(genus, 5), 1),
    mo_genus_new6 = abbreviate_mo(genus, 6),
    mo_genus_new7 = abbreviate_mo(genus, 7),
    mo_genus_new8 = abbreviate_mo(genus, 8),
    mo_genus_new = case_when(
      !is.na(mo_genus_old) ~ mo_genus_old,
      !mo_genus_new5 %in% mo_genus_old ~ mo_genus_new5,
      !mo_genus_new6 %in% mo_genus_old ~ mo_genus_new6,
      !mo_genus_new7 %in% mo_genus_old ~ mo_genus_new7,
      !mo_genus_new8 %in% mo_genus_old ~ mo_genus_new8,
      !mo_genus_new5b %in% mo_genus_old ~ mo_genus_new5b,
      TRUE ~ mo_genus_old
    ),
    mo_duplicated = duplicated(mo_genus_new),
    mo_genus_new = case_when(
      !mo_duplicated ~ mo_genus_new,
      mo_duplicated & mo_genus_new == mo_genus_new5 ~ mo_genus_new6,
      mo_duplicated & mo_genus_new == mo_genus_new6 ~ mo_genus_new7,
      mo_duplicated & mo_genus_new == mo_genus_new7 ~ mo_genus_new8,
      mo_duplicated & mo_genus_new == mo_genus_new8 ~ mo_genus_new5b,
      TRUE ~ NA_character_
    ),
    mo_duplicated = duplicated(mo_genus_new)
  ) %>%
  ungroup()
if (any(mo_genus$mo_duplicated, na.rm = TRUE) | anyNA(mo_genus$mo_genus_new)) stop("Duplicate MO codes for genus!")
# no duplicates *within kingdoms*, so keep the right columns for left joining later
mo_genus <- mo_genus %>%
  select(kingdom, genus, mo_genus = mo_genus_new)

# same for species - keep old where available and create new per kingdom-genus where needed:
mo_species <- taxonomy %>%
  filter(rank == "species") %>%
  distinct(kingdom, genus, species) %>%
  left_join(microorganisms %>%
    filter(rank == "species") %>%
    transmute(mo_species_old = gsub("^[A-Z]+_[A-Z]+_", "", as.character(mo)), kingdom, genus, species) %>%
    filter(mo_species_old %unlike% "-") %>%
    distinct(kingdom, genus, species, .keep_all = TRUE),
  by = c("kingdom", "genus", "species")
  ) %>%
  distinct(kingdom, genus, species, .keep_all = TRUE) %>%
  group_by(kingdom, genus) %>%
  mutate(
    mo_species_new4 = abbreviate_mo(species, 4, hyphen_as_space = TRUE),
    mo_species_new5 = abbreviate_mo(species, 5, hyphen_as_space = TRUE),
    mo_species_new5b = paste0(abbreviate_mo(species, 5, hyphen_as_space = TRUE), 1),
    mo_species_new6 = abbreviate_mo(species, 6, hyphen_as_space = TRUE),
    mo_species_new7 = abbreviate_mo(species, 7, hyphen_as_space = TRUE),
    mo_species_new8 = abbreviate_mo(species, 8, hyphen_as_space = TRUE),
    mo_species_new = case_when(
      !is.na(mo_species_old) ~ mo_species_old,
      !mo_species_new4 %in% mo_species_old ~ mo_species_new4,
      !mo_species_new5 %in% mo_species_old ~ mo_species_new5,
      !mo_species_new6 %in% mo_species_old ~ mo_species_new6,
      !mo_species_new7 %in% mo_species_old ~ mo_species_new7,
      !mo_species_new8 %in% mo_species_old ~ mo_species_new8,
      !mo_species_new5b %in% mo_species_old ~ mo_species_new5b,
      TRUE ~ mo_species_old
    ),
    mo_duplicated = duplicated(mo_species_new),
    mo_species_new = case_when(
      !mo_duplicated ~ mo_species_new,
      mo_duplicated & mo_species_new == mo_species_new4 ~ mo_species_new5,
      mo_duplicated & mo_species_new == mo_species_new5 ~ mo_species_new6,
      mo_duplicated & mo_species_new == mo_species_new6 ~ mo_species_new7,
      mo_duplicated & mo_species_new == mo_species_new7 ~ mo_species_new8,
      mo_duplicated & mo_species_new == mo_species_new8 ~ mo_species_new5b,
      TRUE ~ NA_character_
    ),
    mo_duplicated = duplicated(mo_species_new)
  ) %>%
  ungroup()
if (any(mo_species$mo_duplicated, na.rm = TRUE) | anyNA(mo_species$mo_species_new)) stop("Duplicate MO codes for species!")
# no duplicates *within kingdoms*, so keep the right columns for left joining later
mo_species <- mo_species %>%
  select(kingdom, genus, species, mo_species = mo_species_new)

# same for subspecies - keep old where available and create new per kingdom-genus-species where needed:
mo_subspecies <- taxonomy %>%
  filter(rank == "subspecies") %>%
  distinct(kingdom, genus, species, subspecies) %>%
  left_join(microorganisms %>%
    filter(rank %in% c("subspecies", "subsp.", "infraspecies")) %>%
    transmute(mo_subspecies_old = gsub("^[A-Z]+_[A-Z]+_[A-Z]+_", "", as.character(mo)), kingdom, genus, species, subspecies) %>%
    filter(mo_subspecies_old %unlike% "-") %>%
    distinct(kingdom, genus, species, subspecies, .keep_all = TRUE),
  by = c("kingdom", "genus", "species", "subspecies")
  ) %>%
  distinct(kingdom, genus, species, subspecies, .keep_all = TRUE) %>%
  group_by(kingdom, genus, species) %>%
  mutate(
    mo_subspecies_new4 = abbreviate_mo(subspecies, 4, hyphen_as_space = TRUE),
    mo_subspecies_new5 = abbreviate_mo(subspecies, 5, hyphen_as_space = TRUE),
    mo_subspecies_new5b = paste0(abbreviate_mo(subspecies, 5, hyphen_as_space = TRUE), 1),
    mo_subspecies_new6 = abbreviate_mo(subspecies, 6, hyphen_as_space = TRUE),
    mo_subspecies_new7 = abbreviate_mo(subspecies, 7, hyphen_as_space = TRUE),
    mo_subspecies_new8 = abbreviate_mo(subspecies, 8, hyphen_as_space = TRUE),
    mo_subspecies_new = case_when(
      !is.na(mo_subspecies_old) ~ mo_subspecies_old,
      !mo_subspecies_new4 %in% mo_subspecies_old ~ mo_subspecies_new4,
      !mo_subspecies_new5 %in% mo_subspecies_old ~ mo_subspecies_new5,
      !mo_subspecies_new6 %in% mo_subspecies_old ~ mo_subspecies_new6,
      !mo_subspecies_new7 %in% mo_subspecies_old ~ mo_subspecies_new7,
      !mo_subspecies_new8 %in% mo_subspecies_old ~ mo_subspecies_new8,
      !mo_subspecies_new5b %in% mo_subspecies_old ~ mo_subspecies_new5b,
      TRUE ~ mo_subspecies_old
    ),
    mo_duplicated = duplicated(mo_subspecies_new),
    mo_subspecies_new = case_when(
      !mo_duplicated ~ mo_subspecies_new,
      mo_duplicated & mo_subspecies_new == mo_subspecies_new4 ~ mo_subspecies_new5,
      mo_duplicated & mo_subspecies_new == mo_subspecies_new5 ~ mo_subspecies_new6,
      mo_duplicated & mo_subspecies_new == mo_subspecies_new6 ~ mo_subspecies_new7,
      mo_duplicated & mo_subspecies_new == mo_subspecies_new7 ~ mo_subspecies_new8,
      mo_duplicated & mo_subspecies_new == mo_subspecies_new8 ~ mo_subspecies_new5b,
      TRUE ~ NA_character_
    ),
    mo_duplicated = duplicated(mo_subspecies_new)
  ) %>%
  ungroup()
if (any(mo_subspecies$mo_duplicated, na.rm = TRUE) | anyNA(mo_subspecies$mo_subspecies_new)) stop("Duplicate MO codes for subspecies!")
# no duplicates *within kingdoms*, so keep the right columns for left joining later
mo_subspecies <- mo_subspecies %>%
  select(kingdom, genus, species, subspecies, mo_subspecies = mo_subspecies_new)

# unknowns - manually added
mo_unknown <- AMR::microorganisms %>%
  filter(fullname %like% "unknown") %>%
  transmute(fullname, mo_unknown = as.character(mo))

# apply the new codes!
taxonomy <- taxonomy %>%
  left_join(mo_kingdom, by = "kingdom") %>%
  left_join(mo_phylum, by = c("kingdom", "phylum")) %>%
  left_join(mo_class, by = c("kingdom", "class")) %>%
  left_join(mo_order, by = c("kingdom", "order")) %>%
  left_join(mo_family, by = c("kingdom", "family")) %>%
  left_join(mo_genus, by = c("kingdom", "genus")) %>%
  left_join(mo_species, by = c("kingdom", "genus", "species")) %>%
  left_join(mo_subspecies, by = c("kingdom", "genus", "species", "subspecies")) %>%
  left_join(mo_unknown, by = "fullname") %>%
  mutate(across(starts_with("mo_"), function(x) ifelse(is.na(x), "", x))) %>%
  mutate(
    mo = case_when(
      fullname %like% "unknown" ~ mo_unknown,
      # add special cases for taxons higher than genus
      rank == "kingdom" ~ paste(mo_kingdom, "[KNG]", toupper(kingdom), sep = "_"),
      rank == "phylum" ~ paste(mo_kingdom, mo_phylum, sep = "_"),
      rank == "class" ~ paste(mo_kingdom, mo_class, sep = "_"),
      rank == "order" ~ paste(mo_kingdom, mo_order, sep = "_"),
      rank == "family" ~ paste(mo_kingdom, mo_family, sep = "_"),
      TRUE ~ paste(mo_kingdom, mo_genus, mo_species, mo_subspecies, sep = "_")
    ),
    mo = trimws(gsub("_+$", "", mo)),
    .before = 1
  ) %>%
  select(!starts_with("mo_")) %>%
  arrange(fullname) %>%
  distinct(fullname, .keep_all = TRUE)


# Remove unwanted taxonomic entries from Protoza/Fungi --------------------

# this must be done after the microbial ID generation, since it will otherwise generate a lot of different IDs
taxonomy <- taxonomy %>% 
  filter(
    # Protozoa:
    !(phylum %in% c("Choanozoa", "Mycetozoa") & prevalence == 3),
    # Fungi:
    !(phylum %in% c("Ascomycota", "Zygomycota", "Basidiomycota") & prevalence == 3),
    !(genus %in% c("Phoma", "Leptosphaeria") & rank %in% c("species", "subspecies")), # only genus of this rare fungus, with resp. 1300 and 800 species
    # (leave Alternaria in there, part of human mycobiome and opportunistic pathogen)
    # Animalia:
    !genus %in% c("Lucilia", "Lumbricus"),
    !(genus %in% c("Aedes", "Anopheles") & rank %in% c("species", "subspecies")), # only genus of the many hundreds of mosquitoes species
    kingdom != "Plantae") # this kingdom only contained Curvularia and Hymenolepis, which have coincidental twin names with Fungi

message("\nCongratulations! The new taxonomic table will contain ", format(nrow(taxonomy), big.mark = ","), " rows.\n")


# Add SNOMED CT -----------------------------------------------------------

# we will use Public Health Information Network Vocabulary Access and Distribution System (PHIN VADS)
# as a source, which copies directly from the latest US SNOMED CT version
# - go to https://phinvads.cdc.gov/vads/ViewValueSet.action?oid=2.16.840.1.114222.4.11.1009
# - check that current online version is higher than SNOMED_VERSION$current_version
# - if so, click on 'Download Value Set', choose 'TXT'
snomed <- vroom("data-raw/SNOMED_PHVS_Microorganism_CDC_V12.txt", skip = 3) %>%
  select(1:2) %>%
  setNames(c("snomed", "mo")) %>%
  mutate(snomed = as.character(snomed))

# try to get name of MO
snomed <- snomed %>%
  mutate(mo = gsub("ss. ", "", mo, fixed = TRUE)) %>%
  mutate(fullname = case_when(
    mo %like_case% "[A-Z][a-z]+ [a-z]+ [a-z]{4,} " ~ gsub("(^|.*)([A-Z][a-z]+ [a-z]+ [a-z]{4,}) .*", "\\2", mo),
    mo %like_case% "[A-Z][a-z]+ [a-z]{4,} " ~ gsub("(^|.*)([A-Z][a-z]+ [a-z]{4,}) .*", "\\2", mo),
    mo %like_case% "[A-Z][a-z]+" ~ gsub("(^|.*)([A-Z][a-z]+) .*", "\\2", mo),
    TRUE ~ NA_character_
  )) %>%
  filter(fullname %in% taxonomy$fullname)

message(nrow(snomed), " SNOMED codes will be added to ", n_distinct(snomed$fullname), " microorganisms")

snomed <- snomed %>%
  group_by(fullname) %>%
  summarise(snomed = list(snomed))

taxonomy <- taxonomy %>%
  left_join(snomed, by = "fullname")


# Clean data set ----------------------------------------------------------

# format to tibble and check again for invalid characters
taxonomy <- taxonomy %>%
  select(mo, fullname, status, kingdom:subspecies, rank, ref, source, starts_with("lpsn"), starts_with("gbif"), prevalence, snomed) %>%
  df_remove_nonASCII()

# set class <mo>
class(taxonomy$mo) <- c("mo", "character")

# Moraxella catarrhalis was named Branhamella catarrhalis (Catlin, 1970), but this is unaccepted in clinical microbiology
# we keep them both
taxonomy$status[which(taxonomy$fullname == "Moraxella catarrhalis")]
taxonomy$lpsn_renamed_to[which(taxonomy$fullname == "Moraxella catarrhalis")]
taxonomy$status[which(taxonomy$fullname == "Moraxella catarrhalis")] <- "accepted"
taxonomy$lpsn_renamed_to[which(taxonomy$fullname == "Moraxella catarrhalis")] <- NA_character_

taxonomy <- taxonomy %>% 
  AMR:::dataset_UTF8_to_ASCII()


# Save to package ---------------------------------------------------------

microorganisms <- taxonomy
usethis::use_data(microorganisms, overwrite = TRUE, version = 2, compress = "xz")
rm(microorganisms)

# DON'T FORGET TO UPDATE R/globals.R!


# Test updates ------------------------------------------------------------

# and check: these codes should not be missing (will otherwise throw a unit test error):
AMR::microorganisms.codes %>% filter(!mo %in% taxonomy$mo)
AMR::rsi_translation %>% filter(!mo %in% taxonomy$mo)
AMR::example_isolates %>% filter(!mo %in% taxonomy$mo)
AMR::intrinsic_resistant %>% filter(!mo %in% taxonomy$mo)

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


# run the unit tests
Sys.setenv(NOT_CRAN = "true")
testthat::test_file("tests/testthat/test-data.R")
testthat::test_file("tests/testthat/test-mo.R")
testthat::test_file("tests/testthat/test-mo_property.R")
