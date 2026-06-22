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
# how to conduct AMR data analysis: https://amr-for-r.org              #
# ==================================================================== #


# TODO Giardia zat in protozoa en animalia, dus sorteren op prioriteit kingdom en dan distinct op fullname



# ! THIS SCRIPT REQUIRES AT LEAST 16 GB RAM !
# (at least 12 GB will be used by the R session for the size of the files)

# GBIF:
# 1. Go to https://doi.org/10.48580/dgxjw and find the download link
#    under "Endpoints" and unpack "Taxon.tsv" from it (~2.5 GB)
#    ALSO BE SURE to get the date of release and update R/aa_globals.R later!
# LPSN:
# 2. Go to https://lpsn.dsmz.de/downloads (register first) and download the latest
#    CSV file (~12,5 MB) and rename to "taxonomy.csv"
#    ALSO BE SURE to get the date of release and update R/aa_globals.R later!
# MycoBank:
# 3. Go to https://www.mycobank.org/ and find the download link of all entries
#    (last time https://www.mycobank.org/Images/MBList.zip) and unpack
#    "MBList.xlsx" from it (~120 MB)
#    ALSO BE SURE to get the date of release and update R/aa_globals.R later!
# Bartlett:
# 4. For data about human pathogens, we use Bartlett et al. (2022),
#    https://doi.org/10.1099/mic.0.001269. Their latest supplementary material
#    can be found here: https://github.com/padpadpadpad/bartlett_et_al_2022_human_pathogens.
#    Download their latest xlsx file in the `data` folder and save it to our
#    `data-raw` folder.
# 5. Go to BacDive data base for the oxygen tolerance and cell shape.
#    Go to https://bacdive.dsmz.de/advsearch, filter 'Oxygen tolerance' or
#    'Cell shape' on "*" and click Submit and click on the 'Download tabel as CSV' button
# 6. Set these locations to the paths where the files are:
folder_location <- "~/Downloads/"
file_gbif <- paste0(folder_location, "/xr_latest_dwca/Taxon.tsv")
file_lpsn <- paste0(folder_location, "taxonomy.csv")
file_mycobank <- paste0(folder_location, "MBList-2.xlsx")

file_bartlett <- "data-raw/bartlett_et_al_2022_human_pathogens.xlsx"

file_oxygen_tolerance <- paste0(folder_location, "oxygen.csv")
file_cell_shape <- paste0(folder_location, "shape.csv")

# 4. Run the rest of this script line by line and check everything :)

if (!file.exists(file_gbif)) {
  stop("GBIF file not found")
}
if (!file.exists(file_lpsn)) {
  stop("LPSN file not found")
}
if (!file.exists(file_mycobank)) {
  stop("MB file not found")
}
if (!file.exists(file_bartlett)) {
  stop("Bartlett et al. Excel file not found")
}
if (!file.exists(file_oxygen_tolerance)) {
  stop("BacDive oxygen tolerance CSV file not found")
}
if (!file.exists(file_cell_shape)) {
  stop("BacDive cell shape CSV file not found")
}

library(dplyr)
library(tidyr)
library(vroom) # to import files
library(rvest) # to scrape LPSN website
library(progress) # to show progress bars
library(readxl) # to read the MycoBank and Bartlett Excel files
devtools::load_all(".") # to load the AMR package


# Helper functions --------------------------------------------------------------------------------

get_author_year <- function(ref) {
  # Only keep first author, e.g. transform 'Smith, Jones, 2011' to 'Smith et al., 2011'

  # normalise Unicode to NFC (precomposed) so iconv can transliterate
  authors2 <- stringi::stri_trans_nfc(ref)
  # fix known encoding errors in source data
  authors2 <- gsub("\u0092", "'", authors2)              # Windows-1252 right quote
  authors2 <- gsub("\u0091", "'", authors2)              # Windows-1252 left quote
  authors2 <- gsub("\u04AB", "\u00E7", authors2)         # Cyrillic cedilla to Latin cedilla
  authors2 <- gsub("[\u0400-\u04FF]", "", authors2)      # remove remaining Cyrillic characters
  authors2 <- iconv(authors2, from = "UTF-8", to = "ASCII//TRANSLIT")
  
  authors2 <- gsub(" ?\\(Approved Lists [0-9]+\\) ?", " ", authors2)
  authors2 <- gsub(" +", " ", authors2)
  authors2 <- trimws(authors2)
  # remove leading and trailing brackets
  authors2 <- trimws(gsub("^[(](.*)[)]$", "\\1", authors2))
  # only take part after brackets if there's a name
  authors2 <- if_else(
    grepl(".*[)] [a-zA-Z]+.*", authors2),
    gsub(".*[)] (.*)", "\\1", authors2),
    authors2
  )
  # replace parentheses with emend. to get the latest authors
  # authors2 <- gsub("(", " emend. ", authors2, fixed = TRUE)
  # authors2 <- gsub(")", "", authors2, fixed = TRUE)
  
  # remove any remaining parentheses
  authors2 <- gsub("[()]", "", authors2)
  authors2 <- gsub(" +", " ", authors2)
  authors2 <- trimws(authors2)
  
  # strip emend. and everything after it to retain the combination authority
  authors2 <- gsub(" ?emend[.]?.*", "", authors2)
  
  # get year from last 4 digits
  lastyear <- as.integer(gsub(".*([0-9]{4})$", "\\1", authors2))
  # can never be later than now
  lastyear <- if_else(
    lastyear > as.integer(format(Sys.Date(), "%Y")),
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
  authors <- gsub("[A-Z-][a-z-]?[.]", "", authors, ignore.case = FALSE)
  # remove trailing and leading spaces
  authors <- trimws(authors)
  # strip emend. and everything after it to retain the combination authority
  authors <- gsub(" ?emend[.]?.*", "", authors)
  # only keep first author and replace all others by 'et al'
  authors <- gsub("(,| and| et| &| ex| emend\\.?) .*", " et al.", authors)
  # et al. always with ending dot
  authors <- gsub(" et al\\.?", " et al.", authors)
  authors <- gsub(" ?,$", "", authors)
  # don't start with 'sensu' or 'ehrenb'
  authors <- gsub(
    "^(sensu|Ehrenb.?|corrig.?) ",
    "",
    authors,
    ignore.case = TRUE
  )
  # no initials, only surname
  authors <- trimws(authors)
  authors <- gsub("^([A-Z-][.])+( & ?)?", "", authors, ignore.case = FALSE)
  authors <- gsub("^([A-Z-]+ )+", "", authors, ignore.case = FALSE)
  # remove dots
  authors <- gsub(".", "", authors, fixed = TRUE)
  authors <- gsub("et al", "et al.", authors, fixed = TRUE)
  authors[nchar(authors) <= 1] <- ""
  # combine author and year if year is available
  ref <- if_else(!is.na(lastyear), paste0(authors, ", ", lastyear), authors)
  # fix beginning and ending
  ref <- gsub(", $", "", ref)
  ref <- gsub("^, ", "", ref)
  ref <- gsub("^(emend|et al.,?)", "", ref)
  ref <- trimws(ref)
  ref <- gsub("'", "", ref)

  # a lot start with a lowercase character - fix that
  ref[!grepl("^d[A-Z]", ref)] <- gsub(
    "^([a-z])",
    "\\U\\1",
    ref[!grepl("^d[A-Z]", ref)],
    perl = TRUE
  )
  # specific one for the French that are named dOrbigny
  ref[grepl("^d[A-Z]", ref)] <- gsub("^d", "d'", ref[grepl("^d[A-Z]", ref)])
  ref <- gsub(" +", " ", ref)
  ref <- trimws(ref)
  ref <- gsub("^NA, ?", "", ref)
  ref[ref %in% c("", "NA")] <- NA_character_
  ref
}

# to retrieve LPSN and authors from LPSN website
# e.g., get_lpsn_and_author("genus", "Klebsiella")
get_lpsn_and_author <- function(rank, name) {
  name <- gsub("^Candidatus ", "", name)
  url <- paste0(
    "https://lpsn.dsmz.de/",
    tolower(rank), "/",
    gsub(" ", "-", tolower(name))
  )
  page_txt <- tryCatch({
    resp <- curl::curl_fetch_memory(url)
    if (resp$status_code == 200) {
      read_html(rawToChar(resp$content))
    } else {
      NULL
    }
  }, error = function(e) NULL)
  if (is.null(page_txt)) {
    warning("No LPSN found for ", tolower(rank), " '", name, "'")
    lpsn <- NA_character_
    ref <- NA_character_
    status <- "unknown"
  } else {
    page_txt <- page_txt %>%
      html_element("#detail-page") %>%
      html_text()
    lpsn <- gsub(
      ".*Record number:[\r\n\t ]*([0-9]+).*",
      "\\1",
      page_txt,
      perl = FALSE
    )
    ref <- page_txt %>%
      gsub(".*?Name: (.*[0-9]{4}?).*", "\\1", ., perl = FALSE) %>%
      gsub(name, "", ., fixed = TRUE) %>%
      gsub("^\"?Candidatus ?\"?", "", .) %>%
      trimws()
    status <- trimws(gsub(
      ".*Nomenclatural status:[\r\n\t ]*([a-zA-Z, ]+)[\r\n\t].*",
      "\\1",
      page_txt,
      perl = FALSE
    ))
    if (
      (status %like% "validly published" & status %unlike% "not valid") |
        status %like% "[\r\n\t]"
    ) {
      # we used to take "accepted" for every LPSN record, also candidates. Now only for missing values and explicit accepted ones.
      status <- "accepted"
    } else {
      status <- "not validly published"
    }
  }
  c("lpsn" = lpsn, "ref" = ref, "status" = status)
}

# this will e.g. take the family from the root genus record, and gives all species of that family
get_top_lvl <- function(current, rank, source, rank_target, target) {
  current.bak <- current
  current <- current[target != ""]
  rank <- rank[target != ""]
  source <- source[target != ""]
  target <- target[target != ""]
  if (length(current) == 0) {
    current.bak
  } else if (!rank_target %in% rank) {
    current.bak[1]
  } else {
    if (n_distinct(source) > 1 && "GBIF" %in% source) {
      # prefer LPSN and MycoBank over GBIF, which is often not up-to-date at all
      current <- current[source != "GBIF"]
      rank <- rank[source != "GBIF"]
      source <- source[source != "GBIF"]
    }
    out <- current[rank == rank_target][1]
    if (out %in% c("", NA)) {
      out <- names(sort(
        table(current[which(!current %in% c("", NA))]),
        decreasing = TRUE
      )[1])
      if (is.null(out)) {
        out <- ""
      }
    }
    out
  }
}

# MB/ June 2024: after years still useless, does not contain full taxonomy, e.g. LPSN::request(cred, category = "family") is empty.
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

# Read LPSN data ----------------------------------------------------------------------------------

taxonomy_lpsn.bak <- vroom(file_lpsn, guess_max = 1e5)

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
    status = if_else(is.na(record_lnk), "accepted", "synonym"),
    ref = authors,
    lpsn = as.character(record_no),
    lpsn_parent = NA_character_,
    lpsn_renamed_to = as.character(record_lnk)
  ) %>%
  mutate(source = "LPSN")

# integrity tests
sort(table(taxonomy_lpsn$rank))
sort(table(taxonomy_lpsn$status))
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
  # this will not alter `taxonomy_lpsn` yet
  message("Downloading page ", page, "...", appendLF = TRUE)
  url <- paste0("https://lpsn.dsmz.de/genus?page=", page)
  x <- tryCatch(read_html(url), error = function(e) {
    message("Waiting 10 seconds because of error: ", conditionMessage(e))
    Sys.sleep(10)
    read_html(url)
  })
  x <- x %>%
    # class "main-list" is the main table
    html_element(".main-list") %>%
    # get every list element with a set <id> attribute
    html_elements("li[id]")
  pb <- progress_bar$new(
    total = length(x),
    format = "[:bar] :current/:total :eta"
  )
  for (i in seq_len(length(x))) {
    pb$tick()
    elements <- x[[i]] %>% html_elements("a")
    hrefs <- elements %>% html_attr("href")
    ranks <- hrefs %>% gsub(".*/(.*?)/.*", "\\1", .)
    names <- elements %>%
      html_text() %>%
      gsub('"', "", ., fixed = TRUE)
    # no species, this must be until genus level
    hrefs <- hrefs[ranks != "species"]
    names <- names[ranks != "species"]
    ranks <- ranks[ranks != "species"]

    suppressMessages(
      df <- names %>%
        tibble() %>%
        t() %>%
        as_tibble(.name_repair = "unique") %>%
        setNames(ranks) %>%
        # no candidates please
        filter(genus %unlike% "^(Candidatus|\\[)")
    )

    taxonomy_lpsn_missing <- taxonomy_lpsn_missing %>%
      bind_rows(df)
  }
  message(
    "  => ",
    length(x),
    " entries incl. candidates (cleaned total: ",
    nrow(taxonomy_lpsn_missing),
    ")"
  )
}
taxonomy_lpsn_missing <- taxonomy_lpsn_missing %>% distinct()
saveRDS(taxonomy_lpsn_missing, "data-raw/taxonomy_lpsn_missing.rds")
# taxonomy_lpsn_missing <- readRDS("data-raw/taxonomy_lpsn_missing.rds")

# duplicate genera:
taxonomy_lpsn_missing %>% filter(genus %in% taxonomy_lpsn_missing$genus[duplicated(taxonomy_lpsn_missing$genus)])
# look them up on LPSN, then find out which to keep (the ones validly published under ICNP)
# to remove:
taxonomy_lpsn_missing <- taxonomy_lpsn_missing %>%
  filter(!(genus == "Halalkalibacterium" & family == "Balneolaceae"),
         !(genus == "Pusillimonas" & family == "Oscillospiraceae"),
         !(genus == "Rhodococcus" & family == "Chroococcaceae"))

taxonomy_lpsn <- taxonomy_lpsn %>%
  left_join(taxonomy_lpsn_missing, by = "genus") %>%
  select(domain:family, everything()) %>%
  # remove entries like "[Bacteria, no family]" and "[Bacteria, no class]"
  mutate_all(function(x) if_else(x %like_case% " no ", NA_character_, x))

taxonomy_lpsn.bak2 <- taxonomy_lpsn
# download family taxonomic data (e.g. authors of Enterobacteriaceae) directly from LPSN website using scraping, by using get_lpsn_and_author()
# try it first:
# get_lpsn_and_author("genus", "Escherichia")
# get_lpsn_and_author("family", "Enterobacteriaceae")
pb <- progress_bar$new(
  total = length(unique(taxonomy_lpsn$family)),
  format = "[:bar] :current/:total :eta"
)
for (f in unique(taxonomy_lpsn$family)) {
  pb$tick()
  if (is.na(f)) {
    next
  }
  tax_info <- get_lpsn_and_author("Family", f)
  taxonomy_lpsn <- taxonomy_lpsn %>%
    bind_rows(tibble(
      domain = taxonomy_lpsn$domain[which(taxonomy_lpsn$family == f)[1]],
      kingdom = taxonomy_lpsn$kingdom[which(taxonomy_lpsn$family == f)[1]],
      phylum = taxonomy_lpsn$phylum[which(taxonomy_lpsn$family == f)[1]],
      class = taxonomy_lpsn$class[which(taxonomy_lpsn$family == f)[1]],
      order = taxonomy_lpsn$order[which(taxonomy_lpsn$family == f)[1]],
      family = f,
      rank = "family",
      status = unname(tax_info["status"]),
      source = "LPSN",
      lpsn = unname(tax_info["lpsn"]),
      ref = unname(tax_info["ref"])
    ))
}
# download order directly from LPSN website using scraping
# try it first:
# get_lpsn_and_author("order", "Enterobacterales")
pb <- progress_bar$new(
  total = length(unique(taxonomy_lpsn$order)),
  format = "[:bar] :current/:total :eta"
)
for (o in unique(taxonomy_lpsn$order)) {
  pb$tick()
  if (is.na(o)) {
    next
  }
  tax_info <- get_lpsn_and_author("Order", o)
  taxonomy_lpsn <- taxonomy_lpsn %>%
    bind_rows(tibble(
      domain = taxonomy_lpsn$domain[which(taxonomy_lpsn$order == o)[1]],
      kingdom = taxonomy_lpsn$kingdom[which(taxonomy_lpsn$order == o)[1]],
      phylum = taxonomy_lpsn$phylum[which(taxonomy_lpsn$order == o)[1]],
      class = taxonomy_lpsn$class[which(taxonomy_lpsn$order == o)[1]],
      order = o,
      rank = "order",
      status = unname(tax_info["status"]),
      source = "LPSN",
      lpsn = unname(tax_info["lpsn"]),
      ref = unname(tax_info["ref"])
    ))
}
# download class directly from LPSN website using scraping
# try it first:
# get_lpsn_and_author("class", "Gammaproteobacteria")
pb <- progress_bar$new(
  total = length(unique(taxonomy_lpsn$class)),
  format = "[:bar] :current/:total :eta"
)
for (cc in unique(taxonomy_lpsn$class)) {
  pb$tick()
  if (is.na(cc)) {
    next
  }
  tax_info <- get_lpsn_and_author("Class", cc)
  taxonomy_lpsn <- taxonomy_lpsn %>%
    bind_rows(tibble(
      domain = taxonomy_lpsn$domain[which(taxonomy_lpsn$class == cc)[1]],
      kingdom = taxonomy_lpsn$kingdom[which(taxonomy_lpsn$class == cc)[1]],
      phylum = taxonomy_lpsn$phylum[which(taxonomy_lpsn$class == cc)[1]],
      class = cc,
      rank = "class",
      status = unname(tax_info["status"]),
      source = "LPSN",
      lpsn = unname(tax_info["lpsn"]),
      ref = unname(tax_info["ref"])
    ))
}
# download phylum directly from LPSN website using scraping
# try it first:
# get_lpsn_and_author("phylum", "Pseudomonadota")
pb <- progress_bar$new(
  total = length(unique(taxonomy_lpsn$phylum)),
  format = "[:bar] :current/:total :eta"
)
for (p in unique(taxonomy_lpsn$phylum)) {
  pb$tick()
  if (is.na(p)) {
    next
  }
  tax_info <- get_lpsn_and_author("Phylum", p)
  taxonomy_lpsn <- taxonomy_lpsn %>%
    bind_rows(tibble(
      domain = taxonomy_lpsn$domain[which(taxonomy_lpsn$phylum == p)[1]],
      kingdom = taxonomy_lpsn$kingdom[which(taxonomy_lpsn$phylum == p)[1]],
      phylum = p,
      rank = "phylum",
      status = unname(tax_info["status"]),
      source = "LPSN",
      lpsn = unname(tax_info["lpsn"]),
      ref = unname(tax_info["ref"])
    ))
}
# download kingdom directly from LPSN website using scraping
# try it first:
# get_lpsn_and_author("kingdom", "Pseudomonadati")
pb <- progress_bar$new(
  total = length(unique(taxonomy_lpsn$kingdom)),
  format = "[:bar] :current/:total :eta"
)
for (k in unique(taxonomy_lpsn$kingdom)) {
  pb$tick()
  if (is.na(k)) {
    next
  }
  tax_info <- get_lpsn_and_author("Kingdom", k)
  taxonomy_lpsn <- taxonomy_lpsn %>%
    bind_rows(tibble(
      domain = taxonomy_lpsn$domain[which(taxonomy_lpsn$kingdom == k)[1]],
      kingdom = k,
      rank = "kingdom",
      status = unname(tax_info["status"]),
      source = "LPSN",
      lpsn = unname(tax_info["lpsn"]),
      ref = unname(tax_info["ref"])
    ))
}
# download domain directly from LPSN website using scraping
# try it first:
# get_lpsn_and_author("domain", "Bacteria")
pb <- progress_bar$new(
  total = length(unique(taxonomy_lpsn$domain)),
  format = "[:bar] :current/:total :eta"
)
for (d in unique(taxonomy_lpsn$domain)) {
  pb$tick()
  if (is.na(d)) {
    next
  }
  tax_info <- get_lpsn_and_author("Domain", d)
  taxonomy_lpsn <- taxonomy_lpsn %>%
    bind_rows(tibble(
      domain = d,
      rank = "domain",
      status = unname(tax_info["status"]),
      source = "LPSN",
      lpsn = unname(tax_info["lpsn"]),
      ref = unname(tax_info["ref"])
    ))
}

taxonomy_lpsn <- taxonomy_lpsn %>%
  filter(status != "not validly published")

taxonomy_lpsn <- taxonomy_lpsn %>%
  mutate(ref = get_author_year(ref))

# select final set
taxonomy_lpsn <- taxonomy_lpsn %>%
  select(
    domain,
    kingdom,
    phylum,
    class,
    order,
    family,
    genus,
    species,
    subspecies,
    rank,
    status,
    ref,
    lpsn,
    lpsn_parent,
    lpsn_renamed_to,
    source
  )

# integrity tests
sort(table(taxonomy_lpsn$rank))
# should only be 'accepted' and 'synonym':
sort(table(taxonomy_lpsn$status))
saveRDS(taxonomy_lpsn, "data-raw/taxonomy_lpsn.rds", version = 2)


# Read MycoBank data ------------------------------------------------------------------------------

taxonomy_mycobank <- read_excel(file_mycobank, guess_max = 1e5)
taxonomy_mycobank.bak <- taxonomy_mycobank

taxonomy_mycobank <- taxonomy_mycobank.bak %>%
  mutate(
    clean_name = stringr::str_extract(
      gsub("Current name: ", "", Synonymy, fixed = TRUE),
      "^([A-Z][a-z]+)( [a-z]+)?( [a-z]+[.] [a-z]+)?"
    )
  ) %>%
  transmute(
    mycobank = `MycoBank #`,
    fullname = gsub(" +", " ", `Taxon name`),
    current = clean_name,
    ref = paste0(Authors, ", ", `Year of effective publication`),
    rank = `Rank`,
    status = `Name status`,
    mycobank_renamed_to = taxonomy_mycobank.bak$`MycoBank #`[match(
      clean_name,
      taxonomy_mycobank.bak$`Taxon name`
    )],
    Classification
  ) %>%
  separate(
    Classification,
    sep = ", ",
    into = paste0("tax_", letters[1:20]),
    remove = TRUE
  ) %>%
  mutate(
    rank = case_when(
      fullname %like_case% "^[A-Z][a-z]+ .+ .+" ~ "subsp.",
      rank == "-" & fullname %like_case% "^[A-Z][a-z]+ [a-z-]+$" ~ "sp.",
      rank == "-" &
        paste0(fullname, " ") %in%
          gsub(
            "(^[A-Z][a-z]+ ).*",
            "\\1",
            trimws(taxonomy_mycobank.bak$`Taxon name`),
            perl = TRUE
          ) ~
        "gen.",
      # we take a leap here for the family and order
      rank == "-" & fullname %like% "ceae$" ~ "fam.",
      rank == "-" & fullname %like% "ales$" ~ "ordo",
      # tax_d and tax_e are the subdivision/class columns, kind of; prefer e over d
      rank == "-" & fullname %in% tax_e ~ "cl.",
      rank == "-" & fullname %in% tax_d ~ "cl.",
      TRUE ~ rank
    )
  ) %>%
  filter(
    # remove subkingdom, subfamilies, etc, but keep subspecies
    rank %unlike% "sub" | rank == "subsp.",
    !tolower(rank) %in%
      c("var.", "sect.", "ser.", "tr.", "f.", "race", "stirps", "*", "-"),
    # we also remove orthographic variants here, because our algorithms will give valid results based on misspelling anyway
    !tolower(status) %in%
      c(
        "invalid",
        "deleted",
        "uncertain",
        "illegitimate",
        "orthographic variant",
        "unavailable"
      )
  ) %>%
  mutate(
    mycobank_renamed_to = if_else(
      mycobank_renamed_to == mycobank,
      NA_character_,
      mycobank_renamed_to
    )
  ) %>%
  # only keep 1 entry for the kingdom (regnum)
  filter(fullname == "Fungi" | rank != "regn.") %>%
  # no other kingdoms than fungi
  filter(fullname == "Fungi" | tax_a == "Fungi") %>%
  select_if(function(x) !all(is.na(x)))

taxonomy_mycobank %>% count(rank, sort = TRUE)

table(taxonomy_mycobank$status)
taxonomy_mycobank <- taxonomy_mycobank %>%
  mutate(status = if_else(!is.na(mycobank_renamed_to), "synonym", "accepted"))
table(taxonomy_mycobank$status)

taxonomy_mycobank2 <- taxonomy_mycobank

taxonomy_mycobank <- taxonomy_mycobank2 %>%
  mutate(
    rank = recode_values(
      rank,
      "subsp." ~ "subspecies",
      "sp." ~ "species",
      "gen." ~ "genus",
      "fam." ~ "family",
      "ordo" ~ "order",
      "cl." ~ "class",
      "div." ~ "phylum",
      "regn." ~ "kingdom",
      default = paste0("#", rank)
    )
  )
taxonomy_mycobank %>% count(rank, sort = TRUE)
# some subspecies have wild denotions for the last term, such as 'b.', 'bb.', 'mut.', 'tax.vag.', etc (tens of them)
# we only keep these
taxonomy_mycobank <- taxonomy_mycobank %>%
  filter(
    rank != "subspecies" |
      fullname %like% " (f|f\\.sp|ssp|subf|subsp|subvar|var)[.] "
  ) %>%
  mutate(
    fullname = ifelse(
      rank == "subspecies",
      # take the 1st, 2nd and 4th term
      sapply(strsplit(fullname, " "), function(x) {
        paste(x[c(1, 2, 4)], collapse = " ")
      }),
      fullname
    )
  ) %>%
  # and make fullname distinct again
  arrange(rank, mycobank_renamed_to) %>%
  distinct(fullname, .keep_all = TRUE) %>%
  arrange(fullname)

taxonomy_mycobank %>% count(rank, sort = TRUE)
taxonomy_mycobank %>%
  filter(rank %like% "#") %>%
  count(rank)

taxonomy_mycobank <- taxonomy_mycobank %>%
  filter(rank %unlike% "#")

taxonomy_mycobank3 <- taxonomy_mycobank

taxonomy_mycobank <- taxonomy_mycobank3
# MycoBank just pasted their taxonomy together into 1 field, and now some classes are in the division (tax_b) column, it's horrible
# so we decide based on the fullname and rank column per record
# use this to determine how far to go:
any(
  taxonomy_mycobank$fullname[taxonomy_mycobank$rank == "class"] %in%
    taxonomy_mycobank$tax_e
) # replace tax_e with tax_f, etc
any(
  taxonomy_mycobank$fullname[taxonomy_mycobank$rank == "class"] %in%
    taxonomy_mycobank$tax_f
)
taxonomy_mycobank <- taxonomy_mycobank %>%
  mutate(
    kingdom = "Fungi", # we already filtered everything else, and MycoBank 90157 has '-' as current name...
    phylum = case_when(
      rank == "phylum" ~ fullname,
      tax_b %in%
        taxonomy_mycobank$fullname[taxonomy_mycobank$rank == "phylum"] ~
        tax_b,
      tax_c %in%
        taxonomy_mycobank$fullname[taxonomy_mycobank$rank == "phylum"] ~
        tax_c,
      TRUE ~ ""
    ),
    class = case_when(
      rank == "class" ~ fullname,
      tax_c %in% taxonomy_mycobank$fullname[taxonomy_mycobank$rank == "class"] ~
        tax_c,
      tax_d %in% taxonomy_mycobank$fullname[taxonomy_mycobank$rank == "class"] ~
        tax_d,
      tax_e %in% taxonomy_mycobank$fullname[taxonomy_mycobank$rank == "class"] ~
        tax_e,
      TRUE ~ ""
    ),
    order = case_when(
      rank == "order" ~ fullname,
      tax_c %in% taxonomy_mycobank$fullname[taxonomy_mycobank$rank == "order"] ~
        tax_c,
      tax_d %in% taxonomy_mycobank$fullname[taxonomy_mycobank$rank == "order"] ~
        tax_d,
      tax_e %in% taxonomy_mycobank$fullname[taxonomy_mycobank$rank == "order"] ~
        tax_e,
      tax_f %in% taxonomy_mycobank$fullname[taxonomy_mycobank$rank == "order"] ~
        tax_f,
      tax_g %in% taxonomy_mycobank$fullname[taxonomy_mycobank$rank == "order"] ~
        tax_g,
      TRUE ~ ""
    ),
    family = case_when(
      rank == "family" ~ fullname,
      tax_c %in%
        taxonomy_mycobank$fullname[taxonomy_mycobank$rank == "family"] ~
        tax_c,
      tax_d %in%
        taxonomy_mycobank$fullname[taxonomy_mycobank$rank == "family"] ~
        tax_d,
      tax_e %in%
        taxonomy_mycobank$fullname[taxonomy_mycobank$rank == "family"] ~
        tax_e,
      tax_f %in%
        taxonomy_mycobank$fullname[taxonomy_mycobank$rank == "family"] ~
        tax_f,
      tax_g %in%
        taxonomy_mycobank$fullname[taxonomy_mycobank$rank == "family"] ~
        tax_g,
      tax_h %in%
        taxonomy_mycobank$fullname[taxonomy_mycobank$rank == "family"] ~
        tax_h,
      tax_i %in%
        taxonomy_mycobank$fullname[taxonomy_mycobank$rank == "family"] ~
        tax_i,
      TRUE ~ ""
    ),
    genus = case_when(
      rank == "genus" ~ fullname,
      tax_b %in% taxonomy_mycobank$fullname[taxonomy_mycobank$rank == "genus"] ~
        tax_b,
      tax_c %in% taxonomy_mycobank$fullname[taxonomy_mycobank$rank == "genus"] ~
        tax_c,
      tax_d %in% taxonomy_mycobank$fullname[taxonomy_mycobank$rank == "genus"] ~
        tax_d,
      tax_e %in% taxonomy_mycobank$fullname[taxonomy_mycobank$rank == "genus"] ~
        tax_e,
      tax_f %in% taxonomy_mycobank$fullname[taxonomy_mycobank$rank == "genus"] ~
        tax_f,
      tax_g %in% taxonomy_mycobank$fullname[taxonomy_mycobank$rank == "genus"] ~
        tax_g,
      tax_h %in% taxonomy_mycobank$fullname[taxonomy_mycobank$rank == "genus"] ~
        tax_h,
      tax_i %in% taxonomy_mycobank$fullname[taxonomy_mycobank$rank == "genus"] ~
        tax_i,
      tax_j %in% taxonomy_mycobank$fullname[taxonomy_mycobank$rank == "genus"] ~
        tax_j,
      tax_k %in% taxonomy_mycobank$fullname[taxonomy_mycobank$rank == "genus"] ~
        tax_k,
      TRUE ~ ""
    ),
    species = case_when(
      rank == "species" & fullname %like% " " ~
        gsub(".* (.*)", "\\1", fullname, perl = TRUE),
      rank == "subspecies" & fullname %like% " " ~
        gsub(".* (.*) .*", "\\1", fullname, perl = TRUE),
      TRUE ~ ""
    ),
    subspecies = case_when(
      rank == "subspecies" & fullname %like% " " ~
        gsub(".* .* (.*)", "\\1", fullname, perl = TRUE),
      TRUE ~ ""
    )
  )

# use this to get all the genera with updated names from MO_RELEVANT_GENERA:
genera_overview <- AMR::microorganisms %>%
  filter(genus %in% AMR:::MO_RELEVANT_GENERA) %>%
  pull(fullname) %>%
  mo_current() %>% # renamed to's will be included below
  mo_genus() %>%
  unique() %>%
  sort()

# keep only the relevant ones
include <- taxonomy_mycobank %>%
  filter(
    genus %in% genera_overview | !rank %in% c("genus", "species", "subspecies")
  ) %>%
  filter(!(genus == "" & rank %in% c("genus", "species", "subspecies")))
missing_renamed_to <- include$mycobank_renamed_to[which(
  !is.na(include$mycobank_renamed_to) &
    !include$mycobank_renamed_to %in% include$mycobank
)]
taxonomy_mycobank <- include %>%
  bind_rows(taxonomy_mycobank %>% filter(mycobank %in% missing_renamed_to)) %>%
  arrange(fullname)

# clean up authors and add last columns
taxonomy_mycobank <- taxonomy_mycobank %>%
  mutate(
    domain = kingdom,
    source = "MycoBank",
    mycobank_parent = NA_character_,
    ref = get_author_year(ref)
  )

# select final set
taxonomy_mycobank <- taxonomy_mycobank %>%
  select(
    fullname,
    domain,
    kingdom,
    phylum,
    class,
    order,
    family,
    genus,
    species,
    subspecies,
    rank,
    status,
    ref,
    mycobank,
    mycobank_parent,
    mycobank_renamed_to,
    source
  )

# not all 'renamed to' records are available, some were even just orthographic variants (that have an invalid status)
taxonomy_mycobank$status[
  !is.na(taxonomy_mycobank$mycobank_renamed_to) &
    !taxonomy_mycobank$mycobank_renamed_to %in% taxonomy_mycobank$mycobank
] <- "accepted"
taxonomy_mycobank$mycobank_renamed_to[
  taxonomy_mycobank$status == "accepted"
] <- NA

taxonomy_mycobank %>% count(status)
saveRDS(taxonomy_mycobank, "data-raw/taxonomy_mycobank.rds", version = 2)


# Read GBIF data ----------------------------------------------------------------------------------

taxonomy_gbif.bak <- vroom(file_gbif, guess_max = 5e5)
colnames(taxonomy_gbif.bak) <- gsub(".*:(.*)", "\\1", colnames(taxonomy_gbif.bak))

# include all fungal orders from the mycobank db
include_fungal_orders <- unique(taxonomy_mycobank$order[
  !taxonomy_mycobank$order %in% c("", NA)
])

# check some columns to validate below filters
taxonomy_gbif.bak %>% count(taxonomicStatus, sort = TRUE)
taxonomy_gbif.bak %>% count(taxonRank, sort = TRUE)
taxonomy_gbif.bak %>% count(kingdom, sort = TRUE)

taxonomy_gbif <- taxonomy_gbif.bak %>%
  mutate(
    # strip the authors from the scientific name
    fullname = case_when(
      is.na(scientificNameAuthorship) | scientificNameAuthorship == "" ~ scientificName,
      endsWith(scientificName, scientificNameAuthorship) ~
        trimws(substr(scientificName, 1L, nchar(scientificName) - nchar(scientificNameAuthorship))),
      TRUE ~ scientificName
    )
  ) %>%
  mutate(
    genus = if_else(taxonRank == "genus", fullname, genus),
    family = if_else(taxonRank == "family", fullname, family),
    order = if_else(taxonRank == "order", fullname, order),
    class = if_else(taxonRank == "class", fullname, class),
    phylum = if_else(taxonRank == "phylum", fullname, phylum),
    kingdom = if_else(taxonRank == "kingdom", fullname, kingdom)
  )
taxonomy_gbif0 <- taxonomy_gbif

taxonomy_gbif <- taxonomy_gbif0 %>%
  # immediately filter rows we really never want
  filter(
    # never doubtful status, only accepted and all synonyms, and only ranked items
    !taxonomicStatus %in% c("doubtful", "misapplied", "ambiguous synonym"),
    taxonRank != "unranked",
    #  include these kingdoms (no Chromista)
    kingdom %in%
      c("Archaea", "Bacteria", "Protozoa", na.omit(unique(taxonomy_lpsn$kingdom))) |
      # include all of these fungal orders
      order %in% include_fungal_orders |
      # and all of these important genera (see "data-raw/_pre_commit_checks.R")
      # (they also contain bacteria and protozoa, but these will get higher prevalence scores later on)
      genus %in% AMR:::MO_RELEVANT_GENERA
  ) %>%
  mutate(
    # set the right domain for the prokaryotes based on LPSN
    domain = case_when(
      kingdom %in% na.omit(taxonomy_lpsn$kingdom[taxonomy_lpsn$domain == "Bacteria"]) ~ "Bacteria",
      kingdom %in% na.omit(taxonomy_lpsn$kingdom[taxonomy_lpsn$domain == "Archaea"]) ~ "Archaea",
      TRUE ~ kingdom
    )
  )
taxonomy_gbif <- taxonomy_gbif %>%
  bind_rows(taxonomy_gbif0 %>%
              filter(acceptedNameUsageID %in% taxonomy_gbif$taxonID,
                     !taxonomicStatus %in% c("doubtful", "misapplied", "ambiguous synonym")) %>%
              mutate(status = "synonym"))

taxonomy_gbif <- taxonomy_gbif %>%
  select(
    fullname,
    domain,
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
  )
taxonomy_gbif1 <- taxonomy_gbif

taxonomy_gbif <- taxonomy_gbif1 %>%
  mutate(
    status = if_else(status %like% "accepted", "accepted", "synonym"),
    # checked taxonRank - the "form" and "variety" always have a subspecies
    # see: taxonomy_gbif.bak %>% filter(taxonRank %in% c("form", "variety")) %>% count(taxonRank, is.na(infraspecificEpithet), sort = TRUE)
    rank = if_else(rank %in% c("form", "variety"), "subspecies", rank),
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
  distinct(
    domain,
    kingdom,
    fullname,
    .keep_all = TRUE
  ) %>%
  filter(
    domain %unlike% "[0-9]",
    kingdom %unlike% "[0-9]",
    phylum %unlike% "[0-9]",
    class %unlike% "[0-9]",
    order %unlike% "[0-9]",
    family %unlike% "[0-9]",
    genus %unlike% "[0-9]"
  )

# fix the missing synonym taxonomy:

rank_levels <- c("kingdom", "phylum", "class", "order", "family", "genus", "species", "subspecies")

accepted <- taxonomy_gbif %>%
  filter(!is.na(gbif)) %>%
  select(
    gbif,
    acc_domain = domain,
    acc_kingdom = kingdom,
    acc_phylum = phylum,
    acc_class = class,
    acc_order = order,
    acc_family = family
  ) %>%
  distinct(gbif, .keep_all = TRUE)

taxonomy_gbif <- taxonomy_gbif %>%
  left_join(accepted, by = c("gbif_renamed_to" = "gbif")) %>%
  mutate(
    name_parts = strsplit(fullname, " "),
    rank_num = match(rank, rank_levels, nomatch = NA_integer_),
    
    # Step 1: parse the synonym's own name components from fullname
    kingdom = if_else(is.na(kingdom) & rank == "kingdom",
                      sapply(name_parts, `[`, 1L), kingdom),
    phylum = if_else(is.na(phylum) & rank == "phylum",
                     sapply(name_parts, `[`, 1L), phylum),
    class = if_else(is.na(class) & rank == "class",
                    sapply(name_parts, `[`, 1L), class),
    order = if_else(is.na(order) & rank == "order",
                    sapply(name_parts, `[`, 1L), order),
    family = if_else(is.na(family) & rank == "family",
                     sapply(name_parts, `[`, 1L), family),
    genus = if_else(is.na(genus) & rank %in% c("genus", "species", "subspecies"),
                    sapply(name_parts, `[`, 1L), genus),
    species = if_else(is.na(species) & rank %in% c("species", "subspecies"),
                      sapply(name_parts, `[`, 2L), species),
    subspecies = if_else(is.na(subspecies) & rank == "subspecies",
                         sapply(name_parts, `[`, 3L), subspecies),
    
    # Step 2: fill higher taxonomy from accepted name, only for ranks ABOVE the synonym's rank
    domain  = if_else(is.na(domain) & !is.na(acc_domain), acc_domain, domain),
    kingdom = if_else(is.na(kingdom) & !is.na(acc_kingdom) & rank_num > 1L, acc_kingdom, kingdom),
    phylum  = if_else(is.na(phylum) & !is.na(acc_phylum) & rank_num > 2L, acc_phylum, phylum),
    class   = if_else(is.na(class) & !is.na(acc_class) & rank_num > 3L, acc_class, class),
    order   = if_else(is.na(order) & !is.na(acc_order) & rank_num > 4L, acc_order, order),
    family  = if_else(is.na(family) & !is.na(acc_family) & rank_num > 5L, acc_family, family)
  ) %>%
  select(-starts_with("acc_"), -name_parts, -rank_num)


taxonomy_gbif <- taxonomy_gbif %>%
  mutate(ref = get_author_year(ref)) %>%
  filter(ref != "AmSOD")

# integrity tests
sort(table(taxonomy_gbif$rank))
sort(table(taxonomy_gbif$status))

saveRDS(taxonomy_gbif, "data-raw/taxonomy_gbif.rds", version = 2)


# *** Saved intermediate results (all taxonomies) *** ---------------------------------------------

# *** SaveRDS from above, this does not need to be run:
taxonomy_lpsn <- readRDS("data-raw/taxonomy_lpsn.rds")
taxonomy_mycobank <- readRDS("data-raw/taxonomy_mycobank.rds")
taxonomy_gbif <- readRDS("data-raw/taxonomy_gbif.rds")
# this just allows to always get back to this point by simply loading the files from data-raw/.


# Add full names ----------------------------------------------------------------------------------

taxonomy_gbif <- taxonomy_gbif %>%
  # clean NAs and add fullname
  mutate(
    across(domain:subspecies, function(x) if_else(is.na(x), "", x)),
    fullname = trimws(case_when(
      rank == "family" ~ family,
      rank == "order" ~ order,
      rank == "class" ~ class,
      rank == "phylum" ~ phylum,
      rank == "kingdom" ~ kingdom,
      rank == "domain" ~ domain,
      TRUE ~ paste(genus, species, subspecies) # already trimmed 6 lines up
    )),
    .before = 1
  ) %>%
  # keep only one GBIF taxon ID per full name
  arrange(fullname, gbif) %>%
  distinct(kingdom, rank, fullname, .keep_all = TRUE)

taxonomy_lpsn <- taxonomy_lpsn %>%
  # clean NAs and add fullname
  mutate(
    across(domain:subspecies, function(x) if_else(is.na(x), "", x)),
    fullname = trimws(case_when(
      rank == "family" ~ family,
      rank == "order" ~ order,
      rank == "class" ~ class,
      rank == "phylum" ~ phylum,
      rank == "kingdom" ~ kingdom,
      rank == "domain" ~ domain,
      TRUE ~ paste(genus, species, subspecies) # already trimmed 6 lines up
    )),
    .before = 1
  ) %>%
  # keep only one LPSN record ID per full name
  arrange(fullname, lpsn) %>%
  distinct(kingdom, rank, fullname, .keep_all = TRUE)

taxonomy_mycobank <- taxonomy_mycobank %>%
  # clean NAs and add fullname
  mutate(
    across(domain:subspecies, function(x) if_else(is.na(x), "", x)),
    fullname = trimws(case_when(
      rank == "family" ~ family,
      rank == "order" ~ order,
      rank == "class" ~ class,
      rank == "phylum" ~ phylum,
      rank == "kingdom" ~ kingdom,
      rank == "domain" ~ domain,
      TRUE ~ paste(genus, species, subspecies) # already trimmed 7 lines up
    )),
    .before = 1
  ) %>%
  # keep only one MycoBank record ID per full name
  arrange(fullname, mycobank) %>%
  distinct(kingdom, rank, fullname, .keep_all = TRUE)


# Combine the datasets ----------------------------------------------------------------------------

taxonomy <- taxonomy_lpsn %>%
  filter(!domain %in% c("", NA)) %>%
  # add fungi
  bind_rows(taxonomy_mycobank %>% filter(!domain %in% c("", NA))) %>%
  # add GBIF to the bottom
  bind_rows(taxonomy_gbif %>% filter(!domain %in% c("", NA))) %>%
  # group on unique species
  group_by(domain, fullname) %>%
  # fill the NAs in LPSN/GBIF fields and ref with the other source (so LPSN: 123 and GBIF: NA will become LPSN: 123 and GBIF: 123)
  mutate(across(matches("^(lpsn|mycobank|gbif|ref)"), function(x) {
    rep(x[!is.na(x)][1], length(x))
  })) %>%
  # ungroup again
  ungroup() %>%
  # only keep unique species per domain
  distinct(domain, fullname, .keep_all = TRUE) %>%
  arrange(fullname) %>% 
  select(fullname, everything())

# get missing entries from existing microorganisms data set
taxonomy.old <- AMR::microorganisms %>%
  select(any_of(colnames(taxonomy))) %>%
  # TODO ONLY IN 2026!
  mutate(domain = kingdom) %>% 
  # END OF TODO
  filter(
    !paste(domain, fullname) %in% paste(taxonomy$domain, taxonomy$fullname),
    # these will be added later:
    tolower(source) != "manually added"
  )
taxonomy <- taxonomy %>%
  bind_rows(taxonomy.old) %>%
  arrange(fullname) %>%
  filter(fullname != "")

# fix rank
taxonomy %>% count(rank, sort = TRUE)
taxonomy <- taxonomy %>%
  mutate(
    rank = case_when(
      subspecies != "" ~ "subspecies",
      species != "" ~ "species",
      genus != "" ~ "genus",
      family != "" ~ "family",
      order != "" ~ "order",
      class != "" ~ "class",
      phylum != "" ~ "phylum",
      kingdom != "" ~ "kingdom",
      domain != "" ~ "domain",
      TRUE ~ NA_character_
    )
  )
taxonomy %>% count(rank, sort = TRUE)

# Resolve genera that appear in multiple domains (e.g. Giardia in Protozoa
# and Animalia). Prefer Bacteria > Fungi > Protozoa > Archaea > Chromista >
# Animalia, then overwrite domain-to-family for the entire genus.
taxonomy %>% count(domain, sort = TRUE)
domain_priority <- c(
  "Bacteria" = 1L,
  "Fungi" = 2L,
  "Protozoa" = 3L,
  "Archaea" = 4L,
  "Chromista" = 5L,
  "Animalia" = 6L,
  "Plantae" = 7L
)

# find the best domain per genus from genus-rank records
best_domain <- taxonomy %>%
  filter(rank == "genus", genus != "", domain %in% names(domain_priority)) %>%
  mutate(
    dprio = domain_priority[domain],
    sprio = recode_values(source,
                          "LPSN" ~ 1L,
                          "MycoBank" ~ 2L,
                          "GBIF" ~ 3L,
                          default = 4L)) %>%
  arrange(genus, dprio, sprio) %>%
  distinct(genus, .keep_all = TRUE) %>%
  select(genus, best_domain = domain, best_kingdom = kingdom,
         best_phylum = phylum, best_class = class,
         best_order = order, best_family = family)

taxonomy <- taxonomy %>%
  left_join(best_domain, by = "genus") %>%
  mutate(
    domain  = if_else(genus != "" & !is.na(best_domain), best_domain, domain),
    kingdom = if_else(genus != "" & !is.na(best_kingdom), best_kingdom, kingdom),
    phylum  = if_else(genus != "" & !is.na(best_phylum) & rank %in%
                        c("genus", "species", "subspecies"), best_phylum, phylum),
    class   = if_else(genus != "" & !is.na(best_class) & rank %in%
                        c("genus", "species", "subspecies"), best_class, class),
    order   = if_else(genus != "" & !is.na(best_order) & rank %in%
                        c("genus", "species", "subspecies"), best_order, order),
    family  = if_else(genus != "" & !is.na(best_family) & rank %in%
                        c("genus", "species", "subspecies"), best_family, family)) %>%
  select(-starts_with("best_"))

# check if some genera within kingdoms have multiple families / orders, etc:
taxonomy %>%
  filter(genus != "") %>%
  group_by(kingdom, genus) %>%
  filter(n_distinct(family) > 1) %>%
  View()
# then fix (we're still arranged by domain, source here):
taxonomy <- taxonomy %>%
  group_by(kingdom, genus) %>%
  mutate(family = get_top_lvl(family, rank, source, "genus", genus)) %>%
  group_by(kingdom, family) %>%
  mutate(order = get_top_lvl(order, rank, source, "family", family)) %>%
  group_by(kingdom, order) %>%
  mutate(class = get_top_lvl(class, rank, source, "order", order)) %>%
  group_by(kingdom, class) %>%
  mutate(phylum = get_top_lvl(phylum, rank, source, "class", class)) %>%
  ungroup()

# remove the taxonomy where it must remain empty
taxonomy <- taxonomy %>%
  mutate(
    kingdom = if_else(rank %in% c("domain"), "", kingdom),
    phylum = if_else(rank %in% c("domain", "kingdom"), "", phylum),
    class = if_else(rank %in% c("domain", "kingdom", "phylum"), "", class),
    order = if_else(rank %in% c("domain", "kingdom", "phylum", "class"), "", order),
    family = if_else(
      rank %in% c("domain", "kingdom", "phylum", "class", "order"),
      "",
      family
    ),
    genus = if_else(
      rank %in% c("domain", "kingdom", "phylum", "class", "order", "family"),
      "",
      genus
    ),
    species = if_else(
      rank %in% c("domain", "kingdom", "phylum", "class", "order", "family", "genus"),
      "",
      species
    ),
    subspecies = if_else(
      rank %in%
        c("domain", "kingdom", "phylum", "class", "order", "family", "genus", "species"),
      "",
      subspecies
    )
  )

# recreate fullnames and keep unique - we're still arranged by domain, source here
taxonomy <- taxonomy %>%
  mutate(
    fullname = trimws(case_when(
      rank == "family" ~ family,
      rank == "order" ~ order,
      rank == "class" ~ class,
      rank == "phylum" ~ phylum,
      rank == "kingdom" ~ kingdom,
      rank == "domain" ~ domain,
      TRUE ~ paste(genus, species, subspecies) # already trimmed 7 lines up
    ))) %>% 
  arrange(fullname) %>% 
  distinct(fullname, .keep_all = TRUE)


# *** Save intermediate results (0) *** -----------------------------------------------------------

saveRDS(taxonomy, "data-raw/taxonomy0.rds")
# taxonomy <- readRDS("data-raw/taxonomy0.rds")


# Add missing taxonomic entries and deduplicate ---------------------------------------------------
# Ensure every referenced rank has its own row (domain through species).
# Where possible, enrich with GBIF identifiers.

current_gbif <- taxonomy_gbif.bak %>%
  filter(is.na(acceptedNameUsageID)) %>%
  mutate(
    taxonID = as.character(taxonID),
    parentNameUsageID = as.character(parentNameUsageID),
    domain = case_when(
      kingdom %in% na.omit(taxonomy_lpsn$kingdom[taxonomy_lpsn$domain == "Bacteria"]) ~ "Bacteria",
      kingdom %in% na.omit(taxonomy_lpsn$kingdom[taxonomy_lpsn$domain == "Archaea"]) ~ "Archaea",
      TRUE ~ kingdom
    )
  )

rank_levels <- c("domain", "kingdom", "phylum", "class", "order", "family", "genus")

taxonomy_all_missing <- map_dfr(rank_levels, function(rank_name) {
  candidates <- taxonomy %>%
    filter(.data[[rank_name]] != "") %>%
    distinct(across(all_of(unique(c("domain", rank_levels[seq_len(which(rank_levels == rank_name))]))))) %>%
    mutate(fullname = .data[[rank_name]], rank = rank_name) %>%
    anti_join(
      taxonomy %>% filter(rank == rank_name),
      by = unique(c("domain", setNames(rank_name, rank_name), "rank"))
    )
  
  if (nrow(candidates) == 0) return(tibble())
  
  message("Adding missing: ", rank_name, "... n = ", nrow(candidates))
  
  # enrich from GBIF
  gbif_lookup <- current_gbif %>%
    filter(taxonRank == rank_name) %>%
    select(all_of(unique(c("domain", rank_name))),
           ref = scientificNameAuthorship,
           gbif = taxonID,
           gbif_parent = parentNameUsageID) %>%
    distinct()
  
  candidates %>%
    left_join(gbif_lookup, by = unique(c("domain", rank_name))) %>%
    mutate(
      source = if_else(!is.na(gbif), "GBIF", "manually added"),
      status = if_else(!is.na(gbif), "accepted", "unknown")
    )
})

# species implied by subspecies but missing as a species-rank row
missing_species <- taxonomy %>%
  filter(species != "") %>%
  distinct(domain, genus, species, .keep_all = TRUE) %>%
  select(domain:species) %>%
  mutate(fullname = paste(genus, species), rank = "species") %>%
  anti_join(
    taxonomy %>% filter(rank == "species"),
    by = c("domain", "genus", "species", "rank")
  ) %>%
  left_join(
    current_gbif %>%
      filter(taxonRank == "species") %>%
      select(domain, genus, species = specificEpithet,
             ref = scientificNameAuthorship,
             gbif = taxonID, gbif_parent = parentNameUsageID) %>%
      distinct(),
    by = c("domain", "genus", "species")
  ) %>%
  mutate(
    source = if_else(!is.na(gbif), "GBIF", "manually added"),
    status = if_else(!is.na(gbif), "accepted", "unknown")
  )

taxonomy <- taxonomy %>%
  bind_rows(taxonomy_all_missing, missing_species) %>%
  mutate(
    domain = if_else(is.na(domain) | domain == "", kingdom, domain),
    across(kingdom:subspecies, \(x) if_else(is.na(x), "", x))
  )

# TODO perhaps only for 2026: remove old kingdoms
taxonomy <- taxonomy %>% 
  filter(!(rank == "kingdom" & fullname %in% c("Bacteria", "Archaea")))

# Deduplicate: where the same fullname appears at multiple ranks within a
# domain (e.g. Nitrospira as genus and class), keep the lowest rank as-is
# and disambiguate higher ranks by appending {rank}.
rank_priority <- c("subspecies" = 1L, "species" = 2L, "genus" = 3L,
                   "family" = 4L, "order" = 5L, "class" = 6L,
                   "phylum" = 7L, "kingdom" = 8L, "domain" = 9L)
source_priority <- c("LPSN" = 1L, "MycoBank" = 2L, "GBIF" = 3L, "manually added" = 4L)

taxonomy <- taxonomy %>%
  mutate(rank_index = rank_priority[rank],
         source_index = source_priority[source]) %>%
  arrange(domain, fullname, rank_index, source_index) %>%
  distinct(domain, fullname, rank_index, .keep_all = TRUE) %>% 
  group_by(domain, fullname) %>%
  mutate(fullname = if_else(
    row_number() > 1,
    paste0(fullname, " {", rank, "}"),
    fullname
  )) %>%
  ungroup() %>%
  select(-rank_index, -source_index) %>%
  arrange(domain, fullname, ref) %>%
  distinct(domain, fullname, .keep_all = TRUE) %>%
  filter(domain != "")


# *** Save intermediate results (1) *** -----------------------------------------------------------

saveRDS(taxonomy, "data-raw/taxonomy1.rds")
# taxonomy <- readRDS("data-raw/taxonomy1.rds")


# Get previously manually added entries -----------------------------------------------------------

manually_added <- AMR::microorganisms %>%
  filter(
    tolower(source) == "manually added",
    !paste(kingdom, fullname) %in% paste(taxonomy$kingdom, taxonomy$fullname),
    !rank %in% c("domain", "kingdom", "phylum", "class", "order", "family")
  ) %>%
  select(fullname:subspecies, ref, source, rank)

# Build a lookup for each child->parent rank pair, preferring LPSN > MycoBank > any
source_priority <- c("LPSN" = 1L, "MycoBank" = 2L, "GBIF" = 3L, "manually added" = 4L)

fill_from_parent <- function(data, taxonomy, child_rank, parent_rank) {
  lookup <- taxonomy %>%
    filter(.data[[child_rank]] != "", .data[[parent_rank]] != "") %>%
    mutate(sprio = source_priority[source]) %>%
    arrange(.data[[child_rank]], sprio) %>%
    distinct(.data[[child_rank]], .keep_all = TRUE) %>%
    select(all_of(c(child_rank, parent_rank)))
  
  data %>%
    rows_update(lookup, by = child_rank, unmatched = "ignore")
}

# Walk up the hierarchy: genus->family->order->class->phylum->kingdom
rank_pairs <- tibble(
  child  = c("genus", "family", "order", "class", "phylum"),
  parent = c("family", "order", "class", "phylum", "kingdom")
)

for (i in seq_len(nrow(rank_pairs))) {
  manually_added <- fill_from_parent(
    manually_added, taxonomy,
    rank_pairs$child[i], rank_pairs$parent[i]
  )
}

manually_added <- manually_added %>%
  mutate(
    status = "unknown",
    rank = if_else(fullname %like% "unknown", "(unknown rank)", rank)
  ) %>%
  filter(!fullname %in% taxonomy$fullname)

taxonomy <- taxonomy %>%
  bind_rows(manually_added) %>%
  arrange(fullname)

table(taxonomy$rank, useNA = "always")


# Remove childless higher-rank entries (family through kingdom) -----------------------------------

# Work bottom-up: family first, then order, class, phylum, kingdom.
# After each step, the next rank up may have lost its last child.
for (rank_name in c("family", "order", "class", "phylum", "kingdom")) {
  has_children <- taxonomy %>%
    filter(rank != rank_name, .data[[rank_name]] != "") %>%
    distinct(.data[[rank_name]]) %>%
    pull()
  
  n_before <- nrow(taxonomy)
  taxonomy <- taxonomy %>%
    filter(!(rank == rank_name & !fullname %in% has_children))
  message("Removed ", n_before - nrow(taxonomy), " childless ", rank_name, " entries")
}


# Get LPSN data for records missing from `taxonomy_lpsn` ------------------------------------------

# Weirdly enough, some LPSN records are lacking from the API and the CSV file (i.e., `taxonomy_lpsn`),
# such as family Thiotrichaceae and its order Thiotrichales, or the genus Coleospermum. When running
# get_lpsn_and_author("family", "Thiotrichaceae") you do get a result, vs. taxonomy_lpsn %>% filter(family == "Thiotrichaceae").
# So check every non-LPSN records from the kingdom of Bacteria and add it
lpsn_genera <- taxonomy %>%
  filter(source == "LPSN", rank == "genus") %>%
  pull(genus)
gbif_bacteria <- which(
  taxonomy$domain == "Bacteria" &
    taxonomy$source %in% c("GBIF", "manually added") &
    (taxonomy$rank %in% c("phylum", "class", "order", "family", "genus") |
       (taxonomy$rank == "species" & taxonomy$genus %in% lpsn_genera))
)
message(length(gbif_bacteria))
added <- 0
pb <- progress_bar$new(
  total = length(gbif_bacteria),
  format = "[:bar] :current/:total :eta"
)
for (record in gbif_bacteria) {
  pb$tick()
  Sys.sleep(0.1)
  lpsn <- get_lpsn_and_author(
    rank = taxonomy$rank[record],
    name = taxonomy$fullname[record]
  )
  if (is.na(lpsn["lpsn"])) {
    next
  } else {
    added <- added + 1
    taxonomy$source[record] <- "LPSN"
    taxonomy$lpsn[record] <- unname(lpsn["lpsn"])
    taxonomy$ref[record] <- unname(lpsn["ref"])
    taxonomy$status[record] <- unname(lpsn["status"])
  }
}
warnings()
message(added, " GBIF records altered to latest LPSN")

# parent LPSNs will be added later


# *** Save intermediate results (1b) *** ----------------------------------------------------------

saveRDS(taxonomy, "data-raw/taxonomy1b.rds")
# taxonomy <- readRDS("data-raw/taxonomy1b.rds")


# Clean scientific reference ----------------------------------------------------------------------

taxonomy <- taxonomy %>%
  mutate(ref = get_author_year(ref))


# Get the latest upper taxonomy from LPSN/MycoBank for GBIF data ----------------------------------

# this fix is required after the LPSN import of get_lpsn_and_author()
taxonomy <- taxonomy %>%
  mutate(
    # set the right domain for the prokaryotes based on LPSN
    domain = case_when(
      kingdom %in% na.omit(taxonomy_lpsn$kingdom[taxonomy_lpsn$domain == "Bacteria"]) ~ "Bacteria",
      kingdom %in% na.omit(taxonomy_lpsn$kingdom[taxonomy_lpsn$domain == "Archaea"]) ~ "Archaea",
      TRUE ~ kingdom
    )
  )
taxonomy$domain[taxonomy$domain == "(unknown kingdom)"] <- "(unknown domain)"

# Fix GBIF taxonomy using authoritative source per domain -----------------------------------------

# e.g., phylum above class "Bacilli" was still "Firmicutes" in 2023, should be "Bacillota" per LPSN
rank_pairs <- tibble(
  child  = c("genus", "family", "order", "class", "phylum"),
  parent = c("family", "order",  "class", "phylum", "kingdom")
)

for (d in unique(taxonomy$domain[taxonomy$domain != "(unknown domain)"])) {
  src <- if (d == "Fungi") "MycoBank" else "LPSN"
  message("Fixing GBIF taxonomy for domain ", d, " based on ", src, "...",
          appendLF = FALSE)
  
  for (i in seq_len(nrow(rank_pairs))) {
    child_col  <- rank_pairs$child[i]
    parent_col <- rank_pairs$parent[i]
    
    # build lookup: for each child value in this domain, get the parent from the authoritative source
    lookup <- taxonomy %>%
      filter(
        domain == d,
        source == src,
        .data[[child_col]] != ""
      )
    
    # special case: kingdom lookup excludes "Bacteria"/"Archaea" as kingdom values
    if (parent_col == "kingdom") {
      lookup <- lookup %>%
        filter(!kingdom %in% c("Bacteria", "Archaea"))
    }
    
    lookup <- lookup %>%
      distinct(.data[[child_col]], .keep_all = TRUE) %>%
      select(all_of(c(child_col, parent_col)))
    
    # join and overwrite
    taxonomy <- taxonomy %>%
      rows_update(
        taxonomy %>%
          filter(domain == d, .data[[child_col]] != "") %>%
          select(-all_of(parent_col)) %>%
          left_join(lookup, by = child_col) %>%
          filter(!is.na(.data[[parent_col]])) %>%
          select(fullname, all_of(parent_col)),
        by = "fullname",
        unmatched = "ignore"
      )
  }
  message(" OK.")
}

# Fix unknown kingdoms and rank
taxonomy <- taxonomy %>%
  mutate(
    kingdom = case_when(
      domain %in% c("Bacteria", "Archaea") & phylum %like% "unknown" ~ "(unknown kingdom)",
      .default = kingdom
    ),
    rank = case_when(
      rank %in% c("(unknown rank)", "species group") ~ rank,
      subspecies != "" ~ "subspecies",
      species != ""    ~ "species",
      genus != ""      ~ "genus",
      family != ""     ~ "family",
      order != ""      ~ "order",
      class != ""      ~ "class",
      phylum != ""     ~ "phylum",
      kingdom != ""    ~ "kingdom",
      .default = NA_character_
    )
  )


# *** Save intermediate results (1c) *** ----------------------------------------------------------

saveRDS(taxonomy, "data-raw/taxonomy1c.rds")
# taxonomy <- readRDS("data-raw/taxonomy1c.rds")


# Add parent identifiers --------------------------------------------------------------------------
taxonomy$source[taxonomy$source == "Manually added"] <- "manually added"
taxonomy <- taxonomy %>%
  filter(status %in% c("accepted", "synonym") | source == "manually added")

# Add domain rows
taxonomy <- taxonomy %>%
  bind_rows(
    taxonomy %>%
      filter(
        kingdom == domain,
        rank == "kingdom",
        !domain %in% c("Bacteria", "Archaea"),
        fullname %unlike% "\\{kingdom\\}"
      ) %>%
      mutate(
        kingdom = "", rank = "domain", status = "accepted",
        source = "manually added", ref = NA_character_,
        across(matches("^(lpsn|mycobank|gbif)"), \(x) NA_character_)
      )
  )
taxonomy$rank[taxonomy$kingdom == taxonomy$domain &
                taxonomy$rank == "kingdom" &
                taxonomy$domain %in% c("Bacteria", "Archaea")] <- "domain"
taxonomy$kingdom[taxonomy$kingdom == taxonomy$domain &
                   taxonomy$rank == "domain"] <- ""

taxonomy <- taxonomy %>%
  arrange(fullname, lpsn, mycobank, desc(gbif)) %>%
  distinct(fullname, .keep_all = TRUE) %>%
  arrange(fullname)

# Resolve parent fullname for each record
taxonomy <- taxonomy %>%
  mutate(
    .parent_fullname = case_when(
      rank == "subspecies"              ~ paste(genus, species),
      rank == "species"                 ~ genus,
      rank == "genus"   & family != ""  ~ family,
      rank == "genus"   & order != ""   ~ order,
      rank == "genus"   & class != ""   ~ class,
      rank == "genus"   & phylum != ""  ~ phylum,
      rank == "genus"                   ~ domain,
      rank == "family"  & order != ""   ~ order,
      rank == "family"  & class != ""   ~ class,
      rank == "family"  & phylum != ""  ~ phylum,
      rank == "family"                  ~ domain,
      rank == "order"   & class != ""   ~ class,
      rank == "order"   & phylum != ""  ~ phylum,
      rank == "order"                   ~ domain,
      rank == "class"   & phylum != ""  ~ phylum,
      rank == "class"                   ~ domain,
      rank == "phylum"  & kingdom != "" ~ kingdom,
      rank == "phylum"                  ~ domain,
      rank == "kingdom"                 ~ domain,
      .default = NA_character_
    )
  )

# Build lookup vectors: fullname -> id
lpsn_lookup     <- setNames(taxonomy$lpsn, taxonomy$fullname)
mycobank_lookup <- setNames(taxonomy$mycobank, taxonomy$fullname)
gbif_lookup     <- setNames(taxonomy$gbif, taxonomy$fullname)

taxonomy <- taxonomy %>%
  mutate(
    lpsn_parent     = unname(lpsn_lookup[.parent_fullname]),
    mycobank_parent = unname(mycobank_lookup[.parent_fullname]),
    gbif_parent     = unname(gbif_lookup[.parent_fullname])
  ) %>%
  select(-.parent_fullname)

# Check: these still have no record in our data set
which(!taxonomy$lpsn_parent %in% taxonomy$lpsn)
which(!taxonomy$mycobank_parent %in% taxonomy$mycobank)
which(!taxonomy$gbif_parent %in% taxonomy$gbif)

# the domains now have rank = NA
taxonomy$fullname[is.na(taxonomy$rank)]
taxonomy$rank[is.na(taxonomy$rank)] <- "domain"


# Add prevalence ----------------------------------------------------------------------------------

# this part is required here, because it's needed for filtering on 'relevant' species to keep later on

pathogens <- read_excel(file_bartlett, sheet = "Tab 6 Full List")

# get all established, both old and current taxonomic names
established <- pathogens %>%
  filter(status == "established") %>%
  mutate(fullname = paste(genus, species)) %>%
  pull(fullname) %>%
  c(
    unlist(mo_current(.)),
    unlist(mo_synonyms(., keep_synonyms = FALSE))
  ) %>%
  strsplit(" ", fixed = TRUE) %>%
  sapply(function(x) if (length(x) == 1) x else paste(x[1], x[2])) %>%
  sort() %>%
  unique()

# get all putative, both old and current taxonomic names
putative <- pathogens %>%
  filter(status == "putative") %>%
  mutate(fullname = paste(genus, species)) %>%
  pull(fullname) %>%
  c(
    unlist(mo_current(.)),
    unlist(mo_synonyms(., keep_synonyms = FALSE))
  ) %>%
  strsplit(" ", fixed = TRUE) %>%
  sapply(function(x) if (length(x) == 1) x else paste(x[1], x[2])) %>%
  sort() %>%
  unique()

established <- established[established %unlike% "unknown"]
putative <- putative[putative %unlike% "unknown"]

established_genera <- established %>%
  strsplit(" ", fixed = TRUE) %>%
  sapply(function(x) x[1]) %>%
  sort() %>%
  unique()

putative_genera <- putative %>%
  strsplit(" ", fixed = TRUE) %>%
  sapply(function(x) x[1]) %>%
  sort() %>%
  unique()

nonbacterial_genera <- AMR:::MO_RELEVANT_GENERA %>%
  c(
    unlist(mo_current(.)),
    unlist(mo_synonyms(., keep_synonyms = FALSE))
  ) %>%
  strsplit(" ", fixed = TRUE) %>%
  sapply(function(x) x[1]) %>%
  sort() %>%
  unique()
nonbacterial_genera <- nonbacterial_genera[
  nonbacterial_genera %unlike% "unknown"
]

# update prevalence based on taxonomy (following the recent and thorough work of Bartlett et al., 2022)
# see https://doi.org/10.1099/mic.0.001269
taxonomy <- taxonomy %>%
  mutate(
    prevalence = case_when(
      # genera of the pathogens mentioned in the World Health Organization's (WHO) Priority Pathogen List
      genus %in% MO_WHO_PRIORITY_GENERA ~ 1.0,
      # 'established' means 'have infected at least three persons in three or more references'
      paste(genus, species) %in%
        established &
        rank %in% c("species", "subspecies") ~
        1.15,
      # other genera in the 'established' group
      genus %in% established_genera & rank == "genus" ~ 1.15,

      # 'putative' means 'fewer than three known cases'
      paste(genus, species) %in%
        putative &
        rank %in% c("species", "subspecies") ~
        1.25,
      
      # other genera in the 'putative' group
      genus %in% putative_genera & rank == "genus" ~ 1.25,
      # species and subspecies in 'established' and 'putative' groups
      genus %in%
        c(established_genera, putative_genera) &
        rank %in% c("species", "subspecies") ~
        1.5,
      
      # non-bacterial genera/species/subspecies of clinical relevance
      genus %in% nonbacterial_genera &
        domain != "Bacteria" &
        rank %in% c("genus", "species", "subspecies") ~
        1.25,
      # all others
      TRUE ~ 2.0
    )
  )

table(taxonomy$prevalence, useNA = "always")
# (a lot will be removed further below)


# *** Save intermediate results (2) *** -----------------------------------------------------------

saveRDS(taxonomy, "data-raw/taxonomy2.rds")
# taxonomy <- readRDS("data-raw/taxonomy2.rds")


# Remove unwanted taxonomic entries ---------------------------------------------------------------

part1 <- taxonomy %>%
  filter(
    # keep all we added ourselves:
    source == "manually added" |
      # keep all bacteria anyway, main focus of our package:
      domain == "Bacteria" |
      # and these kingdoms are very small, and also mainly microorganisms:
      domain == "Protozoa" |
      domain == "Archaea" |
      domain == "Chromista" |
      # keep everything from family up, 'ghost' entries will be removed later on
      rank %in% c("domain", "kingdom", "phylum", "class", "order", "family") |
      # other relevant genera to keep:
      genus %in% AMR:::MO_RELEVANT_GENERA |
      # relevant for biotechnology:
      genus %in%
        c(
          "Archaeoglobus",
          "Desulfurococcus",
          "Ferroglobus",
          "Ferroplasma",
          "Halobacterium",
          "Halococcus",
          "Haloferax",
          "Naegleria",
          "Nanoarchaeum",
          "Nitrosopumilus",
          "Nosema",
          "Pleistophora",
          "Pyrobaculum",
          "Pyrococcus",
          "Spraguea",
          "Thelohania",
          "Thermococcus",
          "Thermoplasma",
          "Thermoproteus",
          "Tritrichomonas"
        ) |
      genus %like_case% "Methano" |
      # domain of Protozoa:
      (phylum %in% c("Choanozoa", "Mycetozoa") & prevalence < 2) |
      # Fungi:
      (domain == "Fungi" &
        (!rank %in% c("genus", "species", "subspecies") |
          prevalence < 2 |
          class == "Pichiomycetes")) |
      # Animalia:
      genus %in% c("Lucilia", "Lumbricus") |
      (class == "Insecta" & !rank %in% c("species", "subspecies")) | # keep only genus of insects, not all of their (sub)species
      (genus == "Amoeba" & domain != "Animalia") # keep only in the protozoa, not the animalia
  ) %>%
  # this domain only contained Curvularia and Hymenolepis, which have coincidental twin names with Fungi
  filter(
    domain != "Plantae",
    !(genus %in% c("Aedes", "Anopheles") & rank %in% c("species", "subspecies"))
  )

# now get the parents and old names
part2 <- taxonomy %>%
  filter(
    gbif %in%
      c(
        part1$gbif_parent[!is.na(part1$gbif_parent)],
        part1$gbif_renamed_to[!is.na(part1$gbif_renamed_to)]
      ) |
      mycobank %in%
        c(
          part1$mycobank_parent[!is.na(part1$mycobank_parent)],
          part1$mycobank_renamed_to[!is.na(part1$mycobank_renamed_to)]
        ) |
      lpsn %in%
        c(
          part1$lpsn_parent[!is.na(part1$lpsn_parent)],
          part1$lpsn_renamed_to[!is.na(part1$lpsn_renamed_to)]
        )
  )
parts <- bind_rows(part1, part2)

part3 <- taxonomy %>%
  filter(
    gbif %in%
      c(
        parts$gbif_parent[!is.na(parts$gbif_parent)],
        parts$gbif_renamed_to[!is.na(parts$gbif_renamed_to)]
      ) |
      mycobank %in%
        c(
          parts$mycobank_parent[!is.na(parts$mycobank_parent)],
          parts$mycobank_renamed_to[!is.na(parts$mycobank_renamed_to)]
        ) |
      lpsn %in%
        c(
          parts$lpsn_parent[!is.na(parts$lpsn_parent)],
          parts$lpsn_renamed_to[!is.na(parts$lpsn_renamed_to)]
        )
  )
parts <- bind_rows(part1, part2, part3)

part4 <- taxonomy %>%
  filter(
    gbif %in%
      c(
        parts$gbif_parent[!is.na(parts$gbif_parent)],
        parts$gbif_renamed_to[!is.na(parts$gbif_renamed_to)]
      ) |
      mycobank %in%
        c(
          parts$mycobank_parent[!is.na(parts$mycobank_parent)],
          parts$mycobank_renamed_to[!is.na(parts$mycobank_renamed_to)]
        ) |
      lpsn %in%
        c(
          parts$lpsn_parent[!is.na(parts$lpsn_parent)],
          parts$lpsn_renamed_to[!is.na(parts$lpsn_renamed_to)]
        )
  )
parts <- bind_rows(part1, part2, part3, part4)


# Keep accepted names that any kept synonym points to
accepted_targets <- taxonomy %>%
  filter(
    gbif %in% na.omit(parts$gbif_renamed_to) |
      lpsn %in% na.omit(parts$lpsn_renamed_to) |
      mycobank %in% na.omit(parts$mycobank_renamed_to)
  )

taxonomy <- bind_rows(parts, accepted_targets) %>%
  mutate(
    # here we must prefer in this order: LPSN > MycoBank > GBIF
    source_index = case_when(
      source == "LPSN" ~ 1,
      source == "MycoBank" ~ 2,
      source == "GBIF" ~ 3,
      # manually added:
      TRUE ~ 4
    ),
    # also arrange on rank, otherwise e.g. Cyptococcus comes from Animalia, not Fungi
    domain_index = case_when(
      domain == "Bacteria" ~ 1,
      domain == "Fungi" ~ 2,
      domain == "Protozoa" ~ 3,
      domain == "Archaea" ~ 4,
      domain == "Chromista" ~ 5,
      domain == "Animalia" ~ 6,
      TRUE ~ 7
    )
  ) %>%
  arrange(source_index, domain_index, fullname) %>%
  distinct(fullname, .keep_all = TRUE) %>%
  select(-c(source_index, domain_index))

taxonomy <- taxonomy %>%
  mutate(
    lpsn_renamed_to = if_else(
      !is.na(lpsn_renamed_to) & !lpsn_renamed_to %in% lpsn,
      NA_character_, lpsn_renamed_to
    ),
    mycobank_renamed_to = if_else(
      !is.na(mycobank_renamed_to) & !mycobank_renamed_to %in% mycobank,
      NA_character_, mycobank_renamed_to
    ),
    gbif_renamed_to = if_else(
      !is.na(gbif_renamed_to) & !gbif_renamed_to %in% gbif,
      NA_character_, gbif_renamed_to
    )
  )

# Verify no unchased references remain
unchased <- taxonomy %>%
  filter(
    (!is.na(gbif_renamed_to) & !gbif_renamed_to %in% c(taxonomy$gbif, NA)) |
      (!is.na(lpsn_renamed_to) & !lpsn_renamed_to %in% c(taxonomy$lpsn, NA)) |
      (!is.na(mycobank_renamed_to) & !mycobank_renamed_to %in% c(taxonomy$mycobank, NA))
  )
if (nrow(unchased) > 0) {
  warning(nrow(unchased), " records have renamed_to references not present in taxonomy")
  table(unchased$prevalence)
}

# first make sure that species of a genus cannot be across multiple kingdoms or domains, since otherwise IDs cannot be given correctly
# (e.g., Amoeba has some species in Protozoa and some in Animalia)
# we acknowledge that this may be taxonomically right, but we need the same DOMAIN_GENUS identifier for each record
genus_ref <- taxonomy %>%
  filter(genus != "") %>%
  group_by(genus) %>%
  summarise(
    .ref_domain     = coalesce(first(domain[rank == "genus"]), first(domain)),
    .ref_kingdom    = coalesce(first(kingdom[rank == "genus"]), first(kingdom)),
    .ref_phylum     = coalesce(first(phylum[rank == "genus"]), first(phylum)),
    .ref_class      = coalesce(first(class[rank == "genus"]), first(class)),
    .ref_order      = coalesce(first(order[rank == "genus"]), first(order)),
    .ref_family     = coalesce(first(family[rank == "genus"]), first(family)),
    .ref_prevalence = coalesce(first(prevalence[rank == "genus"]), first(prevalence)),
    .groups = "drop"
  )

taxonomy <- taxonomy %>%
  left_join(genus_ref, by = "genus") %>%
  mutate(
    across_flag = rank != "genus" & genus != "" & fullname %unlike% "unknown",
    domain     = if_else(across_flag, .ref_domain, domain),
    kingdom    = if_else(across_flag, .ref_kingdom, kingdom),
    phylum     = if_else(across_flag, .ref_phylum, phylum),
    class      = if_else(across_flag, .ref_class, class),
    order      = if_else(across_flag, .ref_order, order),
    family     = if_else(across_flag, .ref_family, family),
    prevalence = if_else(across_flag, .ref_prevalence, prevalence)
  ) %>%
  select(-starts_with(".ref_"), -across_flag) %>%
  arrange(fullname)

# update the taxonomic names based on new genus classification
taxonomy <- taxonomy %>%
  # first, fix kingdom based on family level (only when family is not empty)
  group_by(domain, family) %>%
  mutate(
    kingdom = if_else(
      family != "",
      get_top_lvl(kingdom, rank, source, "phylum", phylum),
      kingdom
    ),
    phylum = if_else(
      family != "",
      first(phylum[phylum != "" & rank != "phylum"]),
      phylum
    ),
    class = if_else(
      family != "",
      first(class[class != "" & rank != "class"]),
      class
    ),
    order = if_else(
      family != "",
      first(order[order != "" & rank != "order"]),
      order
    )
  ) %>%
  # next, update all taxonomic layers
  group_by(domain, genus) %>%
  mutate(family = get_top_lvl(family, rank, source, "genus", genus)) %>%
  group_by(domain, family) %>%
  mutate(order = get_top_lvl(order, rank, source, "family", family)) %>%
  group_by(domain, order) %>%
  mutate(class = get_top_lvl(class, rank, source, "order", order)) %>%
  group_by(domain, class) %>%
  mutate(phylum = get_top_lvl(phylum, rank, source, "class", class)) %>%
  ungroup() %>%
  arrange(fullname) %>%
  mutate(across(kingdom:family, ~ if_else(is.na(.x), "", .x)))

# no ghost families, orders, classes, phyla
# (but keep the ghost families of bacteria)
taxonomy <- taxonomy %>%
  group_by(domain, family) %>%
  filter(
    n() > 1 |
      fullname %like% "unknown" |
      rank %in% c("domain", "genus", "species", "subspecies") |
      domain == "Bacteria"
  ) %>%
  group_by(domain, order) %>%
  filter(
    n() > 1 |
      fullname %like% "unknown" |
      rank %in% c("domain", "genus", "species", "subspecies") |
      family != ""
  ) %>%
  group_by(domain, class) %>%
  filter(
    n() > 1 |
      fullname %like% "unknown" |
      rank %in% c("domain", "genus", "species", "subspecies") |
      family != "" |
      order != ""
  ) %>%
  group_by(domain, phylum) %>%
  filter(
    n() > 1 |
      fullname %like% "unknown" |
      rank %in% c("domain", "genus", "species", "subspecies") |
      family != "" |
      order != "" |
      class != ""
  ) %>%
  ungroup()

# fix text in `domain` column where rank is "domain"
taxonomy$domain[taxonomy$rank == "domain"] <- gsub(" \\{domain\\}", "", taxonomy$fullname[taxonomy$rank == "domain"])


# *** Save intermediate results (2b) *** ----------------------------------------------------------

saveRDS(taxonomy, "data-raw/taxonomy2b.rds")
# taxonomy <- readRDS("data-raw/taxonomy2b.rds")


# Add Plasmodium species --------------------------------------------------

# MB 2026-05-09/ Plasmodium is not in the GBIF/COL dataset with its species, so we add it here manually.
plasmodia <- data.frame(
  domain = "Protozoa",
  kingdom = "Protozoa",
  phylum = "Apicomplexa",
  class = "Aconoidasida",
  order = "Haemospororida",
  family = "Plasmodiidae",
  genus = "Plasmodium",
  species = c("", "accipiteris", "achiotense", "achromaticum", "acuminatum", "adleri", "aegyptensis", "aeuminatum", "agamae", "alaudae", "alloelongatum", "anasum", "anomaluri", "arachniformis", "ashfordi", "atheruri", "attenuatum", "audaciosum", "auffenbergi", "aurulentum", "australis", "azurophilum", "balli", "bambusicolai", "basilisci", "beaucournui", "beebei", "beltrani", "berghei", "bertii", "bigueti", "billbrayi", "billcollinsi", "bioccai", "biziurae", "blacklocki", "booliati", "bouillize", "brasilianum", "brodeni", "brumpti", "brygooi", "bubalis", "bucki", "buteonis", "caloti", "capistrani", "caprae", "carmelinoi", "carteri", "cathemerium", "caucasica", "cephalophi", "cercopitheci", "chabaudi", "chiricahuae", "circularis", "circumflexum", "clelandi", "cnemaspi", "cnemidophori", "coatneyi", "coggeshalli", "colombiense", "columbae", "coluzzii", "cordyli", "coturnixi", "coulangesi", "cuculus", "cyclopsi", "cynomolgi", "cynomolgi", "cynomolgi", "cynomolgi", "delichoni", "dherteae", "diminutivum", "diploglossi", "dissanaikei", "dominicana", "dorsti", "draconis", "durae", "egerniae", "elongatum", "eylesi", "fairchildi", "falciparum", "fallax", "fieldi", "fischeri", "floridense", "foleyi", "formosanum", "forresteri", "fragile", "gabaldoni", "gaboni", "gallinaceum", "garnhami", "gemini", "georgesi", "ghadiriani", "giganteum", "ginsburgi", "giovannolai", "girardi", "globularis", "gloriai", "gologoense", "golvani", "gonatodi", "gonderi", "gracilis", "griffithsi", "guangdong", "gundersi", "guyannense", "hegneri", "heischi", "hermani", "heroni", "heteronucleare", "hexamerium", "hispaniolae", "hoionucleophilum", "holaspi", "holti", "homocircumflexum", "homopolare", "huffi", "hydrochaeri", "hylobati", "icipeensis", "iguanae", "incertae", "inopinatum", "intabazwe", "inui", "japonicum", "jeanriouxi", "jefferyi", "jiangi", "josephinae", "joyeuxi", "juxtanucleare", "kachelibaensis", "kadogoi", "kaninii", "kempi", "kentropyxi", "knowlesi", "knowlesi", "knowlesi", "koreafense", "kyaii", "lacertiliae", "lagopi", "lainsoni", "landauae", "lemuris", "lenoblei", "lepidoptiformis", "leucocytica", "lionatum", "lomamiensis", "lophurae", "loveridgei", "lucens", "lutzi", "lygosomae", "mabuiae", "mackerrasae", "mackiei", "maculilabre", "maior", "majus", "malagasi", "malariae", "marginatum", "matutinum", "megaglobularis", "megalotrypa", "melanipherum", "melanoleuca", "merulae", "mexicanum", "michikoa", "minasense", "minuoviride", "modestum", "mohammedi", "morulum", "multiformis", "multivacuolaris", "narayani", "necatrix", "neusticuri", "nucleophilium", "octamerium", "odhiamboi", "odocoilei", "ovale", "ovale", "ovale", "pachysomum", "paddae", "papernai", "parahexamerium", "paranucleophilum", "parvulum", "pedioecetii", "pelaezi", "percygarnhami", "pessoai", "petersi", "pifanoi", "pinotti", "pitheci", "pitmani", "polare", "polymorphum", "praefalciparum", "pulmophilium", "pythonias", "quelea", "reichenowi", "relictum", "reniai", "rhacodactyli", "rhadinurum", "rhodaini", "robinsoni", "rousetti", "rousseloti", "rouxi", "sandoshami", "sapaaensis", "sasai", "saurocaudatum", "scelopori", "schwetzi", "scorzai", "semiovale", "semnopitheci", "sergentorum", "silvaticum", "simium", "simplex", "smirnovi", "snounoui", "stellatum", "stuthionis", "tanzaniae", "tejerai", "telfordi", "tenue", "tomodoni", "torrealbai", "toucani", "traguli", "tranieri", "tribolonti", "tropiduri", "tumbayaensis", "tyrio", "uilenbergi", "uluguruense", "unalis", "uncinatum", "uzungwiense", "vacuolatum", "valkiunasi", "vastator", "vaughani", "vautieri", "venkataramiahii", "vinckei", "vivax", "volans", "voltaicum", "watteni", "wenyoni", "yoelii", "youngi", "zonuriae"),
  subspecies = c("", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "bastianelli", "ceylonensis", "cynomolgi", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "edesoni", "knowlesi", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "curtisi", "wallikeri", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", ""),
  source = "manually added",
  status = "accepted",
  prevalence = 1.25
) %>%
  mutate(rank = ifelse(species == "", "genus",
                       ifelse(subspecies != "", "subspecies",
                              "species")),
         fullname = trimws(paste(genus, species, subspecies)))

taxonomy <- taxonomy %>%
  filter(genus != "Plasmodium") %>%
  bind_rows(plasmodia) %>%
  arrange(fullname)


# Add microbial IDs -------------------------------------------------------------------------------

# (MO codes in the AMR package have the form: DOMAIN_GENUS_SPECIES_SUBSPECIES where all are abbreviated)

# Domain is abbreviated with 1 character, with exceptions for Animalia and Plantae
mo_domain <- taxonomy %>%
  filter(rank == "domain") %>%
  select(domain) %>%
  mutate(
    mo_domain = case_when(
      domain == "Animalia" ~ "AN",
      domain == "Archaea" ~ "A",
      domain == "Bacteria" ~ "B",
      domain == "Chromista" ~ "C",
      domain == "Fungi" ~ "F",
      domain == "Plantae" ~ "PL", # this is actually not part of the scope, having 0 results (2024)
      domain == "Protozoa" ~ "P",
      TRUE ~ ""
    )
  )
mo_domain

# in 2026 we added `domain`, remove this line later and just keep it `AMR::microorganisms`
existing_mo_tbl <- AMR::microorganisms %>% rename(domain = kingdom)


# phylum until family are abbreviated with 8 characters and prefixed with their rank

# Phylum - keep old and fill up for new ones
mo_phylum <- taxonomy %>%
  filter(rank == "phylum") %>%
  distinct(domain, phylum) %>%
  left_join(
    existing_mo_tbl %>%
      filter(rank == "phylum") %>%
      transmute(
        domain,
        phylum = fullname,
        mo_old = gsub("[A-Z]{1,2}_", "", as.character(mo))
      ),
    by = c("domain", "phylum")
  ) %>%
  group_by(domain) %>%
  mutate(
    mo_phylum8 = AMR:::abbreviate_mo(phylum, minlength = 8, prefix = "[PHL]_"),
    mo_phylum9 = AMR:::abbreviate_mo(phylum, minlength = 9, prefix = "[PHL]_"),
    mo_phylum = if_else(!is.na(mo_old), mo_old, mo_phylum8),
    mo_duplicated = duplicated(mo_phylum),
    mo_phylum = if_else(mo_duplicated, mo_phylum9, mo_phylum),
    mo_duplicated = duplicated(mo_phylum)
  ) %>%
  ungroup()
if (any(mo_phylum$mo_duplicated, na.rm = TRUE)) {
  stop("Duplicate MO codes for phylum!")
}
mo_phylum <- mo_phylum %>%
  select(domain, phylum, mo_phylum)

# Class - keep old and fill up for new ones
mo_class <- taxonomy %>%
  filter(rank == "class") %>%
  distinct(domain, class) %>%
  left_join(
    existing_mo_tbl %>%
      filter(rank == "class") %>%
      transmute(
        domain,
        class = fullname,
        mo_old = gsub("[A-Z]{1,2}_", "", as.character(mo))
      ),
    by = c("domain", "class")
  ) %>%
  group_by(domain) %>%
  mutate(
    mo_class8 = AMR:::abbreviate_mo(class, minlength = 8, prefix = "[CLS]_"),
    mo_class9 = AMR:::abbreviate_mo(class, minlength = 9, prefix = "[CLS]_"),
    mo_class = if_else(!is.na(mo_old), mo_old, mo_class8),
    mo_duplicated = duplicated(mo_class),
    mo_class = if_else(mo_duplicated, mo_class9, mo_class),
    mo_duplicated = duplicated(mo_class)
  ) %>%
  ungroup()
if (any(mo_class$mo_duplicated, na.rm = TRUE)) {
  stop("Duplicate MO codes for class!")
}
mo_class <- mo_class %>%
  select(domain, class, mo_class)

# Order - keep old and fill up for new ones
mo_order <- taxonomy %>%
  filter(rank == "order") %>%
  distinct(domain, order) %>%
  left_join(
    existing_mo_tbl %>%
      filter(rank == "order") %>%
      transmute(
        domain,
        order = fullname,
        mo_old = gsub("[A-Z]{1,2}_", "", as.character(mo))
      ),
    by = c("domain", "order")
  ) %>%
  group_by(domain) %>%
  mutate(
    mo_order8 = AMR:::abbreviate_mo(order, minlength = 8, prefix = "[ORD]_"),
    mo_order9 = AMR:::abbreviate_mo(order, minlength = 9, prefix = "[ORD]_"),
    mo_order = if_else(!is.na(mo_old), mo_old, mo_order8),
    mo_duplicated = duplicated(mo_order),
    mo_order = if_else(mo_duplicated, mo_order9, mo_order),
    mo_duplicated = duplicated(mo_order)
  ) %>%
  ungroup()
if (any(mo_order$mo_duplicated, na.rm = TRUE)) {
  stop("Duplicate MO codes for order!")
}
mo_order <- mo_order %>%
  select(domain, order, mo_order)

# Family - keep old and fill up for new ones
mo_family <- taxonomy %>%
  filter(rank == "family") %>%
  distinct(domain, family) %>%
  left_join(
    existing_mo_tbl %>%
      filter(rank == "family") %>%
      transmute(
        domain,
        family = fullname,
        mo_old = gsub("[A-Z]{1,2}_", "", as.character(mo))
      ),
    by = c("domain", "family")
  ) %>%
  group_by(domain) %>%
  mutate(
    mo_family8 = AMR:::abbreviate_mo(family, minlength = 8, prefix = "[FAM]_"),
    mo_family9 = AMR:::abbreviate_mo(family, minlength = 9, prefix = "[FAM]_"),
    mo_family = if_else(!is.na(mo_old), mo_old, mo_family8),
    mo_duplicated = duplicated(mo_family),
    mo_family = if_else(mo_duplicated, mo_family9, mo_family),
    mo_duplicated = duplicated(mo_family)
  ) %>%
  ungroup()
if (any(mo_family$mo_duplicated, na.rm = TRUE)) {
  stop("Duplicate MO codes for family!")
}
mo_family <- mo_family %>%
  select(domain, family, mo_family)

# construct code part for genus - keep old code where available and generate new ones where needed
mo_genus <- taxonomy %>%
  filter(rank == "genus") %>%
  distinct(domain, genus) %>%
  # get available old MO codes
  left_join(
    existing_mo_tbl %>%
      filter(rank == "genus") %>%
      transmute(
        mo_genus_old = gsub("^[A-Z]+_", "", as.character(mo)),
        domain,
        genus
      ) %>%
      distinct(domain, genus, .keep_all = TRUE),
    by = c("domain", "genus")
  ) %>%
  distinct(domain, genus, .keep_all = TRUE) %>%
  # since domain is part of the code, genus abbreviations may be duplicated between kingdoms
  group_by(domain) %>%
  # generate new MO codes for genus and set the right one
  mutate(
    mo_genus_new5 = AMR:::abbreviate_mo(genus, 5),
    mo_genus_new5b = paste0(AMR:::abbreviate_mo(genus, 5), 1),
    mo_genus_new5c = paste0(AMR:::abbreviate_mo(genus, 5), 2),
    mo_genus_new6 = AMR:::abbreviate_mo(genus, 6),
    mo_genus_new7 = AMR:::abbreviate_mo(genus, 7),
    mo_genus_new8 = AMR:::abbreviate_mo(genus, 8),
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
      mo_duplicated & mo_genus_new == mo_genus_old ~ mo_genus_new6,
      mo_duplicated & mo_genus_new == mo_genus_new5 ~ mo_genus_new6,
      mo_duplicated & mo_genus_new == mo_genus_new6 ~ mo_genus_new7,
      mo_duplicated & mo_genus_new == mo_genus_new7 ~ mo_genus_new8,
      mo_duplicated & mo_genus_new == mo_genus_new8 ~ mo_genus_new5b,
      TRUE ~ NA_character_
    ),
    mo_duplicated = duplicated(mo_genus_new),
    mo_genus_new = case_when(
      !mo_duplicated ~ mo_genus_new,
      mo_duplicated & mo_genus_new == mo_genus_new6 ~ mo_genus_new7,
      mo_duplicated & mo_genus_new == mo_genus_new7 ~ mo_genus_new8,
      mo_duplicated & mo_genus_new == mo_genus_new8 ~ mo_genus_new5b,
      mo_duplicated ~ mo_genus_new5c,
      TRUE ~ NA_character_
    ),
    mo_duplicated = duplicated(mo_genus_new)
  ) %>%
  ungroup()
if (any(mo_genus$mo_duplicated, na.rm = TRUE) | anyNA(mo_genus$mo_genus_new)) {
  stop("Duplicate MO codes for genus!")
}
# mo_genus %>% filter(mo_duplicated)
# no duplicates *within kingdoms*, so keep the right columns for left joining later
mo_genus <- mo_genus %>%
  select(domain, genus, mo_genus = mo_genus_new)

# same for species - keep old where available and create new per domain-genus where needed:
mo_species <- taxonomy %>%
  filter(rank == "species") %>%
  distinct(domain, genus, species) %>%
  left_join(
    existing_mo_tbl %>%
      filter(rank == "species") %>%
      transmute(
        mo_species_old = gsub("^[A-Z]+_[A-Z]+_", "", as.character(mo)),
        domain,
        genus,
        species
      ) %>%
      filter(mo_species_old %unlike% "-") %>%
      distinct(domain, genus, species, .keep_all = TRUE),
    by = c("domain", "genus", "species")
  ) %>%
  distinct(domain, genus, species, .keep_all = TRUE) %>%
  group_by(domain, genus) %>%
  mutate(
    mo_species_new4 = AMR:::abbreviate_mo(species, 4, hyphen_as_space = TRUE),
    mo_species_new5 = AMR:::abbreviate_mo(species, 5, hyphen_as_space = TRUE),
    mo_species_new5b = paste0(
      AMR:::abbreviate_mo(species, 5, hyphen_as_space = TRUE),
      1
    ),
    mo_species_new6 = AMR:::abbreviate_mo(species, 6, hyphen_as_space = TRUE),
    mo_species_new7 = AMR:::abbreviate_mo(species, 7, hyphen_as_space = TRUE),
    mo_species_new8 = AMR:::abbreviate_mo(species, 8, hyphen_as_space = TRUE),
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
if (
  any(mo_species$mo_duplicated, na.rm = TRUE) | anyNA(mo_species$mo_species_new)
) {
  stop("Duplicate MO codes for species!")
}
# no duplicates *within kingdoms*, so keep the right columns for left joining later
mo_species <- mo_species %>%
  select(domain, genus, species, mo_species = mo_species_new)

# same for subspecies - keep old where available and create new per domain-genus-species where needed:
mo_subspecies <- taxonomy %>%
  filter(rank == "subspecies") %>%
  distinct(domain, genus, species, subspecies) %>%
  left_join(
    existing_mo_tbl %>%
      filter(rank %in% c("subspecies", "subsp.", "infraspecies")) %>%
      transmute(
        mo_subspecies_old = gsub(
          "^[A-Z]+_[A-Z]+_[A-Z]+_",
          "",
          as.character(mo)
        ),
        domain,
        genus,
        species,
        subspecies
      ) %>%
      filter(mo_subspecies_old %unlike% "-") %>%
      distinct(domain, genus, species, subspecies, .keep_all = TRUE),
    by = c("domain", "genus", "species", "subspecies")
  ) %>%
  distinct(domain, genus, species, subspecies, .keep_all = TRUE) %>%
  group_by(domain, genus, species) %>%
  mutate(
    mo_subspecies_new4 = AMR:::abbreviate_mo(
      subspecies,
      4,
      hyphen_as_space = TRUE
    ),
    mo_subspecies_new5 = AMR:::abbreviate_mo(
      subspecies,
      5,
      hyphen_as_space = TRUE
    ),
    mo_subspecies_new5b = paste0(
      AMR:::abbreviate_mo(subspecies, 5, hyphen_as_space = TRUE),
      1
    ),
    mo_subspecies_new6 = AMR:::abbreviate_mo(
      subspecies,
      6,
      hyphen_as_space = TRUE
    ),
    mo_subspecies_new7 = AMR:::abbreviate_mo(
      subspecies,
      7,
      hyphen_as_space = TRUE
    ),
    mo_subspecies_new8 = AMR:::abbreviate_mo(
      subspecies,
      8,
      hyphen_as_space = TRUE
    ),
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
      mo_duplicated & mo_subspecies_new == mo_subspecies_new4 ~
        mo_subspecies_new5,
      mo_duplicated & mo_subspecies_new == mo_subspecies_new5 ~
        mo_subspecies_new6,
      mo_duplicated & mo_subspecies_new == mo_subspecies_new6 ~
        mo_subspecies_new7,
      mo_duplicated & mo_subspecies_new == mo_subspecies_new7 ~
        mo_subspecies_new8,
      mo_duplicated & mo_subspecies_new == mo_subspecies_new8 ~
        mo_subspecies_new5b,
      TRUE ~ NA_character_
    ),
    mo_duplicated = duplicated(mo_subspecies_new)
  ) %>%
  ungroup()
if (
  any(mo_subspecies$mo_duplicated, na.rm = TRUE) |
    anyNA(mo_subspecies$mo_subspecies_new)
) {
  stop("Duplicate MO codes for subspecies!")
}
# no duplicates *within kingdoms*, so keep the right columns for left joining later
mo_subspecies <- mo_subspecies %>%
  select(domain, genus, species, subspecies, mo_subspecies = mo_subspecies_new)

# unknowns - manually added
mo_unknown <- existing_mo_tbl %>%
  filter(fullname %like% "unknown") %>%
  transmute(fullname, mo_unknown = as.character(mo))

# apply the new codes!
taxonomy1 <- taxonomy
taxonomy <- taxonomy %>%
  left_join(mo_domain, by = "domain") %>%
  left_join(mo_phylum, by = c("domain", "phylum")) %>%
  left_join(mo_class, by = c("domain", "class")) %>%
  left_join(mo_order, by = c("domain", "order")) %>%
  left_join(mo_family, by = c("domain", "family")) %>%
  left_join(mo_genus, by = c("domain", "genus")) %>%
  left_join(mo_species, by = c("domain", "genus", "species")) %>%
  left_join(
    mo_subspecies,
    by = c("domain", "genus", "species", "subspecies")
  ) %>%
  left_join(mo_unknown, by = "fullname") %>%
  mutate(across(starts_with("mo_"), function(x) if_else(is.na(x), "", x))) %>%
  mutate(
    mo = case_when(
      fullname %like% "unknown" ~ mo_unknown,
      # add special cases for taxons higher than genus
      rank == "domain" ~
        paste(mo_domain, "[DMN]", toupper(domain), sep = "_"),
      rank == "kingdom" ~
        paste(mo_domain, "[KNG]", toupper(kingdom), sep = "_"),
      rank == "phylum" ~ paste(mo_domain, mo_phylum, sep = "_"),
      rank == "class" ~ paste(mo_domain, mo_class, sep = "_"),
      rank == "order" ~ paste(mo_domain, mo_order, sep = "_"),
      rank == "family" ~ paste(mo_domain, mo_family, sep = "_"),
      TRUE ~ paste(mo_domain, mo_genus, mo_species, mo_subspecies, sep = "_")
    ),
    mo = trimws(gsub("_+$", "", mo)),
    .before = 1
  ) %>%
  select(!starts_with("mo_")) %>%
  arrange(fullname)

# now check these duplicate fullnames, they should not exist (didn't at least in 2024)
# this example shows that some could have same fullnames: taxonomy %>% filter(class == genus & class != "")
# such as Kapibacteria (both genus and class), but the class last time correctly got the " {class}" suffix in the fullname
taxonomy %>%
  filter(fullname %in% .[duplicated(fullname), "fullname", drop = TRUE]) %>%
  View()

# OTHERWISE uncomment this section:
# # we will prefer domains as Bacteria > Fungi > Protozoa > Archaea > Animalia, so check the ones within 1 domain
# taxonomy %>%
#   filter(fullname %in% .[duplicated(fullname), "fullname", drop = TRUE]) %>%
#   group_by(fullname) %>%
#   filter(n_distinct(domain) == 1)
#
# # fullnames must be unique, we'll keep the most relevant ones only
# taxonomy <- taxonomy %>%
#   mutate(rank_index = case_when(
#     domain == "Bacteria" ~ 1,
#     domain == "Fungi" ~ 2,
#     domain == "Protozoa" ~ 3,
#     domain == "Archaea" ~ 4,
#     domain == "Chromista" ~ 5,
#     domain == "Animalia" ~ 6,
#     TRUE ~ 7
#   )) %>%
#   arrange(fullname, rank_index) %>%
#   distinct(fullname, .keep_all = TRUE) %>%
#   select(-rank_index) %>%
#   filter(mo != "")

# keep the codes from manually added ones
manual_mos <- as.character(existing_mo_tbl$mo)[match(
  taxonomy$fullname[taxonomy$source == "manually added"],
  existing_mo_tbl$fullname
)]
taxonomy$mo[taxonomy$source == "manually added"][!is.na(manual_mos)] <- manual_mos[!is.na(manual_mos)]

# remove the manually added ones that are now ghosts
# MB 2026-05-09/ not needed anymore; 0 rows
# taxonomy %>%
#   filter((mo %like% "__" & source == "manually added")) %>%
#   View()
# taxonomy <- taxonomy %>%
#   filter(!(mo %like% "__" & source == "manually added"))

# this must not exist (2024: GBIF mess with missing genera - remove them):
taxonomy %>%
  filter(mo %like% "__") %>%
  View()
taxonomy <- taxonomy %>% filter(mo %unlike% "__")

# keep unique taxonomy
taxonomy <- bind_rows(
  taxonomy %>% 
    filter(fullname %like% "unknown"),
  taxonomy %>%
    distinct(across(domain:subspecies), .keep_all = TRUE) %>%
    filter(fullname %unlike% "unknown")
)

saveRDS(taxonomy, "data-raw/taxonomy2c.rds")
# taxonomy <- readRDS("data-raw/taxonomy2c.rds")

# Some integrity checks ---------------------------------------------------------------------------

# are mo codes unique?
taxonomy %>%
  filter(mo %in% .[duplicated(mo), "mo", drop = TRUE]) %>%
  mutate(removed_if_distinct_is_applied = duplicated(mo), .before = 1) %>%
  View()
# checked 2026, can all be removed, they are actual duplicates
taxonomy <- taxonomy %>%
  distinct(mo, .keep_all = TRUE)
# Also:
# - There are multiple cases where the suffix " {species}" was added - these are just duplicates.
#   Some are Salmonella errors at LPSN's side - typhimurium is NOT a species of it, species is enterica, serovar is Typhimurium (that we use as subspecies)
# taxonomy %>%
#   filter(mo %in% .[duplicated(mo), "mo", drop = TRUE]) %>%
#   arrange(mo) %>%
#   View()
# keep the firsts
# taxonomy <- taxonomy %>%
#   arrange(mo) %>%
#   distinct(mo, .keep_all = TRUE)

# are fullnames unique?
taxonomy %>%
  arrange(fullname) %>%
  filter(fullname %in% .[duplicated(fullname), "fullname", drop = TRUE]) %>%
  View()
# # we can now solve it like this, but check every year if this is okay
# taxonomy <- taxonomy %>%
#   group_by(fullname) %>%
#   mutate(fullname = if_else(duplicated(fullname), paste0(fullname, " {", rank, "}"), fullname)) %>%
#   ungroup()

# are all GBIFs available?
taxonomy %>%
  filter(
    (!gbif_parent %in% gbif) |
      (!lpsn_parent %in% lpsn) |
      (!mycobank_parent %in% mycobank)
  ) %>%
  count(
    source = case_when(
      !gbif_parent %in% gbif ~ "GBIF",
      !lpsn_parent %in% lpsn ~ "LPSN",
      !mycobank_parent %in% mycobank ~ "MycoBank",
      TRUE ~ "?"
    ),
    rank
  )

# so fix again all parent identifiers (exact same code as earlier)
taxonomy <- taxonomy %>%
  mutate(
    lpsn_parent = case_when(
      rank == "kingdom" ~ lpsn[match(domain, fullname)],
      # in phylum, take parent from kingdom if available, otherwise domain
      rank == "phylum" & kingdom != "" ~ lpsn[match(kingdom, fullname)],
      # in class, take parent from phylum if available, otherwise domain
      rank == "class" & phylum != "" ~ lpsn[match(phylum, fullname)],
      rank == "class" ~ lpsn[match(domain, fullname)],
      # in order, take parent from class if available, otherwise phylum, otherwise domain
      rank == "order" & class != "" ~ lpsn[match(class, fullname)],
      rank == "order" & phylum != "" ~ lpsn[match(phylum, fullname)],
      rank == "order" ~ lpsn[match(domain, fullname)],
      # family
      rank == "family" & order != "" ~ lpsn[match(order, fullname)],
      rank == "family" & class != "" ~ lpsn[match(class, fullname)],
      rank == "family" & phylum != "" ~ lpsn[match(phylum, fullname)],
      rank == "family" ~ lpsn[match(domain, fullname)],
      # genus
      rank == "genus" & family != "" ~ lpsn[match(family, fullname)],
      rank == "genus" & order != "" ~ lpsn[match(order, fullname)],
      rank == "genus" & class != "" ~ lpsn[match(class, fullname)],
      rank == "genus" & phylum != "" ~ lpsn[match(phylum, fullname)],
      rank == "genus" ~ lpsn[match(domain, fullname)],
      # species, always has a genus
      rank == "species" ~ lpsn[match(genus, fullname)],
      # subspecies, always has a genus + species
      rank == "subspecies" ~ lpsn[match(paste(genus, species), fullname)],
      TRUE ~ NA_character_
    ),
    mycobank_parent = case_when(
      rank == "kingdom" ~ mycobank[match(domain, fullname)],
      # in phylum, take parent from kingdom if available, otherwise domain
      rank == "phylum" & kingdom != "" ~ mycobank[match(kingdom, fullname)],
      # class
      rank == "class" & phylum != "" ~ mycobank[match(phylum, fullname)],
      rank == "class" ~ mycobank[match(domain, fullname)],
      # order
      rank == "order" & class != "" ~ mycobank[match(class, fullname)],
      rank == "order" & phylum != "" ~ mycobank[match(phylum, fullname)],
      rank == "order" ~ mycobank[match(domain, fullname)],
      # family
      rank == "family" & order != "" ~ mycobank[match(order, fullname)],
      rank == "family" & class != "" ~ mycobank[match(class, fullname)],
      rank == "family" & phylum != "" ~ mycobank[match(phylum, fullname)],
      rank == "family" ~ mycobank[match(domain, fullname)],
      # genus
      rank == "genus" & family != "" ~ mycobank[match(family, fullname)],
      rank == "genus" & order != "" ~ mycobank[match(order, fullname)],
      rank == "genus" & class != "" ~ mycobank[match(class, fullname)],
      rank == "genus" & phylum != "" ~ mycobank[match(phylum, fullname)],
      rank == "genus" ~ mycobank[match(domain, fullname)],
      # species
      rank == "species" ~ mycobank[match(genus, fullname)],
      # subspecies
      rank == "subspecies" ~ mycobank[match(paste(genus, species), fullname)],
      TRUE ~ NA_character_
    ),
    gbif_parent = case_when(
      rank == "kingdom" ~ gbif[match(domain, fullname)],
      # in phylum, take parent from kingdom if available, otherwise domain
      rank == "phylum" & kingdom != "" ~ gbif[match(kingdom, fullname)],
      # class
      rank == "class" & phylum != "" ~ gbif[match(phylum, fullname)],
      rank == "class" ~ gbif[match(domain, fullname)],
      # order
      rank == "order" & class != "" ~ gbif[match(class, fullname)],
      rank == "order" & phylum != "" ~ gbif[match(phylum, fullname)],
      rank == "order" ~ gbif[match(domain, fullname)],
      # family
      rank == "family" & order != "" ~ gbif[match(order, fullname)],
      rank == "family" & class != "" ~ gbif[match(class, fullname)],
      rank == "family" & phylum != "" ~ gbif[match(phylum, fullname)],
      rank == "family" ~ gbif[match(domain, fullname)],
      # genus
      rank == "genus" & family != "" ~ gbif[match(family, fullname)],
      rank == "genus" & order != "" ~ gbif[match(order, fullname)],
      rank == "genus" & class != "" ~ gbif[match(class, fullname)],
      rank == "genus" & phylum != "" ~ gbif[match(phylum, fullname)],
      rank == "genus" ~ gbif[match(domain, fullname)],
      # species
      rank == "species" ~ gbif[match(genus, fullname)],
      # subspecies
      rank == "subspecies" ~ gbif[match(paste(genus, species), fullname)],
      TRUE ~ NA_character_
    )
  )

# check again
taxonomy %>%
  filter(
    (!gbif_parent %in% gbif) |
      (!lpsn_parent %in% lpsn) |
      (!mycobank_parent %in% mycobank)
  ) %>%
  count(
    source = case_when(
      !gbif_parent %in% gbif ~ "GBIF",
      !lpsn_parent %in% lpsn ~ "LPSN",
      !mycobank_parent %in% mycobank ~ "MycoBank",
      TRUE ~ "?"
    ),
    rank
  )


# *** Save intermediate results (3) *** -----------------------------------------------------------

saveRDS(taxonomy, "data-raw/taxonomy3.rds")
# taxonomy <- readRDS("data-raw/taxonomy3.rds")

# Finalised taxonomy (without additional data) ----------------------------------------------------

message(
  "\nCongratulations! The new taxonomic table will contain ",
  format(nrow(taxonomy), big.mark = " "),
  " rows.\n",
  "This was ",
  format(nrow(AMR::microorganisms), big.mark = " "),
  " rows.\n"
)

# these are the new ones:
taxonomy %>%
  # filter(!paste(kingdom, fullname) %in% paste(AMR::microorganisms$kingdom, AMR::microorganisms$fullname)) %>%
  filter(!fullname %in% AMR::microorganisms$fullname) %>%
  View()
# these are to be removed:
AMR::microorganisms %>%
  filter(!fullname %in% taxonomy$fullname) %>%
  View()
# be extra critical to these:
AMR::microorganisms %>%
  filter(!fullname %in% taxonomy$fullname, source == "LPSN", status == "accepted") %>%
  count(rank)

# Add SNOMED CT -----------------------------------------------------------------------------------

# we will use Public Health Information Network Vocabulary Access and Distribution System (PHIN VADS)
# as a source, which copies directly from the latest US SNOMED CT version
# - go to https://phinvads.cdc.gov/vads/ViewValueSet.action?oid=2.16.840.1.114222.4.11.1009
# - check that current online version is higher than TAXONOMY_VERSION$SNOMED
# - if so, click on 'Download Value Set', choose 'TXT'
snomed <- vroom("data-raw/SNOMED_PHVS_Microorganism_CDC_V12.txt", skip = 3, guess_max = 1e5) %>%
  select(1:2) %>%
  setNames(c("snomed", "mo")) %>%
  mutate(snomed = as.character(snomed))

# try to get name of MO
snomed <- snomed %>%
  mutate(mo = gsub("(ss[.]|subspecies) ", "", mo)) %>%
  mutate(
    fullname = case_when(
      mo %like_case% "[A-Z][a-z]+ [a-z]+ [a-z]{4,} " ~
        gsub("(^|.*)([A-Z][a-z]+ [a-z]+ [a-z]{4,}) .*", "\\2", mo),
      mo %like_case% "[A-Z][a-z]+ [a-z]{4,} " ~
        gsub("(^|.*)([A-Z][a-z]+ [a-z]{4,}) .*", "\\2", mo),
      mo %like_case% "[A-Z][a-z]+" ~
        gsub("(^|.*)([A-Z][a-z]+)( .*|$)", "\\2", mo),
      TRUE ~ NA_character_
    )
  )
snomed <- snomed %>%
  filter(fullname %in% taxonomy$fullname)

message(
  nrow(snomed),
  " SNOMED codes will be added to ",
  n_distinct(snomed$fullname),
  " microorganisms"
)

snomed <- snomed %>%
  group_by(fullname) %>%
  summarise(snomed = list(snomed))

taxonomy <- taxonomy %>%
  left_join(snomed, by = "fullname")


# Add oxygen tolerance (aerobe/anaerobe) ----------------------------------------------------------

bacdive <- vroom::vroom(file_oxygen_tolerance, skip = 2) %>%
  select(species, oxygen = `Oxygen tolerance`)
bacdive <- bacdive %>%
  # fill in missing species from previous rows
  mutate(fullname = if_else(is.na(species), lag(species), species)) %>%
  filter(
    !is.na(species),
    !is.na(oxygen),
    oxygen %unlike% "tolerant",
    species %unlike% "unclassified"
  ) %>%
  select(-species)
bacdive <- bacdive %>%
  # now determine type per species
  group_by(fullname) %>%
  summarise(
    fullname = first(fullname),
    oxygen_tolerance = case_when(
      any(oxygen %like% "facultative") ~ "facultative anaerobe",
      all(oxygen == "microaerophile") ~ "microaerophile",
      all(oxygen %in% c("anaerobe", "obligate anaerobe")) ~ "anaerobe",
      all(oxygen %in% c("anaerobe", "obligate anaerobe", "microaerophile")) ~
        "anaerobe/microaerophile",
      all(oxygen %in% c("aerobe", "obligate aerobe")) ~ "aerobe",
      all(!oxygen %in% c("anaerobe", "obligate anaerobe")) ~ "aerobe",
      all(c("aerobe", "anaerobe") %in% oxygen) ~ "facultative anaerobe",
      TRUE ~ NA_character_
    )
  )
# now find all synonyms and copy them from their current taxonomic names
synonyms <- taxonomy %>%
  filter(status == "synonym") %>%
  transmute(
    mo,
    fullname_old = fullname,
    current = synonym_mo_to_accepted_mo(
      mo,
      fill_in_accepted = FALSE,
      dataset = taxonomy
    )
  ) %>%
  filter(!is.na(current)) %>%
  mutate(fullname = taxonomy$fullname[match(current, taxonomy$mo)]) %>%
  left_join(bacdive, by = "fullname") %>%
  filter(!is.na(oxygen_tolerance)) %>%
  select(fullname, oxygen_tolerance)

bacdive <- bacdive %>%
  bind_rows(synonyms) %>%
  distinct()

bacdive_genus <- bacdive %>%
  mutate(
    oxygen = oxygen_tolerance,
    genus = taxonomy$genus[match(fullname, taxonomy$fullname)]
  ) %>%
  group_by(fullname = genus) %>%
  summarise(
    oxygen_tolerance = case_when(
      any(oxygen == "facultative anaerobe") ~ "facultative anaerobe",
      any(oxygen == "anaerobe/microaerophile") ~ "anaerobe/microaerophile",
      all(oxygen == "microaerophile") ~ "microaerophile",
      all(oxygen == "anaerobe") ~ "anaerobe",
      all(oxygen == "aerobe") ~ "aerobe",
      TRUE ~ "facultative anaerobe"
    )
  )
bacdive <- bacdive %>%
  bind_rows(bacdive_genus) %>%
  arrange(fullname)

bacdive_other <- taxonomy %>%
  filter(
    domain == "Bacteria",
    rank == "species",
    !fullname %in% bacdive$fullname,
    genus %in% bacdive$fullname
  ) %>%
  select(fullname, genus) %>%
  left_join(bacdive, by = c("genus" = "fullname")) %>%
  mutate(
    oxygen_tolerance = if_else(
      oxygen_tolerance %in%
        c("aerobe", "anaerobe", "microaerophile", "anaerobe/microaerophile"),
      oxygen_tolerance,
      paste("likely", oxygen_tolerance)
    )
  ) %>%
  select(fullname, oxygen_tolerance) %>%
  distinct(fullname, .keep_all = TRUE)

bacdive <- bacdive %>%
  bind_rows(bacdive_other) %>%
  arrange(fullname) %>%
  distinct(fullname, .keep_all = TRUE)

taxonomy <- taxonomy %>%
  left_join(bacdive, by = "fullname") %>%
  relocate(oxygen_tolerance, .after = ref)

taxonomy %>% count(oxygen_tolerance)


# Add morphology ---------------------------------------------------------------------

bacdive_shape <- vroom::vroom(file_cell_shape, skip = 2, guess_max = 1e5) %>%
  select(species, shape = `Cell shape`)
bacdive_shape <- bacdive_shape %>%
  # fill in missing species from previous rows
  mutate(fullname = if_else(is.na(species), lag(species), species)) %>%
  filter(
    !is.na(species),
    !is.na(shape),
    species %unlike% "unclassified"
  ) %>%
  select(-species)
bacdive_shape <- bacdive_shape %>%
  # map raw BacDive values to a controlled vocabulary
  mutate(
    shape = case_when(
      shape %in% c("coccus-shaped", "sphere-shaped", "diplococcus-shaped") ~ "cocci",
      shape %in% c("oval-shaped", "ovoid-shaped") ~ "coccobacilli",
      shape %in% c("rod-shaped", "curved-shaped", "vibrio-shaped", "flask-shaped") ~ "rods",
      shape %in% c("spiral-shaped", "helical-shaped") ~ "spirilla",
      shape == "filament-shaped" ~ "filamentous",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(shape)) %>%
  # now determine shape per species by majority vote
  group_by(fullname) %>%
  summarise(
    morphology = names(sort(table(shape), decreasing = TRUE))[1]
  )
# now find all synonyms and copy them from their current taxonomic names
synonyms_shape <- taxonomy %>%
  filter(status == "synonym") %>%
  transmute(
    mo,
    fullname_old = fullname,
    current = synonym_mo_to_accepted_mo(
      mo,
      fill_in_accepted = FALSE,
      dataset = taxonomy
    )
  ) %>%
  filter(!is.na(current)) %>%
  mutate(fullname = taxonomy$fullname[match(current, taxonomy$mo)]) %>%
  left_join(bacdive_shape, by = "fullname") %>%
  filter(!is.na(morphology)) %>%
  select(fullname, morphology)

bacdive_shape <- bacdive_shape %>%
  bind_rows(synonyms_shape) %>%
  distinct()

bacdive_shape_genus <- bacdive_shape %>%
  mutate(
    shape_raw = morphology,
    genus = taxonomy$genus[match(fullname, taxonomy$fullname)]
  ) %>%
  group_by(fullname = genus) %>%
  summarise(
    morphology = names(sort(table(shape_raw), decreasing = TRUE))[1]
  )
bacdive_shape <- bacdive_shape %>%
  bind_rows(bacdive_shape_genus) %>%
  arrange(fullname)

bacdive_shape_other <- taxonomy %>%
  filter(
    domain == "Bacteria",
    rank == "species",
    !fullname %in% bacdive_shape$fullname,
    genus %in% bacdive_shape$fullname
  ) %>%
  select(fullname, genus) %>%
  left_join(bacdive_shape, by = c("genus" = "fullname")) %>%
  mutate(
    morphology = case_when(
      fullname %like% "coccus" ~ "cocci",
      fullname %in% taxonomy$fullname[taxonomy$order %in% c("Enterobacterales", "Caryophanales", "Lactobacillales")] ~ morphology,
      TRUE ~ paste("likely", morphology))
  ) %>%
  select(fullname, morphology) %>%
  distinct(fullname, .keep_all = TRUE)

bacdive_shape <- bacdive_shape %>%
  bind_rows(bacdive_shape_other) %>%
  arrange(fullname) %>%
  distinct(fullname, .keep_all = TRUE)

taxonomy <- taxonomy %>%
  left_join(bacdive_shape, by = "fullname") %>%
  relocate(morphology, .after = oxygen_tolerance)

# Override: genera that are clinically established coccobacilli but where BacDive
# majority vote yields "rods" due to observer disagreement on the rod/oval boundary.
# These genera are universally reported as coccobacilli on Gram stain in clinical
# microbiology practice.
coccobacilli_genera <- c(
  "Acinetobacter", "Aggregatibacter", "Brucella",
  "Gardnerella", "Haemophilus", "Kingella",
  "Moraxella", "Pasteurella"
)
taxonomy <- taxonomy %>%
  mutate(
    morphology = case_when(
      genus %in% coccobacilli_genera & is.na(morphology) ~ "likely coccobacilli",
      genus %in% coccobacilli_genera &
        morphology %in% c("rods", "cocci") ~ "coccobacilli",
      genus %in% coccobacilli_genera &
        morphology %in% c("likely rods", "likely cocci") ~ "likely coccobacilli",
      TRUE ~ morphology
    )
  )

# Spirochaetes: the entire phylum is spirochaete by definition, fill in where missing
taxonomy <- taxonomy %>%
  mutate(
    morphology = case_when(
      phylum %in% c("Spirochaetota", "Spirochaetes") & is.na(morphology) ~ "likely spirilla",
      phylum %in% c("Spirochaetota", "Spirochaetes") &
        morphology %in% c("rods", "likely rods") ~ "spirilla",
      TRUE ~ morphology
    )
  )

taxonomy %>% count(morphology)


# Restore 'synonym' microorganisms to 'accepted' --------------------------------------------------

# If there are some synonyms that need to be corrected to 'accepted', you can do that here.
# Before 2024, we encountered this (but currently, this is all good, no action needed):

# according to LPSN: Stenotrophomonas maltophilia is the correct name if this species is regarded as a separate species (i.e., if its nomenclatural type is not assigned to another species whose name is validly published, legitimate and not rejected and has priority) within a separate genus Stenotrophomonas.
# https://lpsn.dsmz.de/species/stenotrophomonas-maltophilia

"Moraxella catarrhalis" %in% taxonomy$fullname
"Stenotrophomonas maltophilia" %in% taxonomy$fullname
# # all MO's to keep as 'accepted', not as 'synonym':
# to_restore <- c(
#   "Stenotrophomonas maltophilia",
#   "Moraxella catarrhalis"
# )
# taxonomy %>% filter(fullname %in% to_restore) %>% View()
# all(to_restore %in% taxonomy$fullname)
# for (nm in to_restore) {
#   taxonomy$lpsn_renamed_to[which(taxonomy$fullname == nm)] <- NA
#   taxonomy$gbif_renamed_to[which(taxonomy$fullname == nm)] <- NA
#   taxonomy$status[which(taxonomy$fullname == nm)] <- "accepted"
# }


# Add species groups ------------------------------------------------------------------------------

# just before the end, make sure we get the species group from the previous data set
old_groups <- AMR::microorganisms %>%
  filter(rank == "species group") %>%
  mutate(mo = as.character(mo))

# use current taxonomy
groups <- taxonomy[match(old_groups$genus, taxonomy$genus), ]
groups <- groups %>%
  mutate(
    mo = as.character(old_groups$mo),
    fullname = old_groups$fullname,
    species = old_groups$species,
    subspecies = old_groups$subspecies,
    rank = old_groups$rank,
    ref = old_groups$ref,
    status = "accepted",
    source = "manually added",
    lpsn = NA_character_,
    lpsn_renamed_to = NA_character_,
    mycobank = NA_character_,
    mycobank_renamed_to = NA_character_,
    gbif = NA_character_,
    gbif_renamed_to = NA_character_
  ) %>%
  select(
    -c(
      lpsn,
      lpsn_renamed_to,
      mycobank,
      mycobank_renamed_to,
      gbif,
      gbif_renamed_to,
      snomed
    )
  )

taxonomy <- taxonomy %>%
  filter(!mo %in% groups$mo) %>%
  bind_rows(groups)

# we added an MO code, so make sure everything is still unique
any(duplicated(taxonomy$mo))
taxonomy$mo[duplicated(taxonomy$mo)]
any(duplicated(taxonomy$fullname))
taxonomy$fullname[duplicated(taxonomy$fullname)]


# Set unknown ranks -------------------------------------------------------------------------------

# MB 2026-05-18/ not needed anymore, can be removed next year
taxonomy$rank[which(taxonomy$fullname %like% "unknown")] <- "(unknown rank)"


# Fix repetitive elements in MO code --------------------------------------------------------------

# this happened in early 2025 (and May 2026), check that MO codes do not have repeated elements
# fixed it then like this: microorganisms$mo <- gsub("B_SCLLM_CNNM_LNSM_LNSM_LNSM_LNSM", "B_SCLLM_CNNM", microorganisms$mo)
taxonomy %>%
  filter(mo %like% "_.*_.*_.*_") %>%
  View()
# Pattern: capture the prefix, then remove its repetition before the final epithet
taxonomy$mo[taxonomy$mo %like% "_.*_.*_.*_"] <- gsub("^(B_[A-Z0-9]+(?:_[A-Z0-9]+)?)_\\1_([A-Z]+)$",
                                                     "\\1_\\2",
                                                     taxonomy$mo[taxonomy$mo %like% "_.*_.*_.*_"],
                                                     perl = TRUE)
taxonomy %>%
  filter(mo %like% "_.*_.*_.*_") %>%
  View()


# Remove childless non-bacterial genera -----------------------------------------------------------

# this removes all genera that have no species, except for the domain of Bacteria.
taxonomy <- taxonomy %>%
  filter(
    rank != "genus" |
      domain == "Bacteria" |
      genus %in% taxonomy$genus[taxonomy$rank == "species"]
  )

# then remove all childless upper taxonomy caused by this
for (rank_name in c("family", "order", "class", "phylum", "kingdom")) {
  has_children <- taxonomy %>%
    filter(rank != rank_name, .data[[rank_name]] != "") %>%
    distinct(.data[[rank_name]]) %>%
    pull()
  
  n_before <- nrow(taxonomy)
  taxonomy <- taxonomy %>%
    filter(!(rank == rank_name & !fullname %in% has_children))
  message("Removed ", n_before - nrow(taxonomy), " childless ", rank_name, " entries")
}

# Clean-up of nonsense names ----------------------------------------------------------------------

table(unlist(strsplit(microorganisms$fullname, "")))
microorganisms <- microorganisms %>% 
  filter(fullname %unlike% " [0-9]") %>% 
  filter(fullname %unlike% " (i|ii|iii|iv|v|vi|vii|viii|ix|x)$" | genus == "Streptococcus") %>% 
  filter(fullname %unlike% "[.]") %>% 
  mutate(across(c(fullname, domain:subspecies), function(x) gsub("[\"'`]", "", x)))



# Some final checks -------------------------------------------------------------------------------

# all our previously manually added names should be in it
all(
  microorganisms$fullname[microorganisms$source == "manually added"] %in%
    taxonomy$fullname
)
microorganisms$fullname[
  !microorganisms$fullname[microorganisms$source == "manually added"] %in%
    taxonomy$fullname
]

# we added an MO code, so make sure everything is still unique
any(duplicated(taxonomy$mo))
any(duplicated(taxonomy$fullname))


# Update other data sets --------------------------------------------------------------------------

fix_old_mos <- function(dataset, new_ref, drop = FALSE, col = "mo") {
  before <- dataset
  
  mo_names <- AMR::microorganisms$fullname[match(before[[col]], AMR::microorganisms$mo)]
  matches <- new_ref$mo[match(mo_names, new_ref$fullname)]
  after <- before
  after[[col]] <- matches
  class(after[[col]]) <- c("mo", "character")
  if (drop == TRUE) {
    after <- after[!is.na(matches), ]
    message("Dropping ", nrow(before) - nrow(after), " rows")
  }
  message("Updated ", sum(!after[[col]] %in% before[[col]]), " MO codes")
  after
}
save_df <- function(...) {
  usethis::use_data(
    ...,
    overwrite = TRUE,
    version = 2,
    compress = "xz"
  )
}
microorganisms.codes <- fix_old_mos(AMR::microorganisms.codes, taxonomy, drop = TRUE)
if (!identical(microorganisms.codes, AMR::microorganisms.codes)) save_df(microorganisms.codes); rm(microorganisms.codes)

clinical_breakpoints <- fix_old_mos(AMR::clinical_breakpoints, taxonomy)
if (!identical(clinical_breakpoints, AMR::clinical_breakpoints)) save_df(clinical_breakpoints); rm(clinical_breakpoints)

example_isolates <- fix_old_mos(AMR::example_isolates, taxonomy)
if (!identical(example_isolates, AMR::example_isolates)) save_df(example_isolates); rm(example_isolates)

intrinsic_resistant <- fix_old_mos(AMR::intrinsic_resistant, taxonomy)
if (!identical(intrinsic_resistant, AMR::intrinsic_resistant)) save_df(intrinsic_resistant); rm(intrinsic_resistant)

microorganisms.groups <- fix_old_mos(AMR::microorganisms.groups, taxonomy, col = "mo")
microorganisms.groups <- fix_old_mos(microorganisms.groups, taxonomy, col = "mo_group")
if (!identical(microorganisms.groups, AMR::microorganisms.groups)) save_df(microorganisms.groups); rm(microorganisms.groups)


# *** Save to package *** -------------------------------------------------------------------------

# format to tibble and remove non-ASCII characters

saveRDS(taxonomy, "data-raw/taxonomy3b.rds")
# taxonomy <- readRDS("data-raw/taxonomy3b.rds")
taxonomy <- taxonomy %>%
  arrange(fullname) %>%
  select(
    mo,
    fullname,
    status,
    domain:subspecies,
    rank,
    ref,
    oxygen_tolerance,
    morphology,
    source,
    starts_with("lpsn"),
    starts_with("mycobank"),
    starts_with("gbif"),
    prevalence,
    snomed
  ) %>%
  AMR:::dataset_UTF8_to_ASCII()

microorganisms <- taxonomy

# set class <mo>
class(microorganisms$mo) <- c("mo", "character")
usethis::use_data(
  microorganisms,
  overwrite = TRUE,
  version = 2,
  compress = "xz"
)
rm(microorganisms)

# DON'T FORGET TO UPDATE R/_globals.R!

# load new data set
devtools::load_all(".")

anyNA(microorganisms$mo)
anyNA(microorganisms.codes$mo)
anyNA(intrinsic_resistant$mo)
anyNA(clinical_breakpoints$mo)

# load new data sets again
devtools::load_all(".")
source("data-raw/_pre_commit_checks.R")
devtools::load_all(".")


# run the unit tests
Sys.setenv(NOT_CRAN = "true")
testthat::test_file("tests/testthat/test-data.R")
testthat::test_file("tests/testthat/test-mo.R")
testthat::test_file("tests/testthat/test-mo_property.R")
