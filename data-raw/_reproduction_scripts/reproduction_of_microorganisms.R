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
# Center Groningen in the Netherlands, in collaboration with many      #
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


# ! THIS SCRIPT REQUIRES AT LEAST 16 GB RAM !
# (at least 12 GB will be used by the R session for the size of the files)

# GBIF:
# 1. Go to https://doi.org/10.48580/dgxjw and click Download in thd menu.
#    Check Classification and Extended, and download the file (~1.2 GB)
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
#    'Cell shape' on "*" and click Submit and click on the 'Download tabel as CSV' button.
#    Last time these worked: 
#    - wget -O bacdive_oxygen_tolerance.csv 'https://bacdive.dsmz.de/advsearch/csv?fg%5B1%5D%5Bgc%5D=OR&fg%5B1%5D%5Bfl%5D%5B1%5D%5Bfd%5D=Oxygen+tolerance&fg%5B1%5D%5Bfl%5D%5B1%5D%5Bfv%5D=%2A&fg%5B1%5D%5Bfl%5D%5B1%5D%5Bfvd%5D=oxygen_tolerance-oxygen_tol-4'
#    - wget -O bacdive_cell_shape.csv 'https://bacdive.dsmz.de/advsearch/csv?fg%5B1%5D%5Bgc%5D=OR&fg%5B1%5D%5Bfl%5D%5B1%5D%5Bfd%5D=Cell+shape&fg%5B1%5D%5Bfl%5D%5B1%5D%5Bfv%5D=%2A&fg%5B1%5D%5Bfl%5D%5B1%5D%5Bfvd%5D=cell_morphology-cell_shape-2'
# 6. Set these locations to the paths where the files are (or set the environment variable
#    AMR_TAXONOMY_SOURCE_DIR to the folder, e.g. for an unattended run):
folder_location <- Sys.getenv("AMR_TAXONOMY_SOURCE_DIR", unset = "~/Downloads/")
folder_location <- paste0(sub("/+$", "", path.expand(folder_location)), "/")
# (COL is downloaded as a COL Data Package (ColDP), unpack "NameUsage.tsv" and "metadata.yaml" to a folder "COL")
file_gbif <- paste0(folder_location, "COL/NameUsage.tsv")
file_lpsn <- paste0(folder_location, "taxonomy.csv")
file_mycobank <- paste0(folder_location, "MBList.xlsx")

file_bartlett <- "data-raw/bartlett_et_al_2022_human_pathogens.xlsx"

file_oxygen_tolerance <- paste0(folder_location, "bacdive_oxygen_tolerance.csv")
file_cell_shape <- paste0(folder_location, "bacdive_cell_shape.csv")

# 7. Run the rest of this script line by line and check everything :)
#    All review moments use review(): in an interactive session they open View(), otherwise they
#    print the first rows (50, or set options(AMR_build_print_rows = ...)). Set
#    options(AMR_build_view = FALSE) to only print, e.g. for a run with source(..., echo = TRUE)
#    and reading the log afterwards. For an unattended run by Claude, see
#    data-raw/_reproduction_scripts/reproduction_of_microorganisms_CLAUDE.md.

# The source files can only be downloaded by a human (registration, manual export buttons), so this
# script stops immediately if any of them is missing, instead of continuing with old data.
raw_files <- c(
  "GBIF/COL NameUsage.tsv (step 1)" = file_gbif,
  "LPSN taxonomy.csv (step 2)" = file_lpsn,
  "MycoBank MBList.xlsx (step 3)" = file_mycobank,
  "Bartlett et al. Excel file (step 4)" = file_bartlett,
  "BacDive bacdive_oxygen_tolerance.csv (step 5)" = file_oxygen_tolerance,
  "BacDive bacdive_cell_shape.csv (step 5)" = file_cell_shape
)
raw_missing <- raw_files[!file.exists(raw_files)]
if (length(raw_missing) > 0) {
  stop(
    "The following source files are missing and must be downloaded by a human first ",
    "(see steps 1-6 at the top of this script):\n",
    paste0("  - ", names(raw_missing), ": ", raw_missing, collapse = "\n"),
    call. = FALSE
  )
}
# downloaded files of more than 120 days old are probably from the previous update
raw_downloads <- raw_files[raw_files != file_bartlett]
raw_age <- as.numeric(Sys.Date() - as.Date(file.mtime(raw_downloads)), units = "days")
if (any(raw_age > 120)) {
  warning(
    "These source files are older than 120 days, check that they were downloaded for this update: ",
    toString(names(raw_downloads)[raw_age > 120]),
    call. = FALSE
  )
}
rm(raw_missing, raw_downloads, raw_age)

library(dplyr)
library(tidyr)
library(vroom) # to import files
library(rvest) # to scrape LPSN website
library(progress) # to show progress bars
library(readxl) # to read the MycoBank and Bartlett Excel files
devtools::load_all(".") # to load the AMR package

# keep the data set of the last release at hand (not that of the development version, which may contain
# flaws of an earlier build), it is the basis for relevance, earlier manually added taxa and all comparisons
previous_release <- system("git describe --tags --abbrev=0", intern = TRUE)
microorganisms_old <- local({
  tmp <- tempfile(fileext = ".rda")
  system2("git", c("show", paste0(previous_release, ":data/microorganisms.rda")), stdout = tmp)
  env <- new.env()
  load(tmp, envir = env)
  env$microorganisms
})
# releases until v3.0.1 had no `domain`, their kingdom was the domain (Bacteria, Archaea, Fungi, etc.)
if (!"domain" %in% colnames(microorganisms_old)) {
  microorganisms_old <- microorganisms_old %>%
    mutate(domain = kingdom, .before = kingdom)
}
message("Previous data set taken from release ", previous_release, " (", nrow(microorganisms_old), " rows)")
# the other data sets of this package contain codes of the development version, `AMR::microorganisms` is
# needed to translate these, and will be replaced at the end of this script
microorganisms_dev <- AMR::microorganisms

# the current names of `x` according to the development version, but only those in the same domain as `x` (`domain`, or
# else as in the development version): names that are not in that version are matched to the closest name by
# mo_current(), which can be another organism (e.g. in 2026, the puffball genus Lycoperdon for the filarial nematode
# genus Brugia, which was missing in the development version)
current_names_same_domain <- function(x, domain = NULL) {
  current <- as.character(suppressWarnings(suppressMessages(mo_current(x))))
  if (is.null(domain)) {
    domain <- microorganisms_dev$domain[match(x, microorganisms_dev$fullname)]
  }
  domain_current <- microorganisms_dev$domain[match(current, microorganisms_dev$fullname)]
  unique(current[!is.na(domain) & !is.na(domain_current) & domain == domain_current])
}

# the integrity rules of the data set and the MO code registry, shared with the unit tests
source("tests/testthat/helper-microorganisms.R")


# Helper functions --------------------------------------------------------------------------------

# show data that must be reviewed, without blocking a non-interactive run
# with options(AMR_build_review_dir = "some/folder"), every review is also saved as a CSV file in that
# folder, so that a human can review them afterwards (e.g. in a pull request). The file name is based on
# the title only, so running a section again overwrites its own file.
review <- function(x, title = deparse(substitute(x))) {
  message("\n>>> REVIEW: ", title, " (", nrow(x), " rows)")
  review_dir <- getOption("AMR_build_review_dir")
  if (!is.null(review_dir)) {
    dir.create(review_dir, recursive = TRUE, showWarnings = FALSE)
    file_name <- paste0(gsub("[^a-z0-9]+", "_", gsub("^[^a-z0-9]+|[^a-z0-9]+$", "", tolower(title))), ".csv")
    out <- as.data.frame(x)
    # list columns (such as `snomed`) cannot be written to CSV
    out[] <- lapply(out, function(col) if (is.list(col)) vapply(col, toString, character(1)) else col)
    utils::write.csv(out, file.path(review_dir, file_name), row.names = FALSE, na = "")
  }
  if (nrow(x) > 0) {
    if (interactive() && isTRUE(getOption("AMR_build_view", TRUE))) {
      utils::View(x, title)
    } else {
      print(utils::head(as.data.frame(x), getOption("AMR_build_print_rows", 50)))
    }
  }
  invisible(x)
}

# taxonomic names in ASCII and without quotes
ascii_names <- function(x) {
  gsub("[\"'`]", "", stringi::stri_trans_general(x, "Latin-ASCII"))
}

# Microsporidia are fungi (decision by Matthijs S. Berends, 5 October 2026), while COL places them in the Protozoa and
# releases until v3.0.1 had most of them in the Protozoa as well; their codes of earlier releases are translated by name
is_microsporidian <- function(phylum, class) {
  phylum %in% c("Microsporidia", "Microspora", "Chytridiopsidomycota") | class %in% "Microsporea"
}
microsporidia_to_fungi <- function(df) {
  fix <- is_microsporidian(df$phylum, df$class)
  df$domain[fix] <- "Fungi"
  df$kingdom[fix] <- "Fungi"
  df
}

# all names in the sources with their domain, as "domain name" (for GBIF all accepted names, not only the ones
# selected by this script), to check whether a name still exists in a domain
names_in_sources <- function() {
  gbif_accepted <- is.na(taxonomy_gbif.bak$acceptedNameUsageID)
  gbif_domain <- if_else(
    is_microsporidian(taxonomy_gbif.bak$phylum, taxonomy_gbif.bak$class), "Fungi", taxonomy_gbif.bak$kingdom
  )
  unique(c(
    # (the kingdom is the domain for all but the prokaryotes, which are covered by LPSN)
    paste(gbif_domain[gbif_accepted], taxonomy_gbif.bak$scientificName[gbif_accepted]),
    paste(taxonomy_gbif$domain, taxonomy_gbif$fullname),
    paste(taxonomy_mycobank$domain, taxonomy_mycobank$fullname),
    paste(taxonomy_lpsn$domain, trimws(gsub(" +", " ", paste(
      taxonomy_lpsn$genus, coalesce(taxonomy_lpsn$species, ""), coalesce(taxonomy_lpsn$subspecies, "")
    ))))
  ))
}

# the priority of sources, used everywhere in this script
source_priority <- c("LPSN" = 1L, "MycoBank" = 2L, "GBIF" = 3L, "manually added" = 4L, "inferred" = 5L)
source_prio <- function(source) {
  coalesce(unname(source_priority[source]), 6L)
}

# 'inferred' is used during this script for records that are not in any source, but that we must
# create ourselves (e.g. a genus of which only species are available). These are NOT manually curated
# and must therefore follow all the relevance filters. They are renamed to "manually added" at the end.

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
    authors2 %like_case% ".*[)] [a-zA-Z]+.*",
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
  ref[ref %unlike_case% "^d[A-Z]"] <- gsub(
    "^([a-z])",
    "\\U\\1",
    ref[ref %unlike_case% "^d[A-Z]"],
    perl = TRUE
  )
  # specific one for the French that are named dOrbigny
  ref[ref %like_case% "^d[A-Z]"] <- gsub("^d", "d'", ref[ref %like_case% "^d[A-Z]"])
  ref <- gsub(" +", " ", ref)
  ref <- trimws(ref)
  ref <- gsub("^NA, ?", "", ref)
  ref[ref %in% c("", "NA")] <- NA_character_
  ref
}

# to retrieve LPSN and authors from LPSN website
# e.g., get_lpsn_and_author("genus", "Klebsiella")
# results are cached in `lpsn_cache` (and in data-raw/lpsn_scrape_cache.rds using save_lpsn_cache()),
# so an interrupted run does not have to download thousands of pages again
lpsn_cache_file <- "data-raw/lpsn_scrape_cache.rds"
lpsn_cache <- new.env()
if (file.exists(lpsn_cache_file)) {
  # only use a cache of this month, the LPSN is updated continuously
  if (format(file.mtime(lpsn_cache_file), "%Y-%m") == format(Sys.Date(), "%Y-%m")) {
    list2env(readRDS(lpsn_cache_file), envir = lpsn_cache)
  }
}
save_lpsn_cache <- function() {
  saveRDS(as.list(lpsn_cache), lpsn_cache_file, version = 2)
}
get_lpsn_and_author <- function(rank, name, tries = 3) {
  name <- gsub("^Candidatus ", "", name)
  # (subspecies without "subsp.", e.g. /subspecies/clavibacter-michiganensis-nebraskensis)
  name <- gsub(" subsp[.] ", " ", name)
  url <- paste0(
    "https://lpsn.dsmz.de/",
    tolower(rank), "/",
    gsub(" ", "-", tolower(name))
  )
  # (cached results of before October 2026 lack the taxonomic status, so these are retrieved again)
  if (!is.null(lpsn_cache[[url]]) && "taxonomic_status" %in% names(lpsn_cache[[url]])) {
    return(lpsn_cache[[url]])
  }
  page_txt <- NULL
  not_found <- FALSE
  for (try in seq_len(tries)) {
    # be polite to the LPSN server
    Sys.sleep(0.1)
    page_txt <- tryCatch({
      resp <- curl::curl_fetch_memory(url)
      if (resp$status_code == 200) {
        read_html(rawToChar(resp$content))
      } else {
        # 404 is a valid answer, other codes (e.g. 429, 503) are worth a retry
        not_found <- resp$status_code == 404
        NULL
      }
    }, error = function(e) NULL)
    if (!is.null(page_txt) || not_found) {
      break
    }
    Sys.sleep(5 * try)
  }
  if (!is.null(page_txt)) {
    page_txt <- page_txt %>%
      html_element("#detail-page") %>%
      html_text()
  }
  # a page without a record number is no valid LPSN record (otherwise the whole page text would be
  # returned by gsub() below and stored as the LPSN ID)
  if (is.null(page_txt) || is.na(page_txt) || page_txt %unlike_case% "Record number:") {
    if (!not_found) {
      warning("No LPSN found for ", tolower(rank), " '", name, "'", call. = FALSE)
    }
    out <- c(
      "lpsn" = NA_character_, "ref" = NA_character_, "status" = "unknown",
      "taxonomic_status" = NA_character_, "correct_name" = NA_character_
    )
  } else {
    lpsn <- gsub(
      ".*Record number:[\r\n\t ]*([0-9]+).*",
      "\\1",
      page_txt,
      perl = FALSE
    )
    if (lpsn %unlike_case% "^[0-9]+$") {
      lpsn <- NA_character_
    }
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
    # the nomenclatural status above does not tell whether a name is current: a validly published name can be a
    # synonym (e.g. Eubacterium lentum, with the correct name Eggerthella lenta)
    page_flat <- gsub("[\r\n\t ]+", " ", page_txt)
    taxonomic_status <- if (page_flat %like_case% "Taxonomic status:") {
      trimws(sub(".*?Taxonomic status: *([a-z]+( [a-z]+)?).*", "\\1", page_flat))
    } else {
      NA_character_
    }
    correct_name <- NA_character_
    if (taxonomic_status %in% "synonym" && page_flat %like_case% "Correct name: ") {
      # e.g. 'Correct name: Eggerthella lenta (Eggerth 1935) Wade et al. 1999', without authors and "subsp."
      correct_name <- sub(".*?Correct name: \"?([A-Z][a-z]+( (subsp[.] )?[a-z][a-z-]+){0,2}).*", "\\1", page_flat)
      correct_name <- gsub(" subsp[.] ", " ", correct_name)
      status <- "synonym"
    }
    out <- c(
      "lpsn" = lpsn, "ref" = ref, "status" = status,
      "taxonomic_status" = taxonomic_status, "correct_name" = correct_name
    )
  }
  # do not cache network failures, so they will be retried in a next run
  if (!is.null(page_txt) || not_found) {
    lpsn_cache[[url]] <- out
  }
  out
}

# Applies a result of get_lpsn_and_author() to record `i` of `df`. A synonym in LPSN is linked to its correct name if
# that is in `df` (by the first identifier the correct name has, LPSN > MycoBank > GBIF); without its correct name in
# `df` it stays a synonym without a current name, which the build resolves in 'Synonyms without a current name'.
apply_lpsn_result <- function(df, i, lpsn) {
  df$source[i] <- "LPSN"
  df$lpsn[i] <- unname(lpsn["lpsn"])
  df$ref[i] <- unname(lpsn["ref"])
  df$status[i] <- unname(lpsn["status"])
  df[i, c("lpsn_renamed_to", "mycobank_renamed_to", "gbif_renamed_to")] <- NA_character_
  if (lpsn["status"] %in% "synonym") {
    target <- which(df$fullname == lpsn["correct_name"] & df$domain == df$domain[i] & df$status != "synonym")
    if (length(target) == 1) {
      for (id in c("lpsn", "mycobank", "gbif")) {
        if (!is.na(df[[id]][target])) {
          df[[paste0(id, "_renamed_to")]][i] <- df[[id]][target]
          break
        }
      }
    }
  }
  df
}

# get the LPSN record of all unique taxa of a certain rank (with their higher taxonomy),
# used for families until domains below. Binding once at the end is much faster than binding in the loop.
get_lpsn_rank <- function(taxonomy_lpsn, rank) {
  higher_ranks <- c("domain", "kingdom", "phylum", "class", "order", "family")
  higher_ranks <- higher_ranks[seq_len(which(higher_ranks == rank) - 1)]
  taxa <- taxonomy_lpsn %>%
    filter(!is.na(.data[[rank]])) %>%
    distinct(.data[[rank]], .keep_all = TRUE) %>%
    select(all_of(c(higher_ranks, rank)))
  pb <- progress_bar$new(
    total = nrow(taxa),
    format = paste0(rank, " [:bar] :current/:total :eta")
  )
  info <- lapply(taxa[[rank]], function(nm) {
    pb$tick()
    get_lpsn_and_author(tools::toTitleCase(rank), nm)
  })
  save_lpsn_cache()
  taxa %>%
    mutate(
      rank = rank,
      status = vapply(info, function(x) unname(x["status"]), character(1)),
      source = "LPSN",
      lpsn = vapply(info, function(x) unname(x["lpsn"]), character(1)),
      ref = vapply(info, function(x) unname(x["ref"]), character(1))
    )
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
    pick <- function(current, rank) {
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
    out <- ""
    if (n_distinct(source) > 1 && "GBIF" %in% source) {
      # prefer LPSN and MycoBank over GBIF, which is often not up-to-date at all
      out <- pick(current[source != "GBIF"], rank[source != "GBIF"])
    }
    if (out == "") {
      # but never return an empty value if GBIF does have one (e.g. Mucor has no family in MycoBank)
      out <- pick(current, rank)
    }
    out
  }
}

# for a set of records, take the most reliable non-empty value: of the record with the own rank
# (e.g. the family record itself when determining the order of a family), else of the most
# authoritative source, else the most common value
consensus_value <- function(x, rank, source, own_rank) {
  ok <- !x %in% c("", NA)
  if (!any(ok)) {
    return("")
  }
  x <- x[ok]
  rank <- rank[ok]
  source <- source[ok]
  own <- x[rank == own_rank]
  if (length(own) > 0) {
    return(own[order(source_prio(source[rank == own_rank]))][1])
  }
  best_source <- min(source_prio(source))
  x <- x[source_prio(source) == best_source]
  names(sort(table(x), decreasing = TRUE))[1]
}

# all records whose fullname is not in `taxonomy` yet, but which are referred to by other records
# (e.g. a family of which only genera are available) - these will be added as 'inferred' records,
# with identifiers from GBIF where available
add_missing_parents <- function(taxonomy, current_gbif) {
  rank_levels <- c("domain", "kingdom", "phylum", "class", "order", "family", "genus")
  missing_higher <- bind_rows(lapply(rank_levels, function(rank_name) {
    candidates <- taxonomy %>%
      filter(.data[[rank_name]] != "") %>%
      distinct(across(all_of(unique(c("domain", rank_levels[seq_len(which(rank_levels == rank_name))]))))) %>%
      mutate(fullname = .data[[rank_name]], rank = rank_name) %>%
      # placeholders such as "(unknown class)" are no taxa
      filter(fullname %unlike% "^[(]unknown") %>%
      anti_join(
        taxonomy %>% filter(rank == rank_name),
        by = unique(c("domain", setNames(rank_name, rank_name), "rank"))
      ) %>%
      # one record per taxon, inconsistent higher taxonomy will be fixed later
      distinct(domain, fullname, .keep_all = TRUE)
    if (nrow(candidates) == 0) {
      return(tibble())
    }
    message("Adding missing: ", rank_name, "... n = ", nrow(candidates))
    # enrich from GBIF
    gbif_lookup <- current_gbif %>%
      filter(taxonRank == rank_name) %>%
      select(
        all_of(unique(c("domain", rank_name))),
        ref = scientificNameAuthorship,
        gbif = taxonID,
        gbif_parent = parentNameUsageID
      ) %>%
      distinct(across(all_of(unique(c("domain", rank_name)))), .keep_all = TRUE)
    candidates %>%
      left_join(gbif_lookup, by = unique(c("domain", rank_name))) %>%
      mutate(
        source = if_else(!is.na(gbif), "GBIF", "inferred"),
        status = if_else(!is.na(gbif), "accepted", "unknown")
      )
  }))
  # species implied by subspecies but missing as a species-rank row
  missing_species <- taxonomy %>%
    # placeholders such as "(unknown species)" are no taxa
    filter(species != "", genus %unlike% "^[(]unknown", species %unlike% "^[(]unknown") %>%
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
        select(
          domain, genus,
          species = specificEpithet,
          ref = scientificNameAuthorship,
          gbif = taxonID,
          gbif_parent = parentNameUsageID
        ) %>%
        distinct(domain, genus, species, .keep_all = TRUE),
      by = c("domain", "genus", "species")
    ) %>%
    mutate(
      source = if_else(!is.na(gbif), "GBIF", "inferred"),
      status = if_else(!is.na(gbif), "accepted", "unknown")
    )
  if (nrow(missing_species) > 0) {
    message("Adding missing: species... n = ", nrow(missing_species))
  }
  taxonomy %>%
    bind_rows(missing_higher, missing_species) %>%
    mutate(
      domain = if_else(is.na(domain) | domain == "", kingdom, domain),
      across(kingdom:subspecies, \(x) if_else(is.na(x), "", x))
    )
}

# Clinically relevant genera ----------------------------------------------------------------------

# One set of genera that must never get lost, used by all filters below. It contains:
# - our own curated list of non-bacterial genera (MO_RELEVANT_GENERA, see data-raw/_pre_commit_checks.R)
# - the WHO priority genera
# - all genera of the human pathogens in Bartlett et al. (2022)
# - all genera that are used in the other data sets of this package (breakpoints, species groups, codes, etc.), but
#   not of intrinsic_resistant: that is computed from this data set afterwards, and as it lists nearly all bacterial
#   genera, using it here made every build depend on the previous one (decision by Matthijs S. Berends, 7 October 2026)
# - all genera that the interpretive rules name (EUCAST expected phenotypes and expert rules), read from the rules
#   themselves, so that they stay in sync with them
# (until October 2026, also all non-bacterial genera that were relevant in the previous version of this data set, but
# that made every build depend on the previous one; their genera are in MO_RELEVANT_GENERA now, decision by Matthijs
# S. Berends, 7 October 2026)
# - and of all these: their current names and their synonyms (e.g. Candida -> Candidozyma, Nakaseomyces)
pathogens <- read_excel(file_bartlett, sheet = "Tab 6 Full List")

genera_of_mo <- function(x) {
  # these data sets contain codes of the development version, which may not be in the last release
  unique(microorganisms_dev$genus[match(as.character(x), as.character(microorganisms_dev$mo))])
}
genera_in_interpretive_rules <- function() {
  r <- AMR:::INTERPRETIVE_RULES_DF
  r <- r[r$like.is.one_of %in% c("is", "one_of") & r$if_mo_property %in% c("genus", "genus_species", "fullname"), , drop = FALSE]
  nms <- trimws(unlist(strsplit(r$this_value, ",", fixed = TRUE)))
  # group names, such as 'Enterobacter cloacae complex', are taken by the genera of their members
  in_groups <- AMR::microorganisms.groups$mo_name[AMR::microorganisms.groups$mo_group_name %in% nms]
  unique(sub(" .*", "", c(nms[!nms %in% AMR::microorganisms.groups$mo_group_name], in_groups)))
}
relevant_genera <- c(
  AMR:::MO_RELEVANT_GENERA,
  genera_in_interpretive_rules(),
  AMR:::MO_WHO_PRIORITY_GENERA,
  pathogens$genus,
  genera_of_mo(AMR::clinical_breakpoints$mo),
  genera_of_mo(AMR::microorganisms.groups$mo),
  genera_of_mo(AMR::microorganisms.groups$mo_group),
  genera_of_mo(AMR::microorganisms.codes$mo),
  genera_of_mo(AMR::example_isolates$mo)
)
relevant_genera <- relevant_genera[!relevant_genera %in% c("", NA) & relevant_genera %unlike% "unknown"]
relevant_genera <- c(
  relevant_genera,
  # synonyms of the genera themselves (only of names in the development version, as mo_synonyms() would otherwise
  # look up the closest match, e.g. a fungus for a parasite that was missing in that version)
  unlist(mo_synonyms(intersect(relevant_genera, microorganisms_dev$fullname), keep_synonyms = FALSE))
) %>%
  strsplit(" ", fixed = TRUE) %>%
  vapply(function(x) x[1], character(1)) %>%
  unique() %>%
  sort()
relevant_genera <- relevant_genera[!relevant_genera %in% c("", NA) & relevant_genera %unlike% "unknown"]
message(length(relevant_genera), " clinically relevant genera will be protected")
# The current names of all species of these genera, e.g. Candida auris -> Candidozyma auris, are protected as well,
# with the record of their genus. Not their whole genus: in the previous data set, the prevalence of species was set
# per genus, so e.g. every historical Fusarium name was 'relevant', including plant pathogens that are now in
# Cercospora or Colletotrichum, which would make these large plant pathogenic genera 'clinically relevant' entirely
# (in 2026, this tripled the size of the data set to 222,000 records).
relevant_current_species <- microorganisms_old %>%
  filter(genus %in% relevant_genera, rank %in% c("species", "subspecies"))
relevant_current_species <- current_names_same_domain(relevant_current_species$fullname, relevant_current_species$domain)
relevant_current_species <- relevant_current_species[!relevant_current_species %in% c("", NA) &
  relevant_current_species %unlike% "unknown"]
relevant_current_genera <- setdiff(unique(gsub(" .*", "", relevant_current_species)), relevant_genera)
message(
  length(relevant_current_species), " current names of their species will be protected, with ",
  length(relevant_current_genera), " additional genus records"
)

# Genus names that exist in multiple domains (homonyms) cannot always be resolved automatically.
# E.g., Capillaria (nematode vs fungus) and Necator (hookworm vs fungus) would otherwise end up as
# fungi, since Fungi have priority over Animalia. These overrides have the highest priority.
# A report of all unresolved homonyms of relevant genera is shown below in 'Combine the datasets',
# so extend this list there when needed.
genus_domain_override <- c(
  # the labyrinthulids are stramenopiles, not fungi, as all other labyrinthulids in COL (decision by Matthijs S.
  # Berends, 7 October 2026; until v3.0.1, Aplanochytrium was in the Fungi)
  "Aplanochytrium" = "Chromista",
  "Capillaria" = "Animalia",
  "Graphium" = "Fungi",
  "Hymenolepis" = "Animalia",
  "Necator" = "Animalia",
  "Schizophyllum" = "Fungi",
  "Trichophyton" = "Fungi",
  "Wangiella" = "Fungi"
)

# Genera of which the released data sets contained other organisms with the same genus name, that cannot be told apart
# from the right ones (e.g. until v3.0.1, the Graphium butterflies were in the Fungi with a fungal lineage). Of these
# genera, released species and subspecies are only kept if a current source still has their name in this domain, the
# codes of the others are retired (and must be reviewed).
released_homonym_genera <- c(
  "Graphium" = "Fungi"
)

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

# All taxa in the genera of the WHO priority pathogen lists (AMR:::MO_WHO_PRIORITY_GENERA) are always protected, also
# their names that are not validly published under the ICNP, such as Mycobacterium canettii, Mycobacterium orygis
# and Klebsiella quasivariicola, which LPSN has as preferred names (decision by Matthijs S. Berends, 7 October 2026)
is_who_priority_genus <- function(genus) genus %in% AMR:::MO_WHO_PRIORITY_GENERA


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

# duplicate genera (homonyms) would duplicate rows in the join below:
review(
  taxonomy_lpsn_missing %>% filter(genus %in% taxonomy_lpsn_missing$genus[duplicated(taxonomy_lpsn_missing$genus)]),
  "LPSN genera with multiple lineages"
)
# look them up on LPSN, then find out which to keep (the ones validly published under ICNP)
# to remove:
taxonomy_lpsn_missing <- taxonomy_lpsn_missing %>%
  filter(
    !(genus == "Halalkalibacterium" & family == "Balneolaceae"),
    !(genus == "Pusillimonas" & family == "Oscillospiraceae"),
    !(genus == "Rhodococcus" & family == "Chroococcaceae")
  )
# any remaining duplicates: keep the most complete lineage, so that the join below cannot duplicate records
if (anyDuplicated(taxonomy_lpsn_missing$genus) > 0) {
  warning(
    "Genera with multiple LPSN lineages left, keeping the most complete one: ",
    toString(unique(taxonomy_lpsn_missing$genus[duplicated(taxonomy_lpsn_missing$genus)])),
    call. = FALSE
  )
  taxonomy_lpsn_missing <- taxonomy_lpsn_missing %>%
    arrange(genus, rowSums(is.na(across(c(domain, kingdom, phylum, class, order, family))))) %>%
    distinct(genus, .keep_all = TRUE)
}

taxonomy_lpsn <- taxonomy_lpsn %>%
  left_join(taxonomy_lpsn_missing, by = "genus", relationship = "many-to-one") %>%
  select(domain, kingdom, phylum, class, order, family, everything()) %>%
  # remove entries like "[Bacteria, no family]" and "[Bacteria, no class]"
  mutate(across(
    c(domain, kingdom, phylum, class, order, family),
    function(x) if_else(x %like_case% " no ", NA_character_, x)
  ))

taxonomy_lpsn.bak2 <- taxonomy_lpsn
# download higher taxonomic data (e.g. authors of Enterobacteriaceae) directly from LPSN website using
# scraping, by using get_lpsn_and_author(). Try it first:
# get_lpsn_and_author("genus", "Escherichia")
# get_lpsn_and_author("family", "Enterobacteriaceae")
# get_lpsn_and_author("order", "Enterobacterales")
# get_lpsn_and_author("class", "Gammaproteobacteria")
# get_lpsn_and_author("phylum", "Pseudomonadota")
# get_lpsn_and_author("kingdom", "Pseudomonadati")
# get_lpsn_and_author("domain", "Bacteria")
taxonomy_lpsn <- taxonomy_lpsn %>%
  bind_rows(
    get_lpsn_rank(taxonomy_lpsn, "family"),
    get_lpsn_rank(taxonomy_lpsn, "order"),
    get_lpsn_rank(taxonomy_lpsn, "class"),
    get_lpsn_rank(taxonomy_lpsn, "phylum"),
    get_lpsn_rank(taxonomy_lpsn, "kingdom"),
    get_lpsn_rank(taxonomy_lpsn, "domain")
  )
# these could not be downloaded at all (network errors), run the lines above again to retry them:
review(taxonomy_lpsn %>% filter(status == "unknown"), "LPSN higher taxa that could not be downloaded")

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
    # without a year, 'Smith, NA' would otherwise become 'Smith et al.' in get_author_year()
    ref = if_else(
      is.na(`Year of effective publication`),
      as.character(Authors),
      paste0(Authors, ", ", `Year of effective publication`)
    ),
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
      ),
    # names with other characters than Latin letters, spaces and hyphens are no valid names, such as
    # 'Cordyceps *jezoensoides' or 'Pyrenopeziza poae' with Cyrillic letters that look like Latin ones
    stringi::stri_detect_regex(fullname, "^[\\p{Script=Latin} -]+$")
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
  # and make fullname distinct again, preferring accepted names (NA would otherwise be sorted last)
  arrange(rank, !is.na(mycobank_renamed_to), mycobank_renamed_to) %>%
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
# MycoBank pasted their taxonomy together into 1 field (split into tax_a, tax_b, ... above), and the
# position of each rank differs per record, since e.g. subphyla and subclasses are not always present.
# Previously we only looked in fixed columns per rank, which left e.g. Mucor without phylum, class,
# order and family. So now for every rank, we take the first tax_* column with a known name of that rank.
tax_cols <- grep("^tax_", colnames(taxonomy_mycobank), value = TRUE)
first_name_of_rank <- function(df, rank_name) {
  known <- unique(df$fullname[df$rank == rank_name])
  vals <- lapply(df[tax_cols], function(x) if_else(x %in% known, x, NA_character_))
  coalesce(do.call(coalesce, unname(vals)), "")
}
taxonomy_mycobank <- taxonomy_mycobank %>%
  mutate(
    kingdom = "Fungi", # we already filtered everything else, and MycoBank 90157 has '-' as current name...
    phylum = if_else(rank == "phylum", fullname, first_name_of_rank(taxonomy_mycobank, "phylum")),
    class = if_else(rank == "class", fullname, first_name_of_rank(taxonomy_mycobank, "class")),
    order = if_else(rank == "order", fullname, first_name_of_rank(taxonomy_mycobank, "order")),
    family = if_else(rank == "family", fullname, first_name_of_rank(taxonomy_mycobank, "family")),
    genus = case_when(
      rank == "genus" ~ fullname,
      # the genus of a (sub)species is always its first word, never rely on the classification for this
      rank %in% c("species", "subspecies") ~ gsub(" .*", "", fullname),
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
# these should be very few now:
taxonomy_mycobank %>%
  filter(rank %in% c("genus", "species"), family == "") %>%
  count(rank)

# (sub)species without a family in MycoBank get the most common family of their genus
taxonomy_mycobank <- taxonomy_mycobank %>%
  group_by(genus) %>%
  mutate(across(
    c(phylum, class, order, family),
    function(x) {
      if (genus[1] == "" || all(x != "") || all(x == "")) {
        x
      } else {
        if_else(x == "", names(sort(table(x[x != ""]), decreasing = TRUE))[1], x)
      }
    }
  )) %>%
  ungroup()

# keep only the relevant genera (see 'Clinically relevant genera' above), but not the fungal homonyms
# of relevant non-fungal genera, such as Necator (otherwise the hookworm would become a fungus)
genera_overview <- relevant_genera[
  !relevant_genera %in% names(genus_domain_override)[genus_domain_override != "Fungi"]
]
include <- taxonomy_mycobank %>%
  filter(
    genus %in% genera_overview | !rank %in% c("genus", "species", "subspecies") |
      fullname %in% relevant_current_species | (rank == "genus" & fullname %in% relevant_current_genera)
  ) %>%
  filter(!(genus == "" & rank %in% c("genus", "species", "subspecies")))
missing_renamed_to <- include$mycobank_renamed_to[which(
  !is.na(include$mycobank_renamed_to) &
    !include$mycobank_renamed_to %in% include$mycobank
)]
renamed_to_records <- taxonomy_mycobank %>%
  filter(mycobank %in% missing_renamed_to)
taxonomy_mycobank <- include %>%
  bind_rows(
    renamed_to_records,
    # and the genus records of these current names, e.g. a new genus of a renamed relevant species
    taxonomy_mycobank %>%
      filter(rank == "genus", fullname %in% renamed_to_records$genus, !mycobank %in% include$mycobank)
  ) %>%
  distinct(mycobank, .keep_all = TRUE) %>%
  arrange(fullname)
rm(renamed_to_records)

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

# @resume-block gbif_raw after=taxonomy_gbif
# (needed by `current_gbif` when resuming from a checkpoint with run_microorganisms_build.R)
# COL is downloaded as a COL Data Package (ColDP), of which NameUsage.tsv contains all names. Its columns are
# translated to the Darwin Core names that this script was built on. In ColDP, `parentID` is the parent of an
# accepted name, but the accepted name of a synonym (as `acceptedNameUsageID` in Darwin Core). Synonyms have
# no higher taxonomy in either format. The file is not quoted, so quotes must not be interpreted.
taxonomy_gbif.bak <- vroom(
  file_gbif,
  delim = "\t",
  quote = "",
  col_types = cols(.default = col_character()),
  col_select = c(
    "col:ID", "col:parentID", "col:status", "col:scientificName", "col:authorship", "col:rank",
    "col:specificEpithet", "col:infraspecificEpithet",
    "col:kingdom", "col:phylum", "col:class", "col:order", "col:family", "col:genus", "clb:merged"
  )
)
colnames(taxonomy_gbif.bak) <- gsub(".*:(.*)", "\\1", colnames(taxonomy_gbif.bak))
taxonomy_gbif.bak <- taxonomy_gbif.bak %>%
  transmute(
    taxonID = ID,
    parentNameUsageID = if_else(status %in% c("accepted", "provisionally accepted"), parentID, NA_character_),
    acceptedNameUsageID = if_else(status %in% c("accepted", "provisionally accepted"), NA_character_, parentID),
    taxonomicStatus = status,
    taxonRank = rank,
    scientificName,
    scientificNameAuthorship = authorship,
    specificEpithet,
    infraspecificEpithet,
    kingdom,
    phylum,
    class,
    order,
    family,
    genus,
    # TRUE for names that the COL eXtended Release (XR) added programmatically to the curated base release
    xr_merged = merged %in% "true"
  )
# @end-resume-block

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
taxonomy_gbif0 <- taxonomy_gbif %>%
  # names that the COL eXtended Release added programmatically are only used if they are clinically relevant or were
  # part of the last release (decision by Matthijs S. Berends, 5 October 2026): the XR adds many names of which the
  # quality has not been checked by the COL editors (synonyms have no genus, so the genus is taken from the name)
  filter(
    !xr_merged |
      coalesce(genus, sub(" .*", "", fullname)) %in% relevant_genera |
      fullname %in% relevant_current_species |
      fullname %in% microorganisms_old$fullname
  )

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
      # and all relevant genera (see 'Clinically relevant genera' above)
      # (they also contain bacteria and protozoa, but these will get higher prevalence scores later on)
      genus %in% relevant_genera
  ) %>%
  mutate(
    # set the right domain for the prokaryotes based on LPSN
    domain = case_when(
      kingdom %in% na.omit(taxonomy_lpsn$kingdom[taxonomy_lpsn$domain == "Bacteria"]) ~ "Bacteria",
      kingdom %in% na.omit(taxonomy_lpsn$kingdom[taxonomy_lpsn$domain == "Archaea"]) ~ "Archaea",
      TRUE ~ kingdom
    )
  ) %>%
  microsporidia_to_fungi()
# add all synonyms of the included records (most synonyms in GBIF lack their higher taxonomy, so the
# filter above does not catch them)
taxonomy_gbif <- taxonomy_gbif %>%
  bind_rows(
    taxonomy_gbif0 %>%
      filter(
        acceptedNameUsageID %in% taxonomy_gbif$taxonID,
        !taxonID %in% taxonomy_gbif$taxonID,
        !taxonomicStatus %in% c("doubtful", "misapplied", "ambiguous synonym")
      )
  )

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
    acc_family = family,
    acc_genus = genus
  ) %>%
  distinct(gbif, .keep_all = TRUE)

# the higher taxonomy of accepted genera, to be used for synonyms of species in another genus than
# their accepted name (e.g. 'Ascaris cati' is a synonym of 'Toxocara cati', but its family must be the one
# of Ascaris, not of Toxocara)
accepted_genera <- taxonomy_gbif %>%
  filter(rank == "genus", status %like% "accepted", !is.na(genus), !is.na(domain)) %>%
  distinct(domain, genus, .keep_all = TRUE) %>%
  select(
    domain, genus,
    gen_phylum = phylum,
    gen_class = class,
    gen_order = order,
    gen_family = family
  )

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
    
    # Step 2: fill domain and kingdom from the accepted name, only for ranks ABOVE the synonym's rank
    domain  = if_else(is.na(domain) & !is.na(acc_domain), acc_domain, domain),
    kingdom = if_else(is.na(kingdom) & !is.na(acc_kingdom) & rank_num > 1L, acc_kingdom, kingdom),
    # Step 3: phylum to family may only come from the accepted name if that is in the same genus
    # (or if the synonym is a genus or higher), otherwise from the synonym's own genus (see Step 4)
    same_lineage = rank_num <= 6L | (!is.na(acc_genus) & genus == acc_genus),
    phylum  = if_else(is.na(phylum) & !is.na(acc_phylum) & rank_num > 2L & same_lineage, acc_phylum, phylum),
    class   = if_else(is.na(class) & !is.na(acc_class) & rank_num > 3L & same_lineage, acc_class, class),
    order   = if_else(is.na(order) & !is.na(acc_order) & rank_num > 4L & same_lineage, acc_order, order),
    family  = if_else(is.na(family) & !is.na(acc_family) & rank_num > 5L & same_lineage, acc_family, family)
  ) %>%
  # Step 4: species in another genus than their accepted name get the taxonomy of their own genus
  left_join(accepted_genera, by = c("domain", "genus")) %>%
  mutate(
    phylum = if_else(is.na(phylum) & rank_num > 6L, gen_phylum, phylum),
    class = if_else(is.na(class) & rank_num > 6L, gen_class, class),
    order = if_else(is.na(order) & rank_num > 6L, gen_order, order),
    family = if_else(is.na(family) & rank_num > 6L, gen_family, family)
  ) %>%
  select(-starts_with("acc_"), -starts_with("gen_"), -name_parts, -rank_num, -same_lineage)
rm(accepted, accepted_genera)


taxonomy_gbif <- taxonomy_gbif %>%
  mutate(ref = get_author_year(ref)) %>%
  # (records without authors must be kept: `NA != "AmSOD"` is NA, which filter() would drop, as happened in 2026 to
  # e.g. Balamuthia mandrillaris and thousands of other names without authors in COL)
  filter(is.na(ref) | ref != "AmSOD")

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

# names must be ASCII before they are made unique, otherwise e.g. 'Fusarium aloës' and 'Fusarium aloes'
# both end up as 'Fusarium aloes' in the final data set (which happened to 12 names in 2026), see ascii_names()
taxonomy_gbif <- taxonomy_gbif %>%
  mutate(across(domain:subspecies, ascii_names))
taxonomy_lpsn <- taxonomy_lpsn %>%
  mutate(across(domain:subspecies, ascii_names))
taxonomy_mycobank <- taxonomy_mycobank %>%
  mutate(across(domain:subspecies, ascii_names))

# Keeps one record per name (kingdom, rank and full name) of a source. A source can contain the same name more than
# once, e.g. LPSN has the correct name Eggerthella lenta Wade et al. 1999 and the illegitimate homotypic synonym
# Eggerthella lenta Kageyama et al. 1999, and COL has Anisakis simplex Dujardin, 1845 (accepted) and Anisakis simplex
# Rudolphi, 1809 (synonym). The accepted record is kept, and every link to a dropped record is moved to the kept record.
# (Until October 2026, the record with the lowest identifier was kept, which made e.g. Eggerthella lenta, Gordonia
# amarae and Anisakis simplex synonyms without a current name.)
one_record_per_name <- function(df, id) {
  renamed_to <- paste0(id, "_renamed_to")
  df <- df %>%
    arrange(fullname, status != "accepted", .data[[id]]) %>%
    group_by(kingdom, rank, fullname) %>%
    mutate(.kept_id = first(.data[[id]])) %>%
    ungroup()
  dropped <- df %>%
    filter(!is.na(.data[[id]]), .data[[id]] != .kept_id)
  relink <- match(df[[renamed_to]], dropped[[id]])
  df[[renamed_to]][!is.na(relink)] <- dropped$.kept_id[relink[!is.na(relink)]]
  # a record never points to itself
  df[[renamed_to]][which(df[[renamed_to]] == df[[id]])] <- NA_character_
  df %>%
    distinct(kingdom, rank, fullname, .keep_all = TRUE) %>%
    select(-.kept_id)
}

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
  # keep only one GBIF taxon ID per full name, preferably the accepted one
  one_record_per_name("gbif")

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
  # keep only one LPSN record ID per full name, preferably the accepted one
  one_record_per_name("lpsn")

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
  # keep only one MycoBank record ID per full name, preferably the accepted one
  one_record_per_name("mycobank")


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
# (only names that are completely absent now, so that a taxon that moved to another domain, such as
# Trichophyton in 2026, does not come back in its old domain)
source_names <- names_in_sources()
source_names_any_domain <- sub("^[^ ]+ ", "", source_names)
taxonomy.old <- microorganisms_old %>%
  microsporidia_to_fungi() %>%
  select(any_of(colnames(taxonomy))) %>%
  filter(
    !fullname %in% taxonomy$fullname,
    # these will be added later:
    tolower(source) != "manually added",
    # not the (sub)species of which the species now only exists in another domain, as they are other organisms with
    # the same genus name (e.g. the subspecies of the butterfly Graphium sarpedon, which were in the Fungi until v3.0.1)
    !(rank %in% c("species", "subspecies") &
      paste(genus, species) %in% source_names_any_domain &
      !paste(domain, genus, species) %in% source_names),
    # and of the genera in `released_homonym_genera`, only the names that a current source still has in this domain
    !(rank %in% c("species", "subspecies") & genus %in% names(released_homonym_genera) &
      !paste(domain, fullname) %in% source_names)
  ) %>%
  # an old record outside the scope of its source is not that source's record (anymore), such as the microsporidian
  # families Culicosporidae and Janacekiidae, which were in the Protozoa with source MycoBank until v3.0.1
  mutate(source = if_else(
    (source == "MycoBank" & domain != "Fungi") | (source == "LPSN" & !domain %in% c("Bacteria", "Archaea")),
    "manually added",
    source
  ))
taxonomy <- taxonomy %>%
  bind_rows(taxonomy.old) %>%
  arrange(fullname) %>%
  filter(fullname != "")
rm(source_names, source_names_any_domain)

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
# and Animalia, or Trichophyton in Fungi and as a cyanobacterium in GBIF).
# The winning domain per genus is decided in this order:
# 1. Our own overrides in `genus_domain_override` (see 'Clinically relevant genera' above).
# 2. The domain of the genus in the last release, if it was in only one domain there, so that released codes
#    keep denoting the same organisms (see the review 'Homonyms decided by the domain of the last release'),
#    but only if an authoritative source supports that domain, or if no domain has an authoritative source.
# 3. Authoritative sources (LPSN, MycoBank) always beat GBIF and other sources,
#    so GBIF can never overrule LPSN or MycoBank. MycoBank records of a genus without any phylum are no
#    authoritative fungal records (MycoBank also indexes e.g. Plasmodium and Sarcocystis).
# 4. Within LPSN and MycoBank: accepted beats non-accepted, based on the genus-rank
#    record of that source/domain (or its species if the genus record is missing).
#    So a synonymic bacterial genus in LPSN loses from an accepted fungal genus in
#    MycoBank with the same name, and vice versa, unless rule 2 applies (such as for Pirella).
# 5. Source: LPSN > MycoBank > GBIF > other.
# 6. Genus-rank records over species/subspecies records.
# 7. Domain: Bacteria > Fungi > Protozoa > Archaea > Chromista > Animalia >
#    Plantae, which only breaks ties within the same source (e.g. Giardia within GBIF).
# The lineage of the winning domain then overwrites domain-to-family for the entire genus.
# Records of the losing domain are:
# - moved to the winning domain if they are protists that were placed in another kingdom (GBIF places
#   e.g. some Trypanosoma, Giardia and Amoeba species in Animalia), or if they have the same phylum or
#   class as the winning lineage (e.g. microsporidia placed in Protozoa instead of Fungi);
# - removed in all other cases, as they are other organisms with the same genus name (e.g. the
#   Graphium butterflies vs the Graphium moulds, or the Pirella fungi vs the released Pirella bacteria).
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
# the domain of every genus in the last release, if it was in only one domain there
released_genus_domain <- microorganisms_old %>%
  microsporidia_to_fungi() %>%
  filter(rank == "genus", domain %in% names(domain_priority)) %>%
  distinct(genus, domain) %>%
  group_by(genus) %>%
  filter(n() == 1) %>%
  ungroup()
released_genus_domain <- setNames(released_genus_domain$domain, released_genus_domain$genus)
first_nonempty <- function(x) {
  coalesce(x[!x %in% c("", NA)][1], "")
}

genus_candidates <- taxonomy %>%
  filter(
    genus != "",
    rank %in% c("genus", "species", "subspecies"),
    domain %in% names(domain_priority),
    # a source only counts within its own scope, so that a lineage inherited
    # from an older (possibly flawed) version of `microorganisms` cannot vote
    !(source == "LPSN" & !domain %in% c("Bacteria", "Archaea")),
    !(source == "MycoBank" & domain != "Fungi")
  ) %>%
  group_by(genus, source, domain) %>%
  mutate(
    # status of this genus in this source/domain: from the genus-rank record if
    # available, otherwise accepted if any of its species/subspecies is accepted
    genus_status = coalesce(
      first(status[rank == "genus"]),
      if_else(any(status == "accepted"), "accepted", first(status))
    ),
    # MycoBank also indexes names of organisms that it does not classify as fungi, without any phylum, such as
    # Plasmodium Marchiafava & Celli, Sarcocystis and Pseudospirillum: these are no authoritative fungal records
    authoritative = source == "LPSN" | (source == "MycoBank" & any(phylum != ""))
  ) %>%
  group_by(genus) %>%
  mutate(
    released = unname(released_genus_domain[genus]),
    # the released domain only decides if an authoritative source supports it, or if no domain has an
    # authoritative source (so a single GBIF record of a released genus cannot beat e.g. hundreds of MycoBank fungi),
    # and never on the basis of unclassified MycoBank records alone (e.g. the labyrinthulids, which are in Chromista)
    released_decides = !is.na(released) &
      any(domain == released & (authoritative | source != "MycoBank")) &
      (any(authoritative & domain == released) | !any(authoritative))
  ) %>%
  ungroup() %>%
  mutate(
    override = unname(genus_domain_override[genus]),
    oprio = case_when(
      is.na(override) ~ 1L,
      domain == override ~ 0L,
      TRUE ~ 2L
    ),
    # the domain of the genus in the last release, so that its released codes keep denoting the same organisms
    # (e.g. the bacterial Moorella, Serpula and Chainia, now synonyms in LPSN, may not lose to an accepted fungal
    # genus with the same name, as their codes would then be retired)
    relprio = case_when(
      !released_decides ~ 1L,
      domain == released ~ 0L,
      TRUE ~ 2L
    ),
    tprio = if_else(authoritative, 1L, 2L),
    stprio = if_else(tprio == 2L | genus_status %in% "accepted", 1L, 2L),
    sprio = if_else(source == "MycoBank" & !authoritative, source_prio("GBIF") + 1L, source_prio(source)),
    rprio = if_else(rank == "genus", 1L, 2L),
    dprio = unname(domain_priority[domain])
  ) %>%
  arrange(genus, oprio, relprio, tprio, stprio, sprio, rprio, dprio)

# REVIEW: genera of which the domain of the last release won, while the other rules would have chosen another
# domain - check that the released domain is right, and extend `genus_domain_override` if it is not
genus_candidates %>%
  arrange(genus, oprio, tprio, stprio, sprio, rprio, dprio) %>%
  group_by(genus) %>%
  summarise(
    domain_without_release = first(domain),
    released = first(released),
    released_decides = first(released_decides),
    override = first(override)
  ) %>%
  filter(!is.na(released), is.na(override), domain_without_release != released) %>%
  review("Homonyms decided by the domain of the last release")

best_domain <- genus_candidates %>%
  group_by(genus) %>%
  mutate(best_domain = first(domain)) %>%
  # take every rank from the most reliable record that has a value for it, so that an empty family
  # in the winning record (e.g. Mucor in MycoBank) cannot erase the family from another source
  filter(domain == best_domain) %>%
  summarise(
    best_domain = first(best_domain),
    best_kingdom = first_nonempty(kingdom),
    best_phylum = first_nonempty(phylum),
    best_class = first_nonempty(class),
    best_order = first_nonempty(order),
    best_family = first_nonempty(family),
    .groups = "drop"
  )

# REVIEW: genus homonyms of clinically relevant genera, check the winning domain and extend
# `genus_domain_override` at the start of this script if it is wrong
genus_candidates %>%
  filter(genus %in% relevant_genera) %>%
  group_by(genus) %>%
  filter(n_distinct(domain) > 1) %>%
  group_by(genus, domain) %>%
  summarise(
    n = n(),
    sources = paste(sort(unique(source)), collapse = "+"),
    example = first(fullname[rank == "species"]),
    .groups = "drop"
  ) %>%
  left_join(best_domain %>% select(genus, best_domain), by = "genus") %>%
  mutate(in_override_list = genus %in% names(genus_domain_override)) %>%
  review("Homonyms of relevant genera (check `best_domain`)")

taxonomy <- taxonomy %>%
  left_join(best_domain, by = "genus") %>%
  mutate(
    .resolve = case_when(
      genus == "" | is.na(best_domain) | !domain %in% names(domain_priority) ~ "keep",
      domain == best_domain ~ "keep",
      # LPSN and MycoBank records of a losing homonym are other organisms under another
      # nomenclatural code (e.g. the bacterial Pirella species if the fungal Pirella wins)
      source %in% c("LPSN", "MycoBank") ~ "remove",
      # protists placed in another kingdom
      best_domain %in% c("Protozoa", "Chromista") & domain %in% c("Protozoa", "Chromista", "Animalia") ~ "move",
      # same phylum or class, so the same organisms in another kingdom (e.g. microsporidia as Protozoa
      # instead of Fungi, which share the class Microsporea)
      (phylum != "" & phylum == best_phylum) | (class != "" & class == best_class) ~ "move",
      TRUE ~ "remove"
    )
  )
review(
  taxonomy %>% filter(.resolve == "remove") %>% count(genus, domain, best_domain, source, sort = TRUE),
  "Records removed as homonyms of a genus in another domain"
)
# these are other organisms with the same genus name, so they may never be restored as released taxa later on
# (e.g. the Graphium butterflies, which were in the Fungi until v3.0.1), saved for 'Restore released taxa'
saveRDS(
  taxonomy %>%
    filter(.resolve == "remove") %>%
    distinct(fullname, .keep_all = TRUE) %>%
    transmute(fullname, domain, genus, best_domain),
  "data-raw/taxonomy_homonyms_removed.rds",
  version = 2
)
taxonomy <- taxonomy %>%
  filter(.resolve != "remove") %>%
  mutate(
    is_genus_rank = genus != "" & !is.na(best_domain) & rank %in% c("genus", "species", "subspecies"),
    domain  = if_else(genus != "" & !is.na(best_domain), best_domain, domain),
    kingdom = if_else(genus != "" & !is.na(best_kingdom) & best_kingdom != "", best_kingdom, kingdom),
    phylum  = if_else(is_genus_rank & best_phylum != "", best_phylum, phylum),
    class   = if_else(is_genus_rank & best_class != "", best_class, class),
    order   = if_else(is_genus_rank & best_order != "", best_order, order),
    family  = if_else(is_genus_rank & best_family != "", best_family, family)
  ) %>%
  select(-starts_with("best_"), -is_genus_rank, -.resolve)
rm(genus_candidates, best_domain)

# integrity tests: no LPSN record may have left the prokaryotes and no MycoBank
# record may have left the fungi (see issue #309, Trichophyton as cyanobacterium).
# Both must return 0 rows - otherwise check the genus resolution above.
taxonomy %>%
  filter(source == "LPSN", !domain %in% c("Bacteria", "Archaea", "", NA)) %>%
  count(domain, genus, sort = TRUE)
taxonomy %>%
  filter(source == "MycoBank", !domain %in% c("Fungi", "", NA)) %>%
  count(domain, genus, sort = TRUE)
# and these known cross-domain homonyms must end up in the right domain
taxonomy %>%
  filter(
    rank == "genus",
    genus %in% c("Trichophyton", "Pirella", "Tubercularia", "Necator", "Capillaria", "Graphium", "Trypanosoma")
  ) %>%
  select(fullname, domain, kingdom, phylum, family, source)

# check if some genera within kingdoms have multiple families / orders, etc:
taxonomy %>%
  filter(genus != "") %>%
  group_by(kingdom, genus) %>%
  filter(n_distinct(family) > 1) %>%
  ungroup() %>%
  review("Genera with multiple families")
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
  # keep the record of the most authoritative source: LPSN > MycoBank > GBIF > other
  mutate(sprio = source_prio(source)) %>%
  arrange(fullname, sprio) %>%
  distinct(fullname, .keep_all = TRUE) %>%
  select(-sprio)


# *** Save intermediate results (0) *** -----------------------------------------------------------

saveRDS(taxonomy, "data-raw/taxonomy0.rds")
# taxonomy <- readRDS("data-raw/taxonomy0.rds")


# Add missing taxonomic entries and deduplicate ---------------------------------------------------
# Ensure every referenced rank has its own row (domain through species).
# Where possible, enrich with GBIF identifiers.

# @resume-block current_gbif after=taxonomy1
# (needed by add_missing_parents() after filtering, when resuming from a checkpoint)
current_gbif <- taxonomy_gbif.bak %>%
  filter(is.na(acceptedNameUsageID)) %>%
  mutate(
    taxonID = as.character(taxonID),
    parentNameUsageID = as.character(parentNameUsageID),
    domain = case_when(
      kingdom %in% na.omit(taxonomy_lpsn$kingdom[taxonomy_lpsn$domain == "Bacteria"]) ~ "Bacteria",
      kingdom %in% na.omit(taxonomy_lpsn$kingdom[taxonomy_lpsn$domain == "Archaea"]) ~ "Archaea",
      TRUE ~ kingdom
    ),
    # the own rank of a record is not always filled in, so take it from the name (as for `taxonomy_gbif`)
    canonical = case_when(
      is.na(scientificNameAuthorship) | scientificNameAuthorship == "" ~ scientificName,
      endsWith(scientificName, scientificNameAuthorship) ~
        trimws(substr(scientificName, 1L, nchar(scientificName) - nchar(scientificNameAuthorship))),
      TRUE ~ scientificName
    ),
    genus = if_else(taxonRank == "genus", canonical, genus),
    family = if_else(taxonRank == "family", canonical, family),
    order = if_else(taxonRank == "order", canonical, order),
    class = if_else(taxonRank == "class", canonical, class),
    phylum = if_else(taxonRank == "phylum", canonical, phylum)
  ) %>%
  microsporidia_to_fungi() %>%
  select(-canonical)
# @end-resume-block

taxonomy <- add_missing_parents(taxonomy, current_gbif)

# remove the old kingdoms Bacteria and Archaea, these are domains
taxonomy <- taxonomy %>%
  filter(!(rank == "kingdom" & fullname %in% c("Bacteria", "Archaea")))

# Deduplicate: where the same fullname appears at multiple ranks within a
# domain (e.g. Nitrospira as genus and class), keep the lowest rank as-is
# and disambiguate higher ranks by appending {rank}.
rank_priority <- c(
  "subspecies" = 1L, "species" = 2L, "genus" = 3L,
  "family" = 4L, "order" = 5L, "class" = 6L,
  "phylum" = 7L, "kingdom" = 8L, "domain" = 9L
)

taxonomy <- taxonomy %>%
  mutate(
    rank_index = rank_priority[rank],
    source_index = source_prio(source)
  ) %>%
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

# (note: until 2026, also records that were only inferred by this script were labelled 'manually added',
# so these are not all curated entries - they follow the relevance filters in 'Remove unwanted taxonomic entries')
homonyms_removed <- readRDS("data-raw/taxonomy_homonyms_removed.rds")
source_names <- names_in_sources()
manually_added <- microorganisms_old %>%
  microsporidia_to_fungi() %>%
  filter(
    tolower(source) == "manually added",
    !paste(kingdom, fullname) %in% paste(taxonomy$kingdom, taxonomy$fullname),
    !rank %in% c("domain", "kingdom", "phylum", "class", "order", "family"),
    # not the names that were removed above as another organism with the same genus name, unless the name still
    # exists in this domain in the sources (e.g. the Graphium butterflies, which were 'manually added' as Fungi
    # until v3.0.1, while the fungal genus Trichurus is a fungus indeed)
    !(fullname %in% homonyms_removed$fullname & !paste(domain, fullname) %in% source_names),
    # and of the genera in `released_homonym_genera`, only the names that a current source still has in this domain
    !(genus %in% names(released_homonym_genera) & !paste(domain, fullname) %in% source_names)
  ) %>%
  select(fullname:subspecies, ref, source, rank, status)

# Build a lookup for each child->parent rank pair, preferring LPSN > MycoBank > any
fill_from_parent <- function(data, taxonomy, child_rank, parent_rank) {
  lookup <- taxonomy %>%
    filter(.data[[child_rank]] != "", .data[[parent_rank]] != "") %>%
    mutate(sprio = source_prio(source)) %>%
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
    # keep a curated 'accepted' status (such as of Plasmodium), all others are unknown
    status = if_else(status %in% "accepted", status, "unknown"),
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
    taxonomy$source %in% c("GBIF", "manually added", "inferred") &
    (taxonomy$rank %in% c("phylum", "class", "order", "family", "genus") |
       (taxonomy$rank == "species" & taxonomy$genus %in% lpsn_genera))
)
message(length(gbif_bacteria))
added <- 0
pb <- progress_bar$new(
  total = length(gbif_bacteria),
  format = "[:bar] :current/:total :eta"
)
for (i in seq_along(gbif_bacteria)) {
  pb$tick()
  record <- gbif_bacteria[i]
  if (i %% 500 == 0) {
    save_lpsn_cache()
  }
  lpsn <- get_lpsn_and_author(
    rank = taxonomy$rank[record],
    name = taxonomy$fullname[record]
  )
  if (is.na(lpsn["lpsn"])) {
    next
  } else {
    added <- added + 1
    taxonomy <- apply_lpsn_result(taxonomy, record, lpsn)
  }
}
save_lpsn_cache()
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
      # an empty parent in the authoritative source must never erase a known parent from GBIF
      # (e.g. MycoBank has no family for Mucor, GBIF has Mucoraceae)
      filter(!.data[[parent_col]] %in% c("", NA)) %>%
      # prefer the record of the child itself, e.g. the family record for the order of a family
      arrange(rank != child_col) %>%
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
# Names that are not validly published (according to LPSN) or could not be checked are removed, unless
# they are clinically relevant. E.g., Tropheryma whipplei is 'not validly published' under the ICNP, and
# was therefore removed from this data set in 2026.
# For prokaryotes, only the human pathogens of Bartlett et al. (and their genera) are protected, otherwise all
# non-validly published species of e.g. Streptococcus would return.
protected <- (taxonomy$domain %in% c("Bacteria", "Archaea") &
  (paste(taxonomy$genus, taxonomy$species) %in% paste(pathogens$genus, pathogens$species) | is_who_priority_genus(taxonomy$genus) |
    (taxonomy$rank == "genus" & taxonomy$genus %in% pathogens$genus))) |
  (!taxonomy$domain %in% c("Bacteria", "Archaea") &
    (taxonomy$genus %in% relevant_genera | taxonomy$fullname %in% relevant_current_species))
review(
  taxonomy %>% filter(!status %in% c("accepted", "synonym"), protected, !source %in% c("manually added", "inferred")),
  "Not validly published or unknown status, but kept as clinically relevant"
)
taxonomy <- taxonomy %>%
  mutate(status = if_else(!status %in% c("accepted", "synonym") & protected & source != "inferred", "accepted", status)) %>%
  filter(status %in% c("accepted", "synonym") | source %in% c("manually added", "inferred"))
rm(protected)

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
  arrange(fullname, source_prio(source), lpsn, mycobank, desc(gbif)) %>%
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
# (`pathogens` was read at the start of this script)

# @resume-block prevalence after=taxonomy2
# (compute_prevalence() is used again later on, when resuming from a checkpoint)
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
    current_names_same_domain(.),
    unlist(mo_synonyms(intersect(., microorganisms_dev$fullname), keep_synonyms = FALSE))
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
# this is a function, since it must be run again after the higher taxonomy has been harmonised below
compute_prevalence <- function(taxonomy) {
  taxonomy <- taxonomy %>%
    mutate(
      who_genus = genus %in% MO_WHO_PRIORITY_GENERA,
      prevalence = case_when(
        # 'established' means 'have infected at least three persons in three or more references'
        paste(genus, species) %in% established & rank %in% c("species", "subspecies") ~ 1.15,
        # other genera in the 'established' group
        genus %in% established_genera & rank == "genus" ~ 1.15,
        # 'putative' means 'fewer than three known cases'
        paste(genus, species) %in% putative & rank %in% c("species", "subspecies") ~ 1.25,
        # other genera in the 'putative' group
        genus %in% putative_genera & rank == "genus" ~ 1.25,
        # species and subspecies in 'established' and 'putative' groups
        genus %in% c(established_genera, putative_genera) & rank %in% c("species", "subspecies") ~ 1.5,
        # non-bacterial genera/species/subspecies of clinical relevance
        genus %in% nonbacterial_genera & domain != "Bacteria" & rank %in% c("genus", "species", "subspecies") ~ 1.25,
        # all others
        TRUE ~ 2.0
      )
    )
  # a current name is as relevant as its synonyms, so that e.g. Candidozyma auris (renamed from
  # Candida auris in 2024) gets the prevalence of Candida auris, without having to update any list.
  # This uses the prevalence without the WHO rule below, since that rule applies to whole genera:
  # otherwise every species that ever was in e.g. Pseudomonas would become 1.0.
  for (id in c("lpsn", "mycobank", "gbif")) {
    renamed_col <- paste0(id, "_renamed_to")
    if (!all(c(id, renamed_col) %in% colnames(taxonomy))) {
      next
    }
    via_synonyms <- taxonomy %>%
      filter(!is.na(.data[[renamed_col]])) %>%
      group_by(target = .data[[renamed_col]]) %>%
      summarise(syn_prevalence = min(prevalence), .groups = "drop")
    taxonomy <- taxonomy %>%
      left_join(via_synonyms, by = setNames("target", id), na_matches = "never") %>%
      mutate(prevalence = pmin(prevalence, syn_prevalence, na.rm = TRUE)) %>%
      select(-syn_prevalence)
  }
  # genera of the pathogens mentioned in the World Health Organization's (WHO) Priority Pathogen List
  taxonomy <- taxonomy %>%
    mutate(prevalence = if_else(who_genus, 1.0, prevalence)) %>%
    select(-who_genus)
  # and a genus is as relevant as its most relevant accepted species
  genus_prevalence <- taxonomy %>%
    filter(rank == "species", status == "accepted", genus != "") %>%
    group_by(domain, genus) %>%
    summarise(gen_prevalence = min(prevalence), .groups = "drop")
  taxonomy %>%
    left_join(genus_prevalence, by = c("domain", "genus")) %>%
    mutate(prevalence = if_else(rank == "genus", pmin(prevalence, gen_prevalence, na.rm = TRUE), prevalence)) %>%
    select(-gen_prevalence)
}
# @end-resume-block
taxonomy <- compute_prevalence(taxonomy)

table(taxonomy$prevalence, useNA = "always")
# (a lot will be removed further below)


# *** Save intermediate results (2) *** -----------------------------------------------------------

saveRDS(taxonomy, "data-raw/taxonomy2.rds")
# taxonomy <- readRDS("data-raw/taxonomy2.rds")


# Remove unwanted taxonomic entries ---------------------------------------------------------------

part1 <- taxonomy %>%
  filter(
    # keep all unknowns we added ourselves (other 'manually added' records follow the
    # rules below, since until 2026 also records that were only inferred were labelled 'manually added')
    fullname %like% "unknown" |
      # keep all bacteria anyway, main focus of our package:
      domain == "Bacteria" |
      # these domains are very small, and also mainly microorganisms (Protozoa: not anymore since the COL eXtended
      # Release, so they follow the same rules as the Fungi below, decision by Matthijs S. Berends, 5 October 2026):
      (domain == "Protozoa" & (!rank %in% c("genus", "species", "subspecies") | prevalence < 2)) |
      domain == "Archaea" |
      domain == "Chromista" |
      # keep everything from family up, 'ghost' entries will be removed later on
      rank %in% c("domain", "kingdom", "phylum", "class", "order", "family") |
      # other relevant genera to keep (see 'Clinically relevant genera' above):
      genus %in% relevant_genera |
      # and the current names of their species, with the record of their genus:
      fullname %in% relevant_current_species |
      (rank == "genus" & fullname %in% relevant_current_genera) |
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
    # species groups are added again at the end from the previous data set, and would otherwise get the
    # MO code of their genus below (e.g. B_MYCBC for 'Mycobacterium avium-intracellulare complex')
    rank != "species group",
    !(genus %in% c("Aedes", "Anopheles") & rank %in% c("species", "subspecies"))
  )

# now get the parents, the current names of kept synonyms, and the synonyms of kept (sub)species,
# and repeat until nothing is added anymore (previously this was done a fixed number of 3 times)
parts <- part1
repeat {
  ids_wanted <- function(id) {
    na.omit(c(parts[[paste0(id, "_parent")]], parts[[paste0(id, "_renamed_to")]]))
  }
  kept_species <- parts %>% filter(rank %in% c("species", "subspecies"))
  added <- taxonomy %>%
    filter(
      !fullname %in% parts$fullname,
      # parents and current names
      gbif %in% ids_wanted("gbif") |
        mycobank %in% ids_wanted("mycobank") |
        lpsn %in% ids_wanted("lpsn") |
        # old names of kept (sub)species, so that as.mo() can still find e.g. old fungal names
        (status == "synonym" & rank %in% c("species", "subspecies") &
          (gbif_renamed_to %in% na.omit(kept_species$gbif) |
            mycobank_renamed_to %in% na.omit(kept_species$mycobank) |
            lpsn_renamed_to %in% na.omit(kept_species$lpsn)))
    )
  message("Adding ", nrow(added), " parents, current names and synonyms")
  if (nrow(added) == 0) {
    break
  }
  parts <- bind_rows(parts, added)
}
rm(added, kept_species)

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

# No lichens (decision by Matthijs S. Berends, 5 October 2026): these entered as current names of kept synonyms, with all
# their own synonyms, and are not relevant for this package. Lichens are taken as the
# lichen-forming classes and the order Verrucariales, unless their genus is clinically relevant. Synonyms that point to
# a removed lichen are removed as well. Released lichens are restored below, as released taxa are never removed.
lichens <- taxonomy %>%
  filter(
    domain == "Fungi",
    class %in% c("Arthoniomycetes", "Candelariomycetes", "Lecanoromycetes", "Lichinomycetes") | order == "Verrucariales",
    !genus %in% relevant_genera
  )
lichens <- bind_rows(
  lichens,
  taxonomy %>%
    filter(
      status == "synonym",
      (!is.na(gbif_renamed_to) & gbif_renamed_to %in% na.omit(lichens$gbif)) |
        (!is.na(mycobank_renamed_to) & mycobank_renamed_to %in% na.omit(lichens$mycobank))
    )
) %>%
  distinct(fullname, .keep_all = TRUE)
review(
  lichens %>% count(class, order, status, name = "n_records"),
  "Lichens and their synonyms, these will be removed"
)
taxonomy <- taxonomy %>%
  filter(!fullname %in% lichens$fullname)
rm(lichens)

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
    .ref_phylum     = consensus_value(phylum, rank, source, "genus"),
    .ref_class      = consensus_value(class, rank, source, "genus"),
    .ref_order      = consensus_value(order, rank, source, "genus"),
    .ref_family     = consensus_value(family, rank, source, "genus"),
    .groups = "drop"
  )

taxonomy <- taxonomy %>%
  left_join(genus_ref, by = "genus") %>%
  mutate(
    across_flag = rank != "genus" & genus != "" & fullname %unlike% "unknown",
    domain     = if_else(across_flag, .ref_domain, domain),
    kingdom    = if_else(across_flag, .ref_kingdom, kingdom),
    # also fill the genus record itself where it has no value
    phylum     = if_else((across_flag | (rank == "genus" & phylum == "")) & .ref_phylum != "", .ref_phylum, phylum),
    class      = if_else((across_flag | (rank == "genus" & class == "")) & .ref_class != "", .ref_class, class),
    order      = if_else((across_flag | (rank == "genus" & order == "")) & .ref_order != "", .ref_order, order),
    family     = if_else((across_flag | (rank == "genus" & family == "")) & .ref_family != "", .ref_family, family)
  ) %>%
  select(-starts_with(".ref_"), -across_flag) %>%
  arrange(fullname)
# (until 2026, the prevalence of the genus was also copied to all its species here, which made the
# species-level information of Bartlett et al. useless - the prevalence is now recalculated below)

# update the taxonomic names based on new genus classification
taxonomy <- taxonomy %>%
  # first, fix the higher taxonomy of families, preferring the family record itself, then the most
  # authoritative source (previously the first record in alphabetical order was taken, so that one outdated
  # GBIF record could change the phylum of a whole family)
  group_by(domain, family) %>%
  mutate(
    kingdom = if_else(family != "", consensus_value(kingdom, rank, source, "family"), kingdom),
    phylum = if_else(family != "", consensus_value(phylum, rank, source, "family"), phylum),
    class = if_else(family != "", consensus_value(class, rank, source, "family"), class),
    order = if_else(family != "", consensus_value(order, rank, source, "family"), order)
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

# recalculate the prevalence, now that all species have the taxonomy of their genus
taxonomy <- compute_prevalence(taxonomy)

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

# the filters above may have removed parents of records that were kept (e.g. the genus of a kept synonym),
# which would leave these records without a valid MO code later on - so add them again
n_before <- nrow(taxonomy)
taxonomy <- add_missing_parents(taxonomy, current_gbif) %>%
  filter(!(rank == "kingdom" & fullname %in% c("Bacteria", "Archaea")))
message("Added ", nrow(taxonomy) - n_before, " missing parents after filtering")
review(
  taxonomy %>% filter(duplicated(fullname) | duplicated(fullname, fromLast = TRUE)),
  "Duplicate full names after adding missing parents"
)
taxonomy <- taxonomy %>%
  arrange(fullname, source_prio(source)) %>%
  distinct(fullname, .keep_all = TRUE) %>%
  # this also sets the prevalence of the added records
  compute_prevalence()


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

# MB 2026-10-03/ the same applies to Leishmania (in 2026, COL contained only a handful of its species,
# not L. donovani, L. major or L. aethiopica), so add the clinically relevant species that are missing
leishmania_genus <- taxonomy %>% filter(fullname == "Leishmania", rank == "genus")
if (nrow(leishmania_genus) == 1) {
  leishmania <- tibble(
    species = c(
      "aethiopica", "amazonensis", "braziliensis", "donovani", "guyanensis", "infantum", "lainsoni",
      "major", "martiniquensis", "mexicana", "naiffi", "panamensis", "peruviana", "shawi", "tropica",
      "venezuelensis"
    )
  ) %>%
    mutate(
      fullname = paste("Leishmania", species),
      domain = leishmania_genus$domain,
      kingdom = leishmania_genus$kingdom,
      phylum = leishmania_genus$phylum,
      class = leishmania_genus$class,
      order = leishmania_genus$order,
      family = leishmania_genus$family,
      genus = "Leishmania",
      subspecies = "",
      rank = "species",
      source = "manually added",
      status = "accepted",
      prevalence = leishmania_genus$prevalence
    ) %>%
    filter(!fullname %in% taxonomy$fullname)
  message("Adding ", nrow(leishmania), " Leishmania species")
  taxonomy <- taxonomy %>%
    bind_rows(leishmania) %>%
    arrange(fullname)
} else {
  warning("Leishmania genus not found, check this!", call. = FALSE)
}

# MB 2026-10-05/ in 2026, COL contained no species of Brugia at all, so add the species that cause human lymphatic
# filariasis (with the lineage of the genus)
brugia_genus <- taxonomy %>% filter(fullname == "Brugia", rank == "genus")
if (nrow(brugia_genus) == 1) {
  brugia <- tibble(species = c("malayi", "timori")) %>%
    mutate(
      fullname = paste("Brugia", species),
      domain = brugia_genus$domain,
      kingdom = brugia_genus$kingdom,
      phylum = brugia_genus$phylum,
      class = brugia_genus$class,
      order = brugia_genus$order,
      family = brugia_genus$family,
      genus = "Brugia",
      subspecies = "",
      rank = "species",
      source = "manually added",
      status = "accepted",
      prevalence = brugia_genus$prevalence
    ) %>%
    filter(!fullname %in% taxonomy$fullname)
  message("Adding ", nrow(brugia), " Brugia species")
  taxonomy <- taxonomy %>%
    bind_rows(brugia) %>%
    arrange(fullname)
  rm(brugia)
} else {
  warning("Brugia genus not found, check this!", call. = FALSE)
}
rm(brugia_genus)

# MB 2026-10-05/ COL does not contain Cystoisospora belli at all, only Isospora belli, so add it as a synonym of
# Isospora belli, so that as.mo() can find it (decision by Matthijs S. Berends, 5 October 2026)
isospora_belli <- taxonomy %>% filter(fullname == "Isospora belli", rank == "species")
if (nrow(isospora_belli) == 1 && !"Cystoisospora belli" %in% taxonomy$fullname) {
  taxonomy <- taxonomy %>%
    bind_rows(
      isospora_belli %>%
        mutate(
          fullname = "Cystoisospora belli",
          genus = "Cystoisospora",
          status = "synonym",
          source = "manually added",
          ref = "",
          across(any_of(c("lpsn", "lpsn_parent", "lpsn_renamed_to", "mycobank", "mycobank_parent", "mycobank_renamed_to",
                          "gbif", "gbif_parent")), ~NA_character_),
          gbif_renamed_to = isospora_belli$gbif
        )
    ) %>%
    arrange(fullname)
  message("Added Cystoisospora belli as synonym of Isospora belli")
} else if (nrow(isospora_belli) != 1) {
  warning("Isospora belli not found, check this!", call. = FALSE)
}
rm(isospora_belli)


# Clean-up of nonsense names ----------------------------------------------------------------------

# (until 2026, this was applied to `microorganisms` instead of `taxonomy`, so it had no effect, and it was
# done after the MO codes were made - names such as 'Verticillium sbs.' then got a code like F_VRTCL_SBS.)
table(unlist(strsplit(taxonomy$fullname, "")))
n_before <- nrow(taxonomy)
taxonomy <- taxonomy %>%
  filter(fullname %unlike% " [0-9]") %>%
  filter(fullname %unlike% " (i|ii|iii|iv|v|vi|vii|viii|ix|x)$" | genus == "Streptococcus") %>%
  filter(fullname %unlike% "[.]") %>%
  # ASCII names without quotes, also for the records from the previous data set
  mutate(across(c(fullname, domain:subspecies), ascii_names)) %>%
  distinct(fullname, .keep_all = TRUE)
message("Removed ", n_before - nrow(taxonomy), " nonsense names")

# Restore released taxa ---------------------------------------------------------------------------

# A taxon that was part of a release (since v2.0.0) is never removed, since users have stored its MO code
# (decision by Matthijs S. Berends, 3 October 2026). Taxa that are missing now, e.g. because they are not
# validly published or not clinically relevant, are restored from their last release, with their old code and
# their last known status. The only exception: records that were removed as homonyms above, i.e. other organisms
# with the same genus name (such as the Graphium butterflies, which were wrongly in the Fungi until v3.0.1).
# Their codes are retired, with the reason, in data-raw/microorganisms_files/mo_code_retirements.csv, and are
# translated to NA by as.mo() with that reason, so that they never lead to a wrong organism and are never reused.
restore_registry <- read_mo_registry(".")
restore_renames <- read_mo_renames(".")
restore_registry$expected_name <- restore_registry$fullname
renamed <- match(restore_registry$mo, restore_renames$mo)
restore_registry$expected_name[!is.na(renamed)] <- restore_renames$new_name[renamed[!is.na(renamed)]]
missing_released <- restore_registry %>%
  filter(
    # by name and rank, since a name can be used at two ranks (e.g. the genus Kapabacteria and the class
    # 'Kapabacteria {class}'): if only one of them is still present, the other one is missing too (not by domain,
    # since a code of a taxon that moved to another domain, such as P_TXPL for Toxoplasma, is translated by its name)
    !paste(rank, mo_name_without_suffix(expected_name)) %in%
      paste(taxonomy$rank, mo_name_without_suffix(taxonomy$fullname)),
    # species groups are added at the end of this script
    rank != "species group"
  )
homonyms_removed <- if (file.exists("data-raw/taxonomy_homonyms_removed.rds")) {
  readRDS("data-raw/taxonomy_homonyms_removed.rds")
} else {
  tibble(fullname = character(0), domain = character(0), genus = character(0), best_domain = character(0))
}

# codes to retire: released taxa that are other organisms than the genus with their name
retirements_file <- "data-raw/microorganisms_files/mo_code_retirements.csv"
retirements <- utils::read.csv(retirements_file, colClasses = "character", na.strings = character(0))
source_names <- names_in_sources()
new_retirements <- missing_released %>%
  inner_join(
    homonyms_removed %>% select(fullname, removed_domain = domain, removed_genus = genus, best_domain),
    by = c("expected_name" = "fullname")
  ) %>%
  filter(
    !mo %in% retirements$mo,
    # a released taxon in the winning domain of which the name still exists in that domain in the sources is the
    # same organism, not another one (e.g. the fungal genus Trichurus, while 'Trichurus' in Animalia is a misspelling
    # of the whipworm genus Trichuris), so it is restored below instead of retired
    !(domain == best_domain & paste(domain, expected_name) %in% source_names)
  ) %>%
  transmute(
    mo,
    fullname = expected_name,
    reason = paste0(
      "another organism (", removed_domain, ") with the same genus name as ", removed_genus, " (", best_domain, ")"
    ),
    decided_by = "taxonomy build (automatic, homonym of a genus in another domain)",
    date = as.character(Sys.Date())
  )
# released (sub)species of the genera in `released_homonym_genera` that no current source has in their domain cannot be
# told apart from the other organisms with the same genus name, so they are not restored but retired, to be reviewed
new_retirements <- new_retirements %>%
  bind_rows(
    missing_released %>%
      filter(
        genus %in% names(released_homonym_genera),
        rank %in% c("species", "subspecies"),
        domain == unname(released_homonym_genera[genus]),
        !paste(domain, expected_name) %in% source_names,
        !mo %in% c(retirements$mo, new_retirements$mo)
      ) %>%
      transmute(
        mo,
        fullname = expected_name,
        reason = paste0(
          "released in the genus ", genus, " (", domain, "), which also contained other organisms with the same genus ",
          "name until v3.0.1, but no current source has this name in the ", domain
        ),
        decided_by = "taxonomy build (automatic, released_homonym_genera), NEEDS REVIEW",
        date = as.character(Sys.Date())
      )
  )
review(new_retirements, "Released MO codes that are retired, since their taxon is another organism with the same genus name")
if (nrow(new_retirements) > 0) {
  retirements <- bind_rows(retirements, new_retirements) %>% arrange(mo)
  utils::write.csv(retirements, retirements_file, row.names = FALSE, na = "")
}
missing_released <- missing_released %>%
  filter(!mo %in% retirements$mo)

rank_order <- c(
  "subspecies" = 1L, "species" = 2L, "genus" = 3L, "family" = 4L, "order" = 5L, "class" = 6L,
  "phylum" = 7L, "kingdom" = 8L, "domain" = 9L
)

# get the full records from the releases in which these taxa were last present
read_release_data <- function(tag) {
  tmp <- tempfile(fileext = ".rda")
  on.exit(unlink(tmp))
  if (!identical(system2("git", c("show", paste0(tag, ":data/microorganisms.rda")), stdout = tmp), 0L)) {
    stop("Could not read data/microorganisms.rda of release ", tag, call. = FALSE)
  }
  env <- new.env()
  load(tmp, envir = env)
  out <- as.data.frame(env$microorganisms, stringsAsFactors = FALSE)
  out$mo <- as.character(out$mo)
  # until v3.0.1, `kingdom` contained what is now called `domain`
  if (!"domain" %in% colnames(out)) {
    out$domain <- out$kingdom
    out$domain[out$domain == "(unknown kingdom)"] <- "(unknown domain)"
    out$kingdom <- NA_character_
  }
  out
}
released_records <- bind_rows(lapply(unique(missing_released$last_release), function(tag) {
  read_release_data(tag) %>%
    filter(mo %in% missing_released$mo[missing_released$last_release == tag]) %>%
    select(any_of(c(
      "mo", "fullname", "status", "domain", "kingdom", "phylum", "class", "order", "family", "genus", "species",
      "subspecies", "rank", "ref", "lpsn_renamed_to", "mycobank_renamed_to", "gbif_renamed_to", "prevalence"
    ))) %>%
    mutate(across(c(fullname, domain:subspecies), ascii_names))
}))

# a name can be registered in two domains (e.g. the genus Octospora as a fungus until v3.0.1 and as a microsporidian
# until v2.1.1), but names must be unique: only the one of the most recent release is restored, the codes of the older
# ones are retired, as they cannot be restored; whether they denote the same organism cannot be decided automatically
older_homonyms <- released_records %>%
  left_join(missing_released %>% select(mo, last_release), by = "mo") %>%
  mutate(release_version = package_version(sub("^v", "", last_release))) %>%
  group_by(fullname, rank) %>%
  filter(n_distinct(domain) > 1) %>%
  arrange(desc(release_version), .by_group = TRUE) %>%
  mutate(
    newer_domain = first(domain),
    # a microsporidian in both releases is the same organism, now in the Fungi (see is_microsporidian())
    both_microsporidian = all(is_microsporidian(phylum, class))
  ) %>%
  filter(domain != newer_domain) %>%
  ungroup()
# the older records of the same microsporidia are not restored, their codes are translated by name
released_records <- released_records %>%
  filter(!mo %in% older_homonyms$mo[older_homonyms$both_microsporidian]) %>%
  microsporidia_to_fungi()
older_homonyms <- older_homonyms %>%
  filter(!both_microsporidian)
older_homonym_retirements <- older_homonyms %>%
  transmute(
    mo,
    fullname,
    reason = paste0(
      "a more recent release has its name in another domain (", newer_domain, "), which may be the same organism in ",
      "a new classification or another organism, and it could not be restored"
    ),
    decided_by = "taxonomy build (automatic, name in another domain in a more recent release), NEEDS REVIEW",
    date = as.character(Sys.Date())
  ) %>%
  filter(!mo %in% read_mo_retirements(".")$mo)
review(older_homonym_retirements, "Released MO codes that are retired, since a more recent release has their name in another domain")
if (nrow(older_homonym_retirements) > 0) {
  utils::write.csv(
    bind_rows(read_mo_retirements("."), older_homonym_retirements) %>% arrange(mo),
    retirements_file,
    row.names = FALSE, na = ""
  )
}
released_records <- released_records %>%
  filter(!mo %in% older_homonyms$mo)

# a released taxon of which the genus now only exists in another domain cannot be restored consistently: it may be
# the same organism in a new classification, or another organism (homonym), which cannot be decided automatically.
# Its code is retired (as.mo() then gives NA with this reason, and translates it if the name exists elsewhere).
genus_domains <- taxonomy %>%
  filter(rank == "genus") %>%
  distinct(genus, domain)
conflicting <- released_records %>%
  filter(
    genus != "",
    !paste(domain, genus) %in% paste(genus_domains$domain, genus_domains$genus),
    genus %in% genus_domains$genus
  ) %>%
  left_join(genus_domains %>% group_by(genus) %>% summarise(now_in = toString(domain)), by = "genus")
conflict_retirements <- conflicting %>%
  transmute(
    mo,
    fullname,
    reason = paste0("its genus ", genus, " is now only known in another domain (", now_in, ") and could not be restored"),
    decided_by = "taxonomy build (automatic, genus in another domain), NEEDS REVIEW",
    date = as.character(Sys.Date())
  ) %>%
  filter(!mo %in% read_mo_retirements(".")$mo)
review(conflict_retirements, "Released MO codes that are retired, since their genus is now only known in another domain")
if (nrow(conflict_retirements) > 0) {
  utils::write.csv(
    bind_rows(read_mo_retirements("."), conflict_retirements) %>% arrange(mo),
    retirements_file,
    row.names = FALSE, na = ""
  )
}
released_records <- released_records %>%
  filter(!mo %in% conflicting$mo)

# the current kingdom of prokaryotes (until v3.0.1, the kingdom was the domain) and the current higher taxonomy
# of the genus if it still exists
kingdom_of_phylum <- taxonomy %>%
  filter(domain %in% c("Bacteria", "Archaea"), phylum != "", kingdom != "") %>%
  count(domain, phylum, kingdom) %>%
  arrange(desc(n)) %>%
  distinct(domain, phylum, .keep_all = TRUE) %>%
  select(domain, phylum, new_kingdom = kingdom)
lineage_of_genus <- taxonomy %>%
  filter(rank == "genus") %>%
  distinct(domain, genus, .keep_all = TRUE) %>%
  select(domain, genus, g_kingdom = kingdom, g_phylum = phylum, g_class = class, g_order = order, g_family = family)
restored <- released_records %>%
  left_join(kingdom_of_phylum, by = c("domain", "phylum")) %>%
  left_join(lineage_of_genus, by = c("domain", "genus")) %>%
  mutate(
    kingdom = case_when(
      !is.na(g_kingdom) ~ g_kingdom,
      domain %in% c("Bacteria", "Archaea") & !is.na(new_kingdom) ~ new_kingdom,
      domain %in% c("Bacteria", "Archaea") ~ "(unknown kingdom)",
      TRUE ~ coalesce(kingdom, domain)
    ),
    phylum = coalesce(g_phylum, phylum),
    class = coalesce(g_class, class),
    order = coalesce(g_order, order),
    family = coalesce(g_family, family),
    across(kingdom:subspecies, function(x) if_else(is.na(x), "", x)),
    # older releases used placeholders such as "(unknown class)", which would become records themselves
    across(phylum:family, function(x) if_else(x %like% "^[(]unknown", "", x)),
    # older releases also had other status values, such as "not validly published"
    status = if_else(status %in% c("accepted", "synonym", "unknown"), status, "unknown"),
    source = "manually added"
  ) %>%
  select(-mo, -new_kingdom, -starts_with("g_"))
review(
  restored %>% count(domain, rank, status, sort = TRUE),
  "Released taxa that were missing and are restored (by domain, rank and status)"
)
review(
  restored %>% filter(prevalence < 2) %>% select(fullname, domain, rank, status, ref, prevalence),
  "Restored released taxa that are clinically relevant"
)
n_before <- nrow(taxonomy)
taxonomy <- taxonomy %>%
  bind_rows(restored) %>%
  add_missing_parents(current_gbif) %>%
  mutate(fullname = sub(" [{][a-z]+[}]$", "", fullname)) %>%
  # a name at two ranks: keep the lowest rank as-is and append {rank} to the higher one, as in 'Deduplicate' above
  # (the restored records can bring such a name back); also between domains, as names must be unique (e.g. the
  # fungal genus Acantharia and the radiolarian class 'Acantharia {class}', decision by Matthijs S. Berends, 7 October
  # 2026)
  group_by(fullname) %>%
  mutate(fullname = if_else(
    n_distinct(rank) > 1 & rank_order[rank] > min(rank_order[rank]),
    paste0(fullname, " {", rank, "}"),
    fullname
  )) %>%
  ungroup() %>%
  arrange(fullname, source_prio(source)) %>%
  distinct(fullname, .keep_all = TRUE) %>%
  compute_prevalence()
message("Restored ", nrow(restored), " released taxa (", nrow(taxonomy) - n_before, " records including their parents)")
rm(restore_registry, restore_renames, renamed, missing_released, homonyms_removed, retirements, new_retirements, source_names, older_homonyms, older_homonym_retirements,
   released_records, genus_domains, conflicting, conflict_retirements, kingdom_of_phylum, lineage_of_genus, restored,
   rank_order)


# Link orthographic variants without a current name ----------------------------------------------

# A synonym without a current name of which the epithet only differs from an accepted name in the same genus by a
# Latin gender ending (e.g. -us, -a, -um) is an orthographic variant of that name: under the nomenclatural codes,
# such epithets are variants of one name, so two valid species in one genus cannot differ only in this way.
# E.g. 'Nakaseomyces glabrata' -> 'Nakaseomyces glabratus' (from Candida glabrata, renamed in 2022).
epithet_stem <- function(x) {
  sub("(us|um|a|is|e|er|ra|rum)$", "", x)
}
variant_candidates <- taxonomy %>%
  filter(
    status == "synonym",
    is.na(lpsn_renamed_to), is.na(mycobank_renamed_to), is.na(gbif_renamed_to),
    rank %in% c("species", "subspecies")
  ) %>%
  mutate(stem = epithet_stem(if_else(rank == "species", species, subspecies)))
accepted_variants <- taxonomy %>%
  filter(status == "accepted", rank %in% c("species", "subspecies")) %>%
  mutate(stem = epithet_stem(if_else(rank == "species", species, subspecies))) %>%
  select(domain, genus, rank, stem, species_acc = species, target = fullname, t_lpsn = lpsn, t_mycobank = mycobank, t_gbif = gbif)
variant_links <- variant_candidates %>%
  inner_join(accepted_variants, by = c("domain", "genus", "rank", "stem")) %>%
  # for subspecies, the species must be the same variant too
  filter(rank == "species" | epithet_stem(species) == epithet_stem(species_acc)) %>%
  # only unambiguous links, and the target must have an identifier to link to
  group_by(fullname) %>%
  filter(n() == 1) %>%
  ungroup() %>%
  filter(!is.na(t_lpsn) | !is.na(t_mycobank) | !is.na(t_gbif))
review(variant_links %>% select(fullname, target, source), "Orthographic variants linked to their accepted name")
taxonomy <- taxonomy %>%
  left_join(variant_links %>% select(fullname, t_lpsn, t_mycobank, t_gbif), by = "fullname") %>%
  mutate(
    lpsn_renamed_to = coalesce(lpsn_renamed_to, t_lpsn),
    mycobank_renamed_to = coalesce(mycobank_renamed_to, t_mycobank),
    gbif_renamed_to = coalesce(gbif_renamed_to, t_gbif)
  ) %>%
  select(-t_lpsn, -t_mycobank, -t_gbif)
rm(variant_candidates, accepted_variants, variant_links)


# Break cycles of synonyms ------------------------------------------------------------------------

# Sources can contradict each other, e.g. in 2026 GBIF listed Capillidium denaeosporus as synonym of
# Conidiobolus denaeosporum and also the other way around, so that the 'current name' depended on the starting
# point. In a cycle, the statement of the most authoritative source counts: that record stays a synonym, and its
# target becomes accepted.
next_synonym_step <- function(df) {
  by_gbif <- df$fullname[match(df$gbif_renamed_to, df$gbif, incomparables = NA)]
  by_mycobank <- df$fullname[match(df$mycobank_renamed_to, df$mycobank, incomparables = NA)]
  by_lpsn <- df$fullname[match(df$lpsn_renamed_to, df$lpsn, incomparables = NA)]
  out <- coalesce(by_lpsn, by_mycobank, by_gbif)
  out[df$status != "synonym"] <- NA_character_
  setNames(out, df$fullname)
}
synonym_next <- next_synonym_step(taxonomy)
synonym_cycles <- list()
for (start in names(synonym_next)[!is.na(synonym_next)]) {
  path <- start
  repeat {
    nxt <- unname(synonym_next[path[length(path)]])
    if (is.na(nxt)) {
      break
    }
    if (nxt %in% path) {
      cycle <- sort(path[match(nxt, path):length(path)])
      synonym_cycles[[paste(cycle, collapse = " | ")]] <- cycle
      break
    }
    path <- c(path, nxt)
    if (length(path) > 25) {
      break
    }
  }
}
# If multiple records in the cycle come from the most authoritative source, the sources contradict themselves and
# the direction cannot be determined automatically (in 2026, even the authors in `ref` of these GBIF records were
# unreliable). Then ALL records in the cycle become accepted, so that no wrong current name is given, and they are
# shown for a human decision.
# Human decisions for such cycles: the genera that hold the current names (decision by Matthijs S. Berends, 5 October
# 2026; Capillidium for the former Conidiobolus species, Scolecobasidium for the former Ochroconis and Pseudosigmoidea
# species).
synonym_cycle_current_genera <- c("Capillidium", "Scolecobasidium")
cycle_fixes <- bind_rows(lapply(synonym_cycles, function(cycle_members) {
  # (not named `cycle`, since tibble() below would then use its own new column `cycle`)
  members <- taxonomy[match(cycle_members, taxonomy$fullname), , drop = FALSE]
  best <- members[source_prio(members$source) == min(source_prio(members$source)), , drop = FALSE]
  decided <- members$fullname[members$genus %in% synonym_cycle_current_genera]
  if (length(decided) == 1) {
    tibble(
      cycle = paste(cycle_members, collapse = " | "),
      made_accepted = decided,
      reason = "human decision in `synonym_cycle_current_genera`"
    )
  } else if (nrow(best) == 1) {
    tibble(
      cycle = paste(cycle_members, collapse = " | "),
      made_accepted = unname(synonym_next[best$fullname]),
      reason = paste("statement of", best$fullname, "from", best$source, "(most authoritative source)")
    )
  } else {
    tibble(
      cycle = paste(cycle_members, collapse = " | "),
      made_accepted = cycle_members,
      reason = paste("contradicting statements of", toString(unique(best$source)), "- all accepted, NEEDS A HUMAN DECISION")
    )
  }
}))
review(cycle_fixes, "Cycles of synonyms")
if (nrow(cycle_fixes) > 0) {
  taxonomy <- taxonomy %>%
    mutate(
      fix = fullname %in% cycle_fixes$made_accepted,
      status = if_else(fix, "accepted", status),
      lpsn_renamed_to = if_else(fix, NA_character_, lpsn_renamed_to),
      mycobank_renamed_to = if_else(fix, NA_character_, mycobank_renamed_to),
      gbif_renamed_to = if_else(fix, NA_character_, gbif_renamed_to)
    ) %>%
    select(-fix)
}
rm(synonym_next, synonym_cycles, cycle_fixes, synonym_cycle_current_genera)

# Current names that a source has the wrong way around, as accepted name = its synonym in the source. Only for flaws in
# the source data that no rule can resolve, and explain each.
accepted_name_override <- c(
  # in 2026, COL had Enterocytozoon bieneusi as a synonym of 'Encephalitozoon bieneusi', while Enterocytozoon bieneusi
  # is the accepted name of this human pathogen (decision by Matthijs S. Berends, 5 October 2026)
  "Enterocytozoon bieneusi" = "Encephalitozoon bieneusi",
  # the same for Encephalitozoon cuniculi, which COL had as a synonym of 'Nosema cuniculi' (same decision)
  "Encephalitozoon cuniculi" = "Nosema cuniculi",
  # in 2026, the genus Erysipelatoclostridium (Yutin and Galperin 2013, not in LPSN) and its species were accepted in
  # GBIF and earlier releases, next to the current Thomasclavelia (Lawson et al. 2023) with the same species, so they
  # become synonyms (decision by Matthijs S. Berends, 6 October 2026; E. merdavium is retired, see the next block)
  "Thomasclavelia" = "Erysipelatoclostridium",
  "Thomasclavelia cocleata" = "Erysipelatoclostridium cocleatum",
  "Thomasclavelia ramosa" = "Erysipelatoclostridium ramosum",
  "Thomasclavelia saccharogumia" = "Erysipelatoclostridium saccharogumia",
  "Thomasclavelia spiroformis" = "Erysipelatoclostridium spiroforme"
)
for (accepted_name in names(accepted_name_override)) {
  i_acc <- which(taxonomy$fullname == accepted_name)
  i_syn <- which(taxonomy$fullname == accepted_name_override[accepted_name])
  if (length(i_acc) != 1 || length(i_syn) != 1) {
    warning("accepted_name_override: ", accepted_name, " or its synonym not found, check this!", call. = FALSE)
    next
  }
  taxonomy$status[i_acc] <- "accepted"
  taxonomy[i_acc, c("lpsn_renamed_to", "mycobank_renamed_to", "gbif_renamed_to")] <- NA_character_
  taxonomy$status[i_syn] <- "synonym"
  taxonomy[i_syn, c("lpsn_renamed_to", "mycobank_renamed_to", "gbif_renamed_to")] <- NA_character_
  # the synonym points to the accepted name by the identifier that the accepted name has
  for (id in c("lpsn", "mycobank", "gbif")) {
    if (!is.na(taxonomy[[id]][i_acc])) {
      taxonomy[[paste0(id, "_renamed_to")]][i_syn] <- taxonomy[[id]][i_acc]
      break
    }
  }
  # the other synonyms of the former accepted name now point to the accepted name as well
  for (id in c("lpsn", "mycobank", "gbif")) {
    old_id <- taxonomy[[id]][i_syn]
    if (!is.na(old_id)) {
      repoint <- which(taxonomy[[paste0(id, "_renamed_to")]] %in% old_id)
      taxonomy[[paste0(id, "_renamed_to")]][repoint] <- taxonomy[[paste0(id, "_renamed_to")]][i_syn]
    }
  }
}
rm(accepted_name_override, accepted_name, i_acc, i_syn, id, old_id, repoint)

# Released taxa that are retired by a human decision in data-raw/microorganisms_files/mo_code_retirements.csv (the
# automatic retirements of the taxonomy build are handled above). Released taxa are otherwise never removed, but these
# can enter in other ways than the restoration above, e.g. as previously manually added entries.
retired_by_decision <- read_mo_retirements(".") %>%
  filter(decided_by %unlike% "^taxonomy build")
review(
  taxonomy %>% filter(fullname %in% retired_by_decision$fullname) %>% select(fullname, rank, status, source),
  "Released taxa that are retired by a human decision, these will be removed"
)
taxonomy <- taxonomy %>%
  filter(!fullname %in% retired_by_decision$fullname)
rm(retired_by_decision)


# Synonyms without a current name -----------------------------------------------------------------

# A synonym without a current name promises something the data cannot deliver, as as.mo() cannot update it. Most of
# these are manually added entries and released records restored from earlier releases, that took their status from
# older sources: e.g. 'Mycobacterium orygis' was a GBIF synonym without a current name in v3.0.1, while LPSN has it as
# a preferred name that is not validly published. Rules (decision by Matthijs S. Berends, 7 October 2026):
# - prokaryotes are looked up in LPSN: a synonym there is linked to its correct name (if that is in the data set), a
#   correct name becomes accepted, and a name that is not validly published is only kept if it is protected (see
#   'Add parent identifiers', e.g. all taxa of the WHO priority genera), as accepted;
# - all others are removed, and their released codes are retired with the reason in mo_code_retirements.csv, so that
#   they are never given to another taxon and as.mo() translates them to NA with that reason;
# - a record that still has children that are not removed gets the status 'unknown' instead, so that the hierarchy
#   stays complete.
synonym_target <- function(df) {
  # the record of the current name, with the priority of synonym_mo_to_accepted_mo(): LPSN > MycoBank > GBIF
  target <- rep(NA_integer_, nrow(df))
  for (id in c("gbif", "mycobank", "lpsn")) {
    t <- match(df[[paste0(id, "_renamed_to")]], df[[id]], incomparables = NA)
    target[!is.na(t)] <- t[!is.na(t)]
  }
  target
}
synonyms_without_current_name <- function(df) {
  target <- synonym_target(df)
  resolves <- df$status != "synonym"
  for (step in seq_len(10)) {
    resolves <- resolves | (df$status == "synonym" & resolves[target] %in% TRUE)
  }
  df$status == "synonym" & !resolves
}

no_current_name <- which(synonyms_without_current_name(taxonomy))
message(length(no_current_name), " synonyms without a current name")
lpsn_outcome <- rep(NA_character_, nrow(taxonomy))
prokaryotes <- no_current_name[taxonomy$domain[no_current_name] %in% c("Bacteria", "Archaea") &
  taxonomy$rank[no_current_name] %in% c("phylum", "class", "order", "family", "genus", "species", "subspecies")]
for (i in prokaryotes) {
  lpsn <- get_lpsn_and_author(rank = taxonomy$rank[i], name = taxonomy$fullname[i])
  if (is.na(lpsn["lpsn"])) {
    lpsn_outcome[i] <- "not in LPSN"
    next
  }
  taxonomy <- apply_lpsn_result(taxonomy, i, lpsn)
  taxonomy$ref[i] <- get_author_year(taxonomy$ref[i])
  lpsn_outcome[i] <- paste0("LPSN: ", lpsn["taxonomic_status"], if (!is.na(lpsn["correct_name"])) paste0(" of ", lpsn["correct_name"]))
}
save_lpsn_cache()
# names that are not validly published, according to the policy in 'Add parent identifiers'
not_valid <- intersect(prokaryotes, which(taxonomy$status == "not validly published"))
protected_now <- (paste(taxonomy$genus, taxonomy$species) %in% paste(pathogens$genus, pathogens$species) | is_who_priority_genus(taxonomy$genus) |
  (taxonomy$rank == "genus" & taxonomy$genus %in% pathogens$genus))
taxonomy$status[not_valid[protected_now[not_valid]]] <- "accepted"
taxonomy$status[not_valid[!protected_now[not_valid]]] <- "synonym" # (without a current name, so removed below)
rm(not_valid, protected_now)

# what remains is removed, unless it has children that are kept
to_remove <- synonyms_without_current_name(taxonomy)
has_kept_children <- rep(FALSE, nrow(taxonomy))
for (rank_name in c("phylum", "class", "order", "family", "genus", "species")) {
  parents <- which(to_remove & taxonomy$rank == rank_name)
  if (length(parents) == 0) next
  child_rank <- c("phylum", "class", "order", "family", "genus", "species", "subspecies")
  child_rank <- child_rank[seq(which(child_rank == rank_name) + 1, length(child_rank))]
  kept_children <- taxonomy %>%
    filter(!to_remove, rank %in% child_rank) %>%
    distinct(domain, value = .data[[rank_name]])
  if (rank_name == "species") {
    kept_children <- taxonomy %>%
      filter(!to_remove, rank == "subspecies") %>%
      distinct(domain, value = paste(genus, species))
  }
  parent_value <- if (rank_name == "species") taxonomy$fullname[parents] else taxonomy[[rank_name]][parents]
  has_kept_children[parents] <- paste(taxonomy$domain[parents], parent_value) %in% paste(kept_children$domain, kept_children$value)
}
outcome <- taxonomy %>%
  mutate(
    lpsn_outcome = lpsn_outcome,
    outcome = case_when(
      seq_len(n()) %in% no_current_name & !to_remove ~ "kept (current name or protected)",
      to_remove & has_kept_children ~ "status 'unknown' (has children that are kept)",
      to_remove ~ "removed"
    )
  ) %>%
  filter(!is.na(outcome))
review(
  outcome %>% select(fullname, domain, rank, source, prevalence, lpsn_outcome, outcome),
  "Synonyms without a current name and their outcome"
)
taxonomy$status[to_remove & has_kept_children] <- "unknown"
removed <- taxonomy[to_remove & !has_kept_children, , drop = FALSE]
taxonomy <- taxonomy[!(to_remove & !has_kept_children), , drop = FALSE]

# retire the released codes of the removed records
registry_now <- read_mo_registry(".")
retirements_now <- read_mo_retirements(".")
# (by rank and name, not by domain: a taxon can have moved to another domain since its release, such as the
# microsporidia, released in the Protozoa; only names that are not in the data set anymore)
new_retirements <- registry_now %>%
  filter(
    paste(rank, mo_name_without_suffix(fullname)) %in% paste(removed$rank, removed$fullname),
    !mo_name_without_suffix(fullname) %in% taxonomy$fullname,
    !mo %in% retirements_now$mo
  ) %>%
  mutate(
    lpsn_outcome = outcome$lpsn_outcome[match(paste(rank, fullname), paste(outcome$rank, outcome$fullname))],
    reason = paste0(
      "synonym without a current name, in no current source as a current name",
      if_else(is.na(lpsn_outcome), "", paste0(" (", lpsn_outcome, ")"))
    ),
    decided_by = "taxonomy build (automatic, synonym without a current name)",
    date = as.character(Sys.Date())
  ) %>%
  select(mo, fullname, reason, decided_by, date)
review(new_retirements, "Released MO codes that are retired, since their taxon was a synonym without a current name")
if (nrow(new_retirements) > 0) {
  utils::write.csv(
    bind_rows(retirements_now, new_retirements) %>% arrange(mo),
    retirements_file,
    row.names = FALSE,
    na = ""
  )
}
rm(no_current_name, lpsn_outcome, prokaryotes, to_remove, has_kept_children, outcome, removed, registry_now,
   retirements_now, new_retirements)


# Fix genera that are synonyms while they contain accepted species --------------------------------

# e.g. in 2026, MycoBank listed the genus Blastomyces as a synonym, while Blastomyces dermatitidis was accepted
genera_with_accepted_species <- taxonomy %>%
  filter(rank == "species", status == "accepted") %>%
  distinct(domain, genus)
genus_synonyms_to_fix <- taxonomy %>%
  filter(rank == "genus", status == "synonym") %>%
  semi_join(genera_with_accepted_species, by = c("domain", "genus"))
review(
  genus_synonyms_to_fix %>% select(fullname, domain, family, source, lpsn_renamed_to, mycobank_renamed_to, gbif_renamed_to),
  "Synonym genera with accepted species, these will become accepted"
)
taxonomy <- taxonomy %>%
  mutate(
    fix = fullname %in% genus_synonyms_to_fix$fullname,
    status = if_else(fix, "accepted", status),
    lpsn_renamed_to = if_else(fix, NA_character_, lpsn_renamed_to),
    mycobank_renamed_to = if_else(fix, NA_character_, mycobank_renamed_to),
    gbif_renamed_to = if_else(fix, NA_character_, gbif_renamed_to)
  ) %>%
  select(-fix)
rm(genera_with_accepted_species, genus_synonyms_to_fix)


# Fill empty families of genera -------------------------------------------------------------------

# Some genera have no family in their source, while GBIF/COL has one for the same genus in the same order (e.g. in 2026,
# MycoBank had no family for Lacazia, while COL has Ajellomycetaceae). For prokaryotes, LPSN is followed as it is.
gbif_genus_family <- current_gbif %>%
  filter(taxonRank == "genus", !is.na(family), family != "", !domain %in% c("Bacteria", "Archaea")) %>%
  distinct(domain, genus, .keep_all = TRUE) %>%
  select(domain, genus, gbif_order = order, gbif_family = family)
# Some clinically relevant genera have no family in any source (in 2026, COL had no family for these genera). Only
# add a genus here if its family is certain, and explain it.
genus_family_override <- c(
  # the human roundworm, family Ascarididae Baird, 1853
  "Ascaris" = "Ascarididae",
  # the free-living amoeba causing granulomatous amoebic encephalitis, family Balamuthiidae Cavalier-Smith, 2004
  "Balamuthia" = "Balamuthiidae"
)
family_fill <- taxonomy %>%
  filter(rank %in% c("genus", "species", "subspecies"), genus != "", family == "", !domain %in% c("Bacteria", "Archaea")) %>%
  left_join(gbif_genus_family, by = c("domain", "genus")) %>%
  mutate(new_family = case_when(
    genus %in% names(genus_family_override) ~ unname(genus_family_override[genus]),
    !is.na(gbif_family) & (order == "" | order == gbif_order) ~ gbif_family,
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(new_family))
review(
  family_fill %>% filter(rank == "genus") %>% select(fullname, domain, order, new_family, source),
  "Genera without a family that get a family from GBIF or the override list"
)
taxonomy <- taxonomy %>%
  left_join(family_fill %>% select(domain, fullname, new_family), by = c("domain", "fullname")) %>%
  mutate(family = coalesce(new_family, family)) %>%
  select(-new_family) %>%
  # the new families need records of their own
  add_missing_parents(current_gbif) %>%
  arrange(fullname, source_prio(source)) %>%
  distinct(fullname, .keep_all = TRUE) %>%
  compute_prevalence()
rm(gbif_genus_family, genus_family_override, family_fill)


# Add missing parent records ----------------------------------------------------------------------

# Every phylum, class, order and family that a current record refers to must exist as a record in the same domain
# (in October 2026 e.g. the family Plasmodiidae was only in the Chromista, as in COL, while Plasmodium is in the
# Protozoa). A record of that name in another domain that is not used there is moved, otherwise a record is added.
# Genera that are moved to another domain as a whole, as their only source has them in another domain than all
# related taxa (`genus_domain_override` above only chooses between domains in which a genus occurs). Their parent
# records follow below. Decision by Matthijs S. Berends, 7 October 2026:
genus_domain_move <- c(
  # a plasmodial slime mould (Myxomycetes) that is only in MycoBank, while all other Myxomycetes are in Protozoa
  "Dictydiaethalium" = "Protozoa"
)
for (g in names(genus_domain_move)) {
  i <- which(taxonomy$genus == g & taxonomy$rank %in% c("genus", "species", "subspecies"))
  message("Moving ", length(i), " records of ", g, " to ", genus_domain_move[g])
  taxonomy$domain[i] <- genus_domain_move[g]
  taxonomy$kingdom[i] <- genus_domain_move[g]
}

parent_log <- tibble()
for (r in c("family", "order", "class", "phylum")) {
  lower <- c("phylum", "class", "order", "family", "genus", "species", "subspecies")
  lower <- lower[seq(which(lower == r) + 1, length(lower))]
  above <- c("kingdom", "phylum", "class", "order", "family")
  above <- above[seq_len(which(above == r) - 1)]
  children <- taxonomy %>%
    filter(status != "synonym", rank %in% lower, .data[[r]] != "", .data[[r]] %unlike% "^[(]")
  needed <- children %>% distinct(domain, name = .data[[r]])
  missing_parents <- needed %>%
    # (by the rank field, as a record can have a rank suffix, e.g. 'Acantharia {class}')
    filter(!paste(domain, name) %in% paste(taxonomy$domain[taxonomy$rank == r], taxonomy[[r]][taxonomy$rank == r]))
  for (k in seq_len(nrow(missing_parents))) {
    dom <- missing_parents$domain[k]
    nm <- missing_parents$name[k]
    kid <- children %>%
      filter(domain == dom, .data[[r]] == nm) %>%
      arrange(prevalence, fullname) %>%
      slice(1)
    same_name <- which(taxonomy$fullname == nm)
    used_elsewhere <- any(children$domain != dom & children[[r]] == nm)
    if (length(same_name) == 1 && taxonomy$rank[same_name] == r && !used_elsewhere) {
      old_domain <- taxonomy$domain[same_name]
      taxonomy$domain[same_name] <- dom
      for (a in above) taxonomy[[a]][same_name] <- kid[[a]]
      action <- paste("moved from", old_domain)
    } else if (length(same_name) == 0) {
      new_record <- kid
      new_record$fullname <- nm
      new_record$rank <- r
      new_record$status <- "accepted"
      new_record$source <- "manually added"
      new_record$ref <- NA_character_
      for (col in intersect(c("lpsn", "lpsn_parent", "lpsn_renamed_to", "mycobank", "mycobank_parent",
                              "mycobank_renamed_to", "gbif", "gbif_parent", "gbif_renamed_to"), colnames(new_record))) {
        new_record[[col]] <- NA_character_
      }
      for (lf in c(r, lower)) new_record[[lf]] <- ""
      new_record[[r]] <- nm
      taxonomy <- bind_rows(taxonomy, new_record)
      action <- "added"
    } else {
      action <- "not possible, the name is used in another domain or at another rank"
    }
    parent_log <- bind_rows(parent_log, tibble(rank = r, name = nm, domain = dom, action = action))
  }
}
review(parent_log, "Missing parent records that were added or moved")
rm(parent_log, lower, above, children, needed, missing_parents)


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

# @resume-block mo_registry after=taxonomy2c
# (the MO code registry is used by the final checks, when resuming from a checkpoint)
# parts of existing MO codes, split on the underscore (genus codes can contain digits or hyphens, such as
# B_RBCLM1 or B_CND-P, so a regex with only [A-Z] does not work and caused repeated elements like
# B_SCLLM_CNNM_LNSM_LNSM_LNSM in earlier versions)
mo_part <- function(mo, i) {
  vapply(strsplit(as.character(mo), "_", fixed = TRUE), function(x) if (length(x) >= i) x[i] else NA_character_, character(1))
}
# parts of a name, e.g. the species of 'Escherichia coli'
mo_part_of_name <- function(x, i) {
  vapply(strsplit(x, " ", fixed = TRUE), function(x) if (length(x) >= i) x[i] else "", character(1))
}

# Existing codes come from the MO code registry (data-raw/microorganisms_files/), which contains every code of
# every release since v2.0.0, and NOT from the previous data set: codes that only existed in a development
# version carry no weight (decision by Matthijs S. Berends, 3 October 2026). A registered code is never given
# to another taxon, also not when its taxon is not in the new data set anymore, since users store these codes.
mo_registry <- read_mo_registry(".")
mo_renames <- read_mo_renames(".")
if (is.null(mo_registry)) {
  stop("MO code registry not found, see data-raw/microorganisms_files/README.md", call. = FALSE)
}
release_order <- function(release) {
  as.integer(factor(release, levels = unique(release[order(numeric_version(sub("^v", "", release)))])))
}
existing_mo_tbl <- mo_registry %>%
  # approved renames, e.g. a corrected spelling: the code continues with the new name
  left_join(mo_renames %>% select(mo, new_name), by = "mo") %>%
  mutate(
    fullname = coalesce(new_name, fullname),
    species = if_else(!is.na(new_name) & rank %in% c("species", "subspecies"), mo_part_of_name(new_name, 2), species),
    subspecies = if_else(!is.na(new_name) & rank == "subspecies", mo_part_of_name(new_name, 3), subspecies)
  ) %>%
  select(-new_name) %>%
  # if a taxon had multiple codes over time, the code of the most recent release is used
  arrange(desc(release_order(last_release)), desc(release_order(first_release))) %>%
  distinct(domain, rank, fullname, .keep_all = TRUE)

# all codes that have ever been used must never be given to another taxon, also not when their taxon
# is not in the new data set anymore (users may have stored old codes)
reserved_genus_codes <- mo_registry %>%
  filter(rank %in% c("genus", "species", "subspecies")) %>%
  transmute(domain, code = mo_part(mo, 2)) %>%
  distinct()
reserved_species_codes <- mo_registry %>%
  filter(rank %in% c("species", "subspecies")) %>%
  transmute(domain, genus_code = mo_part(mo, 2), code = mo_part(mo, 3)) %>%
  distinct()
# (until October 2026, subspecies codes were not reserved, so a registered code of a removed subspecies could be
# given to another name, e.g. F_CANDD_MELBS_MMBR of Candida melibiosi membranaefaciens to C. m. membranifaciens)
reserved_subspecies_codes <- mo_registry %>%
  filter(rank == "subspecies") %>%
  transmute(domain, genus_code = mo_part(mo, 2), species_code = mo_part(mo, 3), code = mo_part(mo, 4)) %>%
  filter(!is.na(code)) %>%
  distinct()
# @end-resume-block


# phylum until family are abbreviated with 8 characters and prefixed with their rank

# keep the old code where available, and give every new taxon the first candidate code that is not in use
# and was never used before (`reserved`) - if all candidates are taken, a number is added to the first one
assign_codes <- function(old, candidates, reserved = character(0)) {
  used <- unique(c(reserved, old[!is.na(old)]))
  out <- old
  for (i in which(is.na(old))) {
    cand <- unique(candidates[[i]][!is.na(candidates[[i]]) & candidates[[i]] != ""])
    pick <- cand[!cand %in% used][1]
    k <- 1
    while (is.na(pick)) {
      if (!paste0(cand[1], k) %in% used) {
        pick <- paste0(cand[1], k)
      }
      k <- k + 1
    }
    out[i] <- pick
    used <- c(used, pick)
  }
  out
}

# Phylum until family - keep old (from the registry) and fill up for new ones, never with a registered code
higher_rank_codes <- function(rank_name, tag) {
  reserved <- mo_registry %>%
    filter(rank == rank_name) %>%
    transmute(domain, code = sub("^[A-Z]{1,2}_", "", mo))
  out <- taxonomy %>%
    filter(rank == rank_name) %>%
    distinct(domain, name = .data[[rank_name]]) %>%
    left_join(
      existing_mo_tbl %>%
        filter(rank == rank_name) %>%
        transmute(domain, name = mo_name_without_suffix(fullname), mo_old = sub("^[A-Z]{1,2}_", "", mo)) %>%
        distinct(domain, name, .keep_all = TRUE),
      by = c("domain", "name")
    ) %>%
    group_by(domain) %>%
    mutate(
      code = assign_codes(
        old = mo_old,
        candidates = Map(
          c,
          AMR:::abbreviate_mo(name, minlength = 8, prefix = paste0("[", tag, "]_")),
          AMR:::abbreviate_mo(name, minlength = 9, prefix = paste0("[", tag, "]_"))
        ),
        reserved = reserved$code[reserved$domain == cur_group()$domain]
      )
    ) %>%
    ungroup()
  if (anyDuplicated(paste(out$domain, out$code)) > 0 || anyNA(out$code)) {
    stop("Duplicate MO codes for ", rank_name, "!", call. = FALSE)
  }
  out %>%
    select(domain, name, code) %>%
    setNames(c("domain", rank_name, paste0("mo_", rank_name)))
}
mo_phylum <- higher_rank_codes("phylum", "PHL")
mo_class <- higher_rank_codes("class", "CLS")
mo_order <- higher_rank_codes("order", "ORD")
mo_family <- higher_rank_codes("family", "FAM")


# construct code part for genus - keep old code where available and generate new ones where needed
mo_genus <- taxonomy %>%
  filter(rank == "genus") %>%
  arrange(status, fullname) %>% # sort on accepted < synonym < unknown, added 2026-09-03
  distinct(domain, genus) %>%
  # get available old MO codes
  left_join(
    existing_mo_tbl %>%
      filter(rank == "genus") %>%
      transmute(
        mo_genus_old = mo_part(mo, 2),
        domain,
        genus
      ) %>%
      distinct(domain, genus, .keep_all = TRUE),
    by = c("domain", "genus")
  ) %>%
  # since domain is part of the code, genus abbreviations may be duplicated between domains
  group_by(domain) %>%
  mutate(
    mo_genus = assign_codes(
      old = mo_genus_old,
      candidates = Map(
        c,
        AMR:::abbreviate_mo(genus, 5),
        AMR:::abbreviate_mo(genus, 6),
        AMR:::abbreviate_mo(genus, 7),
        AMR:::abbreviate_mo(genus, 8),
        paste0(AMR:::abbreviate_mo(genus, 5), 1),
        paste0(AMR:::abbreviate_mo(genus, 5), 2)
      ),
      reserved = reserved_genus_codes$code[reserved_genus_codes$domain == cur_group()$domain]
    )
  ) %>%
  ungroup()
if (anyDuplicated(paste(mo_genus$domain, mo_genus$mo_genus)) > 0 || anyNA(mo_genus$mo_genus)) {
  stop("Duplicate MO codes for genus!")
}
# no duplicates *within domains*, so keep the right columns for left joining later
mo_genus <- mo_genus %>%
  select(domain, genus, mo_genus)

# same for species - keep old where available and create new per domain-genus where needed:
mo_species <- taxonomy %>%
  filter(rank == "species") %>%
  arrange(status, fullname) %>% # sort on accepted < synonym < unknown, added 2026-09-03
  distinct(domain, genus, species) %>%
  left_join(mo_genus, by = c("domain", "genus")) %>%
  left_join(
    existing_mo_tbl %>%
      filter(rank == "species") %>%
      transmute(
        mo_species_old = mo_part(mo, 3),
        domain,
        genus,
        species
      ) %>%
      # only valid codes (in 2026, some old codes contained quotes, e.g. F_FUSRM_LO"ES)
      filter(mo_species_old %like_case% "^[A-Z0-9]+$") %>%
      distinct(domain, genus, species, .keep_all = TRUE),
    by = c("domain", "genus", "species")
  ) %>%
  group_by(domain, genus) %>%
  mutate(
    mo_species = assign_codes(
      old = mo_species_old,
      candidates = Map(
        c,
        AMR:::abbreviate_mo(species, 4, hyphen_as_space = TRUE),
        AMR:::abbreviate_mo(species, 5, hyphen_as_space = TRUE),
        AMR:::abbreviate_mo(species, 6, hyphen_as_space = TRUE),
        AMR:::abbreviate_mo(species, 7, hyphen_as_space = TRUE),
        AMR:::abbreviate_mo(species, 8, hyphen_as_space = TRUE),
        paste0(AMR:::abbreviate_mo(species, 5, hyphen_as_space = TRUE), 1)
      ),
      reserved = reserved_species_codes$code[
        reserved_species_codes$domain == cur_group()$domain &
          reserved_species_codes$genus_code %in% mo_genus[1]
      ]
    )
  ) %>%
  ungroup()
if (anyDuplicated(paste(mo_species$domain, mo_species$genus, mo_species$mo_species)) > 0 || anyNA(mo_species$mo_species)) {
  stop("Duplicate MO codes for species!")
}
# no duplicates *within genera*, so keep the right columns for left joining later
mo_species <- mo_species %>%
  select(domain, genus, species, mo_species)

# same for subspecies - keep old where available and create new per domain-genus-species where needed:
mo_subspecies <- taxonomy %>%
  filter(rank == "subspecies") %>%
  arrange(status, fullname) %>% # sort on accepted < synonym < unknown, added 2026-09-03
  distinct(domain, genus, species, subspecies) %>%
  left_join(
    existing_mo_tbl %>%
      filter(rank %in% c("subspecies", "subsp.", "infraspecies")) %>%
      transmute(
        mo_subspecies_old = mo_part(mo, 4),
        domain,
        genus,
        species,
        subspecies
      ) %>%
      # only valid codes (in 2026, some old codes contained quotes, e.g. F_FUSRM_LO"ES)
      filter(mo_subspecies_old %like_case% "^[A-Z0-9]+$") %>%
      distinct(domain, genus, species, subspecies, .keep_all = TRUE),
    by = c("domain", "genus", "species", "subspecies")
  ) %>%
  left_join(mo_genus, by = c("domain", "genus")) %>%
  left_join(mo_species, by = c("domain", "genus", "species")) %>%
  group_by(domain, genus, species) %>%
  mutate(
    mo_subspecies = assign_codes(
      old = mo_subspecies_old,
      candidates = Map(
        c,
        AMR:::abbreviate_mo(subspecies, 4, hyphen_as_space = TRUE),
        AMR:::abbreviate_mo(subspecies, 5, hyphen_as_space = TRUE),
        AMR:::abbreviate_mo(subspecies, 6, hyphen_as_space = TRUE),
        AMR:::abbreviate_mo(subspecies, 7, hyphen_as_space = TRUE),
        AMR:::abbreviate_mo(subspecies, 8, hyphen_as_space = TRUE),
        paste0(AMR:::abbreviate_mo(subspecies, 5, hyphen_as_space = TRUE), 1)
      ),
      reserved = reserved_subspecies_codes$code[
        reserved_subspecies_codes$domain == cur_group()$domain &
          reserved_subspecies_codes$genus_code %in% mo_genus[1] &
          reserved_subspecies_codes$species_code %in% mo_species[1]
      ]
    )
  ) %>%
  ungroup()
if (anyDuplicated(paste(mo_subspecies$domain, mo_subspecies$genus, mo_subspecies$species, mo_subspecies$mo_subspecies)) > 0 ||
  anyNA(mo_subspecies$mo_subspecies)) {
  stop("Duplicate MO codes for subspecies!")
}
# no duplicates *within species*, so keep the right columns for left joining later
mo_subspecies <- mo_subspecies %>%
  select(domain, genus, species, subspecies, mo_subspecies)

# unknowns - manually added
mo_unknown <- existing_mo_tbl %>%
  filter(fullname %like% "unknown") %>%
  transmute(fullname, mo_unknown = as.character(mo))
# unknowns that were added after the last release keep their code too, as long as it is not registered
mo_unknown <- mo_unknown %>%
  bind_rows(
    microorganisms_old %>%
      filter(fullname %like% "unknown", !fullname %in% mo_unknown$fullname, !mo %in% mo_registry$mo) %>%
      transmute(fullname, mo_unknown = as.character(mo))
  )

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
  review("Duplicate full names after adding MO codes")

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
# (by domain and name, so that a code can never end up in another domain)
manual_mos <- as.character(existing_mo_tbl$mo)[match(
  paste(taxonomy$domain, taxonomy$fullname)[taxonomy$source == "manually added"],
  paste(existing_mo_tbl$domain, existing_mo_tbl$fullname)
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
# (since 2026, missing parents are added again after filtering, so this should be 0 rows)
taxonomy %>%
  filter(mo %like% "__") %>%
  review("Records without a valid MO code, these will be removed")
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

# Salmonella serovars are stored as 'subspecies' in this data set (e.g. 'Salmonella Typhi'), while LPSN still lists
# their old species names as synonyms (e.g. 'Salmonella typhi'). These would make as.mo("S. typhi") return
# Salmonella enterica instead of Salmonella Typhi, so they are removed (until 2026, they were removed as duplicate
# MO codes, since they got the same code as the serovar).
salmonella_serovars <- tolower(taxonomy$fullname[taxonomy$genus == "Salmonella" & taxonomy$rank == "subspecies"])
taxonomy %>%
  filter(genus == "Salmonella", rank == "species", status != "accepted", tolower(fullname) %in% salmonella_serovars) %>%
  select(mo, fullname, status, source) %>%
  review("Salmonella species synonyms that are serovars in this data set, these will be removed")
taxonomy <- taxonomy %>%
  filter(!(genus == "Salmonella" & rank == "species" & status != "accepted" & tolower(fullname) %in% salmonella_serovars))
rm(salmonella_serovars)

# Some integrity checks ---------------------------------------------------------------------------

# are mo codes unique?
taxonomy %>%
  filter(mo %in% .[duplicated(mo), "mo", drop = TRUE]) %>%
  mutate(removed_if_distinct_is_applied = duplicated(mo), .before = 1) %>%
  review("Duplicate MO codes")
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
  review("Duplicate full names")
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
  format(nrow(microorganisms_old), big.mark = " "),
  " rows.\n"
)

# these are the new ones:
taxonomy %>%
  filter(!fullname %in% microorganisms_old$fullname) %>%
  review("New taxa")
# these are to be removed:
microorganisms_old %>%
  filter(!fullname %in% taxonomy$fullname) %>%
  review("Removed taxa")
# be extra critical to these:
microorganisms_old %>%
  filter(!fullname %in% taxonomy$fullname, source == "LPSN", status == "accepted") %>%
  count(rank)
# and to these: clinically relevant taxa that will be removed
microorganisms_old %>%
  filter(!fullname %in% taxonomy$fullname, prevalence < 2, rank %in% c("genus", "species")) %>%
  select(mo, fullname, status, source, prevalence) %>%
  review("Removed taxa that were clinically relevant (prevalence < 2)")
# and these: taxa that moved to another domain (this changes their MO code)
taxonomy %>%
  select(fullname, domain, mo) %>%
  inner_join(microorganisms_old %>% select(fullname, domain_old = domain, mo_old = mo), by = "fullname") %>%
  filter(domain != domain_old) %>%
  count(domain_old, domain, genus = gsub(" .*", "", fullname), sort = TRUE) %>%
  review("Genera that moved to another domain")
# and the family of the most relevant genera, which should be checked every year
taxonomy %>%
  filter(rank == "genus", prevalence <= 1.25) %>%
  select(fullname, domain, phylum, class, order, family, source, status) %>%
  left_join(microorganisms_old %>% select(fullname, family_old = family), by = "fullname") %>%
  filter(family != family_old | family == "") %>%
  review("Relevant genera with a new or empty family")

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
  # a species can have multiple rows, the next rows have an empty species (lag() only filled one row and
  # was then removed by the filter below anyway, so all those values were lost)
  mutate(fullname = species) %>%
  tidyr::fill(fullname, .direction = "down") %>%
  filter(
    !is.na(fullname),
    !is.na(oxygen),
    oxygen %unlike% "tolerant",
    fullname %unlike% "unclassified"
  ) %>%
  select(-species)
bacdive <- bacdive %>%
  # BacDive has values per strain: a minority of strains reported as facultative anaerobe may not decide for the whole
  # species (in 2026, 6 of the 52 records of the obligate anaerobe Clostridioides difficile were 'facultative anaerobe')
  group_by(fullname) %>%
  filter(oxygen %unlike% "facultative" | mean(oxygen %like% "facultative") >= 0.5) %>%
  ungroup() %>%
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
  # a species can have multiple rows, the next rows have an empty species (lag() only filled one row and
  # was then removed by the filter below anyway, so all those values were lost)
  mutate(fullname = species) %>%
  tidyr::fill(fullname, .direction = "down") %>%
  filter(
    !is.na(fullname),
    !is.na(shape),
    fullname %unlike% "unclassified"
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
# (the cause was found in October 2026: the old species code was extracted with a regex that failed on genus
# codes with digits, such as B_RBCLM1 - this is now done with mo_part(), so this should be 0 rows)
taxonomy %>%
  filter(mo %like% "_.*_.*_.*_") %>%
  review("MO codes with repeated elements")
# Pattern: capture the prefix, then remove its repetition before the final epithet
taxonomy$mo[taxonomy$mo %like% "_.*_.*_.*_"] <- gsub(
  "^([A-Z]{1,2}_[A-Z0-9-]+(?:_[A-Z0-9-]+)?)_\\1_([A-Z0-9-]+)$",
  "\\1_\\2",
  taxonomy$mo[taxonomy$mo %like% "_.*_.*_.*_"],
  perl = TRUE
)
taxonomy %>%
  filter(mo %like% "_.*_.*_.*_") %>%
  review("MO codes with repeated elements after fixing")


# Remove childless non-bacterial genera -----------------------------------------------------------

# this removes all genera that have no species, except for the domain of Bacteria, and except for released taxa
# (a taxon that was part of a release is never removed, see 'Restore released taxa')
released_names <- mo_name_without_suffix(read_mo_registry(".")$fullname)
# nor the current name of a synonym, as that synonym would then have no current name (in October 2026 e.g. the
# genus Monotosporella, the current name of Monosporella)
is_current_name_of_synonym <- function(df) {
  syn <- df$status == "synonym"
  (!is.na(df$lpsn) & df$lpsn %in% df$lpsn_renamed_to[syn]) |
    (!is.na(df$mycobank) & df$mycobank %in% df$mycobank_renamed_to[syn]) |
    (!is.na(df$gbif) & df$gbif %in% df$gbif_renamed_to[syn])
}
taxonomy <- taxonomy %>%
  filter(
    rank != "genus" |
      domain == "Bacteria" |
      genus %in% taxonomy$genus[taxonomy$rank == "species"] |
      mo_name_without_suffix(fullname) %in% released_names |
      is_current_name_of_synonym(taxonomy)
  )

# then remove all childless upper taxonomy caused by this (again except for released taxa)
for (rank_name in c("family", "order", "class", "phylum", "kingdom")) {
  has_children <- taxonomy %>%
    filter(rank != rank_name, .data[[rank_name]] != "") %>%
    distinct(.data[[rank_name]]) %>%
    pull()

  n_before <- nrow(taxonomy)
  taxonomy <- taxonomy %>%
    filter(!(rank == rank_name & !fullname %in% has_children & !mo_name_without_suffix(fullname) %in% released_names &
      !is_current_name_of_synonym(taxonomy)))
  message("Removed ", n_before - nrow(taxonomy), " childless ", rank_name, " entries")
}
rm(released_names, is_current_name_of_synonym)

# records that were only inferred by this script get the same label as before
taxonomy$source[taxonomy$source == "inferred"] <- "manually added"


# Harmonise the higher taxonomy ---------------------------------------------------------------------

# Every record takes its higher taxonomy from the record of its parent, from the top down, so that e.g. a genus is
# always in the same order, class, phylum and kingdom as its family record. Without this, records from different
# sources disagreed (in October 2026 e.g. hundreds of genera had another phylum or order than their family).
harmonise_higher_taxonomy <- function(df) {
  ranks <- c("kingdom", "phylum", "class", "order", "family", "genus")
  steps <- list(
    c(child = "class", parent = "phylum"),
    c(child = "order", parent = "class"),
    c(child = "family", parent = "order"),
    c(child = "genus", parent = "family"),
    c(child = "species", parent = "genus"),
    c(child = "subspecies", parent = "genus")
  )
  target <- rep(NA_integer_, nrow(df))
  for (id in c("gbif", "mycobank", "lpsn")) {
    t <- match(df[[paste0(id, "_renamed_to")]], df[[id]], incomparables = NA)
    target[!is.na(t) & df$status == "synonym"] <- t[!is.na(t) & df$status == "synonym"]
  }
  for (step in steps) {
    parent_rank <- step[["parent"]]
    above <- ranks[seq_len(which(ranks == parent_rank) - 1)]
    # a parent name that is a synonym is replaced by its current name (in October 2026 e.g. genera in the family
    # 'Acetohalobiaceae' or 'Ruminococcaceae', which are synonyms); not the genus, which is part of the name itself
    if (parent_rank %in% c("kingdom", "phylum", "class", "order", "family")) {
      syn_parent <- which(df$rank == parent_rank & df$status == "synonym" & !is.na(target) & df$rank[target] %in% parent_rank)
      outdated <- which(df[[parent_rank]] %in% df$fullname[syn_parent] & df$rank != parent_rank & df$status != "synonym")
      df[[parent_rank]][outdated] <- df$fullname[target[syn_parent]][match(df[[parent_rank]][outdated], df$fullname[syn_parent])]
    }
    # the current record of the parent, or else its synonym record (e.g. the species of the synonym genus Amphithrix)
    parents <- df[df$rank == parent_rank, , drop = FALSE]
    parents <- parents[order(parents$status == "synonym"), , drop = FALSE]
    parents <- parents[!duplicated(paste(parents$domain, parents[[parent_rank]])), , drop = FALSE]
    child <- which(df$rank == step[["child"]] & df[[parent_rank]] != "")
    p <- match(paste(df$domain[child], df[[parent_rank]][child]), paste(parents$domain, parents[[parent_rank]]))
    child <- child[!is.na(p)]
    p <- p[!is.na(p)]
    for (r in above) {
      df[[r]][child] <- parents[[r]][p]
    }
  }
  df
}
taxonomy_before <- taxonomy
taxonomy <- harmonise_higher_taxonomy(taxonomy)
changed <- which(Reduce(`|`, lapply(c("kingdom", "phylum", "class", "order", "family"), function(r) taxonomy[[r]] != taxonomy_before[[r]])))
review(
  bind_cols(
    taxonomy[changed, c("fullname", "rank", "domain")],
    taxonomy_before[changed, c("kingdom", "phylum", "class", "order", "family")] %>% rename_with(~ paste0(.x, "_before")),
    taxonomy[changed, c("kingdom", "phylum", "class", "order", "family")]
  ),
  "Records whose higher taxonomy was harmonised with the record of their parent"
)
rm(taxonomy_before, changed)


# Some final checks -------------------------------------------------------------------------------

# all our previously manually added names should be in it - review the ones that are not
# (until 2026, also records that were only inferred were labelled 'manually added', so many will be gone)
microorganisms_old %>%
  filter(source == "manually added", !fullname %in% taxonomy$fullname) %>%
  select(mo, fullname, rank, status, prevalence) %>%
  review("Previously manually added taxa that are not in the new data set")

# 'renamed to' identifiers of records that were removed after the earlier clean-up (e.g. by deduplication) lead
# nowhere, so remove them again: these synonyms then have no current name, instead of a broken one
n_dangling <- sum(
  (!is.na(taxonomy$lpsn_renamed_to) & !taxonomy$lpsn_renamed_to %in% taxonomy$lpsn) |
    (!is.na(taxonomy$mycobank_renamed_to) & !taxonomy$mycobank_renamed_to %in% taxonomy$mycobank) |
    (!is.na(taxonomy$gbif_renamed_to) & !taxonomy$gbif_renamed_to %in% taxonomy$gbif)
)
taxonomy <- taxonomy %>%
  mutate(
    lpsn_renamed_to = if_else(lpsn_renamed_to %in% lpsn, lpsn_renamed_to, NA_character_),
    mycobank_renamed_to = if_else(mycobank_renamed_to %in% mycobank, mycobank_renamed_to, NA_character_),
    gbif_renamed_to = if_else(gbif_renamed_to %in% gbif, gbif_renamed_to, NA_character_)
  )
message("Removed dangling 'renamed to' identifiers of ", n_dangling, " records")

# An identifier of a source denotes one taxon, so records of the same domain and rank that share one are the same taxon
# under two spellings (in October 2026 e.g. Candida haemulonii and 'Candida haemulonis', restored from earlier releases
# with the same GBIF identifier). Records of another domain or rank that share an identifier are not the same taxon
# (old GBIF identifiers collide, e.g. 439 for both the fungal order Septobasidiales and the nematode order Strongylida),
# these are only listed. One record is kept as the current name: the one with its name in a current source, else for
# prokaryotes the correct name according to LPSN, else the first in alphabetical order (then marked for review). The
# others become its synonyms, so that their codes keep working.
names_in_current_sources <- function() {
  read_names <- function(file) {
    if (!file.exists(file)) {
      return(character(0))
    }
    d <- readRDS(file)
    if ("fullname" %in% colnames(d)) {
      return(d$fullname)
    }
    trimws(gsub(" +", " ", paste(
      ifelse(is.na(d$genus), "", d$genus), ifelse(is.na(d$species), "", d$species), ifelse(is.na(d$subspecies), "", d$subspecies)
    )))
  }
  unique(c(read_names("data-raw/taxonomy_lpsn.rds"), read_names("data-raw/taxonomy_mycobank.rds"), read_names("data-raw/taxonomy_gbif.rds")))
}
source_names_now <- names_in_current_sources()
shared_id_log <- tibble()
for (id in c("lpsn", "mycobank", "gbif")) {
  x <- taxonomy[[id]]
  dup_ids <- unique(x[!is.na(x) & duplicated(x)])
  for (dup_id in dup_ids) {
    rows <- which(taxonomy[[id]] == dup_id)
    if (length(rows) < 2) next
    if (length(unique(paste(taxonomy$domain[rows], taxonomy$rank[rows]))) > 1) {
      # not the same taxon: the identifier is removed from the records that are not in a current source
      not_current <- rows[!taxonomy$fullname[rows] %in% source_names_now]
      if (length(not_current) == length(rows)) not_current <- rows
      taxonomy[[id]][not_current] <- NA_character_
      # (without its identifier, the record is no record of that source anymore)
      taxonomy$source[not_current[is.na(taxonomy$lpsn[not_current]) & is.na(taxonomy$mycobank[not_current]) & is.na(taxonomy$gbif[not_current])]] <- "manually added"
      shared_id_log <- bind_rows(shared_id_log, tibble(
        id = id, identifier = dup_id, kept = NA_character_,
        made_synonym = NA_character_,
        decided_by = paste("other domain or rank, identifier removed from:", paste(taxonomy$fullname[not_current], collapse = "; "))
      ))
      next
    }
    in_source <- taxonomy$fullname[rows] %in% source_names_now
    lpsn_known <- rep(FALSE, length(rows))
    lpsn_results <- NULL
    if (!any(in_source) && all(taxonomy$domain[rows] %in% c("Bacteria", "Archaea"))) {
      # the spelling that LPSN knows (e.g. Lacrimispora indica, not 'Lacrimispora indicum')
      lpsn_results <- lapply(rows, function(i) get_lpsn_and_author(taxonomy$rank[i], taxonomy$fullname[i]))
      lpsn_known <- vapply(lpsn_results, function(l) !is.na(l["lpsn"]), logical(1))
    }
    ord <- order(!in_source, !lpsn_known, taxonomy$status[rows] == "synonym", taxonomy$fullname[rows])
    rows <- rows[ord]
    keep <- rows[1]
    decided_by <- if (any(in_source)) "name in a current source" else if (any(lpsn_known)) "spelling known to LPSN" else "alphabetical order, NEEDS REVIEW"
    for (other in rows[-1]) {
      taxonomy$status[other] <- "synonym"
      taxonomy[other, c("lpsn", "mycobank", "gbif", "lpsn_renamed_to", "mycobank_renamed_to", "gbif_renamed_to")] <- NA_character_
      taxonomy[[paste0(id, "_renamed_to")]][other] <- dup_id
      taxonomy$source[other] <- "manually added"
    }
    # the kept spelling takes its LPSN record, so that it is linked to its correct name if it is a synonym there
    if (!is.null(lpsn_results) && lpsn_known[ord][1]) {
      status_before <- taxonomy$status[keep]
      taxonomy <- apply_lpsn_result(taxonomy, keep, lpsn_results[ord][[1]])
      # (a synonym in LPSN whose correct name is not in the data set keeps its status, as it would have no current name)
      if (taxonomy$status[keep] == "synonym" && is.na(taxonomy$lpsn_renamed_to[keep]) &&
        is.na(taxonomy$mycobank_renamed_to[keep]) && is.na(taxonomy$gbif_renamed_to[keep])) {
        taxonomy$status[keep] <- status_before
      }
      taxonomy[[id]][keep] <- dup_id
      taxonomy$ref[keep] <- get_author_year(taxonomy$ref[keep])
      if (taxonomy$status[keep] == "not validly published") taxonomy$status[keep] <- "accepted"
    }
    shared_id_log <- bind_rows(shared_id_log, tibble(
      id = id, identifier = dup_id, kept = taxonomy$fullname[keep],
      made_synonym = paste(taxonomy$fullname[rows[-1]], collapse = "; "), decided_by = decided_by
    ))
  }
}
save_lpsn_cache()
review(shared_id_log, "Records that shared a source identifier, one kept as the current name")
rm(source_names_now, shared_id_log, x, dup_ids, names_in_current_sources)

# only synonyms point to a current name: the status follows the priority of the sources (LPSN > MycoBank > GBIF), so
# e.g. a name that is current in MycoBank but a synonym in COL is accepted, and the COL pointer is left out (until
# October 2026, 1,734 accepted MycoBank records kept such a pointer)
taxonomy <- taxonomy %>%
  mutate(across(c(lpsn_renamed_to, mycobank_renamed_to, gbif_renamed_to), ~ if_else(status == "synonym", .x, NA_character_)))
# and a parent is always another record of a higher rank (e.g. the GBIF record of Animalia had itself as parent)
rank_level <- c("domain" = 0, "kingdom" = 1, "phylum" = 2, "class" = 3, "order" = 4, "family" = 5, "genus" = 6, "species" = 7, "subspecies" = 8)
for (id in c("lpsn", "mycobank", "gbif")) {
  p <- match(taxonomy[[paste0(id, "_parent")]], taxonomy[[id]], incomparables = NA)
  invalid <- !is.na(p) & (p == seq_len(nrow(taxonomy)) |
    (taxonomy$rank %in% names(rank_level) & taxonomy$rank[p] %in% names(rank_level) &
      rank_level[taxonomy$rank[p]] >= rank_level[taxonomy$rank]))
  message("Removed ", sum(invalid), " invalid ", id, " parent identifiers")
  taxonomy[[paste0(id, "_parent")]][invalid] <- NA_character_
}
rm(rank_level, p, invalid)

# synonyms without a current name are no error, but every one of them is a name that as.mo() cannot update
review(mo_synonyms_without_current_name(taxonomy), "Synonyms without a current name")

# the integrity rules of the data set and the MO code registry (shared with the unit tests in
# tests/testthat/helper-microorganisms.R): these must never fail, as the package and its users rely on them
integrity_issues <- mo_integrity_issues(taxonomy, registry = mo_registry, renames = mo_renames, retirements = read_mo_retirements("."))
if (!is.null(mo_integrity_report(integrity_issues))) {
  stop("The new data set breaks integrity rules, fix the cause in this script:\n",
    mo_integrity_report(integrity_issues),
    call. = FALSE
  )
}


# Update other data sets --------------------------------------------------------------------------

# `current = TRUE` replaces synonyms by their current name (decision by Matthijs S. Berends, 7 October 2026: the data
# sets refer only to current names, except intrinsic_resistant, which lists synonyms on purpose)
fix_old_mos <- function(dataset, new_ref, drop = FALSE, col = "mo", current = TRUE) {
  before <- dataset
  
  # the data sets contain codes of the development version
  mo_names <- microorganisms_dev$fullname[match(before[[col]], microorganisms_dev$mo)]
  # names of the development version that do not exist anymore, with their name in the new data set
  renamed_names <- c(
    # EUCAST and MycoBank use the species name Trichophyton indotineae, not a variety of T. mentagrophytes
    "Trichophyton mentagrophytes indotineae" = "Trichophyton indotineae",
    # in LPSN, the genus Mycobacteroides is a synonym of Mycobacterium again
    "Mycobacteroides stephanolepidis" = "Mycobacterium stephanolepidis"
  )
  mo_names <- if_else(mo_names %in% names(renamed_names), unname(renamed_names[mo_names]), mo_names)
  matches <- new_ref$mo[match(mo_names, new_ref$fullname)]
  if (isTRUE(current)) {
    matches <- AMR:::synonym_mo_to_accepted_mo(matches, fill_in_accepted = TRUE, dataset = new_ref)
  }
  unmatched <- !is.na(before[[col]]) & is.na(matches)
  review(
    tibble(mo = as.character(before[[col]][unmatched]), fullname = mo_names[unmatched]) %>% count(mo, fullname),
    paste("Codes without a match in", deparse(substitute(dataset)), col)
  )
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

# (intrinsic_resistant is not updated here: it is computed from the new data set afterwards, see
# data-raw/_reproduction_scripts/reproduction_of_intrinsic_resistant.R)

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

# check again after the conversion to ASCII
integrity_issues <- mo_integrity_issues(taxonomy, registry = mo_registry, renames = mo_renames, retirements = read_mo_retirements("."))
if (!is.null(mo_integrity_report(integrity_issues))) {
  stop("The data set breaks integrity rules after the conversion to ASCII:\n",
    mo_integrity_report(integrity_issues),
    call. = FALSE
  )
}

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
anyNA(clinical_breakpoints$mo)

# load new data sets again
devtools::load_all(".")
source("data-raw/_pre_commit_checks.R")
devtools::load_all(".")


# recreate the data sets that are computed from this one (species groups, EUCAST breakpoints, intrinsic resistance),
# in the right order, see data-raw/_reproduction_scripts/run_after_microorganisms_build.R
source("data-raw/_reproduction_scripts/run_after_microorganisms_build.R")
devtools::load_all(".")

# run the unit tests
Sys.setenv(NOT_CRAN = "true")
testthat::test_file("tests/testthat/test-data.R")
testthat::test_file("tests/testthat/test-mo.R")
testthat::test_file("tests/testthat/test-mo_property.R")
