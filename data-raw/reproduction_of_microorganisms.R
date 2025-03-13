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
# how to conduct AMR data analysis: https://msberends.github.io/AMR/   #
# ==================================================================== #

# ! THIS SCRIPT REQUIRES AT LEAST 16 GB RAM !
# (at least 12 GB will be used by the R session for the size of the files)

# GBIF:
# 1. Go to https://doi.org/10.15468/39omei and find the download link for the
#    latest GBIF backbone taxonony under "Endpoints" and unpack "Taxon.tsv" from it (~1.0 GB)
#    ALSO BE SURE to get the date of release and update R/aa_globals.R later!
# LPSN:
# 2. Go to https://lpsn.dsmz.de/downloads (register first) and download the latest
#    CSV file (~12,5 MB) and rename to "taxonomy.csv"
#    ALSO BE SURE to get the date of release and update R/aa_globals.R later!
# MycoBank:
# 3. Go to https://www.mycobank.org/ and find the download link of all entries
#    (last time https://www.mycobank.org/images/MBList.zip) and unpack
#    "MBList.xlsx" from it (~120 MB)
#    ALSO BE SURE to get the date of release and update R/aa_globals.R later!
# Bartlett:
# 4. For data about human pathogens, we use Bartlett et al. (2022),
#    https://doi.org/10.1099/mic.0.001269. Their latest supplementary material
#    can be found here: https://github.com/padpadpadpad/bartlett_et_al_2022_human_pathogens.
#    Download their latest xlsx file in the `data` folder and save it to our
#    `data-raw` folder.
# 5. Set these locations to the paths where the files are:
folder_location <- "~/Downloads/"
file_gbif <- paste0(folder_location, "/backbone/Taxon.tsv")
file_lpsn <- paste0(folder_location, "taxonomy.csv")
file_mycobank <- paste0(folder_location, "MBList.xlsx")

file_bartlett <- "data-raw/bartlett_et_al_2022_human_pathogens.xlsx"

# 4. Run the rest of this script line by line and check everything :)

if (!file.exists(file_gbif)) stop("GBIF file not found")
if (!file.exists(file_lpsn)) stop("LPSN file not found")
if (!file.exists(file_mycobank)) stop("MB file not found")
if (!file.exists(file_bartlett)) stop("Bartlett et al. Excel file not found")

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
  
  authors2 <- iconv(ref, from = "UTF-8", to = "ASCII//TRANSLIT")
  authors2 <- gsub(" ?\\(Approved Lists [0-9]+\\) ?", " () ", authors2)
  authors2 <- gsub(" [)(]+ $", "", authors2)
  # remove leading and trailing brackets
  authors2 <- trimws(gsub("^[(](.*)[)]$", "\\1", authors2))
  # only take part after brackets if there's a name
  authors2 <- if_else(grepl(".*[)] [a-zA-Z]+.*", authors2),
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
  lastyear <- if_else(lastyear > as.integer(format(Sys.Date(), "%Y")),
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
  authors <- gsub("^([A-Z-][.])+( & ?)?", "", authors, ignore.case = FALSE)
  authors <- gsub("^([A-Z-]+ )+", "", authors, ignore.case = FALSE)
  # remove dots
  authors <- gsub(".", "", authors, fixed = TRUE)
  authors <- gsub("et al", "et al.", authors, fixed = TRUE)
  authors[nchar(authors) <= 3] <- ""
  # combine author and year if year is available
  ref <- if_else(!is.na(lastyear),
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


# to retrieve LPSN and authors from LPSN website
# e.g., get_lpsn_and_author("genus", "Klebsiella")
get_lpsn_and_author <- function(rank, name) {
  name <- gsub("^Candidatus ", "", name)
  url <- paste0("https://lpsn.dsmz.de/", tolower(rank), "/", tolower(name))
  page_txt <- tryCatch(read_html(url), error = function(e) NULL)
  if (is.null(page_txt)) {
    warning("No LPSN found for ", tolower(rank), " '", name, "'")
    lpsn <- NA_character_
    ref <- NA_character_
    status <- "unknown"
  } else {
    page_txt <- page_txt %>%
      html_element("#detail-page") %>%
      html_text()
    lpsn <- gsub(".*Record number:[\r\n\t ]*([0-9]+).*", "\\1", page_txt, perl = FALSE)
    ref <- page_txt %>%
      gsub(".*?Name: (.*[0-9]{4}?).*", "\\1", ., perl = FALSE) %>%
      gsub(name, "", ., fixed = TRUE) %>%
      gsub("^\"?Candidatus ?\"?", "", .) %>%
      trimws()
    status <- trimws(gsub(".*Nomenclatural status:[\r\n\t ]*([a-zA-Z, ]+)[\r\n\t].*", "\\1", page_txt, perl = FALSE))
    if ((status %like% "validly published" & status %unlike% "not valid") | status %like% "[\r\n\t]") {
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
      out <- names(sort(table(current[which(!current %in% c("", NA))]), decreasing = TRUE)[1])
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
  x <- tryCatch(read_html(url),
                error = function(e) {
                  message("Waiting 10 seconds because of error: ", e$message)
                  Sys.sleep(10)
                  read_html(url)
                })
  x <- x %>%
    # class "main-list" is the main table
    html_element(".main-list") %>%
    # get every list element with a set <id> attribute
    html_elements("li[id]")
  pb <- progress_bar$new(total = length(x), format = "[:bar] :current/:total :eta")
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
    ranks[ranks == "domain"] <- "kingdom"
    
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
  message("  => ", length(x), " entries incl. candidates (cleaned total: ", nrow(taxonomy_lpsn_missing), ")")
}
taxonomy_lpsn_missing <- taxonomy_lpsn_missing %>% distinct()
saveRDS(taxonomy_lpsn_missing, "data-raw/taxonomy_lpsn_missing.rds")
# taxonomy_lpsn_missing <- readRDS("data-raw/taxonomy_lpsn_missing.rds")

# had to pick the right genus/family combination here:
taxonomy_lpsn_missing <- taxonomy_lpsn_missing %>% filter(!(genus == "Pusillimonas" & family == "Oscillospiraceae"))
taxonomy_lpsn.bak2 <- taxonomy_lpsn.bak

taxonomy_lpsn <- taxonomy_lpsn %>%
  left_join(taxonomy_lpsn_missing, by = "genus") %>%
  select(kingdom:family, everything()) %>%
  # remove entries like "[Bacteria, no family]" and "[Bacteria, no class]"
  mutate_all(function(x) if_else(x %like_case% " no ", NA_character_, x))

taxonomy_lpsn.bak2 <- taxonomy_lpsn
# download family directly from LPSN website using scraping, by using get_lpsn_and_author()
# try it first:
# get_lpsn_and_author("genus", "Escherichia")
# get_lpsn_and_author("family", "Enterobacteriaceae")
pb <- progress_bar$new(total = length(unique(taxonomy_lpsn$family)), format = "[:bar] :current/:total :eta")
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
      status = unname(tax_info["status"]),
      source = "LPSN",
      lpsn = unname(tax_info["lpsn"]),
      ref = unname(tax_info["ref"])
    ))
}
# download order directly from LPSN website using scraping
# try it first:
# get_lpsn_and_author("order", "Enterobacterales")
pb <- progress_bar$new(total = length(unique(taxonomy_lpsn$order)), format = "[:bar] :current/:total :eta")
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
      status = unname(tax_info["status"]),
      source = "LPSN",
      lpsn = unname(tax_info["lpsn"]),
      ref = unname(tax_info["ref"])
    ))
}
# download class directly from LPSN website using scraping
# try it first:
# get_lpsn_and_author("class", "Gammaproteobacteria")
pb <- progress_bar$new(total = length(unique(taxonomy_lpsn$class)), format = "[:bar] :current/:total :eta")
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
      status = unname(tax_info["status"]),
      source = "LPSN",
      lpsn = unname(tax_info["lpsn"]),
      ref = unname(tax_info["ref"])
    ))
}
# download phylum directly from LPSN website using scraping
# try it first:
# get_lpsn_and_author("phylum", "Pseudomonadota")
pb <- progress_bar$new(total = length(unique(taxonomy_lpsn$phylum)), format = "[:bar] :current/:total :eta")
for (p in unique(taxonomy_lpsn$phylum)) {
  pb$tick()
  if (is.na(p)) next
  tax_info <- get_lpsn_and_author("Phylum", p)
  taxonomy_lpsn <- taxonomy_lpsn %>%
    bind_rows(tibble(
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
# get_lpsn_and_author("kingdom", "Bacteria")
pb <- progress_bar$new(total = length(unique(taxonomy_lpsn$kingdom)), format = "[:bar] :current/:total :eta")
for (k in unique(taxonomy_lpsn$kingdom)) {
  pb$tick()
  if (is.na(k)) next
  tax_info <- get_lpsn_and_author("Domain", k)
  taxonomy_lpsn <- taxonomy_lpsn %>%
    bind_rows(tibble(
      kingdom = k,
      rank = "kingdom",
      status = unname(tax_info["status"]),
      source = "LPSN",
      lpsn = unname(tax_info["lpsn"]),
      ref = unname(tax_info["ref"])
    ))
}

taxonomy_lpsn <- taxonomy_lpsn %>%
  filter(status != "not validly published")

# integrity tests
sort(table(taxonomy_lpsn$rank))
# should only be 'accepted' and 'synonym':
sort(table(taxonomy_lpsn$status))


# Read MycoBank data ------------------------------------------------------------------------------

taxonomy_mycobank <- read_excel(file_mycobank, guess_max = 1e5)
taxonomy_mycobank.bak <- taxonomy_mycobank

taxonomy_mycobank <- taxonomy_mycobank %>%
  transmute(mycobank = `MycoBank #`,
            fullname = gsub(" +", " ", `Taxon name`),
            current = `Current name.Taxon name`,
            ref = paste0(Authors, ", ", `Year of effective publication`),
            rank = `Rank.Rank name`,
            status = `Name status`,
            mycobank_renamed_to = taxonomy_mycobank.bak$`MycoBank #`[match(`Current name.Taxon name`, taxonomy_mycobank.bak$`Taxon name`)],
            Classification
  ) %>%
  separate(Classification, sep = ", ", into = paste0("tax_", letters[1:20]), remove = TRUE) %>% 
  mutate(rank = case_when(
    fullname %like_case% "^[A-Z][a-z]+ .+ .+" ~ "subsp.",
    rank == "-" & fullname %like_case% "^[A-Z][a-z]+ [a-z-]+$" ~ "sp.",
    rank == "-" & paste0(fullname, " ") %in% gsub("(^[A-Z][a-z]+ ).*", "\\1", trimws(taxonomy_mycobank.bak$`Taxon name`), perl = TRUE) ~ "gen.",
    # we take a leap here for the family and order
    rank == "-" & fullname %like% "ceae$" ~ "fam.",
    rank == "-" & fullname %like% "ales$" ~ "ordo",
    # tax_d and tax_e are the subdivision/class columns, kind of; prefer e over d
    rank == "-" & fullname %in% tax_e ~ "cl.",
    rank == "-" & fullname %in% tax_d ~ "cl.",
    TRUE ~ rank
  )) %>% 
  filter(rank %unlike% "sub" & !tolower(rank) %in% c("var.", "sect.", "ser.", "tr.", "f.", "race", "stirps", "*", "-"),
         # we also remove orthographic variants here, because our algorithms will give valid results based on misspelling
         !tolower(status) %in% c("invalid", "deleted", "uncertain", "illegitimate", "orthographic variant", "unavailable")) %>% 
  mutate(mycobank_renamed_to = if_else(mycobank_renamed_to == mycobank, NA_character_, mycobank_renamed_to)) %>% 
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

taxonomy_mycobank <- taxonomy_mycobank %>% 
  mutate(rank = case_match(rank,
                           "sp." ~ "species",
                           "gen." ~ "genus",
                           "fam." ~ "family",
                           "ordo" ~ "order",
                           "cl." ~ "class",
                           "div." ~ "phylum",
                           "regn." ~ "kingdom",
                           .default = paste0("#", rank)))
taxonomy_mycobank %>% filter(rank %like% "#")
taxonomy_mycobank %>% count(rank, sort = TRUE)

taxonomy_mycobank3 <- taxonomy_mycobank

# MycoBank just pasted their taxonomy together into 1 field, and now some classes are in the division (tax_b) column, it's horrible
# so we decide based on the fullname and rank column per record
# use this to determine how far to go:
any(taxonomy_mycobank$fullname[taxonomy_mycobank$rank == "class"] %in% taxonomy_mycobank$tax_e) # replace tax_e with tax_f, etc
any(taxonomy_mycobank$fullname[taxonomy_mycobank$rank == "class"] %in% taxonomy_mycobank$tax_f)
taxonomy_mycobank <- taxonomy_mycobank %>% 
  mutate(kingdom = "Fungi", # we already filtered everything else, and MycoBank 90157 has '-' as current name...
         phylum = case_when(rank == "phylum" ~ fullname,
                            tax_b %in% taxonomy_mycobank$fullname[taxonomy_mycobank$rank == "phylum"] ~ tax_b,
                            tax_c %in% taxonomy_mycobank$fullname[taxonomy_mycobank$rank == "phylum"] ~ tax_c,
                            TRUE ~ ""),
         class = case_when(rank == "class" ~ fullname,
                           tax_c %in% taxonomy_mycobank$fullname[taxonomy_mycobank$rank == "class"] ~ tax_c,
                           tax_d %in% taxonomy_mycobank$fullname[taxonomy_mycobank$rank == "class"] ~ tax_d,
                           tax_e %in% taxonomy_mycobank$fullname[taxonomy_mycobank$rank == "class"] ~ tax_e,
                           TRUE ~ ""),
         order = case_when(rank == "order" ~ fullname,
                           tax_c %in% taxonomy_mycobank$fullname[taxonomy_mycobank$rank == "order"] ~ tax_c,
                           tax_d %in% taxonomy_mycobank$fullname[taxonomy_mycobank$rank == "order"] ~ tax_d,
                           tax_e %in% taxonomy_mycobank$fullname[taxonomy_mycobank$rank == "order"] ~ tax_e,
                           tax_f %in% taxonomy_mycobank$fullname[taxonomy_mycobank$rank == "order"] ~ tax_f,
                           tax_g %in% taxonomy_mycobank$fullname[taxonomy_mycobank$rank == "order"] ~ tax_g,
                           TRUE ~ ""),
         family = case_when(rank == "family" ~ fullname,
                            tax_c %in% taxonomy_mycobank$fullname[taxonomy_mycobank$rank == "family"] ~ tax_c,
                            tax_d %in% taxonomy_mycobank$fullname[taxonomy_mycobank$rank == "family"] ~ tax_d,
                            tax_e %in% taxonomy_mycobank$fullname[taxonomy_mycobank$rank == "family"] ~ tax_e,
                            tax_f %in% taxonomy_mycobank$fullname[taxonomy_mycobank$rank == "family"] ~ tax_f,
                            tax_g %in% taxonomy_mycobank$fullname[taxonomy_mycobank$rank == "family"] ~ tax_g,
                            tax_h %in% taxonomy_mycobank$fullname[taxonomy_mycobank$rank == "family"] ~ tax_h,
                            tax_i %in% taxonomy_mycobank$fullname[taxonomy_mycobank$rank == "family"] ~ tax_i,
                            TRUE ~ ""),
         genus = case_when(rank == "genus" ~ fullname,
                           tax_b %in% taxonomy_mycobank$fullname[taxonomy_mycobank$rank == "genus"] ~ tax_b,
                           tax_c %in% taxonomy_mycobank$fullname[taxonomy_mycobank$rank == "genus"] ~ tax_c,
                           tax_d %in% taxonomy_mycobank$fullname[taxonomy_mycobank$rank == "genus"] ~ tax_d,
                           tax_e %in% taxonomy_mycobank$fullname[taxonomy_mycobank$rank == "genus"] ~ tax_e,
                           tax_f %in% taxonomy_mycobank$fullname[taxonomy_mycobank$rank == "genus"] ~ tax_f,
                           tax_g %in% taxonomy_mycobank$fullname[taxonomy_mycobank$rank == "genus"] ~ tax_g,
                           tax_h %in% taxonomy_mycobank$fullname[taxonomy_mycobank$rank == "genus"] ~ tax_h,
                           tax_i %in% taxonomy_mycobank$fullname[taxonomy_mycobank$rank == "genus"] ~ tax_i,
                           tax_k %in% taxonomy_mycobank$fullname[taxonomy_mycobank$rank == "genus"] ~ tax_j,
                           tax_k %in% taxonomy_mycobank$fullname[taxonomy_mycobank$rank == "genus"] ~ tax_k,
                           TRUE ~ ""),
         species = case_when(rank == "species" & fullname %like% " " ~ gsub(".* (.*)", "\\1", fullname, perl = TRUE),
                             TRUE ~ "")
  )

# FOR 2025: use this to get all the genera with updated names from MO_RELEVANT_GENERA:
# AMR::microorganisms %>% filter(genus %in% MO_RELEVANT_GENERA) %>% pull(fullname) %>% mo_current() %>% mo_genus() %>% unique() %>% sort()

# keep only the relevant ones
taxonomy_mycobank <- taxonomy_mycobank %>% 
  filter(genus %in% AMR:::MO_RELEVANT_GENERA |
           !rank %in% c("genus", "species")) %>% 
  filter(!(genus == "" & rank %in% c("genus", "species")))

# clean up authors and add last columns
taxonomy_mycobank <- taxonomy_mycobank %>% 
  mutate(subspecies = "",
         source = "MycoBank",
         mycobank_parent = NA_character_)

# select final set
taxonomy_mycobank <- taxonomy_mycobank %>%
  select(fullname, kingdom, phylum, class, order, family, genus, species, subspecies,
         rank, status, ref, mycobank, mycobank_parent, mycobank_renamed_to, source)

# not all 'renamed to' records are available, some were even just orthographic variants (that have an invalid status)
taxonomy_mycobank$status[!is.na(taxonomy_mycobank$mycobank_renamed_to) & !taxonomy_mycobank$mycobank_renamed_to %in% taxonomy_mycobank$mycobank] <- "accepted"
taxonomy_mycobank$mycobank_renamed_to[taxonomy_mycobank$status == "accepted"] <- NA

taxonomy_mycobank %>% count(status)


# Read GBIF data ----------------------------------------------------------------------------------

taxonomy_gbif.bak <- vroom(file_gbif, guess_max = 1e5)

# include all fungal orders from the mycobank db
include_fungal_orders <- unique(taxonomy_mycobank$order[!taxonomy_mycobank$order %in% c("", NA)])

# check some columns to validate below filters
taxonomy_gbif.bak %>% count(taxonomicStatus, sort = TRUE)
taxonomy_gbif.bak %>% count(taxonRank, sort = TRUE)

taxonomy_gbif <- taxonomy_gbif.bak %>%
  # immediately filter rows we really never want
  filter(
    # never doubtful status, only accepted and all synonyms, and only ranked items
    taxonomicStatus != "doubtful",
    taxonRank != "unranked",
    #  include these kingdoms (no Chromista)
    kingdom %in% c("Archaea", "Bacteria", "Protozoa") |
      # include all of these fungal orders
      order %in% include_fungal_orders |
      # and all of these important genera (see "data-raw/_pre_commit_checks.R")
      # (they also contain bacteria and protozoa, but these will get higher prevalence scores later on)
      genus %in% AMR:::MO_RELEVANT_GENERA
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
    status = if_else(status == "accepted", "accepted", "synonym"),
    # checked taxonRank - the "form" and "variety" always have a subspecies, so:
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

taxonomy_gbif


# Save intermediate results -----------------------------------------------------------------------

saveRDS(taxonomy_lpsn, "data-raw/taxonomy_lpsn.rds", version = 2)
saveRDS(taxonomy_mycobank, "data-raw/taxonomy_mycobank.rds", version = 2)
saveRDS(taxonomy_gbif, "data-raw/taxonomy_gbif.rds", version = 2)
# this allows to always get back to this point by simply loading the files from data-raw/.


# Add full names ----------------------------------------------------------------------------------

taxonomy_gbif <- taxonomy_gbif %>%
  # clean NAs and add fullname
  mutate(across(kingdom:subspecies, function(x) if_else(is.na(x), "", x)),
         fullname = trimws(case_when(
           rank == "family" ~ family,
           rank == "order" ~ order,
           rank == "class" ~ class,
           rank == "phylum" ~ phylum,
           rank == "kingdom" ~ kingdom,
           TRUE ~ paste(genus, species, subspecies) # already trimmed 6 lines up
         )), .before = 1
  ) %>%
  # keep only one GBIF taxon ID per full name
  arrange(fullname, gbif) %>%
  distinct(kingdom, rank, fullname, .keep_all = TRUE)

taxonomy_lpsn <- taxonomy_lpsn %>%
  # clean NAs and add fullname
  mutate(across(kingdom:subspecies, function(x) if_else(is.na(x), "", x)),
         fullname = trimws(case_when(
           rank == "family" ~ family,
           rank == "order" ~ order,
           rank == "class" ~ class,
           rank == "phylum" ~ phylum,
           rank == "kingdom" ~ kingdom,
           TRUE ~ paste(genus, species, subspecies) # already trimmed 6 lines up
         )), .before = 1
  ) %>%
  # keep only one LPSN record ID per full name
  arrange(fullname, lpsn) %>%
  distinct(kingdom, rank, fullname, .keep_all = TRUE)

taxonomy_mycobank <- taxonomy_mycobank %>%
  # clean NAs and add fullname
  mutate(across(kingdom:subspecies, function(x) if_else(is.na(x), "", x)),
         fullname = trimws(case_when(
           rank == "family" ~ family,
           rank == "order" ~ order,
           rank == "class" ~ class,
           rank == "phylum" ~ phylum,
           rank == "kingdom" ~ kingdom,
           TRUE ~ paste(genus, species, subspecies) # already trimmed 6 lines up
         )), .before = 1
  ) %>%
  # keep only one MycoBank record ID per full name
  arrange(fullname, mycobank) %>%
  distinct(kingdom, rank, fullname, .keep_all = TRUE)


# Combine the datasets ----------------------------------------------------------------------------

taxonomy <- taxonomy_lpsn %>%
  # add fungi
  bind_rows(taxonomy_mycobank) %>% 
  # add GBIF to the bottom
  bind_rows(taxonomy_gbif) %>% 
  # group on unique species
  group_by(kingdom, fullname) %>%
  # fill the NAs in LPSN/GBIF fields and ref with the other source (so LPSN: 123 and GBIF: NA will become LPSN: 123 and GBIF: 123)
  mutate(across(matches("^(lpsn|mycobank|gbif|ref)"), function(x) rep(x[!is.na(x)][1], length(x)))) %>%
  # ungroup again
  ungroup() %>% 
  # only keep unique species per kingdom
  distinct(kingdom, fullname, .keep_all = TRUE) %>% 
  arrange(fullname)

# get missing entries from existing microorganisms data set
taxonomy.old <- AMR::microorganisms %>%
  select(any_of(colnames(taxonomy))) %>%
  filter(
    !paste(kingdom, fullname) %in% paste(taxonomy$kingdom, taxonomy$fullname),
    # these will be added later:
    source != "manually added")
taxonomy <- taxonomy %>%
  bind_rows(taxonomy.old) %>%
  arrange(fullname) %>%
  filter(fullname != "")

# fix rank
taxonomy %>% count(rank, sort = TRUE)
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
taxonomy %>% count(rank, sort = TRUE)

taxonomy0 <- taxonomy

# at this point, it happens that some genera within kingdoms have multiple families / orders, etc., see here:
taxonomy %>% filter(genus != "") %>% group_by(kingdom, genus) %>% filter(n_distinct(family) > 1) %>% View()
# so make this universal
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
# and remove the taxonomy where it must remain empty
taxonomy <- taxonomy %>% 
  mutate(phylum = if_else(rank %in% c("kingdom"), "", phylum),
         class = if_else(rank %in% c("kingdom", "phylum"), "", class),
         order = if_else(rank %in% c("kingdom", "phylum", "class"), "", order),
         family = if_else(rank %in% c("kingdom", "phylum", "class", "order"), "", family),
         genus = if_else(rank %in% c("kingdom", "phylum", "class", "order", "family"), "", genus),
         species = if_else(rank %in% c("kingdom", "phylum", "class", "order", "family", "genus"), "", species),
         subspecies = if_else(rank %in% c("kingdom", "phylum", "class", "order", "family", "genus", "species"), "", subspecies))


# Save intermediate results (0) -------------------------------------------------------------------

saveRDS(taxonomy, "data-raw/taxonomy0.rds")
# taxonomy <- readRDS("data-raw/taxonomy0.rds")

# Add missing and fix old taxonomic entries -------------------------------------------------------

# this part will make sure that the whole taxonomy of every included species exists, so no missing genera, classes, etc.

current_gbif <- taxonomy_gbif.bak %>%
  filter(is.na(acceptedNameUsageID)) %>%
  mutate(
    taxonID = as.character(taxonID),
    parentNameUsageID = as.character(parentNameUsageID)
  )

# add missing kingdoms
taxonomy_all_missing <- taxonomy %>%
  filter(kingdom != "") %>%
  distinct(kingdom) %>%
  mutate(
    fullname = kingdom,
    rank = "kingdom"
  ) %>%
  filter(!paste(kingdom, rank) %in% paste(taxonomy$kingdom, taxonomy$rank)) %>%
  left_join(
    current_gbif %>%
      select(kingdom, rank = taxonRank, ref = scientificNameAuthorship, gbif = taxonID, gbif_parent = parentNameUsageID),
    by = c("kingdom", "rank")
  ) %>%
  mutate(source = if_else(!is.na(gbif), "GBIF", "manually added"),
         status = if_else(!is.na(gbif), "accepted", "unknown"))

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
      rank = i_name
    ) %>%
    filter(!paste(kingdom, .[[ncol(.) - 2]], rank) %in% paste(taxonomy$kingdom, taxonomy[[i + 1]], taxonomy$rank)) %>%
    # get GBIF identifier where available
    left_join(
      current_gbif %>%
        select(kingdom, all_of(i_name), rank = taxonRank, ref = scientificNameAuthorship, gbif = taxonID, gbif_parent = parentNameUsageID),
      by = c("kingdom", "rank", i_name)
    ) %>%
    mutate(source = if_else(!is.na(gbif), "GBIF", "manually added"),
           status = if_else(!is.na(gbif), "accepted", "unknown"))
  message("n = ", nrow(to_add))
  taxonomy_all_missing <- taxonomy_all_missing %>%
    bind_rows(to_add)
  rm(to_add)
}
taxonomy_all_missing %>% View()

taxonomy <- taxonomy %>%
  bind_rows(taxonomy_all_missing)

# fix for duplicate fullnames within a kingdom (such as Nitrospira which is the name of the genus AND its class)
taxonomy <- taxonomy %>%
  mutate(
    rank_index = case_when(
      rank == "subspecies" ~ 1,
      rank == "species" ~ 2,
      rank == "genus" ~ 3,
      rank == "family" ~ 4,
      rank == "order" ~ 5,
      rank == "class" ~ 6,
      TRUE ~ 7
    ),
    fullname_rank = paste0(fullname, " {", rank, "}")
  ) %>%
  arrange(kingdom, fullname, rank_index) %>%
  group_by(kingdom, fullname) %>%
  mutate(fullname = if_else(row_number() > 1, fullname_rank, fullname)) %>%
  ungroup() %>%
  select(-fullname_rank, -rank_index) %>%
  arrange(fullname)

# now also add missing species that have subspecies (requires combination with genus)
missing_species <- taxonomy %>%
  filter(species != "") %>%
  distinct(kingdom, genus, species, .keep_all = TRUE) %>%
  select(kingdom:species) %>%
  mutate(
    fullname = paste(genus, species),
    rank = "species"
  ) %>%
  filter(!paste(kingdom, genus, species, rank) %in% paste(taxonomy$kingdom, taxonomy$genus, taxonomy$species, taxonomy$rank)) %>%
  # get GBIF identifier where available
  left_join(
    current_gbif %>%
      select(kingdom, genus, species = specificEpithet, rank = taxonRank, ref = scientificNameAuthorship, gbif = taxonID, gbif_parent = parentNameUsageID),
    by = c("kingdom", "rank", "genus", "species")
  ) %>%
  mutate(source = if_else(!is.na(gbif), "GBIF", "manually added"),
         status = if_else(!is.na(gbif), "accepted", "unknown"))

taxonomy <- taxonomy %>%
  bind_rows(missing_species)

# remove NAs from taxonomy again, and keep unique full names
taxonomy <- taxonomy %>%
  mutate(across(kingdom:subspecies, function(x) if_else(is.na(x), "", x))) %>%
  arrange(kingdom, fullname, ref) %>%
  distinct(kingdom, fullname, .keep_all = TRUE) %>%
  filter(kingdom != "")


# Save intermediate results (1) -------------------------------------------------------------------

saveRDS(taxonomy, "data-raw/taxonomy1.rds")
# taxonomy <- readRDS("data-raw/taxonomy1.rds")


# Get previously manually added entries -----------------------------------------------------------

manually_added <- AMR::microorganisms %>%
  filter(source == "manually added",
         !paste(kingdom, fullname) %in% paste(taxonomy$kingdom, taxonomy$fullname),
         !rank %in% c("kingdom", "phylum", "class", "order", "family")) %>%
  select(fullname:subspecies, ref, source, rank)

# get latest taxonomy for those entries
`%if_na%` <- function(x, y) if (is.na(x)) y else x
for (g in unique(manually_added$genus[manually_added$genus != "" & manually_added$genus %in% taxonomy$genus])) {
  manually_added$family[which(manually_added$genus == g)] <-
    taxonomy$family[which(taxonomy$genus == g & !is.na(taxonomy$lpsn))][1] %if_na%
    taxonomy$family[which(taxonomy$genus == g & !is.na(taxonomy$mycobank))][1] %if_na%
    taxonomy$family[which(taxonomy$genus == g)][1]
}
for (f in unique(manually_added$family[manually_added$family != "" & manually_added$family %in% taxonomy$family])) {
  manually_added$order[which(manually_added$family == f)] <- 
    taxonomy$order[which(taxonomy$family == f & !is.na(taxonomy$lpsn))][1] %if_na%
    taxonomy$order[which(taxonomy$family == f & !is.na(taxonomy$mycobank))][1] %if_na%
    taxonomy$order[which(taxonomy$family == f)][1]
}
for (o in unique(manually_added$order[manually_added$order != "" & manually_added$order %in% taxonomy$order])) {
  manually_added$class[which(manually_added$order == o)] <- 
    taxonomy$class[which(taxonomy$order == o & !is.na(taxonomy$lpsn))][1] %if_na%
    taxonomy$class[which(taxonomy$order == o & !is.na(taxonomy$mycobank))][1] %if_na%
    taxonomy$class[which(taxonomy$order == o)][1]
}
for (cc in unique(manually_added$class[manually_added$class != "" & manually_added$class %in% taxonomy$class])) {
  manually_added$phylum[which(manually_added$class == cc)] <- 
    taxonomy$phylum[which(taxonomy$class == cc & !is.na(taxonomy$lpsn))][1] %if_na%
    taxonomy$phylum[which(taxonomy$class == cc & !is.na(taxonomy$mycobank))][1] %if_na%
    taxonomy$phylum[which(taxonomy$class == cc)][1]
}
for (p in unique(manually_added$phylum[manually_added$phylum != "" & manually_added$phylum %in% taxonomy$phylum])) {
  manually_added$kingdom[which(manually_added$phylum == p)] <- 
    taxonomy$kingdom[which(taxonomy$phylum == p & !is.na(taxonomy$lpsn))][1] %if_na%
    taxonomy$kingdom[which(taxonomy$phylum == p & !is.na(taxonomy$mycobank))][1] %if_na%
    taxonomy$kingdom[which(taxonomy$phylum == p)][1]
}

manually_added <- manually_added %>%
  mutate(
    status = "unknown",
    rank = if_else(fullname %like% "unknown", "(unknown rank)", rank)
  )
manually_added %>% View()

# these are now included in the new taxonomy, check them
manually_added %>% filter(fullname %in% taxonomy$fullname)

taxonomy <- taxonomy %>%
  # here also the 'unknowns' are added, such as "(unknown fungus)"
  bind_rows(manually_added) %>%
  arrange(fullname)

table(taxonomy$rank, useNA = "always")


# Get LPSN data for records missing from `taxonomy_lpsn` ------------------------------------------

# Weirdly enough, some LPSN records are lacking from the API and the CSV file (i.e., `taxonomy_lpsn`), 
# such as family Thiotrichaceae and its order Thiotrichales, or the genus Coleospermum. When running
# get_lpsn_and_author("family", "Thiotrichaceae") you do get a result, vs. taxonomy_lpsn %>% filter(family == "Thiotrichaceae").
# So check every non-LPSN records from the kingdom of Bacteria and add it
gbif_bacteria <- which(taxonomy$kingdom == "Bacteria" & taxonomy$source %in% c("GBIF", "manually added") & taxonomy$rank %in% c("phylum", "class", "order", "family", "genus"))
added <- 0
pb <- progress_bar$new(total = length(gbif_bacteria), format = "[:bar] :current/:total :eta")
for (record in gbif_bacteria) {
  pb$tick()
  lpsn <- get_lpsn_and_author(rank = taxonomy$rank[record],
                              name = taxonomy$fullname[record])
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

saveRDS(taxonomy, "data-raw/taxonomy1b.rds")
# taxonomy <- readRDS("data-raw/taxonomy1b.rds")


# Clean scientific reference ----------------------------------------------------------------------

taxonomy <- taxonomy %>%
  mutate(ref = get_author_year(ref))


# Get the latest upper taxonomy from LPSN/MycoBank for GBIF data ---------------------------------------

# we did this for the manually added ones (from the current AMR pkg version), but not for the new GBIF records
# (e.g., phylum above class "Bacilli" was still "Firmicutes" in 2023, should be "Bacillota")
for (k in unique(taxonomy$kingdom[taxonomy$kingdom != ""])) {
  if (k == "Fungi") {
    src <- "MycoBank"
  } else {
    src <- "LPSN"
  }
  message("Fixing GBIF taxonomy for kingdom ", k, " based on ", src, ".", appendLF = FALSE)
  i <- 0
  for (g in unique(taxonomy$genus[taxonomy$genus != "" & taxonomy$kingdom == k & taxonomy$source == src])) {
    i <- i + 1
    if (i %% 50 == 0) message(".", appendLF = FALSE)
    taxonomy$family[which(taxonomy$genus == g & taxonomy$kingdom == k)] <- taxonomy$family[which(taxonomy$genus == g & taxonomy$kingdom == k & taxonomy$source == src)][1]
  }
  for (f in unique(taxonomy$family[taxonomy$family != "" & taxonomy$kingdom == k & taxonomy$source == src])) {
    i <- i + 1
    if (i %% 50 == 0) message(".", appendLF = FALSE)
    taxonomy$order[which(taxonomy$family == f & taxonomy$kingdom == k)] <- taxonomy$order[which(taxonomy$family == f & taxonomy$kingdom == k & taxonomy$source == src)][1]
  }
  for (o in unique(taxonomy$order[taxonomy$order != "" & taxonomy$kingdom == k & taxonomy$source == src])) {
    i <- i + 1
    if (i %% 50 == 0) message(".", appendLF = FALSE)
    taxonomy$class[which(taxonomy$order == o & taxonomy$kingdom == k)] <- taxonomy$class[which(taxonomy$order == o & taxonomy$kingdom == k & taxonomy$source == src)][1]
  }
  for (cc in unique(taxonomy$class[taxonomy$class != "" & taxonomy$kingdom == k & taxonomy$source == src])) {
    i <- i + 1
    if (i %% 50 == 0) message(".", appendLF = FALSE)
    taxonomy$phylum[which(taxonomy$class == cc & taxonomy$kingdom == k)] <- taxonomy$phylum[which(taxonomy$class == cc & taxonomy$kingdom == k & taxonomy$source == src)][1]
  }
  message("OK.")
}


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

saveRDS(taxonomy, "data-raw/taxonomy1c.rds")
# taxonomy <- readRDS("data-raw/taxonomy1c.rds")


# Add parent identifiers --------------------------------------------------------------------------

# requires full name and full taxonomy
taxonomy <- taxonomy %>%
  mutate(
    lpsn_parent = case_when(
      rank == "phylum" ~ lpsn[match(kingdom, fullname)],
      # in class, take parent from phylum if available, otherwise kingdom
      rank == "class" & phylum != "" ~ lpsn[match(phylum, fullname)],
      rank == "class" ~ lpsn[match(kingdom, fullname)],
      # in order, take parent from class if available, otherwise phylum, otherwise kingdom
      rank == "order" & class != "" ~ lpsn[match(class, fullname)],
      rank == "order" & phylum != "" ~ lpsn[match(phylum, fullname)],
      rank == "order" ~ lpsn[match(kingdom, fullname)],
      # family
      rank == "family" & order != "" ~ lpsn[match(order, fullname)],
      rank == "family" & class != "" ~ lpsn[match(class, fullname)],
      rank == "family" & phylum != "" ~ lpsn[match(phylum, fullname)],
      rank == "family" ~ lpsn[match(kingdom, fullname)],
      # genus
      rank == "genus" & family != "" ~ lpsn[match(family, fullname)],
      rank == "genus" & order != "" ~ lpsn[match(order, fullname)],
      rank == "genus" & class != "" ~ lpsn[match(class, fullname)],
      rank == "genus" & phylum != "" ~ lpsn[match(phylum, fullname)],
      rank == "genus" ~ lpsn[match(kingdom, fullname)],
      # species, always has a genus
      rank == "species" ~ lpsn[match(genus, fullname)],
      # subspecies, always has a genus + species
      rank == "subspecies" ~ lpsn[match(paste(genus, species), fullname)],
      TRUE ~ NA_character_),
    mycobank_parent = case_when(
      rank == "phylum" ~ mycobank[match(kingdom, fullname)],
      # class
      rank == "class" & phylum != "" ~ mycobank[match(phylum, fullname)],
      rank == "class" ~ mycobank[match(kingdom, fullname)],
      # order
      rank == "order" & class != "" ~ mycobank[match(class, fullname)],
      rank == "order" & phylum != "" ~ mycobank[match(phylum, fullname)],
      rank == "order" ~ mycobank[match(kingdom, fullname)],
      # family
      rank == "family" & order != "" ~ mycobank[match(order, fullname)],
      rank == "family" & class != "" ~ mycobank[match(class, fullname)],
      rank == "family" & phylum != "" ~ mycobank[match(phylum, fullname)],
      rank == "family" ~ mycobank[match(kingdom, fullname)],
      # genus
      rank == "genus" & family != "" ~ mycobank[match(family, fullname)],
      rank == "genus" & order != "" ~ mycobank[match(order, fullname)],
      rank == "genus" & class != "" ~ mycobank[match(class, fullname)],
      rank == "genus" & phylum != "" ~ mycobank[match(phylum, fullname)],
      rank == "genus" ~ mycobank[match(kingdom, fullname)],
      # species
      rank == "species" ~ mycobank[match(genus, fullname)],
      # subspecies
      rank == "subspecies" ~ mycobank[match(paste(genus, species), fullname)],
      TRUE ~ NA_character_),
    gbif_parent = case_when(
      rank == "phylum" ~ gbif[match(kingdom, fullname)],
      # class
      rank == "class" & phylum != "" ~ gbif[match(phylum, fullname)],
      rank == "class" ~ gbif[match(kingdom, fullname)],
      # order
      rank == "order" & class != "" ~ gbif[match(class, fullname)],
      rank == "order" & phylum != "" ~ gbif[match(phylum, fullname)],
      rank == "order" ~ gbif[match(kingdom, fullname)],
      # family
      rank == "family" & order != "" ~ gbif[match(order, fullname)],
      rank == "family" & class != "" ~ gbif[match(class, fullname)],
      rank == "family" & phylum != "" ~ gbif[match(phylum, fullname)],
      rank == "family" ~ gbif[match(kingdom, fullname)],
      # genus
      rank == "genus" & family != "" ~ gbif[match(family, fullname)],
      rank == "genus" & order != "" ~ gbif[match(order, fullname)],
      rank == "genus" & class != "" ~ gbif[match(class, fullname)],
      rank == "genus" & phylum != "" ~ gbif[match(phylum, fullname)],
      rank == "genus" ~ gbif[match(kingdom, fullname)],
      # species
      rank == "species" ~ gbif[match(genus, fullname)],
      # subspecies
      rank == "subspecies" ~ gbif[match(paste(genus, species), fullname)],
      TRUE ~ NA_character_))

# these still have no record in our data set:
which(!taxonomy$lpsn_parent %in% taxonomy$lpsn)
which(!taxonomy$mycobank_parent %in% taxonomy$mycobank)
which(!taxonomy$gbif_parent %in% taxonomy$gbif)



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
nonbacterial_genera <- nonbacterial_genera[nonbacterial_genera %unlike% "unknown"]

# update prevalence based on taxonomy (following the recent and thorough work of Bartlett et al., 2022)
# see https://doi.org/10.1099/mic.0.001269
taxonomy <- taxonomy %>%
  mutate(prevalence = case_when(
    # genera of the pathogens mentioned in the World Health Organization's (WHO) Priority Pathogen List
    genus %in% MO_WHO_PRIORITY_GENERA ~ 1.0,
    # 'established' means 'have infected at least three persons in three or more references'
    paste(genus, species) %in% established & rank %in% c("species", "subspecies") ~ 1.15,
    # other genera in the 'established' group
    genus %in% established_genera & rank == "genus" ~ 1.15,
    
    # 'putative' means 'fewer than three known cases'
    paste(genus, species) %in% putative & rank %in% c("species", "subspecies") ~ 1.25,
    # other genera in the 'putative' group
    genus %in% putative_genera & rank == "genus" ~ 1.25,
    
    # we keep track of prevalent genera too of non-bacterial species
    genus %in% AMR:::MO_RELEVANT_GENERA & kingdom != "Bacteria" & rank %in% c("genus", "species", "subspecies") ~ 1.25,
    
    # species and subspecies in 'established' and 'putative' groups
    genus %in% c(established_genera, putative_genera) & rank %in% c("species", "subspecies") ~ 1.5,
    # other species from a genus in either group
    genus %in% nonbacterial_genera & rank %in% c("genus", "species", "subspecies") ~ 1.5,
    
    # all others
    TRUE ~ 2.0
  ))

table(taxonomy$prevalence, useNA = "always")
# (a lot will be removed further below)


# Save intermediate results (2) -------------------------------------------------------------------

saveRDS(taxonomy, "data-raw/taxonomy2.rds")
# taxonomy <- readRDS("data-raw/taxonomy2.rds")


# Remove unwanted taxonomic entries ---------------------------------------------------------------

part1 <- taxonomy %>%
  filter(
    # keep all we added ourselves:
    source == "manually added" |
      # keep all bacteria anyway, main focus of our package:
      kingdom == "Bacteria" |
      # and these kingdoms are very small, and also mainly microorganisms:
      kingdom == "Protozoa" |
      kingdom == "Archaea" |
      kingdom == "Chromista" |
      # keep everything from family up, 'ghost' entries will be removed later on
      rank %in% c("kingdom", "phylum", "class", "order", "family") |
      # other relevant genera to keep:
      genus %in% AMR:::MO_RELEVANT_GENERA |
      # relevant for biotechnology:
      genus %in% c("Archaeoglobus", "Desulfurococcus", "Ferroglobus", "Ferroplasma", "Halobacterium", "Halococcus", "Haloferax", "Naegleria", "Nanoarchaeum", "Nitrosopumilus", "Nosema", "Pleistophora", "Pyrobaculum", "Pyrococcus", "Spraguea", "Thelohania", "Thermococcus", "Thermoplasma", "Thermoproteus", "Tritrichomonas") |
      genus %like_case% "Methano" |
      # kingdom of Protozoa:
      (phylum %in% c("Choanozoa", "Mycetozoa") & prevalence < 2) |
      # Fungi:
      (kingdom == "Fungi" & (!rank %in% c("genus", "species", "subspecies") | prevalence < 2 | class == "Pichiomycetes")) |
      # Animalia:
      genus %in% c("Lucilia", "Lumbricus") |
      (class == "Insecta" & !rank %in% c("species", "subspecies")) | # keep only genus of insects, not all of their (sub)species
      (genus == "Amoeba" & kingdom != "Animalia") # keep only in the protozoa, not the animalia
  ) %>%
  # this kingdom only contained Curvularia and Hymenolepis, which have coincidental twin names with Fungi
  filter(kingdom != "Plantae",
         !(genus %in% c("Aedes", "Anopheles") & rank %in% c("species", "subspecies")))

# now get the parents and old names
part2 <- taxonomy %>% 
  filter(gbif %in% c(part1$gbif_parent[!is.na(part1$gbif_parent)], part1$gbif_renamed_to[!is.na(part1$gbif_renamed_to)]) |
           mycobank %in% c(part1$mycobank_parent[!is.na(part1$mycobank_parent)], part1$mycobank_renamed_to[!is.na(part1$mycobank_renamed_to)]) |
           lpsn %in% c(part1$lpsn_parent[!is.na(part1$lpsn_parent)], part1$lpsn_renamed_to[!is.na(part1$lpsn_renamed_to)]))
parts <- bind_rows(part1, part2)

part3 <- taxonomy %>% 
  filter(gbif %in% c(parts$gbif_parent[!is.na(parts$gbif_parent)], parts$gbif_renamed_to[!is.na(parts$gbif_renamed_to)]) |
           mycobank %in% c(parts$mycobank_parent[!is.na(parts$mycobank_parent)], parts$mycobank_renamed_to[!is.na(parts$mycobank_renamed_to)]) |
           lpsn %in% c(parts$lpsn_parent[!is.na(parts$lpsn_parent)], parts$lpsn_renamed_to[!is.na(parts$lpsn_renamed_to)]))
parts <- bind_rows(part1, part2, part3)

part4 <- taxonomy %>% 
  filter(gbif %in% c(parts$gbif_parent[!is.na(parts$gbif_parent)], parts$gbif_renamed_to[!is.na(parts$gbif_renamed_to)]) |
           mycobank %in% c(parts$mycobank_parent[!is.na(parts$mycobank_parent)], parts$mycobank_renamed_to[!is.na(parts$mycobank_renamed_to)]) |
           lpsn %in% c(parts$lpsn_parent[!is.na(parts$lpsn_parent)], parts$lpsn_renamed_to[!is.na(parts$lpsn_renamed_to)]))
parts <- bind_rows(part1, part2, part3, part4)

taxonomy <- bind_rows(part1, part2, part3, part4) %>% 
  mutate(
    # here we must prefer in this order: LPSN > MycoBank > GBIF
    source_index = case_when(source == "LPSN" ~ 1,
                             source == "MycoBank" ~ 2,
                             source == "GBIF" ~ 3,
                             # manually added:
                             TRUE ~ 4),
    # also arrange on rank, otherwise e.g. Cyptococcus comes from Animalia, not Fungi
    rank_index = case_when(
      kingdom == "Bacteria" ~ 1,
      kingdom == "Fungi" ~ 2,
      kingdom == "Protozoa" ~ 3,
      kingdom == "Archaea" ~ 4,
      kingdom == "Chromista" ~ 5,
      kingdom == "Animalia" ~ 6,
      TRUE ~ 7)) %>%
  arrange(source_index, rank_index, fullname) %>% 
  distinct(fullname, .keep_all = TRUE) %>% 
  select(-c(source_index, rank_index))

# some manual updates, remove later
taxonomy$family[taxonomy$genus == "Blastocystis"] <- "Blastocystidae"
taxonomy$ref[which(taxonomy$fullname == "Salmonella Typhimurium")] <- "Castellani et al., 1919"
taxonomy$ref[which(taxonomy$fullname == "Salmonella Paratyphi")] <- "Ezaki et al., 2000"
taxonomy$ref[which(taxonomy$fullname == "Salmonella Typhi")] <- "Warren et al., 1930"

# first make sure that species of a genus cannot be across multiple kingdoms, since otherwise IDs cannot be given correctly
# (e.g., Amoeba has some species in Protozoa and some in Animalia)
# we acknowledge that this may be taxonomically right, but we need the same KINGDOM_GENUS identifier for each record
taxonomy <- taxonomy %>%
  group_by(genus) %>%
  mutate(
    kingdom = if_else(rank != "genus" & genus != "" & fullname %unlike% "unknown", first(kingdom[rank == "genus"]), kingdom),
    phylum = if_else(rank != "genus" & genus != "" & fullname %unlike% "unknown", first(phylum[rank == "genus"]), phylum),
    class = if_else(rank != "genus" & genus != "" & fullname %unlike% "unknown", first(class[rank == "genus"]), class),
    order = if_else(rank != "genus" & genus != "" & fullname %unlike% "unknown", first(order[rank == "genus"]), order),
    family = if_else(rank != "genus" & genus != "" & fullname %unlike% "unknown", first(family[rank == "genus"]), family),
    prevalence = if_else(rank != "genus" & genus != "" & fullname %unlike% "unknown", first(prevalence[rank == "genus"]), prevalence)
  ) %>%
  ungroup() %>% 
  arrange(fullname)

# update the taxonomic names based on new genus classification
taxonomy <- taxonomy %>% 
  # first, fix kingdom based on family level (only when family is not empty)
  group_by(family) %>%
  mutate(kingdom = if_else(family != "", get_top_lvl(kingdom, rank, source, "phylum", phylum), kingdom),
         phylum = if_else(family != "", first(phylum[phylum != "" & rank != "phylum"]), phylum),
         class = if_else(family != "", first(class[class != "" & rank != "class"]), class),
         order = if_else(family != "", first(order[order != "" & rank != "order"]), order)) %>% 
  # next, update all taxonomic layers
  group_by(kingdom, genus) %>% 
  mutate(family = get_top_lvl(family, rank, source, "genus", genus)) %>% 
  group_by(kingdom, family) %>%
  mutate(order = get_top_lvl(order, rank, source, "family", family)) %>% 
  group_by(kingdom, order) %>%
  mutate(class = get_top_lvl(class, rank, source, "order", order)) %>% 
  group_by(kingdom, class) %>%
  mutate(phylum = get_top_lvl(phylum, rank, source, "class", class)) %>% 
  ungroup() %>% 
  arrange(fullname) %>% 
  mutate(across(phylum:family, ~if_else(is.na(.x), "", .x)))

# no ghost families, orders, classes, phyla
# (but keep the ghost families of bacteria)
taxonomy <- taxonomy %>%
  group_by(kingdom, family) %>%
  filter(n() > 1 | fullname %like% "unknown" | rank %in% c("kingdom", "genus", "species", "subspecies") | kingdom == "Bacteria") %>%
  group_by(kingdom, order) %>%
  filter(n() > 1 | fullname %like% "unknown" | rank %in% c("kingdom", "genus", "species", "subspecies") | family != "") %>%
  group_by(kingdom, class) %>%
  filter(n() > 1 | fullname %like% "unknown" | rank %in% c("kingdom", "genus", "species", "subspecies") | family != "" | order != "") %>%
  group_by(kingdom, phylum) %>%
  filter(n() > 1 | fullname %like% "unknown" | rank %in% c("kingdom", "genus", "species", "subspecies") | family != "" | order != "" | class != "") %>%
  ungroup()

saveRDS(taxonomy, "data-raw/taxonomy2b.rds")
# taxonomy <- readRDS("data-raw/taxonomy2b.rds")


# Add microbial IDs -------------------------------------------------------------------------------

# (MO codes in the AMR package have the form: KINGDOM_GENUS_SPECIES_SUBSPECIES where all are abbreviated)

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
    kingdom == "Plantae" ~ "PL", # this is actually not part of the scope, having 0 results (2024)
    kingdom == "Protozoa" ~ "P",
    TRUE ~ ""
  ))

# phylum until family are abbreviated with 8 characters and prefixed with their rank

# Phylum - keep old and fill up for new ones
mo_phylum <- taxonomy %>%
  filter(rank == "phylum") %>%
  distinct(kingdom, phylum) %>%
  left_join(
    AMR::microorganisms %>%
      filter(rank == "phylum") %>%
      transmute(kingdom,
                phylum = fullname,
                mo_old = gsub("[A-Z]{1,2}_", "", as.character(mo))
      ),
    by = c("kingdom", "phylum")
  ) %>%
  group_by(kingdom) %>%
  mutate(
    mo_phylum8 = AMR:::abbreviate_mo(phylum, minlength = 8, prefix = "[PHL]_"),
    mo_phylum9 = AMR:::abbreviate_mo(phylum, minlength = 9, prefix = "[PHL]_"),
    mo_phylum = if_else(!is.na(mo_old), mo_old, mo_phylum8),
    mo_duplicated = duplicated(mo_phylum),
    mo_phylum = if_else(mo_duplicated, mo_phylum9, mo_phylum),
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
  left_join(
    AMR::microorganisms %>%
      filter(rank == "class") %>%
      transmute(kingdom,
                class = fullname,
                mo_old = gsub("[A-Z]{1,2}_", "", as.character(mo))
      ),
    by = c("kingdom", "class")
  ) %>%
  group_by(kingdom) %>%
  mutate(
    mo_class8 = AMR:::abbreviate_mo(class, minlength = 8, prefix = "[CLS]_"),
    mo_class9 = AMR:::abbreviate_mo(class, minlength = 9, prefix = "[CLS]_"),
    mo_class = if_else(!is.na(mo_old), mo_old, mo_class8),
    mo_duplicated = duplicated(mo_class),
    mo_class = if_else(mo_duplicated, mo_class9, mo_class),
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
  left_join(
    AMR::microorganisms %>%
      filter(rank == "order") %>%
      transmute(kingdom,
                order = fullname,
                mo_old = gsub("[A-Z]{1,2}_", "", as.character(mo))
      ),
    by = c("kingdom", "order")
  ) %>%
  group_by(kingdom) %>%
  mutate(
    mo_order8 = AMR:::abbreviate_mo(order, minlength = 8, prefix = "[ORD]_"),
    mo_order9 = AMR:::abbreviate_mo(order, minlength = 9, prefix = "[ORD]_"),
    mo_order = if_else(!is.na(mo_old), mo_old, mo_order8),
    mo_duplicated = duplicated(mo_order),
    mo_order = if_else(mo_duplicated, mo_order9, mo_order),
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
  left_join(
    AMR::microorganisms %>%
      filter(rank == "family") %>%
      transmute(kingdom,
                family = fullname,
                mo_old = gsub("[A-Z]{1,2}_", "", as.character(mo))
      ),
    by = c("kingdom", "family")
  ) %>%
  group_by(kingdom) %>%
  mutate(
    mo_family8 = AMR:::abbreviate_mo(family, minlength = 8, prefix = "[FAM]_"),
    mo_family9 = AMR:::abbreviate_mo(family, minlength = 9, prefix = "[FAM]_"),
    mo_family = if_else(!is.na(mo_old), mo_old, mo_family8),
    mo_duplicated = duplicated(mo_family),
    mo_family = if_else(mo_duplicated, mo_family9, mo_family),
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
  left_join(
    AMR::microorganisms %>%
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
if (any(mo_genus$mo_duplicated, na.rm = TRUE) | anyNA(mo_genus$mo_genus_new)) stop("Duplicate MO codes for genus!")
# mo_genus %>% filter(mo_duplicated)
# no duplicates *within kingdoms*, so keep the right columns for left joining later
mo_genus <- mo_genus %>%
  select(kingdom, genus, mo_genus = mo_genus_new)

# same for species - keep old where available and create new per kingdom-genus where needed:
mo_species <- taxonomy %>%
  filter(rank == "species") %>%
  distinct(kingdom, genus, species) %>%
  left_join(
    AMR::microorganisms %>%
      filter(rank == "species") %>%
      transmute(mo_species_old = gsub("^[A-Z]+_[A-Z]+_", "", as.character(mo)), kingdom, genus, species) %>%
      filter(mo_species_old %unlike% "-") %>%
      distinct(kingdom, genus, species, .keep_all = TRUE),
    by = c("kingdom", "genus", "species")
  ) %>%
  distinct(kingdom, genus, species, .keep_all = TRUE) %>%
  group_by(kingdom, genus) %>%
  mutate(
    mo_species_new4 = AMR:::abbreviate_mo(species, 4, hyphen_as_space = TRUE),
    mo_species_new5 = AMR:::abbreviate_mo(species, 5, hyphen_as_space = TRUE),
    mo_species_new5b = paste0(AMR:::abbreviate_mo(species, 5, hyphen_as_space = TRUE), 1),
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
if (any(mo_species$mo_duplicated, na.rm = TRUE) | anyNA(mo_species$mo_species_new)) stop("Duplicate MO codes for species!")
# no duplicates *within kingdoms*, so keep the right columns for left joining later
mo_species <- mo_species %>%
  select(kingdom, genus, species, mo_species = mo_species_new)

# same for subspecies - keep old where available and create new per kingdom-genus-species where needed:
mo_subspecies <- taxonomy %>%
  filter(rank == "subspecies") %>%
  distinct(kingdom, genus, species, subspecies) %>%
  left_join(
    AMR::microorganisms %>%
      filter(rank %in% c("subspecies", "subsp.", "infraspecies")) %>%
      transmute(mo_subspecies_old = gsub("^[A-Z]+_[A-Z]+_[A-Z]+_", "", as.character(mo)), kingdom, genus, species, subspecies) %>%
      filter(mo_subspecies_old %unlike% "-") %>%
      distinct(kingdom, genus, species, subspecies, .keep_all = TRUE),
    by = c("kingdom", "genus", "species", "subspecies")
  ) %>%
  distinct(kingdom, genus, species, subspecies, .keep_all = TRUE) %>%
  group_by(kingdom, genus, species) %>%
  mutate(
    mo_subspecies_new4 = AMR:::abbreviate_mo(subspecies, 4, hyphen_as_space = TRUE),
    mo_subspecies_new5 = AMR:::abbreviate_mo(subspecies, 5, hyphen_as_space = TRUE),
    mo_subspecies_new5b = paste0(AMR:::abbreviate_mo(subspecies, 5, hyphen_as_space = TRUE), 1),
    mo_subspecies_new6 = AMR:::abbreviate_mo(subspecies, 6, hyphen_as_space = TRUE),
    mo_subspecies_new7 = AMR:::abbreviate_mo(subspecies, 7, hyphen_as_space = TRUE),
    mo_subspecies_new8 = AMR:::abbreviate_mo(subspecies, 8, hyphen_as_space = TRUE),
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
  mutate(across(starts_with("mo_"), function(x) if_else(is.na(x), "", x))) %>%
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
  arrange(fullname)

# now check these duplicate fullnames, they should not exist (didn't at least in 2024)
# this example shows that some could have same fullnames: taxonomy %>% filter(class == genus & class != "")
# such as Kapibacteria (both genus and class), but the class last time correctly got the " {class}" suffix in the fullname
taxonomy %>%
  filter(fullname %in% .[duplicated(fullname), "fullname", drop = TRUE]) %>%
  View()
# OTHERWISE uncomment this section:
# # we will prefer kingdoms as Bacteria > Fungi > Protozoa > Archaea > Animalia, so check the ones within 1 kingdom
# taxonomy %>%
#   filter(fullname %in% .[duplicated(fullname), "fullname", drop = TRUE]) %>%
#   group_by(fullname) %>% 
#   filter(n_distinct(kingdom) == 1)
# 
# # fullnames must be unique, we'll keep the most relevant ones only
# taxonomy <- taxonomy %>%
#   mutate(rank_index = case_when(
#     kingdom == "Bacteria" ~ 1,
#     kingdom == "Fungi" ~ 2,
#     kingdom == "Protozoa" ~ 3,
#     kingdom == "Archaea" ~ 4,
#     kingdom == "Chromista" ~ 5,
#     kingdom == "Animalia" ~ 6,
#     TRUE ~ 7
#   )) %>%
#   arrange(fullname, rank_index) %>%
#   distinct(fullname, .keep_all = TRUE) %>%
#   select(-rank_index) %>%
#   filter(mo != "")

# These were overwritten to Bacteria, fix them
taxonomy$kingdom[taxonomy$mo == "UNKNOWN"] <- "(unknown kingdom)"
taxonomy$kingdom[taxonomy$mo == "F_FUNGUS"] <- "Fungi"
taxonomy$kingdom[taxonomy$mo == "F_YEAST"] <- "Fungi"
taxonomy$kingdom[taxonomy$mo == "P_PROTOZOAN"] <- "Protozoa"


# keep the codes from manually added ones
manual_mos <- as.character(AMR::microorganisms$mo)[match(taxonomy$fullname[taxonomy$source == "manually added"], AMR::microorganisms$fullname)]
taxonomy$mo[taxonomy$source == "manually added"][!is.na(manual_mos)] <- manual_mos[!is.na(manual_mos)]
# fix the kingdoms for manual codes (sometimes species go from Fungi to Chromista or so)
substr(taxonomy$mo[!substr(taxonomy$mo, 1, 1) %in% c("B", "U") & taxonomy$kingdom == "Bacteria" & !is.na(taxonomy$kingdom)], 1, 1) <- "B"
substr(taxonomy$mo[!substr(taxonomy$mo, 1, 1) %in% c("F", "U") & taxonomy$kingdom == "Fungi" & !is.na(taxonomy$kingdom)], 1, 1) <- "F"
substr(taxonomy$mo[!substr(taxonomy$mo, 1, 1) %in% c("C", "U") & taxonomy$kingdom == "Chromista" & !is.na(taxonomy$kingdom)], 1, 1) <- "C"
substr(taxonomy$mo[!substr(taxonomy$mo, 1, 1) %in% c("A", "U") & taxonomy$kingdom == "Archaea" & !is.na(taxonomy$kingdom)], 1, 1) <- "A"
substr(taxonomy$mo[!substr(taxonomy$mo, 1, 1) %in% c("P", "U") & taxonomy$kingdom == "Protozoa" & !is.na(taxonomy$kingdom)], 1, 1) <- "P"
substr(taxonomy$mo[!substr(taxonomy$mo, 1, 2) %in% c("AN", "UN") & taxonomy$kingdom == "Animalia" & !is.na(taxonomy$kingdom)], 1, 2) <- "AN"
substr(taxonomy$mo[!substr(taxonomy$mo, 1, 2) %in% c("PL", "UN") & taxonomy$kingdom == "Plantae" & !is.na(taxonomy$kingdom)], 1, 2) <- "PL"

# remove the manually added ones that are now ghosts
taxonomy %>%
  filter((mo %like% "__" & source == "manually added")) %>% 
  View()
taxonomy <- taxonomy %>%
  filter(!(mo %like% "__" & source == "manually added"))

# this must not exist (2024: GBIF mess with missing genera - remove them):
taxonomy %>%
  filter(mo %like% "__") %>%
  View()
taxonomy <- taxonomy %>% filter(mo %unlike% "__")

saveRDS(taxonomy, "data-raw/taxonomy2c.rds")
# taxonomy <- readRDS("data-raw/taxonomy2c.rds")


# Some integrity checks ---------------------------------------------------------------------------

# are mo codes unique?
taxonomy %>% filter(mo %in% .[duplicated(mo), "mo", drop = TRUE]) %>% arrange(mo) %>% View()
# Nope: (2024)
# - There is one weird entry from MycoBank, this is recorded as class while it's a phylum
fix <- which(taxonomy$fullname == "Chytridiopsidomycota")
taxonomy$class[fix] <- ""
taxonomy$rank[fix] <- "phylum"
taxonomy$mo[fix] <- paste0("F_", AMR:::abbreviate_mo("Chytridiopsidomycota", minlength = 8, prefix = "[PHL]_"))
# Also:
# - There are multiple cases where the suffix " {species}" was added - these are just duplicates.
#   Some are Salmonella errors at LPSN's side - typhimurium is NOT a species of it, species is enterica, serovar is Typhimurium (that we use as subspecies)
taxonomy %>% filter(mo %in% .[duplicated(mo), "mo", drop = TRUE]) %>% arrange(mo) %>% View()
# keep the firsts
taxonomy <- taxonomy %>% arrange(mo) %>% distinct(mo, .keep_all = TRUE)

# are fullnames unique?
taxonomy %>% arrange(fullname) %>% filter(fullname %in% .[duplicated(fullname), "fullname", drop = TRUE]) %>% View()
# # we can now solve it like this, but check every year if this is okay
# taxonomy <- taxonomy %>%
#   group_by(fullname) %>%
#   mutate(fullname = if_else(duplicated(fullname), paste0(fullname, " {", rank, "}"), fullname)) %>%
#   ungroup()


# are all GBIFs available?
taxonomy %>%
  filter((!gbif_parent %in% gbif) | (!lpsn_parent %in% lpsn) | (!mycobank_parent %in% mycobank)) %>%
  count(source = case_when(!gbif_parent %in% gbif ~ "GBIF",
                           !lpsn_parent %in% lpsn ~ "LPSN",
                           !mycobank_parent %in% mycobank ~ "MycoBank",
                           TRUE ~ "?"),
        rank)

# so fix again all parent identifiers (exact same code as earlier)
taxonomy <- taxonomy %>%
  mutate(
    lpsn_parent = case_when(
      rank == "phylum" ~ lpsn[match(kingdom, fullname)],
      # in class, take parent from phylum if available, otherwise kingdom
      rank == "class" & phylum != "" ~ lpsn[match(phylum, fullname)],
      rank == "class" ~ lpsn[match(kingdom, fullname)],
      # in order, take parent from class if available, otherwise phylum, otherwise kingdom
      rank == "order" & class != "" ~ lpsn[match(class, fullname)],
      rank == "order" & phylum != "" ~ lpsn[match(phylum, fullname)],
      rank == "order" ~ lpsn[match(kingdom, fullname)],
      # family
      rank == "family" & order != "" ~ lpsn[match(order, fullname)],
      rank == "family" & class != "" ~ lpsn[match(class, fullname)],
      rank == "family" & phylum != "" ~ lpsn[match(phylum, fullname)],
      rank == "family" ~ lpsn[match(kingdom, fullname)],
      # genus
      rank == "genus" & family != "" ~ lpsn[match(family, fullname)],
      rank == "genus" & order != "" ~ lpsn[match(order, fullname)],
      rank == "genus" & class != "" ~ lpsn[match(class, fullname)],
      rank == "genus" & phylum != "" ~ lpsn[match(phylum, fullname)],
      rank == "genus" ~ lpsn[match(kingdom, fullname)],
      # species, always has a genus
      rank == "species" ~ lpsn[match(genus, fullname)],
      # subspecies, always has a genus + species
      rank == "subspecies" ~ lpsn[match(paste(genus, species), fullname)],
      TRUE ~ NA_character_),
    mycobank_parent = case_when(
      rank == "phylum" ~ mycobank[match(kingdom, fullname)],
      # class
      rank == "class" & phylum != "" ~ mycobank[match(phylum, fullname)],
      rank == "class" ~ mycobank[match(kingdom, fullname)],
      # order
      rank == "order" & class != "" ~ mycobank[match(class, fullname)],
      rank == "order" & phylum != "" ~ mycobank[match(phylum, fullname)],
      rank == "order" ~ mycobank[match(kingdom, fullname)],
      # family
      rank == "family" & order != "" ~ mycobank[match(order, fullname)],
      rank == "family" & class != "" ~ mycobank[match(class, fullname)],
      rank == "family" & phylum != "" ~ mycobank[match(phylum, fullname)],
      rank == "family" ~ mycobank[match(kingdom, fullname)],
      # genus
      rank == "genus" & family != "" ~ mycobank[match(family, fullname)],
      rank == "genus" & order != "" ~ mycobank[match(order, fullname)],
      rank == "genus" & class != "" ~ mycobank[match(class, fullname)],
      rank == "genus" & phylum != "" ~ mycobank[match(phylum, fullname)],
      rank == "genus" ~ mycobank[match(kingdom, fullname)],
      # species
      rank == "species" ~ mycobank[match(genus, fullname)],
      # subspecies
      rank == "subspecies" ~ mycobank[match(paste(genus, species), fullname)],
      TRUE ~ NA_character_),
    gbif_parent = case_when(
      rank == "phylum" ~ gbif[match(kingdom, fullname)],
      # class
      rank == "class" & phylum != "" ~ gbif[match(phylum, fullname)],
      rank == "class" ~ gbif[match(kingdom, fullname)],
      # order
      rank == "order" & class != "" ~ gbif[match(class, fullname)],
      rank == "order" & phylum != "" ~ gbif[match(phylum, fullname)],
      rank == "order" ~ gbif[match(kingdom, fullname)],
      # family
      rank == "family" & order != "" ~ gbif[match(order, fullname)],
      rank == "family" & class != "" ~ gbif[match(class, fullname)],
      rank == "family" & phylum != "" ~ gbif[match(phylum, fullname)],
      rank == "family" ~ gbif[match(kingdom, fullname)],
      # genus
      rank == "genus" & family != "" ~ gbif[match(family, fullname)],
      rank == "genus" & order != "" ~ gbif[match(order, fullname)],
      rank == "genus" & class != "" ~ gbif[match(class, fullname)],
      rank == "genus" & phylum != "" ~ gbif[match(phylum, fullname)],
      rank == "genus" ~ gbif[match(kingdom, fullname)],
      # species
      rank == "species" ~ gbif[match(genus, fullname)],
      # subspecies
      rank == "subspecies" ~ gbif[match(paste(genus, species), fullname)],
      TRUE ~ NA_character_))

# check again
taxonomy %>%
  filter((!gbif_parent %in% gbif) | (!lpsn_parent %in% lpsn) | (!mycobank_parent %in% mycobank)) %>%
  count(source = case_when(!gbif_parent %in% gbif ~ "GBIF",
                           !lpsn_parent %in% lpsn ~ "LPSN",
                           !mycobank_parent %in% mycobank ~ "MycoBank",
                           TRUE ~ "?"),
        rank)


# Save intermediate results (3) -------------------------------------------------------------------

saveRDS(taxonomy, "data-raw/taxonomy3.rds")
# taxonomy <- readRDS("data-raw/taxonomy3.rds")


# Finalised taxonomy (without additional data) ----------------------------------------------------

message(
  "\nCongratulations! The new taxonomic table will contain ", format(nrow(taxonomy), big.mark = " "), " rows.\n",
  "This was ", format(nrow(AMR::microorganisms), big.mark = " "), " rows.\n"
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


# Add SNOMED CT -----------------------------------------------------------------------------------

# we will use Public Health Information Network Vocabulary Access and Distribution System (PHIN VADS)
# as a source, which copies directly from the latest US SNOMED CT version
# - go to https://phinvads.cdc.gov/vads/ViewValueSet.action?oid=2.16.840.1.114222.4.11.1009
# - check that current online version is higher than TAXONOMY_VERSION$SNOMED
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
  ))
snomed <- snomed %>%
  filter(fullname %in% taxonomy$fullname)

message(nrow(snomed), " SNOMED codes will be added to ", n_distinct(snomed$fullname), " microorganisms")

snomed <- snomed %>%
  group_by(fullname) %>%
  summarise(snomed = list(snomed))

taxonomy <- taxonomy %>%
  left_join(snomed, by = "fullname")


# Add oxygen tolerance (aerobe/anaerobe) ----------------------------------------------------------

# We will use the BacDive data base for this:
# - go to https://bacdive.dsmz.de/advsearch a
# - filter 'Oxygen tolerance' on "*" and click Submit
# - click on the 'Download tabel as CSV' button
bacdive <- vroom::vroom("data-raw/bacdive.csv", skip = 2) %>% 
  select(species, oxygen = `Oxygen tolerance`)
bacdive <- bacdive %>% 
  # fill in missing species from previous rows
  mutate(fullname = if_else(is.na(species), lag(species), species)) %>%
  filter(!is.na(species), !is.na(oxygen), oxygen %unlike% "tolerant", species %unlike% "unclassified") %>% 
  select(-species)
bacdive <- bacdive %>% 
  # now determine type per species
  group_by(fullname) %>%
  summarise(fullname = first(fullname),
            oxygen_tolerance = case_when(any(oxygen %like% "facultative") ~ "facultative anaerobe",
                                         all(oxygen == "microaerophile") ~ "microaerophile",
                                         all(oxygen %in% c("anaerobe", "obligate anaerobe")) ~ "anaerobe",
                                         all(oxygen %in% c("anaerobe", "obligate anaerobe", "microaerophile")) ~ "anaerobe/microaerophile",
                                         all(oxygen %in% c("aerobe", "obligate aerobe")) ~ "aerobe",
                                         all(!oxygen %in% c("anaerobe", "obligate anaerobe")) ~ "aerobe",
                                         all(c("aerobe", "anaerobe") %in% oxygen) ~ "facultative anaerobe",
                                         TRUE ~ NA_character_))
# now find all synonyms and copy them from their current taxonomic names
synonyms <- taxonomy %>% 
  filter(status == "synonym") %>% 
  transmute(mo,
            fullname_old = fullname,
            current = synonym_mo_to_accepted_mo(mo, fill_in_accepted = FALSE, dataset = taxonomy)) %>% 
  filter(!is.na(current)) %>% 
  mutate(fullname = taxonomy$fullname[match(current, taxonomy$mo)]) %>% 
  left_join(bacdive, by = "fullname") %>% 
  filter(!is.na(oxygen_tolerance)) %>% 
  select(fullname, oxygen_tolerance)

bacdive <- bacdive %>% 
  bind_rows(synonyms) %>% 
  distinct()

bacdive_genus <- bacdive %>%
  mutate(oxygen = oxygen_tolerance,
         genus = taxonomy$genus[match(fullname, taxonomy$fullname)]) %>% 
  group_by(fullname = genus) %>% 
  summarise(oxygen_tolerance = case_when(any(oxygen == "facultative anaerobe") ~ "facultative anaerobe",
                                         any(oxygen == "anaerobe/microaerophile") ~ "anaerobe/microaerophile",
                                         all(oxygen == "microaerophile") ~ "microaerophile",
                                         all(oxygen == "anaerobe") ~ "anaerobe",
                                         all(oxygen == "aerobe") ~ "aerobe",
                                         TRUE ~ "facultative anaerobe"))
bacdive <- bacdive %>% 
  bind_rows(bacdive_genus) %>% 
  arrange(fullname)

bacdive_other <- taxonomy %>%
  filter(kingdom == "Bacteria", rank == "species", !fullname %in% bacdive$fullname, genus %in% bacdive$fullname) %>%
  select(fullname, genus) %>%
  left_join(bacdive, by = c("genus" = "fullname")) %>%
  mutate(oxygen_tolerance = if_else(oxygen_tolerance %in% c("aerobe", "anaerobe", "microaerophile", "anaerobe/microaerophile"),
                                    oxygen_tolerance,
                                    paste("likely", oxygen_tolerance))) %>% 
  select(fullname, oxygen_tolerance) %>% 
  distinct(fullname, .keep_all = TRUE)

bacdive <- bacdive %>% 
  bind_rows(bacdive_other) %>% 
  arrange(fullname) %>% 
  distinct(fullname, .keep_all = TRUE)

taxonomy <- taxonomy %>%
  left_join(bacdive, by = "fullname") %>% 
  relocate(oxygen_tolerance, .after = ref)


# Restore 'synonym' microorganisms to 'accepted' --------------------------------------------------

# If there are some synonyms that need to be corrected to 'accepted', you can do that here.
# Before 2024, we encountered this (but currently, this is all good, no action needed):

# according to LPSN: Stenotrophomonas maltophilia is the correct name if this species is regarded as a separate species (i.e., if its nomenclatural type is not assigned to another species whose name is validly published, legitimate and not rejected and has priority) within a separate genus Stenotrophomonas.
# https://lpsn.dsmz.de/species/stenotrophomonas-maltophilia

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
old_groups <- microorganisms %>%
  filter(rank == "species group") %>% 
  mutate(mo = as.character(mo))
# fix, remove later
substr(old_groups$mo[!substr(old_groups$mo, 1, 1) %in% c("B", "U") & old_groups$kingdom == "Bacteria" & !is.na(old_groups$kingdom)], 1, 1) <- "B"
substr(old_groups$mo[!substr(old_groups$mo, 1, 1) %in% c("F", "U") & old_groups$kingdom == "Fungi" & !is.na(old_groups$kingdom)], 1, 1) <- "F"
substr(old_groups$mo[!substr(old_groups$mo, 1, 1) %in% c("C", "U") & old_groups$kingdom == "Chromista" & !is.na(old_groups$kingdom)], 1, 1) <- "C"
substr(old_groups$mo[!substr(old_groups$mo, 1, 1) %in% c("A", "U") & old_groups$kingdom == "Archaea" & !is.na(old_groups$kingdom)], 1, 1) <- "A"
substr(old_groups$mo[!substr(old_groups$mo, 1, 1) %in% c("P", "U") & old_groups$kingdom == "Protozoa" & !is.na(old_groups$kingdom)], 1, 1) <- "P"
substr(old_groups$mo[!substr(old_groups$mo, 1, 2) %in% c("AN", "UN") & old_groups$kingdom == "Animalia" & !is.na(old_groups$kingdom)], 1, 2) <- "AN"
substr(old_groups$mo[!substr(old_groups$mo, 1, 2) %in% c("PL", "UN") & old_groups$kingdom == "Plantae" & !is.na(old_groups$kingdom)], 1, 2) <- "PL"

# use current taxonomy
groups <- taxonomy[match(old_groups$genus, taxonomy$genus), ]
groups <- groups %>% 
  mutate(mo = as.character(old_groups$mo),
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
         gbif_renamed_to = NA_character_) %>% 
  select(-c(lpsn, lpsn_renamed_to, mycobank, mycobank_renamed_to, gbif, gbif_renamed_to, snomed))

taxonomy <- taxonomy %>%
  filter(!mo %in% groups$mo) %>% 
  bind_rows(groups)

# we added an MO code, so make sure everything is still unique
any(duplicated(taxonomy$mo))
taxonomy$mo[duplicated(taxonomy$mo)]
any(duplicated(taxonomy$fullname))
taxonomy$fullname[duplicated(taxonomy$fullname)]


# Set unknown ranks -------------------------------------------------------------------------------

taxonomy$rank[which(taxonomy$fullname %like% "unknown")] <- "(unknown rank)"


# Some final checks -------------------------------------------------------------------------------

fix_old_mos <- function(dataset) {
  df <- dataset %>% mutate(mo = as.character(mo))
  df$mo[which(!df$mo %in% taxonomy$mo)] <- df %>% filter(!mo %in% taxonomy$mo) %>% mutate(name = mo_name(mo, keep_synonyms = TRUE), new_mo = taxonomy$mo[match(name, taxonomy$fullname)]) %>% pull(new_mo)
  class(df$mo) <- c("mo", "character")
  df
}

# and check: these codes should not be missing (will otherwise throw a unit test error):
AMR::microorganisms.codes %>% filter(!mo %in% taxonomy$mo)
  # fix:
  microorganisms.codes <- fix_old_mos(AMR::microorganisms.codes)
  usethis::use_data(microorganisms.codes, overwrite = TRUE, version = 2, compress = "xz")
  rm(microorganisms.codes)
  
AMR::clinical_breakpoints %>% filter(!mo %in% taxonomy$mo)

AMR::example_isolates %>% filter(!mo %in% taxonomy$mo)

AMR::intrinsic_resistant %>% filter(!mo %in% taxonomy$mo)
  # fix:
  intrinsic_resistant <- fix_old_mos(AMR::intrinsic_resistant)
  usethis::use_data(intrinsic_resistant, overwrite = TRUE, version = 2, compress = "xz")
  rm(intrinsic_resistant)


# all our previously manually added names should be in it
all(microorganisms$fullname[microorganisms$source == "manually added"] %in% taxonomy$fullname)
microorganisms$fullname[!microorganisms$fullname[microorganisms$source == "manually added"] %in% taxonomy$fullname]

# we added an MO code, so make sure everything is still unique
any(duplicated(taxonomy$mo))
any(duplicated(taxonomy$fullname))


# Save to package ---------------------------------------------------------------------------------

# format to tibble and remove non-ASCII characters

saveRDS(taxonomy, "data-raw/taxonomy3b.rds")
# taxonomy <- readRDS("data-raw/taxonomy3b.rds")
taxonomy <- taxonomy %>%
  arrange(fullname) %>% 
  select(mo, fullname, status, kingdom:subspecies, rank, ref, oxygen_tolerance, source, starts_with("lpsn"), starts_with("mycobank"), starts_with("gbif"), prevalence, snomed) %>%
  df_remove_nonASCII()

microorganisms <- taxonomy

# set class <mo>
class(microorganisms$mo) <- c("mo", "character")
microorganisms <- microorganisms %>% arrange(fullname) %>% df_remove_nonASCII()
usethis::use_data(microorganisms, overwrite = TRUE, version = 2, compress = "xz")
rm(microorganisms)

# DON'T FORGET TO UPDATE R/_globals.R!

# load new data set
devtools::load_all(".")


# Reset previously changed mo codes ---------------------------------------------------------------

if (!identical(clinical_breakpoints$mo, as.mo(clinical_breakpoints$mo, language = NULL, keep_synonyms = FALSE))) {
  clinical_breakpoints$mo <- as.mo(clinical_breakpoints$mo, language = NULL, keep_synonyms = FALSE)
  usethis::use_data(clinical_breakpoints, overwrite = TRUE, version = 2, compress = "xz")
  rm(clinical_breakpoints)
}

if (!identical(microorganisms.codes$mo, as.mo(microorganisms.codes$mo, language = NULL, keep_synonyms = FALSE))) {
  microorganisms.codes <- microorganisms.codes %>% filter(mo %in% microorganisms$mo)
  microorganisms.codes$mo <- as.mo(microorganisms.codes$mo, language = NULL, keep_synonyms = FALSE)
  usethis::use_data(microorganisms.codes, overwrite = TRUE, version = 2, compress = "xz")
  rm(microorganisms.codes)
}

if (!identical(microorganisms.groups$mo_group, as.mo(microorganisms.groups$mo_group, language = NULL, keep_synonyms = FALSE)) ||
    !identical(microorganisms.groups$mo, as.mo(microorganisms.groups$mo, language = NULL, keep_synonyms = FALSE))) {
  microorganisms.groups$mo_group <- as.mo(microorganisms.groups$mo_group_name, language = NULL, keep_synonyms = FALSE)
  microorganisms.groups$mo <- as.mo(microorganisms.groups$mo_name, language = NULL, keep_synonyms = FALSE)
  usethis::use_data(microorganisms.groups, overwrite = TRUE, version = 2)
  rm(microorganisms.groups)
}

if (!identical(example_isolates$mo, as.mo(example_isolates$mo, language = NULL, keep_synonyms = FALSE))) {
  example_isolates$mo <- as.mo(example_isolates$mo, language = NULL, keep_synonyms = FALSE)
  usethis::use_data(example_isolates, overwrite = TRUE, version = 2)
  rm(example_isolates)
}

if (!identical(intrinsic_resistant$mo, as.mo(intrinsic_resistant$mo, language = NULL))) {
  intrinsic_resistant$mo <- as.mo(intrinsic_resistant$mo, language = NULL, keep_synonyms = FALSE)
  intrinsic_resistant <- distinct(intrinsic_resistant)
  usethis::use_data(intrinsic_resistant, overwrite = TRUE, version = 2)
  rm(intrinsic_resistant)
}


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
