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

#' Get Properties of a Microorganism
#'
#' Use these functions to return a specific property of a microorganism based on the latest accepted taxonomy. All input values will be evaluated internally with [as.mo()], which makes it possible to use microbial abbreviations, codes and names as input. See *Examples*.
#' @param x any [character] (vector) that can be coerced to a valid microorganism code with [as.mo()]. Can be left blank for auto-guessing the column containing microorganism codes if used in a data set, see *Examples*.
#' @param property one of the column names of the [microorganisms] data set: `r vector_or(colnames(microorganisms), sort = FALSE, quotes = TRUE)`, or must be `"shortname"`
#' @inheritParams as.mo
#' @param ... other arguments passed on to [as.mo()], such as 'minimum_matching_score', 'ignore_pattern', and 'remove_from_input'
#' @param ab any (vector of) text that can be coerced to a valid antibiotic code with [as.ab()]
#' @param open browse the URL using [`browseURL()`][utils::browseURL()]
#' @details All functions will, at default, keep old taxonomic properties. Please refer to this example, knowing that *Escherichia blattae* was renamed to *Shimwellia blattae* in 2010:
#' - `mo_name("Escherichia blattae")` will return `"Shimwellia blattae"` (with a message about the renaming)
#' - `mo_ref("Escherichia blattae", keep_synonyms = TRUE)` will return `"Burgess et al., 1973"` (with a warning about the renaming)
#' - `mo_ref("Shimwellia blattae", keep_synonyms = FALSE)` will return `"Priest et al., 2010"` (without a message)
#'
#' The short name - [mo_shortname()] - almost always returns the first character of the genus and the full species, like `"E. coli"`. Exceptions are abbreviations of staphylococci (such as *"CoNS"*, Coagulase-Negative Staphylococci) and beta-haemolytic streptococci (such as *"GBS"*, Group B Streptococci). Please bear in mind that e.g. *E. coli* could mean *Escherichia coli* (kingdom of Bacteria) as well as *Entamoeba coli* (kingdom of Protozoa). Returning to the full name will be done using [as.mo()] internally, giving priority to bacteria and human pathogens, i.e. `"E. coli"` will be considered *Escherichia coli*. In other words, `mo_fullname(mo_shortname("Entamoeba coli"))` returns `"Escherichia coli"`.
#'
#' Since the top-level of the taxonomy is sometimes referred to as 'kingdom' and sometimes as 'domain', the functions [mo_kingdom()] and [mo_domain()] return the exact same results.
#'
#' The Gram stain - [mo_gramstain()] - will be determined based on the taxonomic kingdom and phylum. According to Cavalier-Smith (2002, [PMID 11837318](https://pubmed.ncbi.nlm.nih.gov/11837318)), who defined subkingdoms Negibacteria and Posibacteria, only these phyla are Posibacteria: Actinobacteria, Chloroflexi, Firmicutes and Tenericutes. These bacteria are considered Gram-positive, except for members of the class Negativicutes which are Gram-negative. Members of other bacterial phyla are all considered Gram-negative. Species outside the kingdom of Bacteria will return a value `NA`. Functions [mo_is_gram_negative()] and [mo_is_gram_positive()] always return `TRUE` or `FALSE` (except when the input is `NA` or the MO code is `UNKNOWN`), thus always return `FALSE` for species outside the taxonomic kingdom of Bacteria.
#'
#' Determination of yeasts - [mo_is_yeast()] - will be based on the taxonomic kingdom and class. *Budding yeasts* are fungi of the phylum Ascomycetes, class Saccharomycetes (also called Hemiascomycetes). *True yeasts* are aggregated into the underlying order Saccharomycetales. Thus, for all microorganisms that are fungi and member of the taxonomic class Saccharomycetes, the function will return `TRUE`. It returns `FALSE` otherwise (except when the input is `NA` or the MO code is `UNKNOWN`).
#'
#' Intrinsic resistance - [mo_is_intrinsic_resistant()] - will be determined based on the [intrinsic_resistant] data set, which is based on `r format_eucast_version_nr(3.3)`. The [mo_is_intrinsic_resistant()] functions can be vectorised over arguments `x` (input for microorganisms) and over `ab` (input for antibiotics).
#'
#' All output [will be translated][translate] where possible.
#'
#' The function [mo_url()] will return the direct URL to the online database entry, which also shows the scientific reference of the concerned species.
#'
#' SNOMED codes - [mo_snomed()] - are from the version of `r documentation_date(TAXONOMY_VERSION$SNOMED$accessed_date)`. See *Source* and the [microorganisms] data set for more info.
#' @inheritSection mo_matching_score Matching Score for Microorganisms
#' @inheritSection as.mo Source
#' @rdname mo_property
#' @name mo_property
#' @return
#' - An [integer] in case of [mo_year()]
#' - A [list] in case of [mo_taxonomy()] and [mo_info()]
#' - A named [character] in case of [mo_url()]
#' - A [numeric] in case of [mo_snomed()]
#' - A [character] in all other cases
#' @export
#' @seealso Data set [microorganisms]
#' @inheritSection AMR Reference Data Publicly Available
#' @examples
#' # taxonomic tree -----------------------------------------------------------
#' mo_kingdom("Klebsiella pneumoniae")
#' mo_phylum("Klebsiella pneumoniae")
#' mo_class("Klebsiella pneumoniae")
#' mo_order("Klebsiella pneumoniae")
#' mo_family("Klebsiella pneumoniae")
#' mo_genus("Klebsiella pneumoniae")
#' mo_species("Klebsiella pneumoniae")
#' mo_subspecies("Klebsiella pneumoniae")
#'
#' # colloquial properties ----------------------------------------------------
#' mo_name("Klebsiella pneumoniae")
#' mo_fullname("Klebsiella pneumoniae")
#' mo_shortname("Klebsiella pneumoniae")
#'
#' # other properties ---------------------------------------------------------
#' mo_gramstain("Klebsiella pneumoniae")
#' mo_snomed("Klebsiella pneumoniae")
#' mo_type("Klebsiella pneumoniae")
#' mo_rank("Klebsiella pneumoniae")
#' mo_url("Klebsiella pneumoniae")
#' mo_synonyms("Klebsiella pneumoniae")
#'
#' # scientific reference -----------------------------------------------------
#' mo_ref("Klebsiella pneumoniae")
#' mo_authors("Klebsiella pneumoniae")
#' mo_year("Klebsiella pneumoniae")
#' mo_lpsn("Klebsiella pneumoniae")
#' mo_gbif("Klebsiella pneumoniae")
#'
#' # abbreviations known in the field -----------------------------------------
#' mo_genus("MRSA")
#' mo_species("MRSA")
#' mo_shortname("VISA")
#' mo_gramstain("VISA")
#'
#' mo_genus("EHEC")
#' mo_species("EHEC")
#'
#' # known subspecies ---------------------------------------------------------
#' mo_fullname("K. pneu rh")
#' mo_shortname("K. pneu rh")
#'
#' \donttest{
#' # Becker classification, see ?as.mo ----------------------------------------
#' mo_fullname("Staph. epi")
#' mo_fullname("Staph. epi", Becker = TRUE)
#' mo_shortname("Staph. epi")
#' mo_shortname("Staph. epi", Becker = TRUE)
#'
#' # Lancefield classification, see ?as.mo ------------------------------------
#' mo_fullname("S. pyo")
#' mo_fullname("S. pyo", Lancefield = TRUE)
#' mo_shortname("S. pyo")
#' mo_shortname("S. pyo", Lancefield = TRUE)
#'
#'
#' # language support  --------------------------------------------------------
#' mo_gramstain("Klebsiella pneumoniae", language = "de") # German
#' mo_gramstain("Klebsiella pneumoniae", language = "nl") # Dutch
#' mo_gramstain("Klebsiella pneumoniae", language = "es") # Spanish
#' mo_gramstain("Klebsiella pneumoniae", language = "el") # Greek
#' mo_gramstain("Klebsiella pneumoniae", language = "uk") # Ukrainian
#'
#' # mo_type is equal to mo_kingdom, but mo_kingdom will remain official
#' mo_kingdom("Klebsiella pneumoniae")
#' mo_type("Klebsiella pneumoniae")
#' mo_type("Klebsiella pneumoniae", language = "nl")
#'
#' mo_fullname("S. pyogenes", Lancefield = TRUE, language = "de")
#' mo_fullname("S. pyogenes", Lancefield = TRUE, language = "uk")
#'
#'
#' # other --------------------------------------------------------------------
#'
#' mo_is_yeast(c("Candida", "Trichophyton", "Klebsiella"))
#'
#' # gram stains and intrinsic resistance can be used as a filter in dplyr verbs
#' if (require("dplyr")) {
#'   example_isolates %>%
#'     filter(mo_is_gram_positive()) %>% 
#'     count(mo_genus(), count = TRUE)
#' }
#' if (require("dplyr")) {
#'   example_isolates %>%
#'     filter(mo_is_intrinsic_resistant(ab = "vanco")) %>% 
#'     count(mo_genus(), count = TRUE)
#' }
#'
#'
#' # get a list with the complete taxonomy (from kingdom to subspecies)
#' mo_taxonomy("Klebsiella pneumoniae")
#'
#' # get a list with the taxonomy, the authors, Gram-stain,
#' # SNOMED codes, and URL to the online database
#' mo_info("Klebsiella pneumoniae")
#' }
mo_name <- function(x, language = get_AMR_locale(), keep_synonyms = getOption("AMR_keep_synonyms", FALSE), ...) {
  if (missing(x)) {
    # this tries to find the data and an <mo> column
    x <- find_mo_col(fn = "mo_name")
  }
  meet_criteria(x, allow_NA = TRUE)
  language <- validate_language(language)
  meet_criteria(keep_synonyms, allow_class = "logical", has_length = 1)

  translate_into_language(mo_validate(x = x, property = "fullname", language = language, keep_synonyms = keep_synonyms, ...),
    language = language,
    only_unknown = FALSE,
    only_affect_mo_names = TRUE
  )
}

#' @rdname mo_property
#' @export
mo_fullname <- mo_name

#' @rdname mo_property
#' @export
mo_shortname <- function(x, language = get_AMR_locale(), keep_synonyms = getOption("AMR_keep_synonyms", FALSE), ...) {
  if (missing(x)) {
    # this tries to find the data and an <mo> column
    x <- find_mo_col(fn = "mo_shortname")
  }
  meet_criteria(x, allow_NA = TRUE)
  language <- validate_language(language)
  meet_criteria(keep_synonyms, allow_class = "logical", has_length = 1)

  x.mo <- as.mo(x, language = language, keep_synonyms = keep_synonyms, ...)

  metadata <- get_mo_uncertainties()

  replace_empty <- function(x) {
    x[x == ""] <- "spp."
    x
  }

  # get first char of genus and complete species in English
  genera <- mo_genus(x.mo, language = NULL, keep_synonyms = keep_synonyms)
  shortnames <- paste0(substr(genera, 1, 1), ". ", replace_empty(mo_species(x.mo, language = NULL, keep_synonyms = keep_synonyms)))

  # exceptions for where no species is known
  shortnames[shortnames %like% ".[.] spp[.]"] <- genera[shortnames %like% ".[.] spp[.]"]
  # exceptions for staphylococci
  shortnames[shortnames == "S. coagulase-negative"] <- "CoNS"
  shortnames[shortnames == "S. coagulase-positive"] <- "CoPS"
  # exceptions for streptococci: Group A Streptococcus -> GAS
  shortnames[shortnames %like% "S. group [ABCDFGHK]"] <- paste0("G", gsub("S. group ([ABCDFGHK])", "\\1", shortnames[shortnames %like% "S. group [ABCDFGHK]"], perl = TRUE), "S")
  # unknown species etc.
  shortnames[shortnames %like% "unknown"] <- paste0("(", trimws2(gsub("[^a-zA-Z -]", "", shortnames[shortnames %like% "unknown"], perl = TRUE)), ")")

  shortnames[is.na(x.mo)] <- NA_character_
  load_mo_uncertainties(metadata)
  translate_into_language(shortnames, language = language, only_unknown = FALSE, only_affect_mo_names = TRUE)
}



#' @rdname mo_property
#' @export
mo_subspecies <- function(x, language = get_AMR_locale(), keep_synonyms = getOption("AMR_keep_synonyms", FALSE), ...) {
  if (missing(x)) {
    # this tries to find the data and an <mo> column
    x <- find_mo_col(fn = "mo_subspecies")
  }
  meet_criteria(x, allow_NA = TRUE)
  language <- validate_language(language)
  meet_criteria(keep_synonyms, allow_class = "logical", has_length = 1)

  translate_into_language(mo_validate(x = x, property = "subspecies", language = language, keep_synonyms = keep_synonyms, ...), language = language, only_unknown = TRUE)
}

#' @rdname mo_property
#' @export
mo_species <- function(x, language = get_AMR_locale(), keep_synonyms = getOption("AMR_keep_synonyms", FALSE), ...) {
  if (missing(x)) {
    # this tries to find the data and an <mo> column
    x <- find_mo_col(fn = "mo_species")
  }
  meet_criteria(x, allow_NA = TRUE)
  language <- validate_language(language)
  meet_criteria(keep_synonyms, allow_class = "logical", has_length = 1)

  translate_into_language(mo_validate(x = x, property = "species", language = language, keep_synonyms = keep_synonyms, ...), language = language, only_unknown = TRUE)
}

#' @rdname mo_property
#' @export
mo_genus <- function(x, language = get_AMR_locale(), keep_synonyms = getOption("AMR_keep_synonyms", FALSE), ...) {
  if (missing(x)) {
    # this tries to find the data and an <mo> column
    x <- find_mo_col(fn = "mo_genus")
  }
  meet_criteria(x, allow_NA = TRUE)
  language <- validate_language(language)
  meet_criteria(keep_synonyms, allow_class = "logical", has_length = 1)

  translate_into_language(mo_validate(x = x, property = "genus", language = language, keep_synonyms = keep_synonyms, ...), language = language, only_unknown = TRUE)
}

#' @rdname mo_property
#' @export
mo_family <- function(x, language = get_AMR_locale(), keep_synonyms = getOption("AMR_keep_synonyms", FALSE), ...) {
  if (missing(x)) {
    # this tries to find the data and an <mo> column
    x <- find_mo_col(fn = "mo_family")
  }
  meet_criteria(x, allow_NA = TRUE)
  language <- validate_language(language)
  meet_criteria(keep_synonyms, allow_class = "logical", has_length = 1)

  translate_into_language(mo_validate(x = x, property = "family", language = language, keep_synonyms = keep_synonyms, ...), language = language, only_unknown = TRUE)
}

#' @rdname mo_property
#' @export
mo_order <- function(x, language = get_AMR_locale(), keep_synonyms = getOption("AMR_keep_synonyms", FALSE), ...) {
  if (missing(x)) {
    # this tries to find the data and an <mo> column
    x <- find_mo_col(fn = "mo_order")
  }
  meet_criteria(x, allow_NA = TRUE)
  language <- validate_language(language)
  meet_criteria(keep_synonyms, allow_class = "logical", has_length = 1)

  translate_into_language(mo_validate(x = x, property = "order", language = language, keep_synonyms = keep_synonyms, ...), language = language, only_unknown = TRUE)
}

#' @rdname mo_property
#' @export
mo_class <- function(x, language = get_AMR_locale(), keep_synonyms = getOption("AMR_keep_synonyms", FALSE), ...) {
  if (missing(x)) {
    # this tries to find the data and an <mo> column
    x <- find_mo_col(fn = "mo_class")
  }
  meet_criteria(x, allow_NA = TRUE)
  language <- validate_language(language)
  meet_criteria(keep_synonyms, allow_class = "logical", has_length = 1)

  translate_into_language(mo_validate(x = x, property = "class", language = language, keep_synonyms = keep_synonyms, ...), language = language, only_unknown = TRUE)
}

#' @rdname mo_property
#' @export
mo_phylum <- function(x, language = get_AMR_locale(), keep_synonyms = getOption("AMR_keep_synonyms", FALSE), ...) {
  if (missing(x)) {
    # this tries to find the data and an <mo> column
    x <- find_mo_col(fn = "mo_phylum")
  }
  meet_criteria(x, allow_NA = TRUE)
  language <- validate_language(language)
  meet_criteria(keep_synonyms, allow_class = "logical", has_length = 1)

  translate_into_language(mo_validate(x = x, property = "phylum", language = language, keep_synonyms = keep_synonyms, ...), language = language, only_unknown = TRUE)
}

#' @rdname mo_property
#' @export
mo_kingdom <- function(x, language = get_AMR_locale(), keep_synonyms = getOption("AMR_keep_synonyms", FALSE), ...) {
  if (missing(x)) {
    # this tries to find the data and an <mo> column
    x <- find_mo_col(fn = "mo_kingdom")
  }
  meet_criteria(x, allow_NA = TRUE)
  language <- validate_language(language)
  meet_criteria(keep_synonyms, allow_class = "logical", has_length = 1)

  translate_into_language(mo_validate(x = x, property = "kingdom", language = language, keep_synonyms = keep_synonyms, ...), language = language, only_unknown = TRUE)
}

#' @rdname mo_property
#' @export
mo_domain <- mo_kingdom

#' @rdname mo_property
#' @export
mo_type <- function(x, language = get_AMR_locale(), keep_synonyms = getOption("AMR_keep_synonyms", FALSE), ...) {
  if (missing(x)) {
    # this tries to find the data and an <mo> column
    x <- find_mo_col(fn = "mo_type")
  }
  meet_criteria(x, allow_NA = TRUE)
  language <- validate_language(language)
  meet_criteria(keep_synonyms, allow_class = "logical", has_length = 1)

  x.mo <- as.mo(x, language = language, keep_synonyms = keep_synonyms, ...)
  out <- mo_kingdom(x.mo, language = NULL, keep_synonyms = keep_synonyms)
  out[which(mo_is_yeast(x.mo, keep_synonyms = keep_synonyms))] <- "Yeasts"
  translate_into_language(out, language = language, only_unknown = FALSE)
}

#' @rdname mo_property
#' @export
mo_status <- function(x, language = get_AMR_locale(), keep_synonyms = getOption("AMR_keep_synonyms", FALSE), ...) {
  if (missing(x)) {
    # this tries to find the data and an <mo> column
    x <- find_mo_col(fn = "mo_status")
  }
  meet_criteria(x, allow_NA = TRUE)
  language <- validate_language(language)
  meet_criteria(keep_synonyms, allow_class = "logical", has_length = 1)

  translate_into_language(mo_validate(x = x, property = "status", language = language, keep_synonyms = keep_synonyms, ...), language = language, only_unknown = TRUE)
}

#' @rdname mo_property
#' @export
mo_gramstain <- function(x, language = get_AMR_locale(), keep_synonyms = getOption("AMR_keep_synonyms", FALSE), ...) {
  if (missing(x)) {
    # this tries to find the data and an <mo> column
    x <- find_mo_col(fn = "mo_gramstain")
  }
  meet_criteria(x, allow_NA = TRUE)
  language <- validate_language(language)
  meet_criteria(keep_synonyms, allow_class = "logical", has_length = 1)

  x.mo <- as.mo(x, language = language, keep_synonyms = keep_synonyms, ...)
  metadata <- get_mo_uncertainties()

  x <- rep(NA_character_, length(x))
  # make all bacteria Gram negative
  x[mo_kingdom(x.mo, language = NULL, keep_synonyms = keep_synonyms) == "Bacteria"] <- "Gram-negative"
  # overwrite these 4 phyla with Gram-positives
  # Source: https://itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=956097 (Cavalier-Smith, 2002)
  x[(mo_phylum(x.mo, language = NULL, keep_synonyms = keep_synonyms) %in% c(
    "Actinobacteria",
    "Chloroflexi",
    "Firmicutes",
    "Tenericutes",
    "Bacillota" # this one is new! It was renamed from Firmicutes by Gibbons et al., 2021
  ) &
    # but class Negativicutes (of phylum Firmicutes) are Gram-negative!
    mo_class(x.mo, language = NULL, keep_synonyms = keep_synonyms) != "Negativicutes")
  # and of course our own ID for Gram-positives
  | x.mo == "B_GRAMP"] <- "Gram-positive"

  load_mo_uncertainties(metadata)
  translate_into_language(x, language = language, only_unknown = FALSE)
}

#' @rdname mo_property
#' @export
mo_is_gram_negative <- function(x, language = get_AMR_locale(), keep_synonyms = getOption("AMR_keep_synonyms", FALSE), ...) {
  if (missing(x)) {
    # this tries to find the data and an <mo> column
    x <- find_mo_col(fn = "mo_is_gram_negative")
  }
  meet_criteria(x, allow_NA = TRUE)
  language <- validate_language(language)
  meet_criteria(keep_synonyms, allow_class = "logical", has_length = 1)

  x.mo <- as.mo(x, language = language, keep_synonyms = keep_synonyms, ...)
  metadata <- get_mo_uncertainties()
  grams <- mo_gramstain(x.mo, language = NULL, keep_synonyms = keep_synonyms)
  load_mo_uncertainties(metadata)
  out <- grams == "Gram-negative" & !is.na(grams)
  out[x.mo %in% c(NA_character_, "UNKNOWN")] <- NA
  out
}

#' @rdname mo_property
#' @export
mo_is_gram_positive <- function(x, language = get_AMR_locale(), keep_synonyms = getOption("AMR_keep_synonyms", FALSE), ...) {
  if (missing(x)) {
    # this tries to find the data and an <mo> column
    x <- find_mo_col(fn = "mo_is_gram_positive")
  }
  meet_criteria(x, allow_NA = TRUE)
  language <- validate_language(language)
  meet_criteria(keep_synonyms, allow_class = "logical", has_length = 1)

  x.mo <- as.mo(x, language = language, keep_synonyms = keep_synonyms, ...)
  metadata <- get_mo_uncertainties()
  grams <- mo_gramstain(x.mo, language = NULL, keep_synonyms = keep_synonyms)
  load_mo_uncertainties(metadata)
  out <- grams == "Gram-positive" & !is.na(grams)
  out[x.mo %in% c(NA_character_, "UNKNOWN")] <- NA
  out
}

#' @rdname mo_property
#' @export
mo_is_yeast <- function(x, language = get_AMR_locale(), keep_synonyms = getOption("AMR_keep_synonyms", FALSE), ...) {
  if (missing(x)) {
    # this tries to find the data and an <mo> column
    x <- find_mo_col(fn = "mo_is_yeast")
  }
  meet_criteria(x, allow_NA = TRUE)
  language <- validate_language(language)
  meet_criteria(keep_synonyms, allow_class = "logical", has_length = 1)

  x.mo <- as.mo(x, language = language, keep_synonyms = keep_synonyms, ...)
  metadata <- get_mo_uncertainties()

  x.kingdom <- mo_kingdom(x.mo, language = NULL, keep_synonyms = keep_synonyms)
  x.class <- mo_class(x.mo, language = NULL, keep_synonyms = keep_synonyms)

  load_mo_uncertainties(metadata)

  out <- rep(FALSE, length(x))
  out[x.kingdom == "Fungi" & x.class == "Saccharomycetes"] <- TRUE
  out[x.mo %in% c(NA_character_, "UNKNOWN")] <- NA
  out
}

#' @rdname mo_property
#' @export
mo_is_intrinsic_resistant <- function(x, ab, language = get_AMR_locale(), keep_synonyms = getOption("AMR_keep_synonyms", FALSE), ...) {
  if (missing(x)) {
    # this tries to find the data and an <mo> column
    x <- find_mo_col(fn = "mo_is_intrinsic_resistant")
  }
  meet_criteria(x, allow_NA = TRUE)
  meet_criteria(ab, allow_NA = FALSE)
  language <- validate_language(language)
  meet_criteria(keep_synonyms, allow_class = "logical", has_length = 1)

  x <- as.mo(x, language = language, keep_synonyms = keep_synonyms, ...)
  ab <- as.ab(ab, language = NULL, flag_multiple_results = FALSE, info = FALSE)

  if (length(x) == 1 & length(ab) > 1) {
    x <- rep(x, length(ab))
  } else if (length(ab) == 1 & length(x) > 1) {
    ab <- rep(ab, length(x))
  }
  if (length(x) != length(ab)) {
    stop_("length of `x` and `ab` must be equal, or one of them must be of length 1.")
  }

  # show used version number once per session (AMR_env will reload every session)
  if (message_not_thrown_before("mo_is_intrinsic_resistant", "version.mo", entire_session = TRUE)) {
    message_(
      "Determining intrinsic resistance based on ",
      format_eucast_version_nr(3.3, markdown = FALSE), ". ",
      font_red("This note will be shown once per session.")
    )
  }

  # runs against internal vector: INTRINSIC_R (see zzz.R)
  paste(x, ab) %in% INTRINSIC_R
}

#' @rdname mo_property
#' @export
mo_snomed <- function(x, language = get_AMR_locale(), keep_synonyms = getOption("AMR_keep_synonyms", FALSE), ...) {
  if (missing(x)) {
    # this tries to find the data and an <mo> column
    x <- find_mo_col(fn = "mo_snomed")
  }
  meet_criteria(x, allow_NA = TRUE)
  language <- validate_language(language)
  meet_criteria(keep_synonyms, allow_class = "logical", has_length = 1)

  mo_validate(x = x, property = "snomed", language = language, keep_synonyms = keep_synonyms, ...)
}

#' @rdname mo_property
#' @export
mo_ref <- function(x, language = get_AMR_locale(), keep_synonyms = getOption("AMR_keep_synonyms", FALSE), ...) {
  if (missing(x)) {
    # this tries to find the data and an <mo> column
    x <- find_mo_col(fn = "mo_ref")
  }
  meet_criteria(x, allow_NA = TRUE)
  language <- validate_language(language)
  meet_criteria(keep_synonyms, allow_class = "logical", has_length = 1)

  mo_validate(x = x, property = "ref", language = language, keep_synonyms = keep_synonyms, ...)
}

#' @rdname mo_property
#' @export
mo_authors <- function(x, language = get_AMR_locale(), keep_synonyms = getOption("AMR_keep_synonyms", FALSE), ...) {
  if (missing(x)) {
    # this tries to find the data and an <mo> column
    x <- find_mo_col(fn = "mo_authors")
  }
  meet_criteria(x, allow_NA = TRUE)
  language <- validate_language(language)
  meet_criteria(keep_synonyms, allow_class = "logical", has_length = 1)

  x <- mo_validate(x = x, property = "ref", language = language, keep_synonyms = keep_synonyms, ...)
  # remove last 4 digits and presumably the comma and space that preceed them
  x[!is.na(x)] <- gsub(",? ?[0-9]{4}", "", x[!is.na(x)], perl = TRUE)
  suppressWarnings(x)
}

#' @rdname mo_property
#' @export
mo_year <- function(x, language = get_AMR_locale(), keep_synonyms = getOption("AMR_keep_synonyms", FALSE), ...) {
  if (missing(x)) {
    # this tries to find the data and an <mo> column
    x <- find_mo_col(fn = "mo_year")
  }
  meet_criteria(x, allow_NA = TRUE)
  language <- validate_language(language)
  meet_criteria(keep_synonyms, allow_class = "logical", has_length = 1)

  x <- mo_validate(x = x, property = "ref", language = language, keep_synonyms = keep_synonyms, ...)
  # get last 4 digits
  x[!is.na(x)] <- gsub(".*([0-9]{4})$", "\\1", x[!is.na(x)], perl = TRUE)
  suppressWarnings(as.integer(x))
}

#' @rdname mo_property
#' @export
mo_lpsn <- function(x, language = get_AMR_locale(), keep_synonyms = getOption("AMR_keep_synonyms", FALSE), ...) {
  if (missing(x)) {
    # this tries to find the data and an <mo> column
    x <- find_mo_col(fn = "mo_lpsn")
  }
  meet_criteria(x, allow_NA = TRUE)
  language <- validate_language(language)
  meet_criteria(keep_synonyms, allow_class = "logical", has_length = 1)

  mo_validate(x = x, property = "lpsn", language = language, keep_synonyms = keep_synonyms, ...)
}

#' @rdname mo_property
#' @export
mo_gbif <- function(x, language = get_AMR_locale(), keep_synonyms = getOption("AMR_keep_synonyms", FALSE), ...) {
  if (missing(x)) {
    # this tries to find the data and an <mo> column
    x <- find_mo_col(fn = "mo_gbif")
  }
  meet_criteria(x, allow_NA = TRUE)
  language <- validate_language(language)
  meet_criteria(keep_synonyms, allow_class = "logical", has_length = 1)

  mo_validate(x = x, property = "gbif", language = language, keep_synonyms = keep_synonyms, ...)
}

#' @rdname mo_property
#' @export
mo_rank <- function(x, language = get_AMR_locale(), keep_synonyms = getOption("AMR_keep_synonyms", FALSE), ...) {
  if (missing(x)) {
    # this tries to find the data and an <mo> column
    x <- find_mo_col(fn = "mo_rank")
  }
  meet_criteria(x, allow_NA = TRUE)
  language <- validate_language(language)
  meet_criteria(keep_synonyms, allow_class = "logical", has_length = 1)

  mo_validate(x = x, property = "rank", language = language, keep_synonyms = keep_synonyms, ...)
}

#' @rdname mo_property
#' @export
mo_taxonomy <- function(x, language = get_AMR_locale(), keep_synonyms = getOption("AMR_keep_synonyms", FALSE), ...) {
  if (missing(x)) {
    # this tries to find the data and an <mo> column
    x <- find_mo_col(fn = "mo_taxonomy")
  }
  meet_criteria(x, allow_NA = TRUE)
  language <- validate_language(language)
  meet_criteria(keep_synonyms, allow_class = "logical", has_length = 1)

  x <- as.mo(x, language = language, keep_synonyms = keep_synonyms, ...)
  metadata <- get_mo_uncertainties()

  out <- list(
    kingdom = mo_kingdom(x, language = language, keep_synonyms = keep_synonyms),
    phylum = mo_phylum(x, language = language, keep_synonyms = keep_synonyms),
    class = mo_class(x, language = language, keep_synonyms = keep_synonyms),
    order = mo_order(x, language = language, keep_synonyms = keep_synonyms),
    family = mo_family(x, language = language, keep_synonyms = keep_synonyms),
    genus = mo_genus(x, language = language, keep_synonyms = keep_synonyms),
    species = mo_species(x, language = language, keep_synonyms = keep_synonyms),
    subspecies = mo_subspecies(x, language = language, keep_synonyms = keep_synonyms)
  )

  load_mo_uncertainties(metadata)
  out
}

#' @rdname mo_property
#' @export
mo_synonyms <- function(x, language = get_AMR_locale(), keep_synonyms = getOption("AMR_keep_synonyms", FALSE), ...) {
  if (missing(x)) {
    # this tries to find the data and an <mo> column
    x <- find_mo_col(fn = "mo_synonyms")
  }
  meet_criteria(x, allow_NA = TRUE)
  language <- validate_language(language)
  meet_criteria(keep_synonyms, allow_class = "logical", has_length = 1)

  x.mo <- as.mo(x, language = language, keep_synonyms = keep_synonyms, ...)
  metadata <- get_mo_uncertainties()

  syns <- lapply(x.mo, function(y) {
    gbif <- AMR::microorganisms$gbif[match(y, AMR::microorganisms$mo)]
    lpsn <- AMR::microorganisms$lpsn[match(y, AMR::microorganisms$mo)]
    out <- AMR::microorganisms[which(AMR::microorganisms$lpsn_renamed_to == lpsn | AMR::microorganisms$gbif_renamed_to == gbif), "fullname", drop = TRUE]
    if (length(out) == 0) {
      NULL
    } else {
      out
    }
  })

  if (length(syns) > 1) {
    names(syns) <- mo_name(x)
    result <- syns
  } else {
    result <- unlist(syns)
  }

  load_mo_uncertainties(metadata)
  result
}

#' @rdname mo_property
#' @export
mo_info <- function(x, language = get_AMR_locale(), keep_synonyms = getOption("AMR_keep_synonyms", FALSE), ...) {
  if (missing(x)) {
    # this tries to find the data and an <mo> column
    x <- find_mo_col(fn = "mo_info")
  }
  meet_criteria(x, allow_NA = TRUE)
  language <- validate_language(language)
  meet_criteria(keep_synonyms, allow_class = "logical", has_length = 1)

  x <- as.mo(x, language = language, keep_synonyms = keep_synonyms, ...)
  metadata <- get_mo_uncertainties()

  info <- lapply(x, function(y) {
    c(
      mo_taxonomy(y, language = language, keep_synonyms = keep_synonyms),
      list(
        status = mo_status(y, language = language, keep_synonyms = keep_synonyms),
        synonyms = mo_synonyms(y, keep_synonyms = keep_synonyms),
        gramstain = mo_gramstain(y, language = language, keep_synonyms = keep_synonyms),
        url = unname(mo_url(y, open = FALSE, keep_synonyms = keep_synonyms)),
        ref = mo_ref(y, keep_synonyms = keep_synonyms),
        snomed = unlist(mo_snomed(y, keep_synonyms = keep_synonyms))
      )
    )
  })
  if (length(info) > 1) {
    names(info) <- mo_name(x)
    result <- info
  } else {
    result <- info[[1L]]
  }

  load_mo_uncertainties(metadata)
  result
}

#' @rdname mo_property
#' @export
mo_url <- function(x, open = FALSE, language = get_AMR_locale(), keep_synonyms = getOption("AMR_keep_synonyms", FALSE), ...) {
  if (missing(x)) {
    # this tries to find the data and an <mo> column
    x <- find_mo_col(fn = "mo_url")
  }
  meet_criteria(x, allow_NA = TRUE)
  meet_criteria(open, allow_class = "logical", has_length = 1)
  language <- validate_language(language)
  meet_criteria(keep_synonyms, allow_class = "logical", has_length = 1)

  x.mo <- as.mo(x = x, language = language, keep_synonyms = keep_synonyms, ... = ...)
  metadata <- get_mo_uncertainties()

  x.rank <- AMR::microorganisms$rank[match(x.mo, AMR::microorganisms$mo)]
  x.name <- AMR::microorganisms$fullname[match(x.mo, AMR::microorganisms$mo)]
  x.lpsn <- AMR::microorganisms$lpsn[match(x.mo, AMR::microorganisms$mo)]
  x.gbif <- AMR::microorganisms$gbif[match(x.mo, AMR::microorganisms$mo)]

  u <- character(length(x))
  u[!is.na(x.gbif)] <- paste0(TAXONOMY_VERSION$GBIF$url, "/species/", x.gbif[!is.na(x.gbif)])
  # overwrite with LPSN:
  u[!is.na(x.lpsn)] <- paste0(TAXONOMY_VERSION$LPSN$url, "/", x.rank[!is.na(x.lpsn)], "/", gsub(" ", "-", tolower(x.name[!is.na(x.lpsn)]), fixed = TRUE))

  names(u) <- x.name

  if (isTRUE(open)) {
    if (length(u) > 1) {
      warning_("in `mo_url()`: only the first URL will be opened, as `browseURL()` only suports one string.")
    }
    utils::browseURL(u[1L])
  }

  load_mo_uncertainties(metadata)
  u
}


#' @rdname mo_property
#' @export
mo_property <- function(x, property = "fullname", language = get_AMR_locale(), keep_synonyms = getOption("AMR_keep_synonyms", FALSE), ...) {
  if (missing(x)) {
    # this tries to find the data and an <mo> column
    x <- find_mo_col(fn = "mo_property")
  }
  meet_criteria(x, allow_NA = TRUE)
  meet_criteria(property, allow_class = "character", has_length = 1, is_in = colnames(AMR::microorganisms))
  language <- validate_language(language)
  meet_criteria(keep_synonyms, allow_class = "logical", has_length = 1)

  translate_into_language(mo_validate(x = x, property = property, language = language, keep_synonyms = keep_synonyms, ...), language = language, only_unknown = TRUE)
}

mo_validate <- function(x, property, language, keep_synonyms = keep_synonyms, ...) {

  # try to catch an error when inputting an invalid argument
  # so the 'call.' can be set to FALSE
  tryCatch(x[1L] %in% unlist(AMR::microorganisms[1, property, drop = TRUE]),
    error = function(e) stop(e$message, call. = FALSE)
  )

  dots <- list(...)
  Becker <- dots$Becker
  if (is.null(Becker) || property %in% c("kingdom", "phylum", "class", "order", "family", "genus")) {
    Becker <- FALSE
  }
  Lancefield <- dots$Lancefield
  if (is.null(Lancefield) || property %in% c("kingdom", "phylum", "class", "order", "family", "genus")) {
    Lancefield <- FALSE
  }
  has_Becker_or_Lancefield <- Becker %in% c(TRUE, "all") || Lancefield %in% c(TRUE, "all")

  # get microorganisms data set, but remove synonyms if keep_synonyms is FALSE
  mo_data_check <- AMR::microorganisms[which(AMR::microorganisms$status %in% if (isTRUE(keep_synonyms)) c("synonym", "accepted") else "accepted"), , drop = FALSE]

  if (all(x %in% c(mo_data_check$mo, NA)) && !has_Becker_or_Lancefield) {
    # do nothing, just don't run the other if-else's
  } else if (all(x %in% c(unlist(mo_data_check[[property]]), NA)) && !has_Becker_or_Lancefield) {
    # no need to do anything, just return it
    return(x)
  } else {
    # we need to get MO codes now
    x <- replace_old_mo_codes(x, property = property)
    x <- as.mo(x, language = language, keep_synonyms = keep_synonyms, ...)
  }

  # get property reeaaally fast using match()
  x <- AMR::microorganisms[[property]][match(x, AMR::microorganisms$mo)]

  if (property == "mo") {
    return(set_clean_class(x, new_class = c("mo", "character")))
  } else if (property == "snomed") {
    return(sort(as.character(eval(parse(text = x)))))
  } else {
    # everything else is character
    return(as.character(x))
  }
}

find_mo_col <- function(fn) {
  # this function tries to find an mo column in the data the function was called in,
  # which is useful when functions are used within dplyr verbs
  df <- get_current_data(arg_name = "x", call = -3) # will return an error if not found
  mo <- NULL
  try(
    {
      mo <- suppressMessages(search_type_in_df(df, "mo"))
    },
    silent = TRUE
  )
  if (!is.null(df) && !is.null(mo) && is.data.frame(df)) {
    if (message_not_thrown_before(fn = fn)) {
      message_("Using column '", font_bold(mo), "' as input for `", fn, "()`")
    }
    return(df[, mo, drop = TRUE])
  } else {
    stop_("argument `x` is missing and no column with info about microorganisms could be found.", call = -2)
  }
}
