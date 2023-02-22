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

#' Get Properties of a Microorganism
#'
#' Use these functions to return a specific property of a microorganism based on the latest accepted taxonomy. All input values will be evaluated internally with [as.mo()], which makes it possible to use microbial abbreviations, codes and names as input. See *Examples*.
#' @param x any [character] (vector) that can be coerced to a valid microorganism code with [as.mo()]. Can be left blank for auto-guessing the column containing microorganism codes if used in a data set, see *Examples*.
#' @param property one of the column names of the [microorganisms] data set: `r vector_or(colnames(microorganisms), sort = FALSE, quotes = TRUE)`, or must be `"shortname"`
#' @inheritParams as.mo
#' @param ... other arguments passed on to [as.mo()], such as 'minimum_matching_score', 'ignore_pattern', and 'remove_from_input'
#' @param ab any (vector of) text that can be coerced to a valid antibiotic drug code with [as.ab()]
#' @param open browse the URL using [`browseURL()`][utils::browseURL()]
#' @details All functions will, at default, **not** keep old taxonomic properties, as synonyms are automatically replaced with the current taxonomy. Take for example *Enterobacter aerogenes*, which was initially named in 1960 but renamed to *Klebsiella aerogenes* in 2017:
#' - `mo_genus("Enterobacter aerogenes")` will return `"Klebsiella"` (with a note about the renaming)
#' - `mo_genus("Enterobacter aerogenes", keep_synonyms = TRUE)` will return `"Enterobacter"` (with a once-per-session warning that the name is outdated)
#' - `mo_ref("Enterobacter aerogenes")` will return `"Tindall et al., 2017"` (with a note)
#' - `mo_ref("Enterobacter aerogenes", keep_synonyms = TRUE)` will return `"Hormaeche et al., 1960"` (with a warning)
#'
#' The short name ([mo_shortname()]) returns the first character of the genus and the full species, such as `"E. coli"`, for species and subspecies. Exceptions are abbreviations of staphylococci (such as *"CoNS"*, Coagulase-Negative Staphylococci) and beta-haemolytic streptococci (such as *"GBS"*, Group B Streptococci). Please bear in mind that e.g. *E. coli* could mean *Escherichia coli* (kingdom of Bacteria) as well as *Entamoeba coli* (kingdom of Protozoa). Returning to the full name will be done using [as.mo()] internally, giving priority to bacteria and human pathogens, i.e. `"E. coli"` will be considered *Escherichia coli*. As a result, `mo_fullname(mo_shortname("Entamoeba coli"))` returns `"Escherichia coli"`.
#'
#' Since the top-level of the taxonomy is sometimes referred to as 'kingdom' and sometimes as 'domain', the functions [mo_kingdom()] and [mo_domain()] return the exact same results.
#'
#' Determination of human pathogenicity ([mo_pathogenicity()]) is strongly based on Bartlett *et al.* (2022, \doi{10.1099/mic.0.001269}). This function returns a [factor] with the levels *Pathogenic*, *Potentially pathogenic*, *Non-pathogenic*, and *Unknown*.
#'
#' Determination of the Gram stain ([mo_gramstain()]) will be based on the taxonomic kingdom and phylum. Originally, Cavalier-Smith defined the so-called subkingdoms Negibacteria and Posibacteria (2002, [PMID 11837318](https://pubmed.ncbi.nlm.nih.gov/11837318/)), and only considered these phyla as Posibacteria: Actinobacteria, Chloroflexi, Firmicutes, and Tenericutes. These phyla were later renamed to Actinomycetota, Chloroflexota, Bacillota, and Mycoplasmatota (2021, [PMID 34694987](https://pubmed.ncbi.nlm.nih.gov/34694987/)). Bacteria in these phyla are considered Gram-positive in this `AMR` package, except for members of the class Negativicutes (within phylum Bacillota) which are Gram-negative. All other bacteria are considered Gram-negative. Species outside the kingdom of Bacteria will return a value `NA`. Functions [mo_is_gram_negative()] and [mo_is_gram_positive()] always return `TRUE` or `FALSE` (or `NA` when the input is `NA` or the MO code is `UNKNOWN`), thus always return `FALSE` for species outside the taxonomic kingdom of Bacteria.
#'
#' Determination of yeasts ([mo_is_yeast()]) will be based on the taxonomic kingdom and class. *Budding yeasts* are fungi of the phylum Ascomycota, class Saccharomycetes (also called Hemiascomycetes). *True yeasts* are aggregated into the underlying order Saccharomycetales. Thus, for all microorganisms that are member of the taxonomic class Saccharomycetes, the function will return `TRUE`. It returns `FALSE` otherwise (or `NA` when the input is `NA` or the MO code is `UNKNOWN`).
#'
#' Determination of intrinsic resistance ([mo_is_intrinsic_resistant()]) will be based on the [intrinsic_resistant] data set, which is based on `r format_eucast_version_nr(3.3)`. The [mo_is_intrinsic_resistant()] function can be vectorised over both argument `x` (input for microorganisms) and `ab` (input for antibiotics).
#'
#' The function [mo_url()] will return the direct URL to the online database entry, which also shows the scientific reference of the concerned species.
#'
#' SNOMED codes ([mo_snomed()]) are from the version of `r documentation_date(TAXONOMY_VERSION$SNOMED$accessed_date)`. See *Source* and the [microorganisms] data set for more info.
#'
#' Old taxonomic names (so-called 'synonyms') can be retrieved with [mo_synonyms()] (which will have the scientific reference as [name][base::names()]), the current taxonomic name can be retrieved with [mo_current()]. Both functions return full names.
#'
#' All output [will be translated][translate] where possible.
#' @section Matching Score for Microorganisms:
#' This function uses [as.mo()] internally, which uses an advanced algorithm to translate arbitrary user input to valid taxonomy using a so-called matching score. You can read about this public algorithm on the [MO matching score page][mo_matching_score()].
#' @inheritSection as.mo Source
#' @rdname mo_property
#' @name mo_property
#' @return
#' - An [integer] in case of [mo_year()]
#' - An [ordered factor][factor] in case of [mo_pathogenicity()]
#' - A [list] in case of [mo_taxonomy()], [mo_synonyms()], [mo_snomed()] and [mo_info()]
#' - A named [character] in case of [mo_url()]
#' - A [character] in all other cases
#' @export
#' @seealso Data set [microorganisms]
#' @inheritSection AMR Reference Data Publicly Available
#' @examples
#' # taxonomic tree -----------------------------------------------------------
#'
#' mo_kingdom("Klebsiella pneumoniae")
#' mo_phylum("Klebsiella pneumoniae")
#' mo_class("Klebsiella pneumoniae")
#' mo_order("Klebsiella pneumoniae")
#' mo_family("Klebsiella pneumoniae")
#' mo_genus("Klebsiella pneumoniae")
#' mo_species("Klebsiella pneumoniae")
#' mo_subspecies("Klebsiella pneumoniae")
#'
#'
#' # full names and short names -----------------------------------------------
#'
#' mo_name("Klebsiella pneumoniae")
#' mo_fullname("Klebsiella pneumoniae")
#' mo_shortname("Klebsiella pneumoniae")
#'
#'
#' # other properties ---------------------------------------------------------
#'
#' mo_pathogenicity("Klebsiella pneumoniae")
#' mo_gramstain("Klebsiella pneumoniae")
#' mo_snomed("Klebsiella pneumoniae")
#' mo_type("Klebsiella pneumoniae")
#' mo_rank("Klebsiella pneumoniae")
#' mo_url("Klebsiella pneumoniae")
#' mo_is_yeast(c("Candida", "Trichophyton", "Klebsiella"))
#'
#'
#' # scientific reference -----------------------------------------------------
#'
#' mo_ref("Klebsiella aerogenes")
#' mo_authors("Klebsiella aerogenes")
#' mo_year("Klebsiella aerogenes")
#' mo_lpsn("Klebsiella aerogenes")
#' mo_gbif("Klebsiella aerogenes")
#' mo_synonyms("Klebsiella aerogenes")
#' 
#'
#' # abbreviations known in the field -----------------------------------------
#'
#' mo_genus("MRSA")
#' mo_species("MRSA")
#' mo_shortname("VISA")
#' mo_gramstain("VISA")
#'
#' mo_genus("EHEC")
#' mo_species("EIEC")
#' mo_name("UPEC")
#'
#'
#' # known subspecies ---------------------------------------------------------
#'
#' mo_fullname("K. pneu rh")
#' mo_shortname("K. pneu rh")
#'
#' \donttest{
#' # Becker classification, see ?as.mo ----------------------------------------
#'
#' mo_fullname("Staph epidermidis")
#' mo_fullname("Staph epidermidis", Becker = TRUE)
#' mo_shortname("Staph epidermidis")
#' mo_shortname("Staph epidermidis", Becker = TRUE)
#'
#'
#' # Lancefield classification, see ?as.mo ------------------------------------
#'
#' mo_fullname("Strep agalactiae")
#' mo_fullname("Strep agalactiae", Lancefield = TRUE)
#' mo_shortname("Strep agalactiae")
#' mo_shortname("Strep agalactiae", Lancefield = TRUE)
#'
#'
#' # language support  --------------------------------------------------------
#'
#' mo_gramstain("Klebsiella pneumoniae", language = "de") # German
#' mo_gramstain("Klebsiella pneumoniae", language = "nl") # Dutch
#' mo_gramstain("Klebsiella pneumoniae", language = "es") # Spanish
#' mo_gramstain("Klebsiella pneumoniae", language = "el") # Greek
#' mo_gramstain("Klebsiella pneumoniae", language = "uk") # Ukrainian
#'
#' # mo_type is equal to mo_kingdom, but mo_kingdom will remain untranslated
#' mo_kingdom("Klebsiella pneumoniae")
#' mo_type("Klebsiella pneumoniae")
#' mo_kingdom("Klebsiella pneumoniae", language = "zh") # Chinese, no effect
#' mo_type("Klebsiella pneumoniae", language = "zh") # Chinese, translated
#'
#' mo_fullname("S. pyogenes", Lancefield = TRUE, language = "de")
#' mo_fullname("S. pyogenes", Lancefield = TRUE, language = "uk")
#'
#'
#' # other --------------------------------------------------------------------
#'
#' # gram stains and intrinsic resistance can be used as a filter in dplyr verbs
#' if (require("dplyr")) {
#'   example_isolates %>%
#'     filter(mo_is_gram_positive()) %>%
#'     count(mo_genus(), sort = TRUE)
#' }
#' if (require("dplyr")) {
#'   example_isolates %>%
#'     filter(mo_is_intrinsic_resistant(ab = "vanco")) %>%
#'     count(mo_genus(), sort = TRUE)
#' }
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
    # this tries to find the data and an 'mo' column
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
    # this tries to find the data and an 'mo' column
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
  shortnames[shortnames %like_case% "S. Group [ABCDFGHK]"] <- paste0("G", gsub("S. Group ([ABCDFGHK])", "\\1", shortnames[shortnames %like_case% "S. Group [ABCDFGHK]"], perl = TRUE), "S")
  # unknown species etc.
  shortnames[shortnames %like% "unknown"] <- paste0("(", trimws2(gsub("[^a-zA-Z -]", "", shortnames[shortnames %like% "unknown"], perl = TRUE)), ")")

  shortnames[mo_rank(x.mo) %in% c("kingdom", "phylum", "class", "order", "family")] <- mo_name(x.mo, language = NULL, keep_synonyms = keep_synonyms)

  shortnames[is.na(x.mo)] <- NA_character_
  load_mo_uncertainties(metadata)
  translate_into_language(shortnames, language = language, only_unknown = FALSE, only_affect_mo_names = TRUE)
}



#' @rdname mo_property
#' @export
mo_subspecies <- function(x, language = get_AMR_locale(), keep_synonyms = getOption("AMR_keep_synonyms", FALSE), ...) {
  if (missing(x)) {
    # this tries to find the data and an 'mo' column
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
    # this tries to find the data and an 'mo' column
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
    # this tries to find the data and an 'mo' column
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
    # this tries to find the data and an 'mo' column
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
    # this tries to find the data and an 'mo' column
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
    # this tries to find the data and an 'mo' column
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
    # this tries to find the data and an 'mo' column
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
    # this tries to find the data and an 'mo' column
    x <- find_mo_col(fn = "mo_kingdom")
  }
  meet_criteria(x, allow_NA = TRUE)
  language <- validate_language(language)
  meet_criteria(keep_synonyms, allow_class = "logical", has_length = 1)

  translate_into_language(mo_validate(x = x, property = "kingdom", language = language, keep_synonyms = keep_synonyms, ...), language = language, only_unknown = TRUE)
}

#' @rdname mo_property
#' @export
mo_domain <- function(x, language = get_AMR_locale(), keep_synonyms = getOption("AMR_keep_synonyms", FALSE), ...) {
  if (missing(x)) {
    # this tries to find the data and an 'mo' column
    x <- find_mo_col(fn = "mo_domain")
  }
  mo_kingdom(x = x, language = language, keep_synonyms = keep_synonyms, ...)
}

#' @rdname mo_property
#' @export
mo_type <- function(x, language = get_AMR_locale(), keep_synonyms = getOption("AMR_keep_synonyms", FALSE), ...) {
  if (missing(x)) {
    # this tries to find the data and an 'mo' column
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
    # this tries to find the data and an 'mo' column
    x <- find_mo_col(fn = "mo_status")
  }
  meet_criteria(x, allow_NA = TRUE)
  language <- validate_language(language)
  meet_criteria(keep_synonyms, allow_class = "logical", has_length = 1)

  translate_into_language(mo_validate(x = x, property = "status", language = language, keep_synonyms = keep_synonyms, ...), language = language, only_unknown = TRUE)
}

#' @rdname mo_property
#' @export
mo_pathogenicity <- function(x, language = get_AMR_locale(), keep_synonyms = getOption("AMR_keep_synonyms", FALSE), ...) {
  if (missing(x)) {
    # this tries to find the data and an 'mo' column
    x <- find_mo_col(fn = "mo_pathogenicity")
  }
  meet_criteria(x, allow_NA = TRUE)
  language <- validate_language(language)
  meet_criteria(keep_synonyms, allow_class = "logical", has_length = 1)

  add_MO_lookup_to_AMR_env()

  x.mo <- as.mo(x, language = language, keep_synonyms = keep_synonyms, ...)
  metadata <- get_mo_uncertainties()

  prev <- AMR_env$MO_lookup$prevalence[match(x.mo, AMR_env$MO_lookup$mo)]
  kngd <- AMR_env$MO_lookup$kingdom[match(x.mo, AMR_env$MO_lookup$mo)]
  rank <- AMR_env$MO_lookup$rank[match(x.mo, AMR_env$MO_lookup$mo)]

  out <- factor(
    ifelse(prev == 1 & kngd == "Bacteria" & rank != "genus",
      "Pathogenic",
      ifelse(prev < 2 & kngd == "Fungi",
        "Potentially pathogenic",
        ifelse(prev == 2 & kngd == "Bacteria",
          "Non-pathogenic",
          ifelse(kngd == "Bacteria",
            "Potentially pathogenic",
            "Unknown"
          )
        )
      )
    ),
    levels = c("Pathogenic", "Potentially pathogenic", "Non-pathogenic", "Unknown"),
    ordered = TRUE
  )

  load_mo_uncertainties(metadata)
  out
}

#' @rdname mo_property
#' @export
mo_gramstain <- function(x, language = get_AMR_locale(), keep_synonyms = getOption("AMR_keep_synonyms", FALSE), ...) {
  if (missing(x)) {
    # this tries to find the data and an 'mo' column
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
    # no longer in use, does not hurt to keep here:
    "Actinobacteria",
    "Chloroflexi",
    "Firmicutes",
    "Tenericutes",
    "Actinomycetota", # since 2021, old name was Actinobacteria
    "Chloroflexota", # since 2021, old name was Chloroflexi
    "Bacillota", # since 2021, old name was Firmicutes
    "Mycoplasmatota" # since 2021, old name was Tenericutes
  ) &
    # but class Negativicutes (of phylum Bacillota) are Gram-negative!
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
    # this tries to find the data and an 'mo' column
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
    # this tries to find the data and an 'mo' column
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
    # this tries to find the data and an 'mo' column
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
    # this tries to find the data and an 'mo' column
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

  # runs against internal vector: intrinsic_resistant (see zzz.R)
  add_intrinsic_resistance_to_AMR_env()
  paste(x, ab) %in% AMR_env$intrinsic_resistant
}

#' @rdname mo_property
#' @export
mo_snomed <- function(x, language = get_AMR_locale(), keep_synonyms = getOption("AMR_keep_synonyms", FALSE), ...) {
  if (missing(x)) {
    # this tries to find the data and an 'mo' column
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
    # this tries to find the data and an 'mo' column
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
    # this tries to find the data and an 'mo' column
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
    # this tries to find the data and an 'mo' column
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
    # this tries to find the data and an 'mo' column
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
    # this tries to find the data and an 'mo' column
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
    # this tries to find the data and an 'mo' column
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
    # this tries to find the data and an 'mo' column
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
    # this tries to find the data and an 'mo' column
    x <- find_mo_col(fn = "mo_synonyms")
  }
  meet_criteria(x, allow_NA = TRUE)
  language <- validate_language(language)
  meet_criteria(keep_synonyms, allow_class = "logical", has_length = 1)

  add_MO_lookup_to_AMR_env()

  x.mo <- as.mo(x, language = language, keep_synonyms = keep_synonyms, ...)
  metadata <- get_mo_uncertainties()

  syns <- lapply(x.mo, function(y) {
    gbif <- AMR_env$MO_lookup$gbif[match(y, AMR_env$MO_lookup$mo)]
    lpsn <- AMR_env$MO_lookup$lpsn[match(y, AMR_env$MO_lookup$mo)]
    fullname <- AMR_env$MO_lookup[which(AMR_env$MO_lookup$lpsn_renamed_to == lpsn | AMR_env$MO_lookup$gbif_renamed_to == gbif), "fullname", drop = TRUE]
    if (length(fullname) == 0) {
      NULL
    } else {
      ref <- AMR_env$MO_lookup[which(AMR_env$MO_lookup$lpsn_renamed_to == lpsn | AMR_env$MO_lookup$gbif_renamed_to == gbif), "ref", drop = TRUE]
      names(fullname) <- ref
      fullname
    }
  })

  if (length(syns) == 1) {
    syns <- unlist(syns)
  }
  
  load_mo_uncertainties(metadata)
  syns
}

#' @rdname mo_property
#' @export
mo_current <- function(x, language = get_AMR_locale(), ...) {
  meet_criteria(x, allow_NA = TRUE)
  language <- validate_language(language)
  x.mo <- suppressWarnings(as.mo(x, keep_synonyms = TRUE, ...))
  out <- synonym_mo_to_accepted_mo(x.mo, fill_in_accepted = TRUE)
  mo_name(out, language = language)
}

#' @rdname mo_property
#' @export
mo_info <- function(x, language = get_AMR_locale(), keep_synonyms = getOption("AMR_keep_synonyms", FALSE), ...) {
  if (missing(x)) {
    # this tries to find the data and an 'mo' column
    x <- find_mo_col(fn = "mo_info")
  }
  meet_criteria(x, allow_NA = TRUE)
  language <- validate_language(language)
  meet_criteria(keep_synonyms, allow_class = "logical", has_length = 1)

  x <- as.mo(x, language = language, keep_synonyms = keep_synonyms, ...)
  metadata <- get_mo_uncertainties()

  info <- lapply(x, function(y) {
    c(
      list(mo = as.character(x)),
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
    # this tries to find the data and an 'mo' column
    x <- find_mo_col(fn = "mo_url")
  }
  meet_criteria(x, allow_NA = TRUE)
  meet_criteria(open, allow_class = "logical", has_length = 1)
  language <- validate_language(language)
  meet_criteria(keep_synonyms, allow_class = "logical", has_length = 1)

  add_MO_lookup_to_AMR_env()

  x.mo <- as.mo(x = x, language = language, keep_synonyms = keep_synonyms, ... = ...)
  metadata <- get_mo_uncertainties()

  x.rank <- AMR_env$MO_lookup$rank[match(x.mo, AMR_env$MO_lookup$mo)]
  x.name <- AMR_env$MO_lookup$fullname[match(x.mo, AMR_env$MO_lookup$mo)]
  x.lpsn <- AMR_env$MO_lookup$lpsn[match(x.mo, AMR_env$MO_lookup$mo)]
  x.gbif <- AMR_env$MO_lookup$gbif[match(x.mo, AMR_env$MO_lookup$mo)]

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
    # this tries to find the data and an 'mo' column
    x <- find_mo_col(fn = "mo_property")
  }
  meet_criteria(x, allow_NA = TRUE)
  meet_criteria(property, allow_class = "character", has_length = 1, is_in = colnames(AMR::microorganisms))
  language <- validate_language(language)
  meet_criteria(keep_synonyms, allow_class = "logical", has_length = 1)

  translate_into_language(mo_validate(x = x, property = property, language = language, keep_synonyms = keep_synonyms, ...), language = language, only_unknown = TRUE)
}

mo_validate <- function(x, property, language, keep_synonyms = keep_synonyms, ...) {
  add_MO_lookup_to_AMR_env()

  # try to catch an error when inputting an invalid argument
  # so the 'call.' can be set to FALSE
  tryCatch(x[1L] %in% unlist(AMR_env$MO_lookup[1, property, drop = TRUE]),
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
  mo_data_check <- AMR_env$MO_lookup[which(AMR_env$MO_lookup$status %in% if (isTRUE(keep_synonyms)) c("synonym", "accepted") else "accepted"), , drop = FALSE]

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
  if (property == "snomed") {
    x <- lapply(x, function(y) unlist(AMR_env$MO_lookup$snomed[match(y, AMR_env$MO_lookup$mo)]))
  } else {
    x <- AMR_env$MO_lookup[[property]][match(x, AMR_env$MO_lookup$mo)]
  }

  if (property == "mo") {
    return(set_clean_class(x, new_class = c("mo", "character")))
  } else if (property == "snomed") {
    return(x)
  } else if (property == "prevalence") {
    return(as.double(x))
  } else {
    # everything else as character
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
