# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Data Analysis for R                   #
#                                                                      #
# SOURCE                                                               #
# https://github.com/msberends/AMR                                     #
#                                                                      #
# LICENCE                                                              #
# (c) 2018-2022 Berends MS, Luz CF et al.                              #
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

#' Get Properties of an Antibiotic
#'
#' Use these functions to return a specific property of an antibiotic from the [antibiotics] data set. All input values will be evaluated internally with [as.ab()].
#' @inheritSection lifecycle Stable Lifecycle
#' @param x any (vector of) text that can be coerced to a valid antibiotic code with [as.ab()]
#' @param tolower a [logical] to indicate whether the first [character] of every output should be transformed to a lower case [character]. This will lead to e.g. "polymyxin B" and not "polymyxin b".
#' @param property one of the column names of one of the [antibiotics] data set: `vector_or(colnames(antibiotics), sort = FALSE)`.
#' @param language language of the returned text, defaults to system language (see [get_AMR_locale()]) and can also be set with `getOption("AMR_locale")`. Use `language = NULL` or `language = ""` to prevent translation.
#' @param administration way of administration, either `"oral"` or `"iv"`
#' @param open browse the URL using [utils::browseURL()]
#' @param ... in case of [set_ab_names()] and `data` is a [data.frame]: variables to select (supports tidy selection such as `column1:column4`), otherwise other arguments passed on to [as.ab()]
#' @param data a [data.frame] of which the columns need to be renamed, or a [character] vector of column names
#' @param snake_case a [logical] to indicate whether the names should be in so-called [snake case](https://en.wikipedia.org/wiki/Snake_case): in lower case and all spaces/slashes replaced with an underscore (`_`)
#' @param only_first a [logical] to indicate whether only the first ATC code must be returned, with giving preference to J0-codes (i.e., the antimicrobial drug group)
#' @details All output [will be translated][translate] where possible.
#' 
#' The function [ab_url()] will return the direct URL to the official WHO website. A warning will be returned if the required ATC code is not available.
#' 
#' The function [set_ab_names()] is a special column renaming function for [data.frame]s. It renames columns names that resemble antimicrobial drugs. It always makes sure that the new column names are unique. If `property = "atc"` is set, preference is given to ATC codes from the J-group.
#' @inheritSection as.ab Source
#' @rdname ab_property
#' @name ab_property
#' @return 
#' - An [integer] in case of [ab_cid()]
#' - A named [list] in case of [ab_info()] and multiple [ab_atc()]/[ab_synonyms()]/[ab_tradenames()]
#' - A [double] in case of [ab_ddd()]
#' - A [data.frame] in case of [set_ab_names()]
#' - A [character] in all other cases
#' @export
#' @seealso [antibiotics]
#' @inheritSection AMR Reference Data Publicly Available
#' @inheritSection AMR Read more on Our Website!
#' @examples
#' # all properties:
#' ab_name("AMX")       # "Amoxicillin"
#' ab_atc("AMX")        # "J01CA04" (ATC code from the WHO)
#' ab_cid("AMX")        # 33613 (Compound ID from PubChem)
#' ab_synonyms("AMX")   # a list with brand names of amoxicillin
#' ab_tradenames("AMX") # same
#' ab_group("AMX")      # "Beta-lactams/penicillins"
#' ab_atc_group1("AMX") # "Beta-lactam antibacterials, penicillins"
#' ab_atc_group2("AMX") # "Penicillins with extended spectrum"
#' ab_url("AMX")        # link to the official WHO page
#'
#' # smart lowercase tranformation
#' ab_name(x = c("AMC", "PLB"))  # "Amoxicillin/clavulanic acid" "Polymyxin B"
#' ab_name(x = c("AMC", "PLB"),
#'         tolower = TRUE)       # "amoxicillin/clavulanic acid" "polymyxin B"
#'
#' # defined daily doses (DDD)
#' ab_ddd("AMX", "oral")       #  1.5
#' ab_ddd_units("AMX", "oral") # "g"
#' ab_ddd("AMX", "iv")         #  3
#' ab_ddd_units("AMX", "iv")   # "g"
#'
#' ab_info("AMX")       # all properties as a list
#'
#' # all ab_* functions use as.ab() internally, so you can go from 'any' to 'any':
#' ab_atc("AMP")           # ATC code of AMP (ampicillin)
#' ab_group("J01CA01")     # Drug group of ampicillins ATC code
#' ab_loinc("ampicillin")  # LOINC codes of ampicillin
#' ab_name("21066-6")      # "Ampicillin" (using LOINC)
#' ab_name(6249)           # "Ampicillin" (using CID)
#' ab_name("J01CA01")      # "Ampicillin" (using ATC)
#' 
#' # spelling from different languages and dyslexia are no problem
#' ab_atc("ceftriaxon")
#' ab_atc("cephtriaxone")
#' ab_atc("cephthriaxone")
#' ab_atc("seephthriaaksone")
#' 
#' # use set_ab_names() for renaming columns
#' colnames(example_isolates)
#' colnames(set_ab_names(example_isolates))
#' colnames(set_ab_names(example_isolates, NIT:VAN))
#' \donttest{
#' if (require("dplyr")) {
#'   example_isolates %>%
#'     set_ab_names()
#'     
#'   # this does the same:
#'   example_isolates %>%
#'     rename_with(set_ab_names)
#'   
#'   # set_ab_names() works with any AB property:
#'   example_isolates %>%
#'     set_ab_names(property = "atc")
#'     
#'  example_isolates %>%
#'    set_ab_names(where(is.rsi)) %>%
#'    colnames()
#'    
#'  example_isolates %>%
#'    set_ab_names(NIT:VAN) %>%
#'    colnames()
#' }
#' }
ab_name <- function(x, language = get_AMR_locale(), tolower = FALSE, ...) {
  meet_criteria(x, allow_NA = TRUE)
  meet_criteria(language, has_length = 1, is_in = c(LANGUAGES_SUPPORTED, ""), allow_NULL = TRUE, allow_NA = TRUE)
  meet_criteria(tolower, allow_class = "logical", has_length = 1)
  
  x <- translate_into_language(ab_validate(x = x, property = "name", ...), language = language, only_affect_ab_names = TRUE)
  if (tolower == TRUE) {
    # use perl to only transform the first character
    # as we want "polymyxin B", not "polymyxin b"
    x <- gsub("^([A-Z])", "\\L\\1", x, perl = TRUE)
  }
  x
}

#' @rdname ab_property
#' @export
ab_cid <- function(x, ...) {
  meet_criteria(x, allow_NA = TRUE)
  ab_validate(x = x, property = "cid", ...)
}

#' @rdname ab_property
#' @export
ab_synonyms <- function(x, ...) {
  meet_criteria(x, allow_NA = TRUE)
  syns <- ab_validate(x = x, property = "synonyms", ...)
  names(syns) <- x
  if (length(syns) == 1) {
    unname(unlist(syns))
  } else {
    syns
  }
}

#' @rdname ab_property
#' @export
ab_tradenames <- function(x, ...) {
  meet_criteria(x, allow_NA = TRUE)
  ab_synonyms(x, ...)
}

#' @rdname ab_property
#' @export
ab_group <- function(x, language = get_AMR_locale(), ...) {
  meet_criteria(x, allow_NA = TRUE)
  meet_criteria(language, has_length = 1, is_in = c(LANGUAGES_SUPPORTED, ""), allow_NULL = TRUE, allow_NA = TRUE)
  translate_into_language(ab_validate(x = x, property = "group", ...), language = language, only_affect_ab_names = TRUE)
}

#' @rdname ab_property
#' @aliases ATC
#' @export
ab_atc <- function(x, only_first = FALSE, ...) {
  meet_criteria(x, allow_NA = TRUE)
  meet_criteria(only_first, allow_class = "logical", has_length = 1)
  
  atcs <- ab_validate(x = x, property = "atc", ...)
  
  if (only_first == TRUE) {
    atcs <- vapply(FUN.VALUE = character(1),
                   # get only the first ATC code
                   atcs,
                   function(x) {
                     # try to get the J-group
                     if (any(x %like% "^J")) {
                       x[x %like% "^J"][1L]
                     } else {
                       as.character(x[1L])
                     }
                   })
  } else if (length(atcs) == 1) {
    atcs <- unname(unlist(atcs))
  } else {
    names(atcs) <- x
  }
  
  atcs
}

#' @rdname ab_property
#' @export
ab_atc_group1 <- function(x, language = get_AMR_locale(), ...) {
  meet_criteria(x, allow_NA = TRUE)
  meet_criteria(language, has_length = 1, is_in = c(LANGUAGES_SUPPORTED, ""), allow_NULL = TRUE, allow_NA = TRUE)
  translate_into_language(ab_validate(x = x, property = "atc_group1", ...), language = language, only_affect_ab_names = TRUE)
}

#' @rdname ab_property
#' @export
ab_atc_group2 <- function(x, language = get_AMR_locale(), ...) {
  meet_criteria(x, allow_NA = TRUE)
  meet_criteria(language, has_length = 1, is_in = c(LANGUAGES_SUPPORTED, ""), allow_NULL = TRUE, allow_NA = TRUE)
  translate_into_language(ab_validate(x = x, property = "atc_group2", ...), language = language, only_affect_ab_names = TRUE)
}

#' @rdname ab_property
#' @export
ab_loinc <- function(x, ...) {
  meet_criteria(x, allow_NA = TRUE)
  loincs <- ab_validate(x = x, property = "loinc", ...)
  names(loincs) <- x
  if (length(loincs) == 1) {
    unname(unlist(loincs))
  } else {
    loincs
  }
}

#' @rdname ab_property
#' @export
ab_ddd <- function(x, administration = "oral", ...) {
  meet_criteria(x, allow_NA = TRUE)
  meet_criteria(administration, is_in = c("oral", "iv"), has_length = 1)
  
  x <- as.ab(x, ...)
  ddd_prop <- administration
  # old behaviour
  units <- list(...)$units
  if (!is.null(units) && isTRUE(units)) {
    if (message_not_thrown_before("ab_ddd", entire_session = TRUE)) {
      warning_("in `ab_ddd()`: using `ab_ddd(..., units = TRUE)` is deprecated, use `ab_ddd_units()` to retrieve units instead.",
               "This warning will be shown once per session.")
    }
    ddd_prop <- paste0(ddd_prop, "_units")
  } else {
    ddd_prop <- paste0(ddd_prop, "_ddd")
  }
  out <- ab_validate(x = x, property = ddd_prop)
  
  if (any(ab_name(x, language = NULL) %like% "/" & is.na(out))) {
    warning_("in `ab_ddd()`: DDDs of some combined products are available for different dose combinations and not (yet) part of the AMR package.",
             "Please refer to the WHOCC website:\n",
             "www.whocc.no/ddd/list_of_ddds_combined_products/")
  }
  out
}

#' @rdname ab_property
#' @export
ab_ddd_units <- function(x, administration = "oral", ...) {
  meet_criteria(x, allow_NA = TRUE)
  meet_criteria(administration, is_in = c("oral", "iv"), has_length = 1)
  
  x <- as.ab(x, ...)
  if (any(ab_name(x, language = NULL) %like% "/")) {
    warning_("in `ab_ddd_units()`: DDDs of combined products are available for different dose combinations and not (yet) part of the AMR package.",
             "Please refer to the WHOCC website:\n",
             "www.whocc.no/ddd/list_of_ddds_combined_products/")
  }
  
  ddd_prop <- paste0(administration, "_units")
  ab_validate(x = x, property = ddd_prop)
}

#' @rdname ab_property
#' @export
ab_info <- function(x, language = get_AMR_locale(), ...) {
  meet_criteria(x, allow_NA = TRUE)
  meet_criteria(language, has_length = 1, is_in = c(LANGUAGES_SUPPORTED, ""), allow_NULL = TRUE, allow_NA = TRUE)
  
  x <- as.ab(x, ...)
  list(ab = as.character(x),
       cid = ab_cid(x),
       name = ab_name(x, language = language),
       group = ab_group(x, language = language),
       atc = ab_atc(x),
       atc_group1 = ab_atc_group1(x, language = language),
       atc_group2 = ab_atc_group2(x, language = language),
       tradenames = ab_tradenames(x),
       loinc = ab_loinc(x),
       ddd = list(oral = list(amount = ab_ddd(x, administration = "oral"),
                              units = ab_ddd_units(x, administration = "oral")),
                  iv = list(amount = ab_ddd(x, administration = "iv"),
                            units = ab_ddd_units(x, administration = "iv"))))
}


#' @rdname ab_property
#' @export
ab_url <- function(x, open = FALSE, ...) {
  meet_criteria(x, allow_NA = TRUE)
  meet_criteria(open, allow_class = "logical", has_length = 1)
  
  ab <- as.ab(x = x, ...)
  atcs <- ab_atc(ab, only_first = TRUE)
  u <- paste0("https://www.whocc.no/atc_ddd_index/?code=", atcs, "&showdescription=no")
  u[is.na(atcs)] <- NA_character_
  names(u) <- ab_name(ab)
  
  NAs <- ab_name(ab, tolower = TRUE, language = NULL)[!is.na(ab) & is.na(atcs)]
  if (length(NAs) > 0) {
    warning_("in `ab_url()`: no ATC code available for ", vector_and(NAs, quotes = FALSE), ".")
  }
  
  if (open == TRUE) {
    if (length(u) > 1 & !is.na(u[1L])) {
      warning_("in `ab_url()`: only the first URL will be opened, as `browseURL()` only suports one string.")
    }
    if (!is.na(u[1L])) {
      utils::browseURL(u[1L])
    }
  }
  u
}

#' @rdname ab_property
#' @export
ab_property <- function(x, property = "name", language = get_AMR_locale(), ...) {
  meet_criteria(x, allow_NA = TRUE)
  meet_criteria(property, is_in = colnames(antibiotics), has_length = 1)
  meet_criteria(language, is_in = c(LANGUAGES_SUPPORTED, ""), has_length = 1, allow_NULL = TRUE, allow_NA = TRUE)
  translate_into_language(ab_validate(x = x, property = property, ...), language = language)
}

#' @rdname ab_property
#' @aliases ATC
#' @export
set_ab_names <- function(data, ..., property = "name", language = get_AMR_locale(), snake_case = NULL) {
  meet_criteria(data, allow_class = c("data.frame", "character"))
  meet_criteria(property, is_in = colnames(antibiotics), has_length = 1, ignore.case = TRUE)
  meet_criteria(language, has_length = 1, is_in = c(LANGUAGES_SUPPORTED, ""), allow_NULL = TRUE, allow_NA = TRUE)
  meet_criteria(snake_case, allow_class = "logical", has_length = 1, allow_NULL = TRUE)
  
  x_deparsed <- deparse(substitute(data))
  if (length(x_deparsed) > 1 || any(x_deparsed %unlike% "[a-z]+")) {
    x_deparsed <- "your_data"
  }
  
  property <- tolower(property)
  if (is.null(snake_case)) {
    snake_case <- property == "name"
  }
  
  if (is.data.frame(data)) {
    if (tryCatch(length(list(...)) > 0, error = function(e) TRUE)) {
      df <- pm_select(data, ...)
    } else {
      df <- data
    }
    vars <- get_column_abx(df, info = FALSE, only_rsi_columns = FALSE, sort = FALSE, fn = "set_ab_names")
    if (length(vars) == 0) {
      message_("No columns with antibiotic results found for `set_ab_names()`, leaving names unchanged.")
      return(data)
    }
  } else {
    # quickly get antibiotic codes
    vars_ab <- as.ab(data, fast_mode = TRUE)
    vars <- data[!is.na(vars_ab)]
  }
  x <- vapply(FUN.VALUE = character(1),
              ab_property(vars, property = property, language = language),
              function(x) {
                if (property == "atc") {
                  # try to get the J-group
                  if (any(x %like% "^J")) {
                    x[x %like% "^J"][1L]
                  } else {
                    as.character(x[1L])
                  }
                } else {
                  as.character(x[1L])
                }
              },
              USE.NAMES = FALSE)
  if (any(x %in% c("", NA))) {
    warning_("in `set_ab_names()`: no ", property, " found for column(s): ",
             vector_and(vars[x %in% c("", NA)], sort = FALSE))
    x[x %in% c("", NA)] <- vars[x %in% c("", NA)]
  }
  
  if (snake_case == TRUE) {
    x <- tolower(gsub("[^a-zA-Z0-9]+", "_", x))
  }
  
  if (any(duplicated(x))) {
    # very hacky way of adding the index to each duplicate
    # so      "Amoxicillin", "Amoxicillin",   "Amoxicillin"
    # will be "Amoxicillin", "Amoxicillin_2", "Amoxicillin_3"
    invisible(lapply(unique(x),
                     function(u) {
                       dups <- which(x == u)
                       if (length(dups) > 1) {
                         # there are duplicates
                         dup_add_int <- dups[2:length(dups)]
                         x[dup_add_int] <<- paste0(x[dup_add_int], "_", c(2:length(dups)))
                       }
                     }))
  }
  if (is.data.frame(data)) {
    colnames(data)[colnames(data) %in% vars] <- x
    data
  } else {
    data[which(!is.na(vars_ab))] <- x
    data
  }
}

ab_validate <- function(x, property, ...) {
  
  check_dataset_integrity()
  
  if (tryCatch(all(x[!is.na(x)] %in% AB_lookup$ab), error = function(e) FALSE)) {
    # special case for ab_* functions where class is already <ab>
    x <- AB_lookup[match(x, AB_lookup$ab), property, drop = TRUE]
    
  } else {
    # try to catch an error when inputting an invalid argument
    # so the 'call.' can be set to FALSE
    tryCatch(x[1L] %in% antibiotics[1, property],
             error = function(e) stop(e$message, call. = FALSE))
    
    if (!all(x %in% AB_lookup[, property])) {
      x <- as.ab(x, ...)
      x <- AB_lookup[match(x, AB_lookup$ab), property, drop = TRUE]
    }
  }
  
  if (property == "ab") {
    return(set_clean_class(x, new_class = c("ab", "character")))
  } else if (property == "cid") {
    return(as.integer(x))
  } else if (property %like% "ddd") {
    return(as.double(x))
  } else {
    x[is.na(x)] <- NA
    return(x)
  }
}
