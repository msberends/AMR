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

#' Get Properties of an Antiviral Drug
#'
#' Use these functions to return a specific property of an antiviral drug from the [antivirals] data set. All input values will be evaluated internally with [as.av()].
#' @param x any (vector of) text that can be coerced to a valid antiviral drug code with [as.av()]
#' @param tolower a [logical] to indicate whether the first [character] of every output should be transformed to a lower case [character].
#' @param property one of the column names of one of the [antivirals] data set: `vector_or(colnames(antivirals), sort = FALSE)`.
#' @param language language of the returned text, defaults to system language (see [get_AMR_locale()]) and can also be set with `getOption("AMR_locale")`. Use `language = NULL` or `language = ""` to prevent translation.
#' @param administration way of administration, either `"oral"` or `"iv"`
#' @param open browse the URL using [utils::browseURL()]
#' @param ... other arguments passed on to [as.av()]
#' @details All output [will be translated][translate] where possible.
#'
#' The function [av_url()] will return the direct URL to the official WHO website. A warning will be returned if the required ATC code is not available.
#' @inheritSection as.av Source
#' @rdname av_property
#' @name av_property
#' @return
#' - An [integer] in case of [av_cid()]
#' - A named [list] in case of [av_info()] and multiple [av_atc()]/[av_synonyms()]/[av_tradenames()]
#' - A [double] in case of [av_ddd()]
#' - A [character] in all other cases
#' @export
#' @seealso [antivirals]
#' @inheritSection AMR Reference Data Publicly Available
#' @examples
#' # all properties:
#' av_name("ACI")
#' av_atc("ACI")
#' av_cid("ACI")
#' av_synonyms("ACI")
#' av_tradenames("ACI")
#' av_group("ACI")
#' av_url("ACI")
#'
#' # smart lowercase tranformation
#' av_name(x = c("ACI", "VALA"))
#' av_name(x = c("ACI", "VALA"), tolower = TRUE)
#'
#' # defined daily doses (DDD)
#' av_ddd("ACI", "oral")
#' av_ddd_units("ACI", "oral")
#' av_ddd("ACI", "iv")
#' av_ddd_units("ACI", "iv")
#'
#' av_info("ACI") # all properties as a list
#'
#' # all av_* functions use as.av() internally, so you can go from 'any' to 'any':
#' av_atc("ACI")
#' av_group("J05AB01")
#' av_loinc("abacavir")
#' av_name("29113-8")
#' av_name(135398513)
#' av_name("J05AB01")
av_name <- function(x, language = get_AMR_locale(), tolower = FALSE, ...) {
  meet_criteria(x, allow_NA = TRUE)
  language <- validate_language(language)
  meet_criteria(tolower, allow_class = "logical", has_length = 1)
  
  x <- translate_into_language(av_validate(x = x, property = "name", ...), language = language, only_affect_ab_names = TRUE)
  if (tolower == TRUE) {
    # use perl to only transform the first character
    # as we want "polymyxin B", not "polymyxin b"
    x <- gsub("^([A-Z])", "\\L\\1", x, perl = TRUE)
  }
  x
}

#' @rdname av_property
#' @export
av_cid <- function(x, ...) {
  meet_criteria(x, allow_NA = TRUE)
  av_validate(x = x, property = "cid", ...)
}

#' @rdname av_property
#' @export
av_synonyms <- function(x, ...) {
  meet_criteria(x, allow_NA = TRUE)
  syns <- av_validate(x = x, property = "synonyms", ...)
  names(syns) <- x
  if (length(syns) == 1) {
    unname(unlist(syns))
  } else {
    syns
  }
}

#' @rdname av_property
#' @export
av_tradenames <- function(x, ...) {
  meet_criteria(x, allow_NA = TRUE)
  av_synonyms(x, ...)
}

#' @rdname av_property
#' @export
av_group <- function(x, language = get_AMR_locale(), ...) {
  meet_criteria(x, allow_NA = TRUE)
  language <- validate_language(language)
  translate_into_language(av_validate(x = x, property = "atc_group", ...), language = language, only_affect_ab_names = TRUE)
}

#' @rdname av_property
#' @export
av_atc <- function(x, ...) {
  meet_criteria(x, allow_NA = TRUE)
  # ATCs in the antivirals data set are not a list
  av_validate(x = x, property = "atc", ...)
}

#' @rdname av_property
#' @export
av_loinc <- function(x, ...) {
  meet_criteria(x, allow_NA = TRUE)
  loincs <- av_validate(x = x, property = "loinc", ...)
  names(loincs) <- x
  if (length(loincs) == 1) {
    unname(unlist(loincs))
  } else {
    loincs
  }
}

#' @rdname av_property
#' @export
av_ddd <- function(x, administration = "oral", ...) {
  meet_criteria(x, allow_NA = TRUE)
  meet_criteria(administration, is_in = c("oral", "iv"), has_length = 1)
  
  x <- as.av(x, ...)
  ddd_prop <- paste0(administration, "_ddd")
  out <- av_validate(x = x, property = ddd_prop)
  
  if (any(av_name(x, language = NULL) %like% "/" & is.na(out))) {
    warning_(
      "in `av_ddd()`: DDDs of some combined products are available for different dose combinations and not (yet) part of the AMR package.",
      "Please refer to the WHOCC website:\n",
      "www.whocc.no/ddd/list_of_ddds_combined_products/"
    )
  }
  out
}

#' @rdname av_property
#' @export
av_ddd_units <- function(x, administration = "oral", ...) {
  meet_criteria(x, allow_NA = TRUE)
  meet_criteria(administration, is_in = c("oral", "iv"), has_length = 1)
  
  x <- as.av(x, ...)
  ddd_prop <- paste0(administration, "_units")
  out <- av_validate(x = x, property = ddd_prop)
  
  if (any(av_name(x, language = NULL) %like% "/" & is.na(out))) {
    warning_(
      "in `av_ddd_units()`: DDDs of some combined products are available for different dose combinations and not (yet) part of the AMR package.",
      "Please refer to the WHOCC website:\n",
      "www.whocc.no/ddd/list_of_ddds_combined_products/"
    )
  }
  out
}

#' @rdname av_property
#' @export
av_info <- function(x, language = get_AMR_locale(), ...) {
  meet_criteria(x, allow_NA = TRUE)
  language <- validate_language(language)
  
  x <- as.av(x, ...)
  list(
    av = as.character(x),
    cid = av_cid(x),
    name = av_name(x, language = language),
    group = av_group(x, language = language),
    atc = av_atc(x),
    tradenames = av_tradenames(x),
    loinc = av_loinc(x),
    ddd = list(
      oral = list(
        amount = av_ddd(x, administration = "oral"),
        units = av_ddd_units(x, administration = "oral")
      ),
      iv = list(
        amount = av_ddd(x, administration = "iv"),
        units = av_ddd_units(x, administration = "iv")
      )
    )
  )
}


#' @rdname av_property
#' @export
av_url <- function(x, open = FALSE, ...) {
  meet_criteria(x, allow_NA = TRUE)
  meet_criteria(open, allow_class = "logical", has_length = 1)
  
  av <- as.av(x = x, ...)
  atcs <- av_atc(av, only_first = TRUE)
  u <- paste0("https://www.whocc.no/atc_ddd_index/?code=", atcs, "&showdescription=no")
  u[is.na(atcs)] <- NA_character_
  names(u) <- av_name(av)
  
  NAs <- av_name(av, tolower = TRUE, language = NULL)[!is.na(av) & is.na(atcs)]
  if (length(NAs) > 0) {
    warning_("in `av_url()`: no ATC code available for ", vector_and(NAs, quotes = FALSE), ".")
  }
  
  if (open == TRUE) {
    if (length(u) > 1 && !is.na(u[1L])) {
      warning_("in `av_url()`: only the first URL will be opened, as `browseURL()` only suports one string.")
    }
    if (!is.na(u[1L])) {
      utils::browseURL(u[1L])
    }
  }
  u
}

#' @rdname av_property
#' @export
av_property <- function(x, property = "name", language = get_AMR_locale(), ...) {
  meet_criteria(x, allow_NA = TRUE)
  meet_criteria(property, is_in = colnames(AMR::antivirals), has_length = 1)
  meet_criteria(language, is_in = c(LANGUAGES_SUPPORTED, ""), has_length = 1, allow_NULL = TRUE, allow_NA = TRUE)
  translate_into_language(av_validate(x = x, property = property, ...), language = language)
}

av_validate <- function(x, property, ...) {
  if (tryCatch(all(x[!is.na(x)] %in% AMR_env$AV_lookup$av), error = function(e) FALSE)) {
    # special case for av_* functions where class is already 'av'
    x <- AMR_env$AV_lookup[match(x, AMR_env$AV_lookup$av), property, drop = TRUE]
  } else {
    # try to catch an error when inputting an invalid argument
    # so the 'call.' can be set to FALSE
    tryCatch(x[1L] %in% AMR_env$AV_lookup[1, property, drop = TRUE],
             error = function(e) stop(e$message, call. = FALSE)
    )
    
    if (!all(x %in% AMR_env$AV_lookup[, property, drop = TRUE])) {
      x <- as.av(x, ...)
      if (all(is.na(x)) && is.list(AMR_env$AV_lookup[, property, drop = TRUE])) {
        x <- rep(NA_character_, length(x))
      } else {
        x <- AMR_env$AV_lookup[match(x, AMR_env$AV_lookup$av), property, drop = TRUE]
      }
    }
  }
  
  if (property == "av") {
    return(set_clean_class(x, new_class = c("av", "character")))
  } else if (property == "cid") {
    return(as.integer(x))
  } else if (property %like% "ddd") {
    return(as.double(x))
  } else {
    x[is.na(x)] <- NA
    return(x)
  }
}
