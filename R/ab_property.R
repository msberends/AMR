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

#' Get Properties of an Antibiotic
#'
#' Use these functions to return a specific property of an antibiotic from the [antimicrobials] data set. All input values will be evaluated internally with [as.ab()].
#' @param x Any (vector of) text that can be coerced to a valid antibiotic drug code with [as.ab()].
#' @param tolower A [logical] to indicate whether the first [character] of every output should be transformed to a lower case [character]. This will lead to e.g. "polymyxin B" and not "polymyxin b".
#' @param property One of the column names of one of the [antimicrobials] data set: `vector_or(colnames(antimicrobials), sort = FALSE)`.
#' @param language Language of the returned text - the default is the current system language (see [get_AMR_locale()]) and can also be set with the package option [`AMR_locale`][AMR-options]. Use `language = NULL` or `language = ""` to prevent translation.
#' @param administration Way of administration, either `"oral"` or `"iv"`.
#' @param open Browse the URL using [utils::browseURL()].
#' @param ... In case of [set_ab_names()] and `data` is a [data.frame]: columns to select (supports tidy selection such as `column1:column4`), otherwise other arguments passed on to [as.ab()].
#' @param data A [data.frame] of which the columns need to be renamed, or a [character] vector of column names.
#' @param snake_case A [logical] to indicate whether the names should be in so-called [snake case](https://en.wikipedia.org/wiki/Snake_case): in lower case and all spaces/slashes replaced with an underscore (`_`).
#' @param only_first A [logical] to indicate whether only the first ATC code must be returned, with giving preference to J0-codes (i.e., the antimicrobial drug group).
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
#' @seealso [antimicrobials]
#' @inheritSection AMR Download Our Reference Data
#' @examples
#' # all properties:
#' ab_name("AMX")
#' ab_atc("AMX")
#' ab_cid("AMX")
#' ab_synonyms("AMX")
#' ab_tradenames("AMX")
#' ab_group("AMX")
#' ab_atc_group1("AMX")
#' ab_atc_group2("AMX")
#' ab_url("AMX")
#'
#' # smart lowercase transformation
#' ab_name(x = c("AMC", "PLB"))
#' ab_name(x = c("AMC", "PLB"), tolower = TRUE)
#'
#' # defined daily doses (DDD)
#' ab_ddd("AMX", "oral")
#' ab_ddd_units("AMX", "oral")
#' ab_ddd("AMX", "iv")
#' ab_ddd_units("AMX", "iv")
#'
#' ab_info("AMX") # all properties as a list
#'
#' # all ab_* functions use as.ab() internally, so you can go from 'any' to 'any':
#' ab_atc("AMP")
#' ab_group("J01CA01")
#' ab_loinc("ampicillin")
#' ab_name("21066-6")
#' ab_name(6249)
#' ab_name("J01CA01")
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
#'   example_isolates %>%
#'     set_ab_names(where(is.sir)) %>%
#'     colnames()
#'
#'   example_isolates %>%
#'     set_ab_names(NIT:VAN) %>%
#'     colnames()
#' }
#' }
ab_name <- function(x, language = get_AMR_locale(), tolower = FALSE, ...) {
  meet_criteria(x, allow_NA = TRUE)
  language <- validate_language(language)
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
ab_group <- function(x, language = get_AMR_locale(), all_groups = FALSE, ...) {
  meet_criteria(x, allow_NA = TRUE)
  language <- validate_language(language)
  meet_criteria(all_groups, allow_class = "logical", has_length = 1)

  grps <- ab_validate(x = x, property = "group", ...)
  for (i in seq_along(grps)) {
    # take the first match based on ABX_PRIORITY_LIST
    if (all_groups == FALSE) {
      grps[[i]] <- grps[[i]][1]
    } else if (length(grps[[i]]) > 1) {
      grps[[i]] <- grps[[i]][grps[[i]] != "Beta-lactamase inhibitors"] # leave these out if there are other groups
    }
    if (language != "en") {
      grps[[i]] <- translate_into_language(grps[[i]], language = language, only_affect_ab_names = TRUE)
    }
  }
  names(grps) <- x
  if (length(grps) == 1 || all_groups == FALSE) {
    unname(unlist(grps))
  } else {
    grps
  }
}

#' @rdname ab_property
#' @aliases ATC
#' @export
ab_atc <- function(x, only_first = FALSE, ...) {
  meet_criteria(x, allow_NA = TRUE)
  meet_criteria(only_first, allow_class = "logical", has_length = 1)

  atcs <- ab_validate(x = x, property = "atc", ...)

  if (only_first == TRUE) {
    atcs <- vapply(
      FUN.VALUE = character(1),
      # get only the first ATC code
      atcs,
      function(x) {
        # try to get the J-group
        if (any(x %like% "^J")) {
          x[x %like% "^J"][1L]
        } else {
          as.character(x[1L])
        }
      }
    )
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
  language <- validate_language(language)
  translate_into_language(ab_validate(x = x, property = "atc_group1", ...), language = language, only_affect_ab_names = TRUE)
}

#' @rdname ab_property
#' @export
ab_atc_group2 <- function(x, language = get_AMR_locale(), ...) {
  meet_criteria(x, allow_NA = TRUE)
  language <- validate_language(language)
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
  ddd_prop <- paste0(administration, "_ddd")
  out <- ab_validate(x = x, property = ddd_prop)

  if (any(ab_name(x, language = NULL) %like% "/" & is.na(out))) {
    warning_(
      "in `ab_ddd()`: DDDs of some combined products are available for different dose combinations and not (yet) part of the AMR package.",
      "Please refer to the WHOCC website:\n",
      "atcddd.fhi.no/ddd/list_of_ddds_combined_products/"
    )
  }
  out
}

#' @rdname ab_property
#' @export
ab_ddd_units <- function(x, administration = "oral", ...) {
  meet_criteria(x, allow_NA = TRUE)
  meet_criteria(administration, is_in = c("oral", "iv"), has_length = 1)

  x <- as.ab(x, ...)
  ddd_prop <- paste0(administration, "_units")
  out <- ab_validate(x = x, property = ddd_prop)

  if (any(ab_name(x, language = NULL) %like% "/" & is.na(out))) {
    warning_(
      "in `ab_ddd_units()`: DDDs of some combined products are available for different dose combinations and not (yet) part of the AMR package.",
      "Please refer to the WHOCC website:\n",
      "atcddd.fhi.no/ddd/list_of_ddds_combined_products/"
    )
  }
  out
}

#' @rdname ab_property
#' @export
ab_info <- function(x, language = get_AMR_locale(), ...) {
  meet_criteria(x, allow_NA = TRUE)
  language <- validate_language(language)

  x <- as.ab(x, ...)
  list(
    ab = as.character(x),
    cid = ab_cid(x),
    name = ab_name(x, language = language),
    group = ab_group(x, language = language),
    atc = ab_atc(x),
    atc_group1 = ab_atc_group1(x, language = language),
    atc_group2 = ab_atc_group2(x, language = language),
    tradenames = ab_tradenames(x),
    loinc = ab_loinc(x),
    ddd = list(
      oral = list(
        amount = ab_ddd(x, administration = "oral"),
        units = ab_ddd_units(x, administration = "oral")
      ),
      iv = list(
        amount = ab_ddd(x, administration = "iv"),
        units = ab_ddd_units(x, administration = "iv")
      )
    )
  )
}


#' @rdname ab_property
#' @export
ab_url <- function(x, open = FALSE, ...) {
  meet_criteria(x, allow_NA = TRUE)
  meet_criteria(open, allow_class = "logical", has_length = 1)

  ab <- as.ab(x = x, ...)
  atcs <- ab_atc(ab, only_first = TRUE)
  u <- character(length(atcs))
  # veterinary codes
  u[atcs %like% "^Q"] <- paste0("https://atcddd.fhi.no/atcvet/atcvet_index/?code=", atcs[atcs %like% "^Q"], "&showdescription=no")
  u[atcs %unlike% "^Q"] <- paste0("https://atcddd.fhi.no/atc_ddd_index//?code=", atcs[atcs %unlike% "^Q"], "&showdescription=no")
  u[is.na(atcs)] <- NA_character_
  names(u) <- ab_name(ab)

  NAs <- ab_name(ab, tolower = TRUE, language = NULL)[!is.na(ab) & is.na(atcs)]
  if (length(NAs) > 0) {
    warning_("in `ab_url()`: no ATC code available for ", vector_and(NAs, quotes = FALSE), ".")
  }

  if (open == TRUE) {
    if (length(u) > 1 && !is.na(u[1L])) {
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
  meet_criteria(property, is_in = colnames(AMR::antimicrobials), has_length = 1)
  language <- validate_language(language)
  translate_into_language(ab_validate(x = x, property = property, ...), language = language)
}

#' @rdname ab_property
#' @aliases ATC
#' @export
set_ab_names <- function(data, ..., property = "name", language = get_AMR_locale(), snake_case = NULL) {
  meet_criteria(data, allow_class = c("data.frame", "character"))
  meet_criteria(property, is_in = colnames(AMR::antimicrobials), has_length = 1, ignore.case = TRUE)
  language <- validate_language(language)
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
    if (tryCatch(length(c(...)) > 1, error = function(e) TRUE)) {
      df <- tryCatch(suppressWarnings(pm_select(data, ...)),
        error = function(e) {
          data[, c(...), drop = FALSE]
        }
      )
    } else if (tryCatch(is.character(c(...)), error = function(e) FALSE)) {
      df <- data[, c(...), drop = FALSE]
    } else {
      df <- data
    }
    vars <- get_column_abx(df, info = FALSE, only_sir_columns = FALSE, sort = FALSE, fn = "set_ab_names")
    if (length(vars) == 0) {
      message_("No columns with antibiotic results found for `set_ab_names()`, leaving names unchanged.")
      return(data)
    }
  } else {
    # quickly get antibiotic drug codes
    vars_ab <- as.ab(data, fast_mode = TRUE)
    vars <- data[!is.na(vars_ab)]
  }
  x <- vapply(
    FUN.VALUE = character(1),
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
    USE.NAMES = FALSE
  )
  if (any(x %in% c("", NA))) {
    warning_(
      "in `set_ab_names()`: no ", property, " found for column(s): ",
      vector_and(vars[x %in% c("", NA)], sort = FALSE)
    )
    x[x %in% c("", NA)] <- vars[x %in% c("", NA)]
  }

  if (snake_case == TRUE) {
    x <- tolower(gsub("[^a-zA-Z0-9]+", "_", x))
  }

  if (anyDuplicated(x)) {
    # very hacky way of adding the index to each duplicate
    # so      "Amoxicillin", "Amoxicillin",   "Amoxicillin"
    # will be "Amoxicillin", "Amoxicillin_2", "Amoxicillin_3"
    invisible(lapply(
      unique(x),
      function(u) {
        dups <- which(x == u)
        if (length(dups) > 1) {
          # there are duplicates
          dup_add_int <- dups[2:length(dups)]
          x[dup_add_int] <<- paste0(x[dup_add_int], "_", 2:length(dups))
        }
      }
    ))
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
  if (tryCatch(all(x[!is.na(x)] %in% AMR_env$AB_lookup$ab), error = function(e) FALSE)) {
    # special case for ab_* functions where class is already 'ab'
    x <- AMR_env$AB_lookup[match(x, AMR_env$AB_lookup$ab), property, drop = TRUE]
  } else {
    # try to catch an error when inputting an invalid argument
    # so the 'call.' can be set to FALSE
    tryCatch(x[1L] %in% AMR_env$AB_lookup[1, property, drop = TRUE],
      error = function(e) stop(conditionMessage(e), call. = FALSE)
    )

    if (!all(x %in% AMR_env$AB_lookup[, property, drop = TRUE])) {
      x <- as.ab(x, ...)
      if (all(is.na(x)) && is.list(AMR_env$AB_lookup[, property, drop = TRUE])) {
        x <- rep(NA_character_, length(x))
      } else {
        x <- AMR_env$AB_lookup[match(x, AMR_env$AB_lookup$ab), property, drop = TRUE]
      }
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
