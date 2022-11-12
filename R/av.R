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

#' Transform Input to an Antiviral Agent ID
#'
#' Use this function to determine the antiviral agent code of one or more antiviral agents. The data set [antivirals] will be searched for abbreviations, official names and synonyms (brand names).
#' @param x a [character] vector to determine to antiviral agent ID
#' @param flag_multiple_results a [logical] to indicate whether a note should be printed to the console that probably more than one antiviral agent code or name can be retrieved from a single input value.
#' @param info a [logical] to indicate whether a progress bar should be printed, defaults to `TRUE` only in interactive mode
#' @param ... arguments passed on to internal functions
#' @rdname as.av
#' @inheritSection WHOCC WHOCC
#' @details All entries in the [antivirals] data set have three different identifiers: a human readable EARS-Net code (column `ab`, used by ECDC and WHONET), an ATC code (column `atc`, used by WHO), and a CID code (column `cid`, Compound ID, used by PubChem). The data set contains more than 5,000 official brand names from many different countries, as found in PubChem. Not that some drugs contain multiple ATC codes.
#'
#' All these properties will be searched for the user input. The [as.av()] can correct for different forms of misspelling:
#'
#'  * Wrong spelling of drug names (such as "acyclovir"), which corrects for most audible similarities such as f/ph, x/ks, c/z/s, t/th, etc.
#'  * Too few or too many vowels or consonants
#'  * Switching two characters (such as "aycclovir", often the case in clinical data, when doctors typed too fast)
#'  * Digitalised paper records, leaving artefacts like 0/o/O (zero and O's), B/8, n/r, etc.
#'
#' Use the [`av_*`][av_property()] functions to get properties based on the returned antiviral agent ID, see *Examples*.
#'
#' Note: the [as.av()] and [`av_*`][av_property()] functions may use very long regular expression to match brand names of antimicrobial agents. This may fail on some systems.
#' @section Source:
#' World Health Organization (WHO) Collaborating Centre for Drug Statistics Methodology: \url{https://www.whocc.no/atc_ddd_index/}
#'
#' European Commission Public Health PHARMACEUTICALS - COMMUNITY REGISTER: \url{https://ec.europa.eu/health/documents/community-register/html/reg_hum_atc.htm}
#' @aliases ab
#' @return A [character] [vector] with additional class [`ab`]
#' @seealso
#' * [antivirals] for the [data.frame] that is being used to determine ATCs
#' * [av_from_text()] for a function to retrieve antimicrobial drugs from clinical text (from health care records)
#' @inheritSection AMR Reference Data Publicly Available
#' @export
#' @examples
#' # these examples all return "ACI", the ID of aciclovir:
#' as.av("J05AB01")
#' as.av("J 05 AB 01")
#' as.av("Aciclovir")
#' as.av("aciclo")
#' as.av("   aciclo 123")
#' as.av("ACICL")
#' as.av("ACI")
#' as.av("Virorax") # trade name
#' as.av("Zovirax") # trade name
#'
#' as.av("acyklofir") # severe spelling error, yet works
#'
#' # use av_* functions to get a specific properties (see ?av_property);
#' # they use as.av() internally:
#' av_name("J05AB01")
#' av_name("acicl")
as.av <- function(x, flag_multiple_results = TRUE, info = interactive(), ...) {
  meet_criteria(x, allow_class = c("character", "numeric", "integer", "factor"), allow_NA = TRUE)
  meet_criteria(flag_multiple_results, allow_class = "logical", has_length = 1)
  meet_criteria(info, allow_class = "logical", has_length = 1)

  if (is.av(x)) {
    return(x)
  }
  if (all(x %in% c(AMR_env$AV_lookup$av, NA))) {
    # all valid AB codes, but not yet right class
    return(set_clean_class(x,
      new_class = c("av", "character")
    ))
  }

  initial_search <- is.null(list(...)$initial_search)
  already_regex <- isTRUE(list(...)$already_regex)
  fast_mode <- isTRUE(list(...)$fast_mode)

  x_bak <- x
  x <- toupper(x)

  # remove diacritics
  x <- iconv(x, from = "UTF-8", to = "ASCII//TRANSLIT")
  x <- gsub('"', "", x, fixed = TRUE)
  x <- gsub("(specimen|specimen date|specimen_date|spec_date|gender|^dates?$)", "", x, ignore.case = TRUE, perl = TRUE)
  x_bak_clean <- x
  if (already_regex == FALSE) {
    x_bak_clean <- generalise_antibiotic_name(x_bak_clean)
  }

  x <- unique(x_bak_clean) # this means that every x is in fact generalise_antibiotic_name(x)
  x_new <- rep(NA_character_, length(x))
  x_unknown <- character(0)
  x_unknown_ATCs <- character(0)

  note_if_more_than_one_found <- function(found, index, from_text) {
    if (initial_search == TRUE && isTRUE(length(from_text) > 1)) {
      avnames <- av_name(from_text, tolower = TRUE, initial_search = FALSE)
      if (av_name(found[1L], language = NULL) %like% "(clavulanic acid|avibactam)") {
        avnames <- avnames[!avnames %in% c("clavulanic acid", "avibactam")]
      }
      if (length(avnames) > 1) {
        warning_(
          "More than one result was found for item ", index, ": ",
          vector_and(avnames, quotes = FALSE)
        )
      }
    }
    found[1L]
  }

  # Fill in names, AB codes, CID codes and ATC codes directly (`x` is already clean and uppercase)
  known_names <- x %in% AMR_env$AV_lookup$generalised_name
  x_new[known_names] <- AMR_env$AV_lookup$av[match(x[known_names], AMR_env$AV_lookup$generalised_name)]
  known_codes_av <- x %in% AMR_env$AV_lookup$av
  known_codes_atc <- vapply(FUN.VALUE = logical(1), x, function(x_) x_ %in% unlist(AMR_env$AV_lookup$atc), USE.NAMES = FALSE)
  known_codes_cid <- x %in% AMR_env$AV_lookup$cid
  x_new[known_codes_av] <- AMR_env$AV_lookup$av[match(x[known_codes_av], AMR_env$AV_lookup$av)]
  x_new[known_codes_atc] <- AMR_env$AV_lookup$av[vapply(
    FUN.VALUE = integer(1),
    x[known_codes_atc],
    function(x_) {
      which(vapply(
        FUN.VALUE = logical(1),
        AMR_env$AV_lookup$atc,
        function(atc) x_ %in% atc
      ))[1L]
    },
    USE.NAMES = FALSE
  )]
  x_new[known_codes_cid] <- AMR_env$AV_lookup$av[match(x[known_codes_cid], AMR_env$AV_lookup$cid)]
  previously_coerced <- x %in% AMR_env$av_previously_coerced$x
  x_new[previously_coerced & is.na(x_new)] <- AMR_env$av_previously_coerced$av[match(x[is.na(x_new) & x %in% AMR_env$av_previously_coerced$x], AMR_env$av_previously_coerced$x)]
  already_known <- known_names | known_codes_av | known_codes_atc | known_codes_cid | previously_coerced

  # fix for NAs
  x_new[is.na(x)] <- NA
  already_known[is.na(x)] <- FALSE

  if (initial_search == TRUE && sum(already_known) < length(x)) {
    progress <- progress_ticker(n = sum(!already_known), n_min = 25, print = info) # start if n >= 25
    on.exit(close(progress))
  }

  for (i in which(!already_known)) {
    if (initial_search == TRUE) {
      progress$tick()
    }

    if (is.na(x[i]) || is.null(x[i])) {
      next
    }
    if (identical(x[i], "") ||
      # prevent "bacteria" from coercing to TMP, since Bacterial is a brand name of it:
      identical(tolower(x[i]), "bacteria")) {
      x_unknown <- c(x_unknown, x_bak[x[i] == x_bak_clean][1])
      next
    }
    if (x[i] %like_case% "[A-Z][0-9][0-9][A-Z][A-Z][0-9][0-9]") {
      # seems an ATC code, but the available ones are in `already_known`, so:
      x_unknown <- c(x_unknown, x[i])
      x_unknown_ATCs <- c(x_unknown_ATCs, x[i])
      x_new[i] <- NA_character_
      next
    }

    if (fast_mode == FALSE && flag_multiple_results == TRUE && x[i] %like% "[ ]") {
      from_text <- tryCatch(suppressWarnings(av_from_text(x[i], initial_search = FALSE, translate_av = FALSE)[[1]]),
        error = function(e) character(0)
      )
    } else {
      from_text <- character(0)
    }

    # old code for phenoxymethylpenicillin (Peni V)
    if (x[i] == "PNV") {
      x_new[i] <- "PHN"
      next
    }

    # exact LOINC code
    loinc_found <- unlist(lapply(
      AMR_env$AV_lookup$generalised_loinc,
      function(s) x[i] %in% s
    ))
    found <- AMR_env$AV_lookup$av[loinc_found == TRUE]
    if (length(found) > 0) {
      x_new[i] <- note_if_more_than_one_found(found, i, from_text)
      next
    }

    # exact synonym
    synonym_found <- unlist(lapply(
      AMR_env$AV_lookup$generalised_synonyms,
      function(s) x[i] %in% s
    ))
    found <- AMR_env$AV_lookup$av[synonym_found == TRUE]
    if (length(found) > 0) {
      x_new[i] <- note_if_more_than_one_found(found, i, from_text)
      next
    }

    # length of input is quite long, and Levenshtein distance is only max 2
    if (nchar(x[i]) >= 10) {
      levenshtein <- as.double(utils::adist(x[i], AMR_env$AV_lookup$generalised_name))
      if (any(levenshtein <= 2)) {
        found <- AMR_env$AV_lookup$av[which(levenshtein <= 2)]
        x_new[i] <- note_if_more_than_one_found(found, i, from_text)
        next
      }
    }

    # allow characters that resemble others, but only continue when having more than 3 characters
    if (nchar(x[i]) <= 3) {
      x_unknown <- c(x_unknown, x_bak[x[i] == x_bak_clean][1])
      next
    }
    x_spelling <- x[i]
    if (already_regex == FALSE) {
      x_spelling <- gsub("[IY]+", "[IY]+", x_spelling, perl = TRUE)
      x_spelling <- gsub("(C|K|Q|QU|S|Z|X|KS)+", "(C|K|Q|QU|S|Z|X|KS)+", x_spelling, perl = TRUE)
      x_spelling <- gsub("(PH|F|V)+", "(PH|F|V)+", x_spelling, perl = TRUE)
      x_spelling <- gsub("(TH|T)+", "(TH|T)+", x_spelling, perl = TRUE)
      x_spelling <- gsub("A+", "A+", x_spelling, perl = TRUE)
      x_spelling <- gsub("E+", "E+", x_spelling, perl = TRUE)
      x_spelling <- gsub("O+", "O+", x_spelling, perl = TRUE)
      # allow any ending of -in/-ine and -im/-ime
      x_spelling <- gsub("(\\[IY\\]\\+(N|M)|\\[IY\\]\\+(N|M)E\\+?)$", "[IY]+(N|M)E*", x_spelling, perl = TRUE)
      # allow any ending of -ol/-ole
      x_spelling <- gsub("(O\\+L|O\\+LE\\+)$", "O+LE*", x_spelling, perl = TRUE)
      # allow any ending of -on/-one
      x_spelling <- gsub("(O\\+N|O\\+NE\\+)$", "O+NE*", x_spelling, perl = TRUE)
      # replace multiple same characters to single one with '+', like "ll" -> "l+"
      x_spelling <- gsub("(.)\\1+", "\\1+", x_spelling, perl = TRUE)
      # replace spaces and slashes with a possibility on both
      x_spelling <- gsub("[ /]", "( .*|.*/)", x_spelling, perl = TRUE)
      # correct for digital reading text (OCR)
      x_spelling <- gsub("[NRD8B]", "[NRD8B]", x_spelling, perl = TRUE)
      x_spelling <- gsub("(O|0)", "(O|0)+", x_spelling, perl = TRUE)
      x_spelling <- gsub("++", "+", x_spelling, fixed = TRUE)
    }

    # try if name starts with it
    found <- AMR_env$AV_lookup[which(AMR_env$AV_lookup$generalised_name %like% paste0("^", x_spelling)), "av", drop = TRUE]
    if (length(found) > 0) {
      x_new[i] <- note_if_more_than_one_found(found, i, from_text)
      next
    }
    # try if name ends with it
    found <- AMR_env$AV_lookup[which(AMR_env$AV_lookup$generalised_name %like% paste0(x_spelling, "$")), "av", drop = TRUE]
    if (nchar(x[i]) >= 4 && length(found) > 0) {
      x_new[i] <- note_if_more_than_one_found(found, i, from_text)
      next
    }

    # and try if any synonym starts with it
    synonym_found <- unlist(lapply(
      AMR_env$AV_lookup$generalised_synonyms,
      function(s) any(s %like% paste0("^", x_spelling))
    ))
    found <- AMR_env$AV_lookup$av[synonym_found == TRUE]
    if (length(found) > 0) {
      x_new[i] <- note_if_more_than_one_found(found, i, from_text)
      next
    }

    # INITIAL SEARCH - More uncertain results ----

    if (initial_search == TRUE && fast_mode == FALSE) {
      # only run on first try

      # try by removing all spaces
      if (x[i] %like% " ") {
        found <- suppressWarnings(as.av(gsub(" +", "", x[i], perl = TRUE), initial_search = FALSE))
        if (length(found) > 0 && !is.na(found)) {
          x_new[i] <- note_if_more_than_one_found(found, i, from_text)
          next
        }
      }

      # try by removing all spaces and numbers
      if (x[i] %like% " " || x[i] %like% "[0-9]") {
        found <- suppressWarnings(as.av(gsub("[ 0-9]", "", x[i], perl = TRUE), initial_search = FALSE))
        if (length(found) > 0 && !is.na(found)) {
          x_new[i] <- note_if_more_than_one_found(found, i, from_text)
          next
        }
      }

      # transform back from other languages and try again
      x_translated <- paste(lapply(
        strsplit(x[i], "[^A-Z0-9]"),
        function(y) {
          for (i in seq_len(length(y))) {
            for (lang in LANGUAGES_SUPPORTED[LANGUAGES_SUPPORTED != "en"]) {
              y[i] <- ifelse(tolower(y[i]) %in% tolower(TRANSLATIONS[, lang, drop = TRUE]),
                TRANSLATIONS[which(tolower(TRANSLATIONS[, lang, drop = TRUE]) == tolower(y[i]) &
                  !isFALSE(TRANSLATIONS$fixed)), "pattern"],
                y[i]
              )
            }
          }
          generalise_antibiotic_name(y)
        }
      )[[1]],
      collapse = "/"
      )
      x_translated_guess <- suppressWarnings(as.av(x_translated, initial_search = FALSE))
      if (!is.na(x_translated_guess)) {
        x_new[i] <- x_translated_guess
        next
      }

      # now also try to coerce brandname combinations like "Amoxy/clavulanic acid"
      x_translated <- paste(lapply(
        strsplit(x_translated, "[^A-Z0-9 ]"),
        function(y) {
          for (i in seq_len(length(y))) {
            y_name <- suppressWarnings(av_name(y[i], language = NULL, initial_search = FALSE))
            y[i] <- ifelse(!is.na(y_name),
              y_name,
              y[i]
            )
          }
          generalise_antibiotic_name(y)
        }
      )[[1]],
      collapse = "/"
      )
      x_translated_guess <- suppressWarnings(as.av(x_translated, initial_search = FALSE))
      if (!is.na(x_translated_guess)) {
        x_new[i] <- x_translated_guess
        next
      }

      # try by removing all trailing capitals
      if (x[i] %like_case% "[a-z]+[A-Z]+$") {
        found <- suppressWarnings(as.av(gsub("[A-Z]+$", "", x[i], perl = TRUE), initial_search = FALSE))
        if (!is.na(found)) {
          x_new[i] <- note_if_more_than_one_found(found, i, from_text)
          next
        }
      }

      # keep only letters
      found <- suppressWarnings(as.av(gsub("[^A-Z]", "", x[i], perl = TRUE), initial_search = FALSE))
      if (!is.na(found)) {
        x_new[i] <- note_if_more_than_one_found(found, i, from_text)
        next
      }

      # try from a bigger text, like from a health care record, see ?av_from_text
      # already calculated above if flag_multiple_results = TRUE
      if (flag_multiple_results == TRUE) {
        found <- from_text[1L]
      } else {
        found <- tryCatch(suppressWarnings(av_from_text(x[i], initial_search = FALSE, translate_av = FALSE)[[1]][1L]),
          error = function(e) NA_character_
        )
      }
      if (!is.na(found)) {
        x_new[i] <- note_if_more_than_one_found(found, i, from_text)
        next
      }

      # first 5
      found <- suppressWarnings(as.av(substr(x[i], 1, 5), initial_search = FALSE))
      if (!is.na(found)) {
        x_new[i] <- note_if_more_than_one_found(found, i, from_text)
        next
      }

      # make all consonants facultative
      search_str <- gsub("([BCDFGHJKLMNPQRSTVWXZ])", "\\1*", x[i], perl = TRUE)
      found <- suppressWarnings(as.av(search_str, initial_search = FALSE, already_regex = TRUE))
      # keep at least 4 normal characters
      if (nchar(gsub(".\\*", "", search_str, perl = TRUE)) < 4) {
        found <- NA
      }
      if (!is.na(found)) {
        x_new[i] <- note_if_more_than_one_found(found, i, from_text)
        next
      }

      # make all vowels facultative
      search_str <- gsub("([AEIOUY])", "\\1*", x[i], perl = TRUE)
      found <- suppressWarnings(as.av(search_str, initial_search = FALSE, already_regex = TRUE))
      # keep at least 5 normal characters
      if (nchar(gsub(".\\*", "", search_str, perl = TRUE)) < 5) {
        found <- NA
      }
      if (!is.na(found)) {
        x_new[i] <- note_if_more_than_one_found(found, i, from_text)
        next
      }

      # allow misspelling of vowels
      x_spelling <- gsub("A+", "[AEIOU]+", x_spelling, fixed = TRUE)
      x_spelling <- gsub("E+", "[AEIOU]+", x_spelling, fixed = TRUE)
      x_spelling <- gsub("I+", "[AEIOU]+", x_spelling, fixed = TRUE)
      x_spelling <- gsub("O+", "[AEIOU]+", x_spelling, fixed = TRUE)
      x_spelling <- gsub("U+", "[AEIOU]+", x_spelling, fixed = TRUE)
      found <- suppressWarnings(as.av(x_spelling, initial_search = FALSE, already_regex = TRUE))
      if (!is.na(found)) {
        x_new[i] <- note_if_more_than_one_found(found, i, from_text)
        next
      }

      # try with switched character, like "mreopenem"
      for (j in seq_len(nchar(x[i]))) {
        x_switched <- paste0(
          # beginning part:
          substr(x[i], 1, j - 1),
          # here is the switching of 2 characters:
          substr(x[i], j + 1, j + 1),
          substr(x[i], j, j),
          # ending part:
          substr(x[i], j + 2, nchar(x[i]))
        )
        found <- suppressWarnings(as.av(x_switched, initial_search = FALSE))
        if (!is.na(found)) {
          break
        }
      }
      if (!is.na(found)) {
        x_new[i] <- found[1L]
        next
      }
    } # end of initial_search = TRUE

    # not found
    x_unknown <- c(x_unknown, x_bak[x[i] == x_bak_clean][1])
  }

  if (initial_search == TRUE && sum(already_known) < length(x)) {
    close(progress)
  }

  # save to package env to save time for next time
  if (initial_search == TRUE) {
    AMR_env$av_previously_coerced <- AMR_env$av_previously_coerced[which(!AMR_env$av_previously_coerced$x %in% x), , drop = FALSE]
    AMR_env$av_previously_coerced <- unique(rbind(AMR_env$av_previously_coerced,
      data.frame(
        x = x,
        av = x_new,
        x_bak = x_bak[match(x, x_bak_clean)],
        stringsAsFactors = FALSE
      ),
      stringsAsFactors = FALSE
    ))
  }

  # take failed ATC codes apart from rest
  if (length(x_unknown_ATCs) > 0 && fast_mode == FALSE) {
    warning_(
      "in `as.av()`: these ATC codes are not (yet) in the antivirals data set: ",
      vector_and(x_unknown_ATCs), "."
    )
  }
  x_unknown <- x_unknown[!x_unknown %in% x_unknown_ATCs]
  x_unknown <- c(x_unknown,
                 AMR_env$av_previously_coerced$x_bak[which(AMR_env$av_previously_coerced$x %in% x & is.na(AMR_env$av_previously_coerced$av))])
  if (length(x_unknown) > 0 && fast_mode == FALSE) {
    warning_(
      "in `as.av()`: these values could not be coerced to a valid antiviral agent ID: ",
      vector_and(x_unknown), "."
    )
  }

  x_result <- x_new[match(x_bak_clean, x)]
  if (length(x_result) == 0) {
    x_result <- NA_character_
  }

  set_clean_class(x_result,
    new_class = c("av", "character")
  )
}

#' @rdname as.av
#' @export
is.av <- function(x) {
  inherits(x, "av")
}

# will be exported using s3_register() in R/zzz.R
pillar_shaft.av <- function(x, ...) {
  out <- trimws(format(x))
  out[is.na(x)] <- font_na(NA)
  create_pillar_column(out, align = "left", min_width = 4)
}

# will be exported using s3_register() in R/zzz.R
type_sum.av <- function(x, ...) {
  "av"
}

#' @method print av
#' @export
#' @noRd
print.av <- function(x, ...) {
  cat("Class 'av'\n")
  print(as.character(x), quote = FALSE)
}

#' @method as.data.frame av
#' @export
#' @noRd
as.data.frame.av <- function(x, ...) {
  nm <- deparse1(substitute(x))
  if (!"nm" %in% names(list(...))) {
    as.data.frame.vector(as.av(x), ..., nm = nm)
  } else {
    as.data.frame.vector(as.av(x), ...)
  }
}
#' @method [ av
#' @export
#' @noRd
"[.av" <- function(x, ...) {
  y <- NextMethod()
  attributes(y) <- attributes(x)
  y
}
#' @method [[ av
#' @export
#' @noRd
"[[.av" <- function(x, ...) {
  y <- NextMethod()
  attributes(y) <- attributes(x)
  y
}
#' @method [<- av
#' @export
#' @noRd
"[<-.av" <- function(i, j, ..., value) {
  y <- NextMethod()
  attributes(y) <- attributes(i)
  return_after_integrity_check(y, "antimicrobial code", AMR_env$AV_lookup$av)
}
#' @method [[<- av
#' @export
#' @noRd
"[[<-.av" <- function(i, j, ..., value) {
  y <- NextMethod()
  attributes(y) <- attributes(i)
  return_after_integrity_check(y, "antimicrobial code", AMR_env$AV_lookup$av)
}
#' @method c av
#' @export
#' @noRd
c.av <- function(...) {
  x <- list(...)[[1L]]
  y <- NextMethod()
  attributes(y) <- attributes(x)
  return_after_integrity_check(y, "antimicrobial code", AMR_env$AV_lookup$av)
}

#' @method unique av
#' @export
#' @noRd
unique.av <- function(x, incomparables = FALSE, ...) {
  y <- NextMethod()
  attributes(y) <- attributes(x)
  y
}

#' @method rep av
#' @export
#' @noRd
rep.av <- function(x, ...) {
  y <- NextMethod()
  attributes(y) <- attributes(x)
  y
}

get_translate_av <- function(translate_av) {
  translate_av <- as.character(translate_av)[1L]
  if (translate_av %in% c("TRUE", "official")) {
    return("name")
  } else if (translate_av %in% c(NA_character_, "FALSE")) {
    return(FALSE)
  } else {
    translate_av <- tolower(translate_av)
    stop_ifnot(translate_av %in% colnames(AMR::antivirals),
               "invalid value for 'translate_av', this must be a column name of the antivirals data set\n",
               "or TRUE (equals 'name') or FALSE to not translate at all.",
               call = FALSE
    )
    translate_av
  }
}
