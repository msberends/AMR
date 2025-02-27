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

#' Transform Input to an Antibiotic ID
#'
#' Use this function to determine the antibiotic drug code of one or more antibiotics. The data set [antibiotics] will be searched for abbreviations, official names and synonyms (brand names).
#' @param x a [character] vector to determine to antibiotic ID
#' @param flag_multiple_results a [logical] to indicate whether a note should be printed to the console that probably more than one antibiotic drug code or name can be retrieved from a single input value.
#' @param language language to coerce input values from any of the `r length(LANGUAGES_SUPPORTED)` supported languages - default to the system language if supported (see [get_AMR_locale()])
#' @param info a [logical] to indicate whether a progress bar should be printed - the default is `TRUE` only in interactive mode
#' @param ... arguments passed on to internal functions
#' @rdname as.ab
#' @inheritSection WHOCC WHOCC
#' @details All entries in the [antibiotics] data set have three different identifiers: a human readable EARS-Net code (column `ab`, used by ECDC and WHONET), an ATC code (column `atc`, used by WHO), and a CID code (column `cid`, Compound ID, used by PubChem). The data set contains more than 5,000 official brand names from many different countries, as found in PubChem. Not that some drugs contain multiple ATC codes.
#'
#' All these properties will be searched for the user input. The [as.ab()] can correct for different forms of misspelling:
#'
#'  * Wrong spelling of drug names (such as "tobramicin" or "gentamycin"), which corrects for most audible similarities such as f/ph, x/ks, c/z/s, t/th, etc.
#'  * Too few or too many vowels or consonants
#'  * Switching two characters (such as "mreopenem", often the case in clinical data, when doctors typed too fast)
#'  * Digitalised paper records, leaving artefacts like 0/o/O (zero and O's), B/8, n/r, etc.
#'
#' Use the [`ab_*`][ab_property()] functions to get properties based on the returned antibiotic ID, see *Examples*.
#'
#' Note: the [as.ab()] and [`ab_*`][ab_property()] functions may use very long regular expression to match brand names of antimicrobial drugs. This may fail on some systems.
#'
#' You can add your own manual codes to be considered by [as.ab()] and all [`ab_*`][ab_property()] functions, see [add_custom_antimicrobials()].
#' @section Source:
#' World Health Organization (WHO) Collaborating Centre for Drug Statistics Methodology: \url{https://atcddd.fhi.no/atc_ddd_index/}
#'
#' European Commission Public Health PHARMACEUTICALS - COMMUNITY REGISTER: \url{https://ec.europa.eu/health/documents/community-register/html/reg_hum_atc.htm}
#' @aliases ab
#' @return A [character] [vector] with additional class [`ab`]
#' @seealso
#' * [antibiotics] for the [data.frame] that is being used to determine ATCs
#' * [ab_from_text()] for a function to retrieve antimicrobial drugs from clinical text (from health care records)
#' @inheritSection AMR Reference Data Publicly Available
#' @export
#' @examples
#' # these examples all return "ERY", the ID of erythromycin:
#' as.ab("J01FA01")
#' as.ab("J 01 FA 01")
#' as.ab("Erythromycin")
#' as.ab("eryt")
#' as.ab("ERYT")
#' as.ab("ERY")
#' as.ab("eritromicine") # spelled wrong, yet works
#' as.ab("Erythrocin") # trade name
#' as.ab("Romycin") # trade name
#'
#' # spelling from different languages and dyslexia are no problem
#' ab_atc("ceftriaxon")
#' ab_atc("cephtriaxone") # small spelling error
#' ab_atc("cephthriaxone") # or a bit more severe
#' ab_atc("seephthriaaksone") # and even this works
#'
#' # use ab_* functions to get a specific properties (see ?ab_property);
#' # they use as.ab() internally:
#' ab_name("J01FA01")
#' ab_name("eryt")
#'
#' \donttest{
#' if (require("dplyr")) {
#'   # you can quickly rename 'sir' columns using set_ab_names() with dplyr:
#'   example_isolates %>%
#'     set_ab_names(where(is.sir), property = "atc")
#' }
#' }
as.ab <- function(x, flag_multiple_results = TRUE, language = get_AMR_locale(), info = interactive(), ...) {
  meet_criteria(x, allow_class = c("character", "numeric", "integer", "factor"), allow_NA = TRUE)
  meet_criteria(flag_multiple_results, allow_class = "logical", has_length = 1)
  language <- validate_language(language)
  meet_criteria(info, allow_class = "logical", has_length = 1)

  if (is.ab(x) || all(x %in% c(AMR_env$AB_lookup$ab, NA))) {
    # all valid AB codes, but not yet right class or might have additional attributes as AMR selector
    attributes(x) <- NULL
    return(set_clean_class(x,
      new_class = c("ab", "character")
    ))
  }

  already_regex <- isTRUE(list(...)$already_regex)
  fast_mode <- isTRUE(list(...)$fast_mode)

  x_bak <- x
  x <- toupper(x)

  # remove diacritics
  x <- iconv(x, from = "UTF-8", to = "ASCII//TRANSLIT")
  x <- gsub('"', "", x, fixed = TRUE)
  x <- gsub("(specimen|specimen date|specimen_date|spec_date|gender|^dates?$|animal|host($|[a-z]))", "", x, ignore.case = TRUE, perl = TRUE)
  # penicillin is a special case: we call it so, but then most often mean benzylpenicillin
  x[x %like_case% "^PENICILLIN" & x %unlike_case% "[ /+-]"] <- "benzylpenicillin"
  x_bak_clean <- x
  if (already_regex == FALSE) {
    x_bak_clean <- generalise_antibiotic_name(x_bak_clean)
  }

  x <- unique(x_bak_clean) # this means that every x is in fact generalise_antibiotic_name(x)
  x_new <- rep(NA_character_, length(x))
  x_uncertain <- character(0)
  x_unknown <- character(0)
  x_unknown_ATCs <- character(0)

  note_if_more_than_one_found <- function(found, index, from_text) {
    if (isTRUE(length(from_text) > 1)) {
      abnames <- ab_name(from_text, tolower = TRUE)
      if (ab_name(found[1L], language = NULL) %like% "(clavulanic acid|(avi|tazo|mono|vabor)bactam)") {
        abnames <- abnames[!abnames %in% c("clavulanic acid", "avibactam", "tazobactam", "vaborbactam", "monobactam")]
      }
      if (length(abnames) > 1) {
        if (toupper(paste(abnames, collapse = " ")) %in% AMR_env$AB_lookup$generalised_name) {
          # if the found values combined is a valid AB, return that
          found <- AMR_env$AB_lookup$ab[match(toupper(paste(abnames, collapse = " ")), AMR_env$AB_lookup$generalised_name)][1]
        } else {
          message_(
            "More than one result was found for item ", index, ": ",
            vector_and(abnames, quotes = FALSE)
          )
        }
      }
    }
    found[1L]
  }

  # Fill in names, AB codes, CID codes and ATC codes directly (`x` is already clean and uppercase)
  known_names <- x %in% AMR_env$AB_lookup$generalised_name
  x_new[known_names] <- AMR_env$AB_lookup$ab[match(x[known_names], AMR_env$AB_lookup$generalised_name)]
  known_codes_ab <- x %in% AMR_env$AB_lookup$ab
  known_codes_atc <- vapply(FUN.VALUE = logical(1), gsub(" ", "", x), function(x_) x_ %in% unlist(AMR_env$AB_lookup$atc), USE.NAMES = FALSE)
  known_codes_cid <- x %in% AMR_env$AB_lookup$cid
  x_new[known_codes_ab] <- AMR_env$AB_lookup$ab[match(x[known_codes_ab], AMR_env$AB_lookup$ab)]
  x_new[known_codes_atc] <- AMR_env$AB_lookup$ab[vapply(
    FUN.VALUE = integer(1),
    gsub(" ", "", x[known_codes_atc]),
    function(x_) {
      which(vapply(
        FUN.VALUE = logical(1),
        AMR_env$AB_lookup$atc,
        function(atc) x_ %in% atc
      ))[1L]
    },
    USE.NAMES = FALSE
  )]
  x_new[known_codes_cid] <- AMR_env$AB_lookup$ab[match(x[known_codes_cid], AMR_env$AB_lookup$cid)]
  previously_coerced <- x %in% AMR_env$ab_previously_coerced$x
  x_new[previously_coerced & is.na(x_new)] <- AMR_env$ab_previously_coerced$ab[match(x[is.na(x_new) & x %in% AMR_env$ab_previously_coerced$x], AMR_env$ab_previously_coerced$x)]
  prev <- x_bak[which(x[which(previously_coerced)] %in% x_bak_clean)]
  if (any(previously_coerced) && isTRUE(info) && message_not_thrown_before("as.ab", prev, entire_session = TRUE)) {
    message_(
      "Returning previously coerced value", ifelse(length(unique(prev)) > 1, "s", ""),
      " for ", vector_and(prev), ". Run `ab_reset_session()` to reset this. This note will be shown once per session for this input."
    )
  }

  already_known <- known_names | known_codes_ab | known_codes_atc | known_codes_cid | previously_coerced

  # fix for NAs
  x_new[is.na(x)] <- NA
  already_known[is.na(x)] <- FALSE

  if (sum(already_known) < length(x)) {
    progress <- progress_ticker(n = sum(!already_known), n_min = 25, print = info) # start if n >= 25
    on.exit(close(progress))
  }

  for (i in which(!already_known)) {
    progress$tick()

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
      from_text <- tryCatch(suppressWarnings(ab_from_text(x[i], translate_ab = FALSE)[[1]]),
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
      AMR_env$AB_lookup$generalised_loinc,
      function(s) x[i] %in% s
    ))
    found <- AMR_env$AB_lookup$ab[loinc_found == TRUE]
    if (length(found) > 0) {
      x_new[i] <- note_if_more_than_one_found(found, i, from_text)
      next
    }

    # exact synonym
    synonym_found <- unlist(lapply(
      AMR_env$AB_lookup$generalised_synonyms,
      function(s) x[i] %in% s
    ))
    found <- AMR_env$AB_lookup$ab[synonym_found == TRUE]
    if (length(found) > 0) {
      x_new[i] <- note_if_more_than_one_found(found, i, from_text)
      next
    }

    # exact abbreviation
    abbr_found <- unlist(lapply(
      AMR_env$AB_lookup$generalised_abbreviations,
      # require at least 2 characters for abbreviations
      function(s) x[i] %in% s && nchar(x[i]) >= 2
    ))
    found <- AMR_env$AB_lookup$ab[abbr_found == TRUE]
    if (length(found) > 0) {
      x_new[i] <- note_if_more_than_one_found(found, i, from_text)
      next
    }

    # length of input is quite long, and Levenshtein distance is only max 2
    if (nchar(x[i]) >= 10) {
      levenshtein <- as.double(utils::adist(x[i], AMR_env$AB_lookup$generalised_name))
      if (any(levenshtein <= 2)) {
        found <- AMR_env$AB_lookup$ab[which(levenshtein <= 2)]
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
    found <- AMR_env$AB_lookup[which(AMR_env$AB_lookup$generalised_name %like% paste0("^", x_spelling)), "ab", drop = TRUE]
    if (length(found) > 0) {
      x_new[i] <- note_if_more_than_one_found(found, i, from_text)
      next
    }


    # try if name ends with it
    found <- AMR_env$AB_lookup[which(AMR_env$AB_lookup$generalised_name %like% paste0(x_spelling, "$")), "ab", drop = TRUE]
    if (nchar(x[i]) >= 4 && length(found) > 0) {
      x_new[i] <- note_if_more_than_one_found(found, i, from_text)
      next
    }

    # and try if any synonym starts with it
    synonym_found <- unlist(lapply(
      AMR_env$AB_lookup$generalised_synonyms,
      function(s) any(s %like% paste0("^", x_spelling))
    ))
    found <- AMR_env$AB_lookup$ab[synonym_found == TRUE]
    if (length(found) > 0) {
      x_new[i] <- note_if_more_than_one_found(found, i, from_text)
      next
    }

    # More uncertain results ----
    if (fast_mode == FALSE) {
      ab_df <- AMR_env$AB_lookup
      ab_df$length_name <- nchar(ab_df$generalised_name)
      # now retrieve Levensthein distance for name, synonyms, and translated names
      ab_df$lev_name <- as.double(utils::adist(x[i], ab_df$generalised_name,
        ignore.case = FALSE,
        fixed = TRUE,
        costs = c(insertions = 1, deletions = 1, substitutions = 2),
        counts = FALSE
      ))
      ab_df$lev_syn <- vapply(
        FUN.VALUE = double(1),
        ab_df$generalised_synonyms,
        function(y) {
          ifelse(length(y[nchar(y) >= 5]) == 0,
            999,
            min(as.double(utils::adist(x[i], y[nchar(y) >= 5],
              ignore.case = FALSE,
              fixed = TRUE,
              costs = c(insertions = 1, deletions = 1, substitutions = 2),
              counts = FALSE
            )), na.rm = TRUE)
          )
        },
        USE.NAMES = FALSE
      )
      if (!is.null(language) && language != "en") {
        ab_df$trans <- generalise_antibiotic_name(translate_AMR(ab_df$name, language = language))
        ab_df$lev_trans <- as.double(utils::adist(x[i], ab_df$trans,
          ignore.case = FALSE,
          fixed = TRUE,
          costs = c(insertions = 1, deletions = 1, substitutions = 2),
          counts = FALSE
        ))
      } else {
        ab_df$lev_trans <- ab_df$lev_name
      }

      if (any(ab_df$lev_name < 5, na.rm = TRUE)) {
        x_new[i] <- ab_df$ab[order(ab_df$lev_name)][1]
        x_uncertain <- c(x_uncertain, x_bak[x[i] == x_bak_clean][1])
        next
      } else if (any(ab_df$lev_trans < 5, na.rm = TRUE)) {
        x_new[i] <- ab_df$ab[order(ab_df$lev_trans)][1]
        x_uncertain <- c(x_uncertain, x_bak[x[i] == x_bak_clean][1])
        next
      } else if (any(ab_df$lev_syn < 5, na.rm = TRUE)) {
        x_new[i] <- ab_df$ab[order(ab_df$lev_syn)][1]
        x_uncertain <- c(x_uncertain, x_bak[x[i] == x_bak_clean][1])
        next
      } else {
        # then just take name if Levensthein is max 100% of length of name
        ab_df$lev_len_ratio <- ab_df$lev_name / ab_df$length_name
        if (any(ab_df$lev_len_ratio < 1)) {
          ab_df <- ab_df[ab_df$lev_len_ratio < 1, , drop = FALSE]
          x_new[i] <- ab_df$ab[order(ab_df$lev_name)][1]
          x_uncertain <- c(x_uncertain, x_bak[x[i] == x_bak_clean][1])
          next
        }
      }
    }

    # nothing found
    x_unknown <- c(x_unknown, x_bak[x[i] == x_bak_clean][1])
  }

  if (sum(already_known) < length(x)) {
    close(progress)
  }

  # save to package env to save time for next time
  AMR_env$ab_previously_coerced <- AMR_env$ab_previously_coerced[which(!AMR_env$ab_previously_coerced$x %in% x), , drop = FALSE]
  AMR_env$ab_previously_coerced <- unique(rbind_AMR(
    AMR_env$ab_previously_coerced,
    data.frame(
      x = x,
      ab = x_new,
      x_bak = x_bak[match(x, x_bak_clean)],
      stringsAsFactors = FALSE
    )
  ))

  # take failed ATC codes apart from rest
  if (length(x_unknown_ATCs) > 0 && fast_mode == FALSE) {
    warning_(
      "in `as.ab()`: these ATC codes are not (yet) in the antibiotics data set: ",
      vector_and(x_unknown_ATCs), "."
    )
  }

  # Throw note about uncertainties
  x_unknown <- x_unknown[!x_unknown %in% x_unknown_ATCs]
  x_unknown <- c(
    x_unknown,
    AMR_env$ab_previously_coerced$x_bak[which(AMR_env$ab_previously_coerced$x %in% x & is.na(AMR_env$ab_previously_coerced$ab))]
  )
  x_unknown <- x_unknown[!x_unknown %in% c("", NA)]
  if (length(x_unknown) > 0 && fast_mode == FALSE) {
    warning_(
      "in `as.ab()`: these values could not be coerced to a valid antimicrobial ID: ",
      vector_and(x_unknown), "."
    )
  }

  # Throw note about uncertainties
  if (isTRUE(info) && length(x_uncertain) > 0 && fast_mode == FALSE) {
    x_uncertain <- unique(x_uncertain)
    if (message_not_thrown_before("as.ab", "uncertainties", x_bak)) {
      if (length(x_uncertain) <= 3) {
        examples <- vector_and(
          paste0(
            '"', x_uncertain, '" (assumed ',
            ab_name(AMR_env$ab_previously_coerced$ab[which(AMR_env$ab_previously_coerced$x_bak %in% x_uncertain)], language = NULL, tolower = TRUE),
            ", ", AMR_env$ab_previously_coerced$ab[which(AMR_env$ab_previously_coerced$x_bak %in% x_uncertain)], ")"
          ),
          quotes = FALSE
        )
      } else {
        examples <- paste0(nr2char(length(x_uncertain)), " antimicrobials")
      }
      message_(
        "Antimicrobial translation was uncertain for ", examples,
        ". If required, use `add_custom_antimicrobials()` to add custom entries."
      )
    }
  }

  x_result <- x_new[match(x_bak_clean, x)]
  if (length(x_result) == 0) {
    x_result <- NA_character_
  }

  set_clean_class(x_result,
    new_class = c("ab", "character")
  )
}

#' @rdname as.ab
#' @export
is.ab <- function(x) {
  inherits(x, "ab")
}

#' @rdname as.ab
#' @export
ab_reset_session <- function() {
  if (NROW(AMR_env$ab_previously_coerced) > 0) {
    message_("Reset ", nr2char(NROW(AMR_env$ab_previously_coerced)), " previously matched input value", ifelse(NROW(AMR_env$ab_previously_coerced) > 1, "s", ""), ".")
    AMR_env$ab_previously_coerced <- AMR_env$ab_previously_coerced[0, , drop = FALSE]
    AMR_env$mo_uncertainties <- AMR_env$mo_uncertainties[0, , drop = FALSE]
  } else {
    message_("No previously matched input values to reset.")
  }
}

# will be exported using s3_register() in R/zzz.R
pillar_shaft.ab <- function(x, ...) {
  out <- trimws(format(x))
  out[is.na(x)] <- font_na(NA)

  # add the names to the drugs as mouse-over!
  if (tryCatch(isTRUE(getExportedValue("ansi_has_hyperlink_support", ns = asNamespace("cli"))()), error = function(e) FALSE)) {
    out[!is.na(x)] <- font_url(
      url = paste0(x[!is.na(x)], ": ", ab_name(x[!is.na(x)])),
      txt = out[!is.na(x)]
    )
  }

  create_pillar_column(out, align = "left", min_width = 4)
}

# will be exported using s3_register() in R/zzz.R
type_sum.ab <- function(x, ...) {
  "ab"
}

#' @method print ab
#' @export
#' @noRd
print.ab <- function(x, ...) {
  if (!is.null(attributes(x)$amr_selector)) {
    function_name <- attributes(x)$amr_selector
    message_(
      "This 'ab' vector was retrieved using `", function_name, "()`, which should normally be used inside a `dplyr` verb or `data.frame` call, e.g.:\n",
      "  ", AMR_env$bullet_icon, " your_data %>% select(", function_name, "())\n",
      "  ", AMR_env$bullet_icon, " your_data %>% select(column_a, column_b, ", function_name, "())\n",
      "  ", AMR_env$bullet_icon, " your_data %>% filter(any(", function_name, "() == \"R\"))\n",
      "  ", AMR_env$bullet_icon, " your_data[, ", function_name, "()]\n",
      "  ", AMR_env$bullet_icon, " your_data[, c(\"column_a\", \"column_b\", ", function_name, "())]"
    )
  }
  cat("Class 'ab'\n")
  print(as.character(x), quote = FALSE)
}

#' @method as.data.frame ab
#' @export
#' @noRd
as.data.frame.ab <- function(x, ...) {
  nm <- deparse1(substitute(x))
  if (!"nm" %in% names(list(...))) {
    as.data.frame.vector(as.ab(x), ..., nm = nm)
  } else {
    as.data.frame.vector(as.ab(x), ...)
  }
}
#' @method [ ab
#' @export
#' @noRd
"[.ab" <- function(x, ...) {
  y <- NextMethod()
  attributes(y) <- attributes(x)
  y
}
#' @method [[ ab
#' @export
#' @noRd
"[[.ab" <- function(x, ...) {
  y <- NextMethod()
  attributes(y) <- attributes(x)
  y
}
#' @method [<- ab
#' @export
#' @noRd
"[<-.ab" <- function(i, j, ..., value) {
  y <- NextMethod()
  attributes(y) <- attributes(i)
  return_after_integrity_check(y, "antimicrobial drug code", AMR_env$AB_lookup$ab)
}
#' @method [[<- ab
#' @export
#' @noRd
"[[<-.ab" <- function(i, j, ..., value) {
  y <- NextMethod()
  attributes(y) <- attributes(i)
  return_after_integrity_check(y, "antimicrobial drug code", AMR_env$AB_lookup$ab)
}
#' @method c ab
#' @export
#' @noRd
c.ab <- function(...) {
  x <- list(...)[[1L]]
  y <- NextMethod()
  attributes(y) <- attributes(x)
  return_after_integrity_check(y, "antimicrobial drug code", AMR_env$AB_lookup$ab)
}

#' @method unique ab
#' @export
#' @noRd
unique.ab <- function(x, incomparables = FALSE, ...) {
  y <- NextMethod()
  attributes(y) <- attributes(x)
  y
}

#' @method rep ab
#' @export
#' @noRd
rep.ab <- function(x, ...) {
  y <- NextMethod()
  attributes(y) <- attributes(x)
  y
}

generalise_antibiotic_name <- function(x) {
  x <- toupper(x)
  # remove suffices
  x <- gsub("_(MIC|RSI|SIR|DIS[CK])$", "", x, perl = TRUE)
  # remove disk concentrations, like LVX_NM -> LVX
  x <- gsub("_[A-Z]{2}[0-9_.]{0,3}$", "", x, perl = TRUE)
  # keep only max 1 space
  x <- trimws2(gsub(" +", " ", x, perl = TRUE))
  # non-character, space or number should be a slash
  x <- gsub("[^A-Z0-9 -)(]", "/", x, perl = TRUE)
  # correct for 'high level' antibiotics
  x <- trimws(gsub("([^A-Z0-9/ -]+)?(HIGH(.?LE?VE?L)?|[^A-Z0-9/]H[^A-Z0-9]?L)([^A-Z0-9 -]+)?", "-HIGH", x, perl = TRUE))
  x <- trimws(gsub("^(-HIGH)(.*)", "\\2\\1", x, perl = TRUE))
  # remove part between brackets if that's followed by another string
  x <- gsub("(.*)+ [(].*[)]", "\\1", x)
  # spaces around non-characters must be removed: amox + clav -> amox clav
  x <- gsub("(.*[A-Z0-9]) ([^A-Z0-9].*)", "\\1\\2", x, perl = TRUE)
  x <- gsub("(.*[^A-Z0-9]) ([A-Z0-9].*)", "\\1\\2", x, perl = TRUE)
  # remove hyphen after a starting "co"
  x <- gsub("^CO-", "CO", x, perl = TRUE)
  # replace operators with a space
  x <- gsub("(/| AND | WITH | W/|[+]|[-])+", " ", x, perl = TRUE)
  # replace more than 1 space
  x <- trimws(gsub(" +", " ", x, perl = TRUE))
  # move HIGH to end
  x <- trimws(gsub("(.*) HIGH(.*)", "\\1\\2 HIGH", x, perl = TRUE))
  x
}

get_translate_ab <- function(translate_ab) {
  translate_ab <- as.character(translate_ab)[1L]
  if (translate_ab %in% c("TRUE", "official")) {
    return("name")
  } else if (translate_ab %in% c(NA_character_, "FALSE")) {
    return(FALSE)
  } else {
    translate_ab <- tolower(translate_ab)
    stop_ifnot(translate_ab %in% colnames(AMR::antibiotics),
      "invalid value for 'translate_ab', this must be a column name of the antibiotics data set\n",
      "or TRUE (equals 'name') or FALSE to not translate at all.",
      call = FALSE
    )
    translate_ab
  }
}

create_AB_AV_lookup <- function(df) {
  new_df <- df
  new_df$generalised_name <- generalise_antibiotic_name(new_df$name)
  new_df$generalised_synonyms <- lapply(new_df$synonyms, generalise_antibiotic_name)
  if ("abbreviations" %in% colnames(df)) {
    new_df$generalised_abbreviations <- lapply(new_df$abbreviations, generalise_antibiotic_name)
  }
  new_df$generalised_loinc <- lapply(new_df$loinc, generalise_antibiotic_name)
  new_df$generalised_all <- unname(lapply(
    as.list(as.data.frame(
      t(new_df[,
        c(
          colnames(new_df)[colnames(new_df) %in% c("ab", "av", "atc", "cid", "name")],
          colnames(new_df)[colnames(new_df) %like% "generalised"]
        ),
        drop = FALSE
      ]),
      stringsAsFactors = FALSE
    )),
    function(x) {
      x <- generalise_antibiotic_name(unname(unlist(x)))
      x[x != ""]
    }
  ))
  new_df[, colnames(new_df)[colnames(new_df) %like% "^generalised"]]
}
