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

#' Translate Strings from the AMR Package
#'
#' For language-dependent output of AMR functions, like [mo_name()], [mo_gramstain()], [mo_type()] and [ab_name()].
#' @param x text to translate
#' @param language language to choose. Use one of these supported language names or ISO-639-1 codes: `r paste0('"', sapply(LANGUAGES_SUPPORTED_NAMES, function(x) x[[1]]), '" ("' , LANGUAGES_SUPPORTED, '")', collapse = ", ")`.
#' @details The currently `r length(LANGUAGES_SUPPORTED)` supported languages are `r vector_and(sapply(LANGUAGES_SUPPORTED_NAMES, function(x) x[[1]]), quotes = FALSE, sort = FALSE)`. All these languages have translations available for all antimicrobial agents and colloquial microorganism names.
#'
#' Please read about adding or updating a language in [our Wiki](https://github.com/msberends/AMR/wiki/).
#'
#' ## Changing the Default Language
#' The system language will be used at default (as returned by `Sys.getenv("LANG")` or, if `LANG` is not set, [`Sys.getlocale("LC_COLLATE")`][Sys.getlocale()]), if that language is supported. But the language to be used can be overwritten in two ways and will be checked in this order:
#' 
#'   1. Setting the R option `AMR_locale`, either by using `set_AMR_locale()` or by running e.g. `options(AMR_locale = "de")`.
#'   
#'      Note that setting an \R option only works in the same session. Save the command `options(AMR_locale = "(your language)")` to your `.Rprofile` file to apply it for every session.
#'   2. Setting the system variable `LANGUAGE` or `LANG`, e.g. by adding `LANGUAGE="de_DE.utf8"` to your `.Renviron` file in your home directory.
#' 
#' Thus, if the R option `AMR_locale` is set, the system variables `LANGUAGE` and `LANG` will be ignored.
#' @rdname translate
#' @name translate
#' @export
#' @examples
#' # Current settings (based on system language)
#' ab_name("Ciprofloxacin")
#' mo_name("Coagulase-negative Staphylococcus")
#'
#' # setting another language
#' set_AMR_locale("Greek")
#' ab_name("Ciprofloxacin")
#' mo_name("Coagulase-negative Staphylococcus")
#' 
#' set_AMR_locale("Spanish")
#' ab_name("Ciprofloxacin")
#' mo_name("Coagulase-negative Staphylococcus")
#' 
#' # set_AMR_locale() understands endonyms, English exonyms, and ISO-639-1:
#' set_AMR_locale("Deutsch")
#' set_AMR_locale("German")
#' set_AMR_locale("de")
#'
#' # reset to system default
#' reset_AMR_locale()
get_AMR_locale <- function() {
  if (!is.null(getOption("AMR_locale", default = NULL))) {
    return(validate_language(getOption("AMR_locale"), extra_txt = "set with `options(AMR_locale = ...)`"))
  }
  
  lang <- ""
  # now check the LANGUAGE system variable - return it if set
  if (!identical("", Sys.getenv("LANGUAGE"))) {
    lang <- Sys.getenv("LANGUAGE")
  }
  if (!identical("", Sys.getenv("LANG"))) {
    lang <- Sys.getenv("LANG")
  }
  if (lang == "") {
    lang <- Sys.getlocale("LC_COLLATE")
  }
  
  lang <- find_language(lang)
  if (lang != "en" && interactive() && message_not_thrown_before("get_AMR_locale", entire_session = TRUE)) {
    message_("Assuming the ", LANGUAGES_SUPPORTED_NAMES[[lang]]$exonym,  " language (",
             LANGUAGES_SUPPORTED_NAMES[[lang]]$endonym, ") for the AMR package. Change this with `set_AMR_locale()`. ",
             "This note will be shown once per session.")
  }
  lang
}

#' @rdname translate
#' @export
set_AMR_locale <- function(language) {
  language <- validate_language(language)
  options(AMR_locale = language)
  message_("Using the ", LANGUAGES_SUPPORTED_NAMES[[language]]$exonym,  " language (", LANGUAGES_SUPPORTED_NAMES[[language]]$endonym, ") for the AMR package for this session.")
}

#' @rdname translate
#' @export
reset_AMR_locale <- function() {
  options(AMR_locale = NULL)
}

#' @rdname translate
#' @export
translate_AMR <- function(x, language = get_AMR_locale()) {
  translate_into_language(x, language = language)
}


validate_language <- function(language, extra_txt = character(0)) {
  if (trimws(tolower(language)) %in% c("en", "english", "", "false", NA)) {
    return("en")
  }
  lang <- find_language(language, fallback = FALSE)
  stop_ifnot(length(lang) > 0 && lang %in% LANGUAGES_SUPPORTED,
             "unsupported language for AMR package", extra_txt, ": \"", language, "\". Use one of these language names or ISO-639-1 codes: ",
             paste0('"', vapply(FUN.VALUE = character(1), LANGUAGES_SUPPORTED_NAMES, function(x) x[[1]]),
                    '" ("' , LANGUAGES_SUPPORTED, '")', collapse = ", "),
             call = FALSE)
  lang
}

find_language <- function(language, fallback = TRUE) {
  language <- Map(function(l, n, check = language) {
    grepl(paste0("^(", l[1], "|", l[2], "|",
                 n, "(_|$)|", toupper(n), "(_|$))"),
          check,
          ignore.case = FALSE,
          perl = TRUE,
          useBytes = FALSE)
  },
  LANGUAGES_SUPPORTED_NAMES,
  LANGUAGES_SUPPORTED,
  USE.NAMES = TRUE)
  language <- names(which(language == TRUE))
  if (isTRUE(fallback) && length(language) == 0) {
    # other language -> set to English
    language <- "en"
  }
  language
}

# translate strings based on inst/translations.tsv
translate_into_language <- function(from,
                                    language = get_AMR_locale(), 
                                    only_unknown = FALSE,
                                    only_affect_ab_names = FALSE,
                                    only_affect_mo_names = FALSE) {
  
  if (is.null(language)) {
    return(from)
  }
  if (language %in% c("en", "", NA)) {
    return(from)
  }
  
  df_trans <- TRANSLATIONS # internal data file
  from.bak <- from
  from_unique <- unique(from)
  from_unique_translated <- from_unique
  
  # get ISO-639-1 of language
  lang <- validate_language(language)
  # only keep lines where translation is available for this language
  df_trans <- df_trans[which(!is.na(df_trans[, lang, drop = TRUE])), , drop = FALSE]
  # and where the original string is not equal to the string in the target language
  df_trans <- df_trans[which(df_trans[, "pattern", drop = TRUE] != df_trans[, lang, drop = TRUE]), , drop = FALSE]
  if (only_unknown == TRUE) {
    df_trans <- subset(df_trans, pattern %like% "unknown")
  }
  if (only_affect_ab_names == TRUE) {
    df_trans <- subset(df_trans, affect_ab_name == TRUE)
  }
  if (only_affect_mo_names == TRUE) {
    df_trans <- subset(df_trans, affect_mo_name == TRUE)
  }
  if (NROW(df_trans) == 0) {
    return(from)
  }
  
  # default: case sensitive if value if 'case_sensitive' is missing:
  df_trans$case_sensitive[is.na(df_trans$case_sensitive)] <- TRUE
  # default: not using regular expressions if 'regular_expr' is missing:
  df_trans$regular_expr[is.na(df_trans$regular_expr)] <- FALSE
  
  # check if text to look for is in one of the patterns
  any_form_in_patterns <- tryCatch(
    any(from_unique %like% paste0("(", paste(gsub(" +\\(.*", "", df_trans$pattern), collapse = "|"), ")")),
    error = function(e) {
      warning_("Translation not possible. Please open an issue on GitHub (https://github.com/msberends/AMR/issues).")
      return(FALSE)
    })
  
  if (NROW(df_trans) == 0 | !any_form_in_patterns) {
    return(from)
  }
  
  lapply(seq_len(nrow(df_trans)), 
         function(i) from_unique_translated <<- gsub(pattern = df_trans$pattern[i],
                                                     replacement = df_trans[i, lang, drop = TRUE],
                                                     x = from_unique_translated,
                                                     ignore.case = !df_trans$case_sensitive[i] & df_trans$regular_expr[i], 
                                                     fixed = !df_trans$regular_expr[i],
                                                     perl = df_trans$regular_expr[i]))
  
  # force UTF-8 for diacritics
  from_unique_translated <- enc2utf8(from_unique_translated)
  
  # a kind of left join to get all results back
  from_unique_translated[match(from.bak, from_unique)]
}
