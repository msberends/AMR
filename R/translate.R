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
#' @inheritSection lifecycle Stable Lifecycle
#' @param x text to translate
#' @param lang language to choose. Use one of these supported language names or ISO-639-1 codes: `r paste0('"', names(LANGUAGES_SUPPORTED), '" ("' , LANGUAGES_SUPPORTED, '")', collapse = ", ")`.
#' @details The currently `r length(LANGUAGES_SUPPORTED)` supported languages are `r vector_and(names(LANGUAGES_SUPPORTED), quotes = FALSE, sort = FALSE)`. All these languages have translations available for all antimicrobial agents and colloquial microorganism names.
#'
#' Please suggest your own translations [by creating a new issue on our repository](https://github.com/msberends/AMR/issues/new?title=Translations). Strings will be translated to foreign languages if they are defined in [this repository file](https://github.com/msberends/AMR/blob/main/data-raw/translations.tsv). This file will be read by all functions where a translated output can be desired, like all [`mo_*`][mo_property()] functions (such as [mo_name()], [mo_gramstain()], [mo_type()], etc.) and [`ab_*`][ab_property()] functions (such as [ab_name()], [ab_group()], etc.). 
#'
#' ## Changing the Default Language
#' The system language will be used at default (as returned by `Sys.getenv("LANG")` or, if `LANG` is not set, [Sys.getlocale("LC_COLLATE")]), if that language is supported. But the language to be used can be overwritten in two ways and will be checked in this order:
#' 
#'   1. Setting the R option `AMR_locale`, either by using `set_AMR_locale()` or by running e.g. `options(AMR_locale = "de")`.
#'   
#'      Note that setting an \R option only works in the same session. Save the command `options(AMR_locale = "(your language)")` to your `.Rprofile` file to apply it for every session.
#'   2. Setting the system variable `LANGUAGE` or `LANG`, e.g. by adding `LANGUAGE="de_DE.utf8"` to your `.Renviron` file in your home directory.
#' 
#' Thus, if the R option `AMR_locale` is set, the system variables `LANGUAGE` and `LANG` will be ignored.
#' @inheritSection AMR Read more on Our Website!
#' @rdname translate
#' @name translate
#' @export
#' @examples
#' # Current settings
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
    message_("Assuming the ", names(LANGUAGES_SUPPORTED)[LANGUAGES_SUPPORTED == lang],
             " language for the AMR package. Change this with `set_AMR_locale()`. ",
             "This note will be shown once per session.")
  }
  lang
}

#' @rdname translate
#' @export
set_AMR_locale <- function(lang) {
  lang <- validate_language(lang)
  options(AMR_locale = lang)
  message_("Using the ", names(LANGUAGES_SUPPORTED)[LANGUAGES_SUPPORTED == lang],  " language for the AMR package for this session.")
}

#' @rdname translate
#' @export
reset_AMR_locale <- function() {
  options(AMR_locale = NULL)
  message_("Language for the AMR package reset to English for this session.")
}

#' @rdname translate
#' @export
translate_AMR <- function(x, language = get_AMR_locale()) {
  translate_into_language(x, language = language)
}

validate_language <- function(language, extra_txt = character(0)) {
  language.bak <- language
  language <- LANGUAGES_SUPPORTED[which(tolower(language) == LANGUAGES_SUPPORTED | tolower(names(LANGUAGES_SUPPORTED)) == tolower(language))][1]
  stop_ifnot(language %in% LANGUAGES_SUPPORTED,
             "unsupported language for AMR package", extra_txt, ": \"", language.bak, "\". Use one of these language names or ISO-639-1 codes: ",
             paste0('"', names(LANGUAGES_SUPPORTED), '" ("' , LANGUAGES_SUPPORTED, '")', collapse = ", "),
             call = FALSE)
  unname(language)
}

find_language <- function(lang) {
  # grepl() with ignore.case = FALSE is 8x faster than %like_case%
  if (grepl("^(English|en_|EN_)", lang, ignore.case = FALSE, perl = TRUE)) {
    # as first option to optimise speed
    "en"
  } else if (grepl("^(Chinese|zh_|ZH_)", lang, ignore.case = FALSE, perl = TRUE)) {
    "zh"
  } else if (grepl("^(Danish|Dansk|da_|DA_)", lang, ignore.case = FALSE, perl = TRUE)) {
    "da"
  } else if (grepl("^(Dutch|Nederlands|nl_|NL_)", lang, ignore.case = FALSE, perl = TRUE)) {
    "nl"
  } else if (grepl("^(French|Fran.+ais|fr_|FR_)", lang, ignore.case = FALSE, perl = TRUE)) {
    "fr"
  } else if (grepl("^(German|Deutsch|de_|DE_)", lang, ignore.case = FALSE, perl = TRUE)) {
    "de"
  } else if (grepl("^(Greek|el_|EL_)", lang, ignore.case = FALSE, perl = TRUE)) {
    "el"
  } else if (grepl("^(Italian|Italiano|it_|IT_)", lang, ignore.case = FALSE, perl = TRUE)) {
    "it"
  } else if (grepl("^(Japanese|ja_|JA_)", lang, ignore.case = FALSE, perl = TRUE)) {
    "ja"
  } else if (grepl("^(Polish|polsk|pl_|PL_)", lang, ignore.case = FALSE, perl = TRUE)) {
    "pl"
  } else if (grepl("^(Portuguese|Portugu.+s|pt_|PT_)", lang, ignore.case = FALSE, perl = TRUE)) {
    "pt"
  } else if (grepl("^(Russian|pycc|ru_|RU_)", lang, ignore.case = FALSE, perl = TRUE)) {
    "ru"
  } else if (grepl("^(Spanish|Espa.+ol|es_|ES_)", lang, ignore.case = FALSE, perl = TRUE)) {
    "es"
  } else if (grepl("^(Swedish|Svenskt|sv_|SV_)", lang, ignore.case = FALSE, perl = TRUE)) {
    "sv"
  } else if (grepl("^(Ukrainian|uk_|UK)", lang, ignore.case = FALSE, perl = TRUE)) {
    "uk"
  } else {
    # other language -> set to English
    "en"
  }
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
  
  # name of language used
  language.bak <- language
  language <- LANGUAGES_SUPPORTED[which(tolower(language) == LANGUAGES_SUPPORTED | tolower(names(LANGUAGES_SUPPORTED)) == tolower(language))][1]
  
  stop_ifnot(language %in% LANGUAGES_SUPPORTED,
             "unsupported language: \"", language.bak, "\" - use one of these language names or ISO-639-1 codes: ",
             paste0('"', names(LANGUAGES_SUPPORTED), '" ("' , LANGUAGES_SUPPORTED, '")', collapse = ", "),
             call = FALSE)
  
  # only keep lines where translation is available for this language
  df_trans <- df_trans[which(!is.na(df_trans[, language, drop = TRUE])), , drop = FALSE]
  # and where the original string is not equal to the string in the target language
  df_trans <- df_trans[which(df_trans[, "pattern", drop = TRUE] != df_trans[, language, drop = TRUE]), , drop = FALSE]
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
                                                     replacement = df_trans[i, language, drop = TRUE],
                                                     x = from_unique_translated,
                                                     ignore.case = !df_trans$case_sensitive[i] & df_trans$regular_expr[i], 
                                                     fixed = !df_trans$regular_expr[i],
                                                     perl = df_trans$regular_expr[i]))
  
  # force UTF-8 for diacritics
  from_unique_translated <- enc2utf8(from_unique_translated)
  
  # a kind of left join to get all results back
  from_unique_translated[match(from.bak, from_unique)]
}
