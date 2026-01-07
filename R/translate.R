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

#' Translate Strings from the AMR Package
#'
#' For language-dependent output of `AMR` functions, such as [mo_name()], [mo_gramstain()], [mo_type()] and [ab_name()].
#' @param x Text to translate.
#' @param language Language to choose. Use one of these supported language names or [ISO 639-1 codes](https://en.wikipedia.org/wiki/ISO_639-1): `r vector_or(paste0(sapply(LANGUAGES_SUPPORTED_NAMES, function(x) x[[1]]), " (" , LANGUAGES_SUPPORTED, ")"), quotes = FALSE, sort = FALSE)`.
#' @details The currently `r length(LANGUAGES_SUPPORTED)` supported languages are `r vector_and(paste0(sapply(LANGUAGES_SUPPORTED_NAMES, function(x) x[[1]]), " (" , LANGUAGES_SUPPORTED, ")"), quotes = FALSE, sort = FALSE)`. All these languages have translations available for all antimicrobial drugs and colloquial microorganism names.
#'
#' To permanently silence the once-per-session language note on a non-English operating system, you can set the package option [`AMR_locale`][AMR-options] in your `.Rprofile` file like this:
#'
#' ```r
#' # Open .Rprofile file
#' utils::file.edit("~/.Rprofile")
#'
#' # Then add e.g. Italian support to that file using:
#' options(AMR_locale = "Italian")
#' ```
#'
#' And then save the file.
#'
#' Please read about adding or updating a language in [our Wiki](https://github.com/msberends/AMR/wiki/).
#'
#' ### Changing the Default Language
#' The system language will be used at default (as returned by `Sys.getenv("LANG")` or, if `LANG` is not set, [`Sys.getlocale("LC_COLLATE")`][Sys.getlocale()]), if that language is supported. But the language to be used can be overwritten in two ways and will be checked in this order:
#'
#'   1. Setting the package option [`AMR_locale`][AMR-options], either by using e.g. `set_AMR_locale("German")` or by running e.g. `options(AMR_locale = "German")`.
#'
#'      Note that setting an \R option only works in the same session. Save the command `options(AMR_locale = "(your language)")` to your `.Rprofile` file to apply it for every session. Run `utils::file.edit("~/.Rprofile")` to edit your `.Rprofile` file.
#'   2. Setting the system variable `LANGUAGE` or `LANG`, e.g. by adding `LANGUAGE="de_DE.utf8"` to your `.Renviron` file in your home directory.
#'
#' Thus, if the package option [`AMR_locale`][AMR-options] is set, the system variables `LANGUAGE` and `LANG` will be ignored.
#' @rdname translate
#' @name translate
#' @export
#' @examples
#' # Current settings (based on system language)
#' ab_name("Ciprofloxacin")
#' mo_name("Coagulase-negative Staphylococcus (CoNS)")
#'
#' # setting another language
#' set_AMR_locale("Dutch")
#' ab_name("Ciprofloxacin")
#' mo_name("Coagulase-negative Staphylococcus (CoNS)")
#'
#' # setting yet another language
#' set_AMR_locale("German")
#' ab_name("Ciprofloxacin")
#' mo_name("Coagulase-negative Staphylococcus (CoNS)")
#'
#' # set_AMR_locale() understands endonyms, English exonyms, and ISO 639-1:
#' set_AMR_locale("Deutsch")
#' set_AMR_locale("German")
#' set_AMR_locale("de")
#' ab_name("amox/clav")
#'
#' # reset to system default
#' reset_AMR_locale()
#' ab_name("amox/clav")
get_AMR_locale <- function() {
  # a message for this will be thrown in translate_into_language() if outcome is non-English
  if (!is.null(getOption("AMR_locale"))) {
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
  find_language(lang)
}

#' @rdname translate
#' @export
set_AMR_locale <- function(language) {
  language <- validate_language(language)
  options(AMR_locale = language)
  if (interactive() || identical(Sys.getenv("IN_PKGDOWN"), "true")) {
    # show which language to use now
    message_(
      "Using ", LANGUAGES_SUPPORTED_NAMES[[language]]$exonym,
      ifelse(language != "en",
        paste0(" (", LANGUAGES_SUPPORTED_NAMES[[language]]$endonym, ")"),
        ""
      ),
      " for the AMR package for this session."
    )
  }
}

#' @rdname translate
#' @export
reset_AMR_locale <- function() {
  options(AMR_locale = NULL)
  if (interactive() || identical(Sys.getenv("IN_PKGDOWN"), "true")) {
    # show which language to use now
    language <- suppressMessages(get_AMR_locale())
    message_("Using the ", LANGUAGES_SUPPORTED_NAMES[[language]]$exonym, " language (", LANGUAGES_SUPPORTED_NAMES[[language]]$endonym, ") for the AMR package for this session.")
  }
}

#' @rdname translate
#' @export
translate_AMR <- function(x, language = get_AMR_locale()) {
  translate_into_language(x,
    language = language,
    only_unknown = FALSE,
    only_affect_ab_names = FALSE,
    only_affect_mo_names = FALSE
  )
}


validate_language <- function(language, extra_txt = character(0)) {
  if (length(language) == 0 || isTRUE(trimws2(tolower(language[1])) %in% c("en", "english", "", "false", NA))) {
    return("en")
  } else if (language[1] %in% LANGUAGES_SUPPORTED) {
    return(language[1])
  }
  lang <- find_language(language[1], fallback = FALSE)
  stop_ifnot(length(lang) > 0 && lang %in% LANGUAGES_SUPPORTED,
    "unsupported language for AMR package", extra_txt, ": \"", language, "\". Use one of these language names or ISO 639-1 codes: ",
    paste0('"', vapply(FUN.VALUE = character(1), LANGUAGES_SUPPORTED_NAMES, function(x) x[[1]]),
      '" ("', LANGUAGES_SUPPORTED, '")',
      collapse = ", "
    ),
    call = FALSE
  )
  lang
}

find_language <- function(language, fallback = TRUE) {
  language <- Map(LANGUAGES_SUPPORTED_NAMES,
    LANGUAGES_SUPPORTED,
    f = function(l, n, check = language) {
      grepl(
        paste0(
          "^(", l[1], "|", l[2], "|",
          n, "(_|$)|", toupper(n), "(_|$))"
        ),
        check,
        ignore.case = TRUE,
        perl = TRUE,
        useBytes = FALSE
      )
    },
    USE.NAMES = TRUE
  )
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
  # get ISO 639-1 of language
  lang <- validate_language(language)
  if (lang == "en") {
    # don' translate
    return(from)
  }

  df_trans <- TRANSLATIONS # internal data file
  from.bak <- from
  from_unique <- unique(from)
  from_split_combined <- function(vec) {
    sapply(vec, function(x) {
      if (grepl("/", x, fixed = TRUE)) {
        parts <- strsplit(x, "/", fixed = TRUE)[[1]]
        # Translate each part separately
        translated_parts <- translate_into_language(
          parts,
          language = lang,
          only_unknown = only_unknown,
          only_affect_ab_names = only_affect_ab_names,
          only_affect_mo_names = only_affect_mo_names
        )
        paste(translated_parts, collapse = "/")
      } else {
        x
      }
    }, USE.NAMES = FALSE)
  }
  from_unique_translated <- from_split_combined(from_unique)

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
      warning_("Translation not possible. Please create an issue at ", font_url("https://github.com/msberends/AMR/issues"), ". Many thanks!")
      return(FALSE)
    }
  )

  if (NROW(df_trans) == 0 | !any_form_in_patterns) {
    return(from)
  }

  if (only_affect_ab_names == TRUE) {
    df_trans$pattern[df_trans$regular_expr == TRUE] <- paste0(df_trans$pattern[df_trans$regular_expr == TRUE], "$")
    df_trans$pattern[df_trans$regular_expr == TRUE] <- gsub("$$", "$", df_trans$pattern[df_trans$regular_expr == TRUE], fixed = TRUE)
  }

  df_trans_regex <- df_trans[which(df_trans$regular_expr == TRUE), ]
  # regex part
  lapply(
    # starting with longest pattern, since more general translations are shorter, such as 'Group'
    order(nchar(df_trans_regex$pattern), decreasing = TRUE),
    function(i) {
      from_unique_translated <<- gsub(
        pattern = df_trans_regex$pattern[i],
        replacement = df_trans_regex[i, lang, drop = TRUE],
        x = from_unique_translated,
        ignore.case = !df_trans_regex$case_sensitive[i],
        fixed = FALSE,
        perl = TRUE
      )
    }
  )
  # non-regex part
  translate_tokens <- function(tokens) {
    patterns <- df_trans$pattern[df_trans$regular_expr == FALSE]
    replacements <- df_trans[[lang]][df_trans$regular_expr == FALSE]
    matches <- match(tokens, patterns)
    tokens[!is.na(matches)] <- replacements[matches[!is.na(matches)]]
    tokens
  }
  from_unique_translated <- vapply(
    FUN.VALUE = character(1),
    USE.NAMES = FALSE,
    from_unique_translated,
    function(x) {
      delimiters <- "[ /()]"
      split_regex <- paste0("(?<=", delimiters, ")|(?=", delimiters, ")")
      tokens <- strsplit(x, split_regex, perl = TRUE)[[1]]
      tokens <- translate_tokens(tokens)
      paste(tokens, collapse = "")
    }
  )

  # force UTF-8 for diacritics
  from_unique_translated <- enc2utf8(from_unique_translated)

  # a kind of left join to get all results back
  out <- from_unique_translated[match(from.bak, from_unique)]

  if (!identical(from.bak, out) && get_AMR_locale() == lang && is.null(getOption("AMR_locale", default = NULL)) && message_not_thrown_before("translation", entire_session = TRUE) && interactive()) {
    message(word_wrap(
      "Assuming the ", LANGUAGES_SUPPORTED_NAMES[[lang]]$exonym, " language (",
      LANGUAGES_SUPPORTED_NAMES[[lang]]$endonym, ") for the AMR package. See `set_AMR_locale()` to change this or to silence this once-per-session note.",
      add_fn = list(font_blue), as_note = TRUE
    ))
  }

  out
}
