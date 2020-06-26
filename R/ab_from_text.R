# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Analysis                              #
#                                                                      #
# SOURCE                                                               #
# https://gitlab.com/msberends/AMR                                     #
#                                                                      #
# LICENCE                                                              #
# (c) 2018-2020 Berends MS, Luz CF et al.                              #
#                                                                      #
# This R package is free software; you can freely use and distribute   #
# it for both personal and commercial purposes under the terms of the  #
# GNU General Public License version 2.0 (GNU GPL-2), as published by  #
# the Free Software Foundation.                                        #
#                                                                      #
# We created this package for both routine data analysis and academic  #
# research and it was publicly released in the hope that it will be    #
# useful, but it comes WITHOUT ANY WARRANTY OR LIABILITY.              #
# Visit our website for more info: https://msberends.gitlab.io/AMR.    #
# ==================================================================== #

#' Retrieve antimicrobial drugs from clinical text
#' 
#' Use this function on e.g. clinical texts from health care records. It returns a [list] with all antimicrobial drugs found in the texts.
#' @param text text to analyse
#' @param collapse character to pass on to `paste(..., collapse = ...)` to only return one character per element of `text`, see Examples
#' @param translate_ab a column name of the [antibiotics] data set to translate the antibiotic abbreviations to, using [ab_property()]. Defaults to `FALSE`. Using `TRUE` is equal to using "name".
#' @param ... parameters passed on to [as.ab()]
#' @details Without using `collapse`, this function will return a [list]. This can be convenient to use e.g. inside a `mutate()`):\cr
#' `df %>% mutate(abx = ab_from_text(clinical_text))` 
#' 
#' The returned AB codes can be transformed to official names, groups, etc. with all [ab_property()] functions like [ab_name()] and [ab_group()], or by using the `translate_ab` parameter.
#' 
#' With using `collapse`, this function will return a [character]:\cr
#' `df %>% mutate(abx = ab_from_text(clinical_text, collapse = "|"))` 
#' 
#' This function is also internally used by [as.ab()], although it then only returns the first hit and will throw a note if more results could have been returned.
#' @export
#' @return A [list], or a [character] if `collapse` is not `NULL`
#' @examples 
#' # mind the bad spelling of amoxicillin in this line, 
#' # straight from a true health care record:
#' ab_from_text("28/03/2020 regular amoxicilliin 500mg po tds")
#' 
#' ab_from_text("administered amoxi/clav and cipro")
#' ab_from_text("administered amoxi/clav and cipro", collapse = ", ")
#' 
#' # if you want to know which antibiotic groups were administered, check it:
#' abx <- ab_from_text("administered amoxi/clav and cipro")
#' ab_group(abx[[1]])
#' 
#' if (require(dplyr)) {
#'   tibble(clinical_text = c("given cipro and mero",
#'                            "started on doxy today")) %>% 
#'     mutate(abx = ab_from_text(clinical_text),
#'            abx2 = ab_from_text(clinical_text,
#'                                collapse = "|"),
#'            abx3 = ab_from_text(clinical_text,
#'                                collapse = "|",
#'                                translate_ab = "name"))
#' 
#' }
ab_from_text <- function(text, collapse = NULL, translate_ab = FALSE, ...) {
  
  text <- tolower(as.character(text))
  translate_ab <- get_translate_ab(translate_ab)
  
  abbr <- unlist(antibiotics$abbreviations)
  abbr <- abbr[nchar(abbr) >= 4]
  names <- substr(antibiotics$name, 1, 5)
  synonyms <- unlist(antibiotics$synonyms)
  synonyms <- synonyms[nchar(synonyms) >= 4]
  to_regex <- function(x) {
    paste0("^(",
           paste0(unique(gsub("[^a-z0-9]", ".*", sort(tolower(x)))), collapse = "|"),
           ").*")
  }
  
  text_split_all <- strsplit(text, "[ ;.,:/\\|-]")
  result <- lapply(text_split_all, function(text_split) {
    suppressWarnings(
      out <- as.ab(unique(c(text_split[grep(to_regex(abbr), text_split)],
                     text_split[grep(to_regex(names), text_split)],
                     # regular expression must not be too long, so split synonyms in two:
                     text_split[grep(to_regex(synonyms[c(1:0.5 * length(synonyms))]), text_split)],
                     text_split[grep(to_regex(synonyms[c(0.5 * length(synonyms):length(synonyms))]), text_split)])),
            ...))
    out <- out[!is.na(out)]
    if (length(out) == 0) {
      as.ab(NA)
    } else {
      if (!isFALSE(translate_ab)) {
        out <- ab_property(out, property = translate_ab, initial = FALSE)
      }
      out
    }
  })

  if (!is.null(collapse)) {
    result <- sapply(result, function(x) paste0(x, collapse = collapse))
  }

  result
}
