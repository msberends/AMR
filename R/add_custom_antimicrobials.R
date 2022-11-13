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

#' Add Custom Antimicrobials to This Package
#'
#' With [add_custom_antimicrobials()] you can add your own custom antimicrobial drug codes to the `AMR` package.
#' @param x a [data.frame] resembling the [antibiotics] data set, at least containing columns "ab" and "name"
#' @details Due to how \R works, the [add_custom_antimicrobials()] function has to be run in every \R session - added antimicrobials are not stored between sessions and are thus lost when \R is exited. It is possible to save the antimicrobial additions to your `.Rprofile` file to circumvent this, although this requires to load the `AMR` package at every start-up:
#'
#' ```r
#' # Open .Rprofile file
#' utils::file.edit("~/.Rprofile")
#'
#' # Add custom antibiotic drug codes:
#' library(AMR)
#' add_custom_antimicrobials(
#'   data.frame(ab = "TESTAB",
#'              name = "Test Antibiotic",
#'              group = "Test Group")
#' )
#' ```
#'
#' Use [clear_custom_antimicrobials()] to clear the previously added antimicrobials.
#' @rdname add_custom_antimicrobials
#' @export
#' @examples
#' \donttest{
#'
#' # returns NA and throws a warning (which is now suppressed):
#' suppressWarnings(
#'   as.ab("testab")
#' )
#'
#' # now add a custom entry - it will be considered by as.ab() and
#' # all ab_*() functions
#' add_custom_antimicrobials(
#'   data.frame(
#'     ab = "TESTAB",
#'     name = "Test Antibiotic",
#'     # you can add any property present in the
#'     # 'antibiotics' data set, such as 'group':
#'     group = "Test Group"
#'   )
#' )
#'
#' # "testab" is now a new antibiotic:
#' as.ab("testab")
#' ab_name("testab")
#' ab_group("testab")
#'
#' ab_info("testab")
#'
#'
#' # Add Co-fluampicil, which is one of the many J01CR50 codes, see
#' # https://www.whocc.no/ddd/list_of_ddds_combined_products/
#' add_custom_antimicrobials(
#'   data.frame(
#'     ab = "COFLU",
#'     name = "Co-fluampicil",
#'     atc = "J01CR50",
#'     group = "Beta-lactams/penicillines"
#'   )
#' )
#' ab_atc("Co-fluampicil")
#' ab_name("J01CR50")
#'
#' # even antibiotic selectors work
#' x <- data.frame(
#'   random_column = "some value",
#'   coflu = as.rsi("S"),
#'   ampicillin = as.rsi("R")
#' )
#' x
#' x[, betalactams()]
#' }
add_custom_antimicrobials <- function(x) {
  meet_criteria(x, allow_class = "data.frame")
  stop_ifnot(
    all(c("ab", "name") %in% colnames(x)),
    "`x` must contain columns \"ab\" and \"name\"."
  )
  stop_if(
    any(x$ab %in% AMR_env$AB_lookup$ab),
    "Antimicrobial drug code(s) ", vector_and(x$ab[x$ab %in% AMR_env$AB_lookup$ab]), " already exist in the internal `antibiotics` data set."
  )

  x <- x[, colnames(AMR_env$AB_lookup)[colnames(AMR_env$AB_lookup) %in% colnames(x)], drop = FALSE]
  x$generalised_name <- generalise_antibiotic_name(x$name)
  x$generalised_all <- as.list(x$generalised_name)
  if ("atc" %in% colnames(x)) {
    x$atc <- as.list(x$atc)
  }
  if ("loinc" %in% colnames(x)) {
    x$loinc <- as.list(x$loinc)
  }
  AMR_env$custom_ab_codes <- c(AMR_env$custom_ab_codes, x$ab)
  class(AMR_env$AB_lookup$ab) <- "character"

  new_df <- AMR_env$AB_lookup[0, , drop = FALSE][seq_len(NROW(x)), , drop = FALSE]
  rownames(new_df) <- NULL
  list_cols <- vapply(FUN.VALUE = logical(1), new_df, is.list)
  for (l in which(list_cols)) {
    # prevent binding NULLs in lists, replace with NA
    new_df[, l] <- as.list(NA_character_)
  }
  for (col in colnames(x)) {
    # assign new values
    new_df[, col] <- x[, col, drop = TRUE]
  }
  AMR_env$AB_lookup <- unique(rbind(AMR_env$AB_lookup, new_df))
  class(AMR_env$AB_lookup$ab) <- c("ab", "character")
  message_("Added ", nr2char(nrow(x)), " record", ifelse(nrow(x) > 1, "s", ""), " to the internal `antibiotics` data set.")
}

#' @rdname add_custom_antimicrobials
#' @export
clear_custom_antimicrobials <- function() {
  AMR_env$AB_lookup <- create_AB_lookup()
  AMR_env$custom_ab_codes <- character(0)
  message_("Custom antimicrobials cleared.")
}
