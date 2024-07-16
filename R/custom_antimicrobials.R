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

#' Add Custom Antimicrobials
#'
#' With [add_custom_antimicrobials()] you can add your own custom antimicrobial drug names and codes.
#' @param x a [data.frame] resembling the [antibiotics] data set, at least containing columns "ab" and "name"
#' @details **Important:** Due to how \R works, the [add_custom_antimicrobials()] function has to be run in every \R session - added antimicrobials are not stored between sessions and are thus lost when \R is exited.
#'
#' There are two ways to circumvent this and automate the process of adding antimicrobials:
#'
#' **Method 1:** Using the [package option][AMR-options] [`AMR_custom_ab`][AMR-options], which is the preferred method. To use this method:
#'
#'    1. Create a data set in the structure of the [antibiotics] data set (containing at the very least columns "ab" and "name") and save it with [saveRDS()] to a location of choice, e.g. `"~/my_custom_ab.rds"`, or any remote location.
#'
#'    2. Set the file location to the [package option][AMR-options] [`AMR_custom_ab`][AMR-options]: `options(AMR_custom_ab = "~/my_custom_ab.rds")`. This can even be a remote file location, such as an https URL. Since options are not saved between \R sessions, it is best to save this option to the `.Rprofile` file so that it will be loaded on start-up of \R. To do this, open the `.Rprofile` file using e.g. `utils::file.edit("~/.Rprofile")`, add this text and save the file:
#'
#'       ```r
#'       # Add custom antimicrobial codes:
#'       options(AMR_custom_ab = "~/my_custom_ab.rds")
#'       ```
#'
#'       Upon package load, this file will be loaded and run through the [add_custom_antimicrobials()] function.
#'
#' **Method 2:** Loading the antimicrobial additions directly from your `.Rprofile` file. Note that the definitions will be stored in a user-specific \R file, which is a suboptimal workflow. To use this method:
#'
#'    1. Edit the `.Rprofile` file using e.g. `utils::file.edit("~/.Rprofile")`.
#'
#'    2. Add a text like below and save the file:
#'
#'       ```r
#'        # Add custom antibiotic drug codes:
#'        AMR::add_custom_antimicrobials(
#'          data.frame(ab = "TESTAB",
#'                     name = "Test Antibiotic",
#'                     group = "Test Group")
#'        )
#'       ```
#'
#' Use [clear_custom_antimicrobials()] to clear the previously added antimicrobials.
#' @seealso [add_custom_microorganisms()] to add custom microorganisms.
#' @rdname add_custom_antimicrobials
#' @export
#' @examples
#' \donttest{
#'
#' # returns NA and throws a warning (which is suppressed here):
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
#' # https://atcddd.fhi.no/ddd/list_of_ddds_combined_products/
#' add_custom_antimicrobials(
#'   data.frame(
#'     ab = "COFLU",
#'     name = "Co-fluampicil",
#'     atc = "J01CR50",
#'     group = "Beta-lactams/penicillins"
#'   )
#' )
#' ab_atc("Co-fluampicil")
#' ab_name("J01CR50")
#'
#' # even antibiotic selectors work
#' x <- data.frame(
#'   random_column = "some value",
#'   coflu = as.sir("S"),
#'   ampicillin = as.sir("R")
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
  # remove any extra class/type, such as grouped tbl, or data.table:
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  # keep only columns available in the antibiotics data set
  x <- x[, colnames(AMR_env$AB_lookup)[colnames(AMR_env$AB_lookup) %in% colnames(x)], drop = FALSE]
  x$generalised_name <- generalise_antibiotic_name(x$name)
  x$generalised_all <- as.list(x$generalised_name)
  for (col in colnames(x)) {
    if (is.list(AMR_env$AB_lookup[, col, drop = TRUE]) & !is.list(x[, col, drop = TRUE])) {
      x[, col] <- as.list(x[, col, drop = TRUE])
    }
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
  AMR_env$AB_lookup <- unique(rbind_AMR(AMR_env$AB_lookup, new_df))

  AMR_env$ab_previously_coerced <- AMR_env$ab_previously_coerced[which(!AMR_env$ab_previously_coerced$ab %in% x$ab), , drop = FALSE]
  class(AMR_env$AB_lookup$ab) <- c("ab", "character")
  message_("Added ", nr2char(nrow(x)), " record", ifelse(nrow(x) > 1, "s", ""), " to the internal `antibiotics` data set.")
}

#' @rdname add_custom_antimicrobials
#' @export
clear_custom_antimicrobials <- function() {
  n <- nrow(AMR_env$AB_lookup)
  AMR_env$AB_lookup <- cbind(AMR::antibiotics, AB_LOOKUP)
  n2 <- nrow(AMR_env$AB_lookup)
  AMR_env$custom_ab_codes <- character(0)
  AMR_env$ab_previously_coerced <- AMR_env$ab_previously_coerced[which(AMR_env$ab_previously_coerced$ab %in% AMR_env$AB_lookup$ab), , drop = FALSE]
  message_("Cleared ", nr2char(n - n2), " custom record", ifelse(n - n2 > 1, "s", ""), " from the internal `antibiotics` data set.")
}
