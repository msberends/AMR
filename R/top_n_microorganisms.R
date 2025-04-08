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
# how to conduct AMR data analysis: https://amr-for-r.org/             #
# ==================================================================== #

#' Filter Top *n* Microorganisms
#'
#' This function filters a data set to include only the top *n* microorganisms based on a specified property, such as taxonomic family or genus. For example, it can filter a data set to the top 3 species, or to any species in the top 5 genera, or to the top 3 species in each of the top 5 genera.
#' @param x A data frame containing microbial data.
#' @param n An integer specifying the maximum number of unique values of the `property` to include in the output.
#' @param property A character string indicating the microorganism property to use for filtering. Must be one of the column names of the [microorganisms] data set: `r vector_or(colnames(microorganisms), sort = FALSE, quotes = TRUE)`. If `NULL`, the raw values from `col_mo` will be used without transformation.
#' @param n_for_each An optional integer specifying the maximum number of rows to retain for each value of the selected property. If `NULL`, all rows within the top *n* groups will be included.
#' @param col_mo A character string indicating the column in `x` that contains microorganism names or codes. Defaults to the first column of class [`mo`]. Values will be coerced using [as.mo()].
#' @param ... Additional arguments passed on to [mo_property()] when `property` is not `NULL`.
#' @details This function is useful for preprocessing data before creating [antibiograms][antibiogram()] or other analyses that require focused subsets of microbial data. For example, it can filter a data set to only include isolates from the top 10 species.
#' @export
#' @seealso [mo_property()], [as.mo()], [antibiogram()]
#' @examples
#' # filter to the top 3 species:
#' top_n_microorganisms(example_isolates,
#'   n = 3
#' )
#'
#' # filter to any species in the top 5 genera:
#' top_n_microorganisms(example_isolates,
#'   n = 5, property = "genus"
#' )
#'
#' # filter to the top 3 species in each of the top 5 genera:
#' top_n_microorganisms(example_isolates,
#'   n = 5, property = "genus", n_for_each = 3
#' )
top_n_microorganisms <- function(x, n, property = "fullname", n_for_each = NULL, col_mo = NULL, ...) {
  meet_criteria(x, allow_class = "data.frame") # also checks dimensions to be >0
  meet_criteria(n, allow_class = c("numeric", "integer"), has_length = 1, is_finite = TRUE, is_positive = TRUE)
  meet_criteria(property, allow_class = "character", has_length = 1, is_in = colnames(AMR::microorganisms))
  meet_criteria(n_for_each, allow_class = c("numeric", "integer"), has_length = 1, is_finite = TRUE, is_positive = TRUE, allow_NULL = TRUE)
  meet_criteria(col_mo, allow_class = "character", has_length = 1, allow_NULL = TRUE, is_in = colnames(x))
  if (is.null(col_mo)) {
    col_mo <- search_type_in_df(x = x, type = "mo", info = TRUE)
    stop_if(is.null(col_mo), "`col_mo` must be set")
  }

  x.bak <- x

  x[, col_mo] <- as.mo(x[, col_mo, drop = TRUE], keep_synonyms = TRUE)

  if (is.null(property)) {
    x$prop_val <- x[[col_mo]]
  } else {
    x$prop_val <- mo_property(x[[col_mo]], property = property, ...)
  }
  counts <- sort(table(x$prop_val), decreasing = TRUE)

  n <- as.integer(n)
  if (length(counts) < n) {
    n <- length(counts)
  }
  count_values <- names(counts)[seq_len(n)]
  filtered_rows <- which(x$prop_val %in% count_values)

  if (!is.null(n_for_each)) {
    n_for_each <- as.integer(n_for_each)
    filtered_x <- x[filtered_rows, , drop = FALSE]
    filtered_rows <- do.call(
      c,
      lapply(split(filtered_x, filtered_x$prop_val), function(group) {
        top_values <- names(sort(table(group[[col_mo]]), decreasing = TRUE)[seq_len(n_for_each)])
        top_values <- top_values[!is.na(top_values)]
        which(x[[col_mo]] %in% top_values)
      })
    )
  }

  x.bak[filtered_rows, , drop = FALSE]
}
