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

#' Calculate the Mean AMR Distance
#'
#' Calculates a normalised mean for antimicrobial resistance between multiple observations, to help to identify similar isolates without comparing antibiograms by hand.
#' @param x A vector of class [sir][as.sir()], [mic][as.mic()] or [disk][as.disk()], or a [data.frame] containing columns of any of these classes.
#' @param ... Variables to select. Supports [tidyselect language][tidyselect::language] (such as `column1:column4` and `where(is.mic)`), and can thus also be [antimicrobial selectors][amr_selector()].
#' @param combine_SI A [logical] to indicate whether all values of S, SDD, and I must be merged into one, so the input only consists of S+I vs. R (susceptible vs. resistant) - the default is `TRUE`.
#' @details The mean AMR distance is effectively [the Z-score](https://en.wikipedia.org/wiki/Standard_score); a normalised numeric value to compare AMR test results which can help to identify similar isolates, without comparing antibiograms by hand.
#'
#' MIC values (see [as.mic()]) are transformed with [log2()] first; their distance is thus calculated as `(log2(x) - mean(log2(x))) / sd(log2(x))`.
#'
#' SIR values (see [as.sir()]) are transformed using `"S"` = 1, `"I"` = 2, and `"R"` = 3. If `combine_SI` is `TRUE` (default), the `"I"` will be considered to be 1.
#'
#' For data sets, the mean AMR distance will be calculated per column, after which the mean per row will be returned, see *Examples*.
#'
#' Use [amr_distance_from_row()] to subtract distances from the distance of one row, see *Examples*.
#' @section Interpretation:
#' Isolates with distances less than 0.01 difference from each other should be considered similar. Differences lower than 0.025 should be considered suspicious.
#' @export
#' @examples
#' sir <- random_sir(10)
#' sir
#' mean_amr_distance(sir)
#'
#' mic <- random_mic(10)
#' mic
#' mean_amr_distance(mic)
#' # equal to the Z-score of their log2:
#' (log2(mic) - mean(log2(mic))) / sd(log2(mic))
#'
#' disk <- random_disk(10)
#' disk
#' mean_amr_distance(disk)
#'
#' y <- data.frame(
#'   id = LETTERS[1:10],
#'   amox = random_sir(10, ab = "amox", mo = "Escherichia coli"),
#'   cipr = random_disk(10, ab = "cipr", mo = "Escherichia coli"),
#'   gent = random_mic(10, ab = "gent", mo = "Escherichia coli"),
#'   tobr = random_mic(10, ab = "tobr", mo = "Escherichia coli")
#' )
#' y
#' mean_amr_distance(y)
#' y$amr_distance <- mean_amr_distance(y, is.mic(y))
#' y[order(y$amr_distance), ]
#'
#' if (require("dplyr")) {
#'   y %>%
#'     mutate(
#'       amr_distance = mean_amr_distance(y),
#'       check_id_C = amr_distance_from_row(amr_distance, id == "C")
#'     ) %>%
#'     arrange(check_id_C)
#' }
#' if (require("dplyr")) {
#'   # support for groups
#'   example_isolates %>%
#'     filter(mo_genus() == "Enterococcus" & mo_species() != "") %>%
#'     select(mo, TCY, carbapenems()) %>%
#'     group_by(mo) %>%
#'     mutate(dist = mean_amr_distance(.)) %>%
#'     arrange(mo, dist)
#' }
mean_amr_distance <- function(x, ...) {
  UseMethod("mean_amr_distance")
}

#' @noRd
#' @export
mean_amr_distance.default <- function(x, ...) {
  x <- as.double(x)
  # calculate z-score
  (x - mean(x, na.rm = TRUE)) / stats::sd(x, na.rm = TRUE)
}

#' @noRd
#' @export
mean_amr_distance.mic <- function(x, ...) {
  mean_amr_distance(log2(x))
}

#' @noRd
#' @export
mean_amr_distance.disk <- function(x, ...) {
  mean_amr_distance(as.double(x))
}

#' @rdname mean_amr_distance
#' @export
mean_amr_distance.sir <- function(x, ..., combine_SI = TRUE) {
  meet_criteria(combine_SI, allow_class = "logical", has_length = 1, .call_depth = -1)
  if (isTRUE(combine_SI)) {
    x[x %in% c("I", "SDD")] <- "S"
  }
  mean_amr_distance(as.double(x))
}

#' @rdname mean_amr_distance
#' @export
mean_amr_distance.data.frame <- function(x, ..., combine_SI = TRUE) {
  meet_criteria(combine_SI, allow_class = "logical", has_length = 1, .call_depth = -1)
  df <- x
  if (is_null_or_grouped_tbl(df)) {
    df <- get_current_data("x", -2)
  }
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  if (tryCatch(length(list(...)) > 0, error = function(e) TRUE)) {
    out <- tryCatch(suppressWarnings(c(...)), error = function(e) NULL)
    if (!is.null(out)) {
      df <- df[, out, drop = FALSE]
    } else {
      df <- pm_select(df, ...)
    }
  }
  df_classes <- colnames(df)[vapply(FUN.VALUE = logical(1), df, function(x) is.disk(x) | is.mic(x) | is.disk(x), USE.NAMES = FALSE)]
  df_antibiotics <- unname(get_column_abx(df, info = FALSE))
  df <- df[, colnames(df)[colnames(df) %in% union(df_classes, df_antibiotics)], drop = FALSE]

  stop_if(ncol(df) < 2,
    "data set must contain at least two variables",
    call = -2
  )
  if (message_not_thrown_before("mean_amr_distance", "groups")) {
    message_("Calculating mean AMR distance based on columns ", vector_and(colnames(df), sort = FALSE))
  }

  res <- vapply(
    FUN.VALUE = double(nrow(df)),
    df,
    mean_amr_distance,
    combine_SI = combine_SI
  )
  if (is.null(dim(res))) {
    if (all(is.na(res))) {
      return(NA_real_)
    } else {
      return(mean(res, na.rm = TRUE))
    }
  }
  res <- rowMeans(res, na.rm = TRUE)
  res[is.infinite(res) | is.nan(res)] <- 0
  res
}

#' @rdname mean_amr_distance
#' @param amr_distance The outcome of [mean_amr_distance()].
#' @param row An index, such as a row number.
#' @export
amr_distance_from_row <- function(amr_distance, row) {
  meet_criteria(amr_distance, allow_class = "numeric", is_finite = TRUE)
  meet_criteria(row, allow_class = c("logical", "numeric", "integer"))
  if (is.logical(row)) {
    row <- which(row)
  }
  abs(amr_distance[row] - amr_distance)
}
