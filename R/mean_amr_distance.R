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

#' Mean AMR Distance
#'
#' This function calculates a normalised mean for antimicrobial resistance between multiple observations.
#' @param x a vector of class [rsi][as.rsi()], [rsi][as.rsi()] or [rsi][as.rsi()], or a [data.frame] containing columns of any of these classes
#' @param ... variables to select (supports tidy selection such as `column1:column4` and [`where(is.mic)`][tidyselect::language]), and can thus also be [antibiotic selectors][ab_selector()]
#' @param combine_SI 	a [logical] to indicate whether all values of S and I must be merged into one, so the input only consists of S+I vs. R (susceptible vs. resistant), defaults to `TRUE`
#' @details The mean AMR distance is a normalised numeric value to compare AMR test results and can help to identify similar isolates, without comparing antibiograms by hand. For common numeric data this distance is equal to [Z scores](https://en.wikipedia.org/wiki/Standard_score) (the number of standard deviations from the mean).
#'
#' MIC values (see [as.mic()]) are transformed with [log2()] first; their distance is calculated as `(log2(x) - mean(log2(x))) / sd(log2(x))`.
#'
#' R/SI values (see [as.rsi()]) are transformed using `"S"` = 1, `"I"` = 2, and `"R"` = 3. If `combine_SI` is `TRUE` (default), the `"I"` will be considered to be 1.
#'
#' For data sets, the mean AMR distance will be calculated per variable, after which the mean of all columns will returned per row (using [rowMeans()]), see *Examples*.
#'
#' Use [mean_distance_from_row()] to subtract distances from the distance of one row, see *Examples*.
#' @section Interpretation:
#' Isolates with distances less than 0.01 difference from each other should be considered similar. Differences lower than 0.025 should be considered suspicious.
#' @export
#' @examples
#' x <- random_mic(10)
#' x
#' mean_amr_distance(x)
#'
#' y <- data.frame(
#'   id = LETTERS[1:10],
#'   amox = random_mic(10, ab = "amox", mo = "Escherichia coli"),
#'   cipr = random_mic(10, ab = "cipr", mo = "Escherichia coli"),
#'   gent = random_mic(10, ab = "gent", mo = "Escherichia coli"),
#'   tobr = random_mic(10, ab = "tobr", mo = "Escherichia coli")
#' )
#' y
#' mean_amr_distance(y)
#' y$amr_distance <- mean_amr_distance(y, where(is.mic))
#' y[order(y$amr_distance), ]
#'
#' if (require("dplyr")) {
#'   y %>%
#'     mutate(
#'       amr_distance = mean_amr_distance(., where(is.mic)),
#'       check_id_C = mean_distance_from_row(amr_distance, id == "C")
#'     ) %>%
#'     arrange(check_id_C)
#' }
#' if (require("dplyr")) {
#'   # support for groups
#'   example_isolates %>%
#'     filter(mo_genus() == "Enterococcus" & mo_species() != "") %>%
#'     select(mo, TCY, carbapenems()) %>%
#'     group_by(mo) %>%
#'     mutate(d = mean_amr_distance(., where(is.rsi))) %>%
#'     arrange(mo, d)
#' }
mean_amr_distance <- function(x, ...) {
  UseMethod("mean_amr_distance")
}

#' @rdname mean_amr_distance
#' @export
mean_amr_distance.default <- function(x, ...) {
  x <- as.double(x)
  (x - mean(x, na.rm = TRUE)) / stats::sd(x, na.rm = TRUE)
}

#' @rdname mean_amr_distance
#' @export
mean_amr_distance.mic <- function(x, ...) {
  mean_amr_distance(log2(x))
}

#' @rdname mean_amr_distance
#' @export
mean_amr_distance.disk <- function(x, ...) {
  mean_amr_distance(as.double(x))
}

#' @rdname mean_amr_distance
#' @export
mean_amr_distance.rsi <- function(x, ..., combine_SI = TRUE) {
  meet_criteria(combine_SI, allow_class = "logical", has_length = 1, .call_depth = -1)
  if (isTRUE(combine_SI)) {
    x[x == "I"] <- "S"
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
  if (tryCatch(length(list(...)) > 0, error = function(e) TRUE)) {
    out <- tryCatch(suppressWarnings(c(...)), error = function(e) NULL)
    if (!is.null(out)) {
      df <- df[, out, drop = FALSE]
    } else {
      df <- pm_select(df, ...)
    }
  }
  stop_if(ncol(df) < 2,
    "data set must contain at least two variables",
    call = -2
  )
  if (message_not_thrown_before("mean_amr_distance", "groups")) {
    message_("Calculating mean AMR distance based on columns ", vector_and(colnames(df)))
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
  res[is.infinite(res)] <- 0
  res
}

#' @rdname mean_amr_distance
#' @param mean_distance the outcome of [mean_amr_distance()]
#' @param row an index, such as a row number
#' @export
mean_distance_from_row <- function(mean_distance, row) {
  meet_criteria(mean_distance, allow_class = c("double", "numeric"), is_finite = TRUE)
  meet_criteria(row, allow_class = c("logical", "double", "numeric"))
  if (is.logical(row)) {
    row <- which(row)
  }
  abs(mean_distance[row] - mean_distance)
}
