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

#' Principal Component Analysis (for AMR)
#'
#' Performs a principal component analysis (PCA) based on a data set with automatic determination for afterwards plotting the groups and labels, and automatic filtering on only suitable (i.e. non-empty and numeric) variables.
#' @param x A [data.frame] containing [numeric] columns.
#' @param ... Columns of `x` to be selected for PCA, can be unquoted since it supports quasiquotation.
#' @inheritParams stats::prcomp
#' @details The [pca()] function takes a [data.frame] as input and performs the actual PCA with the \R function [prcomp()].
#'
#' The result of the [pca()] function is a [prcomp] object, with an additional attribute `non_numeric_cols` which is a vector with the column names of all columns that do not contain [numeric] values. These are probably the groups and labels, and will be used by [ggplot_pca()].
#' @return An object of classes [pca] and [prcomp]
#' @importFrom stats prcomp
#' @export
#' @examples
#' # `example_isolates` is a data set available in the AMR package.
#' # See ?example_isolates.
#'
#' \donttest{
#' if (require("dplyr")) {
#'   # calculate the resistance per group first
#'   resistance_data <- example_isolates %>%
#'     group_by(
#'       order = mo_order(mo), # group on anything, like order
#'       genus = mo_genus(mo)
#'     ) %>% #   and genus as we do here;
#'     filter(n() >= 30) %>% # filter on only 30 results per group
#'     summarise_if(is.sir, resistance) # then get resistance of all drugs
#'
#'   # now conduct PCA for certain antimicrobial drugs
#'   pca_result <- resistance_data %>%
#'     pca(AMC, CXM, CTX, CAZ, GEN, TOB, TMP, SXT)
#'
#'   pca_result
#'   summary(pca_result)
#'   # old base R plotting method:
#'   biplot(pca_result)
#' }
#'
#' # new ggplot2 plotting method using this package:
#' if (require("dplyr") && require("ggplot2")) {
#'     ggplot_pca(pca_result)
#' }
#' if (require("dplyr") && require("ggplot2")) {
#'     ggplot_pca(pca_result) +
#'       scale_colour_viridis_d() +
#'       labs(title = "Title here")
#' }
#' }
pca <- function(x,
                ...,
                retx = TRUE,
                center = TRUE,
                scale. = TRUE,
                tol = NULL,
                rank. = NULL) {
  meet_criteria(x, allow_class = "data.frame")
  meet_criteria(retx, allow_class = "logical", has_length = 1)
  meet_criteria(center, allow_class = "logical", has_length = 1)
  meet_criteria(scale., allow_class = "logical", has_length = 1)
  meet_criteria(tol, allow_class = "numeric", has_length = 1, allow_NULL = TRUE)
  meet_criteria(rank., allow_class = "numeric", has_length = 1, allow_NULL = TRUE)

  # unset data.table, tibble, etc.
  # also removes groups made by dplyr::group_by
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  x.bak <- x

  # defuse R expressions, this replaces rlang::enquos()
  dots <- substitute(list(...))
  if (length(dots) > 1) {
    new_list <- list(0)
    for (i in seq_len(length(dots) - 1)) {
      new_list[[i]] <- tryCatch(eval(dots[[i + 1]], envir = x),
        error = function(e) stop(conditionMessage(e), call. = FALSE)
      )
      if (length(new_list[[i]]) == 1) {
        if (is.character(new_list[[i]]) && new_list[[i]] %in% colnames(x)) {
          # this is to support quoted variables: df %pm>% pca("mycol1", "mycol2")
          new_list[[i]] <- x[, new_list[[i]]]
        } else {
          # remove item - it's an argument like `center`
          new_list[[i]] <- NULL
        }
      }
    }

    x <- as.data.frame(new_list, stringsAsFactors = FALSE)
    if (any(vapply(FUN.VALUE = logical(1), x, function(y) !is.numeric(y)))) {
      warning_("in `pca()`: be sure to first calculate the resistance (or susceptibility) of variables with antimicrobial test results, since PCA works with numeric variables only. See Examples in `?pca`.", call = FALSE)
    }

    # set column names
    tryCatch(colnames(x) <- as.character(dots)[2:length(dots)],
      error = function(e) warning("column names could not be set")
    )

    # keep only numeric columns
    x <- x[, vapply(FUN.VALUE = logical(1), x, function(y) is.numeric(y)), drop = FALSE]
    # bind the data set with the non-numeric columns
    x <- cbind(x.bak[, vapply(FUN.VALUE = logical(1), x.bak, function(y) !is.numeric(y) & !all(is.na(y))), drop = FALSE], x)
  }

  x <- pm_ungroup(x) # would otherwise select the grouping vars
  x <- x[rowSums(is.na(x)) == 0, ] # remove columns containing NAs

  pca_data <- x[, which(vapply(FUN.VALUE = logical(1), x, function(x) is.numeric(x))), drop = FALSE]

  message_(
    "Columns selected for PCA: ", vector_and(font_bold(colnames(pca_data), collapse = NULL), quotes = TRUE),
    ". Total observations available: ", nrow(pca_data), "."
  )

  if (getRversion() < "3.4.0") {
    # stats::prcomp prior to 3.4.0 does not have the 'rank.' argument
    pca_model <- prcomp(pca_data, retx = retx, center = center, scale. = scale., tol = tol)
  } else {
    pca_model <- prcomp(pca_data, retx = retx, center = center, scale. = scale., tol = tol, rank. = rank.)
  }
  groups <- x[, vapply(FUN.VALUE = logical(1), x, function(y) !is.numeric(y) & !all(is.na(y))), drop = FALSE]
  rownames(groups) <- NULL
  attr(pca_model, "non_numeric_cols") <- groups
  class(pca_model) <- c("pca", class(pca_model))
  pca_model
}

#' @method print pca
#' @export
#' @noRd
print.pca <- function(x, ...) {
  a <- attributes(x)$non_numeric_cols
  if (!is.null(a)) {
    print_pca_group(a)
    class(x) <- class(x)[class(x) != "pca"]
  }
  print(x, ...)
}

#' @method summary pca
#' @export
#' @noRd
summary.pca <- function(object, ...) {
  a <- attributes(object)$non_numeric_cols
  if (!is.null(a)) {
    print_pca_group(a)
    class(object) <- class(object)[class(object) != "pca"]
  }
  summary(object, ...)
}

print_pca_group <- function(a) {
  grps <- sort(unique(a[, 1, drop = TRUE]))
  cat("Groups (n=", length(grps), ", named as '", colnames(a)[1], "'):\n", sep = "")
  print(grps)
  cat("\n")
}
