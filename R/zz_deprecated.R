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

#' Deprecated Functions, Arguments, or Datasets
#'
#' These objects are so-called '[Deprecated]'. **They will be removed in a future version of this package.** Using these will give a warning with the name of the alternative object it has been replaced by (if there is one).
#' @keywords internal
#' @name AMR-deprecated
#' @rdname AMR-deprecated
NULL

.amr_deprecation_warn <- function() {
  deprecation_warning(old = "antibiotics", new = "antimicrobials", is_dataset = TRUE)
  invisible(NULL)
}
#' @export
`[.deprecated_amr_dataset` <- function(x, ...) {
  .amr_deprecation_warn()
  NextMethod("[")
}
#' @export
`[[.deprecated_amr_dataset` <- function(x, ...) {
  .amr_deprecation_warn()
  NextMethod("[[")
}
#' @export
`$.deprecated_amr_dataset` <- function(x, name) {
  .amr_deprecation_warn()
  NextMethod("$")
}
#' @export
print.deprecated_amr_dataset <- function(x, ...) {
  .amr_deprecation_warn()
  NextMethod("print")
}
#' @export
as.data.frame.deprecated_amr_dataset <- function(x, ...) {
  .amr_deprecation_warn()
  NextMethod("as.data.frame")
}

# REMEMBER to search for `deprecation_warning` in the package code to find all instances.
# currently deprecated arguments at least:
# - `antibiotics` in `antibiogram()`
# - `converse_capped_values` in `as.sir()`

#' @rdname AMR-deprecated
#' @export
ab_class <- function(...) {
  deprecation_warning("ab_class", "amr_class", is_function = TRUE)
  amr_class(...)
}

#' @rdname AMR-deprecated
#' @export
ab_selector <- function(...) {
  deprecation_warning("ab_selector", "amr_selector", is_function = TRUE)
  amr_selector(...)
}

## Helper function ----

deprecation_warning <- function(old = NULL, new = NULL, fn = NULL, extra_msg = NULL, is_function = FALSE, is_dataset = FALSE, is_argument = FALSE) {
  if (is.null(old)) {
    warning_(extra_msg)
  } else if (message_not_thrown_before("deprecation", old, new, entire_session = TRUE)) {
    env <- paste0("deprecated_", old)
    if (!env %in% names(AMR_env)) {
      AMR_env[[paste0("deprecated_", old)]] <- 1
      if (isTRUE(is_function)) {
        old <- paste0(old, "()")
        if (!is.null(new)) {
          new <- paste0(new, "()")
        }
        type <- "function"
      } else if (isTRUE(is_dataset)) {
        type <- "dataset"
      } else if (isTRUE(is_argument)) {
        type <- "argument"
        if (is.null(fn)) {
          stop("Set 'fn' in deprecation_warning()")
        }
      } else {
        stop("Set either 'is_function', 'is_dataset', or 'is_argument' to TRUE in deprecation_warning()")
      }
      warning_(
        ifelse(is.null(new),
          paste0("The `", old, "` ", type, " is deprecated"),
          ifelse(type == "dataset",
            paste0("The `", old, "` ", type, " has been renamed to `", new, "`"),
            ifelse(type == "argument",
              paste0("The `", old, "` ", type, " in `", fn, "()` has been replaced with `", new, "`: `", fn, "(", new, " = ...)`"),
              paste0("The `", old, "` ", type, " has been replaced with `", new, "`")
            )
          )
        ),
        ifelse(type == "dataset",
          ". The old name will be removed in future version, so please update your code.",
          ifelse(type == "argument",
            ". While the old argument still works, it will be removed in a future version, so please update your code.",
            " and will be removed in a future version, see `?AMR-deprecated`."
          )
        ),
        ifelse(!is.null(extra_msg),
          paste0(" ", extra_msg),
          ""
        ),
        "\nThis warning will be shown once per session."
      )
    }
  }
}
