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

#' Join [microorganisms] to a Data Set
#'
#' Join the data set [microorganisms] easily to an existing data set or to a [character] vector.
#' @rdname join
#' @name join
#' @aliases join inner_join
#' @param x existing data set to join, or [character] vector. In case of a [character] vector, the resulting [data.frame] will contain a column 'x' with these values.
#' @param by a variable to join by - if left empty will search for a column with class [`mo`] (created with [as.mo()]) or will be `"mo"` if that column name exists in `x`, could otherwise be a column name of `x` with values that exist in `microorganisms$mo` (such as `by = "bacteria_id"`), or another column in [microorganisms] (but then it should be named, like `by = c("bacteria_id" = "fullname")`)
#' @param suffix if there are non-joined duplicate variables in `x` and `y`, these suffixes will be added to the output to disambiguate them. Should be a [character] vector of length 2.
#' @param ... ignored, only in place to allow future extensions
#' @details **Note:** As opposed to the `join()` functions of `dplyr`, [character] vectors are supported and at default existing columns will get a suffix `"2"` and the newly joined columns will not get a suffix.
#'
#' If the `dplyr` package is installed, their join functions will be used. Otherwise, the much slower [merge()] and [interaction()] functions from base \R will be used.
#' @return a [data.frame]
#' @export
#' @examples
#' left_join_microorganisms(as.mo("K. pneumoniae"))
#' left_join_microorganisms("B_KLBSL_PNMN")
#'
#' df <- data.frame(
#'   date = seq(
#'     from = as.Date("2018-01-01"),
#'     to = as.Date("2018-01-07"),
#'     by = 1
#'   ),
#'   bacteria = as.mo(c(
#'     "S. aureus", "MRSA", "MSSA", "STAAUR",
#'     "E. coli", "E. coli", "E. coli"
#'   )),
#'   stringsAsFactors = FALSE
#' )
#' colnames(df)
#'
#' df_joined <- left_join_microorganisms(df, "bacteria")
#' colnames(df_joined)
#'
#' \donttest{
#' if (require("dplyr")) {
#'   example_isolates %>%
#'     left_join_microorganisms() %>%
#'     colnames()
#' }
#' }
inner_join_microorganisms <- function(x, by = NULL, suffix = c("2", ""), ...) {
  meet_criteria(x, allow_class = c("data.frame", "character"))
  meet_criteria(by, allow_class = "character", allow_NULL = TRUE)
  meet_criteria(suffix, allow_class = "character", has_length = 2)

  join_microorganisms(type = "inner_join", x = x, by = by, suffix = suffix, ...)
}

#' @rdname join
#' @export
left_join_microorganisms <- function(x, by = NULL, suffix = c("2", ""), ...) {
  meet_criteria(x, allow_class = c("data.frame", "character"))
  meet_criteria(by, allow_class = "character", allow_NULL = TRUE)
  meet_criteria(suffix, allow_class = "character", has_length = 2)

  join_microorganisms(type = "left_join", x = x, by = by, suffix = suffix, ...)
}

#' @rdname join
#' @export
right_join_microorganisms <- function(x, by = NULL, suffix = c("2", ""), ...) {
  meet_criteria(x, allow_class = c("data.frame", "character"))
  meet_criteria(by, allow_class = "character", allow_NULL = TRUE)
  meet_criteria(suffix, allow_class = "character", has_length = 2)

  join_microorganisms(type = "right_join", x = x, by = by, suffix = suffix, ...)
}

#' @rdname join
#' @export
full_join_microorganisms <- function(x, by = NULL, suffix = c("2", ""), ...) {
  meet_criteria(x, allow_class = c("data.frame", "character"))
  meet_criteria(by, allow_class = "character", allow_NULL = TRUE)
  meet_criteria(suffix, allow_class = "character", has_length = 2)

  join_microorganisms(type = "full_join", x = x, by = by, suffix = suffix, ...)
}

#' @rdname join
#' @export
semi_join_microorganisms <- function(x, by = NULL, ...) {
  meet_criteria(x, allow_class = c("data.frame", "character"))
  meet_criteria(by, allow_class = "character", allow_NULL = TRUE)

  join_microorganisms(type = "semi_join", x = x, by = by, ...)
}

#' @rdname join
#' @export
anti_join_microorganisms <- function(x, by = NULL, ...) {
  meet_criteria(x, allow_class = c("data.frame", "character"))
  meet_criteria(by, allow_class = "character", allow_NULL = TRUE)

  join_microorganisms(type = "anti_join", x = x, by = by, ...)
}

join_microorganisms <- function(type, x, by, suffix, ...) {
  if (!is.data.frame(x)) {
    if (pkg_is_available("tibble", also_load = FALSE)) {
      x <- import_fn("tibble", "tibble")(mo = x)
    } else {
      x <- data.frame(mo = x, stringsAsFactors = FALSE)
    }
    by <- "mo"
  }
  x.bak <- x
  if (is.null(by)) {
    by <- search_type_in_df(x, "mo", info = FALSE)
    if (is.null(by) && NCOL(x) == 1) {
      by <- colnames(x)[1L]
    } else {
      stop_if(is.null(by), "no column with microorganism names or codes found, set this column with `by`", call = -2)
    }
    message_('Joining, by = "', by, '"', add_fn = font_black, as_note = FALSE) # message same as dplyr::join functions
  }
  if (!all(x[, by, drop = TRUE] %in% AMR_env$MO_lookup$mo, na.rm = TRUE)) {
    x$join.mo <- as.mo(x[, by, drop = TRUE])
    by <- c("join.mo" = "mo")
  } else {
    x[, by] <- as.mo(x[, by, drop = TRUE])
  }

  if (is.null(names(by))) {
    # will always be joined to microorganisms$mo, so add name to that
    by <- stats::setNames("mo", by)
  }

  # use dplyr if available - it's much faster than poorman alternatives
  dplyr_join <- import_fn(name = type, pkg = "dplyr", error_on_fail = FALSE)
  if (!is.null(dplyr_join)) {
    join_fn <- dplyr_join
  } else {
    # otherwise use poorman, see R/aa_helper_pm_functions.R
    join_fn <- get(paste0("pm_", type), envir = asNamespace("AMR"))
  }
  if (type %like% "full|left|right|inner") {
    joined <- join_fn(x = x, y = AMR::microorganisms, by = by, suffix = suffix, ...)
  } else {
    joined <- join_fn(x = x, y = AMR::microorganisms, by = by, ...)
  }

  if ("join.mo" %in% colnames(joined)) {
    if ("mo" %in% colnames(joined)) {
      ind_mo <- which(colnames(joined) %in% c("mo", "join.mo"))
      colnames(joined)[ind_mo[1L]] <- paste0("mo", suffix[1L])
      colnames(joined)[ind_mo[2L]] <- paste0("mo", suffix[2L])
    } else {
      colnames(joined)[colnames(joined) == "join.mo"] <- "mo"
    }
  }

  if (type %like% "full|left|right|inner" && NROW(joined) > NROW(x)) {
    warning_("in `", type, "_microorganisms()`: the newly joined data set contains ", nrow(joined) - nrow(x), " rows more than the number of rows of `x`.")
  }

  as_original_data_class(joined, class(x.bak))
}
