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

#' Join a table with [microorganisms]
#'
#' Join the data set [microorganisms] easily to an existing table or character vector.
#' @inheritSection lifecycle Stable lifecycle
#' @rdname join
#' @name join
#' @aliases join inner_join
#' @param x existing table to join, or character vector
#' @param by a variable to join by - if left empty will search for a column with class [`mo`] (created with [as.mo()]) or will be `"mo"` if that column name exists in `x`, could otherwise be a column name of `x` with values that exist in `microorganisms$mo` (like `by = "bacteria_id"`), or another column in [microorganisms] (but then it should be named, like `by = c("my_genus_species" = "fullname")`)
#' @param suffix if there are non-joined duplicate variables in `x` and `y`, these suffixes will be added to the output to disambiguate them. Should be a character vector of length 2.
#' @param ... other parameters to pass on to [dplyr::join()]
#' @details **Note:** As opposed to the [dplyr::join()] functions of `dplyr`, [`character`] vectors are supported and at default existing columns will get a suffix `"2"` and the newly joined columns will not get a suffix. See [dplyr::join()] for more information.
#' @inheritSection AMR Read more on our website!
#' @export
#' @examples
#' left_join_microorganisms(as.mo("K. pneumoniae"))
#' left_join_microorganisms("B_KLBSL_PNE")
#'
#' library(dplyr)
#' example_isolates %>% left_join_microorganisms()
#'
#' df <- data.frame(date = seq(from = as.Date("2018-01-01"),
#'                             to = as.Date("2018-01-07"),
#'                             by = 1),
#'                  bacteria = as.mo(c("S. aureus", "MRSA", "MSSA", "STAAUR",
#'                                     "E. coli", "E. coli", "E. coli")),
#'                  stringsAsFactors = FALSE)
#' colnames(df)
#' df_joined <- left_join_microorganisms(df, "bacteria")
#' colnames(df_joined)
inner_join_microorganisms <- function(x, by = NULL, suffix = c("2", ""), ...) {
  checked <- joins_check_df(x, by)
  x <- checked$x
  by <- checked$by
  join <- suppressWarnings(
    dplyr::inner_join(x = x, y = microorganisms, by = by, suffix = suffix, ...)
  )
  if (nrow(join) > nrow(x)) {
    warning("The newly joined tbl contains ", nrow(join) - nrow(x), " rows more that its original.")
  }
  join
}

#' @rdname join
#' @export
left_join_microorganisms <- function(x, by = NULL, suffix = c("2", ""), ...) {
  checked <- joins_check_df(x, by)
  x <- checked$x
  by <- checked$by
  join <- suppressWarnings(
    dplyr::left_join(x = x, y = microorganisms, by = by, suffix = suffix, ...)
  )
  if (nrow(join) > nrow(x)) {
    warning("The newly joined tbl contains ", nrow(join) - nrow(x), " rows more that its original.")
  }
  join
}

#' @rdname join
#' @export
right_join_microorganisms <- function(x, by = NULL, suffix = c("2", ""), ...) {
  checked <- joins_check_df(x, by)
  x <- checked$x
  by <- checked$by
  join <- suppressWarnings(
    dplyr::right_join(x = x, y = microorganisms, by = by, suffix = suffix, ...)
  )
  if (nrow(join) > nrow(x)) {
    warning("The newly joined tbl contains ", nrow(join) - nrow(x), " rows more that its original.")
  }
  join
}

#' @rdname join
#' @export
full_join_microorganisms <- function(x, by = NULL, suffix = c("2", ""), ...) {
  checked <- joins_check_df(x, by)
  x <- checked$x
  by <- checked$by
  join <- suppressWarnings(
    dplyr::full_join(x = x, y = microorganisms, by = by, suffix = suffix, ...)
  )
  if (nrow(join) > nrow(x)) {
    warning("The newly joined tbl contains ", nrow(join) - nrow(x), " rows more that its original.")
  }
  join
}

#' @rdname join
#' @export
semi_join_microorganisms <- function(x, by = NULL, ...) {
  checked <- joins_check_df(x, by)
  x <- checked$x
  by <- checked$by
  suppressWarnings(
    dplyr::semi_join(x = x, y = microorganisms, by = by, ...)
  )
}

#' @rdname join
#' @export
anti_join_microorganisms <- function(x, by = NULL, ...) {
  checked <- joins_check_df(x, by)
  x <- checked$x
  by <- checked$by
  suppressWarnings(
    dplyr::anti_join(x = x, y = microorganisms, by = by, ...)
  )
}

joins_check_df <- function(x, by) {
  if (!any(class(x) %in% c("data.frame", "matrix"))) {
    x <- data.frame(mo = as.character(x), stringsAsFactors = FALSE)
    if (is.null(by)) {
      by <- "mo"
    }
  }
  if (is.null(by)) {
    # search for column with class `mo` and return first one found
    by <- colnames(x)[lapply(x, is.mo) == TRUE][1]
    if (is.na(by)) {
      if ("mo" %in% colnames(x)) {
        by <- "mo"
      } else {
        stop("Cannot join - no column found with name or class  `mo`.", call. = FALSE)
      }
    }
    message('Joining, by = "', by, '"') # message same as dplyr::join functions
  }
  if (is.null(names(by))) {
    joinby <- colnames(microorganisms)[1]
    names(joinby) <- by
  } else {
    joinby <- by
  }
  list(x = x,
       by = joinby)
}
