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

#' Class 'disk'
#'
#' This transforms a vector to a new class [`disk`], which is a growth zone size (around an antibiotic disk) in millimeters between 6 and 50.
#' @inheritSection lifecycle Stable lifecycle
#' @rdname as.disk
#' @param x vector
#' @param na.rm a logical indicating whether missing values should be removed
#' @details Interpret disk values as RSI values with [as.rsi()]. It supports guidelines from EUCAST and CLSI.
#' @return An [`integer`] with additional new class [`disk`]
#' @aliases disk
#' @export
#' @seealso [as.rsi()]
#' @inheritSection AMR Read more on our website!
#' @examples
#' # interpret disk values
#' as.rsi(x = 12,
#'        mo = as.mo("S. pneumoniae"),
#'        ab = "AMX",
#'        guideline = "EUCAST")
#' as.rsi(x = 12,
#'        mo = as.mo("S. pneumoniae"),
#'        ab = "AMX",
#'        guideline = "CLSI")
as.disk <- function(x, na.rm = FALSE) {
  if (is.disk(x)) {
    x
  } else {
    x <- x %>% unlist()
    if (na.rm == TRUE) {
      x <- x[!is.na(x)]
    }
    x.bak <- x

    na_before <- length(x[is.na(x)])

    # force it to be integer
    x <- suppressWarnings(as.integer(x))

    # disks can never be less than 6 mm (size of smallest disk) or more than 50 mm
    x[x < 6 | x > 50] <- NA_integer_
    na_after <- length(x[is.na(x)])

    if (na_before != na_after) {
      list_missing <- x.bak[is.na(x) & !is.na(x.bak)] %>%
        unique() %>%
        sort()
      list_missing <- paste0('"', list_missing, '"', collapse = ", ")
      warning(na_after - na_before, " results truncated (",
              round(((na_after - na_before) / length(x)) * 100),
              "%) that were invalid disk zones: ",
              list_missing, call. = FALSE)
    }

    class(x) <- c("disk", "integer")
    x
  }
}

#' @rdname as.disk
#' @export
#' @importFrom dplyr %>%
is.disk <- function(x) {
  inherits(x, "disk")
}

#' @exportMethod print.disk
#' @export
#' @noRd
print.disk <- function(x, ...) {
  cat("Class 'disk'\n")
  print(as.integer(x), quote = FALSE)
}

#' @importFrom pillar type_sum
#' @export
type_sum.disk <- function(x) {
  "disk"
}

#' @importFrom pillar pillar_shaft
#' @export
pillar_shaft.disk <- function(x, ...) {
  out <- trimws(format(x))
  out[is.na(x)] <- pillar::style_na(NA)
  pillar::new_pillar_shaft_simple(out, align = "right", min_width = 3)
}
