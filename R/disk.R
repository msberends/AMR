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
#' This transforms a vector to a new class [`disk`], which is a growth zone size (around an antibiotic disk) in millimetres between 6 and 50.
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
#' # transform existing disk zones to the `disk` class
#' library(dplyr)
#' df <- data.frame(microorganism = "E. coli",
#'                  AMP = 20,
#'                  CIP = 14,
#'                  GEN = 18,
#'                  TOB = 16)
#' df <- df %>% mutate_at(vars(AMP:TOB), as.disk)
#' df
#' 
#' # interpret disk values, see ?as.rsi
#' as.rsi(x = as.disk(18),
#'        mo = "Strep pneu",  # `mo` will be coerced with as.mo()
#'        ab = "ampicillin",  # and `ab` with as.ab()
#'        guideline = "EUCAST")
#'        
#' as.rsi(df)
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

    class(x) <- "disk"
    x
  }
}

all_valid_disks <- function(x) {
  x_disk <- suppressWarnings(as.disk(x[!is.na(x)]))
  !any(is.na(x_disk)) & !all(is.na(x))
}

#' @rdname as.disk
#' @export
#' @importFrom dplyr %>%
is.disk <- function(x) {
  inherits(x, "disk")
}

#' @exportMethod as.data.frame.disk
#' @export
#' @noRd
as.data.frame.disk <- function(x, ...) {
  # same as as.data.frame.integer but with removed stringsAsFactors, since it will be class "disk"
  nm <- paste(deparse(substitute(x), width.cutoff = 500L),
              collapse = " ")
  if (!"nm" %in% names(list(...))) {
    as.data.frame.vector(x, ..., nm = nm)
  } else {
    as.data.frame.vector(x, ...)
  }
}

#' @exportMethod print.disk
#' @export
#' @noRd
print.disk <- function(x, ...) {
  cat("Class 'disk'\n")
  print(as.integer(x), quote = FALSE)
}

#' @importFrom pillar pillar_shaft
#' @export
pillar_shaft.disk <- function(x, ...) {
  out <- trimws(format(x))
  out[is.na(x)] <- pillar::style_na(NA)
  pillar::new_pillar_shaft_simple(out, align = "right", min_width = 3)
}

#' @importFrom vctrs vec_ptype_abbr
#' @export
vec_ptype_abbr.disk <- function(x, ...) {
  "disk"
}

#' @importFrom vctrs vec_ptype_full
#' @export
vec_ptype_full.disk <- function(x, ...) {
  "disk"
}

#' @exportMethod [.disk
#' @export
#' @noRd
"[.disk" <- function(x, ...) {
  y <- NextMethod()
  attributes(y) <- attributes(x)
  y
}
#' @exportMethod [[.disk
#' @export
#' @noRd
"[[.disk" <- function(x, ...) {
  y <- NextMethod()
  attributes(y) <- attributes(x)
  y
}
#' @exportMethod [<-.disk
#' @export
#' @noRd
"[<-.disk" <- function(i, j, ..., value) {
  y <- NextMethod()
  attributes(y) <- attributes(i)
  y
}
#' @exportMethod [[<-.disk
#' @export
#' @noRd
"[[<-.disk" <- function(i, j, ..., value) {
  y <- NextMethod()
  attributes(y) <- attributes(i)
  y
}
#' @exportMethod c.disk
#' @export
#' @noRd
c.disk <- function(x, ...) {
  y <- NextMethod()
  attributes(y) <- attributes(x)
  y
}
