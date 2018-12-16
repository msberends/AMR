# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Analysis                              #
#                                                                      #
# AUTHORS                                                              #
# Berends MS (m.s.berends@umcg.nl), Luz CF (c.f.luz@umcg.nl)           #
#                                                                      #
# LICENCE                                                              #
# This package is free software; you can redistribute it and/or modify #
# it under the terms of the GNU General Public License version 2.0,    #
# as published by the Free Software Foundation.                        #
#                                                                      #
# This R package is distributed in the hope that it will be useful,    #
# but WITHOUT ANY WARRANTY; without even the implied warranty of       #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        #
# GNU General Public License version 2.0 for more details.             #
# ==================================================================== #

#' Deprecated functions
#'
#' These functions are \link{Deprecated}. They will be removed in a future release. Using the functions will give a warning with the name of the function it has been replaced by.
#' @export
#' @keywords internal
#' @name AMR-deprecated
#' @rdname AMR-deprecated
as.bactid <- function(...) {
  .Deprecated("as.mo", package = "AMR")
  as.mo(...)
}

#' @rdname AMR-deprecated
#' @export
is.bactid <- function(...) {
  .Deprecated(new = "is.mo", package = "AMR")
  is.mo(...)
}

#' @rdname AMR-deprecated
#' @export
guess_bactid <- function(...) {
  .Deprecated(new = "guess_mo", package = "AMR")
  guess_mo(...)
}

#' @exportMethod print.bactid
#' @export
#' @noRd
print.bactid <- function(x, ...) {
  cat("Class 'bactid'\n")
  print.default(as.character(x), quote = FALSE)
}

#' @exportMethod as.data.frame.bactid
#' @export
#' @noRd
as.data.frame.bactid <- function (x, ...) {
  # same as as.data.frame.character but with removed stringsAsFactors
  nm <- paste(deparse(substitute(x), width.cutoff = 500L),
              collapse = " ")
  if (!"nm" %in% names(list(...))) {
    as.data.frame.vector(x, ..., nm = nm)
  } else {
    as.data.frame.vector(x, ...)
  }
}

#' @exportMethod pull.bactid
#' @export
#' @importFrom dplyr pull
#' @noRd
pull.bactid <- function(.data, ...) {
  pull(as.data.frame(.data), ...)
}

#' @rdname AMR-deprecated
#' @export
ratio <- function(x, ratio) {
  .Deprecated(package = "AMR")

  if (!all(is.numeric(x))) {
    stop('`x` must be a vector of numeric values.')
  }
  if (length(ratio) == 1) {
    if (ratio %like% '^([0-9]+([.][0-9]+)?[-,:])+[0-9]+([.][0-9]+)?$') {
      # support for "1:2:1", "1-2-1", "1,2,1" and even "1.75:2:1.5"
      ratio <- ratio %>% strsplit("[-,:]") %>% unlist() %>% as.double()
    } else {
      stop('Invalid `ratio`: ', ratio, '.')
    }
  }
  if (length(x) != 1 & length(x) != length(ratio)) {
    stop('`x` and `ratio` must be of same size.')
  }
  sum(x, na.rm = TRUE) * (ratio / sum(ratio, na.rm = TRUE))
}
