# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Analysis                              #
#                                                                      #
# SOURCE                                                               #
# https://gitlab.com/msberends/AMR                                     #
#                                                                      #
# LICENCE                                                              #
# (c) 2019 Berends MS (m.s.berends@umcg.nl), Luz CF (c.f.luz@umcg.nl)  #
#                                                                      #
# This R package is free software; you can freely use and distribute   #
# it for both personal and commercial purposes under the terms of the  #
# GNU General Public License version 2.0 (GNU GPL-2), as published by  #
# the Free Software Foundation.                                        #
#                                                                      #
# This R package was created for academic research and was publicly    #
# released in the hope that it will be useful, but it comes WITHOUT    #
# ANY WARRANTY OR LIABILITY.                                           #
# Visit our website for more info: https://msberends.gitab.io/AMR.     #
# ==================================================================== #

#' Deprecated functions
#'
#' These functions are so-called '\link{Deprecated}'. They will be removed in a future release. Using the functions will give a warning with the name of the function it has been replaced by (if there is one).
#' @inheritSection AMR Read more on our website!
#' @export
#' @keywords internal
#' @name AMR-deprecated
#' @rdname AMR-deprecated
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

#' @rdname AMR-deprecated
#' @export
guess_mo <- function(...) {
  .Deprecated(new = "as.mo", package = "AMR")
  as.mo(...)
}

#' @rdname AMR-deprecated
#' @export
ab_property <- function(...) {
  .Deprecated(new = "atc_property", package = "AMR")
  atc_property(...)
}

#' @rdname AMR-deprecated
#' @export
ab_atc <- function(...) {
  .Deprecated(new = "as.atc", package = "AMR")
  as.atc(...)
}

#' @rdname AMR-deprecated
#' @export
ab_official <- function(...) {
  .Deprecated(new = "atc_official", package = "AMR")
  atc_official(...)
}

#' @rdname AMR-deprecated
#' @export
ab_name <- function(...) {
  .Deprecated(new = "atc_name", package = "AMR")
  atc_name(...)
}

#' @rdname AMR-deprecated
#' @export
ab_trivial_nl <- function(...) {
  .Deprecated(new = "atc_trivial_nl", package = "AMR")
  atc_trivial_nl(...)
}

#' @rdname AMR-deprecated
#' @export
ab_certe <- function(...) {
  .Deprecated(new = "atc_certe", package = "AMR")
  atc_certe(...)
}

#' @rdname AMR-deprecated
#' @export
ab_umcg <- function(...) {
  .Deprecated(new = "atc_umcg", package = "AMR")
  atc_umcg(...)
}

#' @rdname AMR-deprecated
#' @export
ab_tradenames <- function(...) {
  .Deprecated(new = "atc_tradenames", package = "AMR")
  atc_tradenames(...)
}
