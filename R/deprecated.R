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
# Visit our website for more info: https://msberends.gitlab.io/AMR.    #
# ==================================================================== #

#' Deprecated functions
#'
#' These functions are so-called '\link{Deprecated}'. They will be removed in a future release. Using the functions will give a warning with the name of the function it has been replaced by (if there is one).
#' @inheritSection AMR Read more on our website!
#' @export
#' @keywords internal
#' @name AMR-deprecated
#' @rdname AMR-deprecated
p.symbol <- function(...) {
  .Deprecated("p_symbol", package = "AMR")
  AMR::p_symbol(...)
}

#' @rdname AMR-deprecated
#' @export
portion_R <- function(...) {
  .Deprecated("resistance", package = "AMR")
  proportion_R(...)
}

#' @rdname AMR-deprecated
#' @export
portion_IR <- function(...) {
  .Deprecated("proportion_IR", package = "AMR")
  proportion_IR(...)
}

#' @rdname AMR-deprecated
#' @export
portion_I <- function(...) {
  .Deprecated("proportion_I", package = "AMR")
  proportion_I(...)
}

#' @rdname AMR-deprecated
#' @export
portion_SI <- function(...) {
  .Deprecated("susceptibility", package = "AMR")
  proportion_SI(...)
}

#' @rdname AMR-deprecated
#' @export
portion_S <- function(...) {
  .Deprecated("proportion_S", package = "AMR")
  proportion_S(...)
}
