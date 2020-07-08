# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Analysis                              #
#                                                                      #
# SOURCE                                                               #
# https://github.com/msberends/AMR                                     #
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
# Visit our website for more info: https://msberends.github.io/AMR.    #
# ==================================================================== #

#' Deprecated functions
#'
#' These functions are so-called '[Deprecated]'. They will be removed in a future release. Using the functions will give a warning with the name of the function it has been replaced by (if there is one).
#' @inheritSection lifecycle Retired lifecycle
#' @inheritSection AMR Read more on our website!
#' @export
#' @keywords internal
#' @name AMR-deprecated
#' @export
portion_R <- function(...) {
  .Deprecated("resistance()", package = "AMR")
  proportion_R(...)
}

#' @rdname AMR-deprecated
#' @export
portion_IR <- function(...) {
  .Deprecated("proportion_IR()", package = "AMR")
  proportion_IR(...)
}

#' @rdname AMR-deprecated
#' @export
portion_I <- function(...) {
  .Deprecated("proportion_I()", package = "AMR")
  proportion_I(...)
}

#' @rdname AMR-deprecated
#' @export
portion_SI <- function(...) {
  .Deprecated("susceptibility()", package = "AMR")
  proportion_SI(...)
}

#' @rdname AMR-deprecated
#' @export
portion_S <- function(...) {
  .Deprecated("proportion_S()", package = "AMR")
  proportion_S(...)
}

#' @rdname AMR-deprecated
#' @export
portion_df <- function(...) {
  .Deprecated("proportion_df()", package = "AMR")
  proportion_df(...)
}
