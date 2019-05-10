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

#' Transform to ATC code
#'
#' Use this function to determine the ATC code of one or more antibiotics. The data set \code{\link{antibiotics}} will be searched for abbreviations, official names and trade names.
#' @param x character vector to determine \code{ATC} code
#' @rdname as.atc
#' @aliases atc
#' @keywords atc
#' @inheritSection WHOCC WHOCC
#' @export
#' @importFrom dplyr %>% filter slice pull
#' @details Use the \code{\link{ab_property}} functions to get properties based on the returned ATC code, see Examples.
#'
#' In the ATC classification system, the active substances are classified in a hierarchy with five different levels.  The system has fourteen main anatomical/pharmacological groups or 1st levels. Each ATC main group is divided into 2nd levels which could be either pharmacological or therapeutic groups.  The 3rd and 4th levels are chemical, pharmacological or therapeutic subgroups and the 5th level is the chemical substance.  The 2nd, 3rd and 4th levels are often used to identify pharmacological subgroups when that is considered more appropriate than therapeutic or chemical subgroups.
#'   Source: \url{https://www.whocc.no/atc/structure_and_principles/}
#' @return Character (vector) with class \code{"act"}. Unknown values will return \code{NA}.
#' @seealso \code{\link{antibiotics}} for the dataframe that is being used to determine ATCs.
#' @inheritSection AMR Read more on our website!
#' @examples
#' # These examples all return "J01FA01", the ATC code of Erythromycin:
#' as.atc("J01FA01")
#' as.atc("Erythromycin")
#' as.atc("eryt")
#' as.atc("   eryt 123")
#' as.atc("ERYT")
#' as.atc("ERY")
as.atc <- function(x) {
  ab_atc(x)
}

#' @rdname as.atc
#' @export
is.atc <- function(x) {
  identical(class(x), "atc")
}

#' @exportMethod print.atc
#' @export
#' @noRd
print.atc <- function(x, ...) {
  cat("Class 'atc'\n")
  print.default(as.character(x), quote = FALSE)
}

#' @exportMethod as.data.frame.atc
#' @export
#' @noRd
as.data.frame.atc <- function (x, ...) {
  # same as as.data.frame.character but with removed stringsAsFactors
  nm <- paste(deparse(substitute(x), width.cutoff = 500L),
              collapse = " ")
  if (!"nm" %in% names(list(...))) {
    as.data.frame.vector(x, ..., nm = nm)
  } else {
    as.data.frame.vector(x, ...)
  }
}

#' @exportMethod pull.atc
#' @export
#' @importFrom dplyr pull
#' @noRd
pull.atc <- function(.data, ...) {
  pull(as.data.frame(.data), ...)
}
