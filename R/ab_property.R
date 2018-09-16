# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Analysis                              #
#                                                                      #
# AUTHORS                                                              #
# Berends MS (m.s.berends@umcg.nl), Luz CF (c.f.luz@umcg.nl)           #
#                                                                      #
# LICENCE                                                              #
# This program is free software; you can redistribute it and/or modify #
# it under the terms of the GNU General Public License version 2.0,    #
# as published by the Free Software Foundation.                        #
#                                                                      #
# This program is distributed in the hope that it will be useful,      #
# but WITHOUT ANY WARRANTY; without even the implied warranty of       #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        #
# GNU General Public License for more details.                         #
# ==================================================================== #

#' Property of an antibiotic
#'
#' Use these functions to return a specific property of an antibiotic from the \code{\link{antibiotics}} data set, based on their ATC code. Get such a code with \code{\link{as.atc}}.
#' @param x a (vector of a) valid \code{\link{atc}} code or any text that can be coerced to a valid atc with \code{\link{as.atc}}
#' @param property one of the column names of one of the \code{\link{antibiotics}} data set, like \code{"atc"} and \code{"official"}
#' @param language language of the returned text, defaults to the systems language. Either one of \code{"en"} (English) or \code{"nl"} (Dutch).
#' @rdname ab_property
#' @return A vector of values. In case of \code{ab_tradenames}, if \code{x} is of length one, a vector will be returned. Otherwise a \code{\link{list}}, with \code{x} as names.
#' @export
#' @importFrom dplyr %>% left_join pull
#' @seealso \code{\link{antibiotics}}
#' @examples
#' ab_atc("amcl")         # J01CR02
#' ab_name("amcl")        # Amoxicillin and beta-lactamase inhibitor
#' ab_name("amcl", "nl")  # Amoxicilline met enzymremmer
#' ab_trivial_nl("amcl")  # Amoxicilline/clavulaanzuur
#' ab_certe("amcl")       # amcl
#' ab_umcg("amcl")        # AMCL
ab_property <- function(x, property = 'official') {
  property <- property[1]
  if (!property %in% colnames(antibiotics)) {
    stop("invalid property: ", property, " - use a column name of the `antibiotics` data set")
  }
  if (!is.atc(x)) {
    x <- as.atc(x) # this will give a warning if x cannot be coerced
  }
  suppressWarnings(
    data.frame(atc = x, stringsAsFactors = FALSE) %>%
      left_join(AMR::antibiotics, by = "atc") %>%
      pull(property)
  )
}

#' @rdname ab_property
#' @export
ab_atc <- function(x) {
  as.character(as.atc(x))
}

#' @rdname ab_property
#' @export
ab_official <- function(x, language = NULL) {

  if (is.null(language)) {
    language <- Sys.locale()
  } else {
    language <- tolower(language[1])
  }
  if (language %in% c("en", "")) {
    ab_property(x, "official")
  } else if (language == "nl") {
    ab_property(x, "official_nl")
  } else {
    stop("Unsupported language: '", language, "' - use one of: 'en', 'nl'", call. = FALSE)
  }
}

#' @rdname ab_property
#' @export
ab_name <- ab_official

#' @rdname ab_property
#' @export
ab_trivial_nl <- function(x) {
  ab_property(x, "trivial_nl")
}

#' @rdname ab_property
#' @export
ab_certe <- function(x) {
  ab_property(x, "certe")
}

#' @rdname ab_property
#' @export
ab_umcg <- function(x) {
  ab_property(x, "umcg")
}

#' @rdname ab_property
#' @export
ab_tradenames <- function(x) {
  res <- ab_property(x, "trade_name")
  res <- strsplit(res, "|", fixed = TRUE)
  if (length(x) == 1) {
    res <- unlist(res)
  } else {
    names(res) <- x
  }
  res
}
