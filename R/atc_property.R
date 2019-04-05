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

#' Property of an antibiotic
#'
#' Use these functions to return a specific property of an antibiotic from the \code{\link{antibiotics}} data set, based on their ATC code. Get such a code with \code{\link{as.atc}}.
#' @param x a (vector of a) valid \code{\link{atc}} code or any text that can be coerced to a valid atc with \code{\link{as.atc}}
#' @param property one of the column names of one of the \code{\link{antibiotics}} data set, like \code{"atc"} and \code{"official"}
#' @param language language of the returned text, defaults to English (\code{"en"}) and can be set with \code{\link{getOption}("AMR_locale")}. Either one of \code{"en"} (English) or \code{"nl"} (Dutch).
#' @rdname atc_property
#' @return A vector of values. In case of \code{atc_tradenames}, if \code{x} is of length one, a vector will be returned. Otherwise a \code{\link{list}}, with \code{x} as names.
#' @export
#' @importFrom dplyr %>% left_join pull
#' @seealso \code{\link{antibiotics}}
#' @inheritSection AMR Read more on our website!
#' @examples
#' as.atc("amcl")         # J01CR02
#' atc_name("amcl")        # Amoxicillin and beta-lactamase inhibitor
#' atc_name("amcl", "nl")  # Amoxicilline met enzymremmer
#' atc_trivial_nl("amcl")  # Amoxicilline/clavulaanzuur
#' atc_certe("amcl")       # amcl
#' atc_umcg("amcl")        # AMCL
atc_property <- function(x, property = 'official') {
  property <- property[1]
  if (!property %in% colnames(AMR::antibiotics)) {
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

#' @rdname atc_property
#' @export
atc_official <- function(x, language = NULL) {

  if (is.null(language)) {
    language <- getOption("AMR_locale", default = "en")[1L]
  } else {
    language <- tolower(language[1])
  }
  if (language %in% c("en", "")) {
    atc_property(x, "official")
  } else if (language == "nl") {
    atc_property(x, "official_nl")
  } else {
    stop("Unsupported language: '", language, "' - use one of: 'en', 'nl'", call. = FALSE)
  }
}

#' @rdname atc_property
#' @export
atc_name <- atc_official

#' @rdname atc_property
#' @export
atc_trivial_nl <- function(x) {
  atc_property(x, "trivial_nl")
}

#' @rdname atc_property
#' @export
atc_certe <- function(x) {
  atc_property(x, "certe")
}

#' @rdname atc_property
#' @export
atc_umcg <- function(x) {
  atc_property(x, "umcg")
}

#' @rdname atc_property
#' @export
atc_tradenames <- function(x) {
  res <- atc_property(x, "trade_name")
  res <- strsplit(res, "|", fixed = TRUE)
  if (length(x) == 1) {
    res <- unlist(res)
  } else {
    names(res) <- x
  }
  res
}
