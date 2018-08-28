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

#' Property of a microorganism
#'
#' Use these functions to return a specific property of a microorganism from the \code{\link{microorganisms}} data set, based on their \code{bactid}. Get such an ID with \code{\link{as.bactid}}.
#' @param x a (vector of a) valid \code{\link{bactid}} or any text that can be coerced to a valid bactid with \code{\link{as.bactid}}
#' @param property one of the column names of one of the \code{\link{microorganisms}} data set, like \code{"bactid"}, \code{"bactsys"}, \code{"family"}, \code{"genus"}, \code{"species"}, \code{"fullname"}, \code{"gramstain"} and \code{"aerobic"}
#' @rdname mo_property
#' @export
#' @importFrom dplyr %>% left_join pull
#' @seealso \code{\link{microorganisms}}
#' @examples
#' # All properties
#' mo_family("E. coli")       # Enterobacteriaceae
#' mo_genus("E. coli")        # Escherichia
#' mo_species("E. coli")      # coli
#' mo_subspecies("E. coli")   # <NA>
#' mo_fullname("E. coli")     # Escherichia coli
#' mo_type("E. coli")         # Bacteria
#' mo_gramstain("E. coli")    # Negative rods
#' mo_aerobic("E. coli")      # TRUE
#' mo_type_nl("E. coli")      # Bacterie
#' mo_gramstain_nl("E. coli") # Negatieve staven
#'
#'
#' # Abbreviations known in the field
#' mo_genus("EHEC")           # Escherichia
#' mo_species("EHEC")         # coli
#' mo_subspecies("EHEC")      # EHEC
#' mo_fullname("EHEC")        # Escherichia coli (EHEC)
#'
#' mo_genus("MRSA")           # Staphylococcus
#' mo_species("MRSA")         # aureus
#' mo_gramstain("MRSA")       # Positive cocci
#'
#' mo_genus("VISA")           # Staphylococcus
#' mo_species("VISA")         # aureus
#'
#'
#' # Known subspecies
#' mo_genus("doylei")         # Campylobacter
#' mo_species("doylei")       # jejuni
#' mo_fullname("doylei")      # Campylobacter jejuni (doylei)
#'
#'
#' # Anaerobic bacteria
#' mo_genus("B. fragilis")    # Bacteroides
#' mo_species("B. fragilis")  # fragilis
#' mo_aerobic("B. fragilis")  # FALSE
mo_property <- function(x, property = 'fullname') {
  property <- property[1]
  if (!property %in% colnames(microorganisms)) {
    stop("invalid property: ", property, " - use a column name of `microorganisms`")
  }
  if (!is.bactid(x)) {
    x <- as.bactid(x) # this will give a warning if x cannot be coerced
  }
  suppressWarnings(
    data.frame(bactid = x, stringsAsFactors = FALSE) %>%
      left_join(AMR::microorganisms, by = "bactid") %>%
      pull(property)
  )
}

#' @rdname mo_property
#' @export
mo_family <- function(x) {
  mo_property(x, "family")
}

#' @rdname mo_property
#' @export
mo_genus <- function(x) {
  mo_property(x, "genus")
}

#' @rdname mo_property
#' @export
mo_species <- function(x) {
  mo_property(x, "species")
}

#' @rdname mo_property
#' @export
mo_subspecies <- function(x) {
  mo_property(x, "subspecies")
}

#' @rdname mo_property
#' @export
mo_fullname <- function(x) {
  mo_property(x, "fullname")
}

#' @rdname mo_property
#' @export
mo_type <- function(x) {
  mo_property(x, "type")
}

#' @rdname mo_property
#' @export
mo_gramstain <- function(x) {
  mo_property(x, "gramstain")
}

#' @rdname mo_property
#' @export
mo_aerobic <- function(x) {
  mo_property(x, "aerobic")
}

#' @rdname mo_property
#' @export
mo_type_nl <- function(x) {
  mo_property(x, "type_nl")
}

#' @rdname mo_property
#' @export
mo_gramstain_nl <- function(x) {
  mo_property(x, "gramstain_nl")
}
