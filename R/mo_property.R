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
#' Use these functions to return a specific property of a microorganism from the \code{\link{microorganisms}} data set. All input values will be evaluated internally with \code{\link{as.mo}}.
#' @param x any (vector of) text that can be coerced to a valid microorganism code with \code{\link{as.mo}}
#' @param property one of the column names of one of the \code{\link{microorganisms}} data set, like \code{"mo"}, \code{"bactsys"}, \code{"family"}, \code{"genus"}, \code{"species"}, \code{"fullname"}, \code{"gramstain"} and \code{"aerobic"}
#' @inheritParams as.mo
#' @param language language of the returned text, either one of \code{"en"} (English), \code{"de"} (German) or \code{"nl"} (Dutch)
#' @source
#' [1] Becker K \emph{et al.} \strong{Coagulase-Negative Staphylococci}. 2014. Clin Microbiol Rev. 27(4): 870–926. \url{https://dx.doi.org/10.1128/CMR.00109-13}
#'
#' [2] Lancefield RC \strong{A serological differentiation of human and other groups of hemolytic streptococci}. 1933. J Exp Med. 57(4): 571–95. \url{https://dx.doi.org/10.1084/jem.57.4.571}
#' @rdname mo_property
#' @export
#' @importFrom dplyr %>% left_join pull
#' @seealso \code{\link{microorganisms}}
#' @examples
#' # All properties
#' mo_family("E. coli")          # "Enterobacteriaceae"
#' mo_genus("E. coli")           # "Escherichia"
#' mo_species("E. coli")         # "coli"
#' mo_subspecies("E. coli")      # <NA>
#' mo_fullname("E. coli")        # "Escherichia coli"
#' mo_type("E. coli")            # "Bacteria"
#' mo_gramstain("E. coli")       # "Negative rods"
#' mo_aerobic("E. coli")         # TRUE
#'
#' # language support for Spanish, German and Dutch
#' mo_type("E. coli", "es")      # "Bakteria"
#' mo_type("E. coli", "de")      # "Bakterien"
#' mo_type("E. coli", "nl")      # "Bacterie"
#' mo_gramstain("E. coli", "es") # "Bacilos negativos"
#' mo_gramstain("E. coli", "de") # "Negative Staebchen"
#' mo_gramstain("E. coli", "nl") # "Negatieve staven"
#'
#'
#' # Abbreviations known in the field
#' mo_genus("MRSA")              # "Staphylococcus"
#' mo_species("MRSA")            # "aureus"
#' mo_gramstain("MRSA")          # "Positive cocci"
#'
#' mo_genus("VISA")              # "Staphylococcus"
#' mo_species("VISA")            # "aureus"
#'
#'
#' # Known subspecies
#' mo_genus("EHEC")              # "Escherichia"
#' mo_species("EHEC")            # "coli"
#' mo_subspecies("EHEC")         # "EHEC"
#' mo_fullname("EHEC")           # "Escherichia coli (EHEC)"
#'
#' mo_genus("doylei")            # "Campylobacter"
#' mo_species("doylei")          # "jejuni"
#' mo_fullname("doylei")         # "Campylobacter jejuni (doylei)"
#'
#' mo_fullname("K. pneu rh")     # "Klebsiella pneumoniae (rhinoscleromatis)"
#'
#'
#' # Anaerobic bacteria
#' mo_genus("B. fragilis")       # "Bacteroides"
#' mo_species("B. fragilis")     # "fragilis"
#' mo_aerobic("B. fragilis")     # FALSE
#'
#'
#' # Becker classification, see ?as.mo
#' mo_fullname("S. epidermidis")                 # "Staphylococcus epidermidis"
#' mo_fullname("S. epidermidis", Becker = TRUE)  # "Coagulase Negative Staphylococcus (CoNS)"
#'
#' # Lancefield classification, see ?as.mo
#' mo_fullname("S. pyogenes")                    # "Streptococcus pyogenes"
#' mo_fullname("S. pyogenes", Lancefield = TRUE) # "Streptococcus group A"
mo_property <- function(x, property = 'fullname', Becker = FALSE, Lancefield = FALSE) {
  property <- tolower(property[1])
  if (!property %in% colnames(microorganisms)) {
    stop("invalid property: ", property, " - use a column name of the `microorganisms` data set")
  }
  x <- as.mo(x = x, Becker = Becker, Lancefield = Lancefield) # this will give a warning if x cannot be coerced
  suppressWarnings(
    data.frame(mo = x, stringsAsFactors = FALSE) %>%
      left_join(AMR::microorganisms, by = "mo") %>%
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
mo_species <- function(x, Becker = FALSE, Lancefield = FALSE) {
  mo_property(x, "species", Becker = Becker, Lancefield = Lancefield)
}

#' @rdname mo_property
#' @export
mo_subspecies <- function(x, Becker = FALSE, Lancefield = FALSE) {
  mo_property(x, "subspecies", Becker = Becker, Lancefield = Lancefield)
}

#' @rdname mo_property
#' @export
mo_fullname <- function(x, Becker = FALSE, Lancefield = FALSE) {
  mo_property(x, "fullname", Becker = Becker, Lancefield = Lancefield)
}

#' @rdname mo_property
#' @export
mo_type <- function(x, language = "en") {
  mo_property(x, paste0("type", checklang(language)))
}

#' @rdname mo_property
#' @export
mo_gramstain <- function(x, language = "en") {
  mo_property(x, paste0("gramstain", checklang(language)))
}

#' @rdname mo_property
#' @export
mo_aerobic <- function(x) {
  mo_property(x, "aerobic")
}

checklang <- function(language) {
  language <- tolower(language[1])
  supported <- c("en", "de", "nl", "es")
  if (!language %in% c(NULL, "", supported)) {
    stop("invalid language: ", language, " - use one of ", paste0("'", sort(supported), "'", collapse = ", "), call. = FALSE)
  }
  if (language %in% c(NULL, "", "en")) {
    ""
  } else {
    paste0("_", language)
  }
}
