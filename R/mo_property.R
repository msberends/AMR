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
#' @param language language of the returned text, defaults to the systems language. Either one of \code{"en"} (English), \code{"de"} (German), \code{"nl"} (Dutch), \code{"es"} (Spanish) or \code{"pt"} (Portuguese).
#' @source
#' [1] Becker K \emph{et al.} \strong{Coagulase-Negative Staphylococci}. 2014. Clin Microbiol Rev. 27(4): 870–926. \url{https://dx.doi.org/10.1128/CMR.00109-13}
#'
#' [2] Lancefield RC \strong{A serological differentiation of human and other groups of hemolytic streptococci}. 1933. J Exp Med. 57(4): 571–95. \url{https://dx.doi.org/10.1084/jem.57.4.571}
#'
#' [3] Integrated Taxonomic Information System (ITIS) on-line database, \url{https://www.itis.gov}.
#' @rdname mo_property
#' @name mo_property
#' @return A logical (in case of \code{mo_aerobic}), a list (in case of \code{mo_taxonomy}), a character otherwise
#' @export
#' @importFrom dplyr %>% left_join pull
#' @seealso \code{\link{microorganisms}}
#' @examples
#' # All properties
#' mo_phylum("E. coli")          # "Proteobacteria"
#' mo_class("E. coli")           # "Gammaproteobacteria"
#' mo_order("E. coli")           # "Enterobacteriales"
#' mo_family("E. coli")          # "Enterobacteriaceae"
#' mo_genus("E. coli")           # "Escherichia"
#' mo_species("E. coli")         # "coli"
#' mo_subspecies("E. coli")      # ""
#' mo_fullname("E. coli")        # "Escherichia coli"
#' mo_shortname("E. coli")       # "E. coli"
#' mo_type("E. coli")            # "Bacteria"
#' mo_gramstain("E. coli")       # "Negative rods"
#' mo_aerobic("E. coli")         # TRUE
#'
#'
#' # Abbreviations known in the field
#' mo_genus("MRSA")              # "Staphylococcus"
#' mo_species("MRSA")            # "aureus"
#' mo_shortname("MRSA")          # "S. aureus"
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
#' mo_shortname("EHEC")          # "E. coli"
#'
#' mo_genus("doylei")            # "Campylobacter"
#' mo_species("doylei")          # "jejuni"
#' mo_fullname("doylei")         # "Campylobacter jejuni (doylei)"
#'
#' mo_fullname("K. pneu rh")     # "Klebsiella pneumoniae (rhinoscleromatis)"
#' mo_shortname("K. pneu rh")    # "K. pneumoniae"
#'
#'
#' # Anaerobic bacteria
#' mo_genus("B. fragilis")       # "Bacteroides"
#' mo_species("B. fragilis")     # "fragilis"
#' mo_aerobic("B. fragilis")     # FALSE
#'
#'
#' # Becker classification, see ?as.mo
#' mo_fullname("S. epi")                     # "Staphylococcus epidermidis"
#' mo_fullname("S. epi", Becker = TRUE)      # "Coagulase Negative Staphylococcus (CoNS)"
#' mo_shortname("S. epi")                    # "S. epidermidis"
#' mo_shortname("S. epi", Becker = TRUE)     # "CoNS"
#'
#' # Lancefield classification, see ?as.mo
#' mo_fullname("S. pyo")                     # "Streptococcus pyogenes"
#' mo_fullname("S. pyo", Lancefield = TRUE)  # "Streptococcus group A"
#' mo_shortname("S. pyo")                    # "S. pyogenes"
#' mo_shortname("S. pyo", Lancefield = TRUE) # "GAS"
#'
#'
#' # Language support for German, Dutch, Spanish and Portuguese
#' mo_type("E. coli", language = "de")       # "Bakterium"
#' mo_type("E. coli", language = "nl")       # "Bacterie"
#' mo_type("E. coli", language = "es")       # "Bakteria"
#' mo_gramstain("E. coli", language = "de")  # "Negative Staebchen"
#' mo_gramstain("E. coli", language = "nl")  # "Negatieve staven"
#' mo_gramstain("E. coli", language = "es")  # "Bacilos negativos"
#' mo_gramstain("Giardia", language = "pt")  # "Parasitas"
#'
#' mo_fullname("S. pyogenes",
#'             Lancefield = TRUE,
#'             language = "de")              # "Streptococcus Gruppe A"
#' mo_fullname("S. pyogenes",
#'             Lancefield = TRUE,
#'             language = "nl")              # "Streptococcus groep A"
#'
#'
#' # Complete taxonomy up to Phylum, returns a list
#' mo_taxonomy("E. coli")
mo_fullname <- function(x, Becker = FALSE, Lancefield = FALSE, language = NULL) {
  mo_property(x, "fullname", Becker = Becker, Lancefield = Lancefield, language = language)
}

#' @rdname mo_property
#' @export
mo_shortname <- function(x, Becker = FALSE, Lancefield = FALSE, language = NULL) {
  if (Becker %in% c(TRUE, "all") | Lancefield == TRUE) {
    res1 <- as.mo(x)
    res2 <- suppressWarnings(as.mo(x, Becker = Becker, Lancefield = Lancefield))
    res2_fullname <- mo_fullname(res2)
    res2_fullname[res2_fullname %like% "\\(CoNS\\)"] <- "CoNS"
    res2_fullname[res2_fullname %like% "\\(CoPS\\)"] <- "CoPS"
    res2_fullname <- gsub("Streptococcus (group|Gruppe|gruppe|groep|grupo|gruppo|groupe) (.)",
                          "G\\2S",
                          res2_fullname) # turn "Streptococcus group A" and "Streptococcus grupo A" to "GAS"
    res2_fullname_vector <- res2_fullname[res2_fullname == mo_fullname(x)]
    res2_fullname[res2_fullname == mo_fullname(x)] <- paste0(substr(mo_genus(res2_fullname_vector), 1, 1),
                                                             ". ",
                                                             suppressWarnings(mo_species(res2_fullname_vector)))
    if (sum(res1 == res2, na.rm = TRUE) > 0) {
      res1[res1 == res2] <- paste0(substr(mo_genus(res1[res1 == res2]), 1, 1),
                                   ". ",
                                   suppressWarnings(mo_species(res1[res1 == res2])))
    }
    res1[res1 != res2] <- res2_fullname
    result <- as.character(res1)
  } else {
    # return G. species
    result <- paste0(substr(mo_genus(x), 1, 1), ". ", suppressWarnings(mo_species(x)))
  }
  result[result %in% c(". ", "(. ")] <- ""
  mo_translate(result, language = language)
}

#' @rdname mo_property
#' @export
mo_subspecies <- function(x, Becker = FALSE, Lancefield = FALSE, language = NULL) {
  mo_property(x, "subspecies", Becker = Becker, Lancefield = Lancefield, language = language)
}

#' @rdname mo_property
#' @export
mo_species <- function(x, Becker = FALSE, Lancefield = FALSE, language = NULL) {
  mo_property(x, "species", Becker = Becker, Lancefield = Lancefield, language = language)
}

#' @rdname mo_property
#' @export
mo_genus <- function(x, language = NULL) {
  mo_property(x, "genus", language = language)
}

#' @rdname mo_property
#' @export
mo_family <- function(x) {
  mo_property(x, "family")
}

#' @rdname mo_property
#' @export
mo_order <- function(x) {
  mo_property(x, "order")
}

#' @rdname mo_property
#' @export
mo_class <- function(x) {
  mo_property(x, "class")
}

#' @rdname mo_property
#' @export
mo_phylum <- function(x) {
  mo_property(x, "phylum")
}

#' @rdname mo_property
#' @export
mo_type <- function(x, language = NULL) {
  mo_property(x, "type", language = language)
}

#' @rdname mo_property
#' @export
mo_gramstain <- function(x, language = NULL) {
  mo_property(x, "gramstain", language = language)
}

#' @rdname mo_property
#' @export
mo_aerobic <- function(x) {
  mo_property(x, "aerobic")
}

#' @rdname mo_property
#' @export
mo_property <- function(x, property = 'fullname', Becker = FALSE, Lancefield = FALSE, language = NULL) {
  property <- tolower(property[1])
  if (!property %in% colnames(AMR::microorganisms)) {
    stop("invalid property: ", property, " - use a column name of the `microorganisms` data set")
  }
  result1 <- as.mo(x = x, Becker = Becker, Lancefield = Lancefield) # this will give a warning if x cannot be coerced
  result2 <- suppressWarnings(
    data.frame(mo = result1, stringsAsFactors = FALSE) %>%
      left_join(AMR::microorganisms, by = "mo") %>%
      pull(property)
  )
  if (property != "aerobic") {
    # will else not retain `logical` class
    result2[x %in% c("", NA) | result2 %in% c("", NA, "(no MO)")] <- ""
    result2 <- mo_translate(result2, language = language)
  }
  result2
}

#' @rdname mo_property
#' @export
mo_taxonomy <- function(x) {
  x <- as.mo(x)
  base::list(phylum = mo_phylum(x),
             class = mo_class(x),
             order = mo_order(x),
             family = mo_family(x),
             genus = mo_genus(x),
             species = mo_species(x),
             subspecies = mo_subspecies(x))
}

#' @importFrom dplyr %>% case_when
mo_translate <- function(x, language) {
  if (is.null(language)) {
    language <- Sys.locale()
  } else {
    language <- tolower(language[1])
  }
  if (language %in% c("en", "")) {
    return(x)
  }

  supported <- c("en", "de", "nl", "es", "pt", "it", "fr")
  if (!language %in% supported) {
    stop("Unsupported language: '", language, "' - use one of: ", paste0("'", sort(supported), "'", collapse = ", "), call. = FALSE)
  }

  case_when(
    # German
    language == "de" ~ x %>%
      gsub("Coagulase Negative Staphylococcus","Koagulase-negative Staphylococcus", ., fixed = TRUE) %>%
      gsub("Coagulase Positive Staphylococcus","Koagulase-positive Staphylococcus", ., fixed = TRUE) %>%
      gsub("Beta-haemolytic Streptococcus",    "Beta-h\u00e4molytischer Streptococcus", ., fixed = TRUE) %>%
      gsub("(no MO)",          "(kein MO)", ., fixed = TRUE) %>%
      gsub("Negative rods",    "Negative St\u00e4bchen", ., fixed = TRUE) %>%
      gsub("Negative cocci",   "Negative Kokken", ., fixed = TRUE) %>%
      gsub("Positive rods",    "Positive St\u00e4bchen", ., fixed = TRUE) %>%
      gsub("Positive cocci",   "Positive Kokken", ., fixed = TRUE) %>%
      gsub("Parasites",        "Parasiten", ., fixed = TRUE) %>%
      gsub("Fungi and yeasts", "Pilze und Hefen", ., fixed = TRUE) %>%
      gsub("Bacteria",         "Bakterium", ., fixed = TRUE) %>%
      gsub("Fungus/yeast",     "Pilz/Hefe", ., fixed = TRUE) %>%
      gsub("Parasite",         "Parasit", ., fixed = TRUE) %>%
      gsub("biogroup",         "Biogruppe", ., fixed = TRUE) %>%
      gsub("biotype",          "Biotyp", ., fixed = TRUE) %>%
      gsub("vegetative",       "vegetativ", ., fixed = TRUE) %>%
      gsub("([([ ]*?)group",   "\\1Gruppe", .) %>%
      gsub("([([ ]*?)Group",   "\\1Gruppe", .),

    # Dutch
    language == "nl" ~ x %>%
      gsub("Coagulase Negative Staphylococcus","Coagulase-negatieve Staphylococcus", ., fixed = TRUE) %>%
      gsub("Coagulase Positive Staphylococcus","Coagulase-positieve Staphylococcus", ., fixed = TRUE) %>%
      gsub("Beta-haemolytic Streptococcus",    "Beta-hemolytische Streptococcus", ., fixed = TRUE) %>%
      gsub("(no MO)",          "(geen MO)", ., fixed = TRUE) %>%
      gsub("Negative rods",    "Negatieve staven", ., fixed = TRUE) %>%
      gsub("Negative cocci",   "Negatieve kokken", ., fixed = TRUE) %>%
      gsub("Positive rods",    "Positieve staven", ., fixed = TRUE) %>%
      gsub("Positive cocci",   "Positieve kokken", ., fixed = TRUE) %>%
      gsub("Parasites",        "Parasieten", ., fixed = TRUE) %>%
      gsub("Fungi and yeasts", "Schimmels en gisten", ., fixed = TRUE) %>%
      gsub("Bacteria",         "Bacterie", ., fixed = TRUE) %>%
      gsub("Fungus/yeast",     "Schimmel/gist", ., fixed = TRUE) %>%
      gsub("Parasite",         "Parasiet", ., fixed = TRUE) %>%
      gsub("biogroup",         "biogroep", ., fixed = TRUE) %>%
      # gsub("biotype",          "biotype", ., fixed = TRUE) %>%
      gsub("vegetative",       "vegetatief", ., fixed = TRUE) %>%
      gsub("([([ ]*?)group",   "\\1groep", .) %>%
      gsub("([([ ]*?)Group",   "\\1Groep", .),

    # Spanish
    language == "es" ~ x %>%
      gsub("Coagulase Negative Staphylococcus","Staphylococcus coagulasa negativo", ., fixed = TRUE) %>%
      gsub("Coagulase Positive Staphylococcus","Staphylococcus coagulasa positivo", ., fixed = TRUE) %>%
      gsub("Beta-haemolytic Streptococcus",    "Streptococcus Beta-hemol\u00edtico", ., fixed = TRUE) %>%
      gsub("(no MO)",          "(sin MO)", ., fixed = TRUE) %>%
      gsub("Negative rods",    "Bacilos negativos", ., fixed = TRUE) %>%
      gsub("Negative cocci",   "Cocos negativos", ., fixed = TRUE) %>%
      gsub("Positive rods",    "Bacilos positivos", ., fixed = TRUE) %>%
      gsub("Positive cocci",   "Cocos positivos", ., fixed = TRUE) %>%
      gsub("Parasites",        "Par\u00e1sitos", ., fixed = TRUE) %>%
      gsub("Fungi and yeasts", "Hongos y levaduras", ., fixed = TRUE) %>%
      # gsub("Bacteria",         "Bacteria", ., fixed = TRUE) %>%
      gsub("Fungus/yeast",     "Hongo/levadura", ., fixed = TRUE) %>%
      gsub("Parasite",         "Par\u00e1sito", ., fixed = TRUE) %>%
      gsub("biogroup",         "biogrupo", ., fixed = TRUE) %>%
      gsub("biotype",          "biotipo", ., fixed = TRUE) %>%
      gsub("vegetative",       "vegetativo", ., fixed = TRUE) %>%
      gsub("([([ ]*?)group",   "\\1grupo", .) %>%
      gsub("([([ ]*?)Group",   "\\1Grupo", .),

    # Portuguese
    language == "pt" ~ x %>%
      gsub("Coagulase Negative Staphylococcus","Staphylococcus coagulase negativo", ., fixed = TRUE) %>%
      gsub("Coagulase Positive Staphylococcus","Staphylococcus coagulase positivo", ., fixed = TRUE) %>%
      gsub("Beta-haemolytic Streptococcus",    "Streptococcus Beta-hemol\u00edtico", ., fixed = TRUE) %>%
      gsub("(no MO)",          "(sem MO)", ., fixed = TRUE) %>%
      gsub("Negative rods",    "Bacilos negativos", ., fixed = TRUE) %>%
      gsub("Negative cocci",   "Cocos negativos", ., fixed = TRUE) %>%
      gsub("Positive rods",    "Bacilos positivos", ., fixed = TRUE) %>%
      gsub("Positive cocci",   "Cocos positivos", ., fixed = TRUE) %>%
      gsub("Parasites",        "Parasitas", ., fixed = TRUE) %>%
      gsub("Fungi and yeasts", "Cogumelos e leveduras", ., fixed = TRUE) %>%
      gsub("Bacteria",         "Bact\u00e9ria", ., fixed = TRUE) %>%
      gsub("Fungus/yeast",     "Cogumelo/levedura", ., fixed = TRUE) %>%
      gsub("Parasite",         "Parasita", ., fixed = TRUE) %>%
      gsub("biogroup",         "biogrupo", ., fixed = TRUE) %>%
      gsub("biotype",          "bi\u00f3tipo", ., fixed = TRUE) %>%
      gsub("vegetative",       "vegetativo", ., fixed = TRUE) %>%
      gsub("([([ ]*?)group",   "\\1grupo", .) %>%
      gsub("([([ ]*?)Group",   "\\1Grupo", .),

    # Italian
    language == "it" ~ x %>%
      gsub("Coagulase Negative Staphylococcus","Staphylococcus negativo coagulasi", ., fixed = TRUE) %>%
      gsub("Coagulase Positive Staphylococcus","Staphylococcus positivo coagulasi", ., fixed = TRUE) %>%
      gsub("Beta-haemolytic Streptococcus",    "Streptococcus Beta-emolitico", ., fixed = TRUE) %>%
      gsub("(no MO)",          "(non MO)", ., fixed = TRUE) %>%
      gsub("Negative rods",    "Bastoncini Gram-negativi", ., fixed = TRUE) %>%
      gsub("Negative cocci",   "Cocchi Gram-negativi", ., fixed = TRUE) %>%
      gsub("Positive rods",    "Bastoncini Gram-positivi", ., fixed = TRUE) %>%
      gsub("Positive cocci",   "Cocchi Gram-positivi", ., fixed = TRUE) %>%
      gsub("Parasites",        "Parassiti", ., fixed = TRUE) %>%
      gsub("Fungi and yeasts", "Funghi e lieviti", ., fixed = TRUE) %>%
      gsub("Bacteria",         "Batterio", ., fixed = TRUE) %>%
      gsub("Fungus/yeast",     "Fungo/lievito", ., fixed = TRUE) %>%
      gsub("Parasite",         "Parassita", ., fixed = TRUE) %>%
      gsub("biogroup",         "biogruppo", ., fixed = TRUE) %>%
      gsub("biotype",          "biotipo", ., fixed = TRUE) %>%
      gsub("vegetative",       "vegetativo", ., fixed = TRUE) %>%
      gsub("([([ ]*?)group",   "\\1gruppo", .) %>%
      gsub("([([ ]*?)Group",   "\\1Gruppo", .),

    # French
    language == "fr" ~ x %>%
      gsub("Coagulase Negative Staphylococcus","Staphylococcus \u00e0 coagulase n\u00e9gative", ., fixed = TRUE) %>%
      gsub("Coagulase Positive Staphylococcus","Staphylococcus \u00e0 coagulase positif", ., fixed = TRUE) %>%
      gsub("Beta-haemolytic Streptococcus",    "Streptococcus B\u00eata-h\u00e9molytique", ., fixed = TRUE) %>%
      gsub("(no MO)",          "(pas MO)", ., fixed = TRUE) %>%
      gsub("Negative rods",    "Bacilles n\u00e9gatif", ., fixed = TRUE) %>%
      gsub("Negative cocci",   "Cocci n\u00e9gatif", ., fixed = TRUE) %>%
      gsub("Positive rods",    "Bacilles positif", ., fixed = TRUE) %>%
      gsub("Positive cocci",   "Cocci positif", ., fixed = TRUE) %>%
      # gsub("Parasites",        "Parasites", ., fixed = TRUE) %>%
      gsub("Fungi and yeasts", "Champignons et levures", ., fixed = TRUE) %>%
      gsub("Bacteria",         "Bact\u00e9rie", ., fixed = TRUE) %>%
      gsub("Fungus/yeast",     "Champignon/levure", ., fixed = TRUE) %>%
      # gsub("Parasite",         "Parasite", ., fixed = TRUE) %>%
      gsub("biogroup",         "biogroupe", ., fixed = TRUE) %>%
      # gsub("biotype",          "biotype", ., fixed = TRUE) %>%
      gsub("vegetative",       "v\u00e9g\u00e9tatif", ., fixed = TRUE) %>%
      gsub("([([ ]*?)group",   "\\1groupe", .) %>%
      gsub("([([ ]*?)Group",   "\\1Groupe", .)

  )

}
