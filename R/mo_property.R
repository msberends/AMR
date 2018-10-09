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
#' @param property one of the column names of one of the \code{\link{microorganisms}} data set or \code{"shortname"}
#' @param language language of the returned text, defaults to English (\code{"en"}) and can also be set with \code{\link{getOption}("AMR_locale")}. Either one of \code{"en"} (English), \code{"de"} (German), \code{"nl"} (Dutch), \code{"es"} (Spanish) or \code{"pt"} (Portuguese).
#' @param ... other parameters passed on to \code{\link{as.mo}}
#' @inheritSection as.mo ITIS
#' @inheritSection as.mo Source
#' @rdname mo_property
#' @name mo_property
#' @return A \code{list} (in case of \code{mo_taxonomy}) or a \code{character} otherwise
#' @export
#' @seealso \code{\link{microorganisms}}
#' @examples
#' # All properties of Escherichia coli
#' mo_subkingdom("E. coli")      # "Negibacteria"
#' mo_phylum("E. coli")          # "Proteobacteria"
#' mo_class("E. coli")           # "Gammaproteobacteria"
#' mo_order("E. coli")           # "Enterobacteriales"
#' mo_family("E. coli")          # "Enterobacteriaceae"
#' mo_genus("E. coli")           # "Escherichia"
#' mo_species("E. coli")         # "coli"
#' mo_subspecies("E. coli")      # NA
#' mo_fullname("E. coli")        # "Escherichia coli"
#' mo_shortname("E. coli")       # "E. coli"
#' mo_gramstain("E. coli")       # "Gram negative"
#' mo_TSN("E. coli")             # 285
#' mo_type("E. coli")            # "Bacteria"
#' mo_ref("E. coli")             # "Castellani and Chalmers, 1919"
#'
#'
#' # Abbreviations known in the field
#' mo_genus("MRSA")              # "Staphylococcus"
#' mo_species("MRSA")            # "aureus"
#' mo_shortname("MRSA")          # "S. aureus"
#' mo_gramstain("MRSA")          # "Gram positive"
#'
#' mo_genus("VISA")              # "Staphylococcus"
#' mo_species("VISA")            # "aureus"
#'
#'
#' # Known subspecies
#' mo_genus("doylei")            # "Campylobacter"
#' mo_species("doylei")          # "jejuni"
#' mo_fullname("doylei")         # "Campylobacter jejuni doylei"
#'
#' mo_fullname("K. pneu rh")     # "Klebsiella pneumoniae rhinoscleromatis"
#' mo_shortname("K. pneu rh")    # "K. pneumoniae"
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
#' mo_gramstain("E. coli", language = "de")  # "Gramnegativ"
#' mo_gramstain("E. coli", language = "nl")  # "Gram-negatief"
#' mo_gramstain("E. coli", language = "es")  # "Gram negativo"
#'
#' mo_fullname("S. pyogenes",
#'             Lancefield = TRUE,
#'             language = "de")              # "Streptococcus Gruppe A"
#' mo_fullname("S. pyogenes",
#'             Lancefield = TRUE,
#'             language = "nl")              # "Streptococcus groep A"
#'
#'
#' # Complete taxonomy up to Subkingdom, returns a list
#' mo_taxonomy("E. coli")
mo_fullname <- function(x, language = NULL, ...) {
  x <- mo_validate(x = x, property = "fullname", ...)
  mo_translate(x, language = language)
}

#' @rdname mo_property
#' @importFrom dplyr %>% left_join mutate pull
#' @export
mo_shortname <- function(x, language = NULL, ...) {
  dots <- list(...)
  Becker <- dots$Becker
  if (is.null(Becker)) {
    Becker <- FALSE
  }
  Lancefield <- dots$Lancefield
  if (is.null(Lancefield)) {
    Lancefield <- FALSE
  }
  if (Becker %in% c(TRUE, "all") | Lancefield == TRUE) {
    res1 <- AMR::as.mo(x, Becker = FALSE, Lancefield = FALSE, reference_df = dots$reference_df)
    res2 <- suppressWarnings(AMR::as.mo(res1, ...))
    res2_fullname <- mo_fullname(res2)
    res2_fullname[res2_fullname %like% "\\(CoNS\\)"] <- "CoNS"
    res2_fullname[res2_fullname %like% "\\(CoPS\\)"] <- "CoPS"
    res2_fullname <- gsub("Streptococcus (group|Gruppe|gruppe|groep|grupo|gruppo|groupe) (.)",
                          "G\\2S",
                          res2_fullname) # turn "Streptococcus group A" and "Streptococcus grupo A" to "GAS"
    res2_fullname_vector <- res2_fullname[res2_fullname == mo_fullname(res1)]
    res2_fullname[res2_fullname == mo_fullname(res1)] <- paste0(substr(mo_genus(res2_fullname_vector), 1, 1),
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
    x <- AMR::as.mo(x, ...)
    suppressWarnings(
      result <- data.frame(mo = x) %>%
        left_join(AMR::microorganisms, by = "mo") %>%
        mutate(shortname = ifelse(!is.na(genus) & !is.na(species), paste0(substr(genus, 1, 1), ". ", species), NA_character_)) %>%
        pull(shortname)
    )
  }
  mo_translate(result, language = language)
}

#' @rdname mo_property
#' @export
mo_subspecies <- function(x, language = NULL, ...) {
  mo_translate(mo_validate(x = x, property = "subspecies", ...), language = language)
}

#' @rdname mo_property
#' @export
mo_species <- function(x, language = NULL, ...) {
  mo_translate(mo_validate(x = x, property = "species", ...), language = language)
}

#' @rdname mo_property
#' @export
mo_genus <- function(x, language = NULL, ...) {
  mo_translate(mo_validate(x = x, property = "genus", ...), language = language)
}

#' @rdname mo_property
#' @export
mo_family <- function(x, ...) {
  mo_validate(x = x, property = "family", ...)
}

#' @rdname mo_property
#' @export
mo_order <- function(x, ...) {
  mo_validate(x = x, property = "order", ...)
}

#' @rdname mo_property
#' @export
mo_class <- function(x, ...) {
  mo_validate(x = x, property = "class", ...)
}

#' @rdname mo_property
#' @export
mo_phylum <- function(x, ...) {
  mo_validate(x = x, property = "phylum", ...)
}

#' @rdname mo_property
#' @export
mo_subkingdom <- function(x, ...) {
  mo_validate(x = x, property = "subkingdom", ...)
}

#' @rdname mo_property
#' @export
mo_ref <- function(x, ...) {
  mo_validate(x = x, property = "ref", ...)
}

#' @rdname mo_property
#' @export
mo_type <- function(x, language = NULL, ...) {
  mo_translate(mo_validate(x = x, property = "type", ...), language = language)
}

#' @rdname mo_property
#' @export
mo_TSN  <- function(x, ...) {
  mo_validate(x = x, property = "tsn", ...)
}

#' @rdname mo_property
#' @export
mo_gramstain <- function(x, language = NULL, ...) {
  mo_translate(mo_validate(x = x, property = "gramstain", ...), language = language)
}

#' @rdname mo_property
#' @importFrom data.table data.table as.data.table setkey
#' @export
mo_property <- function(x, property = 'fullname', language = NULL, ...) {
  if (length(property) != 1L) {
    stop("'property' must be of length 1.")
  }
  if (!property %in% colnames(AMR::microorganisms)) {
    stop("invalid property: '", property, "' - use a column name of the `microorganisms` data set")
  }

  mo_translate(mo_validate(x = x, property = property, ...), language = language)
}

#' @rdname mo_property
#' @export
mo_taxonomy <- function(x, ...) {
  x <- AMR::as.mo(x, ...)
  base::list(subkingdom = mo_subkingdom(x),
             phylum = mo_phylum(x),
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
    language <- getOption("AMR_locale", default = "en")[1L]
  } else {
    language <- tolower(language[1L])
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
      gsub("Gram negative",    "Gramnegativ", ., fixed = TRUE) %>%
      gsub("Gram positive",    "Grampositiv", ., fixed = TRUE) %>%
      gsub("Bacteria",         "Bakterien", ., fixed = TRUE) %>%
      gsub("Fungi",            "Hefen/Pilze", ., fixed = TRUE) %>%
      gsub("Protozoa",         "Protozoen", ., fixed = TRUE) %>%
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
      gsub("Gram negative",    "Gram-negatief", ., fixed = TRUE) %>%
      gsub("Gram positive",    "Gram-positief", ., fixed = TRUE) %>%
      gsub("Bacteria",         "Bacteri\u00ebn", ., fixed = TRUE) %>%
      gsub("Fungi",            "Schimmels/gisten", ., fixed = TRUE) %>%
      gsub("Protozoa",         "protozo\u00ebn", ., fixed = TRUE) %>%
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
      gsub("Gram negative",    "Gram negativo", ., fixed = TRUE) %>%
      gsub("Gram positive",    "Gram positivo", ., fixed = TRUE) %>%
      gsub("Bacteria",         "Bacterias", ., fixed = TRUE) %>%
      gsub("Fungi",            "Hongos", ., fixed = TRUE) %>%
      gsub("Protozoa",         "Protozoarios", ., fixed = TRUE) %>%
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
      gsub("Gram negative",    "Gram negativo", ., fixed = TRUE) %>%
      gsub("Gram positive",    "Gram positivo", ., fixed = TRUE) %>%
      gsub("Bacteria",         "Bact\u00e9rias", ., fixed = TRUE) %>%
      gsub("Fungi",            "Fungos", ., fixed = TRUE) %>%
      gsub("Protozoa",         "Protozo\u00e1rios", ., fixed = TRUE) %>%
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
      gsub("Gram negative",    "Gram negativo", ., fixed = TRUE) %>%
      gsub("Gram positive",    "Gram positivo", ., fixed = TRUE) %>%
      gsub("Bacteria",         "Batteri", ., fixed = TRUE) %>%
      gsub("Fungi",            "Fungo", ., fixed = TRUE) %>%
      gsub("Protozoa",         "Protozoi", ., fixed = TRUE) %>%
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
      gsub("Gram negative",    "Gram n\u00e9gatif", ., fixed = TRUE) %>%
      gsub("Gram positive",    "Gram positif", ., fixed = TRUE) %>%
      gsub("Bacteria",         "Bact\u00e9ries", ., fixed = TRUE) %>%
      gsub("Fungi",            "Champignons", ., fixed = TRUE) %>%
      gsub("Protozoa",         "Protozoaires", ., fixed = TRUE) %>%
      gsub("biogroup",         "biogroupe", ., fixed = TRUE) %>%
      # gsub("biotype",          "biotype", ., fixed = TRUE) %>%
      gsub("vegetative",       "v\u00e9g\u00e9tatif", ., fixed = TRUE) %>%
      gsub("([([ ]*?)group",   "\\1groupe", .) %>%
      gsub("([([ ]*?)Group",   "\\1Groupe", .)

  )

}

mo_validate <- function(x, property, ...) {

  dots <- list(...)
  Becker <- dots$Becker
  if (is.null(Becker)) {
    Becker <- FALSE
  }
  Lancefield <- dots$Lancefield
  if (is.null(Lancefield)) {
    Lancefield <- FALSE
  }

  if (!all(x %in% AMR::microorganisms[, property]) | Becker %in% c(TRUE, "all") | Lancefield == TRUE) {
    exec_as.mo(x, property = property, ...)
  } else {
    x
  }
}
