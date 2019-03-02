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

#' Property of a microorganism
#'
#' Use these functions to return a specific property of a microorganism from the \code{\link{microorganisms}} data set. All input values will be evaluated internally with \code{\link{as.mo}}.
#' @param x any (vector of) text that can be coerced to a valid microorganism code with \code{\link{as.mo}}
#' @param property one of the column names of one of the \code{\link{microorganisms}} data set or \code{"shortname"}
#' @param language language of the returned text, defaults to system language (see \code{\link{get_locale}}) and can also be set with \code{\link{getOption}("AMR_locale")}. Use \code{language = NULL} or \code{language = ""} to prevent translation.
#' @param ... other parameters passed on to \code{\link{as.mo}}
#' @param open browse the URL using \code{\link[utils]{browseURL}}
#' @details All functions will return the most recently known taxonomic property according to the Catalogue of Life, except for \code{mo_ref}, \code{mo_authors} and \code{mo_year}. This leads to the following results:
#' \itemize{
#'   \item{\code{mo_fullname("Chlamydia psittaci")} will return \code{"Chlamydophila psittaci"} (with a warning about the renaming)}
#'   \item{\code{mo_ref("Chlamydia psittaci")} will return \code{"Page, 1968"} (with a warning about the renaming)}
#'   \item{\code{mo_ref("Chlamydophila psittaci")} will return \code{"Everett et al., 1999"} (without a warning)}
#' }
#'
#' The Gram stain - \code{mo_gramstain()} - will be determined on the taxonomic kingdom and phylum. According to Cavalier-Smith (2002) who defined subkingdoms Negibacteria and Posibacteria, only these phyla are Posibacteria: Actinobacteria, Chloroflexi, Firmicutes and Tenericutes (ref: \url{https://itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=956097}). These bacteria are considered Gram positive - all other bacteria are considered Gram negative. Species outside the kingdom of Bacteria will return a value \code{NA}.
#'
#' The function \code{mo_url()} will return the direct URL to the species in the Catalogue of Life.
#' @inheritSection get_locale Supported languages
#' @inheritSection catalogue_of_life Catalogue of Life
#' @inheritSection as.mo Source
#' @rdname mo_property
#' @name mo_property
#' @return \itemize{
#'   \item{An \code{integer} in case of \code{mo_year}}
#'   \item{A \code{list} in case of \code{mo_taxonomy}}
#'   \item{A named \code{character} in case of \code{mo_url}}
#'   \item{A \code{character} in all other cases}
#' }
#' @export
#' @seealso \code{\link{microorganisms}}
#' @inheritSection AMR Read more on our website!
#' @examples
#' ## taxonomic tree
#' mo_kingdom("E. coli")         # "Bacteria"
#' mo_phylum("E. coli")          # "Proteobacteria"
#' mo_class("E. coli")           # "Gammaproteobacteria"
#' mo_order("E. coli")           # "Enterobacteriales"
#' mo_family("E. coli")          # "Enterobacteriaceae"
#' mo_genus("E. coli")           # "Escherichia"
#' mo_species("E. coli")         # "coli"
#' mo_subspecies("E. coli")      # ""
#'
#' ## colloquial properties
#' mo_fullname("E. coli")        # "Escherichia coli"
#' mo_shortname("E. coli")       # "E. coli"
#'
#' ## other properties
#' mo_gramstain("E. coli")       # "Gram negative"
#' mo_type("E. coli")            # "Bacteria" (equal to kingdom)
#' mo_rank("E. coli")            # "species"
#' mo_url("E. coli")             # get the direct url to the Catalogue of Life
#'
#' ## scientific reference
#' mo_ref("E. coli")             # "Castellani et al., 1919"
#' mo_authors("E. coli")         # "Castellani et al."
#' mo_year("E. coli")            # 1919
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
#' mo_shortname("S. pyo", Lancefield = TRUE) # "GAS" ('Group A streptococci')
#'
#'
#' # language support for German, Dutch, Spanish, Portuguese, Italian and French
#' mo_gramstain("E. coli", language = "de")  # "Gramnegativ"
#' mo_gramstain("E. coli", language = "nl")  # "Gram-negatief"
#' mo_gramstain("E. coli", language = "es")  # "Gram negativo"
#'
#' # mo_type is equal to mo_kingdom, but mo_kingdom will remain official
#' mo_kingdom("E. coli")                     # "Bacteria" on a German system
#' mo_type("E. coli")                        # "Bakterien" on a German system
#' mo_type("E. coli")                        # "Bacteria" on an English system
#'
#' mo_fullname("S. pyogenes",
#'             Lancefield = TRUE,
#'             language = "de")              # "Streptococcus Gruppe A"
#' mo_fullname("S. pyogenes",
#'             Lancefield = TRUE,
#'             language = "nl")              # "Streptococcus groep A"
#'
#'
#' # get a list with the complete taxonomy (kingdom to subspecies)
#' mo_taxonomy("E. coli")
mo_fullname <- function(x, language = get_locale(), ...) {
  x <- mo_validate(x = x, property = "fullname", ...)
  mo_translate(x, language = language)
}

#' @rdname mo_property
#' @importFrom dplyr %>% left_join mutate pull
#' @export
mo_shortname <- function(x, language = get_locale(), ...) {
  dots <- list(...)
  Becker <- dots$Becker
  if (is.null(Becker)) {
    Becker <- FALSE
  }
  Lancefield <- dots$Lancefield
  if (is.null(Lancefield)) {
    Lancefield <- FALSE
  }

  shorten <- function(x) {
    # easiest: no transformations needed
    x <- mo_fullname(x, language = "en")
    # shorten for the ones that have a space: shorten first word and write out second word
    shorten_these <- x %like% " " & !x %like% "Streptococcus group "
    x[shorten_these] <- paste0(substr(x[shorten_these], 1, 1),
                              ". ",
                              x[shorten_these] %>%
                                strsplit(" ", fixed = TRUE) %>%
                                unlist() %>%
                                .[2])
    x
  }

  if (isFALSE(Becker) & isFALSE(Lancefield)) {
    result <- shorten(x)

  } else {
    # get result without transformations
    res1 <- AMR::as.mo(x, Becker = FALSE, Lancefield = FALSE, reference_df = dots$reference_df)
    # and result with transformations
    res2 <- suppressWarnings(AMR::as.mo(res1, ...))
    if (res1 == res2
        & !res1 %like% "^B_STRPT_GR") {
      result <- shorten(x)
    } else {
      res2_fullname <- mo_fullname(res2, language = language)
      res2_fullname[res2_fullname %like% " \\(CoNS\\)"] <- "CoNS"
      res2_fullname[res2_fullname %like% " \\(CoPS\\)"] <- "CoPS"
      res2_fullname[res2_fullname %like% " \\(KNS\\)"] <- "KNS"
      res2_fullname[res2_fullname %like% " \\(KPS\\)"] <- "KPS"
      res2_fullname[res2_fullname %like% " \\(CNS\\)"] <- "CNS"
      res2_fullname[res2_fullname %like% " \\(CPS\\)"] <- "CPS"
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
    }
  }

  mo_translate(result, language = language)
}

#' @rdname mo_property
#' @export
mo_subspecies <- function(x, language = get_locale(), ...) {
  mo_translate(mo_validate(x = x, property = "subspecies", ...), language = language)
}

#' @rdname mo_property
#' @export
mo_species <- function(x, language = get_locale(), ...) {
  mo_translate(mo_validate(x = x, property = "species", ...), language = language)
}

#' @rdname mo_property
#' @export
mo_genus <- function(x, language = get_locale(), ...) {
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
mo_kingdom <- function(x, ...) {
  mo_validate(x = x, property = "kingdom", ...)
}

#' @rdname mo_property
#' @export
mo_type <- function(x, language = get_locale(), ...) {
  mo_translate(mo_validate(x = x, property = "kingdom", ...), language = language)
}

#' @rdname mo_property
#' @export
mo_gramstain <- function(x, language = get_locale(), ...) {
  x.bak <- x
  x.mo <- as.mo(x, ...)
  x.phylum <- mo_phylum(x.mo)
  x[x.phylum %in% c("Actinobacteria",
                    "Chloroflexi",
                    "Firmicutes",
                    "Tenericutes")] <- "Gram positive"
  x[x != "Gram positive"] <- "Gram negative"
  x[mo_kingdom(x.mo) != "Bacteria"] <- NA_character_
  x[x.mo == "B_GRAMP"] <- "Gram positive"
  x[x.mo == "B_GRAMN"] <- "Gram negative"

  mo_translate(x, language = language)
}

#' @rdname mo_property
#' @export
mo_ref <- function(x, ...) {
  mo_validate(x = x, property = "ref", ...)
}

#' @rdname mo_property
#' @export
mo_authors <- function(x, ...) {
  x <- mo_validate(x = x, property = "ref", ...)
  # remove last 4 digits and presumably the comma and space that preceeds them
  x[!is.na(x)] <- gsub(",? ?[0-9]{4}", "", x[!is.na(x)])
  suppressWarnings(x)
}

#' @rdname mo_property
#' @export
mo_year <- function(x, ...) {
  x <- mo_validate(x = x, property = "ref", ...)
  # get last 4 digits
  x[!is.na(x)] <- gsub(".*([0-9]{4})$", "\\1", x[!is.na(x)])
  suppressWarnings(as.integer(x))
}

#' @rdname mo_property
#' @export
mo_rank <- function(x, ...) {
  mo_validate(x = x, property = "rank", ...)
}

#' @rdname mo_property
#' @export
mo_taxonomy <- function(x, ...) {
  x <- AMR::as.mo(x, ...)
  base::list(kingdom = mo_kingdom(x),
             phylum = mo_phylum(x),
             class = mo_class(x),
             order = mo_order(x),
             family = mo_family(x),
             genus = mo_genus(x),
             species = mo_species(x),
             subspecies = mo_subspecies(x))
}

#' @rdname mo_property
#' @importFrom utils browseURL
#' @export
mo_url <- function(x, open = FALSE, ...) {
  u <- mo_validate(x = x, property = "species_id", ...)
  u[u != ""] <- paste0(catalogue_of_life$url, "/details/species/id/", u)
  names(u) <- mo_fullname(x = x, ... =  ...)
  if (open == TRUE) {
    if (length(u) > 1) {
      warning("only the first URL will be opened, as `browseURL` only suports one string.")
    }
    browseURL(u[1L])
  }
  u
}


#' @rdname mo_property
#' @importFrom data.table data.table as.data.table setkey
#' @export
mo_property <- function(x, property = 'fullname', language = get_locale(), ...) {
  if (length(property) != 1L) {
    stop("'property' must be of length 1.")
  }
  if (!property %in% colnames(AMR::microorganisms)) {
    stop("invalid property: '", property, "' - use a column name of the `microorganisms` data set")
  }

  mo_translate(mo_validate(x = x, property = property, ...), language = language)
}

#' @importFrom dplyr %>% case_when
mo_translate <- function(x, language) {
  if (is.null(language)) {
    return(x)
  }
  if (language %in% c("en", "")) {
    return(x)
  }

  supported <- c("en", "de", "nl", "es", "pt", "it", "fr")
  if (!language %in% supported) {
    stop("Unsupported language: '", language, "' - use one of: ", paste0("'", sort(supported), "'", collapse = ", "), call. = FALSE)
  }

  x_tobetranslated <- grepl(x = x,
                            pattern = "(Coagulase Negative Staphylococcus|Coagulase Positive Staphylococcus|Beta-haemolytic Streptococcus|unknown Gram negatives|unknown Gram positives|unknown name|unknown kingdom|unknown phylum|unknown class|unknown order|unknown family|unknown genus|unknown species|unknown subspecies|unknown rank|CoNS|CoPS|Gram negative|Gram positive|Bacteria|Fungi|Protozoa|biogroup|biotype|vegetative|group|Group)")

  if (sum(x_tobetranslated, na.rm = TRUE) == 0) {
    return(x)
  }

  # only translate the ones that need translation
  x[x_tobetranslated] <- case_when(
    # German
    language == "de" ~ x[x_tobetranslated] %>%
      gsub("Coagulase Negative Staphylococcus","Koagulase-negative Staphylococcus", ., fixed = TRUE) %>%
      gsub("Coagulase Positive Staphylococcus","Koagulase-positive Staphylococcus", ., fixed = TRUE) %>%
      gsub("Beta-haemolytic Streptococcus",    "Beta-h\u00e4molytischer Streptococcus", ., fixed = TRUE) %>%
      gsub("unknown Gram negatives",           "unbekannte Gramnegativen", ., fixed = TRUE) %>%
      gsub("unknown Gram positives",           "unbekannte Grampositiven", ., fixed = TRUE) %>%
      gsub("unknown name",                     "unbekannte Name", ., fixed = TRUE) %>%
      gsub("unknown kingdom",                  "unbekanntes Reich", ., fixed = TRUE) %>%
      gsub("unknown phylum",                   "unbekannter Stamm", ., fixed = TRUE) %>%
      gsub("unknown class",                    "unbekannte Klasse", ., fixed = TRUE) %>%
      gsub("unknown order",                    "unbekannte Ordnung", ., fixed = TRUE) %>%
      gsub("unknown family",                   "unbekannte Familie", ., fixed = TRUE) %>%
      gsub("unknown genus",                    "unbekannte Gattung", ., fixed = TRUE) %>%
      gsub("unknown species",                  "unbekannte Art", ., fixed = TRUE) %>%
      gsub("unknown subspecies",               "unbekannte Unterart", ., fixed = TRUE) %>%
      gsub("unknown rank",                     "unbekannter Rang", ., fixed = TRUE) %>%
      gsub("(CoNS)",           "(KNS)", ., fixed = TRUE) %>%
      gsub("(CoPS)",           "(KPS)", ., fixed = TRUE) %>%
      gsub("Gram negative",    "Gramnegativ", ., fixed = TRUE) %>%
      gsub("Gram positive",    "Grampositiv", ., fixed = TRUE) %>%
      gsub("Bacteria",         "Bakterien", ., fixed = TRUE) %>%
      gsub("Fungi",            "Hefen/Pilze", ., fixed = TRUE) %>%
      gsub("Protozoa",         "Protozoen", ., fixed = TRUE) %>%
      gsub("biogroup",         "Biogruppe", ., fixed = TRUE) %>%
      gsub("biotype",          "Biotyp", ., fixed = TRUE) %>%
      gsub("vegetative",       "vegetativ", ., fixed = TRUE) %>%
      gsub("([([ ]*?)group",   "\\1Gruppe", .) %>%
      gsub("([([ ]*?)Group",   "\\1Gruppe", .) %>%
      iconv(to = "UTF-8"),

    # Dutch
    language == "nl" ~ x[x_tobetranslated] %>%
      gsub("Coagulase Negative Staphylococcus","Coagulase-negatieve Staphylococcus", ., fixed = TRUE) %>%
      gsub("Coagulase Positive Staphylococcus","Coagulase-positieve Staphylococcus", ., fixed = TRUE) %>%
      gsub("Beta-haemolytic Streptococcus",    "Beta-hemolytische Streptococcus", ., fixed = TRUE) %>%
      gsub("unknown Gram negatives",           "onbekende Gram-negatieven", ., fixed = TRUE) %>%
      gsub("unknown Gram positives",           "onbekende Gram-positieven", ., fixed = TRUE) %>%
      gsub("unknown name",                     "onbekende naam", ., fixed = TRUE) %>%
      gsub("unknown kingdom",                  "onbekend koninkrijk", ., fixed = TRUE) %>%
      gsub("unknown phylum",                   "onbekende fylum", ., fixed = TRUE) %>%
      gsub("unknown class",                    "onbekende klasse", ., fixed = TRUE) %>%
      gsub("unknown order",                    "onbekende orde", ., fixed = TRUE) %>%
      gsub("unknown family",                   "onbekende familie", ., fixed = TRUE) %>%
      gsub("unknown genus",                    "onbekend geslacht", ., fixed = TRUE) %>%
      gsub("unknown species",                  "onbekende soort", ., fixed = TRUE) %>%
      gsub("unknown subspecies",               "onbekende ondersoort", ., fixed = TRUE) %>%
      gsub("unknown rank",                     "onbekende rang", ., fixed = TRUE) %>%
      gsub("(CoNS)",           "(CNS)", ., fixed = TRUE) %>%
      gsub("(CoPS)",           "(CPS)", ., fixed = TRUE) %>%
      gsub("Gram negative",    "Gram-negatief", ., fixed = TRUE) %>%
      gsub("Gram positive",    "Gram-positief", ., fixed = TRUE) %>%
      gsub("Bacteria",         "Bacteri\u00ebn", ., fixed = TRUE) %>%
      gsub("Fungi",            "Schimmels/gisten", ., fixed = TRUE) %>%
      gsub("Protozoa",         "protozo\u00ebn", ., fixed = TRUE) %>%
      gsub("biogroup",         "biogroep", ., fixed = TRUE) %>%
      # gsub("biotype",          "biotype", ., fixed = TRUE) %>%
      gsub("vegetative",       "vegetatief", ., fixed = TRUE) %>%
      gsub("([([ ]*?)group",   "\\1groep", .) %>%
      gsub("([([ ]*?)Group",   "\\1Groep", .) %>%
      iconv(to = "UTF-8"),

    # Spanish
    language == "es" ~ x[x_tobetranslated] %>%
      gsub("Coagulase Negative Staphylococcus","Staphylococcus coagulasa negativo", ., fixed = TRUE) %>%
      gsub("Coagulase Positive Staphylococcus","Staphylococcus coagulasa positivo", ., fixed = TRUE) %>%
      gsub("Beta-haemolytic Streptococcus",    "Streptococcus Beta-hemol\u00edtico", ., fixed = TRUE) %>%
      gsub("unknown Gram negatives",           "Gram negativos desconocidos", ., fixed = TRUE) %>%
      gsub("unknown Gram positives",           "Gram positivos desconocidos", ., fixed = TRUE) %>%
      gsub("unknown name",                     "nombre desconocido", ., fixed = TRUE) %>%
      gsub("unknown kingdom",                  "reino desconocido", ., fixed = TRUE) %>%
      gsub("unknown phylum",                   "filo desconocido", ., fixed = TRUE) %>%
      gsub("unknown class",                    "clase desconocida", ., fixed = TRUE) %>%
      gsub("unknown order",                    "orden desconocido", ., fixed = TRUE) %>%
      gsub("unknown family",                   "familia desconocida", ., fixed = TRUE) %>%
      gsub("unknown genus",                    "g\u00e9nero desconocido", ., fixed = TRUE) %>%
      gsub("unknown species",                  "especie desconocida", ., fixed = TRUE) %>%
      gsub("unknown subspecies",               "subespecie desconocida", ., fixed = TRUE) %>%
      gsub("unknown rank",                     "rango desconocido", ., fixed = TRUE) %>%
      gsub("Gram negative",    "Gram negativo", ., fixed = TRUE) %>%
      gsub("Gram positive",    "Gram positivo", ., fixed = TRUE) %>%
      gsub("Bacteria",         "Bacterias", ., fixed = TRUE) %>%
      gsub("Fungi",            "Hongos", ., fixed = TRUE) %>%
      gsub("Protozoa",         "Protozoarios", ., fixed = TRUE) %>%
      gsub("biogroup",         "biogrupo", ., fixed = TRUE) %>%
      gsub("biotype",          "biotipo", ., fixed = TRUE) %>%
      gsub("vegetative",       "vegetativo", ., fixed = TRUE) %>%
      gsub("([([ ]*?)group",   "\\1grupo", .) %>%
      gsub("([([ ]*?)Group",   "\\1Grupo", .) %>%
      iconv(to = "UTF-8"),

    # Italian
    language == "it" ~ x[x_tobetranslated] %>%
      gsub("Coagulase Negative Staphylococcus","Staphylococcus negativo coagulasi", ., fixed = TRUE) %>%
      gsub("Coagulase Positive Staphylococcus","Staphylococcus positivo coagulasi", ., fixed = TRUE) %>%
      gsub("Beta-haemolytic Streptococcus",    "Streptococcus Beta-emolitico", ., fixed = TRUE) %>%
      gsub("unknown Gram negatives",           "Gram negativi sconosciuti", ., fixed = TRUE) %>%
      gsub("unknown Gram positives",           "Gram positivi sconosciuti", ., fixed = TRUE) %>%
      gsub("unknown name",                     "nome sconosciuto", ., fixed = TRUE) %>%
      gsub("unknown kingdom",                  "regno sconosciuto", ., fixed = TRUE) %>%
      gsub("unknown phylum",                   "phylum sconosciuto", ., fixed = TRUE) %>%
      gsub("unknown class",                    "classe sconosciuta", ., fixed = TRUE) %>%
      gsub("unknown order",                    "ordine sconosciuto", ., fixed = TRUE) %>%
      gsub("unknown family",                   "famiglia sconosciuta", ., fixed = TRUE) %>%
      gsub("unknown genus",                    "genere sconosciuto", ., fixed = TRUE) %>%
      gsub("unknown species",                  "specie sconosciute", ., fixed = TRUE) %>%
      gsub("unknown subspecies",               "sottospecie sconosciute", ., fixed = TRUE) %>%
      gsub("unknown rank",                     "grado sconosciuto", ., fixed = TRUE) %>%
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
    language == "fr" ~ x[x_tobetranslated] %>%
      gsub("Coagulase Negative Staphylococcus","Staphylococcus \u00e0 coagulase n\u00e9gative", ., fixed = TRUE) %>%
      gsub("Coagulase Positive Staphylococcus","Staphylococcus \u00e0 coagulase positif", ., fixed = TRUE) %>%
      gsub("Beta-haemolytic Streptococcus",    "Streptococcus B\u00eata-h\u00e9molytique", ., fixed = TRUE) %>%
      gsub("unknown Gram negatives",           "Gram n\u00e9gatifs inconnus", ., fixed = TRUE) %>%
      gsub("unknown Gram positives",           "Gram positifs inconnus", ., fixed = TRUE) %>%
      gsub("unknown name",                     "nom inconnu", ., fixed = TRUE) %>%
      gsub("unknown kingdom",                  "r\u00e8gme inconnu", ., fixed = TRUE) %>%
      gsub("unknown phylum",                   "embranchement inconnu", ., fixed = TRUE) %>%
      gsub("unknown class",                    "classe inconnue", ., fixed = TRUE) %>%
      gsub("unknown order",                    "ordre inconnu", ., fixed = TRUE) %>%
      gsub("unknown family",                   "famille inconnue", ., fixed = TRUE) %>%
      gsub("unknown genus",                    "genre inconnu", ., fixed = TRUE) %>%
      gsub("unknown species",                  "esp\u00e8ce inconnue", ., fixed = TRUE) %>%
      gsub("unknown subspecies",               "sous-esp\u00e8ce inconnue", ., fixed = TRUE) %>%
      gsub("unknown rank",                     "rang inconnu", ., fixed = TRUE) %>%
      gsub("Gram negative",    "Gram n\u00e9gatif", ., fixed = TRUE) %>%
      gsub("Gram positive",    "Gram positif", ., fixed = TRUE) %>%
      gsub("Bacteria",         "Bact\u00e9ries", ., fixed = TRUE) %>%
      gsub("Fungi",            "Champignons", ., fixed = TRUE) %>%
      gsub("Protozoa",         "Protozoaires", ., fixed = TRUE) %>%
      gsub("biogroup",         "biogroupe", ., fixed = TRUE) %>%
      # gsub("biotype",          "biotype", ., fixed = TRUE) %>%
      gsub("vegetative",       "v\u00e9g\u00e9tatif", ., fixed = TRUE) %>%
      gsub("([([ ]*?)group",   "\\1groupe", .) %>%
      gsub("([([ ]*?)Group",   "\\1Groupe", .) %>%
      iconv(to = "UTF-8"),

    # Portuguese
    language == "pt" ~ x[x_tobetranslated] %>%
      gsub("Coagulase Negative Staphylococcus","Staphylococcus coagulase negativo", ., fixed = TRUE) %>%
      gsub("Coagulase Positive Staphylococcus","Staphylococcus coagulase positivo", ., fixed = TRUE) %>%
      gsub("Beta-haemolytic Streptococcus",    "Streptococcus Beta-hemol\u00edtico", ., fixed = TRUE) %>%
      gsub("unknown Gram negatives",           "Gram negativos desconhecidos", ., fixed = TRUE) %>%
      gsub("unknown Gram positives",           "Gram positivos desconhecidos", ., fixed = TRUE) %>%
      gsub("unknown name",                     "nome desconhecido", ., fixed = TRUE) %>%
      gsub("unknown kingdom",                  "reino desconhecido", ., fixed = TRUE) %>%
      gsub("unknown phylum",                   "filo desconhecido", ., fixed = TRUE) %>%
      gsub("unknown class",                    "classe desconhecida", ., fixed = TRUE) %>%
      gsub("unknown order",                    "ordem desconhecido", ., fixed = TRUE) %>%
      gsub("unknown family",                   "fam\u00edlia desconhecida", ., fixed = TRUE) %>%
      gsub("unknown genus",                    "g\u00eanero desconhecido", ., fixed = TRUE) %>%
      gsub("unknown species",                  "esp\u00e9cies desconhecida", ., fixed = TRUE) %>%
      gsub("unknown subspecies",               "subesp\u00e9cies desconhecida", ., fixed = TRUE) %>%
      gsub("unknown rank",                     "classifica\u00e7\u00e3o desconhecido", ., fixed = TRUE) %>%
      gsub("Gram negative",    "Gram negativo", ., fixed = TRUE) %>%
      gsub("Gram positive",    "Gram positivo", ., fixed = TRUE) %>%
      gsub("Bacteria",         "Bact\u00e9rias", ., fixed = TRUE) %>%
      gsub("Fungi",            "Fungos", ., fixed = TRUE) %>%
      gsub("Protozoa",         "Protozo\u00e1rios", ., fixed = TRUE) %>%
      gsub("biogroup",         "biogrupo", ., fixed = TRUE) %>%
      gsub("biotype",          "bi\u00f3tipo", ., fixed = TRUE) %>%
      gsub("vegetative",       "vegetativo", ., fixed = TRUE) %>%
      gsub("([([ ]*?)group",   "\\1grupo", .) %>%
      gsub("([([ ]*?)Group",   "\\1Grupo", .) %>%
      iconv(to = "UTF-8"))

  x

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

  if (!"AMR" %in% base::.packages()) {
    library("AMR")
    # check onLoad() in R/zzz.R: data tables are created there.
  }

  if (!all(x %in% microorganisms[, property])
      | Becker %in% c(TRUE, "all")
      | Lancefield %in% c(TRUE, "all")) {
    exec_as.mo(x, property = property, ...)
  } else {
    if (property == "mo") {
      return(structure(x, class = "mo"))
    } else {
      return(x)
    }
  }
}
