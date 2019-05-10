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

#' Property of a microorganism
#'
#' Use these functions to return a specific property of a microorganism from the \code{\link{microorganisms}} data set. All input values will be evaluated internally with \code{\link{as.mo}}.
#' @param x any (vector of) text that can be coerced to a valid microorganism code with \code{\link{as.mo}}
#' @param property one of the column names of one of the \code{\link{microorganisms}} data set or \code{"shortname"}
#' @param language language of the returned text, defaults to system language (see \code{\link{get_locale}}) and can also be set with \code{\link{getOption}("AMR_locale")}. Use \code{language = NULL} or \code{language = ""} to prevent translation.
#' @param ... other parameters passed on to \code{\link{as.mo}}
#' @param open browse the URL using \code{\link[utils]{browseURL}()}
#' @details All functions will return the most recently known taxonomic property according to the Catalogue of Life, except for \code{mo_ref}, \code{mo_authors} and \code{mo_year}. This leads to the following results:
#' \itemize{
#'   \item{\code{mo_fullname("Chlamydia psittaci")} will return \code{"Chlamydophila psittaci"} (with a warning about the renaming)}
#'   \item{\code{mo_ref("Chlamydia psittaci")} will return \code{"Page, 1968"} (with a warning about the renaming)}
#'   \item{\code{mo_ref("Chlamydophila psittaci")} will return \code{"Everett et al., 1999"} (without a warning)}
#' }
#'
#' The Gram stain - \code{mo_gramstain()} - will be determined on the taxonomic kingdom and phylum. According to Cavalier-Smith (2002) who defined subkingdoms Negibacteria and Posibacteria, only these phyla are Posibacteria: Actinobacteria, Chloroflexi, Firmicutes and Tenericutes. These bacteria are considered Gram positive - all other bacteria are considered Gram negative. Species outside the kingdom of Bacteria will return a value \code{NA}.
#'
#' All output will be \link{translate}d where possible.
#'
#' The function \code{mo_url()} will return the direct URL to the online database entry, which also shows the scientific reference of the concerned species.
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
#' mo_url("E. coli")             # get the direct url to the online database entry
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
#' mo_fullname("S. epi", Becker = TRUE)      # "Coagulase-negative Staphylococcus (CoNS)"
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
#' # get a list with the complete taxonomy (from kingdom to subspecies)
#' mo_taxonomy("E. coli")
mo_fullname <- function(x, language = get_locale(), ...) {
  x <- mo_validate(x = x, property = "fullname", ...)
  t(x, language = language)
}

#' @rdname mo_property
#' @importFrom dplyr %>% mutate pull
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

  t(result, language = language)
}

#' @rdname mo_property
#' @export
mo_subspecies <- function(x, language = get_locale(), ...) {
  t(mo_validate(x = x, property = "subspecies", ...), language = language)
}

#' @rdname mo_property
#' @export
mo_species <- function(x, language = get_locale(), ...) {
  t(mo_validate(x = x, property = "species", ...), language = language)
}

#' @rdname mo_property
#' @export
mo_genus <- function(x, language = get_locale(), ...) {
  t(mo_validate(x = x, property = "genus", ...), language = language)
}

#' @rdname mo_property
#' @export
mo_family <- function(x, language = get_locale(), ...) {
  t(mo_validate(x = x, property = "family", ...), language = language)
}

#' @rdname mo_property
#' @export
mo_order <- function(x, language = get_locale(), ...) {
  t(mo_validate(x = x, property = "order", ...), language = language)
}

#' @rdname mo_property
#' @export
mo_class <- function(x, language = get_locale(), ...) {
  t(mo_validate(x = x, property = "class", ...), language = language)
}

#' @rdname mo_property
#' @export
mo_phylum <- function(x, language = get_locale(), ...) {
  t(mo_validate(x = x, property = "phylum", ...), language = language)
}

#' @rdname mo_property
#' @export
mo_kingdom <- function(x, language = get_locale(), ...) {
  if (all(x %in% AMR::microorganisms$kingdom)) {
    return(x)
  }
  x <- as.mo(x, ...)
  kngdm <- mo_validate(x = x, property = "kingdom", ...)
  if (language != "en") {
    kngdm[x == "UNKNOWN"] <- t(kngdm[x == "UNKNOWN"], language = language)
  }
  kngdm
}

#' @rdname mo_property
#' @export
mo_type <- function(x, language = get_locale(), ...) {
  t(mo_validate(x = x, property = "kingdom", ...), language = language)
}

#' @rdname mo_property
#' @export
mo_gramstain <- function(x, language = get_locale(), ...) {
  x.mo <- as.mo(x, ...)
  x.phylum <- mo_phylum(x.mo, language = "en")
  x[x.phylum %in% c("Actinobacteria",
                    "Chloroflexi",
                    "Firmicutes",
                    "Tenericutes")] <- "Gram positive"
  x[x != "Gram positive"] <- "Gram negative"
  x[mo_kingdom(x.mo, language = "en") != "Bacteria"] <- NA_character_
  x[x.mo == "B_GRAMP"] <- "Gram positive"
  x[x.mo == "B_GRAMN"] <- "Gram negative"

  t(x, language = language)
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
mo_taxonomy <- function(x, language = get_locale(),  ...) {
  x <- AMR::as.mo(x, ...)
  base::list(kingdom = mo_kingdom(x, language = language),
             phylum = mo_phylum(x, language = language),
             class = mo_class(x, language = language),
             order = mo_order(x, language = language),
             family = mo_family(x, language = language),
             genus = mo_genus(x, language = language),
             species = mo_species(x, language = language),
             subspecies = mo_subspecies(x, language = language))
}

#' @rdname mo_property
#' @importFrom utils browseURL
#' @importFrom dplyr %>% left_join select mutate case_when
#' @export
mo_url <- function(x, open = FALSE, ...) {
  mo <- AMR::as.mo(x = x, ... = ...)
  df <- data.frame(mo, stringsAsFactors = FALSE) %>%
    left_join(select(AMR::microorganisms, mo, source, species_id), by = "mo") %>%
    mutate(url = case_when(source == "CoL" ~
                             paste0(gsub("{year}", catalogue_of_life$year, catalogue_of_life$url_CoL, fixed = TRUE), "details/species/id/", species_id),
                           source == "DSMZ" ~
                             paste0(catalogue_of_life$url_DSMZ, "?bnu_no=", species_id, "#", species_id),
                           TRUE ~
                             NA_character_))

  u <- df$url
  names(u) <- mo_fullname(mo)
  if (open == TRUE) {
    if (length(u) > 1) {
      warning("only the first URL will be opened, as `browseURL()` only suports one string.")
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

  t(mo_validate(x = x, property = property, ...), language = language)
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

  # try to catch an error when inputting an invalid parameter
  # so the 'call.' can be set to FALSE
  tryCatch(x[1L] %in% AMR::microorganisms[1, property],
           error = function(e) stop(e$message, call. = FALSE))

  if (!all(x %in% AMR::microorganisms[, property])
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
