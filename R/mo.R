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

#' Transform to microorganism ID
#'
#' Use this function to determine a valid microorganism ID (\code{mo}). Determination is done using intelligent rules and the complete taxonomic kingdoms Bacteria, Chromista, Protozoa, Archaea and most microbial species from the kingdom Fungi (see Source). The input can be almost anything: a full name (like \code{"Staphylococcus aureus"}), an abbreviated name (like \code{"S. aureus"}), an abbreviation known in the field (like \code{"MRSA"}), or just a genus. Please see Examples.
#' @param x a character vector or a \code{data.frame} with one or two columns
#' @param Becker a logical to indicate whether \emph{Staphylococci} should be categorised into coagulase-negative \emph{Staphylococci} ("CoNS") and coagulase-positive \emph{Staphylococci} ("CoPS") instead of their own species, according to Karsten Becker \emph{et al.} [1,2]. Note that this does not include species that were newly named after these publications, like \emph{S. caeli}.
#'
#'   This excludes \emph{Staphylococcus aureus} at default, use \code{Becker = "all"} to also categorise \emph{S. aureus} as "CoPS".
#' @param Lancefield a logical to indicate whether beta-haemolytic \emph{Streptococci} should be categorised into Lancefield groups instead of their own species, according to Rebecca C. Lancefield [3]. These \emph{Streptococci} will be categorised in their first group, e.g. \emph{Streptococcus dysgalactiae} will be group C, although officially it was also categorised into groups G and L.
#'
#'   This excludes \emph{Enterococci} at default (who are in group D), use \code{Lancefield = "all"} to also categorise all \emph{Enterococci} as group D.
#' @param allow_uncertain a logical (\code{TRUE} or \code{FALSE}) or a value between 0 and 3 to indicate whether the input should be checked for less possible results, see Details
#' @param reference_df a \code{data.frame} to use for extra reference when translating \code{x} to a valid \code{mo}. See \code{\link{set_mo_source}} and \code{\link{get_mo_source}} to automate the usage of your own codes (e.g. used in your analysis or organisation).
#' @param ... other parameters passed on to functions
#' @rdname as.mo
#' @aliases mo
#' @keywords mo Becker becker Lancefield lancefield guess
#' @details
#' \strong{General info} \cr
#' A microbial ID from this package (class: \code{mo}) typically looks like these examples:\cr
#' \preformatted{
#'   Code              Full name
#'   ---------------   --------------------------------------
#'   B_KLBSL           Klebsiella
#'   B_KLBSL_PNE       Klebsiella pneumoniae
#'   B_KLBSL_PNE_RHI   Klebsiella pneumoniae rhinoscleromatis
#'   |   |    |   |
#'   |   |    |   |
#'   |   |    |    ----> subspecies, a 3-4 letter acronym
#'   |   |     ----> species, a 3-4 letter acronym
#'   |    ----> genus, a 5-7 letter acronym, mostly without vowels
#'    ----> taxonomic kingdom: A (Archaea), AN (Animalia), B (Bacteria),
#'                             C (Chromista), F (Fungi), P (Protozoa) or
#'                             PL (Plantae)
#' }
#'
#' Values that cannot be coered will be considered 'unknown' and will get the MO code \code{UNKNOWN}.
#'
#' Use the \code{\link{mo_property}_*} functions to get properties based on the returned code, see Examples.
#'
#' The algorithm uses data from the Catalogue of Life (see below) and from one other source (see \code{?microorganisms}).
#'
# \strong{Self-learning algoritm} \cr
# The \code{as.mo()} function gains experience from previously determined microbial IDs and learns from it. This drastically improves both speed and reliability. Use \code{clean_mo_history()} to reset the algorithms. Only experience from your current \code{AMR} package version is used. This is done because in the future the taxonomic tree (which is included in this package) may change for any organism and it consequently has to rebuild its knowledge.
#
# Usually, any guess after the first try runs 80-95\% faster than the first try.
#
# For now, learning only works per session. If R is closed or terminated, the algorithms reset. This will probably be resolved in a next version.
#
#' \strong{Intelligent rules} \cr
#' This function uses intelligent rules to help getting fast and logical results. It tries to find matches in this order:
#' \itemize{
#'   \item{Valid MO codes and full names: it first searches in already valid MO code and known genus/species combinations}
#'   \item{Human pathogenic prevalence: it first searches in more prevalent microorganisms, then less prevalent ones (see \emph{Microbial prevalence of pathogens in humans} below)}
#'   \item{Taxonomic kingdom: it first searches in Bacteria/Chromista, then Fungi, then Protozoa}
#'   \item{Breakdown of input values: from here it starts to breakdown input values to find possible matches}
#' }
#'
#' A couple of effects because of these rules:
#' \itemize{
#'   \item{\code{"E. coli"} will return the ID of \emph{Escherichia coli} and not \emph{Entamoeba coli}, although the latter would alphabetically come first}
#'   \item{\code{"H. influenzae"} will return the ID of \emph{Haemophilus influenzae} and not \emph{Haematobacter influenzae} for the same reason}
#'   \item{Something like \code{"stau"} or \code{"S aur"} will return the ID of \emph{Staphylococcus aureus} and not \emph{Staphylococcus auricularis}}
#' }
#' This means that looking up human pathogenic microorganisms takes less time than looking up human non-pathogenic microorganisms.
#'
#' \strong{Uncertain results} \cr
#' The algorithm can additionally use three different levels of uncertainty to guess valid results. The default is \code{allow_uncertain = TRUE}, which is equal to uncertainty level 2. Using \code{allow_uncertain = FALSE} will skip all of these additional rules:
#' \itemize{
#'   \item{(uncertainty level 1): It tries to look for only matching genera, previously accepted (but now invalid) taxonomic names and misspelled input}
#'   \item{(uncertainty level 2): It removed parts between brackets, strips off words from the end one by one and re-evaluates the input with all previous rules}
#'   \item{(uncertainty level 3): It strips off words from the start one by one and tries any part of the name}
#' }
#'
#' You can also use e.g. \code{as.mo(..., allow_uncertain = 1)} to only allow up to level 1 uncertainty.
#'
#' Examples:
#' \itemize{
#'   \item{\code{"Streptococcus group B (known as S. agalactiae)"}. The text between brackets will be removed and a warning will be thrown that the result \emph{Streptococcus group B} (\code{B_STRPT_GRB}) needs review.}
#'   \item{\code{"S. aureus - please mind: MRSA"}. The last word will be stripped, after which the function will try to find a match. If it does not, the second last word will be stripped, etc. Again, a warning will be thrown that the result \emph{Staphylococcus aureus} (\code{B_STPHY_AUR}) needs review.}
#'   \item{\code{"Fluoroquinolone-resistant Neisseria gonorrhoeae"}. The first word will be stripped, after which the function will try to find a match. A warning will be thrown that the result \emph{Neisseria gonorrhoeae} (\code{B_NESSR_GON}) needs review.}
#' }
#'
#' Use \code{mo_failures()} to get a vector with all values that could not be coerced to a valid value.
#'
#' Use \code{mo_uncertainties()} to get a data.frame with all values that were coerced to a valid value, but with uncertainty.
#'
#' Use \code{mo_renamed()} to get a vector with all values that could be coerced based on an old, previously accepted taxonomic name.
#'
#' \strong{Microbial prevalence of pathogens in humans} \cr
#' The intelligent rules take into account microbial prevalence of pathogens in humans. It uses three groups and all (sub)species are in only one group. These groups are:
#' \itemize{
#'   \item{1 (most prevalent): class is Gammaproteobacteria \strong{or} genus is one of: \emph{Enterococcus}, \emph{Staphylococcus}, \emph{Streptococcus}.}
#'   \item{2: phylum is one of: Proteobacteria, Firmicutes, Actinobacteria, Sarcomastigophora \strong{or} genus is one of: \emph{Aspergillus}, \emph{Bacteroides}, \emph{Candida}, \emph{Capnocytophaga}, \emph{Chryseobacterium}, \emph{Cryptococcus}, \emph{Elisabethkingia}, \emph{Flavobacterium}, \emph{Fusobacterium}, \emph{Giardia}, \emph{Leptotrichia}, \emph{Mycoplasma}, \emph{Prevotella}, \emph{Rhodotorula}, \emph{Treponema}, \emph{Trichophyton}, \emph{Ureaplasma}.}
#'   \item{3 (least prevalent): all others.}
#' }
#'
#' Group 1 contains all common Gram positives and Gram negatives, like all Enterobacteriaceae and e.g. \emph{Pseudomonas} and \emph{Legionella}.
#'
#' Group 2 probably contains less microbial pathogens; all other members of phyla that were found in humans in the Northern Netherlands between 2001 and 2018.
#' @inheritSection catalogue_of_life Catalogue of Life
#  (source as a section here, so it can be inherited by other man pages:)
#' @section Source:
#' [1] Becker K \emph{et al.} \strong{Coagulase-Negative Staphylococci}. 2014. Clin Microbiol Rev. 27(4): 870–926. \url{https://dx.doi.org/10.1128/CMR.00109-13}
#'
#' [2] Becker K \emph{et al.} \strong{Implications of identifying the recently defined members of the \emph{S. aureus} complex, \emph{S. argenteus} and \emph{S. schweitzeri}: A position paper of members of the ESCMID Study Group for staphylococci and Staphylococcal Diseases (ESGS).} 2019. Clin Microbiol Infect. \url{https://doi.org/10.1016/j.cmi.2019.02.028}
#'
#' [3] Lancefield RC \strong{A serological differentiation of human and other groups of hemolytic streptococci}. 1933. J Exp Med. 57(4): 571–95. \url{https://dx.doi.org/10.1084/jem.57.4.571}
#'
#' [4] Catalogue of Life: Annual Checklist (public online taxonomic database), \url{http://www.catalogueoflife.org} (check included annual version with \code{\link{catalogue_of_life_version}()}).
#' @export
#' @return Character (vector) with class \code{"mo"}
#' @seealso \code{\link{microorganisms}} for the \code{data.frame} that is being used to determine ID's. \cr
#' The \code{\link{mo_property}} functions (like \code{\link{mo_genus}}, \code{\link{mo_gramstain}}) to get properties based on the returned code.
#' @inheritSection AMR Read more on our website!
#' @importFrom dplyr %>% pull left_join
#' @examples
#' \donttest{
#' # These examples all return "B_STPHY_AUR", the ID of S. aureus:
#' as.mo("sau") # WHONET code
#' as.mo("stau")
#' as.mo("STAU")
#' as.mo("staaur")
#' as.mo("S. aureus")
#' as.mo("S aureus")
#' as.mo("Staphylococcus aureus")
#' as.mo("Staphylococcus aureus (MRSA)")
#' as.mo("Sthafilokkockus aaureuz") # handles incorrect spelling
#' as.mo("MRSA")   # Methicillin Resistant S. aureus
#' as.mo("VISA")   # Vancomycin Intermediate S. aureus
#' as.mo("VRSA")   # Vancomycin Resistant S. aureus
#' as.mo(22242419) # Catalogue of Life ID
#'
#' # Dyslexia is no problem - these all work:
#' as.mo("Ureaplasma urealyticum")
#' as.mo("Ureaplasma urealyticus")
#' as.mo("Ureaplasmium urealytica")
#' as.mo("Ureaplazma urealitycium")
#'
#' as.mo("Streptococcus group A")
#' as.mo("GAS") # Group A Streptococci
#' as.mo("GBS") # Group B Streptococci
#'
#' as.mo("S. epidermidis")                 # will remain species: B_STPHY_EPI
#' as.mo("S. epidermidis", Becker = TRUE)  # will not remain species: B_STPHY_CNS
#'
#' as.mo("S. pyogenes")                    # will remain species: B_STRPT_PYO
#' as.mo("S. pyogenes", Lancefield = TRUE) # will not remain species: B_STRPT_GRA
#'
#' # All mo_* functions use as.mo() internally too (see ?mo_property):
#' mo_genus("E. coli")           # returns "Escherichia"
#' mo_gramstain("E. coli")       # returns "Gram negative"#'
#' }
#' \dontrun{
#' df$mo <- as.mo(df$microorganism_name)
#'
#' # the select function of tidyverse is also supported:
#' library(dplyr)
#' df$mo <- df %>%
#'   select(microorganism_name) %>%
#'   as.mo()
#'
#' # and can even contain 2 columns, which is convenient for genus/species combinations:
#' df$mo <- df %>%
#'   select(genus, species) %>%
#'   as.mo()
#' # although this works easier and does the same:
#' df <- df %>%
#'   mutate(mo = as.mo(paste(genus, species)))
#' }
as.mo <- function(x, Becker = FALSE, Lancefield = FALSE, allow_uncertain = TRUE, reference_df = get_mo_source(), ...) {
  if (!"AMR" %in% base::.packages()) {
    require("AMR")
    # check onLoad() in R/zzz.R: data tables are created there.
  }

  # WHONET: xxx = no growth
  x[tolower(as.character(paste0(x, ""))) %in% c("", "xxx", "na", "nan")] <- NA_character_

  uncertainty_level <- translate_allow_uncertain(allow_uncertain)
  # mo_hist <- get_mo_history(x, uncertainty_level, force = isTRUE(list(...)$force_mo_history))

  if (mo_source_isvalid(reference_df)
      & isFALSE(Becker)
      & isFALSE(Lancefield)
      & !is.null(reference_df)
      & all(x %in% reference_df[,1][[1]])) {

    # has valid own reference_df
    # (data.table not faster here)
    reference_df <- reference_df %>% filter(!is.na(mo))
    # keep only first two columns, second must be mo
    if (colnames(reference_df)[1] == "mo") {
      reference_df <- reference_df[, c(2, 1)]
    } else {
      reference_df <- reference_df[, c(1, 2)]
    }
    colnames(reference_df)[1] <- "x"
    # remove factors, just keep characters
    suppressWarnings(
      reference_df[] <- lapply(reference_df, as.character)
    )
    suppressWarnings(
      y <- data.frame(x = x, stringsAsFactors = FALSE) %>%
        left_join(reference_df, by = "x") %>%
        pull("mo")
    )

  } else if (all(x %in% AMR::microorganisms$mo)
             & isFALSE(Becker)
             & isFALSE(Lancefield)) {
    y <- x

    # } else if (!any(is.na(mo_hist))
    #            & isFALSE(Becker)
    #            & isFALSE(Lancefield)) {
    #   # check previously found results
    #   y <- mo_hist

  } else if (all(tolower(x) %in% microorganismsDT$fullname_lower)
             & isFALSE(Becker)
             & isFALSE(Lancefield)) {
    # we need special treatment for very prevalent full names, they are likely! (case insensitive)
    # e.g. as.mo("Staphylococcus aureus")
    y <- microorganismsDT[prevalence == 1][data.table(fullname_lower = tolower(x)),
                                           on = "fullname_lower",
                                           "mo"][[1]]
    if (any(is.na(y))) {
      y[is.na(y)] <- microorganismsDT[prevalence == 2][data.table(fullname_lower = tolower(x[is.na(y)])),
                                                       on = "fullname_lower",
                                                       "mo"][[1]]
    }
    if (any(is.na(y))) {
      y[is.na(y)] <- microorganismsDT[prevalence == 3][data.table(fullname_lower = tolower(x[is.na(y)])),
                                                       on = "fullname_lower",
                                                       "mo"][[1]]
    }
    # save them to history
    set_mo_history(x, y, 0, force = isTRUE(list(...)$force_mo_history))

  } else {
    # will be checked for mo class in validation and uses exec_as.mo internally if necessary
    y <- mo_validate(x = x, property = "mo",
                     Becker = Becker, Lancefield = Lancefield,
                     allow_uncertain = uncertainty_level, reference_df = reference_df,
                     force_mo_history = isTRUE(list(...)$force_mo_history),
                     ...)
  }


  structure(.Data = y, class = "mo")
}

#' @rdname as.mo
#' @export
is.mo <- function(x) {
  identical(class(x), "mo")
}

#' @importFrom dplyr %>% pull left_join n_distinct progress_estimated filter distinct
#' @importFrom data.table data.table as.data.table setkey
#' @importFrom crayon magenta red blue silver italic
# param property a column name of AMR::microorganisms
# param initial_search logical - is FALSE when coming from uncertain tries, which uses exec_as.mo internally too
# param force_mo_history logical - whether found result must be saved with set_mo_history (default FALSE on non-interactive sessions)
# param debug logical - show different lookup texts while searching
exec_as.mo <- function(x,
                       Becker = FALSE,
                       Lancefield = FALSE,
                       allow_uncertain = TRUE,
                       reference_df = get_mo_source(),
                       property = "mo",
                       initial_search = TRUE,
                       force_mo_history = FALSE,
                       debug = FALSE) {

  if (!"AMR" %in% base::.packages()) {
    require("AMR")
    # check onLoad() in R/zzz.R: data tables are created there.
  }

  # WHONET: xxx = no growth
  x[tolower(as.character(paste0(x, ""))) %in% c("", "xxx", "na", "nan")] <- NA_character_

  if (initial_search == TRUE) {
    options(mo_failures = NULL)
    options(mo_uncertainties = NULL)
    options(mo_renamed = NULL)
  }

  if (NCOL(x) == 2) {
    # support tidyverse selection like: df %>% select(colA, colB)
    # paste these columns together
    x_vector <- vector("character", NROW(x))
    for (i in 1:NROW(x)) {
      x_vector[i] <- paste(pull(x[i,], 1), pull(x[i,], 2), sep = " ")
    }
    x <- x_vector
  } else {
    if (NCOL(x) > 2) {
      stop('`x` can be 2 columns at most', call. = FALSE)
    }
    x[is.null(x)] <- NA

    # support tidyverse selection like: df %>% select(colA)
    if (!is.vector(x) & !is.null(dim(x))) {
      x <- pull(x, 1)
    }
  }

  notes <- character(0)
  uncertainties <- data.frame(input = character(0),
                              fullname = character(0),
                              mo = character(0))
  failures <- character(0)
  uncertainty_level <- translate_allow_uncertain(allow_uncertain)

  x_input <- x
  # already strip leading and trailing spaces
  x <- trimws(x, which = "both")
  # only check the uniques, which is way faster
  x <- unique(x)
  # remove empty values (to later fill them in again with NAs)
  # ("xxx" is WHONET code for 'no growth')
  x <- x[!is.na(x)
         & !is.null(x)
         & !identical(x, "")
         & !identical(x, "xxx")]

  # conversion of old MO codes from v0.5.0 (ITIS) to later versions (Catalogue of Life)
  if (any(x %like% "^[BFP]_[A-Z]{3,7}") & !all(x %in% microorganisms$mo)) {
    leftpart <- gsub("^([BFP]_[A-Z]{3,7}).*", "\\1", x)
    if (any(leftpart %in% names(mo_codes_v0.5.0))) {
      rightpart <- gsub("^[BFP]_[A-Z]{3,7}(.*)", "\\1", x)
      leftpart <- mo_codes_v0.5.0[leftpart]
      x[!is.na(leftpart)] <- paste0(leftpart[!is.na(leftpart)], rightpart[!is.na(leftpart)])
    }
    # now check if some are still old
    still_old <- x[x %in% names(mo_codes_v0.5.0)]
    if (length(still_old) > 0) {
      x[x %in% names(mo_codes_v0.5.0)] <- data.frame(old = still_old, stringsAsFactors = FALSE) %>%
        left_join(data.frame(old = names(mo_codes_v0.5.0),
                             new = mo_codes_v0.5.0,
                             stringsAsFactors = FALSE), by = "old") %>%
        # if they couldn't be found, replace them with the old ones again,
        # so they will throw a warning in the end
        mutate(new = ifelse(is.na(new), old, new)) %>%
        pull(new)
    }
  }

  # defined df to check for
  if (!is.null(reference_df)) {
    if (!mo_source_isvalid(reference_df)) {
      stop("`reference_df` must contain a column `mo` with values from the 'microorganisms' data set.", call. = FALSE)
    }
    reference_df <- reference_df %>% filter(!is.na(mo))
    # keep only first two columns, second must be mo
    if (colnames(reference_df)[1] == "mo") {
      reference_df <- reference_df[, c(2, 1)]
    } else {
      reference_df <- reference_df[, c(1, 2)]
    }
    colnames(reference_df)[1] <- "x"
    # remove factors, just keep characters
    suppressWarnings(
      reference_df[] <- lapply(reference_df, as.character)
    )
  }

  # all empty
  if (all(identical(trimws(x_input), "") | is.na(x_input) | length(x) == 0)) {
    if (property == "mo") {
      return(structure(rep(NA_character_, length(x_input)),
                       class = "mo"))
    } else {
      return(rep(NA_character_, length(x_input)))
    }

  } else if (all(x %in% reference_df[, 1][[1]])) {
    # all in reference df
    colnames(reference_df)[1] <- "x"
    suppressWarnings(
      x <- data.frame(x = x, stringsAsFactors = FALSE) %>%
        left_join(reference_df, by = "x") %>%
        left_join(AMR::microorganisms, by = "mo") %>%
        pull(property)
    )

  } else if (all(x %in% AMR::microorganisms$mo)) {
    # existing mo codes when not looking for property "mo", like mo_genus("B_ESCHR_COL")
    y <- microorganismsDT[prevalence == 1][data.table(mo = x), on = "mo", ..property][[1]]
    if (any(is.na(y))) {
      y[is.na(y)] <- microorganismsDT[prevalence == 2][data.table(mo = x[is.na(y)]),
                                                       on = "mo",
                                                       ..property][[1]]
    }
    if (any(is.na(y))) {
      y[is.na(y)] <- microorganismsDT[prevalence == 3][data.table(mo = x[is.na(y)]),
                                                       on = "mo",
                                                       ..property][[1]]
    }
    x <- y

  } else if (all(x %in% read_mo_history(uncertainty_level,
                                        force = force_mo_history)$x)) {
    # previously found code
    x <- microorganismsDT[data.table(mo = get_mo_history(x,
                                                         uncertainty_level,
                                                         force = force_mo_history)),
                          on = "mo", ..property][[1]]

  } else if (all(tolower(x) %in% microorganismsDT$fullname_lower)) {
    # we need special treatment for very prevalent full names, they are likely!
    # e.g. as.mo("Staphylococcus aureus")
    y <- microorganismsDT[prevalence == 1][data.table(fullname_lower = tolower(x)), on = "fullname_lower", ..property][[1]]
    if (any(is.na(y))) {
      y[is.na(y)] <- microorganismsDT[prevalence == 2][data.table(fullname_lower = tolower(x[is.na(y)])),
                                                       on = "fullname_lower",
                                                       ..property][[1]]
    }
    if (any(is.na(y))) {
      y[is.na(y)] <- microorganismsDT[prevalence == 3][data.table(fullname_lower = tolower(x[is.na(y)])),
                                                       on = "fullname_lower",
                                                       ..property][[1]]
    }
    x <- y

  } else if (all(toupper(x) %in% AMR::microorganisms.codes$code)) {
    # commonly used MO codes
    y <- as.data.table(AMR::microorganisms.codes)[data.table(code = toupper(x)), on = "code", ]
    # save them to history
    set_mo_history(x, y$mo, 0, force = force_mo_history)

    x <- microorganismsDT[data.table(mo = y[["mo"]]), on = "mo", ..property][[1]]

  } else if (!all(x %in% AMR::microorganisms[, property])) {

    strip_whitespace <- function(x) {
      # all whitespaces (tab, new lines, etc.) should be one space
      # and spaces before and after should be omitted
      trimws(gsub("[\\s]+", " ", x, perl = TRUE), which = "both")
    }

    x <- strip_whitespace(x)
    x_backup <- x

    # remove spp and species
    x <- gsub(" +(spp.?|ssp.?|sp.? |ss ?.?|subsp.?|subspecies|biovar |serovar |species)", " ", x_backup, ignore.case = TRUE)
    x <- strip_whitespace(x)

    x_backup_without_spp <- x
    x_species <- paste(x, "species")
    # translate to English for supported languages of mo_property
    x <- gsub("(gruppe|groep|grupo|gruppo|groupe)", "group", x, ignore.case = TRUE)
    x <- gsub("(hefe|gist|gisten|levadura|lievito|fermento|levure)[a-z]*", "yeast", x, ignore.case = TRUE)
    x <- gsub("(schimmels?|mofo|molde|stampo|moisissure)[a-z]*", "fungus", x, ignore.case = TRUE)
    # remove non-text in case of "E. coli" except dots and spaces
    x <- gsub("[^.a-zA-Z0-9/ \\-]+", "", x)
    # replace minus by a space
    x <- gsub("-+", " ", x)
    # replace hemolytic by haemolytic
    x <- gsub("ha?emoly", "haemoly", x)
    # place minus back in streptococci
    x <- gsub("(alpha|beta|gamma).?ha?emoly", "\\1-haemoly", x)
    # remove genus as first word
    x <- gsub("^Genus ", "", x)
    # allow characters that resemble others ----
    if (initial_search == FALSE) {
      x <- tolower(x)
      x <- gsub("[iy]+", "[iy]+", x)
      x <- gsub("(c|k|q|qu|s|z|x|ks)+", "(c|k|q|qu|s|z|x|ks)+", x)
      x <- gsub("(ph|f|v)+", "(ph|f|v)+", x)
      x <- gsub("(th|t)+", "(th|t)+", x)
      x <- gsub("a+", "a+", x)
      x <- gsub("u+", "u+", x)
      # allow any ending of -um, -us, -ium, -icum, -ius, -icus, -ica and -a (needs perl for the negative backward lookup):
      x <- gsub("(u\\+\\(c\\|k\\|q\\|qu\\+\\|s\\|z\\|x\\|ks\\)\\+)(?![a-z])",
                "(u[s|m]|[iy][ck]?u[ms]|[iy]?[ck]?a)", x, ignore.case = TRUE, perl = TRUE)
      x <- gsub("(\\[iy\\]\\+\\(c\\|k\\|q\\|qu\\+\\|s\\|z\\|x\\|ks\\)\\+a\\+)(?![a-z])",
                "(u[s|m]|[iy][ck]?u[ms]|[iy]?[ck]?a)", x, ignore.case = TRUE, perl = TRUE)
      x <- gsub("(\\[iy\\]\\+u\\+m)(?![a-z])",
                "(u[s|m]|[iy][ck]?u[ms]|[iy]?[ck]?a)", x, ignore.case = TRUE, perl = TRUE)
      x <- gsub("e+", "e+", x, ignore.case = TRUE)
      x <- gsub("o+", "o+", x, ignore.case = TRUE)
      x <- gsub("(.)\\1+", "\\1+", x)
      # allow ending in -en or -us
      x <- gsub("e\\+n(?![a-z[])", "(e+n|u+(c|k|q|qu|s|z|x|ks)+)", x, ignore.case = TRUE, perl = TRUE)
    }
    x <- strip_whitespace(x)

    x_trimmed <- x
    x_trimmed_species <- paste(x_trimmed, "species")
    x_trimmed_without_group <- gsub(" gro.u.p$", "", x_trimmed, ignore.case = TRUE)
    # remove last part from "-" or "/"
    x_trimmed_without_group <- gsub("(.*)[-/].*", "\\1", x_trimmed_without_group)
    # replace space and dot by regex sign
    x_withspaces <- gsub("[ .]+", ".* ", x)
    x <- gsub("[ .]+", ".*", x)
    # add start en stop regex
    x <- paste0('^', x, '$')
    x_withspaces_start_only <- paste0('^', x_withspaces)
    x_withspaces_end_only <- paste0(x_withspaces, '$')
    x_withspaces_start_end <- paste0('^', x_withspaces, '$')

    if (isTRUE(debug)) {
      cat(paste0('x                       "', x, '"\n'))
      cat(paste0('x_species               "', x_species, '"\n'))
      cat(paste0('x_withspaces_start_only "', x_withspaces_start_only, '"\n'))
      cat(paste0('x_withspaces_end_only   "', x_withspaces_end_only, '"\n'))
      cat(paste0('x_withspaces_start_end  "', x_withspaces_start_end, '"\n'))
      cat(paste0('x_backup                "', x_backup, '"\n'))
      cat(paste0('x_backup_without_spp    "', x_backup_without_spp, '"\n'))
      cat(paste0('x_trimmed               "', x_trimmed, '"\n'))
      cat(paste0('x_trimmed_species       "', x_trimmed_species, '"\n'))
      cat(paste0('x_trimmed_without_group "', x_trimmed_without_group, '"\n'))
    }

    progress <- progress_estimated(n = length(x), min_time = 3)

    for (i in 1:length(x)) {

      progress$tick()$print()

      if (initial_search == TRUE) {
        found <- microorganismsDT[mo == get_mo_history(x_backup[i],
                                                       uncertainty_level,
                                                       force = force_mo_history),
                                  ..property][[1]]
        # previously found result
        if (length(found) > 0) {
          x[i] <- found[1L]
          next
        }
      }

      found <- microorganismsDT[mo == toupper(x_backup[i]), ..property][[1]]
      # is a valid MO code
      if (length(found) > 0) {
        x[i] <- found[1L]
        next
      }

      found <- microorganismsDT[fullname_lower %in% tolower(c(x_backup[i], x_backup_without_spp[i])), ..property][[1]]
      # most probable: is exact match in fullname
      if (length(found) > 0) {
        x[i] <- found[1L]
        if (initial_search == TRUE) {
          set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history)
        }
        next
      }

      found <- microorganismsDT[col_id == x_backup[i], ..property][[1]]
      # is a valid Catalogue of Life ID
      if (NROW(found) > 0) {
        x[i] <- found[1L]
        if (initial_search == TRUE) {
          set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history)
        }
        next
      }


      # WHONET: xxx = no growth
      if (tolower(as.character(paste0(x_backup_without_spp[i], ""))) %in% c("", "xxx", "na", "nan")) {
        x[i] <- NA_character_
        next
      }

      if (tolower(x_backup_without_spp[i]) %in% c("other", "none", "unknown")) {
        # empty and nonsense values, ignore without warning
        x[i] <- microorganismsDT[mo == "UNKNOWN", ..property][[1]]
        if (initial_search == TRUE) {
          set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history)
        }
        next
      }

      # check for very small input, but ignore the O antigens of E. coli
      if (nchar(gsub("[^a-zA-Z]", "", x_trimmed[i])) < 3
          & !x_backup_without_spp[i] %like% "O?(26|103|104|104|111|121|145|157)") {
        # check if search term was like "A. species", then return first genus found with ^A
        # if (x_backup[i] %like% "[a-z]+ species" | x_backup[i] %like% "[a-z] spp[.]?") {
        #   # get mo code of first hit
        #   found <- microorganismsDT[fullname %like% x_withspaces_start_only[i], mo]
        #   if (length(found) > 0) {
        #     mo_code <- found[1L] %>% strsplit("_") %>% unlist() %>% .[1:2] %>% paste(collapse = "_")
        #     found <- microorganismsDT[mo == mo_code, ..property][[1]]
        #     # return first genus that begins with x_trimmed, e.g. when "E. spp."
        #     if (length(found) > 0) {
        #       x[i] <- found[1L]
        #       if (initial_search == TRUE) {
        #         set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history)
        #       }
        #       next
        #     }
        #   }
        # }
        # fewer than 3 chars and not looked for species, add as failure
        x[i] <- microorganismsDT[mo == "UNKNOWN", ..property][[1]]
        if (initial_search == TRUE) {
          failures <- c(failures, x_backup[i])
          set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history)
        }
        next
      }

      if (x_backup_without_spp[i] %like% "virus") {
        # there is no fullname like virus, so don't try to coerce it
        x[i] <- microorganismsDT[mo == "UNKNOWN", ..property][[1]]
        if (initial_search == TRUE) {
          failures <- c(failures, x_backup[i])
          set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history)
        }
        next
      }

      # translate known trivial abbreviations to genus + species ----
      if (!is.na(x_trimmed[i])) {
        if (toupper(x_backup_without_spp[i]) %in% c('MRSA', 'MSSA', 'VISA', 'VRSA')) {
          x[i] <- microorganismsDT[mo == 'B_STPHY_AUR', ..property][[1]][1L]
          if (initial_search == TRUE) {
            set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history)
          }
          next
        }
        if (toupper(x_backup_without_spp[i]) %in% c('MRSE', 'MSSE')) {
          x[i] <- microorganismsDT[mo == 'B_STPHY_EPI', ..property][[1]][1L]
          if (initial_search == TRUE) {
            set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history)
          }
          next
        }
        if (toupper(x_backup_without_spp[i]) == "VRE"
            | x_backup_without_spp[i] %like% '(enterococci|enterokok|enterococo)[a-z]*?$')  {
          x[i] <- microorganismsDT[mo == 'B_ENTRC', ..property][[1]][1L]
          if (initial_search == TRUE) {
            set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history)
          }
          next
        }
        # support for:
        # - AIEC (Adherent-Invasive E. coli)
        # - ATEC (Atypical Entero-pathogenic E. coli)
        # - DAEC (Diffusely Adhering E. coli)
        # - EAEC (Entero-Aggresive E. coli)
        # - EHEC (Entero-Haemorrhagic E. coli)
        # - EIEC (Entero-Invasive E. coli)
        # - EPEC (Entero-Pathogenic E. coli)
        # - ETEC (Entero-Toxigenic E. coli)
        # - NMEC (Neonatal Meningitis‐causing E. coli)
        # - STEC (Shiga-toxin producing E. coli)
        # - UPEC (Uropathogenic E. coli)
        if (toupper(x_backup_without_spp[i]) %in% c("AIEC", "ATEC", "DAEC", "EAEC", "EHEC", "EIEC", "EPEC", "ETEC", "NMEC", "STEC", "UPEC")
            | x_backup_without_spp[i] %like% "O?(26|103|104|104|111|121|145|157)") {
          x[i] <- microorganismsDT[mo == 'B_ESCHR_COL', ..property][[1]][1L]
          if (initial_search == TRUE) {
            set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history)
          }
          next
        }
        if (toupper(x_backup_without_spp[i]) == 'MRPA') {
          # multi resistant P. aeruginosa
          x[i] <- microorganismsDT[mo == 'B_PSDMN_AER', ..property][[1]][1L]
          if (initial_search == TRUE) {
            set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history)
          }
          next
        }
        if (toupper(x_backup_without_spp[i]) == 'CRS'
            | toupper(x_backup_without_spp[i]) == 'CRSM') {
          # co-trim resistant S. maltophilia
          x[i] <- microorganismsDT[mo == 'B_STNTR_MAL', ..property][[1]][1L]
          if (initial_search == TRUE) {
            set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history)
          }
          next
        }
        if (toupper(x_backup_without_spp[i]) %in% c('PISP', 'PRSP', 'VISP', 'VRSP')) {
          # peni I, peni R, vanco I, vanco R: S. pneumoniae
          x[i] <- microorganismsDT[mo == 'B_STRPT_PNE', ..property][[1]][1L]
          if (initial_search == TRUE) {
            set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history)
          }
          next
        }
        if (x_backup_without_spp[i] %like% '^G[ABCDFGHK]S$') {
          # Streptococci, like GBS = Group B Streptococci (B_STRPT_GRB)
          x[i] <- microorganismsDT[mo == gsub("G([ABCDFGHK])S", "B_STRPT_GR\\1", x_backup_without_spp[i], ignore.case = TRUE), ..property][[1]][1L]
          if (initial_search == TRUE) {
            set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history)
          }
          next
        }
        if (x_backup_without_spp[i] %like% '(streptococ|streptokok).* [ABCDFGHK]$') {
          # Streptococci in different languages, like "estreptococos grupo B"
          x[i] <- microorganismsDT[mo == gsub(".*(streptococ|streptokok|estreptococ).* ([ABCDFGHK])$", "B_STRPT_GR\\2", x_backup_without_spp[i], ignore.case = TRUE), ..property][[1]][1L]
          if (initial_search == TRUE) {
            set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history)
          }
          next
        }
        if (x_backup_without_spp[i] %like% 'group [ABCDFGHK] (streptococ|streptokok|estreptococ)') {
          # Streptococci in different languages, like "Group A Streptococci"
          x[i] <- microorganismsDT[mo == gsub(".*group ([ABCDFGHK]) (streptococ|streptokok|estreptococ).*", "B_STRPT_GR\\1", x_backup_without_spp[i], ignore.case = TRUE), ..property][[1]][1L]
          if (initial_search == TRUE) {
            set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history)
          }
          next
        }
        if (x_backup_without_spp[i] %like% 'haemoly.*strept') {
          # Haemolytic streptococci in different languages
          x[i] <- microorganismsDT[mo == 'B_STRPT_HAE', ..property][[1]][1L]
          if (initial_search == TRUE) {
            set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history)
          }
          next
        }
        # CoNS/CoPS in different languages (support for German, Dutch, Spanish, Portuguese) ----
        if (x_backup_without_spp[i] %like% '[ck]oagulas[ea] negatie?[vf]'
            | x_trimmed[i] %like% '[ck]oagulas[ea] negatie?[vf]'
            | x_backup_without_spp[i] %like% '[ck]o?ns[^a-z]?$') {
          # coerce S. coagulase negative
          x[i] <- microorganismsDT[mo == 'B_STPHY_CNS', ..property][[1]][1L]
          if (initial_search == TRUE) {
            set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history)
          }
          next
        }
        if (x_backup_without_spp[i] %like% '[ck]oagulas[ea] positie?[vf]'
            | x_trimmed[i] %like% '[ck]oagulas[ea] positie?[vf]'
            | x_backup_without_spp[i] %like% '[ck]o?ps[^a-z]?$') {
          # coerce S. coagulase positive
          x[i] <- microorganismsDT[mo == 'B_STPHY_CPS', ..property][[1]][1L]
          if (initial_search == TRUE) {
            set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history)
          }
          next
        }
        if (x_backup_without_spp[i] %like% 'gram[ -]?neg.*'
            | x_backup_without_spp[i] %like% 'negatie?[vf]'
            | x_trimmed[i] %like% 'gram[ -]?neg.*') {
          # coerce Gram negatives
          x[i] <- microorganismsDT[mo == 'B_GRAMN', ..property][[1]][1L]
          if (initial_search == TRUE) {
            set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history)
          }
          next
        }
        if (x_backup_without_spp[i] %like% 'gram[ -]?pos.*'
            | x_backup_without_spp[i] %like% 'positie?[vf]'
            | x_trimmed[i] %like% 'gram[ -]?pos.*') {
          # coerce Gram positives
          x[i] <- microorganismsDT[mo == 'B_GRAMP', ..property][[1]][1L]
          if (initial_search == TRUE) {
            set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history)
          }
          next
        }
        if (x_backup_without_spp[i] %like% "salmonella [a-z]+ ?.*") {
          if (x_backup_without_spp[i] %like% "Salmonella group") {
            # Salmonella Group A to Z, just return S. species for now
            x[i] <- microorganismsDT[mo == 'B_SLMNL', ..property][[1]][1L]
            if (initial_search == TRUE) {
              set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history)
            }
          } else if (grepl("[sS]almonella [A-Z][a-z]+ ?.*", x_backup_without_spp[i], ignore.case = FALSE)) {
            # Salmonella with capital letter species like "Salmonella Goettingen" - they're all S. enterica
            x[i] <- microorganismsDT[mo == 'B_SLMNL_ENT', ..property][[1]][1L]
            if (initial_search == TRUE) {
              set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history)
            }
            uncertainties <- rbind(uncertainties,
                                   data.frame(uncertainty = 1,
                                              input = x_backup_without_spp[i],
                                              fullname = microorganismsDT[mo == "B_SLMNL_ENT", fullname][[1]],
                                              mo = "B_SLMNL_ENT"))
          }
          next
        }
        
        # trivial names known to the field:
        if ("meningococcus" %like% x_trimmed[i]) {
          # coerce S. coagulase positive
          x[i] <- microorganismsDT[mo == 'B_NESSR_MEN', ..property][[1]][1L]
          if (initial_search == TRUE) {
            set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history)
          }
          next
        }
        if ("gonococcus" %like% x_trimmed[i]) {
          # coerce S. coagulase positive
          x[i] <- microorganismsDT[mo == 'B_NESSR_GON', ..property][[1]][1L]
          if (initial_search == TRUE) {
            set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history)
          }
          next
        }
        if ("pneumococcus" %like% x_trimmed[i]) {
          # coerce S. coagulase positive
          x[i] <- microorganismsDT[mo == 'B_STRPT_PNE', ..property][[1]][1L]
          if (initial_search == TRUE) {
            set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history)
          }
          next
        }
      }

      # FIRST TRY FULLNAMES AND CODES ----
      # if only genus is available, return only genus
      if (all(!c(x[i], x_trimmed[i]) %like% " ")) {
        found <- microorganismsDT[fullname_lower %in% tolower(c(x_species[i], x_trimmed_species[i])), ..property][[1]]
        if (length(found) > 0) {
          x[i] <- found[1L]
          if (initial_search == TRUE) {
            set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history)
          }
          next
        }
        if (nchar(x_backup_without_spp[i]) >= 6) {
          found <- microorganismsDT[fullname_lower %like% paste0("^", unregex(x_backup_without_spp[i]), "[a-z]+"), ..property][[1]]
          if (length(found) > 0) {
            x[i] <- found[1L]
            if (initial_search == TRUE) {
              set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history)
            }
            next
          }
        }
        # rest of genus only is in allow_uncertain part.
      }

      # TRY OTHER SOURCES ----
      # WHONET and other common LIS codes
      if (toupper(x_backup[i]) %in% AMR::microorganisms.codes[, 1]) {
        mo_found <- AMR::microorganisms.codes[toupper(x_backup[i]) == AMR::microorganisms.codes[, 1], "mo"][1L]
        if (length(mo_found) > 0) {
          x[i] <- microorganismsDT[mo == mo_found, ..property][[1]][1L]
          if (initial_search == TRUE) {
            set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history)
          }
          next
        }
      }
      if (!is.null(reference_df)) {
        # self-defined reference
        if (x_backup[i] %in% reference_df[, 1]) {
          ref_mo <- reference_df[reference_df[, 1] == x_backup[i], "mo"]
          if (ref_mo %in% microorganismsDT[, mo]) {
            x[i] <- microorganismsDT[mo == ref_mo, ..property][[1]][1L]
            next
          } else {
            warning("Value '", x_backup[i], "' was found in reference_df, but '", ref_mo, "' is not a valid MO code.", call. = FALSE)
          }
        }
      }

      # allow no codes less than 4 characters long, was already checked for WHONET above
      if (nchar(x_backup_without_spp[i]) < 4) {
        x[i] <- microorganismsDT[mo == "UNKNOWN", ..property][[1]]
        if (initial_search == TRUE) {
          failures <- c(failures, x_backup[i])
          set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history)
        }
        next
      }

      check_per_prevalence <- function(data_to_check,
                                       a.x_backup,
                                       b.x_trimmed,
                                       c.x_trimmed_without_group,
                                       d.x_withspaces_start_end,
                                       e.x_withspaces_start_only,
                                       f.x_withspaces_end_only,
                                       g.x_backup_without_spp) {

        # try probable: trimmed version of fullname ----
        found <- data_to_check[fullname_lower %in% tolower(g.x_backup_without_spp), ..property][[1]]
        if (length(found) > 0) {
          return(found[1L])
        }

        # try any match keeping spaces ----
        found <- data_to_check[fullname %like% d.x_withspaces_start_end, ..property][[1]]
        if (length(found) > 0 & nchar(g.x_backup_without_spp) >= 6) {
          return(found[1L])
        }

        # try any match keeping spaces, not ending with $ ----
        found <- data_to_check[fullname %like% paste0(trimws(e.x_withspaces_start_only), " "), ..property][[1]]
        if (length(found) > 0) {
          return(found[1L])
        }
        found <- data_to_check[fullname %like% e.x_withspaces_start_only, ..property][[1]]
        if (length(found) > 0 & nchar(g.x_backup_without_spp) >= 6) {
          return(found[1L])
        }

        # try any match keeping spaces, not start with ^ ----
        found <- data_to_check[fullname %like% paste0(" ", trimws(f.x_withspaces_end_only)), ..property][[1]]
        if (length(found) > 0) {
          return(found[1L])
        }

        # try a trimmed version
        found <- data_to_check[fullname_lower %like% b.x_trimmed
                               | fullname_lower %like% c.x_trimmed_without_group, ..property][[1]]
        if (length(found) > 0 & nchar(g.x_backup_without_spp) >= 6) {
          return(found[1L])
        }


        # try splitting of characters in the middle and then find ID ----
        # only when text length is 6 or lower
        # like esco = E. coli, klpn = K. pneumoniae, stau = S. aureus, staaur = S. aureus
        if (nchar(g.x_backup_without_spp) <= 6) {
          x_length <- nchar(g.x_backup_without_spp)
          x_split <- paste0("^",
                            g.x_backup_without_spp %>% substr(1, x_length / 2),
                            '.* ',
                            g.x_backup_without_spp %>% substr((x_length / 2) + 1, x_length))
          found <- data_to_check[fullname %like% x_split, ..property][[1]]
          if (length(found) > 0) {
            return(found[1L])
          }
        }

        # try fullname without start and without nchar limit of >= 6 ----
        # like "K. pneu rhino" >> "Klebsiella pneumoniae (rhinoscleromatis)" = KLEPNERH
        found <- data_to_check[fullname %like% e.x_withspaces_start_only, ..property][[1]]
        if (length(found) > 0) {
          return(found[1L])
        }

        # didn't found any
        return(NA_character_)
      }

      # FIRST TRY VERY PREVALENT IN HUMAN INFECTIONS ----
      x[i] <- check_per_prevalence(data_to_check = microorganismsDT[prevalence == 1],
                                   a.x_backup = x_backup[i],
                                   b.x_trimmed = x_trimmed[i],
                                   c.x_trimmed_without_group = x_trimmed_without_group[i],
                                   d.x_withspaces_start_end = x_withspaces_start_end[i],
                                   e.x_withspaces_start_only = x_withspaces_start_only[i],
                                   f.x_withspaces_end_only = x_withspaces_end_only[i],
                                   g.x_backup_without_spp =  x_backup_without_spp[i])
      if (!empty_result(x[i])) {
        if (initial_search == TRUE) {
          set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history)
        }
        next
      }
      # THEN TRY PREVALENT IN HUMAN INFECTIONS ----
      x[i] <- check_per_prevalence(data_to_check = microorganismsDT[prevalence == 2],
                                   a.x_backup = x_backup[i],
                                   b.x_trimmed = x_trimmed[i],
                                   c.x_trimmed_without_group = x_trimmed_without_group[i],
                                   d.x_withspaces_start_end = x_withspaces_start_end[i],
                                   e.x_withspaces_start_only = x_withspaces_start_only[i],
                                   f.x_withspaces_end_only = x_withspaces_end_only[i],
                                   g.x_backup_without_spp =  x_backup_without_spp[i])
      if (!empty_result(x[i])) {
        if (initial_search == TRUE) {
          set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history)
        }
        next
      }
      # THEN UNPREVALENT IN HUMAN INFECTIONS ----
      x[i] <- check_per_prevalence(data_to_check = microorganismsDT[prevalence == 3],
                                   a.x_backup = x_backup[i],
                                   b.x_trimmed = x_trimmed[i],
                                   c.x_trimmed_without_group = x_trimmed_without_group[i],
                                   d.x_withspaces_start_end = x_withspaces_start_end[i],
                                   e.x_withspaces_start_only = x_withspaces_start_only[i],
                                   f.x_withspaces_end_only = x_withspaces_end_only[i],
                                   g.x_backup_without_spp =  x_backup_without_spp[i])
      if (!empty_result(x[i])) {
        if (initial_search == TRUE) {
          set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history)
        }
        next
      }

      # MISCELLANEOUS ----

      # look for old taxonomic names ----
      found <- microorganisms.oldDT[fullname_lower == tolower(x_backup[i])
                                    | fullname %like% x_withspaces_start_end[i],]
      if (NROW(found) > 0) {
        col_id_new <- found[1, col_id_new]
        # when property is "ref" (which is the case in mo_ref, mo_authors and mo_year), return the old value, so:
        # mo_ref("Chlamydia psittaci") = "Page, 1968" (with warning)
        # mo_ref("Chlamydophila psittaci") = "Everett et al., 1999"
        if (property == "ref") {
          x[i] <- found[1, ref]
        } else {
          x[i] <- microorganismsDT[col_id == found[1, col_id_new], ..property][[1]]
        }
        was_renamed(name_old = found[1, fullname],
                    name_new = microorganismsDT[col_id == found[1, col_id_new], fullname],
                    ref_old = found[1, ref],
                    ref_new = microorganismsDT[col_id == found[1, col_id_new], ref],
                    mo = microorganismsDT[col_id == found[1, col_id_new], mo])
        if (initial_search == TRUE) {
          set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history)
        }
        next
      }

      # check for uncertain results ----
      uncertain_fn <- function(a.x_backup,
                               b.x_trimmed,
                               c.x_withspaces_start_end,
                               d.x_withspaces_start_only,
                               f.x_withspaces_end_only,
                               g.x_backup_without_spp) {

        if (uncertainty_level == 0) {
          # do not allow uncertainties
          return(NA_character_)
        }

        if (uncertainty_level >= 1) {
          # (1) look again for old taxonomic names, now for G. species ----
          if (isTRUE(debug)) {
            cat("\n[UNCERTAINLY LEVEL 1] (1) look again for old taxonomic names, now for G. species\n")
          }
          found <- microorganisms.oldDT[fullname %like% c.x_withspaces_start_end
                                        | fullname %like% d.x_withspaces_start_only]
          if (NROW(found) > 0 & nchar(g.x_backup_without_spp) >= 6) {
            if (property == "ref") {
              # when property is "ref" (which is the case in mo_ref, mo_authors and mo_year), return the old value, so:
              # mo_ref("Chlamydia psittaci) = "Page, 1968" (with warning)
              # mo_ref("Chlamydophila psittaci) = "Everett et al., 1999"
              x <- found[1, ref]
            } else {
              x <- microorganismsDT[col_id == found[1, col_id_new], ..property][[1]]
            }
            was_renamed(name_old = found[1, fullname],
                        name_new = microorganismsDT[col_id == found[1, col_id_new], fullname],
                        ref_old = found[1, ref],
                        ref_new = microorganismsDT[col_id == found[1, col_id_new], ref],
                        mo = microorganismsDT[col_id == found[1, col_id_new], mo])
            uncertainties <<- rbind(uncertainties,
                                    data.frame(uncertainty = 1,
                                               input = a.x_backup,
                                               fullname = found[1, fullname],
                                               mo = paste("CoL", found[1, col_id])))
            if (initial_search == TRUE) {
              set_mo_history(a.x_backup, get_mo_code(x, property), 1, force = force_mo_history)
            }
            return(x)
          }

          # (2) Try with misspelled input ----
          # just rerun with initial_search = FALSE will used the extensive regex part above
          if (isTRUE(debug)) {
            cat("\n[UNCERTAINLY LEVEL 1] (2) Try with misspelled input\n")
          }
          found <- suppressMessages(suppressWarnings(exec_as.mo(a.x_backup, initial_search = FALSE, allow_uncertain = FALSE, debug = debug)))
          if (!empty_result(found)) {
            found_result <- found
            found <- microorganismsDT[mo == found, ..property][[1]]
            uncertainties <<- rbind(uncertainties,
                                    data.frame(uncertainty = 1,
                                               input = a.x_backup,
                                               fullname = microorganismsDT[mo == found_result[1L], fullname][[1]],
                                               mo = found_result[1L]))
            if (initial_search == TRUE) {
              set_mo_history(a.x_backup, get_mo_code(found[1L], property), 1, force = force_mo_history)
            }
            return(found[1L])
          }
        }

        if (uncertainty_level >= 2) {

          # (3) look for genus only, part of name ----
          if (isTRUE(debug)) {
            cat("\n[UNCERTAINLY LEVEL 2] (3) look for genus only, part of name\n")
          }
          if (nchar(g.x_backup_without_spp) > 4 & !b.x_trimmed %like% " ") {
            if (!grepl("^[A-Z][a-z]+", b.x_trimmed, ignore.case = FALSE)) {
              # not when input is like Genustext, because then Neospora would lead to Actinokineospora
              found <- microorganismsDT[fullname_lower %like% paste(b.x_trimmed, "species"), ..property][[1]]
              if (length(found) > 0) {
                x[i] <- found[1L]
                uncertainties <<- rbind(uncertainties,
                                        data.frame(uncertainty = 2,
                                                   input = a.x_backup,
                                                   fullname = microorganismsDT[mo == found[1L], fullname][[1]],
                                                   mo = found[1L]))
                if (initial_search == TRUE) {
                  set_mo_history(a.x_backup, get_mo_code(x, property), 2, force = force_mo_history)
                }
                return(x)
              }
            }
          }

          # (4) strip values between brackets ----
          if (isTRUE(debug)) {
            cat("\n[UNCERTAINLY LEVEL 2] (4) strip values between brackets\n")
          }
          a.x_backup_stripped <- gsub("( *[(].*[)] *)", " ", a.x_backup)
          a.x_backup_stripped <- trimws(gsub(" +", " ", a.x_backup_stripped))
          found <- suppressMessages(suppressWarnings(exec_as.mo(a.x_backup_stripped, initial_search = FALSE, allow_uncertain = FALSE, debug = debug)))
          if (!empty_result(found) & nchar(g.x_backup_without_spp) >= 6) {
            found_result <- found
            found <- microorganismsDT[mo == found, ..property][[1]]
            uncertainties <<- rbind(uncertainties,
                                    data.frame(uncertainty = 2,
                                               input = a.x_backup,
                                               fullname = microorganismsDT[mo == found_result[1L], fullname][[1]],
                                               mo = found_result[1L]))
            if (initial_search == TRUE) {
              set_mo_history(a.x_backup, get_mo_code(found[1L], property), 2, force = force_mo_history)
            }
            return(found[1L])
          }

          # (5a) try to strip off half an element from end and check the remains ----
          if (isTRUE(debug)) {
            cat("\n[UNCERTAINLY LEVEL 2] (5a) try to strip off half an element from end and check the remains\n")
          }
          x_strip <- a.x_backup %>% strsplit(" ") %>% unlist()
          if (length(x_strip) > 1) {
            for (i in 1:(length(x_strip) - 1)) {
              lastword <- x_strip[length(x_strip) - i + 1]
              lastword_half <- substr(lastword, 1, as.integer(nchar(lastword) / 2))
              # remove last half of the second term
              x_strip_collapsed <- paste(c(x_strip[1:(length(x_strip) - i)], lastword_half), collapse = " ")
              if (nchar(x_strip_collapsed) >= 4) {
                found <- suppressMessages(suppressWarnings(exec_as.mo(x_strip_collapsed, initial_search = FALSE, allow_uncertain = FALSE, debug = debug)))
                if (!empty_result(found)) {
                  found_result <- found
                  found <- microorganismsDT[mo == found, ..property][[1]]
                  uncertainties <<- rbind(uncertainties,
                                          data.frame(uncertainty = 2,
                                                     input = a.x_backup,
                                                     fullname = microorganismsDT[mo == found_result[1L], fullname][[1]],
                                                     mo = found_result[1L]))
                  if (initial_search == TRUE) {
                    set_mo_history(a.x_backup, get_mo_code(found[1L], property), 2, force = force_mo_history)
                  }
                  return(found[1L])
                }
              }
            }
          }
          # (5b) try to strip off one element from end and check the remains ----
          if (isTRUE(debug)) {
            cat("\n[UNCERTAINLY LEVEL 2] (5b) try to strip off one element from end and check the remains\n")
          }
          if (length(x_strip) > 1) {
            for (i in 1:(length(x_strip) - 1)) {
              x_strip_collapsed <- paste(x_strip[1:(length(x_strip) - i)], collapse = " ")
              if (nchar(x_strip_collapsed) >= 4) {
                found <- suppressMessages(suppressWarnings(exec_as.mo(x_strip_collapsed, initial_search = FALSE, allow_uncertain = FALSE, debug = debug)))
                if (!empty_result(found)) {
                  found_result <- found
                  found <- microorganismsDT[mo == found, ..property][[1]]
                  uncertainties <<- rbind(uncertainties,
                                          data.frame(uncertainty = 2,
                                                     input = a.x_backup,
                                                     fullname = microorganismsDT[mo == found_result[1L], fullname][[1]],
                                                     mo = found_result[1L]))
                  if (initial_search == TRUE) {
                    set_mo_history(a.x_backup, get_mo_code(found[1L], property), 2, force = force_mo_history)
                  }
                  return(found[1L])
                }
              }
            }
          }
          # (5c) check for unknown yeasts/fungi ----
          if (isTRUE(debug)) {
            cat("\n[UNCERTAINLY LEVEL 2] (5b) check for unknown yeasts/fungi\n")
          }
          if (b.x_trimmed %like% "yeast") {
            found <- "F_YEAST"
            found_result <- found
            found <- microorganismsDT[mo == found, ..property][[1]]
            uncertainties <<- rbind(uncertainties,
                                    data.frame(uncertainty = 2,
                                               input = a.x_backup,
                                               fullname = microorganismsDT[mo == found_result[1L], fullname][[1]],
                                               mo = found_result[1L]))
            if (initial_search == TRUE) {
              set_mo_history(a.x_backup, get_mo_code(found[1L], property), 2, force = force_mo_history)
            }
            return(found[1L])
          }
          if (b.x_trimmed %like% "fungus") {
            found <- "F_FUNGUS"
            found_result <- found
            found <- microorganismsDT[mo == found, ..property][[1]]
            uncertainties <<- rbind(uncertainties,
                                    data.frame(uncertainty = 2,
                                               input = a.x_backup,
                                               fullname = microorganismsDT[mo == found_result[1L], fullname][[1]],
                                               mo = found_result[1L]))
            if (initial_search == TRUE) {
              set_mo_history(a.x_backup, get_mo_code(found[1L], property), 2, force = force_mo_history)
            }
            return(found[1L])
          }
          # (6) try to strip off one element from start and check the remains (only allow >= 2-part name outcome) ----
          if (isTRUE(debug)) {
            cat("\n[UNCERTAINLY LEVEL 2] (6) try to strip off one element from start and check the remains (only allow >= 2-part name outcome)\n")
          }
          x_strip <- a.x_backup %>% strsplit(" ") %>% unlist()
          if (length(x_strip) > 1 & nchar(g.x_backup_without_spp) >= 6) {
            for (i in 2:(length(x_strip))) {
              x_strip_collapsed <- paste(x_strip[i:length(x_strip)], collapse = " ")
              found <- suppressMessages(suppressWarnings(exec_as.mo(x_strip_collapsed, initial_search = FALSE, allow_uncertain = FALSE, debug = debug)))
              if (!empty_result(found)) {
                found_result <- found
                found <- microorganismsDT[mo == found_result[1L], ..property][[1]]
                # uncertainty level 2 only if searched part contains a space (otherwise it will be found with lvl 3)
                if (x_strip_collapsed %like% " ") {
                  uncertainties <<- rbind(uncertainties,
                                          data.frame(uncertainty = 2,
                                                     input = a.x_backup,
                                                     fullname = microorganismsDT[mo == found_result[1L], fullname][[1]],
                                                     mo = found_result[1L]))
                  if (initial_search == TRUE) {
                    set_mo_history(a.x_backup, get_mo_code(found[1L], property), 2, force = force_mo_history)
                  }
                  return(found[1L])
                }
              }
            }
          }
        }

        if (uncertainty_level >= 3) {
          # (7) try to strip off one element from start and check the remains ----
          if (isTRUE(debug)) {
            cat("\n[UNCERTAINLY LEVEL 3] (7) try to strip off one element from start and check the remains\n")
          }
          x_strip <- a.x_backup %>% strsplit(" ") %>% unlist()
          if (length(x_strip) > 1 & nchar(g.x_backup_without_spp) >= 6) {
            for (i in 2:(length(x_strip))) {
              x_strip_collapsed <- paste(x_strip[i:length(x_strip)], collapse = " ")
              found <- suppressMessages(suppressWarnings(exec_as.mo(x_strip_collapsed, initial_search = FALSE, allow_uncertain = FALSE, debug = debug)))
              if (!empty_result(found)) {
                found_result <- found
                found <- microorganismsDT[mo == found, ..property][[1]]
                uncertainties <<- rbind(uncertainties,
                                        data.frame(uncertainty = 3,
                                                   input = a.x_backup,
                                                   fullname = microorganismsDT[mo == found_result[1L], fullname][[1]],
                                                   mo = found_result[1L]))
                if (initial_search == TRUE) {
                  set_mo_history(a.x_backup, get_mo_code(found[1L], property), 3, force = force_mo_history)
                }
                return(found[1L])
              }
            }
          }

          # (8) part of a name (very unlikely match) ----
          if (isTRUE(debug)) {
            cat("\n[UNCERTAINLY LEVEL 3] (8) part of a name (very unlikely match)\n")
          }
          found <- microorganismsDT[fullname %like% f.x_withspaces_end_only]
          if (nrow(found) > 0) {
            found_result <- found[["mo"]]
            if (!empty_result(found_result) & nchar(g.x_backup_without_spp) >= 6) {
              found <- microorganismsDT[mo == found_result[1L], ..property][[1]]
              uncertainties <<- rbind(uncertainties,
                                      data.frame(uncertainty = 3,
                                                 input = a.x_backup,
                                                 fullname = microorganismsDT[mo == found_result[1L], fullname][[1]],
                                                 mo = found_result[1L]))
              if (initial_search == TRUE) {
                set_mo_history(a.x_backup, get_mo_code(found[1L], property), 3, force = force_mo_history)
              }
              return(found[1L])
            }
          }
        }

        # didn't found in uncertain results too
        return(NA_character_)
      }
      x[i] <- uncertain_fn(x_backup[i],
                           x_trimmed[i],
                           x_withspaces_start_end[i],
                           x_withspaces_start_only[i],
                           x_withspaces_end_only[i],
                           x_backup_without_spp[i])
      if (!empty_result(x[i])) {
        # no set_mo_history here - it is already set in uncertain_fn()
        next
      }

      # no results found: make them UNKNOWN ----
      x[i] <- microorganismsDT[mo == "UNKNOWN", ..property][[1]]
      if (initial_search == TRUE) {
        failures <- c(failures, x_backup[i])
        set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history)
      }
    }
  }

  # handling failures ----
  failures <- failures[!failures %in% c(NA, NULL, NaN)]
  if (length(failures) > 0 & initial_search == TRUE) {
    options(mo_failures = sort(unique(failures)))
    plural <- c("value", "it", "was")
    if (n_distinct(failures) > 1) {
      plural <- c("values", "them", "were")
    }
    total_failures <- length(x_input[as.character(x_input) %in% as.character(failures) & !x_input %in% c(NA, NULL, NaN)])
    total_n <- length(x_input[!x_input %in% c(NA, NULL, NaN)])
    msg <- paste0(nr2char(n_distinct(failures)), " unique ", plural[1],
                  " (covering ", percent(total_failures / total_n, round = 1, force_zero = TRUE),
                  ") could not be coerced and ", plural[3], " considered 'unknown'")
    if (n_distinct(failures) <= 10) {
      msg <- paste0(msg, ": ", paste('"', unique(failures), '"', sep = "", collapse = ', '))
    }
    msg <- paste0(msg,  ". Use mo_failures() to review ", plural[2], ". Edit the `allow_uncertain` parameter if needed (see ?as.mo).")
    warning(red(msg),
            call. = FALSE,
            immediate. = TRUE) # thus will always be shown, even if >= warnings
  }
  # handling uncertainties ----
  if (NROW(uncertainties) > 0 & initial_search == TRUE) {
    options(mo_uncertainties = as.list(distinct(uncertainties, input, .keep_all = TRUE)))

    plural <- c("", "it")
    if (NROW(uncertainties) > 1) {
      plural <- c("s", "them")
    }
    msg <- paste0("\nResult", plural[1], " of ", nr2char(NROW(uncertainties)), " value", plural[1],
                  " was guessed with uncertainty. Use mo_uncertainties() to review ", plural[2], ".")
    warning(red(msg),
            call. = FALSE,
            immediate. = TRUE) # thus will always be shown, even if >= warnings
  }

  # Becker ----
  if (Becker == TRUE | Becker == "all") {
    # See Source. It's this figure:
    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4187637/figure/F3/
    MOs_staph <- microorganismsDT[genus == "Staphylococcus"]
    setkey(MOs_staph, species)
    CoNS <- MOs_staph[species %in% c("arlettae", "auricularis", "capitis",
                                     "caprae", "carnosus", "chromogenes", "cohnii", "condimenti",
                                     "devriesei", "epidermidis", "equorum", "felis",
                                     "fleurettii", "gallinarum", "haemolyticus",
                                     "hominis", "jettensis", "kloosii", "lentus",
                                     "lugdunensis", "massiliensis", "microti",
                                     "muscae", "nepalensis", "pasteuri", "petrasii",
                                     "pettenkoferi", "piscifermentans", "rostri",
                                     "saccharolyticus", "saprophyticus", "sciuri",
                                     "stepanovicii", "simulans", "succinus",
                                     "vitulinus", "warneri", "xylosus")
                      | (species == "schleiferi" & subspecies %in% c("schleiferi", "")), ..property][[1]]
    CoPS <- MOs_staph[species %in% c("simiae", "agnetis",
                                     "delphini", "lutrae",
                                     "hyicus", "intermedius",
                                     "pseudintermedius", "pseudointermedius",
                                     "schweitzeri", "argenteus")
                      | (species == "schleiferi" & subspecies == "coagulans"), ..property][[1]]

    # warn when species found that are not in Becker (2014, PMID 25278577) and Becker (2019, PMID 30872103)
    post_Becker <- c("argensis", "caeli", "cornubiensis", "edaphicus")
    if (any(x %in% MOs_staph[species %in% post_Becker, ..property][[1]])) {

      warning("Becker ", italic("et al."), " (2014, 2019) does not contain these species named after their publication: ",
              italic(paste("S.",
                           sort(mo_species(unique(x[x %in% MOs_staph[species %in% post_Becker, ..property][[1]]]))),
                           collapse = ", ")),
              ".",
              call. = FALSE,
              immediate. = TRUE)
    }

    x[x %in% CoNS] <- microorganismsDT[mo == 'B_STPHY_CNS', ..property][[1]][1L]
    x[x %in% CoPS] <- microorganismsDT[mo == 'B_STPHY_CPS', ..property][[1]][1L]
    if (Becker == "all") {
      x[x %in% microorganismsDT[mo %like% '^B_STPHY_AUR', ..property][[1]]] <- microorganismsDT[mo == 'B_STPHY_CPS', ..property][[1]][1L]
    }
  }

  # Lancefield ----
  if (Lancefield == TRUE | Lancefield == "all") {
    # group A - S. pyogenes
    x[x == microorganismsDT[mo == 'B_STRPT_PYO', ..property][[1]][1L]] <- microorganismsDT[mo == 'B_STRPT_GRA', ..property][[1]][1L]
    # group B - S. agalactiae
    x[x == microorganismsDT[mo == 'B_STRPT_AGA', ..property][[1]][1L]] <- microorganismsDT[mo == 'B_STRPT_GRB', ..property][[1]][1L]
    # group C
    S_groupC <- microorganismsDT %>% filter(genus == "Streptococcus",
                                            species %in% c("equisimilis", "equi",
                                                           "zooepidemicus", "dysgalactiae")) %>%
      pull(property)
    x[x %in% S_groupC] <- microorganismsDT[mo == 'B_STRPT_GRC', ..property][[1]][1L]
    if (Lancefield == "all") {
      # all Enterococci
      x[x %like% "^(Enterococcus|B_ENTRC)"] <- microorganismsDT[mo == 'B_STRPT_GRD', ..property][[1]][1L]
    }
    # group F - S. anginosus
    x[x == microorganismsDT[mo == 'B_STRPT_ANG', ..property][[1]][1L]] <- microorganismsDT[mo == 'B_STRPT_GRF', ..property][[1]][1L]
    # group H - S. sanguinis
    x[x == microorganismsDT[mo == 'B_STRPT_SAN', ..property][[1]][1L]] <- microorganismsDT[mo == 'B_STRPT_GRH', ..property][[1]][1L]
    # group K - S. salivarius
    x[x == microorganismsDT[mo == 'B_STRPT_SAL', ..property][[1]][1L]] <- microorganismsDT[mo == 'B_STRPT_GRK', ..property][[1]][1L]
  }

  # Wrap up ----------------------------------------------------------------

  # comply to x, which is also unique and without empty values
  x_input_unique_nonempty <- unique(x_input[!is.na(x_input)
                                            & !is.null(x_input)
                                            & !identical(x_input, "")
                                            & !identical(x_input, "xxx")])

  # left join the found results to the original input values (x_input)
  df_found <- data.frame(input = as.character(x_input_unique_nonempty),
                         found = as.character(x),
                         stringsAsFactors = FALSE)
  df_input <- data.frame(input = as.character(x_input),
                         stringsAsFactors = FALSE)

  suppressWarnings(
    x <- df_input %>%
      left_join(df_found,
                by = "input") %>%
      pull(found)
  )

  if (property == "mo") {
    class(x) <- "mo"
  }

  if (length(mo_renamed()) > 0) {
    print(mo_renamed())
  }

  x
}

empty_result <- function(x) {
  all(x %in% c(NA, "UNKNOWN"))
}

#' @importFrom crayon italic
was_renamed <- function(name_old, name_new, ref_old = "", ref_new = "", mo = "") {
  if (!is.na(ref_old)) {
    ref_old <- paste0(" (", ref_old, ")")
  } else {
    ref_old <- ""
  }
  if (!is.na(ref_new)) {
    ref_new <- paste0(" (", ref_new, ")")
  } else {
    ref_new <- ""
  }
  if (!is.na(mo)) {
    mo <- paste0(" (", mo, ")")
  } else {
    mo <- ""
  }
  old_values <- paste0(italic(name_old), ref_old)
  old_values <- gsub("et al.", italic("et al."), old_values)
  new_values <- paste0(italic(name_new), ref_new, mo)
  new_values <- gsub("et al.", italic("et al."), new_values)

  names(new_values) <- old_values
  total <- c(getOption("mo_renamed"), new_values)
  options(mo_renamed = total[order(names(total))])
}

#' @exportMethod print.mo
#' @export
#' @noRd
print.mo <- function(x, ...) {
  cat("Class 'mo'\n")
  x_names <- names(x)
  x <- as.character(x)
  names(x) <- x_names
  print.default(x, quote = FALSE)
}

#' @exportMethod summary.mo
#' @importFrom dplyr n_distinct
#' @importFrom clean freq top_freq
#' @export
#' @noRd
summary.mo <- function(object, ...) {
  # unique and top 1-3
  x <- object
  top_3 <- unname(top_freq(freq(x), 3))
  c("Class" = "mo",
    "<NA>" = length(x[is.na(x)]),
    "Unique" = n_distinct(x[!is.na(x)]),
    "#1" = top_3[1],
    "#2" = top_3[2],
    "#3" = top_3[3])
}

#' @exportMethod as.data.frame.mo
#' @export
#' @noRd
as.data.frame.mo <- function(x, ...) {
  # same as as.data.frame.character but with removed stringsAsFactors, since it will be class "mo"
  nm <- paste(deparse(substitute(x), width.cutoff = 500L),
              collapse = " ")
  if (!"nm" %in% names(list(...))) {
    as.data.frame.vector(x, ..., nm = nm)
  } else {
    as.data.frame.vector(x, ...)
  }
}

#' @exportMethod pull.mo
#' @export
#' @importFrom dplyr pull
#' @noRd
pull.mo <- function(.data, ...) {
  pull(as.data.frame(.data), ...)
}

#' @rdname as.mo
#' @export
mo_failures <- function() {
  getOption("mo_failures")
}

#' @rdname as.mo
#' @importFrom crayon italic
#' @export
mo_uncertainties <- function() {
  if (is.null(getOption("mo_uncertainties"))) {
    return(NULL)
  }
  structure(.Data = as.data.frame(getOption("mo_uncertainties"), stringsAsFactors = FALSE),
            class = c("mo_uncertainties", "data.frame"))
}

#' @exportMethod print.mo_uncertainties
#' @importFrom crayon green yellow red white bgGreen bgYellow bgRed
#' @export
#' @noRd
print.mo_uncertainties <- function(x, ...) {
  if (NROW(x) == 0) {
    return(NULL)
  }
  cat(paste0(bold(nr2char(nrow(x)), paste0("unique result", ifelse(nrow(x) > 1, "s", ""), " guessed with uncertainty:")),
             "\n(1 = ", green("renamed/misspelled"),
             ", 2 = ", yellow("uncertain"),
             ", 3 = ", red("very uncertain"), ")\n"))

  msg <- ""
  for (i in 1:nrow(x)) {
    if (x[i, "uncertainty"] == 1) {
      colour1 <- green
      colour2 <- function(...) bgGreen(white(...))
    } else if (x[i, "uncertainty"] == 2) {
      colour1 <- yellow
      colour2 <- bgYellow
    } else {
      colour1 <- red
      colour2 <- function(...) bgRed(white(...))
    }
    msg <- paste(msg,
                 paste0(colour2(paste0(" [", x[i, "uncertainty"], "] ")), ' "', x[i, "input"], '" -> ',
                        colour1(paste0(italic(x[i, "fullname"]), " (", x[i, "mo"], ")"))),
                 sep = "\n")
  }
  cat(msg)
}

#' @rdname as.mo
#' @importFrom crayon strip_style
#' @export
mo_renamed <- function() {
  items <- getOption("mo_renamed")
  if (is.null(items)) {
    return(NULL)
  }

  items <- strip_style(items)
  names(items) <- strip_style(names(items))
  structure(.Data = items,
             class = c("mo_renamed", "character"))
}

#' @exportMethod print.mo_renamed
#' @importFrom crayon blue
#' @export
#' @noRd
print.mo_renamed <- function(x, ...) {
  items <- getOption("mo_renamed")
  base::message(blue(paste("NOTE:", names(items), "was renamed", items, collapse = "\n"), collapse = "\n"))
}

nr2char <- function(x) {
  if (x %in% c(1:10)) {
    v <- c("one" = 1, "two" = 2, "three" = 3, "four" = 4, "five" = 5,
           "six" = 6, "seven" = 7, "eight" = 8, "nine" = 9, "ten" = 10)
    names(v[x])
  } else {
    x
  }
}

unregex <- function(x) {
  gsub("[^a-zA-Z0-9 -]", "", x)
}

get_mo_code <- function(x, property) {
  # don't use right now
  return(NULL)

  if (property == "mo") {
    unique(x)
  } else {
    AMR::microorganisms[base::which(AMR::microorganisms[, property] %in% x),]$mo
  }
}

translate_allow_uncertain <- function(allow_uncertain) {
  if (isTRUE(allow_uncertain)) {
    # default to uncertainty level 2
    allow_uncertain <- 2
  } else {
    allow_uncertain <- as.integer(allow_uncertain)
    if (!allow_uncertain %in% c(0:3)) {
      stop("`allow_uncertain` must be a number between 0 (none) and 3 (all), or TRUE (= 2) or FALSE (= 0).", call. = FALSE)
    }
  }
  allow_uncertain
}

get_mo_failures_uncertainties_renamed <- function() {
  list(failures = getOption("mo_failures"),
       uncertainties = getOption("mo_uncertainties"),
       renamed = getOption("mo_renamed"))
}

load_mo_failures_uncertainties_renamed <- function(metadata) {
  options("mo_failures" = metadata$failures)
  options("mo_uncertainties" = metadata$uncertainties)
  options("mo_renamed" = metadata$renamed)
}
