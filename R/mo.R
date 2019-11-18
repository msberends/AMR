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
#' @param allow_uncertain a number between 0 (or "none") and 3 (or "all"), or TRUE (= 2) or FALSE (= 0) to indicate whether the input should be checked for less probable results, see Details
#' @param reference_df a \code{data.frame} to use for extra reference when translating \code{x} to a valid \code{mo}. See \code{\link{set_mo_source}} and \code{\link{get_mo_source}} to automate the usage of your own codes (e.g. used in your analysis or organisation).
#' @param ... other parameters passed on to functions
#' @rdname as.mo
#' @aliases mo
#' @keywords mo Becker becker Lancefield lancefield guess
#' @details
#' \strong{General info} \cr
#' A microorganism ID from this package (class: \code{mo}) typically looks like these examples:\cr
#' \preformatted{
#'   Code               Full name
#'   ---------------    --------------------------------------
#'   B_KLBSL            Klebsiella
#'   B_KLBSL_PNMN       Klebsiella pneumoniae
#'   B_KLBSL_PNMN_RHNS  Klebsiella pneumoniae rhinoscleromatis
#'   |   |    |    |
#'   |   |    |    |
#'   |   |    |     ---> subspecies, a 4-5 letter acronym
#'   |   |     ----> species, a 4-5 letter acronym
#'   |    ----> genus, a 5-7 letter acronym
#'    ----> taxonomic kingdom: A (Archaea), AN (Animalia), B (Bacteria),
#'                             C (Chromista), F (Fungi), P (Protozoa)
#' }
#'
#' Values that cannot be coered will be considered 'unknown' and will get the MO code \code{UNKNOWN}.
#'
#' Use the \code{\link{mo_property}_*} functions to get properties based on the returned code, see Examples.
#'
#' The algorithm uses data from the Catalogue of Life (see below) and from one other source (see \code{\link{microorganisms}}).
#'
#' The \code{as.mo()} function uses several coercion rules for fast and logical results. It assesses the input matching criteria in the following order:

#' \itemize{
#'   \item{Human pathogenic prevalence: the function  starts with more prevalent microorganisms, followed by less prevalent ones;}
#'   \item{Taxonomic kingdom: the function starts with determining Bacteria, then Fungi, then Protozoa, then others;}
#'   \item{Breakdown of input values to identify possible matches.}
#' }
#'
#' This will lead to the effect that e.g. \code{"E. coli"} (a highly prevalent microorganism found in humans) will return the microbial ID of \emph{Escherichia coli} and not \emph{Entamoeba coli} (a less prevalent microorganism in humans), although the latter would alphabetically come first. 
#' 
#' \strong{Coping with uncertain results} \cr
#' In addition, the \code{as.mo()} function can differentiate four levels of uncertainty to guess valid results: 
#' 
#' \itemize{
#'   \item{Uncertainty level 0: no additional rules are applied;}
#'   \item{Uncertainty level 1: allow previously accepted (but now invalid) taxonomic names and minor spelling errors;}
#'   \item{Uncertainty level 2: allow all of level 1, strip values between brackets, inverse the words of the input, strip off text elements from the end keeping at least two elements;}
#'   \item{Uncertainty level 3: allow all of level 1 and 2, strip off text elements from the end, allow any part of a taxonomic name.}
#' }
#' 
#' This leads to e.g.:
#' 
#' \itemize{
#'   \item{\code{"Streptococcus group B (known as S. agalactiae)"}. The text between brackets will be removed and a warning will be thrown that the result \emph{Streptococcus group B} (\code{B_STRPT_GRPB}) needs review.}
#'   \item{\code{"S. aureus - please mind: MRSA"}. The last word will be stripped, after which the function will try to find a match. If it does not, the second last word will be stripped, etc. Again, a warning will be thrown that the result \emph{Staphylococcus aureus} (\code{B_STPHY_AURS}) needs review.}
#'   \item{\code{"Fluoroquinolone-resistant Neisseria gonorrhoeae"}. The first word will be stripped, after which the function will try to find a match. A warning will be thrown that the result \emph{Neisseria gonorrhoeae} (\code{B_NESSR_GNRR}) needs review.}
#' }
#'
#' The level of uncertainty can be set using the argument \code{allow_uncertain}. The default is \code{allow_uncertain = TRUE}, which is equal to uncertainty level 2. Using \code{allow_uncertain = FALSE} is equal to uncertainty level 0 and will skip all rules. You can also use e.g. \code{as.mo(..., allow_uncertain = 1)} to only allow up to level 1 uncertainty.
#' 
#' There are three helper functions that can be run after then \code{as.mo()} function:
#' \itemize{
#'   \item{Use \code{mo_uncertainties()} to get a \code{data.frame} with all values that were coerced to a valid value, but with uncertainty. The output contains a score, that is calculated as \code{(n - 0.5 * L) / n}, where \emph{n} is the number of characters of the returned full name of the microorganism, and \emph{L} is the \href{https://en.wikipedia.org/wiki/Levenshtein_distance}{Levenshtein distance} between that full name and the user input.}
#'   \item{Use \code{mo_failures()} to get a vector with all values that could not be coerced to a valid value.}
#'   \item{Use \code{mo_renamed()} to get a \code{data.frame} with all values that could be coerced based on an old, previously accepted taxonomic name.}
#' }   
#'
#' \strong{Microbial prevalence of pathogens in humans} \cr
#' The intelligent rules consider the prevalence of microorganisms in humans grouped into three groups, which is available as the \code{prevalence} columns in the \code{\link{microorganisms}} and \code{\link{microorganisms.old}} data sets. The grouping into prevalence groups is based on experience from several microbiological laboratories in the Netherlands in conjunction with international reports on pathogen prevalence.
#' 
#' Group 1 (most prevalent microorganisms) consists of all microorganisms where the taxonomic class is Gammaproteobacteria or where the taxonomic genus is  \emph{Enterococcus}, \emph{Staphylococcus} or \emph{Streptococcus}. This group consequently contains all common Gram-negative bacteria, such as \emph{Pseudomonas} and \emph{Legionella} and all species within the order Enterobacteriales. 
#' 
#' Group 2 consists of all microorganisms where the taxonomic phylum is Proteobacteria, Firmicutes, Actinobacteria or Sarcomastigophora, or where the taxonomic genus is \emph{Aspergillus}, \emph{Bacteroides}, \emph{Candida}, \emph{Capnocytophaga}, \emph{Chryseobacterium}, \emph{Cryptococcus}, \emph{Elisabethkingia}, \emph{Flavobacterium}, \emph{Fusobacterium}, \emph{Giardia}, \emph{Leptotrichia}, \emph{Mycoplasma}, \emph{Prevotella}, \emph{Rhodotorula}, \emph{Treponema}, \emph{Trichophyton} or \emph{Ureaplasma}. 
#' 
#' Group 3 (least prevalent microorganisms) consists of all other microorganisms.
#' 
#' \strong{Self-learning algorithm} \cr
#' The \code{as.mo()} function gains experience from previously determined microorganism IDs and learns from it. This drastically improves both speed and reliability. Use \code{clear_mo_history()} to reset the algorithms. Only experience from your current \code{AMR} package version is used. This is done because in the future the taxonomic tree (which is included in this package) may change for any organism and it consequently has to rebuild its knowledge.
#'
#' Usually, any guess after the first try runs 80-95\% faster than the first try.
#'
# \emph{For now, learning only works per session. If R is closed or terminated, the algorithms reset. This might be resolved in a future version.}
#' This resets with every update of this \code{AMR} package since results are saved to your local package library folder.
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
#' # These examples all return "B_STPHY_AURS", the ID of S. aureus:
#' as.mo("sau") # WHONET code
#' as.mo("stau")
#' as.mo("STAU")
#' as.mo("staaur")
#' as.mo("S. aureus")
#' as.mo("S aureus")
#' as.mo("Staphylococcus aureus")
#' as.mo("Staphylococcus aureus (MRSA)")
#' as.mo("Zthafilokkoockus oureuz") # handles incorrect spelling
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
#' as.mo("S. epidermidis")                 # will remain species: B_STPHY_EPDR
#' as.mo("S. epidermidis", Becker = TRUE)  # will not remain species: B_STPHY_CONS
#'
#' as.mo("S. pyogenes")                    # will remain species: B_STRPT_PYGN
#' as.mo("S. pyogenes", Lancefield = TRUE) # will not remain species: B_STRPT_GRPA
#'
#' # All mo_* functions use as.mo() internally too (see ?mo_property):
#' mo_genus("E. coli")           # returns "Escherichia"
#' mo_gramstain("E. coli")       # returns "Gram negative"
#' 
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
  mo_hist <- get_mo_history(x, 
                            uncertainty_level,
                            force = isTRUE(list(...)$force_mo_history),
                            disable = isTRUE(list(...)$disable_mo_history))
  
  if (mo_source_isvalid(reference_df)
      & isFALSE(Becker)
      & isFALSE(Lancefield)
      & !is.null(reference_df)
      & all(x %in% reference_df[, 1][[1]])) {
    
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
    
  } else if (all(x %in% microorganismsDT$mo)
             & isFALSE(Becker)
             & isFALSE(Lancefield)) {
    y <- x
    
  } else if (!any(is.na(mo_hist))
             & isFALSE(Becker)
             & isFALSE(Lancefield)) {
    # check previously found results
    y <- mo_hist
    
  } else if (all(tolower(x) %in% microorganismsDT$fullname_lower)
             & isFALSE(Becker)
             & isFALSE(Lancefield)) {
    # we need special treatment for very prevalent full names, they are likely! (case insensitive)
    # e.g. as.mo("Staphylococcus aureus")
    y <- data.frame(fullname_lower = tolower(x), 
                    stringsAsFactors = FALSE) %>% 
      left_join(microorganismsDT, by = "fullname_lower") %>% 
      pull(mo)
    
    # don't save valid fullnames to history (i.e. values that are in microorganisms$fullname)
    
  } else {
    # will be checked for mo class in validation and uses exec_as.mo internally if necessary
    y <- mo_validate(x = x, property = "mo",
                     Becker = Becker, Lancefield = Lancefield,
                     allow_uncertain = uncertainty_level, reference_df = reference_df,
                     ...)
  }
  
  
  to_class_mo(y)
}

to_class_mo <- function(x) {
  structure(.Data = x, class = "mo")
}

#' @rdname as.mo
#' @export
is.mo <- function(x) {
  identical(class(x), class(to_class_mo(x)))
}

#' @importFrom dplyr %>% pull left_join n_distinct progress_estimated filter distinct
#' @importFrom data.table data.table as.data.table setkey
#' @importFrom crayon magenta red blue silver italic
#' @importFrom cleaner percentage
# param property a column name of AMR::microorganisms
# param initial_search logical - is FALSE when coming from uncertain tries, which uses exec_as.mo internally too
# param dyslexia_mode logical - also check for characters that resemble others
# param force_mo_history logical - whether found result must be saved with set_mo_history (default FALSE on non-interactive sessions)
# param disable_mo_history logical - whether set_mo_history and get_mo_history should be ignored
# param debug logical - show different lookup texts while searching
# param reference_data_to_use data.frame - the data set to check for
exec_as.mo <- function(x,
                       Becker = FALSE,
                       Lancefield = FALSE,
                       allow_uncertain = TRUE,
                       reference_df = get_mo_source(),
                       property = "mo",
                       initial_search = TRUE,
                       dyslexia_mode = FALSE,
                       force_mo_history = FALSE,
                       disable_mo_history = FALSE,
                       debug = FALSE,
                       reference_data_to_use = microorganismsDT) {

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
  options(mo_renamed_last_run = NULL)
  
  if (NCOL(x) == 2) {
    # support tidyverse selection like: df %>% select(colA, colB)
    # paste these columns together
    x_vector <- vector("character", NROW(x))
    for (i in seq_len(NROW(x))) {
      x_vector[i] <- paste(pull(x[i, ], 1), pull(x[i, ], 2), sep = " ")
    }
    x <- x_vector
  } else {
    if (NCOL(x) > 2) {
      stop("`x` can be 2 columns at most", call. = FALSE)
    }
    x[is.null(x)] <- NA
    
    # support tidyverse selection like: df %>% select(colA)
    if (!is.vector(x) & !is.null(dim(x))) {
      x <- pull(x, 1)
    }
  }
  
  uncertainties <- data.frame(uncertainty = integer(0),
                              input = character(0),
                              fullname = character(0),
                              renamed_to = character(0),
                              mo = character(0), 
                              stringsAsFactors = FALSE)
  failures <- character(0)
  uncertainty_level <- translate_allow_uncertain(allow_uncertain)
  old_mo_warning <- FALSE
  
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
  if (any(x %like_case% "^[BFP]_[A-Z]{3,7}") & !all(x %in% microorganisms$mo)) {
    leftpart <- gsub("^([BFP]_[A-Z]{3,7}).*", "\\1", x)
    if (any(leftpart %in% names(mo_codes_v0.5.0))) {
      old_mo_warning <- TRUE
      rightpart <- gsub("^[BFP]_[A-Z]{3,7}(.*)", "\\1", x)
      leftpart <- mo_codes_v0.5.0[leftpart]
      x[!is.na(leftpart)] <- paste0(leftpart[!is.na(leftpart)], rightpart[!is.na(leftpart)])
    }
    # now check if some are still old
    still_old <- x[x %in% names(mo_codes_v0.5.0)]
    if (length(still_old) > 0) {
      old_mo_warning <- TRUE
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
      return(to_class_mo(rep(NA_character_, length(x_input))))
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
    
  } else if (all(x %in% reference_data_to_use$mo)) {
    # existing mo codes when not looking for property "mo", like mo_genus("B_ESCHR_COL")
    y <- reference_data_to_use[prevalence == 1][data.table(mo = x), on = "mo", ..property][[1]]
    if (any(is.na(y))) {
      y[is.na(y)] <- reference_data_to_use[prevalence == 2][data.table(mo = x[is.na(y)]),
                                                            on = "mo",
                                                            ..property][[1]]
    }
    if (any(is.na(y))) {
      y[is.na(y)] <- reference_data_to_use[prevalence == 3][data.table(mo = x[is.na(y)]),
                                                            on = "mo",
                                                            ..property][[1]]
    }
    x <- y
    
  } else if (all(toupper(x) %in% read_mo_history(uncertainty_level,
                                                 force = force_mo_history,
                                                 disable = disable_mo_history)$x)) {
    
    # previously found code
    x <- data.frame(mo = get_mo_history(x,
                                        uncertainty_level,
                                        force = force_mo_history,
                                        disable = disable_mo_history), 
                    stringsAsFactors = FALSE) %>% 
      left_join(AMR::microorganisms, by = "mo") %>% 
      pull(property)
    
  } else if (all(tolower(x) %in% reference_data_to_use$fullname_lower)) {
    # we need special treatment for very prevalent full names, they are likely!
    # e.g. as.mo("Staphylococcus aureus")
    y <- reference_data_to_use[prevalence == 1][data.table(fullname_lower = tolower(x)), on = "fullname_lower", ..property][[1]]
    if (any(is.na(y))) {
      y[is.na(y)] <- reference_data_to_use[prevalence == 2][data.table(fullname_lower = tolower(x[is.na(y)])),
                                                            on = "fullname_lower",
                                                            ..property][[1]]
    }
    if (any(is.na(y))) {
      y[is.na(y)] <- reference_data_to_use[prevalence == 3][data.table(fullname_lower = tolower(x[is.na(y)])),
                                                            on = "fullname_lower",
                                                            ..property][[1]]
    }
    x <- y
    
  } else if (all(toupper(x) %in% AMR::microorganisms.codes$code)) {
    # commonly used MO codes
    y <- as.data.table(AMR::microorganisms.codes)[data.table(code = toupper(x)), on = "code", ]
    # save them to history
    set_mo_history(x, y$mo, 0, force = force_mo_history, disable = disable_mo_history)
    
    x <- reference_data_to_use[data.table(mo = y[["mo"]]), on = "mo", ..property][[1]]
    
  } else if (all(x %in% microorganisms.translation$mo_old)) {
    # is an old mo code, used in previous versions of this package
    old_mo_warning <- TRUE
    y <- as.data.table(microorganisms.translation)[data.table(mo_old = x), on = "mo_old", "mo_new"][[1]]
    y <- reference_data_to_use[data.table(mo = y), on = "mo", ..property][[1]]
    # don't save to history, as all items are already in microorganisms.translation
    x <- y
    
  } else if (!all(x %in% AMR::microorganisms[, property])) {
    
    strip_whitespace <- function(x, dyslexia_mode) {
      # all whitespaces (tab, new lines, etc.) should be one space
      # and spaces before and after should be omitted
      trimmed <- trimws(gsub("[\\s]+", " ", x, perl = TRUE), which = "both")
      # also, make sure the trailing and leading characters are a-z or 0-9
      # in case of non-regex
      if (dyslexia_mode == FALSE) {
        trimmed <- gsub("^[^a-zA-Z0-9)(]+", "", trimmed)
        trimmed <- gsub("[^a-zA-Z0-9)(]+$", "", trimmed)
      }
      trimmed
    }
    
    x_backup_untouched <- x
    x <- strip_whitespace(x, dyslexia_mode)
    x_backup <- x
    
    # from here on case-insensitive
    x <- tolower(x)
    
    x_backup[grepl("^(fungus|fungi)$", x)] <- "F_FUNGUS" # will otherwise become the kingdom
    
    # remove spp and species
    x <- gsub(" +(spp.?|ssp.?|sp.? |ss ?.?|subsp.?|subspecies|biovar |serovar |species)", " ", x)
    x <- gsub("(spp.?|subsp.?|subspecies|biovar|serovar|species)", "", x)
    x <- strip_whitespace(x, dyslexia_mode)
    
    x_backup_without_spp <- x
    x_species <- paste(x, "species")
    # translate to English for supported languages of mo_property
    x <- gsub("(gruppe|groep|grupo|gruppo|groupe)", "group", x)
    # no groups and complexes as ending
    x <- gsub("(complex|group)$", "", x)
    x <- gsub("((an)?aero+b)[a-z]*", "", x)
    x <- gsub("^atyp[a-z]*", "", x)
    x <- gsub("(vergroen)[a-z]*", "viridans", x)
    x <- gsub("[a-z]*diff?erent[a-z]*", "", x)
    x <- gsub("(hefe|gist|gisten|levadura|lievito|fermento|levure)[a-z]*", "yeast", x)
    x <- gsub("(schimmels?|mofo|molde|stampo|moisissure|fungi)[a-z]*", "fungus", x)
    x <- gsub("fungus[ph|f]rya", "fungiphrya", x)
    # remove non-text in case of "E. coli" except dots and spaces
    x <- trimws(gsub("[^.a-zA-Z0-9/ \\-]+", " ", x))
    # replace minus by a space
    x <- gsub("-+", " ", x)
    # replace hemolytic by haemolytic
    x <- gsub("ha?emoly", "haemoly", x)
    # place minus back in streptococci
    x <- gsub("(alpha|beta|gamma).?ha?emoly", "\\1-haemoly", x)
    # remove genus as first word
    x <- gsub("^genus ", "", x)
    # remove 'uncertain'-like texts
    x <- trimws(gsub("(uncertain|susp[ie]c[a-z]+|verdacht)", "", x))
    # allow characters that resemble others = dyslexia_mode ----
    if (dyslexia_mode == TRUE) {
      x <- tolower(x)
      x <- gsub("[iy]+", "[iy]+", x)
      x <- gsub("(c|k|q|qu|s|z|x|ks)+", "(c|k|q|qu|s|z|x|ks)+", x)
      x <- gsub("(ph|hp|f|v)+", "(ph|hp|f|v)+", x)
      x <- gsub("(th|ht|t)+", "(th|ht|t)+", x)
      x <- gsub("a+", "a+", x)
      x <- gsub("u+", "u+", x)
      # allow any ending of -um, -us, -ium, -icum, -ius, -icus, -ica and -a (needs perl for the negative backward lookup):
      x <- gsub("(u\\+\\(c\\|k\\|q\\|qu\\+\\|s\\|z\\|x\\|ks\\)\\+)(?![a-z])",
                "(u[s|m]|[iy][ck]?u[ms]|[iy]?[ck]?a)", x, perl = TRUE)
      x <- gsub("(\\[iy\\]\\+\\(c\\|k\\|q\\|qu\\+\\|s\\|z\\|x\\|ks\\)\\+a\\+)(?![a-z])",
                "(u[s|m]|[iy][ck]?u[ms]|[iy]?[ck]?a)", x, perl = TRUE)
      x <- gsub("(\\[iy\\]\\+u\\+m)(?![a-z])",
                "(u[s|m]|[iy][ck]?u[ms]|[iy]?[ck]?a)", x, perl = TRUE)
      x <- gsub("e+", "e+", x)
      x <- gsub("o+", "o+", x)
      x <- gsub("(.)\\1+", "\\1+", x)
      # allow multiplication of all other consonants
      x <- gsub("([bdgjlnrw]+)", "\\1+", x)
      # allow ending in -en or -us
      x <- gsub("e\\+n(?![a-z[])", "(e+n|u+(c|k|q|qu|s|z|x|ks)+)", x, perl = TRUE)
      # if the input is longer than 10 characters, allow any forgotten consonant between all characters, as some might just have forgotten one...
      # this will allow "Pasteurella damatis" to be correctly read as "Pasteurella dagmatis".
      consonants <- paste(letters[!letters %in% c("a", "e", "i", "o", "u")], collapse = "")
      x[nchar(x_backup_without_spp) > 10] <- gsub("[+]", paste0("+[", consonants, "]?"), x[nchar(x_backup_without_spp) > 10])
      # allow au and ou after all these regex implementations
      x <- gsub("a+[bcdfghjklmnpqrstvwxyz]?u+[bcdfghjklmnpqrstvwxyz]?", "(a+u+|o+u+)[bcdfghjklmnpqrstvwxyz]?", x, fixed = TRUE)
      x <- gsub("o+[bcdfghjklmnpqrstvwxyz]?u+[bcdfghjklmnpqrstvwxyz]?", "(a+u+|o+u+)[bcdfghjklmnpqrstvwxyz]?", x, fixed = TRUE)
    }
    x <- strip_whitespace(x, dyslexia_mode)
    # make sure to remove regex overkill (will lead to errors)
    x <- gsub("++", "+", x, fixed = TRUE)
    x <- gsub("?+", "?", x, fixed = TRUE)
    
    x_trimmed <- x
    x_trimmed_species <- paste(x_trimmed, "species")
    x_trimmed_without_group <- gsub(" gro.u.p$", "", x_trimmed)
    # remove last part from "-" or "/"
    x_trimmed_without_group <- gsub("(.*)[-/].*", "\\1", x_trimmed_without_group)
    # replace space and dot by regex sign
    x_withspaces <- gsub("[ .]+", ".* ", x)
    x <- gsub("[ .]+", ".*", x)
    # add start en stop regex
    x <- paste0("^", x, "$")
    
    x_withspaces_start_only <- paste0("^", x_withspaces)
    x_withspaces_end_only <- paste0(x_withspaces, "$")
    x_withspaces_start_end <- paste0("^", x_withspaces, "$")
    
    if (isTRUE(debug)) {
      cat(paste0(blue("x"), '                       "', x, '"\n'))
      cat(paste0(blue("x_species"), '               "', x_species, '"\n'))
      cat(paste0(blue("x_withspaces_start_only"), ' "', x_withspaces_start_only, '"\n'))
      cat(paste0(blue("x_withspaces_end_only"), '   "', x_withspaces_end_only, '"\n'))
      cat(paste0(blue("x_withspaces_start_end"), '  "', x_withspaces_start_end, '"\n'))
      cat(paste0(blue("x_backup"), '                "', x_backup, '"\n'))
      cat(paste0(blue("x_backup_without_spp"), '    "', x_backup_without_spp, '"\n'))
      cat(paste0(blue("x_trimmed"), '               "', x_trimmed, '"\n'))
      cat(paste0(blue("x_trimmed_species"), '       "', x_trimmed_species, '"\n'))
      cat(paste0(blue("x_trimmed_without_group"), ' "', x_trimmed_without_group, '"\n'))
    }
    
    progress <- progress_estimated(n = length(x), min_time = 3)

    for (i in seq_len(length(x))) {
      
      progress$tick()$print()
      
      mo_hist <- get_mo_history(x_backup[i], uncertainty_level, force = force_mo_history, disable = disable_mo_history)
      if (initial_search == TRUE & !any(is.na(mo_hist))) {
        # previously found code
        found <- data.frame(mo = mo_hist, 
                            stringsAsFactors = FALSE) %>% 
          left_join(reference_data_to_use, by = "mo") %>% 
          pull(property)
        if (length(found) > 0) {
          x[i] <- found[1L]
          next
        }
      }
      
      if (x_backup[i] %like_case% "\\(unknown [a-z]+\\)") {
        x[i] <- "UNKNOWN"
        next
      }
      
      found <- reference_data_to_use[mo == toupper(x_backup[i]), ..property][[1]]
      # is a valid MO code
      if (length(found) > 0) {
        x[i] <- found[1L]
        next
      }
      
      if (x_backup[i] %in% microorganisms.translation$mo_old) {
        # is an old mo code, used in previous versions of this package
        old_mo_warning <- TRUE
        found <- reference_data_to_use[mo == microorganisms.translation[which(microorganisms.translation$mo_old == x_backup[i]), "mo_new"], ..property][[1]]
        if (length(found) > 0) {
          x[i] <- found[1L]
          # don't save to history, as all items are already in microorganisms.translation
          next
        }
      }

      if (toupper(x_backup_untouched[i]) %in% microorganisms.codes$code) {
        # is a WHONET code, like "HA-"
        found <- microorganismsDT[mo == microorganisms.codes[which(microorganisms.codes$code == toupper(x_backup_untouched[i])), "mo"][1L], ..property][[1]]
        if (length(found) > 0) {
          x[i] <- found[1L]
          # don't save to history, as all items are already in microorganisms.codes
          next
        }
      }
      
      found <- reference_data_to_use[fullname_lower %in% tolower(c(x_backup[i], x_backup_without_spp[i])), ..property][[1]]
      # most probable: is exact match in fullname
      if (length(found) > 0) {
        x[i] <- found[1L]
        if (initial_search == TRUE) {
          # don't save valid fullnames to history (i.e. values that are in microorganisms$fullname)
        }
        next
      }
      
      found <- reference_data_to_use[col_id == x_backup[i], ..property][[1]]
      # is a valid Catalogue of Life ID
      if (NROW(found) > 0) {
        x[i] <- found[1L]
        if (initial_search == TRUE) {
          set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history, disable = disable_mo_history)
        }
        next
      }
      
      # WHONET and other common LIS codes
      if (any(toupper(c(x_backup[i], x_backup_without_spp[i])) %in% AMR::microorganisms.codes$code)) {
        mo_found <- AMR::microorganisms.codes[which(AMR::microorganisms.codes$code %in% toupper(c(x_backup[i], x_backup_without_spp[i]))), "mo"][1L]
        if (length(mo_found) > 0) {
          x[i] <- microorganismsDT[mo == mo_found, ..property][[1]][1L]
          if (initial_search == TRUE) {
            set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history, disable = disable_mo_history)
          }
          next
        }
      }
      if (!is.null(reference_df)) {
        # self-defined reference
        if (x_backup[i] %in% reference_df[, 1]) {
          ref_mo <- reference_df[reference_df[, 1] == x_backup[i], "mo"][[1L]]
          if (ref_mo %in% microorganismsDT[, mo]) {
            x[i] <- microorganismsDT[mo == ref_mo, ..property][[1]][1L]
            next
          } else {
            warning("Value '", x_backup[i], "' was found in reference_df, but '", ref_mo, "' is not a valid MO code.", call. = FALSE)
          }
        }
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
          set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history, disable = disable_mo_history)
        }
        next
      }
      
      # check for very small input, but ignore the O antigens of E. coli
      if (nchar(gsub("[^a-zA-Z]", "", x_trimmed[i])) < 3
          & !x_backup_without_spp[i] %like_case% "[Oo]?(26|103|104|104|111|121|145|157)") {
        # fewer than 3 chars and not looked for species, add as failure
        x[i] <- microorganismsDT[mo == "UNKNOWN", ..property][[1]]
        if (initial_search == TRUE) {
          failures <- c(failures, x_backup[i])
          set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history, disable = disable_mo_history)
        }
        next
      }
      
      if (x_backup_without_spp[i] %like_case% "virus") {
        # there is no fullname like virus, so don't try to coerce it
        x[i] <- NA_character_
        next
      }
      
      # translate known trivial abbreviations to genus + species ----
      if (!is.na(x_trimmed[i])) {
        if (toupper(x_backup_without_spp[i]) %in% c("MRSA", "MSSA", "VISA", "VRSA")
            | x_backup_without_spp[i] %like_case% " (mrsa|mssa|visa|vrsa) ") {
          x[i] <- microorganismsDT[mo == "B_STPHY_AURS", ..property][[1]][1L]
          if (initial_search == TRUE) {
            set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history, disable = disable_mo_history)
          }
          next
        }
        if (toupper(x_backup_without_spp[i]) %in% c("MRSE", "MSSE")
            | x_backup_without_spp[i] %like_case% " (mrse|msse) ") {
          x[i] <- microorganismsDT[mo == "B_STPHY_EPDR", ..property][[1]][1L]
          if (initial_search == TRUE) {
            set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history, disable = disable_mo_history)
          }
          next
        }
        if (toupper(x_backup_without_spp[i]) == "VRE"
            | x_backup_without_spp[i] %like_case% " vre "
            | x_backup_without_spp[i] %like_case% "(enterococci|enterokok|enterococo)[a-z]*?$")  {
          x[i] <- microorganismsDT[mo == "B_ENTRC", ..property][[1]][1L]
          if (initial_search == TRUE) {
            set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history, disable = disable_mo_history)
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
            # also support O-antigens of E. coli: O26, O103, O104, O111, O121, O145, O157
            | x_backup_without_spp[i] %like_case% "o?(26|103|104|111|121|145|157)") {
          x[i] <- microorganismsDT[mo == "B_ESCHR_COLI", ..property][[1]][1L]
          if (initial_search == TRUE) {
            set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history, disable = disable_mo_history)
          }
          next
        }
        if (toupper(x_backup_without_spp[i]) == "MRPA"
            | x_backup_without_spp[i] %like_case% " mrpa ") {
          # multi resistant P. aeruginosa
          x[i] <- microorganismsDT[mo == "B_PSDMN_ARGN", ..property][[1]][1L]
          if (initial_search == TRUE) {
            set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history, disable = disable_mo_history)
          }
          next
        }
        if (toupper(x_backup_without_spp[i]) == "CRSM") {
          # co-trim resistant S. maltophilia
          x[i] <- microorganismsDT[mo == "B_STNTR_MLTP", ..property][[1]][1L]
          if (initial_search == TRUE) {
            set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history, disable = disable_mo_history)
          }
          next
        }
        if (toupper(x_backup_without_spp[i]) %in% c("PISP", "PRSP", "VISP", "VRSP")
            | x_backup_without_spp[i] %like_case% " (pisp|prsp|visp|vrsp) ") {
          # peni I, peni R, vanco I, vanco R: S. pneumoniae
          x[i] <- microorganismsDT[mo == "B_STRPT_PNMN", ..property][[1]][1L]
          if (initial_search == TRUE) {
            set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history, disable = disable_mo_history)
          }
          next
        }
        if (x_backup_without_spp[i] %like_case% "^g[abcdfghk]s$") {
          # Streptococci, like GBS = Group B Streptococci (B_STRPT_GRPB)
          x[i] <- microorganismsDT[mo == toupper(gsub("g([abcdfghk])s", "B_STRPT_GRP\\1", x_backup_without_spp[i])), ..property][[1]][1L]
          if (initial_search == TRUE) {
            set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history, disable = disable_mo_history)
          }
          next
        }
        if (x_backup_without_spp[i] %like_case% "(streptococ|streptokok).* [abcdfghk]$") {
          # Streptococci in different languages, like "estreptococos grupo B"
          x[i] <- microorganismsDT[mo == toupper(gsub(".*(streptococ|streptokok|estreptococ).* ([abcdfghk])$", "B_STRPT_GRP\\2", x_backup_without_spp[i])), ..property][[1]][1L]
          if (initial_search == TRUE) {
            set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history, disable = disable_mo_history)
          }
          next
        }
        if (x_backup_without_spp[i] %like_case% "group [abcdfghk] (streptococ|streptokok|estreptococ)") {
          # Streptococci in different languages, like "Group A Streptococci"
          x[i] <- microorganismsDT[mo == toupper(gsub(".*group ([abcdfghk]) (streptococ|streptokok|estreptococ).*", "B_STRPT_GRP\\1", x_backup_without_spp[i])), ..property][[1]][1L]
          if (initial_search == TRUE) {
            set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history, disable = disable_mo_history)
          }
          next
        }
        if (x_backup_without_spp[i] %like_case% "haemoly.*strept") {
          # Haemolytic streptococci in different languages
          x[i] <- microorganismsDT[mo == "B_STRPT_HAEM", ..property][[1]][1L]
          if (initial_search == TRUE) {
            set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history, disable = disable_mo_history)
          }
          next
        }
        # CoNS/CoPS in different languages (support for German, Dutch, Spanish, Portuguese) ----
        if (x_backup_without_spp[i] %like_case% "[ck]oagulas[ea] negatie?[vf]"
            | x_trimmed[i] %like_case% "[ck]oagulas[ea] negatie?[vf]"
            | x_backup_without_spp[i] %like_case% "[ck]o?ns[^a-z]?$") {
          # coerce S. coagulase negative
          x[i] <- microorganismsDT[mo == "B_STPHY_CONS", ..property][[1]][1L]
          if (initial_search == TRUE) {
            set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history, disable = disable_mo_history)
          }
          next
        }
        if (x_backup_without_spp[i] %like_case% "[ck]oagulas[ea] positie?[vf]"
            | x_trimmed[i] %like_case% "[ck]oagulas[ea] positie?[vf]"
            | x_backup_without_spp[i] %like_case% "[ck]o?ps[^a-z]?$") {
          # coerce S. coagulase positive
          x[i] <- microorganismsDT[mo == "B_STPHY_COPS", ..property][[1]][1L]
          if (initial_search == TRUE) {
            set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history, disable = disable_mo_history)
          }
          next
        }
        # streptococcal groups: milleri and viridans
        if (x_trimmed[i] %like_case% "strepto.* mil+er+i"
            | x_backup_without_spp[i] %like_case% "strepto.* mil+er+i"
            | x_backup_without_spp[i] %like_case% "mgs[^a-z]?$") {
          # Milleri Group Streptococcus (MGS)
          x[i] <- microorganismsDT[mo == "B_STRPT_MILL", ..property][[1]][1L]
          if (initial_search == TRUE) {
            set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history, disable = disable_mo_history)
          }
          next
        }
        if (x_trimmed[i] %like_case% "strepto.* viridans"
            | x_backup_without_spp[i] %like_case% "strepto.* viridans"
            | x_backup_without_spp[i] %like_case% "vgs[^a-z]?$") {
          # Viridans Group Streptococcus (VGS)
          x[i] <- microorganismsDT[mo == "B_STRPT_VIRI", ..property][[1]][1L]
          if (initial_search == TRUE) {
            set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history, disable = disable_mo_history)
          }
          next
        }
        if (x_backup_without_spp[i] %like_case% "gram[ -]?neg.*"
            | x_backup_without_spp[i] %like_case% "negatie?[vf]"
            | x_trimmed[i] %like_case% "gram[ -]?neg.*") {
          # coerce Gram negatives
          x[i] <- microorganismsDT[mo == "B_GRAMN", ..property][[1]][1L]
          if (initial_search == TRUE) {
            set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history, disable = disable_mo_history)
          }
          next
        }
        if (x_backup_without_spp[i] %like_case% "gram[ -]?pos.*"
            | x_backup_without_spp[i] %like_case% "positie?[vf]"
            | x_trimmed[i] %like_case% "gram[ -]?pos.*") {
          # coerce Gram positives
          x[i] <- microorganismsDT[mo == "B_GRAMP", ..property][[1]][1L]
          if (initial_search == TRUE) {
            set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history, disable = disable_mo_history)
          }
          next
        }
        if (x_backup_without_spp[i] %like_case% "mycoba[ck]teri.[nm]?$") {
          # coerce Gram positives
          x[i] <- microorganismsDT[mo == "B_MYCBC", ..property][[1]][1L]
          if (initial_search == TRUE) {
            set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history, disable = disable_mo_history)
          }
          next
        }
        
        if (x_backup_without_spp[i] %like_case% "salmonella [a-z]+ ?.*") {
          if (x_backup_without_spp[i] %like_case% "salmonella group") {
            # Salmonella Group A to Z, just return S. species for now
            x[i] <- microorganismsDT[mo == "B_SLMNL", ..property][[1]][1L]
            if (initial_search == TRUE) {
              set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history, disable = disable_mo_history)
            }
            next
          } else if (grepl("[sS]almonella [A-Z][a-z]+ ?.*", x_backup[i], ignore.case = FALSE)) {
            # Salmonella with capital letter species like "Salmonella Goettingen" - they're all S. enterica
            x[i] <- microorganismsDT[mo == "B_SLMNL_ENTR", ..property][[1]][1L]
            if (initial_search == TRUE) {
              set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history, disable = disable_mo_history)
            }
            uncertainties <- rbind(uncertainties,
                                   format_uncertainty_as_df(uncertainty_level = 1,
                                                            input = x_backup[i],
                                                            result_mo = "B_SLMNL_ENTR"))
            next
          }
        }
        
        # trivial names known to the field:
        if ("meningococcus" %like_case% x_trimmed[i]) {
          # coerce Neisseria meningitidis
          x[i] <- microorganismsDT[mo == "B_NESSR_MNNG", ..property][[1]][1L]
          if (initial_search == TRUE) {
            set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history, disable = disable_mo_history)
          }
          next
        }
        if ("gonococcus" %like_case% x_trimmed[i]) {
          # coerce Neisseria gonorrhoeae
          x[i] <- microorganismsDT[mo == "B_NESSR_GNRR", ..property][[1]][1L]
          if (initial_search == TRUE) {
            set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history, disable = disable_mo_history)
          }
          next
        }
        if ("pneumococcus" %like_case% x_trimmed[i]) {
          # coerce Streptococcus penumoniae
          x[i] <- microorganismsDT[mo == "B_STRPT_PNMN", ..property][[1]][1L]
          if (initial_search == TRUE) {
            set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history, disable = disable_mo_history)
          }
          next
        }
      }
      
      # NOW RUN THROUGH DIFFERENT PREVALENCE LEVELS
      check_per_prevalence <- function(data_to_check,
                                       data.old_to_check,
                                       a.x_backup,
                                       b.x_trimmed,
                                       c.x_trimmed_without_group,
                                       d.x_withspaces_start_end,
                                       e.x_withspaces_start_only,
                                       f.x_withspaces_end_only,
                                       g.x_backup_without_spp,
                                       h.x_species,
                                       i.x_trimmed_species) {
        
        # FIRST TRY FULLNAMES AND CODES ----
        # if only genus is available, return only genus
        
        if (all(!c(x[i], b.x_trimmed) %like_case% " ")) {
          found <- data_to_check[fullname_lower %in% c(h.x_species, i.x_trimmed_species), ..property][[1]]
          if (length(found) > 0) {
            x[i] <- found[1L]
            if (initial_search == TRUE) {
              set_mo_history(a.x_backup, get_mo_code(x[i], property), 0, force = force_mo_history, disable = disable_mo_history)
            }
            return(x[i])
          }
          if (nchar(g.x_backup_without_spp) >= 6) {
            found <- data_to_check[fullname_lower %like_case% paste0("^", unregex(g.x_backup_without_spp), "[a-z]+"), ..property][[1]]
            if (length(found) > 0) {
              x[i] <- found[1L]
              if (initial_search == TRUE) {
                set_mo_history(a.x_backup, get_mo_code(x[i], property), 0, force = force_mo_history, disable = disable_mo_history)
              }
              return(x[i])
            }
          }
          # rest of genus only is in allow_uncertain part.
        }
        
        # allow no codes less than 4 characters long, was already checked for WHONET earlier
        if (nchar(g.x_backup_without_spp) < 4) {
          x[i] <- microorganismsDT[mo == "UNKNOWN", ..property][[1]]
          if (initial_search == TRUE) {
            failures <- c(failures, a.x_backup)
            set_mo_history(a.x_backup, get_mo_code(x[i], property), 0, force = force_mo_history, disable = disable_mo_history)
          }
          return(x[i])
        }
        
        # try probable: trimmed version of fullname ----
        found <- data_to_check[fullname_lower %in% tolower(g.x_backup_without_spp), ..property][[1]]
        if (length(found) > 0) {
          return(found[1L])
        }
        
        # try any match keeping spaces ----
        found <- data_to_check[fullname_lower %like_case% d.x_withspaces_start_end, ..property][[1]]
        if (length(found) > 0 & nchar(g.x_backup_without_spp) >= 6) {
          return(found[1L])
        }
        
        # try any match keeping spaces, not ending with $ ----
        found <- data_to_check[fullname_lower %like_case% paste0(trimws(e.x_withspaces_start_only), " "), ..property][[1]]
        if (length(found) > 0) {
          return(found[1L])
        }
        found <- data_to_check[fullname_lower %like_case% e.x_withspaces_start_only, ..property][[1]]
        if (length(found) > 0 & nchar(g.x_backup_without_spp) >= 6) {
          return(found[1L])
        }
        
        # try any match keeping spaces, not start with ^ ----
        found <- data_to_check[fullname_lower %like_case% paste0(" ", trimws(f.x_withspaces_end_only)), ..property][[1]]
        if (length(found) > 0) {
          return(found[1L])
        }
        
        # try a trimmed version
        found <- data_to_check[fullname_lower %like_case% b.x_trimmed
                               | fullname_lower %like_case% c.x_trimmed_without_group, ..property][[1]]
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
                            ".* ",
                            g.x_backup_without_spp %>% substr((x_length / 2) + 1, x_length))
          found <- data_to_check[fullname_lower %like_case% x_split, ..property][[1]]
          if (length(found) > 0) {
            return(found[1L])
          }
        }
        
        # try fullname without start and without nchar limit of >= 6 ----
        # like "K. pneu rhino" >> "Klebsiella pneumoniae (rhinoscleromatis)" = KLEPNERH
        found <- data_to_check[fullname_lower %like_case% e.x_withspaces_start_only, ..property][[1]]
        if (length(found) > 0) {
          return(found[1L])
        }
        
        # MISCELLANEOUS ----
        
        # look for old taxonomic names ----
        # wait until prevalence == 2 to run the old taxonomic results on both prevalence == 1 and prevalence == 2
        found <- data.old_to_check[fullname_lower == tolower(a.x_backup)
                                   | fullname_lower %like_case% d.x_withspaces_start_end, ]
        if (NROW(found) > 0) {
          col_id_new <- found[1, col_id_new]
          # when property is "ref" (which is the case in mo_ref, mo_authors and mo_year), return the old value, so:
          # mo_ref() of "Chlamydia psittaci" will be "Page, 1968" (with warning)
          # mo_ref() of "Chlamydophila psittaci" will be "Everett et al., 1999"
          if (property == "ref") {
            x[i] <- found[1, ref]
          } else {
            x[i] <- microorganismsDT[col_id == found[1, col_id_new], ..property][[1]]
          }
          options(mo_renamed_last_run = found[1, fullname])
          was_renamed(name_old = found[1, fullname],
                      name_new = microorganismsDT[col_id == found[1, col_id_new], fullname],
                      ref_old = found[1, ref],
                      ref_new = microorganismsDT[col_id == found[1, col_id_new], ref],
                      mo = microorganismsDT[col_id == found[1, col_id_new], mo])
          # no set history on renames
          return(x[i])
        }
        
        # check for uncertain results ----
        uncertain_fn <- function(a.x_backup,
                                 b.x_trimmed,
                                 d.x_withspaces_start_end,
                                 e.x_withspaces_start_only,
                                 f.x_withspaces_end_only,
                                 g.x_backup_without_spp,
                                 uncertain.reference_data_to_use) {
          
          if (uncertainty_level == 0) {
            # do not allow uncertainties
            return(NA_character_)
          }
          
          # UNCERTAINTY LEVEL 1 ----
          if (uncertainty_level >= 1) {
            now_checks_for_uncertainty_level <- 1
            
            # (1) look again for old taxonomic names, now for G. species ----
            if (isTRUE(debug)) {
              cat("\n[ UNCERTAINTY LEVEL", now_checks_for_uncertainty_level, "] (1) look again for old taxonomic names, now for G. species\n")
            }
            if (isTRUE(debug)) {
              message("Running '", d.x_withspaces_start_end, "' and '", e.x_withspaces_start_only, "'")
            }
            found <- data.old_to_check[fullname_lower %like_case% d.x_withspaces_start_end
                                       | fullname_lower %like_case% e.x_withspaces_start_only]
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
              options(mo_renamed_last_run = found[1, fullname])
              uncertainties <<- rbind(uncertainties,
                                      format_uncertainty_as_df(uncertainty_level = now_checks_for_uncertainty_level,
                                                               input = a.x_backup,
                                                               result_mo = microorganismsDT[col_id == found[1, col_id_new], mo]))
              # no set history on renames
              return(x)
            }
            
            # (2) Try with misspelled input ----
            # just rerun with dyslexia_mode = TRUE will used the extensive regex part above
            if (isTRUE(debug)) {
              cat("\n[ UNCERTAINTY LEVEL", now_checks_for_uncertainty_level, "] (2) Try with misspelled input\n")
            }
            if (isTRUE(debug)) {
              message("Running '", a.x_backup, "'")
            }
            # first try without dyslexia mode
            found <- suppressMessages(suppressWarnings(exec_as.mo(a.x_backup, initial_search = FALSE, dyslexia_mode = FALSE, allow_uncertain = FALSE, debug = debug, reference_data_to_use = uncertain.reference_data_to_use)))
            if (empty_result(found)) {
              # then with dyslexia mode
              found <- suppressMessages(suppressWarnings(exec_as.mo(a.x_backup, initial_search = FALSE, dyslexia_mode = TRUE, allow_uncertain = FALSE, debug = debug, reference_data_to_use = uncertain.reference_data_to_use)))
            }
            if (!empty_result(found)) {
              found_result <- found
              found <- reference_data_to_use[mo == found, ..property][[1]]
              uncertainties <<- rbind(uncertainties,
                                      format_uncertainty_as_df(uncertainty_level = now_checks_for_uncertainty_level,
                                                               input = a.x_backup,
                                                               result_mo = found_result[1L]))
              if (initial_search == TRUE) {
                set_mo_history(a.x_backup, get_mo_code(found[1L], property), 1, force = force_mo_history, disable = disable_mo_history)
              }
              return(found[1L])
            }
          }
          
          # UNCERTAINTY LEVEL 2 ----
          if (uncertainty_level >= 2) {
            now_checks_for_uncertainty_level <- 2
            
            # (3) look for genus only, part of name ----
            if (isTRUE(debug)) {
              cat("\n[ UNCERTAINTY LEVEL", now_checks_for_uncertainty_level, "] (3) look for genus only, part of name\n")
            }
            if (nchar(g.x_backup_without_spp) > 4 & !b.x_trimmed %like_case% " ") {
              if (!grepl("^[A-Z][a-z]+", b.x_trimmed, ignore.case = FALSE)) {
                if (isTRUE(debug)) {
                  message("Running '", paste(b.x_trimmed, "species"), "'")
                }
                # not when input is like Genustext, because then Neospora would lead to Actinokineospora
                found <- uncertain.reference_data_to_use[fullname_lower %like_case% paste(b.x_trimmed, "species"), ..property][[1]]
                if (length(found) > 0) {
                  x[i] <- found[1L]
                  uncertainties <<- rbind(uncertainties,
                                          format_uncertainty_as_df(uncertainty_level = now_checks_for_uncertainty_level,
                                                                   input = a.x_backup,
                                                                   result_mo = found_result[1L]))
                  if (initial_search == TRUE) {
                    set_mo_history(a.x_backup, get_mo_code(x, property), 2, force = force_mo_history, disable = disable_mo_history)
                  }
                  return(x)
                }
              }
            }
            
            # (4) strip values between brackets ----
            if (isTRUE(debug)) {
              cat("\n[ UNCERTAINTY LEVEL", now_checks_for_uncertainty_level, "] (4) strip values between brackets\n")
            }
            a.x_backup_stripped <- gsub("( *[(].*[)] *)", " ", a.x_backup)
            a.x_backup_stripped <- trimws(gsub(" +", " ", a.x_backup_stripped))
            if (isTRUE(debug)) {
              message("Running '", a.x_backup_stripped, "'")
            }
            # first try without dyslexia mode
            found <- suppressMessages(suppressWarnings(exec_as.mo(a.x_backup_stripped, initial_search = FALSE, dyslexia_mode = FALSE, allow_uncertain = FALSE, debug = debug, reference_data_to_use = uncertain.reference_data_to_use)))
            if (empty_result(found)) {
              # then with dyslexia mode
              found <- suppressMessages(suppressWarnings(exec_as.mo(a.x_backup_stripped, initial_search = FALSE, dyslexia_mode = TRUE, allow_uncertain = FALSE, debug = debug, reference_data_to_use = uncertain.reference_data_to_use)))
            }
            if (!empty_result(found) & nchar(g.x_backup_without_spp) >= 6) {
              found_result <- found
              found <- reference_data_to_use[mo == found, ..property][[1]]
              uncertainties <<- rbind(uncertainties,
                                      format_uncertainty_as_df(uncertainty_level = now_checks_for_uncertainty_level,
                                                               input = a.x_backup,
                                                               result_mo = found_result[1L]))
              if (initial_search == TRUE) {
                set_mo_history(a.x_backup, get_mo_code(found[1L], property), 2, force = force_mo_history, disable = disable_mo_history)
              }
              return(found[1L])
            }
            
            # (5) inverse input ----
            if (isTRUE(debug)) {
              cat("\n[ UNCERTAINTY LEVEL", now_checks_for_uncertainty_level, "] (5) inverse input\n")
            }
            a.x_backup_inversed <- paste(rev(unlist(strsplit(a.x_backup, split = " "))), collapse = " ")
            if (isTRUE(debug)) {
              message("Running '", a.x_backup_inversed, "'")
            }
            # first try without dyslexia mode
            found <- suppressMessages(suppressWarnings(exec_as.mo(a.x_backup_inversed, initial_search = FALSE, dyslexia_mode = FALSE, allow_uncertain = FALSE, debug = debug, reference_data_to_use = uncertain.reference_data_to_use)))
            if (empty_result(found)) {
              # then with dyslexia mode
              found <- suppressMessages(suppressWarnings(exec_as.mo(a.x_backup_inversed, initial_search = FALSE, dyslexia_mode = TRUE, allow_uncertain = FALSE, debug = debug, reference_data_to_use = uncertain.reference_data_to_use)))
            }
            if (!empty_result(found) & nchar(g.x_backup_without_spp) >= 6) {
              found_result <- found
              found <- reference_data_to_use[mo == found, ..property][[1]]
              uncertainties <<- rbind(uncertainties,
                                      format_uncertainty_as_df(uncertainty_level = now_checks_for_uncertainty_level,
                                                               input = a.x_backup,
                                                               result_mo = found_result[1L]))
              if (initial_search == TRUE) {
                set_mo_history(a.x_backup, get_mo_code(found[1L], property), 2, force = force_mo_history, disable = disable_mo_history)
              }
              return(found[1L])
            }
            
            # (6) try to strip off half an element from end and check the remains ----
            if (isTRUE(debug)) {
              cat("\n[ UNCERTAINTY LEVEL", now_checks_for_uncertainty_level, "] (6) try to strip off half an element from end and check the remains\n")
            }
            x_strip <- a.x_backup %>% strsplit(" ") %>% unlist()
            if (length(x_strip) > 1) {
              for (i in seq_len(length(x_strip) - 1)) {
                lastword <- x_strip[length(x_strip) - i + 1]
                lastword_half <- substr(lastword, 1, as.integer(nchar(lastword) / 2))
                # remove last half of the second term
                x_strip_collapsed <- paste(c(x_strip[seq_len(length(x_strip) - i)], lastword_half), collapse = " ")
                if (nchar(x_strip_collapsed) >= 4 & nchar(lastword_half) > 2) {
                  if (isTRUE(debug)) {
                    message("Running '", x_strip_collapsed, "'")
                  }
                  # first try without dyslexia mode
                  found <- suppressMessages(suppressWarnings(exec_as.mo(x_strip_collapsed, initial_search = FALSE, dyslexia_mode = FALSE, allow_uncertain = FALSE, debug = debug, reference_data_to_use = uncertain.reference_data_to_use)))
                  if (empty_result(found)) {
                    # then with dyslexia mode
                    found <- suppressMessages(suppressWarnings(exec_as.mo(x_strip_collapsed, initial_search = FALSE, dyslexia_mode = TRUE, allow_uncertain = FALSE, debug = debug, reference_data_to_use = uncertain.reference_data_to_use)))
                  }
                  if (!empty_result(found)) {
                    found_result <- found
                    found <- reference_data_to_use[mo == found, ..property][[1]]
                    uncertainties <<- rbind(uncertainties,
                                            format_uncertainty_as_df(uncertainty_level = now_checks_for_uncertainty_level,
                                                                     input = a.x_backup,
                                                                     result_mo = found_result[1L]))
                    if (initial_search == TRUE) {
                      set_mo_history(a.x_backup, get_mo_code(found[1L], property), 2, force = force_mo_history, disable = disable_mo_history)
                    }
                    return(found[1L])
                  }
                }
              }
            }
            # (7) try to strip off one element from end and check the remains ----
            if (isTRUE(debug)) {
              cat("\n[ UNCERTAINTY LEVEL", now_checks_for_uncertainty_level, "] (7) try to strip off one element from end and check the remains\n")
            }
            if (length(x_strip) > 1) {
              for (i in seq_len(length(x_strip) - 1)) {
                x_strip_collapsed <- paste(x_strip[seq_len(length(x_strip) - i)], collapse = " ")
                if (nchar(x_strip_collapsed) >= 6) {
                  if (isTRUE(debug)) {
                    message("Running '", x_strip_collapsed, "'")
                  }
                  # first try without dyslexia mode
                  found <- suppressMessages(suppressWarnings(exec_as.mo(x_strip_collapsed, initial_search = FALSE, dyslexia_mode = FALSE, allow_uncertain = FALSE, debug = debug, reference_data_to_use = uncertain.reference_data_to_use)))
                  if (empty_result(found)) {
                    # then with dyslexia mode
                    found <- suppressMessages(suppressWarnings(exec_as.mo(x_strip_collapsed, initial_search = FALSE, dyslexia_mode = TRUE, allow_uncertain = FALSE, debug = debug, reference_data_to_use = uncertain.reference_data_to_use)))
                  }
                  if (!empty_result(found)) {
                    found_result <- found
                    found <- reference_data_to_use[mo == found, ..property][[1]]
                    uncertainties <<- rbind(uncertainties,
                                            format_uncertainty_as_df(uncertainty_level = now_checks_for_uncertainty_level,
                                                                     input = a.x_backup,
                                                                     result_mo = found_result[1L]))
                    if (initial_search == TRUE) {
                      set_mo_history(a.x_backup, get_mo_code(found[1L], property), 2, force = force_mo_history, disable = disable_mo_history)
                    }
                    return(found[1L])
                  }
                }
              }
            }
            # (8) check for unknown yeasts/fungi ----
            if (isTRUE(debug)) {
              cat("\n[ UNCERTAINTY LEVEL", now_checks_for_uncertainty_level, "] (8) check for unknown yeasts/fungi\n")
            }
            if (b.x_trimmed %like_case% "yeast") {
              found <- "F_YEAST"
              found_result <- found
              found <- microorganismsDT[mo == found, ..property][[1]]
              uncertainties <<- rbind(uncertainties,
                                      format_uncertainty_as_df(uncertainty_level = now_checks_for_uncertainty_level,
                                                               input = a.x_backup,
                                                               result_mo = found_result[1L]))
              if (initial_search == TRUE) {
                set_mo_history(a.x_backup, get_mo_code(found[1L], property), 2, force = force_mo_history, disable = disable_mo_history)
              }
              return(found[1L])
            }
            if (b.x_trimmed %like_case% "(fungus|fungi)" & !b.x_trimmed %like_case% "fungiphrya") {
              found <- "F_FUNGUS"
              found_result <- found
              found <- microorganismsDT[mo == found, ..property][[1]]
              uncertainties <<- rbind(uncertainties,
                                      format_uncertainty_as_df(uncertainty_level = now_checks_for_uncertainty_level,
                                                               input = a.x_backup,
                                                               result_mo = found_result[1L]))
              if (initial_search == TRUE) {
                set_mo_history(a.x_backup, get_mo_code(found[1L], property), 2, force = force_mo_history, disable = disable_mo_history)
              }
              return(found[1L])
            }
            # (9) try to strip off one element from start and check the remains (only allow >= 2-part name outcome) ----
            if (isTRUE(debug)) {
              cat("\n[ UNCERTAINTY LEVEL", now_checks_for_uncertainty_level, "] (9) try to strip off one element from start and check the remains (only allow >= 2-part name outcome)\n")
            }
            x_strip <- a.x_backup %>% strsplit(" ") %>% unlist()
            if (length(x_strip) > 1 & nchar(g.x_backup_without_spp) >= 6) {
              for (i in 2:(length(x_strip))) {
                x_strip_collapsed <- paste(x_strip[i:length(x_strip)], collapse = " ")
                if (isTRUE(debug)) {
                  message("Running '", x_strip_collapsed, "'")
                }
                # first try without dyslexia mode
                found <- suppressMessages(suppressWarnings(exec_as.mo(x_strip_collapsed, initial_search = FALSE, dyslexia_mode = FALSE, allow_uncertain = FALSE, debug = debug, reference_data_to_use = uncertain.reference_data_to_use)))
                if (empty_result(found)) {
                  # then with dyslexia mode
                  found <- suppressMessages(suppressWarnings(exec_as.mo(x_strip_collapsed, initial_search = FALSE, dyslexia_mode = TRUE, allow_uncertain = FALSE, debug = debug, reference_data_to_use = uncertain.reference_data_to_use)))
                }
                if (!empty_result(found)) {
                  found_result <- found
                  found <- reference_data_to_use[mo == found_result[1L], ..property][[1]]
                  # uncertainty level 2 only if searched part contains a space (otherwise it will be found with lvl 3)
                  if (x_strip_collapsed %like_case% " ") {
                    uncertainties <<- rbind(uncertainties,
                                            format_uncertainty_as_df(uncertainty_level = now_checks_for_uncertainty_level,
                                                                     input = a.x_backup,
                                                                     result_mo = found_result[1L]))
                    if (initial_search == TRUE) {
                      set_mo_history(a.x_backup, get_mo_code(found[1L], property), 2, force = force_mo_history, disable = disable_mo_history)
                    }
                    return(found[1L])
                  }
                }
              }
            }
          }
          
          # UNCERTAINTY LEVEL 3 ----
          if (uncertainty_level >= 3) {
            now_checks_for_uncertainty_level <- 3
            
            # (10) try to strip off one element from start and check the remains (any text size) ----
            if (isTRUE(debug)) {
              cat("\n[ UNCERTAINTY LEVEL", now_checks_for_uncertainty_level, "] (10) try to strip off one element from start and check the remains (any text size)\n")
            }
            x_strip <- a.x_backup %>% strsplit(" ") %>% unlist()
            if (length(x_strip) > 1 & nchar(g.x_backup_without_spp) >= 6) {
              for (i in 2:(length(x_strip))) {
                x_strip_collapsed <- paste(x_strip[i:length(x_strip)], collapse = " ")
                if (isTRUE(debug)) {
                  message("Running '", x_strip_collapsed, "'")
                }
                # first try without dyslexia mode
                found <- suppressMessages(suppressWarnings(exec_as.mo(x_strip_collapsed, initial_search = FALSE, dyslexia_mode = FALSE, allow_uncertain = FALSE, debug = debug, reference_data_to_use = uncertain.reference_data_to_use)))
                if (empty_result(found)) {
                  # then with dyslexia mode
                  found <- suppressMessages(suppressWarnings(exec_as.mo(x_strip_collapsed, initial_search = FALSE, dyslexia_mode = TRUE, allow_uncertain = FALSE, debug = debug, reference_data_to_use = uncertain.reference_data_to_use)))
                }
                if (!empty_result(found)) {
                  found_result <- found
                  found <- reference_data_to_use[mo == found, ..property][[1]]
                  uncertainties <<- rbind(uncertainties,
                                          format_uncertainty_as_df(uncertainty_level = now_checks_for_uncertainty_level,
                                                                   input = a.x_backup,
                                                                   result_mo = found_result[1L]))
                  if (initial_search == TRUE) {
                    set_mo_history(a.x_backup, get_mo_code(found[1L], property), 3, force = force_mo_history, disable = disable_mo_history)
                  }
                  return(found[1L])
                }
              }
            }
            # (11) try to strip off one element from end and check the remains (any text size) ----
            # (this is in fact 7 but without nchar limit of >=6)
            if (isTRUE(debug)) {
              cat("\n[ UNCERTAINTY LEVEL", now_checks_for_uncertainty_level, "] (11) try to strip off one element from end and check the remains (any text size)\n")
            }
            if (length(x_strip) > 1) {
              for (i in seq_len(length(x_strip) - 1)) {
                x_strip_collapsed <- paste(x_strip[seq_len(length(x_strip) - i)], collapse = " ")
                if (isTRUE(debug)) {
                  message("Running '", x_strip_collapsed, "'")
                }
                # first try without dyslexia mode
                found <- suppressMessages(suppressWarnings(exec_as.mo(x_strip_collapsed, initial_search = FALSE, dyslexia_mode = FALSE, allow_uncertain = FALSE, debug = debug, reference_data_to_use = uncertain.reference_data_to_use)))
                if (empty_result(found)) {
                  # then with dyslexia mode
                  found <- suppressMessages(suppressWarnings(exec_as.mo(x_strip_collapsed, initial_search = FALSE, dyslexia_mode = TRUE, allow_uncertain = FALSE, debug = debug, reference_data_to_use = uncertain.reference_data_to_use)))
                }
                if (!empty_result(found)) {
                  found_result <- found
                  found <- reference_data_to_use[mo == found, ..property][[1]]
                  uncertainties <<- rbind(uncertainties,
                                          format_uncertainty_as_df(uncertainty_level = now_checks_for_uncertainty_level,
                                                                   input = a.x_backup,
                                                                   result_mo = found_result[1L]))
                  if (initial_search == TRUE) {
                    set_mo_history(a.x_backup, get_mo_code(found[1L], property), 2, force = force_mo_history, disable = disable_mo_history)
                  }
                  return(found[1L])
                }
              }
            }
            
            # (12) part of a name (very unlikely match) ----
            if (isTRUE(debug)) {
              cat("\n[ UNCERTAINTY LEVEL", now_checks_for_uncertainty_level, "] (12) part of a name (very unlikely match)\n")
            }
            if (isTRUE(debug)) {
              message("Running '", f.x_withspaces_end_only, "'")
            }
            found <- reference_data_to_use[fullname_lower %like_case% f.x_withspaces_end_only]
            if (nrow(found) > 0) {
              found_result <- found[["mo"]]
              if (!empty_result(found_result) & nchar(g.x_backup_without_spp) >= 6) {
                found <- reference_data_to_use[mo == found_result[1L], ..property][[1]]
                uncertainties <<- rbind(uncertainties,
                                        format_uncertainty_as_df(uncertainty_level = now_checks_for_uncertainty_level,
                                                                 input = a.x_backup,
                                                                 result_mo = found_result[1L]))
                if (initial_search == TRUE) {
                  set_mo_history(a.x_backup, get_mo_code(found[1L], property), 3, force = force_mo_history, disable = disable_mo_history)
                }
                return(found[1L])
              }
            }
          }
          
          # didn't found in uncertain results too
          return(NA_character_)
        }
        
        # uncertain results
        # wait until prevalence == 2 to run the uncertain results on both prevalence == 1 and prevalence == 2
        if (nrow(data_to_check) == nrow(microorganismsDT[prevalence == 2])) {
          x[i] <- uncertain_fn(a.x_backup = a.x_backup, 
                               b.x_trimmed = b.x_trimmed,
                               d.x_withspaces_start_end = d.x_withspaces_start_end,
                               e.x_withspaces_start_only = e.x_withspaces_start_only, 
                               f.x_withspaces_end_only = f.x_withspaces_end_only,
                               g.x_backup_without_spp = g.x_backup_without_spp,
                               uncertain.reference_data_to_use = microorganismsDT[prevalence %in% c(1, 2)])
          if (!empty_result(x[i])) {
            # no set_mo_history here - it is already set in uncertain_fn()
            return(x[i])
          }
        } else if (nrow(data_to_check) == nrow(microorganismsDT[prevalence == 3])) {
          x[i] <- uncertain_fn(a.x_backup = a.x_backup, 
                               b.x_trimmed = b.x_trimmed,
                               d.x_withspaces_start_end = d.x_withspaces_start_end,
                               e.x_withspaces_start_only = e.x_withspaces_start_only, 
                               f.x_withspaces_end_only = f.x_withspaces_end_only,
                               g.x_backup_without_spp = g.x_backup_without_spp,
                               uncertain.reference_data_to_use = microorganismsDT[prevalence == 3])
          if (!empty_result(x[i])) {
            # no set_mo_history here - it is already set in uncertain_fn()
            return(x[i])
          }
        }
        
        # didn't found any
        return(NA_character_)
      }
      
      # FIRST TRY VERY PREVALENT IN HUMAN INFECTIONS ----
      x[i] <- check_per_prevalence(data_to_check = reference_data_to_use[prevalence == 1],
                                   data.old_to_check = microorganisms.oldDT[prevalence == 1],
                                   a.x_backup = x_backup[i],
                                   b.x_trimmed = x_trimmed[i],
                                   c.x_trimmed_without_group = x_trimmed_without_group[i],
                                   d.x_withspaces_start_end = x_withspaces_start_end[i],
                                   e.x_withspaces_start_only = x_withspaces_start_only[i],
                                   f.x_withspaces_end_only = x_withspaces_end_only[i],
                                   g.x_backup_without_spp = x_backup_without_spp[i],
                                   h.x_species = x_species[i],
                                   i.x_trimmed_species = x_trimmed_species[i])
      if (!empty_result(x[i])) {
        if (initial_search == TRUE) {
          set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history, disable = disable_mo_history)
        }
        next
      }
      
      # THEN TRY PREVALENT IN HUMAN INFECTIONS ----
      x[i] <- check_per_prevalence(data_to_check = reference_data_to_use[prevalence == 2],
                                   data.old_to_check = microorganisms.oldDT[prevalence %in% c(2, 3)], # run all other old MOs the second time,
                                   # otherwise e.g. mo_ref("Chlamydia psittaci") doesn't work correctly
                                   a.x_backup = x_backup[i],
                                   b.x_trimmed = x_trimmed[i],
                                   c.x_trimmed_without_group = x_trimmed_without_group[i],
                                   d.x_withspaces_start_end = x_withspaces_start_end[i],
                                   e.x_withspaces_start_only = x_withspaces_start_only[i],
                                   f.x_withspaces_end_only = x_withspaces_end_only[i],
                                   g.x_backup_without_spp = x_backup_without_spp[i],
                                   h.x_species = x_species[i],
                                   i.x_trimmed_species = x_trimmed_species[i])
      if (!empty_result(x[i])) {
        if (initial_search == TRUE) {
          set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history, disable = disable_mo_history)
        }
        next
      }
      
      # THEN UNPREVALENT IN HUMAN INFECTIONS ----
      x[i] <- check_per_prevalence(data_to_check = reference_data_to_use[prevalence == 3],
                                   data.old_to_check = microorganisms.oldDT[prevalence == 999],
                                   a.x_backup = x_backup[i],
                                   b.x_trimmed = x_trimmed[i],
                                   c.x_trimmed_without_group = x_trimmed_without_group[i],
                                   d.x_withspaces_start_end = x_withspaces_start_end[i],
                                   e.x_withspaces_start_only = x_withspaces_start_only[i],
                                   f.x_withspaces_end_only = x_withspaces_end_only[i],
                                   g.x_backup_without_spp = x_backup_without_spp[i],
                                   h.x_species = x_species[i],
                                   i.x_trimmed_species = x_trimmed_species[i])
      if (!empty_result(x[i])) {
        if (initial_search == TRUE) {
          set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history, disable = disable_mo_history)
        }
        next
      }
      
      # no results found: make them UNKNOWN ----
      x[i] <- microorganismsDT[mo == "UNKNOWN", ..property][[1]]
      if (initial_search == TRUE) {
        failures <- c(failures, x_backup[i])
        set_mo_history(x_backup[i], get_mo_code(x[i], property), 0, force = force_mo_history, disable = disable_mo_history)
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
                  " (covering ", percentage(total_failures / total_n),
                  ") could not be coerced and ", plural[3], " considered 'unknown'")
    if (n_distinct(failures) <= 10) {
      msg <- paste0(msg, ": ", paste('"', unique(failures), '"', sep = "", collapse = ", "))
    }
    msg <- paste0(msg,  ". Use mo_failures() to review ", plural[2], ". Edit the `allow_uncertain` parameter if needed (see ?as.mo).")
    warning(red(msg),
            call. = FALSE,
            immediate. = TRUE) # thus will always be shown, even if >= warnings
  }
  # handling uncertainties ----
  if (NROW(uncertainties) > 0 & initial_search == TRUE) {
    options(mo_uncertainties = as.list(distinct(uncertainties, input, .keep_all = TRUE)))
    
    plural <- c("", "it", "was")
    if (NROW(uncertainties) > 1) {
      plural <- c("s", "them", "were")
    }
    msg <- paste0("\nResult", plural[1], " of ", nr2char(NROW(uncertainties)), " value", plural[1],
                  " ", plural[3], " guessed with uncertainty. Use mo_uncertainties() to review ", plural[2], ".")
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
    
    x[x %in% CoNS] <- microorganismsDT[mo == "B_STPHY_CONS", ..property][[1]][1L]
    x[x %in% CoPS] <- microorganismsDT[mo == "B_STPHY_COPS", ..property][[1]][1L]
    if (Becker == "all") {
      x[x %in% microorganismsDT[mo %like_case% "^B_STPHY_AURS", ..property][[1]]] <- microorganismsDT[mo == "B_STPHY_COPS", ..property][[1]][1L]
    }
  }
  
  # Lancefield ----
  if (Lancefield == TRUE | Lancefield == "all") {
    # group A - S. pyogenes
    x[x == microorganismsDT[mo == "B_STRPT_PYGN", ..property][[1]][1L]] <- microorganismsDT[mo == "B_STRPT_GRPA", ..property][[1]][1L]
    # group B - S. agalactiae
    x[x == microorganismsDT[mo == "B_STRPT_AGLC", ..property][[1]][1L]] <- microorganismsDT[mo == "B_STRPT_GRPB", ..property][[1]][1L]
    # group C
    S_groupC <- microorganismsDT %>% filter(genus == "Streptococcus",
                                            species %in% c("equisimilis", "equi",
                                                           "zooepidemicus", "dysgalactiae")) %>%
      pull(property)
    x[x %in% S_groupC] <- microorganismsDT[mo == "B_STRPT_GRPC", ..property][[1]][1L]
    if (Lancefield == "all") {
      # all Enterococci
      x[x %like% "^(Enterococcus|B_ENTRC)"] <- microorganismsDT[mo == "B_STRPT_GRPD", ..property][[1]][1L]
    }
    # group F - S. anginosus
    x[x == microorganismsDT[mo == "B_STRPT_ANGN", ..property][[1]][1L]] <- microorganismsDT[mo == "B_STRPT_GRPF", ..property][[1]][1L]
    # group H - S. sanguinis
    x[x == microorganismsDT[mo == "B_STRPT_SNGN", ..property][[1]][1L]] <- microorganismsDT[mo == "B_STRPT_GRPH", ..property][[1]][1L]
    # group K - S. salivarius
    x[x == microorganismsDT[mo == "B_STRPT_SLVR", ..property][[1]][1L]] <- microorganismsDT[mo == "B_STRPT_GRPK", ..property][[1]][1L]
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
    x <- to_class_mo(x)
  }
  
  if (length(mo_renamed()) > 0) {
    print(mo_renamed())
  }
  
  if (old_mo_warning == TRUE & property != "mo") {
    warning("The input contained old microorganism IDs from previous versions of this package.\nPlease use `as.mo()` on these old IDs to transform them to the new format.\nSUPPORT FOR THIS WILL BE DROPPED IN A FUTURE VERSION.", call. = FALSE)
  }
  
  x
}

empty_result <- function(x) {
  all(x %in% c(NA, "UNKNOWN"))
}

#' @importFrom crayon italic
was_renamed <- function(name_old, name_new, ref_old = "", ref_new = "", mo = "") {
  newly_set <- data.frame(old_name = name_old, 
                          old_ref = ref_old,
                          new_name = name_new,
                          new_ref = ref_new,
                          mo = mo, 
                          stringsAsFactors = FALSE)
  already_set <- getOption("mo_renamed")
  if (!is.null(already_set)) {
    options(mo_renamed = rbind(already_set, newly_set))
  } else {
    options(mo_renamed = newly_set)
  }
}

format_uncertainty_as_df <- function(uncertainty_level,
                                     input,
                                     result_mo) {
  if (!is.null(getOption("mo_renamed_last_run", default = NULL))) {
    # was found as a renamed mo
    df <- data.frame(uncertainty = uncertainty_level,
                     input = input,
                     fullname = getOption("mo_renamed_last_run"),
                     renamed_to = microorganismsDT[mo == result_mo, fullname][[1]],
                     mo = result_mo,
                     stringsAsFactors = FALSE)
    options(mo_renamed_last_run = NULL)
  } else {
    df <- data.frame(uncertainty = uncertainty_level,
                     input = input,
                     fullname = microorganismsDT[mo == result_mo, fullname][[1]],
                     renamed_to = NA_character_,
                     mo = result_mo,
                     stringsAsFactors = FALSE)
  }
  df
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

#' @importFrom pillar type_sum
#' @export
type_sum.mo <- function(x) {
  "mo"
}

#' @importFrom pillar pillar_shaft
#' @export
pillar_shaft.mo <- function(x, ...) {
  out <- format(x)
  # grey out the kingdom (part until first "_")
  out[!is.na(x)] <- gsub("^([A-Z]+_)(.*)", paste0(pillar::style_subtle("\\1"), "\\2"), out[!is.na(x)])
  # and grey out every _
  out[!is.na(x)] <- gsub("_", pillar::style_subtle("_"), out[!is.na(x)])
  
  # markup NA and UNKNOWN
  out[is.na(x)] <- pillar::style_na("  NA")
  out[x == "UNKNOWN"] <- pillar::style_na("  UNKNOWN")
  
  # make it always fit exactly
  pillar::new_pillar_shaft_simple(out, align = "left", width = max(nchar(x)))
}

#' @exportMethod summary.mo
#' @importFrom dplyr n_distinct
#' @importFrom cleaner freq top_freq
#' @export
#' @noRd
summary.mo <- function(object, ...) {
  # unique and top 1-3
  x <- as.mo(object)
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

#' @exportMethod [.mo
#' @export
#' @noRd
"[.mo" <- function(x, ...) {
  y <- NextMethod()
  attributes(y) <- attributes(x)
  y
}
#' @exportMethod [[.mo
#' @export
#' @noRd
"[[.mo" <- function(x, ...) {
  y <- NextMethod()
  attributes(y) <- attributes(x)
  y
}
#' @exportMethod [<-.mo
#' @export
#' @noRd
"[<-.mo" <- function(i, j, ..., value) {
  y <- NextMethod()
  attributes(y) <- attributes(i)
  class_integrity_check(y, "microorganism code", c(as.character(AMR::microorganisms$mo), as.character(microorganisms.translation$mo_old)))
}
#' @exportMethod [[<-.mo
#' @export
#' @noRd
"[[<-.mo" <- function(i, j, ..., value) {
  y <- NextMethod()
  attributes(y) <- attributes(i)
  class_integrity_check(y, "microorganism code", c(as.character(AMR::microorganisms$mo), as.character(microorganisms.translation$mo_old)))
}
#' @exportMethod c.mo
#' @export
#' @noRd
c.mo <- function(x, ...) {
  y <- NextMethod()
  attributes(y) <- attributes(x)
  class_integrity_check(y, "microorganism code", c(as.character(AMR::microorganisms$mo), as.character(microorganisms.translation$mo_old)))
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
#' @importFrom crayon green yellow red white black bgGreen bgYellow bgRed
#' @importFrom cleaner percentage
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
  for (i in seq_len(nrow(x))) {
    if (x[i, "uncertainty"] == 1) {
      colour1 <- green
      colour2 <- function(...) bgGreen(white(...))
    } else if (x[i, "uncertainty"] == 2) {
      colour1 <- yellow
      colour2 <- function(...) bgYellow(black(...))
    } else {
      colour1 <- red
      colour2 <- function(...) bgRed(white(...))
    }
    msg <- paste(msg,
                 paste0(colour2(paste0(" [", x[i, "uncertainty"], "] ")), ' "', x[i, "input"], '" -> ',
                        colour1(paste0(italic(x[i, "fullname"]),
                                       ifelse(!is.na(x[i, "renamed_to"]), paste(", renamed to", italic(x[i, "renamed_to"])), ""),
                                       " (", x[i, "mo"],
                                       ", score: ", percentage(levenshtein_fraction(x[i, "input"], x[i, "fullname"]), digits = 1),
                                       ")"))),
                 sep = "\n")
  }
  cat(msg)
}

#' @rdname as.mo
#' @importFrom dplyr distinct
#' @export
mo_renamed <- function() {
  items <- getOption("mo_renamed")
  if (is.null(items)) {
    items <- data.frame()
  } else {
    items <- distinct(items, old_name, .keep_all = TRUE)
  }
  structure(.Data = items,
            class = c("mo_renamed", "data.frame"))
}

#' @exportMethod print.mo_renamed
#' @importFrom crayon blue italic
#' @export
#' @noRd
print.mo_renamed <- function(x, ...) {
  if (NROW(x) == 0) {
    return(invisible())
  }
  for (i in seq_len(nrow(x))) {
    message(blue(paste0("NOTE: ", 
                        italic(x$old_name[i]), ifelse(x$old_ref[i] %in% c("", NA), "", 
                                                      paste0(" (",  gsub("et al.", italic("et al."), x$old_ref[i]), ")")),
                        " was renamed ", 
                        italic(x$new_name[i]), ifelse(x$new_ref[i] %in% c("", NA), "", 
                                                      paste0(" (",  gsub("et al.", italic("et al."), x$new_ref[i]), ")")),
                        " [", x$mo[i], "]")))
  }
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
  if (property == "mo") {
    unique(x)
  } else {
    microorganismsDT[get(property) == x, "mo"][[1]]
  }
}

translate_allow_uncertain <- function(allow_uncertain) {
  if (isTRUE(allow_uncertain)) {
    # default to uncertainty level 2
    allow_uncertain <- 2
  } else {
    allow_uncertain[tolower(allow_uncertain) == "none"] <- 0
    allow_uncertain[tolower(allow_uncertain) == "all"] <- 3
    allow_uncertain <- as.integer(allow_uncertain)
    if (!allow_uncertain %in% c(0:3)) {
      stop('`allow_uncertain` must be a number between 0 (or "none") and 3 (or "all"), or TRUE (= 2) or FALSE (= 0).', call. = FALSE)
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

#' @importFrom utils adist
levenshtein_fraction <- function(input, output) {
  levenshtein <- double(length = length(input))
  for (i in seq_len(length(input))) {
    # determine levenshtein distance, but maximise to nchar of output
    levenshtein[i] <- base::min(base::as.double(adist(input[i], output[i], ignore.case = TRUE)),
                                base::nchar(output[i]))
  }
  # self-made score between 0 and 1 (for % certainty, so 0 means huge distance, 1 means no distance)
  (base::nchar(output) - 0.5 * levenshtein) / nchar(output)
}
