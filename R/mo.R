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

#' Transform to microorganism ID
#'
#' Use this function to determine a valid microorganism ID (\code{mo}). Determination is done using Artificial Intelligence (AI) and the complete taxonomic kingdoms \emph{Bacteria}, \emph{Fungi} and \emph{Protozoa} (see Source), so the input can be almost anything: a full name (like \code{"Staphylococcus aureus"}), an abbreviated name (like \code{"S. aureus"}), an abbreviation known in the field (like \code{"MRSA"}), or just a genus. You could also \code{\link{select}} a genus and species column, zie Examples.
#' @param x a character vector or a \code{data.frame} with one or two columns
#' @param Becker a logical to indicate whether \emph{Staphylococci} should be categorised into Coagulase Negative \emph{Staphylococci} ("CoNS") and Coagulase Positive \emph{Staphylococci} ("CoPS") instead of their own species, according to Karsten Becker \emph{et al.} [1].
#'
#'   This excludes \emph{Staphylococcus aureus} at default, use \code{Becker = "all"} to also categorise \emph{S. aureus} as "CoPS".
#' @param Lancefield a logical to indicate whether beta-haemolytic \emph{Streptococci} should be categorised into Lancefield groups instead of their own species, according to Rebecca C. Lancefield [2]. These \emph{Streptococci} will be categorised in their first group, e.g. \emph{Streptococcus dysgalactiae} will be group C, although officially it was also categorised into groups G and L.
#'
#'   This excludes \emph{Enterococci} at default (who are in group D), use \code{Lancefield = "all"} to also categorise all \emph{Enterococci} as group D.
#' @param allow_uncertain a logical to indicate whether the input should be checked for less possible results, see Details
#' @param reference_df a \code{data.frame} to use for extra reference when translating \code{x} to a valid \code{mo}. See \code{\link{set_mo_source}} and \code{\link{get_mo_source}} to automate the usage of your own codes (e.g. used in your analysis or organisation).
#' @rdname as.mo
#' @aliases mo
#' @keywords mo Becker becker Lancefield lancefield guess
#' @details
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
#'    ----> taxonomic kingdom: A (Archaea), B (Bacteria), C (Chromista),
#'                             F (Fungi), P (Protozoa) or V (Viruses)
#' }
#'
#' Use the \code{\link{mo_property}} functions to get properties based on the returned code, see Examples.
#'
#' This function uses Artificial Intelligence (AI) to help getting fast and logical results. It tries to find matches in this order:
#' \itemize{
#'   \item{Taxonomic kingdom: it first searches in Bacteria, then Fungi, then Protozoa}
#'   \item{Human pathogenic prevalence: it first searches in more prevalent microorganisms, then less prevalent ones (see section \emph{Microbial prevalence of pathogens in humans})}
#'   \item{Valid MO codes and full names: it first searches in already valid MO code and known genus/species combinations}
#'   \item{Breakdown of input values: from here it starts to breakdown input values to find possible matches}
#' }
#'
#' A couple of effects because of these rules:
#' \itemize{
#'   \item{\code{"E. coli"} will return the ID of \emph{Escherichia coli} and not \emph{Entamoeba coli}, although the latter would alphabetically come first}
#'   \item{\code{"H. influenzae"} will return the ID of \emph{Haemophilus influenzae} and not \emph{Haematobacter influenzae} for the same reason}
#'   \item{Something like \code{"p aer"} will return the ID of \emph{Pseudomonas aeruginosa} and not \emph{Pasteurella aerogenes}}
#'   \item{Something like \code{"stau"} or \code{"S aur"} will return the ID of \emph{Staphylococcus aureus} and not \emph{Staphylococcus auricularis}}
#' }
#' This means that looking up human pathogenic microorganisms takes less time than looking up human \strong{non}-pathogenic microorganisms.
#'
#' \strong{UNCERTAIN RESULTS} \cr
#' When using \code{allow_uncertain = TRUE} (which is the default setting), it will use additional rules if all previous AI rules failed to get valid results. These are:
#' \itemize{
#'   \item{It tries to look for previously accepted (but now invalid) taxonomic names}
#'   \item{It strips off values between brackets and the brackets itself, and re-evaluates the input with all previous rules}
#'   \item{It strips off words from the end one by one and re-evaluates the input with all previous rules}
#'   \item{It strips off words from the start one by one and re-evaluates the input with all previous rules}
#'   \item{It tries to look for some manual changes which are not yet published to the Catalogue of Life (like \emph{Propionibacterium} not yet being \emph{Cutibacterium})}
#' }
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
#' Use \code{mo_uncertainties()} to get a vector with all values that were coerced to a valid value, but with uncertainty.
#'
#' Use \code{mo_renamed()} to get a vector with all values that could be coerced based on an old, previously accepted taxonomic name.
#'
#' @section Microbial prevalence of pathogens in humans:
#' The artificial intelligence takes into account microbial prevalence of pathogens in humans. It uses three groups and every (sub)species is in the group it matches first. These groups are:
#' \itemize{
#'   \item{1 (most prevalent): class is Gammaproteobacteria \strong{or} genus is one of: \emph{Enterococcus}, \emph{Staphylococcus}, \emph{Streptococcus}.}
#'   \item{2: phylum is one of: Proteobacteria, Firmicutes, Actinobacteria, Sarcomastigophora \strong{or} genus is one of: \emph{Aspergillus}, \emph{Bacteroides}, \emph{Candida}, \emph{Capnocytophaga}, \emph{Chryseobacterium}, \emph{Cryptococcus}, \emph{Elisabethkingia}, \emph{Flavobacterium}, \emph{Fusobacterium}, \emph{Giardia}, \emph{Leptotrichia}, \emph{Mycoplasma}, \emph{Prevotella}, \emph{Rhodotorula}, \emph{Treponema}, \emph{Trichophyton}, \emph{Ureaplasma}.}
#'   \item{3 (least prevalent): all others.}
#' }
#'
#' Group 1 contains all common Gram negatives, like all Enterobacteriaceae and e.g. \emph{Pseudomonas} and \emph{Legionella}.
#'
#' Group 2 probably contains all microbial pathogens ever found in humans.
#' @inheritSection catalogue_of_life Catalogue of Life
#  (source as a section, so it can be inherited by other man pages)
#' @section Source:
#' [1] Becker K \emph{et al.} \strong{Coagulase-Negative Staphylococci}. 2014. Clin Microbiol Rev. 27(4): 870–926. \url{https://dx.doi.org/10.1128/CMR.00109-13}
#'
#' [2] Lancefield RC \strong{A serological differentiation of human and other groups of hemolytic streptococci}. 1933. J Exp Med. 57(4): 571–95. \url{https://dx.doi.org/10.1084/jem.57.4.571}
#'
#' [3] Catalogue of Life: Annual Checklist (public online database), \url{www.catalogueoflife.org}.
#' @export
#' @return Character (vector) with class \code{"mo"}. Unknown values will return \code{NA}.
#' @seealso \code{\link{microorganisms}} for the \code{data.frame} that is being used to determine ID's. \cr
#' The \code{\link{mo_property}} functions (like \code{\link{mo_genus}}, \code{\link{mo_gramstain}}) to get properties based on the returned code.
#' @inheritSection AMR Read more on our website!
#' @examples
#' # These examples all return "B_STPHY_AUR", the ID of S. aureus:
#' as.mo("stau")
#' as.mo("STAU")
#' as.mo("staaur")
#' as.mo("S. aureus")
#' as.mo("S aureus")
#' as.mo("Staphylococcus aureus")
#' as.mo("Staphylococcus aureus (MRSA)")
#' as.mo("MRSA") # Methicillin Resistant S. aureus
#' as.mo("VISA") # Vancomycin Intermediate S. aureus
#' as.mo("VRSA") # Vancomycin Resistant S. aureus
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
#' # Use mo_* functions to get a specific property based on `mo`
#' Ecoli <- as.mo("E. coli")     # returns `B_ESCHR_COL`
#' mo_genus(Ecoli)               # returns "Escherichia"
#' mo_gramstain(Ecoli)           # returns "Gram negative"
#' # but it uses as.mo internally too, so you could also just use:
#' mo_genus("E. coli")           # returns "Escherichia"
#'
#'
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
as.mo <- function(x, Becker = FALSE, Lancefield = FALSE, allow_uncertain = TRUE, reference_df = get_mo_source()) {
  # will be checked for mo class in validation
  mo_validate(x = x, property = "mo",
              Becker = Becker, Lancefield = Lancefield,
              allow_uncertain = allow_uncertain, reference_df = reference_df)
}

#' @rdname as.mo
#' @export
is.mo <- function(x) {
  identical(class(x), "mo")
}

#' @importFrom dplyr %>% pull left_join n_distinct progress_estimated filter
#' @importFrom data.table data.table as.data.table setkey
#' @importFrom crayon magenta red blue silver italic has_color
exec_as.mo <- function(x, Becker = FALSE, Lancefield = FALSE,
                       allow_uncertain = TRUE, reference_df = get_mo_source(),
                       property = "mo", clear_options = TRUE) {

  if (!"AMR" %in% base::.packages()) {
    library("AMR")
    # check onLoad() in R/zzz.R: data tables are created there.
  }

  if (clear_options == TRUE) {
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
  uncertainties <- character(0)
  failures <- character(0)
  x_input <- x
  # already strip leading and trailing spaces
  x <- trimws(x, which = "both")
  # only check the uniques, which is way faster
  x <- unique(x)
  # remove empty values (to later fill them in again with NAs)
  x <- x[!is.na(x) & !is.null(x) & !identical(x, "")]

  # conversion of old MO codes from v0.5.0 (ITIS) to later versions (Catalogue of Life)
  if (any(x %like% "^[BFP]_[A-Z]{3,7}") & !all(x %in% microorganisms$mo)) {
    leftpart <- gsub("^([BFP]_[A-Z]{3,7}).*", "\\1", x)
    if (any(leftpart %in% names(mo_codes_v0.5.0))) {
      rightpart <- gsub("^[BFP]_[A-Z]{3,7}(.*)", "\\1", x)
      leftpart <- mo_codes_v0.5.0[leftpart]
      x[!is.na(leftpart)] <- paste0(leftpart[!is.na(leftpart)], rightpart[!is.na(leftpart)])
    }
  }

  # defined df to check for
  if (!is.null(reference_df)) {
    if (!is.data.frame(reference_df) | NCOL(reference_df) < 2) {
      stop('`reference_df` must be a data.frame with at least two columns.', call. = FALSE)
    }
    if (!"mo" %in% colnames(reference_df)) {
      stop("`reference_df` must contain a column `mo` with values from the 'microorganisms' data set.", call. = FALSE)
    }
    reference_df <- reference_df %>% filter(!is.na(mo))
    # # remove factors, just keep characters
    suppressWarnings(
      reference_df[] <- lapply(reference_df, as.character)
    )
  }

  # all empty
  if (all(identical(trimws(x_input), "") | is.na(x_input))) {
    if (property == "mo") {
      return(structure(rep(NA_character_, length(x_input)),
                       class = "mo"))
    } else {
      return(rep(NA_character_, length(x_input)))
    }

  } else if (all(x %in% reference_df[, 1])
             & all(reference_df[, "mo"] %in% AMR::microorganisms$mo)) {
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
      y[is.na(y)] <- microorganismsDT[prevalence == 2][data.table(mo = x[is.na(y)]), on = "mo", ..property][[1]]
    }
    if (any(is.na(y))) {
      y[is.na(y)] <- microorganismsDT[prevalence == 3][data.table(mo = x[is.na(y)]), on = "mo", ..property][[1]]
    }
    x <- y

  } else if (all(x %in% AMR::microorganisms$fullname)) {
    # we need special treatment for very prevalent full names, they are likely!
    # e.g. as.mo("Staphylococcus aureus")
    y <- microorganismsDT[prevalence == 1][data.table(fullname = x), on = "fullname", ..property][[1]]
    if (any(is.na(y))) {
      y[is.na(y)] <- microorganismsDT[prevalence == 2][data.table(fullname = x[is.na(y)]), on = "fullname", ..property][[1]]
    }
    if (any(is.na(y))) {
      y[is.na(y)] <- microorganismsDT[prevalence == 3][data.table(fullname = x[is.na(y)]), on = "fullname", ..property][[1]]
    }
    x <- y

  } else if (all(toupper(x) %in% AMR::microorganisms.codes$code)) {
    # commonly used MO codes
    y <- as.data.table(AMR::microorganisms.codes)[data.table(code = toupper(x)), on = "code", ]
    x <- microorganismsDT[data.table(mo = y[["mo"]]), on = "mo", ..property][[1]]

  } else if (!all(x %in% AMR::microorganisms[, property])) {

    x_backup <- x

    # remove spp and species
    x <- trimws(gsub(" +(spp.?|ssp.?|sp.? |ss ?.?|subsp.?|subspecies|biovar |serovar |species)", " ", x_backup, ignore.case = TRUE), which = "both")
    x_species <- paste(x, "species")
    # translate to English for supported languages of mo_property
    x <- gsub("(Gruppe|gruppe|groep|grupo|gruppo|groupe)", "group", x, ignore.case = TRUE)
    # remove 'empty' genus and species values
    x <- gsub("(no MO)", "", x, fixed = TRUE)
    # remove non-text in case of "E. coli" except dots and spaces
    x <- gsub("[^.a-zA-Z0-9/ \\-]+", "", x)
    # replace minus by a space
    x <- gsub("-+", " ", x)
    # replace hemolytic by haemolytic
    x <- gsub("ha?emoly", "haemoly", x)
    # place minus back in streptococci
    x <- gsub("(alpha|beta|gamma) ha?emoly", "\\1-haemoly", x)
    # remove genus as first word
    x <- gsub("^Genus ", "", x)

    # but spaces before and after should be omitted
    x <- trimws(x, which = "both")
    x_trimmed <- x
    x_trimmed_species <- paste(x_trimmed, "species")
    x_trimmed_without_group <- gsub(" group$", "", x_trimmed, ignore.case = TRUE)
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

    # cat(paste0('x                       "', x, '"\n'))
    # cat(paste0('x_species               "', x_species, '"\n'))
    # cat(paste0('x_withspaces_start_only "', x_withspaces_start_only, '"\n'))
    # cat(paste0('x_withspaces_end_only   "', x_withspaces_end_only, '"\n'))
    # cat(paste0('x_withspaces_start_end  "', x_withspaces_start_end, '"\n'))
    # cat(paste0('x_backup                "', x_backup, '"\n'))
    # cat(paste0('x_trimmed               "', x_trimmed, '"\n'))
    # cat(paste0('x_trimmed_species       "', x_trimmed_species, '"\n'))
    # cat(paste0('x_trimmed_without_group "', x_trimmed_without_group, '"\n'))

    progress <- progress_estimated(n = length(x), min_time = 3)

    for (i in 1:length(x)) {

      progress$tick()$print()

      found <- microorganismsDT[mo == toupper(x_backup[i]), ..property][[1]]
      # is a valid MO code
      if (length(found) > 0) {
        x[i] <- found[1L]
        next
      }

      if (tolower(x_trimmed[i]) %in% c("", "xxx", "other", "none", "unknown")) {
        # empty and nonsense values, ignore without warning ("xxx" is WHONET code for 'no growth')
        x[i] <- NA_character_
        next
      }

      if (nchar(gsub("[^a-zA-Z]", "", x_trimmed[i])) < 3) {
        # check if search term was like "A. species", then return first genus found with ^A
        if (x_backup[i] %like% "[a-z]+ species" | x_backup[i] %like% "[a-z] spp[.]?") {
          # get mo code of first hit
          found <- microorganismsDT[fullname %like% x_withspaces_start_only[i], mo]
          if (length(found) > 0) {
            mo_code <- found[1L] %>% strsplit("_") %>% unlist() %>% .[1:2] %>% paste(collapse = "_")
            found <- microorganismsDT[mo == mo_code, ..property][[1]]
            # return first genus that begins with x_trimmed, e.g. when "E. spp."
            if (length(found) > 0) {
              x[i] <- found[1L]
              next
            }
          }
        }
        # fewer than 3 chars and not looked for species, add as failure
        x[i] <- NA_character_
        failures <- c(failures, x_backup[i])
        next
      }

      if (x_trimmed[i] %like% "virus") {
        # there is no fullname like virus, so don't try to coerce it
        x[i] <- NA_character_
        failures <- c(failures, x_backup[i])
        next
      }

      # translate known trivial abbreviations to genus + species ----
      if (!is.na(x_trimmed[i])) {
        if (toupper(x_trimmed[i]) %in% c('MRSA', 'MSSA', 'VISA', 'VRSA')) {
          x[i] <- microorganismsDT[mo == 'B_STPHY_AUR', ..property][[1]][1L]
          next
        }
        if (toupper(x_trimmed[i]) %in% c('MRSE', 'MSSE')) {
          x[i] <- microorganismsDT[mo == 'B_STPHY_EPI', ..property][[1]][1L]
          next
        }
        if (toupper(x_trimmed[i]) == "VRE"
            | x_trimmed[i] %like% '(enterococci|enterokok|enterococo)[a-z]*?$')  {
          x[i] <- microorganismsDT[mo == 'B_ENTRC', ..property][[1]][1L]
          next
        }
        if (toupper(x_trimmed[i]) %in% c("EHEC", "EPEC", "EIEC", "STEC", "ATEC")) {
          x[i] <- microorganismsDT[mo == 'B_ESCHR_COL', ..property][[1]][1L]
          next
        }
        if (toupper(x_trimmed[i]) == 'MRPA') {
          # multi resistant P. aeruginosa
          x[i] <- microorganismsDT[mo == 'B_PSDMN_AER', ..property][[1]][1L]
          next
        }
        if (toupper(x_trimmed[i]) == 'CRS'
            | toupper(x_trimmed[i]) == 'CRSM') {
          # co-trim resistant S. maltophilia
          x[i] <- microorganismsDT[mo == 'B_STNTR_MAL', ..property][[1]][1L]
          next
        }
        if (toupper(x_trimmed[i]) %in% c('PISP', 'PRSP', 'VISP', 'VRSP')) {
          # peni I, peni R, vanco I, vanco R: S. pneumoniae
          x[i] <- microorganismsDT[mo == 'B_STRPT_PNE', ..property][[1]][1L]
          next
        }
        if (toupper(x_trimmed[i]) %like% '^G[ABCDFGHK]S$') {
          # Streptococci, like GBS = Group B Streptococci (B_STRPT_GRB)
          x[i] <- microorganismsDT[mo == gsub("G([ABCDFGHK])S", "B_STRPT_GR\\1", x_trimmed[i], ignore.case = TRUE), ..property][[1]][1L]
          next
        }
        if (toupper(x_trimmed[i]) %like% '(streptococc|streptokok).* [ABCDFGHK]$') {
          # Streptococci in different languages, like "estreptococos grupo B"
          x[i] <- microorganismsDT[mo == gsub(".*(streptococ|streptokok|estreptococ).* ([ABCDFGHK])$", "B_STRPT_GR\\2", x_trimmed[i], ignore.case = TRUE), ..property][[1]][1L]
          next
        }
        if (toupper(x_trimmed[i]) %like% 'group [ABCDFGHK] (streptococ|streptokok|estreptococ)') {
          # Streptococci in different languages, like "Group A Streptococci"
          x[i] <- microorganismsDT[mo == gsub(".*group ([ABCDFGHK]) (streptococ|streptokok|estreptococ).*", "B_STRPT_GR\\1", x_trimmed[i], ignore.case = TRUE), ..property][[1]][1L]
          next
        }
        # CoNS/CoPS in different languages (support for German, Dutch, Spanish, Portuguese) ----
        if (tolower(x[i]) %like% '[ck]oagulas[ea] negatie?[vf]'
            | tolower(x_trimmed[i]) %like% '[ck]oagulas[ea] negatie?[vf]'
            | tolower(x[i]) %like% '[ck]o?ns[^a-z]?$') {
          # coerce S. coagulase negative
          x[i] <- microorganismsDT[mo == 'B_STPHY_CNS', ..property][[1]][1L]
          next
        }
        if (tolower(x[i]) %like% '[ck]oagulas[ea] positie?[vf]'
            | tolower(x_trimmed[i]) %like% '[ck]oagulas[ea] positie?[vf]'
            | tolower(x[i]) %like% '[ck]o?ps[^a-z]?$') {
          # coerce S. coagulase positive
          x[i] <- microorganismsDT[mo == 'B_STPHY_CPS', ..property][[1]][1L]
          next
        }
        if (tolower(x[i]) %like% 'gram[ -]?neg.*'
            | tolower(x_trimmed[i]) %like% 'gram[ -]?neg.*') {
          # coerce S. coagulase positive
          x[i] <- microorganismsDT[mo == 'B_GRAMN', ..property][[1]][1L]
          next
        }
        if (tolower(x[i]) %like% 'gram[ -]?pos.*'
            | tolower(x_trimmed[i]) %like% 'gram[ -]?pos.*') {
          # coerce S. coagulase positive
          x[i] <- microorganismsDT[mo == 'B_GRAMP', ..property][[1]][1L]
          next
        }
        if (grepl("[sS]almonella [A-Z][a-z]+ ?.*", x_trimmed[i])) {
          if (x_trimmed[i] %like% "Salmonella group") {
            # Salmonella Group A to Z, just return S. species for now
            x[i] <- microorganismsDT[mo == 'B_SLMNL', ..property][[1]][1L]
            notes <- c(notes,
                       magenta(paste0("Note: ",
                                      italic("Salmonella"), " ", trimws(gsub("Salmonella", "", x_trimmed[i])),
                                      " was considered ",
                                      italic("Salmonella species"),
                                      " (B_SLMNL)")))
          } else {
            # Salmonella with capital letter species like "Salmonella Goettingen" - they're all S. enterica
            x[i] <- microorganismsDT[mo == 'B_SLMNL_ENT', ..property][[1]][1L]
            notes <- c(notes,
                       magenta(paste0("Note: ",
                                      italic("Salmonella"), " ", trimws(gsub("Salmonella", "", x_trimmed[i])),
                                      " was considered a subspecies of ",
                                      italic("Salmonella enterica"),
                                      " (B_SLMNL_ENT)")))
          }
          next
        }
      }

      # FIRST TRY FULLNAMES AND CODES
      # if only genus is available, return only genus
      if (all(!c(x[i], x_trimmed[i]) %like% " ")) {
        found <- microorganismsDT[tolower(fullname) %in% tolower(c(x_species[i], x_trimmed_species[i])), ..property][[1]]
        if (length(found) > 0) {
          x[i] <- found[1L]
          next
        }
        if (nchar(x_trimmed[i]) >= 6) {
          found <- microorganismsDT[tolower(fullname) %like% paste0(x_withspaces_start_only[i], "[a-z]+ species"), ..property][[1]]
          if (length(found) > 0) {
            x[i] <- found[1L]
            next
          }
        }
        # rest of genus only is in allow_uncertain part.
      }

      # TRY OTHER SOURCES ----
      if (toupper(x_backup[i]) %in% AMR::microorganisms.codes[, 1]) {
        mo_found <- AMR::microorganisms.codes[toupper(x_backup[i]) == AMR::microorganisms.codes[, 1], "mo"][1L]
        if (length(mo_found) > 0) {
          x[i] <- microorganismsDT[mo == mo_found, ..property][[1]][1L]
          next
        }
      }
      if (!is.null(reference_df)) {
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

      check_per_prevalence <- function(data_to_check,
                                       a.x_backup,
                                       b.x_trimmed,
                                       c.x_trimmed_without_group,
                                       d.x_withspaces_start_end,
                                       e.x_withspaces_start_only,
                                       f.x_withspaces_end_only) {

        found <- data_to_check[tolower(fullname) %in% tolower(c(a.x_backup, b.x_trimmed)), ..property][[1]]
        # most probable: is exact match in fullname
        if (length(found) > 0) {
          return(found[1L])
        }

        found <- data_to_check[tolower(fullname) == tolower(c.x_trimmed_without_group), ..property][[1]]
        if (length(found) > 0) {
          return(found[1L])
        }

        # try any match keeping spaces ----
        found <- data_to_check[fullname %like% d.x_withspaces_start_end, ..property][[1]]
        if (length(found) > 0 & nchar(b.x_trimmed) >= 6) {
          return(found[1L])
        }

        # try any match keeping spaces, not ending with $ ----
        found <- data_to_check[fullname %like% paste0(trimws(e.x_withspaces_start_only), " "), ..property][[1]]
        if (length(found) > 0) {
          return(found[1L])
        }
        found <- data_to_check[fullname %like% e.x_withspaces_start_only, ..property][[1]]
        if (length(found) > 0 & nchar(b.x_trimmed) >= 6) {
          return(found[1L])
        }

        # try any match keeping spaces, not start with ^ ----
        found <- data_to_check[fullname %like% paste0(" ", trimws(f.x_withspaces_end_only)), ..property][[1]]
        if (length(found) > 0) {
          return(found[1L])
        }
        found <- data_to_check[fullname %like% f.x_withspaces_end_only, ..property][[1]]
        if (length(found) > 0 & nchar(b.x_trimmed) >= 6) {
          return(found[1L])
        }

        # try splitting of characters in the middle and then find ID ----
        # only when text length is 6 or lower
        # like esco = E. coli, klpn = K. pneumoniae, stau = S. aureus, staaur = S. aureus
        if (nchar(b.x_trimmed) <= 6) {
          x_length <- nchar(b.x_trimmed)
          x_split <- paste0("^",
                            b.x_trimmed %>% substr(1, x_length / 2),
                            '.* ',
                            b.x_trimmed %>% substr((x_length / 2) + 1, x_length))
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
                                   f.x_withspaces_end_only = x_withspaces_end_only[i])
      if (!is.na(x[i])) {
        next
      }
      # THEN TRY PREVALENT IN HUMAN INFECTIONS ----
      x[i] <- check_per_prevalence(data_to_check = microorganismsDT[prevalence == 2],
                                   a.x_backup = x_backup[i],
                                   b.x_trimmed = x_trimmed[i],
                                   c.x_trimmed_without_group = x_trimmed_without_group[i],
                                   d.x_withspaces_start_end = x_withspaces_start_end[i],
                                   e.x_withspaces_start_only = x_withspaces_start_only[i],
                                   f.x_withspaces_end_only = x_withspaces_end_only[i])
      if (!is.na(x[i])) {
        next
      }
      # THEN UNPREVALENT IN HUMAN INFECTIONS ----
      x[i] <- check_per_prevalence(data_to_check = microorganismsDT[prevalence == 3],
                                   a.x_backup = x_backup[i],
                                   b.x_trimmed = x_trimmed[i],
                                   c.x_trimmed_without_group = x_trimmed_without_group[i],
                                   d.x_withspaces_start_end = x_withspaces_start_end[i],
                                   e.x_withspaces_start_only = x_withspaces_start_only[i],
                                   f.x_withspaces_end_only = x_withspaces_end_only[i])
      if (!is.na(x[i])) {
        next
      }

      # MISCELLANEOUS ----

      # look for old taxonomic names ----
      found <- microorganisms.oldDT[tolower(fullname) == tolower(x_backup[i])
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
        next
      }

      # check for uncertain results ----
      if (allow_uncertain == TRUE) {

        uncertain_fn <- function(a.x_backup, b.x_trimmed, c.x_withspaces_start_end, d.x_withspaces_start_only) {

          # (1) look for genus only, part of name ----
          if (nchar(b.x_trimmed) > 4 & !b.x_trimmed %like% " ") {
            if (!grepl("^[A-Z][a-z]+", b.x_trimmed, ignore.case = FALSE)) {
              # not when input is like Genustext, because then Neospora would lead to Actinokineospora
              found <- microorganismsDT[tolower(fullname) %like% paste(b.x_trimmed, "species"), ..property][[1]]
              if (length(found) > 0) {
                x[i] <- found[1L]
                uncertainties <<- c(uncertainties,
                                    paste0("'", a.x_backup, "' >> ", microorganismsDT[mo == found[1L], fullname][[1]], " (", found[1L], ")"))
                return(x)
              }
            }
          }

          # (2) look again for old taxonomic names, now for G. species ----
          found <- microorganisms.oldDT[fullname %like% c.x_withspaces_start_end
                                        | fullname %like% d.x_withspaces_start_only]
          if (NROW(found) > 0 & nchar(b.x_trimmed) >= 6) {
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
            uncertainties <<- c(uncertainties,
                                paste0("'", a.x_backup, "' >> ", found[1, fullname], " (Catalogue of Life ID ", found[1, col_id], ")"))
            return(x)
          }

          # (3) strip values between brackets ----
          a.x_backup_stripped <- gsub("( [(].*[)])", "", a.x_backup)
          a.x_backup_stripped <- trimws(gsub("  ", " ", a.x_backup_stripped, fixed = TRUE))
          found <- suppressMessages(suppressWarnings(exec_as.mo(a.x_backup_stripped, clear_options = FALSE, allow_uncertain = FALSE)))
          if (!is.na(found) & nchar(b.x_trimmed) >= 6) {
            found_result <- found
            found <- microorganismsDT[mo == found, ..property][[1]]
            uncertainties <<- c(uncertainties,
                                paste0("'", a.x_backup, "' >> ", microorganismsDT[mo == found_result[1L], fullname][[1]], " (", found_result[1L], ")"))
            return(found[1L])
          }

          # (4) try to strip off one element from end and check the remains ----
          x_strip <- a.x_backup %>% strsplit(" ") %>% unlist()
          if (length(x_strip) > 1 & nchar(b.x_trimmed) >= 6) {
            for (i in 1:(length(x_strip) - 1)) {
              x_strip_collapsed <- paste(x_strip[1:(length(x_strip) - i)], collapse = " ")
              found <- suppressMessages(suppressWarnings(exec_as.mo(x_strip_collapsed, clear_options = FALSE, allow_uncertain = FALSE)))
              if (!is.na(found)) {
                found_result <- found
                found <- microorganismsDT[mo == found, ..property][[1]]
                uncertainties <<- c(uncertainties,
                                    paste0("'", a.x_backup, "' >> ", microorganismsDT[mo == found_result[1L], fullname][[1]], " (", found_result[1L], ")"))
                return(found[1L])
              }
            }
          }

          # (5) try to strip off one element from start and check the remains ----
          x_strip <- a.x_backup %>% strsplit(" ") %>% unlist()
          if (length(x_strip) > 1 & nchar(b.x_trimmed) >= 6) {
            for (i in 2:(length(x_strip))) {
              x_strip_collapsed <- paste(x_strip[i:length(x_strip)], collapse = " ")
              found <- suppressMessages(suppressWarnings(exec_as.mo(x_strip_collapsed, clear_options = FALSE, allow_uncertain = FALSE)))
              if (!is.na(found)) {
                found_result <- found
                found <- microorganismsDT[mo == found, ..property][[1]]
                uncertainties <<- c(uncertainties,
                                    paste0("'", a.x_backup, "' >> ", microorganismsDT[mo == found_result[1L], fullname][[1]], " (", found_result[1L], ")"))
                return(found[1L])
              }
            }
          }

          # (6) not yet implemented taxonomic changes in Catalogue of Life ----
          found <- suppressMessages(suppressWarnings(exec_as.mo(TEMPORARY_TAXONOMY(b.x_trimmed), clear_options = FALSE, allow_uncertain = FALSE)))
          if (!is.na(found)) {
            found_result <- found
            found <- microorganismsDT[mo == found, ..property][[1]]
            warning(silver(paste0('Guessed with uncertainty: "',
                                  a.x_backup, '" >> ', italic(microorganismsDT[mo == found_result[1L], fullname][[1]]), " (", found_result[1L], ")")),
                    call. = FALSE, immediate. = FALSE)
            uncertainties <<- c(uncertainties,
                                paste0('"', a.x_backup, '" >> ', microorganismsDT[mo == found_result[1L], fullname][[1]], " (", found_result[1L], ")"))
            return(found[1L])
          }

          # didn't found in uncertain results too
          return(NA_character_)
        }

        x[i] <- uncertain_fn(x_backup[i], x_trimmed[i], x_withspaces_start_end[i], x_withspaces_start_only[i])
        if (!is.na(x[i])) {
          next
        }
      }

      # not found ----
      x[i] <- NA_character_
      failures <- c(failures, x_backup[i])
    }
  }

  # failures
  failures <- failures[!failures %in% c(NA, NULL, NaN)]
  if (length(failures) > 0) {
    options(mo_failures = sort(unique(failures)))
    plural <- c("value", "it")
    if (n_distinct(failures) > 1) {
      plural <- c("values", "them")
    }
    total_failures <- length(x_input[x_input %in% failures & !x_input %in% c(NA, NULL, NaN)])
    total_n <- length(x_input[!x_input %in% c(NA, NULL, NaN)])
    msg <- paste0("\n", n_distinct(failures), " unique ", plural[1],
                  " (^= ", percent(total_failures / total_n, round = 1, force_zero = TRUE),
                  ") could not be coerced to a valid MO code")
    if (n_distinct(failures) <= 10) {
      msg <- paste0(msg, ": ", paste('"', unique(failures), '"', sep = "", collapse = ', '))
    }
    msg <- paste0(msg,  ". Use mo_failures() to review ", plural[2], ".")
    warning(red(msg),
            call. = FALSE,
            immediate. = TRUE) # thus will always be shown, even if >= warnings
  }
  # uncertainties
  if (length(uncertainties) > 0) {
    options(mo_uncertainties = sort(unique(uncertainties)))
    plural <- c("value", "it")
    if (n_distinct(failures) > 1) {
      plural <- c("values", "them")
    }
    msg <- paste0("\nResults of ", n_distinct(uncertainties), " input ", plural[1],
                  " guessed with uncertainty. Use mo_uncertainties() to review ", plural[2], ".")
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
                                     "caprae", "carnosus", "cohnii", "condimenti",
                                     "devriesei", "epidermidis", "equorum",
                                     "fleurettii", "gallinarum", "haemolyticus",
                                     "hominis", "jettensis", "kloosii", "lentus",
                                     "lugdunensis", "massiliensis", "microti",
                                     "muscae", "nepalensis", "pasteuri", "petrasii",
                                     "pettenkoferi", "piscifermentans", "rostri",
                                     "saccharolyticus", "saprophyticus", "sciuri",
                                     "stepanovicii", "simulans", "succinus",
                                     "vitulinus", "warneri", "xylosus"), ..property][[1]]
    CoPS <- MOs_staph[species %in% c("simiae", "agnetis", "chromogenes",
                                     "delphini", "felis", "lutrae",
                                     "hyicus", "intermedius",
                                     "pseudintermedius", "pseudointermedius",
                                     "schleiferi"), ..property][[1]]
    x[x %in% CoNS] <- microorganismsDT[mo == 'B_STPHY_CNS', ..property][[1]][1L]
    x[x %in% CoPS] <- microorganismsDT[mo == 'B_STPHY_CPS', ..property][[1]][1L]
    if (Becker == "all") {
      x[x == microorganismsDT[mo == 'B_STPHY_AUR', ..property][[1]][1L]] <- microorganismsDT[mo == 'B_STPHY_CPS', ..property][[1]][1L]
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
  x_input_unique_nonempty <- unique(x_input[!is.na(x_input) & !is.null(x_input) & !identical(x_input, "")])

  # left join the found results to the original input values (x_input)
  df_found <- data.frame(input = as.character(x_input_unique_nonempty),
                         found = as.character(x),
                         stringsAsFactors = FALSE)
  df_input <- data.frame(input = as.character(x_input),
                         stringsAsFactors = FALSE)
  x <- df_input %>%
    left_join(df_found,
              by = "input") %>%
    pull(found)

  if (property == "mo") {
    class(x) <- "mo"
  }

  if (length(mo_renamed()) > 0) {
    if (has_color()) {
      notes <- getOption("mo_renamed")
    } else {
      notes <- mo_renamed()
    }
    notes <- sort(notes)
    for (i in 1:length(notes)) {
      base::message(blue(paste("Note:", notes[i])))
    }
  }

  x
}

TEMPORARY_TAXONOMY <- function(x) {
  x[x %like% 'Cutibacterium'] <- gsub('Cutibacterium', 'Propionibacterium', x[x %like% 'Cutibacterium'])
  x
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
  msg <- paste0(italic(name_old), ref_old, " was renamed ", italic(name_new), ref_new, mo)
  msg <- gsub("et al.", italic("et al."), msg)
  options(mo_renamed = sort(msg))
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
#' @export
#' @noRd
summary.mo <- function(object, ...) {
  # unique and top 1-3
  x <- object
  top_3 <- unname(top_freq(freq(x), 3))
  c("Class" = "mo",
    "<NA>" = length(x[is.na(x)]),
    "Unique" = dplyr::n_distinct(x[!is.na(x)]),
    "#1" = top_3[1],
    "#2" = top_3[2],
    "#3" = top_3[3])
}

#' @exportMethod as.data.frame.mo
#' @export
#' @noRd
as.data.frame.mo <- function (x, ...) {
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
#' @export
mo_uncertainties <- function() {
  getOption("mo_uncertainties")
}

#' @rdname as.mo
#' @export
mo_renamed <- function() {
  strip_style(gsub("was renamed", ">>", getOption("mo_renamed"), fixed = TRUE))
}
