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
#' @param allow_uncertain a logical to indicate whether empty results should be checked for only a part of the input string. When results are found, a warning will be given about the uncertainty and the result.
#' @param reference_df a \code{data.frame} to use for extra reference when translating \code{x} to a valid \code{mo}. The first column can be any microbial name, code or ID (used in your analysis or organisation), the second column must be a valid \code{mo} as found in the \code{\link{microorganisms}} data set.
#' @rdname as.mo
#' @aliases mo
#' @keywords mo Becker becker Lancefield lancefield guess
#' @details
#' A microbial ID from this package (class: \code{mo}) typically looks like these examples:\cr
#' \preformatted{
#'   Code                Full name
#'   ---------------     --------------------------------------
#'   B_KLBSL             Klebsiella
#'   B_KLBSL_PNE         Klebsiella pneumoniae
#'   B_KLBSL_PNE_RHI     Klebsiella pneumoniae rhinoscleromatis
#'   |   |    |   |
#'   |   |    |   |
#'   |   |    |    ----> subspecies, a 3-4 letter acronym
#'   |   |     ----> species, a 3-4 letter acronym
#'   |    ----> genus, a 5-7 letter acronym, mostly without vowels
#'    ----> taxonomic kingdom, either B (Bacteria), F (Fungi) or P (Protozoa)
#' }
#'
#' Use the \code{\link{mo_property}} functions to get properties based on the returned code, see Examples.
#'
#' This function uses Artificial Intelligence (AI) to help getting fast and logical results. It tries to find matches in this order:
#' \itemize{
#'   \item{Taxonomic kingdom: it first searches in bacteria, then fungi, then protozoa}
#'   \item{Human pathogenic prevalence: it first searches in more prevalent microorganisms, then less prevalent ones}
#'   \item{Valid MO codes and full names: it first searches in already valid MO code and genus/species combinations}
#'   \item{Breakdown of input values: from here it starts to breakdown input values to find possible matches}
#' }
#'
#' A couple of effects because of these rules
#' \itemize{
#'   \item{\code{"E. coli"} will return the ID of \emph{Escherichia coli} and not \emph{Entamoeba coli}, although the latter would alphabetically come first}
#'   \item{\code{"H. influenzae"} will return the ID of \emph{Haemophilus influenzae} and not \emph{Haematobacter influenzae} for the same reason}
#'   \item{Something like \code{"p aer"} will return the ID of \emph{Pseudomonas aeruginosa} and not \emph{Pasteurella aerogenes}}
#'   \item{Something like \code{"stau"} or \code{"S aur"} will return the ID of \emph{Staphylococcus aureus} and not \emph{Staphylococcus auricularis}}
#' }
#' This means that looking up human pathogenic microorganisms takes less time than looking up human \strong{non}-pathogenic microorganisms.
#'
#' \code{guess_mo} is an alias of \code{as.mo}.
#' @section ITIS:
#' \if{html}{\figure{itis_logo.jpg}{options: height=60px style=margin-bottom:5px} \cr}
#' This package contains the \strong{complete microbial taxonomic data} (with all eight taxonomic ranks - from kingdom to subspecies) from the publicly available Integrated Taxonomic Information System (ITIS, \url{https://www.itis.gov}).
#'
#' All (sub)species from \strong{the taxonomic kingdoms Bacteria, Fungi and Protozoa are included in this package}, as well as all previously accepted names known to ITIS. Furthermore, the responsible authors and year of publication are available. This allows users to use authoritative taxonomic information for their data analysis on any microorganism, not only human pathogens. It also helps to quickly determine the Gram stain of bacteria, since all bacteria are classified into subkingdom Negibacteria or Posibacteria.
#'
#' ITIS is a partnership of U.S., Canadian, and Mexican agencies and taxonomic specialists [3].
#'
#  (source as a section, so it can be inherited by other man pages)
#' @section Source:
#' [1] Becker K \emph{et al.} \strong{Coagulase-Negative Staphylococci}. 2014. Clin Microbiol Rev. 27(4): 870–926. \url{https://dx.doi.org/10.1128/CMR.00109-13}
#'
#' [2] Lancefield RC \strong{A serological differentiation of human and other groups of hemolytic streptococci}. 1933. J Exp Med. 57(4): 571–95. \url{https://dx.doi.org/10.1084/jem.57.4.571}
#'
#' [3] Integrated Taxonomic Information System (ITIS). Retrieved September 2018. \url{http://www.itis.gov}
#' @export
#' @return Character (vector) with class \code{"mo"}. Unknown values will return \code{NA}.
#' @seealso \code{\link{microorganisms}} for the \code{data.frame} with ITIS content that is being used to determine ID's. \cr
#' The \code{\link{mo_property}} functions (like \code{\link{mo_genus}}, \code{\link{mo_gramstain}}) to get properties based on the returned code.
#' @examples
#' # These examples all return "B_STPHY_AUR", the ID of S. aureus:
#' as.mo("stau")
#' as.mo("STAU")
#' as.mo("staaur")
#' as.mo("S. aureus")
#' as.mo("S aureus")
#' as.mo("Staphylococcus aureus")
#' as.mo("MRSA") # Methicillin Resistant S. aureus
#' as.mo("VISA") # Vancomycin Intermediate S. aureus
#' as.mo("VRSA") # Vancomycin Resistant S. aureus
#' as.mo(369)    # Search on TSN (Taxonomic Serial Number), a unique identifier
#'               # for the Integrated Taxonomic Information System (ITIS)
#'
#' as.mo("Streptococcus group A")
#' as.mo("GAS") # Group A Streptococci
#' as.mo("GBS") # Group B Streptococci
#'
#' # guess_mo is an alias of as.mo and works the same
#' guess_mo("S. epidermidis")                 # will remain species: B_STPHY_EPI
#' guess_mo("S. epidermidis", Becker = TRUE)  # will not remain species: B_STPHY_CNS
#'
#' guess_mo("S. pyogenes")                    # will remain species: B_STRPTC_PYO
#' guess_mo("S. pyogenes", Lancefield = TRUE) # will not remain species: B_STRPTC_GRA
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
#'   guess_mo()
#'
#' # and can even contain 2 columns, which is convenient for genus/species combinations:
#' df$mo <- df %>%
#'   select(genus, species) %>%
#'   guess_mo()
#'
#' # same result:
#' df <- df %>%
#'   mutate(mo = guess_mo(paste(genus, species)))
#' }
as.mo <- function(x, Becker = FALSE, Lancefield = FALSE, allow_uncertain = FALSE, reference_df = NULL) {
  structure(mo_validate(x = x, property = "mo",
                        Becker = Becker, Lancefield = Lancefield,
                        allow_uncertain = allow_uncertain, reference_df = reference_df),
            class = "mo")
}

#' @rdname as.mo
#' @export
is.mo <- function(x) {
  # bactid for older releases
  # remove when is.bactid will be removed
  identical(class(x), "mo") | identical(class(x), "bactid")
}

#' @rdname as.mo
#' @export
guess_mo <- as.mo

#' @importFrom dplyr %>% pull left_join
#' @importFrom data.table data.table as.data.table setkey
exec_as.mo <- function(x, Becker = FALSE, Lancefield = FALSE, allow_uncertain = FALSE, reference_df = NULL, property = "mo") {

  # These data.tables are available as data sets when the AMR package is loaded:
  #   microorganismsDT        # this one is sorted by kingdom (B<F<P), prevalence, TSN
  #   microorganisms.prevDT   # same as microorganismsDT, but with prevalence != 9999
  #   microorganisms.unprevDT # same as microorganismsDT, but with prevalence == 9999
  #   microorganisms.oldDT    # old taxonomic names, sorted by name (genus+species), TSN

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

  failures <- character(0)
  x_input <- x
  # only check the uniques, which is way faster
  x <- unique(x)
  # remove empty values (to later fill them in again)
  x <- x[!is.na(x) & !is.null(x) & !identical(x, "")]

  # defined df to check for
  if (!is.null(reference_df)) {
    if (!is.data.frame(reference_df) | NCOL(reference_df) < 2) {
      stop('`reference_df` must be a data.frame with at least two columns.', call. = FALSE)
    }
    # remove factors, just keep characters
    suppressWarnings(
      reference_df[] <- lapply(reference_df, as.character)
    )
  }

  if (all(x %in% microorganismsDT[["mo"]])) {
    # existing mo codes when not looking for property "mo", like mo_genus("B_ESCHR_COL")
    x <- microorganismsDT[data.table(mo = x), on = "mo", ..property][[1]]
  } else if (!is.null(reference_df)
             & all(x %in% reference_df[, 1])
             & all(reference_df[, 2] %in% microorganismsDT[["mo"]])) {
    # manually defined reference
    colnames(reference_df)[1] <- "x"
    colnames(reference_df)[2] <- "mo"
    suppressWarnings(
      x <- data.frame(x = x, stringsAsFactors = FALSE) %>%
        left_join(reference_df, by = "x") %>%
        left_join(AMR::microorganisms, by = "mo") %>%
        pull(property)
    )
  } else if (all(toupper(x) %in% AMR::microorganisms.certe[, "certe"])) {
    # old Certe codes
    y <- as.data.table(AMR::microorganisms.certe)[data.table(certe = toupper(x)), on = "certe", ]
    x <- microorganismsDT[data.table(mo = y[["mo"]]), on = "mo", ..property][[1]]

  } else if (!all(x %in% microorganismsDT[[property]])) {

    x_backup <- trimws(x, which = "both")
    x_species <- paste(x_backup, "species")
    # translate to English for supported languages of mo_property
    x <- gsub("(Gruppe|gruppe|groep|grupo|gruppo|groupe)", "group", x)
    # remove 'empty' genus and species values
    x <- gsub("(no MO)", "", x, fixed = TRUE)
    # remove dots and other non-text in case of "E. coli" except spaces
    x <- gsub("[^a-zA-Z0-9/ \\-]+", "", x)
    # but spaces before and after should be omitted
    x <- trimws(x, which = "both")
    x_trimmed <- x
    x_trimmed_species <- paste(x_trimmed, "species")
    # replace space by regex sign
    x_withspaces <- gsub(" ", ".* ", x, fixed = TRUE)
    x <- gsub(" ", ".*", x, fixed = TRUE)
    # add start en stop regex
    x <- paste0('^', x, '$')
    x_withspaces_start <- paste0('^', x_withspaces)
    x_withspaces <- paste0('^', x_withspaces, '$')

    # cat(paste0('x                  "', x, '"\n'))
    # cat(paste0('x_species          "', x_species, '"\n'))
    # cat(paste0('x_withspaces_start "', x_withspaces_start, '"\n'))
    # cat(paste0('x_withspaces       "', x_withspaces, '"\n'))
    # cat(paste0('x_backup           "', x_backup, '"\n'))
    # cat(paste0('x_trimmed          "', x_trimmed, '"\n'))
    # cat(paste0('x_trimmed_species  "', x_trimmed_species, '"\n'))

    for (i in 1:length(x)) {
      if (identical(x_trimmed[i], "")) {
        # empty values
        x[i] <- NA_character_
        next
      }
      if (nchar(x_trimmed[i]) < 3) {
        # fewer than 3 chars, add as failure
        x[i] <- NA_character_
        failures <- c(failures, x_backup[i])
        next
      }

      # translate known trivial abbreviations to genus + species ----
      if (!is.na(x_trimmed[i])) {
        if (toupper(x_trimmed[i]) == 'MRSA'
            | toupper(x_trimmed[i]) == 'MSSA'
            | toupper(x_trimmed[i]) == 'VISA'
            | toupper(x_trimmed[i]) == 'VRSA') {
          x[i] <- microorganismsDT[mo == 'B_STPHY_AUR', ..property][[1]][1L]
          next
        }
        if (toupper(x_trimmed[i]) == 'MRSE'
            | toupper(x_trimmed[i]) == 'MSSE') {
          x[i] <- microorganismsDT[mo == 'B_STPHY_EPI', ..property][[1]][1L]
          next
        }
        if (toupper(x_trimmed[i]) == 'VRE') {
          x[i] <- microorganismsDT[mo == 'B_ENTRC', ..property][[1]][1L]
          next
        }
        if (toupper(x_trimmed[i]) == 'MRPA') {
          # multi resistant P. aeruginosa
          x[i] <- microorganismsDT[mo == 'B_PDMNS_AER', ..property][[1]][1L]
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
          x[i] <- microorganismsDT[mo == 'B_STRPTC_PNE', ..property][[1]][1L]
          next
        }
        if (toupper(x_trimmed[i]) %like% '^G[ABCDFGHK]S$') {
          x[i] <- microorganismsDT[mo == gsub("G([ABCDFGHK])S", "B_STRPTC_GR\\1", x_trimmed[i]), ..property][[1]][1L]
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
        if (tolower(x[i]) %like% '^gram[ -]+nega.*'
            | tolower(x_trimmed[i]) %like% '^gram[ -]+nega.*') {
          # coerce S. coagulase positive
          x[i] <- microorganismsDT[mo == 'B_GRAMN', ..property][[1]][1L]
          next
        }
        if (tolower(x[i]) %like% '^gram[ -]+posi.*'
            | tolower(x_trimmed[i]) %like% '^gram[ -]+posi.*') {
          # coerce S. coagulase positive
          x[i] <- microorganismsDT[mo == 'B_GRAMP', ..property][[1]][1L]
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
        if (nchar(x_trimmed[i]) > 4) {
          # not when abbr is esco, stau, klpn, etc.
          found <- microorganismsDT[tolower(fullname) %like% gsub(" ", ".*", x_trimmed_species[i], fixed = TRUE), ..property][[1]]
          if (length(found) > 0) {
            x[i] <- found[1L]
            next
          }
        }
      }

      # TRY OTHER SOURCES ----
      if (x_backup[i] %in% AMR::microorganisms.certe$certe) {
        x[i] <- microorganismsDT[mo == AMR::microorganisms.certe[AMR::microorganisms.certe[, 1] == x_backup[i], 2], ..property][[1]][1L]
        # x[i] <- exec_as.mo(x = AMR::microorganisms.certe[AMR::microorganisms.certe$certe == x_backup[i], "mo"],
        #                    property = property)
        # next
      }
      if (x_backup[i] %in% AMR::microorganisms.umcg[, 1]) {
        mo_umcg <- AMR::microorganisms.umcg[AMR::microorganisms.umcg[, 1] == x_backup[i], 2]
        mo_found <- AMR::microorganisms.certe[AMR::microorganisms.certe[, 1] == mo_umcg, 2]
        if (length(mo_found) == 0) {
          # not found
          x[i] <- NA_character_
          failures <- c(failures, x_backup[i])
        } else {
          x[i] <- microorganismsDT[mo == mo_found, ..property][[1]][1L]
        }
        next
      }
      if (x_backup[i] %in% reference_df[, 1]) {
        ref_mo <- reference_df[reference_df[, 1] == x_backup[i], 2]
        if (ref_mo %in% microorganismsDT[, mo]) {
          x[i] <- microorganismsDT[mo == ref_mo, ..property][[1]][1L]
          next
        } else {
          warning("Value '", x_backup[i], "' was found in reference_df, but '", ref_mo, "' is not a valid MO code.", call. = FALSE)
        }
      }

      # TRY FIRST THOUSAND MOST PREVALENT IN HUMAN INFECTIONS ----
      found <- microorganisms.prevDT[tolower(fullname) %in% tolower(c(x_backup[i], x_trimmed[i])), ..property][[1]]
      # most probable: is exact match in fullname
      if (length(found) > 0) {
        x[i] <- found[1L]
        next
      }
      found <- microorganisms.prevDT[tsn == x_trimmed[i], ..property][[1]]
      # is a valid TSN
      if (length(found) > 0) {
        x[i] <- found[1L]
        next
      }
      found <- microorganisms.prevDT[mo == toupper(x_backup[i]), ..property][[1]]
      # is a valid mo
      if (length(found) > 0) {
        x[i] <- found[1L]
        next
      }

      # try any match keeping spaces ----
      found <- microorganisms.prevDT[fullname %like% x_withspaces[i], ..property][[1]]
      if (length(found) > 0) {
        x[i] <- found[1L]
        next
      }

      # try any match keeping spaces, not ending with $ ----
      found <- microorganisms.prevDT[fullname %like% x_withspaces_start[i], ..property][[1]]
      if (length(found) > 0) {
        x[i] <- found[1L]
        next
      }

      # try any match diregarding spaces ----
      found <- microorganisms.prevDT[fullname %like% x[i], ..property][[1]]
      if (length(found) > 0) {
        x[i] <- found[1L]
        next
      }


      # try splitting of characters in the middle and then find ID ----
      # only when text length is 6 or lower
      # like esco = E. coli, klpn = K. pneumoniae, stau = S. aureus, staaur = S. aureus
      if (nchar(x_trimmed[i]) <= 6) {
        x_length <- nchar(x_trimmed[i])
        x[i] <- paste0(x_trimmed[i] %>% substr(1, x_length / 2),
                       '.* ',
                       x_trimmed[i] %>% substr((x_length / 2) + 1, x_length))
        found <- microorganisms.prevDT[fullname %like% paste0('^', x[i]), ..property][[1]]
        if (length(found) > 0) {
          x[i] <- found[1L]
          next
        }
      }

      # try fullname without start and stop regex, to also find subspecies ----
      # like "K. pneu rhino" -> "Klebsiella pneumoniae (rhinoscleromatis)" = KLEPNERH
      found <- microorganisms.prevDT[fullname %like% x_withspaces_start[i], ..property][[1]]
      if (length(found) > 0) {
        x[i] <- found[1L]
        next
      }

      # THEN TRY ALL OTHERS ----
      found <- microorganisms.unprevDT[tolower(fullname) == tolower(x_backup[i]), ..property][[1]]
      # most probable: is exact match in fullname
      if (length(found) > 0) {
        x[i] <- found[1L]
        next
      }
      found <- microorganisms.unprevDT[tolower(fullname) == tolower(x_trimmed[i]), ..property][[1]]
      # most probable: is exact match in fullname
      if (length(found) > 0) {
        x[i] <- found[1L]
        next
      }
      found <- microorganisms.unprevDT[tsn == x_trimmed[i], ..property][[1]]
      # is a valid TSN
      if (length(found) > 0) {
        x[i] <- found[1L]
        next
      }
      found <- microorganisms.unprevDT[mo == toupper(x_backup[i]), ..property][[1]]
      # is a valid mo
      if (length(found) > 0) {
        x[i] <- found[1L]
        next
      }

      # try any match keeping spaces ----
      found <- microorganisms.unprevDT[fullname %like% x_withspaces[i], ..property][[1]]
      if (length(found) > 0) {
        x[i] <- found[1L]
        next
      }

      # try any match keeping spaces, not ending with $ ----
      found <- microorganisms.unprevDT[fullname %like% x_withspaces_start[i], ..property][[1]]
      if (length(found) > 0) {
        x[i] <- found[1L]
        next
      }

      # try any match diregarding spaces ----
      found <- microorganisms.unprevDT[fullname %like% x[i], ..property][[1]]
      if (length(found) > 0) {
        x[i] <- found[1L]
        next
      }

      # try splitting of characters in the middle and then find ID ----
      # only when text length is 6 or lower
      # like esco = E. coli, klpn = K. pneumoniae, stau = S. aureus, staaur = S. aureus
      if (nchar(x_trimmed[i]) <= 6) {
        x_length <- nchar(x_trimmed[i])
        x[i] <- paste0(x_trimmed[i] %>% substr(1, x_length / 2),
                       '.* ',
                       x_trimmed[i] %>% substr((x_length / 2) + 1, x_length))
        found <- microorganisms.unprevDT[fullname %like% paste0('^', x[i]), ..property][[1]]
        if (length(found) > 0) {
          x[i] <- found[1L]
          next
        }
      }

      # try fullname without start and stop regex, to also find subspecies ----
      # like "K. pneu rhino" -> "Klebsiella pneumoniae (rhinoscleromatis)" = KLEPNERH
      found <- microorganisms.unprevDT[fullname %like% x_withspaces_start[i], ..property][[1]]
      if (length(found) > 0) {
        x[i] <- found[1L]
        next
      }

      # MISCELLANEOUS ----

      # look for old taxonomic names ----
      found <- microorganisms.oldDT[tolower(name) == tolower(x_backup[i])
                                    | tsn == x_trimmed[i]
                                    | name %like% x_withspaces[i],]
      if (NROW(found) > 0) {
        # when property is "ref" (which is the case in mo_ref, mo_authors and mo_year), return the old value, so:
        # mo_ref("Chlamydia psittaci) = "Page, 1968" (with warning)
        # mo_ref("Chlamydophila psittaci) = "Everett et al., 1999"
        if (property == "ref") {
          x[i] <- found[1, ref]
        } else {
          x[i] <- microorganismsDT[tsn == found[1, tsn_new], ..property][[1]]
        }
        renamed_note(name_old = found[1, name],
                     name_new = microorganismsDT[tsn == found[1, tsn_new], fullname],
                     ref_old = found[1, ref],
                     ref_new = microorganismsDT[tsn == found[1, tsn_new], ref])
        next
      }

      # check for uncertain results ----
      if (allow_uncertain == TRUE) {
        # (1) look again for old taxonomic names, now for G. species ----
        found <- microorganisms.oldDT[name %like% x_withspaces[i]
                                      | name %like% x_withspaces_start[i]
                                      | name %like% x[i],]
        if (NROW(found) > 0) {
          if (property == "ref") {
            x[i] <- found[1, ref]
          } else {
            x[i] <- microorganismsDT[tsn == found[1, tsn_new], ..property][[1]]
          }
          warning("Uncertain interpretation: '",
                  x_backup[i], "' -> '", found[1, name], "'",
                  call. = FALSE, immediate. = TRUE)
          renamed_note(name_old = found[1, name],
                       name_new = microorganismsDT[tsn == found[1, tsn_new], fullname],
                       ref_old = found[1, ref],
                       ref_new = microorganismsDT[tsn == found[1, tsn_new], ref])
          next
        }

        # (2) try to strip off one element and check the remains
        x_strip <- x_backup[i] %>% strsplit(" ") %>% unlist()
        x_strip <- x_strip[1:length(x_strip) - 1]
        x[i] <- suppressWarnings(suppressMessages(as.mo(x_strip)))
        if (!is.na(x[i])) {
          warning("Uncertain interpretation: '",
                  x_backup[i], "' -> '", microorganismsDT[mo == x[i], fullname], "' (", x[i], ")",
                  call. = FALSE, immediate. = TRUE)
          next
        }
      }

      # not found ----
      x[i] <- NA_character_
      failures <- c(failures, x_backup[i])

    }
  }

  failures <- failures[!failures %in% c(NA, NULL, NaN)]
  if (length(failures) > 0) {
    warning("These ", length(failures) , " values could not be coerced to a valid MO code: ",
            paste('"', unique(failures), '"', sep = "", collapse = ', '),
            ".",
            call. = FALSE)
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
    x[x == microorganismsDT[mo == 'B_STRPTC_PYO', ..property][[1]][1L]] <- microorganismsDT[mo == 'B_STRPTC_GRA', ..property][[1]][1L]
    # group B - S. agalactiae
    x[x == microorganismsDT[mo == 'B_STRPTC_AGA', ..property][[1]][1L]] <- microorganismsDT[mo == 'B_STRPTC_GRB', ..property][[1]][1L]
    # group C
    S_groupC <- microorganismsDT %>% filter(genus == "Streptococcus",
                                            species %in% c("equisimilis", "equi",
                                                           "zooepidemicus", "dysgalactiae")) %>%
      pull(property)
    x[x %in% S_groupC] <- microorganismsDT[mo == 'B_STRPTC_GRC', ..property][[1]][1L]
    if (Lancefield == "all") {
      # all Enterococci
      x[x %like% "^(Enterococcus|B_ENTRC)"] <- microorganismsDT[mo == 'B_STRPTC_GRD', ..property][[1]][1L]
    }
    # group F - S. anginosus
    x[x == microorganismsDT[mo == 'B_STRPTC_ANG', ..property][[1]][1L]] <- microorganismsDT[mo == 'B_STRPTC_GRF', ..property][[1]][1L]
    # group H - S. sanguinis
    x[x == microorganismsDT[mo == 'B_STRPTC_SAN', ..property][[1]][1L]] <- microorganismsDT[mo == 'B_STRPTC_GRH', ..property][[1]][1L]
    # group K - S. salivarius
    x[x == microorganismsDT[mo == 'B_STRPTC_SAL', ..property][[1]][1L]] <- microorganismsDT[mo == 'B_STRPTC_GRK', ..property][[1]][1L]
  }

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
  } else if (property == "tsn") {
    x <- as.integer(x)
  }

  x
}

renamed_note <- function(name_old, name_new, ref_old = "", ref_new = "") {
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
  base::message(paste0("Note: '", name_old, "'", ref_old, " was renamed '", name_new, "'", ref_new))
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
