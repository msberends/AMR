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

#' Determine multidrug-resistant organisms (MDRO)
#'
#' Determine which isolates are multidrug-resistant organisms (MDRO) according to country-specific guidelines.
#' @param tbl table with antibiotic columns, like e.g. \code{amox} and \code{amcl}
#' @param country country to determine guidelines. Should be a code from the \href{https://en.wikipedia.org/wiki/ISO_3166-1_alpha-2#Officially_assigned_code_elements}{list of ISO 3166-1 alpha-2 country codes}. Case-insensitive. Currently supported are \code{de} (Germany) and \code{nl} (the Netherlands).
#' @param col_bactid column name of the bacteria ID in \code{tbl} - values of this column should be present in \code{microorganisms$bactid}, see \code{\link{microorganisms}}
#' @param info print progress
#' @param aminoglycosides,quinolones,carbapenems character vector with column names of antibiotics
#' @param ceftazidime,piperacillin,trimethoprim_sulfa,penicillin,vancomycin column names of antibiotics
#' @param ... parameters that are passed on to \code{MDR}
#' @return Ordered factor with values \code{Positive}, \code{Unconfirmed}, \code{Negative}.
#' @rdname MDRO
#' @export
MDRO <- function(tbl,
                 country,
                 col_bactid = 'bactid',
                 info = TRUE,
                 aminoglycosides = c('gent', 'tobr', 'kana'),
                 quinolones = c('cipr', 'norf'),
                 carbapenems = c('imip', 'mero', 'erta'),
                 ceftazidime = 'cfta',
                 piperacillin = 'pita',
                 trimethoprim_sulfa = 'trsu',
                 penicillin = 'peni',
                 vancomycin = 'vanc') {

  # strip whitespaces
  country <- trimws(country)
  if (length(country) > 1) {
    stop('`country` must be a length one character string.', call. = FALSE)
  }
  if (!country %like% '^[a-z]{2}$') {
    stop('This is not a valid ISO 3166-1 alpha-2 country code: "', country, '". Please see ?MDRO.', call. = FALSE)
  }

  # create list and make country code case-independent
  guideline <- list(country = list(code = tolower(country)))

  # support per country
  if (guideline$country$code == 'de') {
    guideline$country$name <- 'Germany'
    guideline$name <- ''
    guideline$version <- ''
    guideline$source <- ''
  } else if (guideline$country$code == 'nl') {
    guideline$country$name <- 'The Netherlands'
    guideline$name <- 'WIP-Richtlijn BRMO'
    guideline$version <- 'Revision of December 2017'
    guideline$source <- 'https://www.rivm.nl/Documenten_en_publicaties/Professioneel_Praktisch/Richtlijnen/Infectieziekten/WIP_Richtlijnen/WIP_Richtlijnen/Ziekenhuizen/WIP_richtlijn_BRMO_Bijzonder_Resistente_Micro_Organismen_ZKH'
  # add here more countries like this:
  # } else if (country$code == 'AA') {
  #   country$name <- 'country name'
  } else {
    stop('This country code is currently unsupported: ', guideline$country$code, call. = FALSE)
  }

  # Console colours
  # source: http://www.tldp.org/HOWTO/Bash-Prompt-HOWTO/x329.html
  ANSI_red <- "\033[31m"
  ANSI_blue <- "\033[34m"
  ANSI_reset <- "\033[0m"

  if (info == TRUE) {
    cat("Determining Highly Resistant Microorganisms (MDRO), according to:\n",
        "Guideline: ", ANSI_red, guideline$name, ", ", guideline$version, ANSI_reset, "\n",
        "Country  : ", ANSI_red, guideline$country$name, ANSI_reset, "\n",
        "Source   : ", ANSI_blue, guideline$source, ANSI_reset, "\n",
        "\n", sep = "")
  }

  # join microorganisms
  tbl <- tbl %>% left_join_microorganisms(col_bactid)

  tbl$MDRO <- 1

  if (guideline$country$code == 'nl') {
    # BRMO; Bijzonder Resistente Micro-Organismen
    aminoglycosides <- aminoglycosides[aminoglycosides %in% colnames(tbl)]
    quinolones <- quinolones[quinolones %in% colnames(tbl)]
    carbapenems <- carbapenems[carbapenems %in% colnames(tbl)]
    if (!ceftazidime %in% colnames(tbl)) { ceftazidime <- NA }
    if (!piperacillin %in% colnames(tbl)) { piperacillin <- NA }
    if (!trimethoprim_sulfa %in% colnames(tbl)) { trimethoprim_sulfa <- NA }
    if (!penicillin %in% colnames(tbl)) { penicillin <- NA }
    if (!vancomycin %in% colnames(tbl)) { vancomycin <- NA }

    # Table 1
    tbl[which(
      tbl$family == 'Enterobacteriaceae'
      & rowSums(tbl[, aminoglycosides] == 'R', na.rm = TRUE) >= 1
      & rowSums(tbl[, quinolones] == 'R', na.rm = TRUE) >= 1
    ), 'MDRO'] <- 4
    tbl[which(
      tbl$family == 'Enterobacteriaceae'
      & rowSums(tbl[, carbapenems] == 'R', na.rm = TRUE) >= 1
    ), 'MDRO'] <- 3
    # rest is negative
    tbl[which(
      tbl$family == 'Enterobacteriaceae'
      & tbl$MDRO == 1
    ), 'MDRO'] <- 2

    # Table 2
    tbl[which(
      tbl$genus == 'Acinetobacter'
      & rowSums(tbl[, carbapenems] == 'R', na.rm = TRUE) >= 1
    ), 'MDRO'] <- 3
    tbl[which(
      tbl$genus == 'Acinetobacter'
      & rowSums(tbl[, aminoglycosides] == 'R', na.rm = TRUE) >= 1
      & rowSums(tbl[, quinolones] == 'R', na.rm = TRUE) >= 1
    ), 'MDRO'] <- 4
    # rest of Acinetobacter is negative
    tbl[which(
      tbl$genus == 'Acinetobacter'
      & tbl$MDRO == 1
    ), 'MDRO'] <- 2

    tbl[which(
      tbl$fullname %like% 'Stenotrophomonas maltophilia'
      & tbl[, trimethoprim_sulfa] == 'R'
    ), 'MDRO'] <- 4
    # rest of Stenotrophomonas is negative
    tbl[which(
      tbl$fullname %like% 'Stenotrophomonas maltophilia'
      & tbl$MDRO == 1
    ), 'MDRO'] <- 2

    tbl[which(
      tbl$fullname %like% 'Pseudomonas aeruginosa'
      & sum(rowSums(tbl[, carbapenems] == 'R', na.rm = TRUE) >= 1,
           rowSums(tbl[, aminoglycosides] == 'R', na.rm = TRUE) >= 1,
           rowSums(tbl[, quinolones] == 'R', na.rm = TRUE) >= 1,
           tbl[, ceftazidime] == 'R',
           tbl[, piperacillin] == 'R') >= 3
    ), 'MDRO'] <- 4
    # rest of Pseudomonas is negative
    tbl[which(
      tbl$fullname %like% 'Pseudomonas aeruginosa'
      & tbl$MDRO == 1
    ), 'MDRO'] <- 2

    # Table 3
    tbl[which(
      tbl$fullname %like% 'Streptococcus pneumoniae'
      & tbl[, penicillin] == 'R'
    ), 'MDRO'] <- 4
    tbl[which(
      tbl$fullname %like% 'Streptococcus pneumoniae'
      & tbl[, vancomycin] == 'R'
    ), 'MDRO'] <- 4
    # rest of Streptococcus pneumoniae is negative
    tbl[which(
      tbl$fullname %like% 'Streptococcus pneumoniae'
      & tbl$MDRO == 1
    ), 'MDRO'] <- 2

    tbl[which(
      tbl$fullname %like% 'Enterococcus faecium'
      &  rowSums(tbl[, c(penicillin, vancomycin)] == 'R', na.rm = TRUE) >= 1
    ), 'MDRO'] <- 4
    # rest of Enterococcus faecium is negative
    tbl[which(
      tbl$fullname %like% 'Enterococcus faecium'
      & tbl$MDRO == 1
    ), 'MDRO'] <- 2
  }

  factor(x =  tbl$MDRO,
         levels = c(1:4),
         labels = c('Unknown', 'Negative', 'Unconfirmed', 'Positive'),
         ordered = TRUE)
}

#' @rdname MDRO
#' @export
BRMO <- function(tbl, country = "nl", ...) {
  MDRO(tbl = tbl, country = country, ...)
}

#' @rdname MDRO
#' @export
MRGN <- function(tbl, country = "de", ...) {
  MDRO(tbl = tbl, country = country, ...)
}
