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
#' @param country country code to determine guidelines. EUCAST rules will be used when left empty, see Details. Should be or a code from the \href{https://en.wikipedia.org/wiki/ISO_3166-1_alpha-2#Officially_assigned_code_elements}{list of ISO 3166-1 alpha-2 country codes}. Case-insensitive. Currently supported are \code{de} (Germany) and \code{nl} (the Netherlands).
#' @param col_bactid column name of the bacteria ID in \code{tbl} - values of this column should be present in \code{microorganisms$bactid}, see \code{\link{microorganisms}}
#' @param info print progress
#' @param amcl,amik,amox,ampi,azit,aztr,cefa,cfra,cfep,cfot,cfox,cfta,cftr,cfur,chlo,cipr,clar,clin,clox,coli,czol,dapt,doxy,erta,eryt,fosf,fusi,gent,imip,kana,levo,linc,line,mero,metr,mino,moxi,nali,neom,neti,nitr,novo,norf,oflo,peni,pita,poly,qida,rifa,roxi,siso,teic,tetr,tica,tige,tobr,trim,trsu,vanc column names of antibiotics. column names of antibiotics
#' @param ... parameters that are passed on to methods
#' @details When \code{country} will be left blank, guidelines will be taken from EUCAST Expert Rules Version 3.1 "Intrinsic Resistance and Exceptional Phenotypes Tables" (\url{http://www.eucast.org/fileadmin/src/media/PDFs/EUCAST_files/Expert_Rules/Expert_rules_intrinsic_exceptional_V3.1.pdf}).
#' @return Ordered factor with levels \code{Unknown < Negative < Unconfirmed < Positive}.
#' @rdname MDRO
#' @export
#' @examples
#' library(dplyr)
#'
#' septic_patients %>%
#'   mutate(EUCAST = MDRO(.),
#'          BRMO = MDRO(., "nl"))
MDRO <- function(tbl,
                 country = NULL,
                 col_bactid = 'bactid',
                 info = TRUE,
                 amcl = 'amcl',
                 amik = 'amik',
                 amox = 'amox',
                 ampi = 'ampi',
                 azit = 'azit',
                 aztr = 'aztr',
                 cefa = 'cefa',
                 cfra = 'cfra',
                 cfep = 'cfep',
                 cfot = 'cfot',
                 cfox = 'cfox',
                 cfta = 'cfta',
                 cftr = 'cftr',
                 cfur = 'cfur',
                 chlo = 'chlo',
                 cipr = 'cipr',
                 clar = 'clar',
                 clin = 'clin',
                 clox = 'clox',
                 coli = 'coli',
                 czol = 'czol',
                 dapt = 'dapt',
                 doxy = 'doxy',
                 erta = 'erta',
                 eryt = 'eryt',
                 fosf = 'fosf',
                 fusi = 'fusi',
                 gent = 'gent',
                 imip = 'imip',
                 kana = 'kana',
                 levo = 'levo',
                 linc = 'linc',
                 line = 'line',
                 mero = 'mero',
                 metr = 'metr',
                 mino = 'mino',
                 moxi = 'moxi',
                 nali = 'nali',
                 neom = 'neom',
                 neti = 'neti',
                 nitr = 'nitr',
                 novo = 'novo',
                 norf = 'norf',
                 oflo = 'oflo',
                 peni = 'peni',
                 pita = 'pita',
                 poly = 'poly',
                 qida = 'qida',
                 rifa = 'rifa',
                 roxi = 'roxi',
                 siso = 'siso',
                 teic = 'teic',
                 tetr = 'tetr',
                 tica = 'tica',
                 tige = 'tige',
                 tobr = 'tobr',
                 trim = 'trim',
                 trsu = 'trsu',
                 vanc = 'vanc') {

  if (!col_bactid %in% colnames(tbl)) {
    stop('Column ', col_bactid, ' not found.', call. = FALSE)
  }

  # strip whitespaces
  if (length(country) > 1) {
    stop('`country` must be a length one character string.', call. = FALSE)
  }

  if (is.null(country)) {
    country <- 'EUCAST'
  }
  country <- trimws(country)
  if (country != 'EUCAST' & !country %like% '^[a-z]{2}$') {
    stop('This is not a valid ISO 3166-1 alpha-2 country code: "', country, '". Please see ?MDRO.', call. = FALSE)
  }

  # create list and make country code case-independent
  guideline <- list(country = list(code = tolower(country)))

  if (guideline$country$code == 'eucast') {
    guideline$country$name <- '(European guidelines)'
    guideline$name <- 'EUCAST Expert Rules, "Intrinsic Resistance and Exceptional Phenotypes Tables"'
    guideline$version <- 'Version 3.1'
    guideline$source <- 'http://www.eucast.org/fileadmin/src/media/PDFs/EUCAST_files/Expert_Rules/Expert_rules_intrinsic_exceptional_V3.1.pdf'
    # support per country:
  } else if (guideline$country$code == 'de') {
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
    cat("Determining multidrug-resistant organisms (MDRO), according to:\n",
        "Guideline: ", ANSI_red, guideline$name, ", ", guideline$version, ANSI_reset, "\n",
        "Country  : ", ANSI_red, guideline$country$name, ANSI_reset, "\n",
        "Source   : ", ANSI_blue, guideline$source, ANSI_reset, "\n",
        "\n", sep = "")
  }

  # check columns
  col.list <- c(amcl, amik, amox, ampi, azit, aztr, cefa, cfra, cfep, cfot,
                cfox, cfta, cftr, cfur, chlo, cipr, clar, clin, clox, coli,
                czol, dapt, doxy, erta, eryt, fosf, fusi, gent, imip, kana,
                levo, linc, line, mero, metr, mino, moxi, nali, neom, neti, nitr,
                novo, norf, oflo, peni, pita, poly, qida, rifa, roxi, siso,
                teic, tetr, tica, tige, tobr, trim, trsu, vanc)
  col.list <- check_available_columns(tbl = tbl, col.list = col.list, info = info)
  amcl <- col.list[amcl]
  amik <- col.list[amik]
  amox <- col.list[amox]
  ampi <- col.list[ampi]
  azit <- col.list[azit]
  aztr <- col.list[aztr]
  cefa <- col.list[cefa]
  cfra <- col.list[cfra]
  cfep <- col.list[cfep]
  cfot <- col.list[cfot]
  cfox <- col.list[cfox]
  cfta <- col.list[cfta]
  cftr <- col.list[cftr]
  cfur <- col.list[cfur]
  chlo <- col.list[chlo]
  cipr <- col.list[cipr]
  clar <- col.list[clar]
  clin <- col.list[clin]
  clox <- col.list[clox]
  coli <- col.list[coli]
  czol <- col.list[czol]
  dapt <- col.list[dapt]
  doxy <- col.list[doxy]
  erta <- col.list[erta]
  eryt <- col.list[eryt]
  fosf <- col.list[fosf]
  fusi <- col.list[fusi]
  gent <- col.list[gent]
  imip <- col.list[imip]
  kana <- col.list[kana]
  levo <- col.list[levo]
  linc <- col.list[linc]
  line <- col.list[line]
  mero <- col.list[mero]
  metr <- col.list[metr]
  mino <- col.list[mino]
  moxi <- col.list[moxi]
  nali <- col.list[nali]
  neom <- col.list[neom]
  neti <- col.list[neti]
  nitr <- col.list[nitr]
  novo <- col.list[novo]
  norf <- col.list[norf]
  oflo <- col.list[oflo]
  peni <- col.list[peni]
  pita <- col.list[pita]
  poly <- col.list[poly]
  qida <- col.list[qida]
  rifa <- col.list[rifa]
  roxi <- col.list[roxi]
  siso <- col.list[siso]
  teic <- col.list[teic]
  tetr <- col.list[tetr]
  tica <- col.list[tica]
  tige <- col.list[tige]
  tobr <- col.list[tobr]
  trim <- col.list[trim]
  trsu <- col.list[trsu]
  vanc <- col.list[vanc]

  # antibiotic classes
  aminoglycosides <- c(tobr, gent, kana, neom, neti, siso)
  tetracyclines <- c(doxy, mino, tetr) # since EUCAST v3.1 tige(cycline) is set apart
  polymyxins <- c(poly, coli)
  macrolides <- c(eryt, azit, roxi, clar) # since EUCAST v3.1 clin(damycin) is set apart
  glycopeptides <- c(vanc, teic)
  streptogramins <- qida # should officially also be pristinamycin and quinupristin/dalfopristin
  cephalosporins <- c(cfep, cfot, cfox, cfra, cfta, cftr, cfur, czol)
  cephalosporins_3rd <- c(cfot, cftr, cfta)
  carbapenems <- c(erta, imip, mero)
  aminopenicillins <- c(ampi, amox)
  ureidopenicillins <- pita # should officially also be azlo and mezlo
  fluoroquinolones <- c(oflo, cipr, norf, levo, moxi)

  # helper function for editing the table
  trans_tbl <- function(to, rows, cols) {
    cols <- cols[!is.na(cols)]
    if (length(rows) > 0 & length(cols) > 0) {
      col_filter <- which(tbl[, cols] == 'R')
      rows <- rows[rows %in% col_filter]
      tbl[rows, 'MDRO'] <<- to
    }
  }

  # join microorganisms
  tbl <- tbl %>% left_join_microorganisms(col_bactid)

  tbl$MDRO <- NA_integer_

  if (guideline$country$code == 'eucast') {
    # EUCAST ------------------------------------------------------------------
    # Table 5
    trans_tbl(4,
              which(tbl$family == 'Enterobacteriaceae'
                    | tbl$fullname %like% '^Pseudomonas aeruginosa'
                    | tbl$genus == 'Acinetobacter'),
              coli)
    trans_tbl(4,
              which(tbl$fullname %like% '^Salmonella Typhi'),
              c(carbapenems, fluoroquinolones))
    trans_tbl(4,
              which(tbl$fullname %like% '^Haemophilus influenzae'),
              c(cephalosporins_3rd, carbapenems, fluoroquinolones))
    trans_tbl(4,
              which(tbl$fullname %like% '^Moraxella catarrhalis'),
              c(cephalosporins_3rd, fluoroquinolones))
    trans_tbl(4,
              which(tbl$fullname %like% '^Neisseria meningitidis'),
              c(cephalosporins_3rd, fluoroquinolones))
    trans_tbl(4,
              which(tbl$fullname %like% '^Neisseria gonorrhoeae'),
              azit)
    # Table 6
    trans_tbl(4,
              which(tbl$fullname %like% '^Staphylococcus (aureus|epidermidis|coagulase negatief|hominis|haemolyticus|intermedius|pseudointermedius)'),
              c(vanc, teic, dapt, line, qida, tige))
    trans_tbl(4,
              which(tbl$genus == 'Corynebacterium'),
              c(vanc, teic, dapt, line, qida, tige))
    trans_tbl(4,
              which(tbl$fullname %like% '^Streptococcus pneumoniae'),
              c(carbapenems, vanc, teic, dapt, line, qida, tige, rifa))
    trans_tbl(4, # Sr. groups A/B/C/G
              which(tbl$fullname %like% '^Streptococcus (pyogenes|agalactiae|equisimilis|equi|zooepidemicus|dysgalactiae|anginosus)'),
              c(peni, cephalosporins, vanc, teic, dapt, line, qida, tige))
    trans_tbl(4,
              which(tbl$genus == 'Enterococcus'),
              c(dapt, line, tige, teic))
    trans_tbl(4,
              which(tbl$fullname %like% '^Enterococcus faecalis'),
              c(ampi, amox))
    # Table 7
    trans_tbl(4,
              which(tbl$genus == 'Bacteroides'),
              metr)
    trans_tbl(4,
              which(tbl$fullname %like% '^Clostridium difficile'),
              c(metr, vanc))
  }

  if (guideline$country$code == 'de') {
    # Germany -----------------------------------------------------------------
    stop("We are still working on German guidelines in this beta version.", call. = FALSE)
  }

  if (guideline$country$code == 'nl') {
    # Netherlands -------------------------------------------------------------
    aminoglycosides <- aminoglycosides[!is.na(aminoglycosides)]
    fluoroquinolones <- fluoroquinolones[!is.na(fluoroquinolones)]
    carbapenems <- carbapenems[!is.na(carbapenems)]

    # Table 1
    tbl[which(
      tbl$family == 'Enterobacteriaceae'
      & rowSums(tbl[, aminoglycosides] == 'R', na.rm = TRUE) >= 1
      & rowSums(tbl[, fluoroquinolones] == 'R', na.rm = TRUE) >= 1
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
      & rowSums(tbl[, fluoroquinolones] == 'R', na.rm = TRUE) >= 1
    ), 'MDRO'] <- 4
    # rest of Acinetobacter is negative
    tbl[which(
      tbl$genus == 'Acinetobacter'
      & tbl$MDRO == 1
    ), 'MDRO'] <- 2

    tbl[which(
      tbl$fullname %like% 'Stenotrophomonas maltophilia'
      & tbl[, trsu] == 'R'
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
           rowSums(tbl[, fluoroquinolones] == 'R', na.rm = TRUE) >= 1,
           tbl[, cfta] == 'R',
           tbl[, pita] == 'R') >= 3
    ), 'MDRO'] <- 4
    # rest of Pseudomonas is negative
    tbl[which(
      tbl$fullname %like% 'Pseudomonas aeruginosa'
      & tbl$MDRO == 1
    ), 'MDRO'] <- 2

    # Table 3
    tbl[which(
      tbl$fullname %like% 'Streptococcus pneumoniae'
      & tbl[, peni] == 'R'
    ), 'MDRO'] <- 4
    tbl[which(
      tbl$fullname %like% 'Streptococcus pneumoniae'
      & tbl[, vanc] == 'R'
    ), 'MDRO'] <- 4
    # rest of Streptococcus pneumoniae is negative
    tbl[which(
      tbl$fullname %like% 'Streptococcus pneumoniae'
      & tbl$MDRO == 1
    ), 'MDRO'] <- 2

    tbl[which(
      tbl$fullname %like% 'Enterococcus faecium'
      &  rowSums(tbl[, c(peni, vanc)] == 'R', na.rm = TRUE) >= 1
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
  MDRO(tbl = tbl, country = "nl", ...)
}

#' @rdname MDRO
#' @export
MRGN <- function(tbl, country = "de", ...) {
  MDRO(tbl = tbl, country = "de", ...)
}

#' @rdname MDRO
#' @export
EUCAST_exceptional_phenotypes <- function(tbl, country = "EUCAST", ...) {
  MDRO(tbl = tbl, country = "EUCAST", ...)
}
