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

#' Determine multidrug-resistant organisms (MDRO)
#'
#' Determine which isolates are multidrug-resistant organisms (MDRO) according to country-specific guidelines.
#' @param tbl table with antibiotic columns, like e.g. \code{amox} and \code{amcl}
#' @param country country code to determine guidelines. EUCAST rules will be used when left empty, see Details. Should be or a code from the \href{https://en.wikipedia.org/wiki/ISO_3166-1_alpha-2#Officially_assigned_code_elements}{list of ISO 3166-1 alpha-2 country codes}. Case-insensitive. Currently supported are \code{de} (Germany) and \code{nl} (the Netherlands).
#' @param info print progress
#' @inheritParams eucast_rules
#' @param metr column name of an antibiotic, see Antibiotics
#' @param verbose print additional info: missing antibiotic columns per parameter
#' @param ... parameters that are passed on to methods
#' @inheritSection eucast_rules Antibiotics
#' @details When \code{country} will be left blank, guidelines will be taken from EUCAST Expert Rules Version 3.1 "Intrinsic Resistance and Exceptional Phenotypes Tables" (\href{http://www.eucast.org/fileadmin/src/media/PDFs/EUCAST_files/Expert_Rules/Expert_rules_intrinsic_exceptional_V3.1.pdf}{link}).
#' @return Ordered factor with levels \code{Negative < Positive, unconfirmed < Positive}.
#' @rdname mdro
#' @importFrom dplyr %>%
#' @importFrom crayon red blue bold
#' @export
#' @inheritSection AMR Read more on our website!
#' @examples
#' library(dplyr)
#'
#' septic_patients %>%
#'   mutate(EUCAST = mdro(.),
#'          BRMO = brmo(.))
mdro <- function(tbl,
                 country = NULL,
                 col_mo = NULL,
                 info = TRUE,
                 amcl = guess_ab_col(),
                 amik = guess_ab_col(),
                 amox = guess_ab_col(),
                 ampi = guess_ab_col(),
                 azit = guess_ab_col(),
                 aztr = guess_ab_col(),
                 cefa = guess_ab_col(),
                 cfra = guess_ab_col(),
                 cfep = guess_ab_col(),
                 cfot = guess_ab_col(),
                 cfox = guess_ab_col(),
                 cfta = guess_ab_col(),
                 cftr = guess_ab_col(),
                 cfur = guess_ab_col(),
                 chlo = guess_ab_col(),
                 cipr = guess_ab_col(),
                 clar = guess_ab_col(),
                 clin = guess_ab_col(),
                 clox = guess_ab_col(),
                 coli = guess_ab_col(),
                 czol = guess_ab_col(),
                 dapt = guess_ab_col(),
                 doxy = guess_ab_col(),
                 erta = guess_ab_col(),
                 eryt = guess_ab_col(),
                 fosf = guess_ab_col(),
                 fusi = guess_ab_col(),
                 gent = guess_ab_col(),
                 imip = guess_ab_col(),
                 kana = guess_ab_col(),
                 levo = guess_ab_col(),
                 linc = guess_ab_col(),
                 line = guess_ab_col(),
                 mero = guess_ab_col(),
                 metr = guess_ab_col(),
                 mino = guess_ab_col(),
                 moxi = guess_ab_col(),
                 nali = guess_ab_col(),
                 neom = guess_ab_col(),
                 neti = guess_ab_col(),
                 nitr = guess_ab_col(),
                 novo = guess_ab_col(),
                 norf = guess_ab_col(),
                 oflo = guess_ab_col(),
                 peni = guess_ab_col(),
                 pipe = guess_ab_col(),
                 pita = guess_ab_col(),
                 poly = guess_ab_col(),
                 qida = guess_ab_col(),
                 rifa = guess_ab_col(),
                 roxi = guess_ab_col(),
                 siso = guess_ab_col(),
                 teic = guess_ab_col(),
                 tetr = guess_ab_col(),
                 tica = guess_ab_col(),
                 tige = guess_ab_col(),
                 tobr = guess_ab_col(),
                 trim = guess_ab_col(),
                 trsu = guess_ab_col(),
                 vanc = guess_ab_col(),
                 verbose = FALSE) {

  if (!is.data.frame(tbl)) {
    stop("`tbl` must be a data frame.", call. = FALSE)
  }

  # try to find columns based on type
  # -- mo
  if (is.null(col_mo)) {
    col_mo <- search_type_in_df(tbl = tbl, type = "mo")
  }
  if (is.null(col_mo)) {
    stop("`col_mo` must be set.", call. = FALSE)
  }

  # strip whitespaces
  if (length(country) > 1) {
    stop('`country` must be a length one character string.', call. = FALSE)
  }

  if (is.null(country)) {
    country <- 'EUCAST'
  }
  country <- trimws(country)
  if (tolower(country) != 'eucast' & !country %like% '^[a-z]{2}$') {
    stop('This is not a valid ISO 3166-1 alpha-2 country code: "', country, '". Please see ?mdro.', call. = FALSE)
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
    guideline$version <- 'Revision as of December 2017'
    guideline$source <- 'https://www.rivm.nl/Documenten_en_publicaties/Professioneel_Praktisch/Richtlijnen/Infectieziekten/WIP_Richtlijnen/WIP_Richtlijnen/Ziekenhuizen/WIP_richtlijn_BRMO_Bijzonder_Resistente_Micro_Organismen_ZKH'
  # add here more countries like this:
  # } else if (country$code == 'xx') {
  #   country$name <- 'country name'
  } else {
    stop('This country code is currently unsupported: ', guideline$country$code, call. = FALSE)
  }

  if (info == TRUE) {
    cat("Determining multidrug-resistant organisms (MDRO), according to:\n",
        "Guideline: ", red(paste0(guideline$name, ", ", guideline$version, "\n")),
        "Country  : ", red(paste0(guideline$country$name, "\n")),
        "Source   : ", blue(paste0(guideline$source, "\n")),
        "\n", sep = "")
  }

  # check columns
  if (identical(amcl, as.name("guess_ab_col"))) { amcl <- guess_ab_col(tbl, "amcl", verbose = verbose) }
  if (identical(amik, as.name("guess_ab_col"))) { amik <- guess_ab_col(tbl, "amik", verbose = verbose) }
  if (identical(amox, as.name("guess_ab_col"))) { amox <- guess_ab_col(tbl, "amox", verbose = verbose) }
  if (identical(ampi, as.name("guess_ab_col"))) { ampi <- guess_ab_col(tbl, "ampi", verbose = verbose) }
  if (identical(azit, as.name("guess_ab_col"))) { azit <- guess_ab_col(tbl, "azit", verbose = verbose) }
  if (identical(aztr, as.name("guess_ab_col"))) { aztr <- guess_ab_col(tbl, "aztr", verbose = verbose) }
  if (identical(cefa, as.name("guess_ab_col"))) { cefa <- guess_ab_col(tbl, "cefa", verbose = verbose) }
  if (identical(cfra, as.name("guess_ab_col"))) { cfra <- guess_ab_col(tbl, "cfra", verbose = verbose) }
  if (identical(cfep, as.name("guess_ab_col"))) { cfep <- guess_ab_col(tbl, "cfep", verbose = verbose) }
  if (identical(cfot, as.name("guess_ab_col"))) { cfot <- guess_ab_col(tbl, "cfot", verbose = verbose) }
  if (identical(cfox, as.name("guess_ab_col"))) { cfox <- guess_ab_col(tbl, "cfox", verbose = verbose) }
  if (identical(cfta, as.name("guess_ab_col"))) { cfta <- guess_ab_col(tbl, "cfta", verbose = verbose) }
  if (identical(cftr, as.name("guess_ab_col"))) { cftr <- guess_ab_col(tbl, "cftr", verbose = verbose) }
  if (identical(cfur, as.name("guess_ab_col"))) { cfur <- guess_ab_col(tbl, "cfur", verbose = verbose) }
  if (identical(chlo, as.name("guess_ab_col"))) { chlo <- guess_ab_col(tbl, "chlo", verbose = verbose) }
  if (identical(cipr, as.name("guess_ab_col"))) { cipr <- guess_ab_col(tbl, "cipr", verbose = verbose) }
  if (identical(clar, as.name("guess_ab_col"))) { clar <- guess_ab_col(tbl, "clar", verbose = verbose) }
  if (identical(clin, as.name("guess_ab_col"))) { clin <- guess_ab_col(tbl, "clin", verbose = verbose) }
  if (identical(clox, as.name("guess_ab_col"))) { clox <- guess_ab_col(tbl, "clox", verbose = verbose) }
  if (identical(coli, as.name("guess_ab_col"))) { coli <- guess_ab_col(tbl, "coli", verbose = verbose) }
  if (identical(czol, as.name("guess_ab_col"))) { czol <- guess_ab_col(tbl, "czol", verbose = verbose) }
  if (identical(dapt, as.name("guess_ab_col"))) { dapt <- guess_ab_col(tbl, "dapt", verbose = verbose) }
  if (identical(doxy, as.name("guess_ab_col"))) { doxy <- guess_ab_col(tbl, "doxy", verbose = verbose) }
  if (identical(erta, as.name("guess_ab_col"))) { erta <- guess_ab_col(tbl, "erta", verbose = verbose) }
  if (identical(eryt, as.name("guess_ab_col"))) { eryt <- guess_ab_col(tbl, "eryt", verbose = verbose) }
  if (identical(fosf, as.name("guess_ab_col"))) { fosf <- guess_ab_col(tbl, "fosf", verbose = verbose) }
  if (identical(fusi, as.name("guess_ab_col"))) { fusi <- guess_ab_col(tbl, "fusi", verbose = verbose) }
  if (identical(gent, as.name("guess_ab_col"))) { gent <- guess_ab_col(tbl, "gent", verbose = verbose) }
  if (identical(imip, as.name("guess_ab_col"))) { imip <- guess_ab_col(tbl, "imip", verbose = verbose) }
  if (identical(kana, as.name("guess_ab_col"))) { kana <- guess_ab_col(tbl, "kana", verbose = verbose) }
  if (identical(levo, as.name("guess_ab_col"))) { levo <- guess_ab_col(tbl, "levo", verbose = verbose) }
  if (identical(linc, as.name("guess_ab_col"))) { linc <- guess_ab_col(tbl, "linc", verbose = verbose) }
  if (identical(line, as.name("guess_ab_col"))) { line <- guess_ab_col(tbl, "line", verbose = verbose) }
  if (identical(mero, as.name("guess_ab_col"))) { mero <- guess_ab_col(tbl, "mero", verbose = verbose) }
  if (identical(metr, as.name("guess_ab_col"))) { metr <- guess_ab_col(tbl, "metr", verbose = verbose) }
  if (identical(mino, as.name("guess_ab_col"))) { mino <- guess_ab_col(tbl, "mino", verbose = verbose) }
  if (identical(moxi, as.name("guess_ab_col"))) { moxi <- guess_ab_col(tbl, "moxi", verbose = verbose) }
  if (identical(nali, as.name("guess_ab_col"))) { nali <- guess_ab_col(tbl, "nali", verbose = verbose) }
  if (identical(neom, as.name("guess_ab_col"))) { neom <- guess_ab_col(tbl, "neom", verbose = verbose) }
  if (identical(neti, as.name("guess_ab_col"))) { neti <- guess_ab_col(tbl, "neti", verbose = verbose) }
  if (identical(nitr, as.name("guess_ab_col"))) { nitr <- guess_ab_col(tbl, "nitr", verbose = verbose) }
  if (identical(novo, as.name("guess_ab_col"))) { novo <- guess_ab_col(tbl, "novo", verbose = verbose) }
  if (identical(norf, as.name("guess_ab_col"))) { norf <- guess_ab_col(tbl, "norf", verbose = verbose) }
  if (identical(oflo, as.name("guess_ab_col"))) { oflo <- guess_ab_col(tbl, "oflo", verbose = verbose) }
  if (identical(peni, as.name("guess_ab_col"))) { peni <- guess_ab_col(tbl, "peni", verbose = verbose) }
  if (identical(pipe, as.name("guess_ab_col"))) { pipe <- guess_ab_col(tbl, "pipe", verbose = verbose) }
  if (identical(pita, as.name("guess_ab_col"))) { pita <- guess_ab_col(tbl, "pita", verbose = verbose) }
  if (identical(poly, as.name("guess_ab_col"))) { poly <- guess_ab_col(tbl, "poly", verbose = verbose) }
  if (identical(qida, as.name("guess_ab_col"))) { qida <- guess_ab_col(tbl, "qida", verbose = verbose) }
  if (identical(rifa, as.name("guess_ab_col"))) { rifa <- guess_ab_col(tbl, "rifa", verbose = verbose) }
  if (identical(roxi, as.name("guess_ab_col"))) { roxi <- guess_ab_col(tbl, "roxi", verbose = verbose) }
  if (identical(siso, as.name("guess_ab_col"))) { siso <- guess_ab_col(tbl, "siso", verbose = verbose) }
  if (identical(teic, as.name("guess_ab_col"))) { teic <- guess_ab_col(tbl, "teic", verbose = verbose) }
  if (identical(tetr, as.name("guess_ab_col"))) { tetr <- guess_ab_col(tbl, "tetr", verbose = verbose) }
  if (identical(tica, as.name("guess_ab_col"))) { tica <- guess_ab_col(tbl, "tica", verbose = verbose) }
  if (identical(tige, as.name("guess_ab_col"))) { tige <- guess_ab_col(tbl, "tige", verbose = verbose) }
  if (identical(tobr, as.name("guess_ab_col"))) { tobr <- guess_ab_col(tbl, "tobr", verbose = verbose) }
  if (identical(trim, as.name("guess_ab_col"))) { trim <- guess_ab_col(tbl, "trim", verbose = verbose) }
  if (identical(trsu, as.name("guess_ab_col"))) { trsu <- guess_ab_col(tbl, "trsu", verbose = verbose) }
  if (identical(vanc, as.name("guess_ab_col"))) { vanc <- guess_ab_col(tbl, "vanc", verbose = verbose) }
  col.list <- c(amcl, amik, amox, ampi, azit, aztr, cefa, cfra, cfep, cfot,
                cfox, cfta, cftr, cfur, chlo, cipr, clar, clin, clox, coli,
                czol, dapt, doxy, erta, eryt, fosf, fusi, gent, imip, kana,
                levo, linc, line, mero, metr, mino, moxi, nali, neom, neti,
                nitr, novo, norf, oflo, peni, pipe, pita, poly, qida, rifa,
                roxi, siso, teic, tetr, tica, tige, tobr, trim, trsu, vanc)
  if (length(col.list) < 60) {
    warning('Some columns do not exist -- THIS MAY STRONGLY INFLUENCE THE OUTCOME.',
            immediate. = TRUE,
            call. = FALSE)
  }
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
  pipe <- col.list[pipe]
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

  ab_missing <- function(ab) {
    isTRUE(ab %in% c(NULL, NA)) | length(ab) == 0
  }

  # antibiotic classes
  aminoglycosides <- c(tobr, gent) # can also be kana but that one is often intrinsic R
  cephalosporins <- c(cfep, cfot, cfox, cfra, cfta, cftr, cfur, czol)
  cephalosporins_3rd <- c(cfot, cftr, cfta)
  carbapenems <- c(erta, imip, mero)
  fluoroquinolones <- c(oflo, cipr, levo, moxi)

  # helper function for editing the table
  trans_tbl <- function(to, rows, cols, any_all) {
    cols <- cols[!ab_missing(cols)]
    if (length(rows) > 0 & length(cols) > 0) {
      if (any_all == "any") {
        col_filter <- which(tbl[, cols] == 'R')
      } else if (any_all == "all") {
        col_filter <- tbl %>%
          mutate(index = 1:nrow(.)) %>%
          filter_at(vars(cols), all_vars(. == "R")) %>%
          pull((index))
      }
      rows <- rows[rows %in% col_filter]
      tbl[rows, 'MDRO'] <<- to
    }
  }

  tbl <- tbl %>%
    mutate_at(vars(col_mo), as.mo) %>%
    # join to microorganisms data set
    left_join_microorganisms(by = col_mo) %>%
    # add unconfirmed to where genus is available
    mutate(MDRO = ifelse(!is.na(genus), 1, NA_integer_))

  if (guideline$country$code == 'eucast') {
    # EUCAST ------------------------------------------------------------------
    # Table 5
    trans_tbl(3,
              which(tbl$family == 'Enterobacteriaceae'
                    | tbl$fullname %like% '^Pseudomonas aeruginosa'
                    | tbl$genus == 'Acinetobacter'),
              coli,
              "all")
    trans_tbl(3,
              which(tbl$fullname %like% '^Salmonella Typhi'),
              c(carbapenems, fluoroquinolones),
              "any")
    trans_tbl(3,
              which(tbl$fullname %like% '^Haemophilus influenzae'),
              c(cephalosporins_3rd, carbapenems, fluoroquinolones),
              "any")
    trans_tbl(3,
              which(tbl$fullname %like% '^Moraxella catarrhalis'),
              c(cephalosporins_3rd, fluoroquinolones),
              "any")
    trans_tbl(3,
              which(tbl$fullname %like% '^Neisseria meningitidis'),
              c(cephalosporins_3rd, fluoroquinolones),
              "any")
    trans_tbl(3,
              which(tbl$fullname %like% '^Neisseria gonorrhoeae'),
              azit,
              "any")
    # Table 6
    trans_tbl(3,
              which(tbl$fullname %like% '^Staphylococcus (aureus|epidermidis|coagulase negatief|hominis|haemolyticus|intermedius|pseudointermedius)'),
              c(vanc, teic, dapt, line, qida, tige),
              "any")
    trans_tbl(3,
              which(tbl$genus == 'Corynebacterium'),
              c(vanc, teic, dapt, line, qida, tige),
              "any")
    trans_tbl(3,
              which(tbl$fullname %like% '^Streptococcus pneumoniae'),
              c(carbapenems, vanc, teic, dapt, line, qida, tige, rifa),
              "any")
    trans_tbl(3, # Sr. groups A/B/C/G
              which(tbl$fullname %like% '^Streptococcus (pyogenes|agalactiae|equisimilis|equi|zooepidemicus|dysgalactiae|anginosus)'),
              c(peni, cephalosporins, vanc, teic, dapt, line, qida, tige),
              "any")
    trans_tbl(3,
              which(tbl$genus == 'Enterococcus'),
              c(dapt, line, tige, teic),
              "any")
    trans_tbl(3,
              which(tbl$fullname %like% '^Enterococcus faecalis'),
              c(ampi, amox),
              "any")
    # Table 7
    trans_tbl(3,
              which(tbl$genus == 'Bacteroides'),
              metr,
              "any")
    trans_tbl(3,
              which(tbl$fullname %like% '^Clostridium difficile'),
              c(metr, vanc),
              "any")
  }

  if (guideline$country$code == 'de') {
    # Germany -----------------------------------------------------------------
    stop("We are still working on German guidelines in this beta version.", call. = FALSE)
  }

  if (guideline$country$code == 'nl') {
    # Netherlands -------------------------------------------------------------
    aminoglycosides <- aminoglycosides[!ab_missing(aminoglycosides)]
    fluoroquinolones <- fluoroquinolones[!ab_missing(fluoroquinolones)]
    carbapenems <- carbapenems[!ab_missing(carbapenems)]

    # Table 1
    trans_tbl(3,
              which(tbl$family == 'Enterobacteriaceae'),
              c(aminoglycosides, fluoroquinolones),
              "all")

    trans_tbl(2,
              which(tbl$family == 'Enterobacteriaceae'),
              c(carbapenems),
              "any")

    # Table 2
    trans_tbl(2,
              which(tbl$genus == 'Acinetobacter'),
              c(carbapenems),
              "any")
    trans_tbl(3,
              which(tbl$genus == 'Acinetobacter'),
              c(aminoglycosides, fluoroquinolones),
              "all")

    trans_tbl(3,
              which(tbl$fullname %like% '^Stenotrophomonas maltophilia'),
              trsu,
              "all")

    if (!ab_missing(mero) & !ab_missing(imip)
        & !ab_missing(gent) & !ab_missing(tobr)
        & !ab_missing(cipr)
        & !ab_missing(cfta)
        & !ab_missing(pita) ) {
      tbl <- tbl %>% mutate(
        psae = 0,
        psae = ifelse(mero == "R" | imip == "R", psae + 1, psae),
        psae = ifelse(gent == "R" & tobr == "R", psae + 1, psae),
        psae = ifelse(cipr == "R", psae + 1, psae),
        psae = ifelse(cfta == "R", psae + 1, psae),
        psae = ifelse(pita == "R", psae + 1, psae),
        psae = ifelse(is.na(psae), 0, psae)
      )
    } else {
      tbl$psae <- 0
    }
    tbl[which(
      tbl$fullname %like% 'Pseudomonas aeruginosa'
      & tbl$psae >= 3
    ), 'MDRO'] <- 3

    # Table 3
    trans_tbl(3,
              which(tbl$fullname %like% 'Streptococcus pneumoniae'),
              peni,
              "all")
    trans_tbl(3,
              which(tbl$fullname %like% 'Streptococcus pneumoniae'),
              vanc,
              "all")
    trans_tbl(3,
              which(tbl$fullname %like% 'Enterococcus faecium'),
              c(peni, vanc),
              "all")
  }

  factor(x =  tbl$MDRO,
         levels = 1:3,
         labels = c('Negative', 'Positive, unconfirmed', 'Positive'),
         ordered = TRUE)
}

#' @rdname mdro
#' @export
brmo <- function(..., country = "nl") {
  mdro(..., country = "nl")
}

#' @rdname mdro
#' @export
mrgn <- function(tbl, country = "de", ...) {
  mdro(tbl = tbl, country = "de", ...)
}

#' @rdname mdro
#' @export
eucast_exceptional_phenotypes <- function(tbl, country = "EUCAST", ...) {
  mdro(tbl = tbl, country = "EUCAST", ...)
}
