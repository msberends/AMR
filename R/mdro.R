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

#' Determine multidrug-resistant organisms (MDRO)
#'
#' Determine which isolates are multidrug-resistant organisms (MDRO) according to country-specific guidelines.
#' @param tbl table with antibiotic columns, like e.g. \code{amox} and \code{amcl}
#' @param country country code to determine guidelines. EUCAST rules will be used when left empty, see Details. Should be or a code from the \href{https://en.wikipedia.org/wiki/ISO_3166-1_alpha-2#Officially_assigned_code_elements}{list of ISO 3166-1 alpha-2 country codes}. Case-insensitive. Currently supported are \code{de} (Germany) and \code{nl} (the Netherlands).
#' @param info print progress
#' @inheritParams eucast_rules
#' @param metr column name of an antibiotic, see Antibiotics
#' @param ... parameters that are passed on to methods
#' @inheritSection eucast_rules Antibiotics
#' @details When \code{country} will be left blank, guidelines will be taken from EUCAST Expert Rules Version 3.1 "Intrinsic Resistance and Exceptional Phenotypes Tables" (\url{http://www.eucast.org/fileadmin/src/media/PDFs/EUCAST_files/Expert_Rules/Expert_rules_intrinsic_exceptional_V3.1.pdf}).
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
                 amcl = guess_ab(),
                 amik = guess_ab(),
                 amox = guess_ab(),
                 ampi = guess_ab(),
                 azit = guess_ab(),
                 aztr = guess_ab(),
                 cefa = guess_ab(),
                 cfra = guess_ab(),
                 cfep = guess_ab(),
                 cfot = guess_ab(),
                 cfox = guess_ab(),
                 cfta = guess_ab(),
                 cftr = guess_ab(),
                 cfur = guess_ab(),
                 chlo = guess_ab(),
                 cipr = guess_ab(),
                 clar = guess_ab(),
                 clin = guess_ab(),
                 clox = guess_ab(),
                 coli = guess_ab(),
                 czol = guess_ab(),
                 dapt = guess_ab(),
                 doxy = guess_ab(),
                 erta = guess_ab(),
                 eryt = guess_ab(),
                 fosf = guess_ab(),
                 fusi = guess_ab(),
                 gent = guess_ab(),
                 imip = guess_ab(),
                 kana = guess_ab(),
                 levo = guess_ab(),
                 linc = guess_ab(),
                 line = guess_ab(),
                 mero = guess_ab(),
                 metr = guess_ab(),
                 mino = guess_ab(),
                 moxi = guess_ab(),
                 nali = guess_ab(),
                 neom = guess_ab(),
                 neti = guess_ab(),
                 nitr = guess_ab(),
                 novo = guess_ab(),
                 norf = guess_ab(),
                 oflo = guess_ab(),
                 peni = guess_ab(),
                 pipe = guess_ab(),
                 pita = guess_ab(),
                 poly = guess_ab(),
                 qida = guess_ab(),
                 rifa = guess_ab(),
                 roxi = guess_ab(),
                 siso = guess_ab(),
                 teic = guess_ab(),
                 tetr = guess_ab(),
                 tica = guess_ab(),
                 tige = guess_ab(),
                 tobr = guess_ab(),
                 trim = guess_ab(),
                 trsu = guess_ab(),
                 vanc = guess_ab()) {

  if (!is.data.frame(tbl)) {
    stop("`tbl` must be a data frame.", call. = FALSE)
  }

  # try to find columns based on type
  # -- mo
  if (is.null(col_mo) & "mo" %in% lapply(tbl, class)) {
    col_mo <- colnames(tbl)[lapply(tbl, class) == "mo"][1]
    message(blue(paste0("NOTE: Using column `", bold(col_mo), "` as input for `col_mo`.")))
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
    guideline$version <- 'Revision of December 2017'
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
  if (identical(amcl, as.name("guess_ab"))) { amcl <- guess_ab(tbl, "amcl", verbose = info) }
  if (identical(amik, as.name("guess_ab"))) { amik <- guess_ab(tbl, "amik", verbose = info) }
  if (identical(amox, as.name("guess_ab"))) { amox <- guess_ab(tbl, "amox", verbose = info) }
  if (identical(ampi, as.name("guess_ab"))) { ampi <- guess_ab(tbl, "ampi", verbose = info) }
  if (identical(azit, as.name("guess_ab"))) { azit <- guess_ab(tbl, "azit", verbose = info) }
  if (identical(aztr, as.name("guess_ab"))) { aztr <- guess_ab(tbl, "aztr", verbose = info) }
  if (identical(cefa, as.name("guess_ab"))) { cefa <- guess_ab(tbl, "cefa", verbose = info) }
  if (identical(cfra, as.name("guess_ab"))) { cfra <- guess_ab(tbl, "cfra", verbose = info) }
  if (identical(cfep, as.name("guess_ab"))) { cfep <- guess_ab(tbl, "cfep", verbose = info) }
  if (identical(cfot, as.name("guess_ab"))) { cfot <- guess_ab(tbl, "cfot", verbose = info) }
  if (identical(cfox, as.name("guess_ab"))) { cfox <- guess_ab(tbl, "cfox", verbose = info) }
  if (identical(cfta, as.name("guess_ab"))) { cfta <- guess_ab(tbl, "cfta", verbose = info) }
  if (identical(cftr, as.name("guess_ab"))) { cftr <- guess_ab(tbl, "cftr", verbose = info) }
  if (identical(cfur, as.name("guess_ab"))) { cfur <- guess_ab(tbl, "cfur", verbose = info) }
  if (identical(chlo, as.name("guess_ab"))) { chlo <- guess_ab(tbl, "chlo", verbose = info) }
  if (identical(cipr, as.name("guess_ab"))) { cipr <- guess_ab(tbl, "cipr", verbose = info) }
  if (identical(clar, as.name("guess_ab"))) { clar <- guess_ab(tbl, "clar", verbose = info) }
  if (identical(clin, as.name("guess_ab"))) { clin <- guess_ab(tbl, "clin", verbose = info) }
  if (identical(clox, as.name("guess_ab"))) { clox <- guess_ab(tbl, "clox", verbose = info) }
  if (identical(coli, as.name("guess_ab"))) { coli <- guess_ab(tbl, "coli", verbose = info) }
  if (identical(czol, as.name("guess_ab"))) { czol <- guess_ab(tbl, "czol", verbose = info) }
  if (identical(dapt, as.name("guess_ab"))) { dapt <- guess_ab(tbl, "dapt", verbose = info) }
  if (identical(doxy, as.name("guess_ab"))) { doxy <- guess_ab(tbl, "doxy", verbose = info) }
  if (identical(erta, as.name("guess_ab"))) { erta <- guess_ab(tbl, "erta", verbose = info) }
  if (identical(eryt, as.name("guess_ab"))) { eryt <- guess_ab(tbl, "eryt", verbose = info) }
  if (identical(fosf, as.name("guess_ab"))) { fosf <- guess_ab(tbl, "fosf", verbose = info) }
  if (identical(fusi, as.name("guess_ab"))) { fusi <- guess_ab(tbl, "fusi", verbose = info) }
  if (identical(gent, as.name("guess_ab"))) { gent <- guess_ab(tbl, "gent", verbose = info) }
  if (identical(imip, as.name("guess_ab"))) { imip <- guess_ab(tbl, "imip", verbose = info) }
  if (identical(kana, as.name("guess_ab"))) { kana <- guess_ab(tbl, "kana", verbose = info) }
  if (identical(levo, as.name("guess_ab"))) { levo <- guess_ab(tbl, "levo", verbose = info) }
  if (identical(linc, as.name("guess_ab"))) { linc <- guess_ab(tbl, "linc", verbose = info) }
  if (identical(line, as.name("guess_ab"))) { line <- guess_ab(tbl, "line", verbose = info) }
  if (identical(mero, as.name("guess_ab"))) { mero <- guess_ab(tbl, "mero", verbose = info) }
  if (identical(metr, as.name("guess_ab"))) { metr <- guess_ab(tbl, "metr", verbose = info) }
  if (identical(mino, as.name("guess_ab"))) { mino <- guess_ab(tbl, "mino", verbose = info) }
  if (identical(moxi, as.name("guess_ab"))) { moxi <- guess_ab(tbl, "moxi", verbose = info) }
  if (identical(nali, as.name("guess_ab"))) { nali <- guess_ab(tbl, "nali", verbose = info) }
  if (identical(neom, as.name("guess_ab"))) { neom <- guess_ab(tbl, "neom", verbose = info) }
  if (identical(neti, as.name("guess_ab"))) { neti <- guess_ab(tbl, "neti", verbose = info) }
  if (identical(nitr, as.name("guess_ab"))) { nitr <- guess_ab(tbl, "nitr", verbose = info) }
  if (identical(novo, as.name("guess_ab"))) { novo <- guess_ab(tbl, "novo", verbose = info) }
  if (identical(norf, as.name("guess_ab"))) { norf <- guess_ab(tbl, "norf", verbose = info) }
  if (identical(oflo, as.name("guess_ab"))) { oflo <- guess_ab(tbl, "oflo", verbose = info) }
  if (identical(peni, as.name("guess_ab"))) { peni <- guess_ab(tbl, "peni", verbose = info) }
  if (identical(pipe, as.name("guess_ab"))) { pipe <- guess_ab(tbl, "pipe", verbose = info) }
  if (identical(pita, as.name("guess_ab"))) { pita <- guess_ab(tbl, "pita", verbose = info) }
  if (identical(poly, as.name("guess_ab"))) { poly <- guess_ab(tbl, "poly", verbose = info) }
  if (identical(qida, as.name("guess_ab"))) { qida <- guess_ab(tbl, "qida", verbose = info) }
  if (identical(rifa, as.name("guess_ab"))) { rifa <- guess_ab(tbl, "rifa", verbose = info) }
  if (identical(roxi, as.name("guess_ab"))) { roxi <- guess_ab(tbl, "roxi", verbose = info) }
  if (identical(siso, as.name("guess_ab"))) { siso <- guess_ab(tbl, "siso", verbose = info) }
  if (identical(teic, as.name("guess_ab"))) { teic <- guess_ab(tbl, "teic", verbose = info) }
  if (identical(tetr, as.name("guess_ab"))) { tetr <- guess_ab(tbl, "tetr", verbose = info) }
  if (identical(tica, as.name("guess_ab"))) { tica <- guess_ab(tbl, "tica", verbose = info) }
  if (identical(tige, as.name("guess_ab"))) { tige <- guess_ab(tbl, "tige", verbose = info) }
  if (identical(tobr, as.name("guess_ab"))) { tobr <- guess_ab(tbl, "tobr", verbose = info) }
  if (identical(trim, as.name("guess_ab"))) { trim <- guess_ab(tbl, "trim", verbose = info) }
  if (identical(trsu, as.name("guess_ab"))) { trsu <- guess_ab(tbl, "trsu", verbose = info) }
  if (identical(vanc, as.name("guess_ab"))) { vanc <- guess_ab(tbl, "vanc", verbose = info) }
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

  # antibiotic classes
  aminoglycosides <- c(tobr, gent) # can also be kana but that one is often intrinsic R
  cephalosporins <- c(cfep, cfot, cfox, cfra, cfta, cftr, cfur, czol)
  cephalosporins_3rd <- c(cfot, cftr, cfta)
  carbapenems <- c(erta, imip, mero)
  fluoroquinolones <- c(oflo, cipr, levo, moxi)

  # helper function for editing the table
  trans_tbl <- function(to, rows, cols, any_all) {
    cols <- cols[!is.na(cols)]
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
    aminoglycosides <- aminoglycosides[!is.na(aminoglycosides)]
    fluoroquinolones <- fluoroquinolones[!is.na(fluoroquinolones)]
    carbapenems <- carbapenems[!is.na(carbapenems)]

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

    if (!is.na(mero) & !is.na(imip)
        & !is.na(gent) & !is.na(tobr)
        & !is.na(cipr)
        & !is.na(cfta)
        & !is.na(pita) ) {
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
