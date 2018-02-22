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

#' EUCAST expert rules
#'
#' Apply expert rules (like intrinsic resistance), as defined by the European Committee on Antimicrobial Susceptibility Testing (EUCAST, \url{http://eucast.org}), see \emph{Source}.
#' @param tbl table with antibiotic columns, like e.g. \code{amox} and \code{amcl}
#' @param col_bactcode column name of the bacteria ID in \code{tbl} - should also be present in \code{bactlist$bactid}, see \code{\link{bactlist}}.
#' @param info print progress
#' @param amcl,amik,amox,ampi,azit,aztr,cefa,cfra,cfep,cfot,cfox,cfta,cftr,cfur,chlo,cipr,clar,clin,clox,coli,czol,dapt,doxy,erta,eryt,fosf,fusi,gent,imip,kana,levo,linc,line,mero,mino,moxi,nali,neom,neti,nitr,novo,norf,oflo,peni,pita,poly,qida,rifa,roxi,siso,teic,tetr,tica,tige,tobr,trim,trsu,vanc column names of antibiotics. Use \code{NA} to skip a column, like \code{tica = NA}. Non-existing column will be skipped.
#' @param ... parameters that are passed on to \code{EUCAST_rules}
#' @rdname EUCAST
#' @export
#' @importFrom dplyr %>% left_join select
#' @return table with edited variables of antibiotics.
#' @source
#'   EUCAST Expert Rules Version 2.0: \cr
#'   Leclercq et al. \strong{EUCAST expert rules in antimicrobial susceptibility testing.} \emph{Clin Microbiol Infect.} 2013;19(2):141-60. \cr
#'   \url{https://doi.org/10.1111/j.1469-0691.2011.03703.x} \cr
#'   \cr
#'   EUCAST Expert Rules Version 3.1: \cr
#'   \url{http://www.eucast.org/expert_rules_and_intrinsic_resistance}
#' @examples
#' a <- data.frame(bactid = c("STAAUR", "ESCCOL", "KLEPNE", "PSEAER"), 
#'                 vanc = "-",
#'                 amox = "-",
#'                 coli = "-",
#'                 cfta = "-",
#'                 cfur = "-",
#'                 stringsAsFactors = FALSE)
#' a
#' 
#' b <- EUCAST_rules(a)
#' b
EUCAST_rules <- function(tbl,
                         col_bactcode = 'bactid',
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
  
  if (!col_bactcode %in% colnames(tbl)) {
    stop('Column ', col_bactcode, ' not found.')
  }
  
  # kolommen controleren
  col.list <- c(amcl, amik, amox, ampi, azit, aztr, cefa, cfra, cfep,
                cfot, cfox, cfta, cftr, cfur, cipr, clar, clin, clox, coli, czol,
                dapt, doxy, erta, eryt, fusi, gent, imip, kana, levo, linc, line,
                mero, mino, moxi, nali, neom, neti, nitr, novo, norf, oflo, peni,
                pita, poly, qida, rifa, roxi, siso, teic, tetr, tica, tige, tobr,
                trim, trsu, vanc)
  col.list <- col.list[!is.na(col.list)]
  if (!all(col.list %in% colnames(tbl))) {
    if (info == TRUE) {
      cat('\n')
    }
    if (info == TRUE) {
      warning('These columns do not exist and will be ignored:\n',
              col.list[!(col.list %in% colnames(tbl))] %>% toString(),
              immediate. = TRUE,
              call. = FALSE)
    }
    if (!amcl %in% colnames(tbl)) { amcl <- NA }
    if (!amik %in% colnames(tbl)) { amik <- NA }
    if (!amox %in% colnames(tbl)) { amox <- NA }
    if (!ampi %in% colnames(tbl)) { ampi <- NA }
    if (!azit %in% colnames(tbl)) { azit <- NA }
    if (!aztr %in% colnames(tbl)) { aztr <- NA }
    if (!cefa %in% colnames(tbl)) { cefa <- NA }
    if (!cfra %in% colnames(tbl)) { cfra <- NA }
    if (!cfep %in% colnames(tbl)) { cfep <- NA }
    if (!cfot %in% colnames(tbl)) { cfot <- NA }
    if (!cfox %in% colnames(tbl)) { cfox <- NA }
    if (!cfta %in% colnames(tbl)) { cfta <- NA }
    if (!cftr %in% colnames(tbl)) { cftr <- NA }
    if (!cfur %in% colnames(tbl)) { cfur <- NA }
    if (!chlo %in% colnames(tbl)) { chlo <- NA }
    if (!cipr %in% colnames(tbl)) { cipr <- NA }
    if (!clar %in% colnames(tbl)) { clar <- NA }
    if (!clin %in% colnames(tbl)) { clin <- NA }
    if (!clox %in% colnames(tbl)) { clox <- NA }
    if (!coli %in% colnames(tbl)) { coli <- NA }
    if (!czol %in% colnames(tbl)) { czol <- NA }
    if (!dapt %in% colnames(tbl)) { dapt <- NA }
    if (!doxy %in% colnames(tbl)) { doxy <- NA }
    if (!erta %in% colnames(tbl)) { erta <- NA }
    if (!eryt %in% colnames(tbl)) { eryt <- NA }
    if (!fosf %in% colnames(tbl)) { fosf <- NA }
    if (!fusi %in% colnames(tbl)) { fusi <- NA }
    if (!gent %in% colnames(tbl)) { gent <- NA }
    if (!imip %in% colnames(tbl)) { imip <- NA }
    if (!kana %in% colnames(tbl)) { kana <- NA }
    if (!levo %in% colnames(tbl)) { levo <- NA }
    if (!linc %in% colnames(tbl)) { linc <- NA }
    if (!line %in% colnames(tbl)) { line <- NA }
    if (!mero %in% colnames(tbl)) { mero <- NA }
    if (!mino %in% colnames(tbl)) { mino <- NA }
    if (!moxi %in% colnames(tbl)) { moxi <- NA }
    if (!nali %in% colnames(tbl)) { nali <- NA }
    if (!neom %in% colnames(tbl)) { neom <- NA }
    if (!neti %in% colnames(tbl)) { neti <- NA }
    if (!nitr %in% colnames(tbl)) { nitr <- NA }
    if (!novo %in% colnames(tbl)) { novo <- NA }
    if (!norf %in% colnames(tbl)) { norf <- NA }
    if (!oflo %in% colnames(tbl)) { oflo <- NA }
    if (!peni %in% colnames(tbl)) { peni <- NA }
    if (!pita %in% colnames(tbl)) { pita <- NA }
    if (!poly %in% colnames(tbl)) { poly <- NA }
    if (!qida %in% colnames(tbl)) { qida <- NA }
    if (!rifa %in% colnames(tbl)) { rifa <- NA }
    if (!roxi %in% colnames(tbl)) { roxi <- NA }
    if (!siso %in% colnames(tbl)) { siso <- NA }
    if (!teic %in% colnames(tbl)) { teic <- NA }
    if (!tetr %in% colnames(tbl)) { tetr <- NA }
    if (!tica %in% colnames(tbl)) { tica <- NA }
    if (!tige %in% colnames(tbl)) { tige <- NA }
    if (!tobr %in% colnames(tbl)) { tobr <- NA }
    if (!trim %in% colnames(tbl)) { trim <- NA }
    if (!trsu %in% colnames(tbl)) { trsu <- NA }
    if (!vanc %in% colnames(tbl)) { vanc <- NA }
  }
  
  total <- 0
  
  # functie voor uitvoeren
  edit_rsi <- function(to, rows, cols) {
    #voortgang$tick()$print()
    cols <- cols[!is.na(cols)]
    if (length(rows) > 0 & length(cols) > 0) {
      tbl[rows, cols] <<- to
      total <<- total + (length(rows) * length(cols))
    }
  }
  
  # bactlist aan vastknopen (bestaande kolommen krijgen extra suffix)
  joinby <- colnames(AMR::bactlist)[1]
  names(joinby) <- col_bactcode
  tbl <- tbl %>% left_join(y = AMR::bactlist, by = joinby, suffix = c("_tempbactlist", ""))
  
  # antibioticagroepen
  aminoglycosiden <- c(tobr, gent, kana, neom, neti, siso)
  tetracyclines <- c(doxy, mino, tetr) # sinds EUCAST v3.1 is tige(cycline) apart
  polymyxines <- c(poly, coli)
  macroliden <- c(eryt, azit, roxi, clar) # sinds EUCAST v3.1 is clinda apart
  glycopeptiden <- c(vanc, teic)
  streptogramines <- qida # eigenlijk pristinamycine en quinupristine/dalfopristine
  cefalosporines <- c(cfep, cfot, cfox, cfra, cfta, cftr, cfur, czol)
  carbapenems <- c(erta, imip, mero)
  aminopenicillines <- c(ampi, amox)
  ureidopenicillines <- pita # eigenlijk ook azlo en mezlo
  fluorochinolonen <- c(oflo, cipr, norf, levo, moxi)
  
  if (info == TRUE) {
    cat('\nApplying EUCAST expert rules on',
        tbl[!is.na(tbl$genus),] %>% nrow(),
        'isolates according to "EUCAST Expert Rules Version 3.1"\n\n')
  }
  
  # Table 1: Intrinsic resistance in Enterobacteriaceae ----
  if (info == TRUE) {
    cat('...Table 1: Intrinsic resistance in Enterobacteriaceae\n')
  }
  #voortgang <- progress_estimated(17)
  # Intrisiek R voor groep
  edit_rsi(to = 'R',
           rows = which(tbl$family == 'Enterobacteriaceae'),
           cols = c(peni, glycopeptiden, fusi, macroliden, linc, streptogramines, rifa, dapt, line))
  # Citrobacter
  edit_rsi(to = 'R',
           rows = which(tbl$fullname %like% '^Citrobacter (koseri|amalonaticus|sedlakii|farmeri|rodentium)'),
           cols = c(aminopenicillines, tica))
  edit_rsi(to = 'R',
           rows = which(tbl$fullname %like% '^Citrobacter (freundii|braakii|murliniae|werkmanii|youngae)'),
           cols = c(aminopenicillines, amcl, czol, cfox))
  # Enterobacter
  edit_rsi(to = 'R',
           rows = which(tbl$fullname %like% '^Enterobacter cloacae'),
           cols = c(aminopenicillines, amcl, czol, cfox))
  edit_rsi(to = 'R',
           rows = which(tbl$fullname %like% '^Enterobacter aerogenes'),
           cols = c(aminopenicillines, amcl, czol, cfox))
  # Escherichia
  edit_rsi(to = 'R',
           rows = which(tbl$fullname %like% '^Escherichia hermanni'),
           cols = c(aminopenicillines, tica))
  # Hafnia
  edit_rsi(to = 'R',
           rows = which(tbl$fullname %like% '^Hafnia alvei'),
           cols = c(aminopenicillines, amcl, czol, cfox))
  # Klebsiella
  edit_rsi(to = 'R',
           rows = which(tbl$fullname %like% '^Klebsiella'),
           cols = c(aminopenicillines, tica))
  # Morganella / Proteus
  edit_rsi(to = 'R',
           rows = which(tbl$fullname %like% '^Morganella morganii'),
           cols = c(aminopenicillines, amcl, czol, tetracyclines, polymyxines, nitr))
  edit_rsi(to = 'R',
           rows = which(tbl$fullname %like% '^Proteus mirabilis'),
           cols = c(tetracyclines, tige, polymyxines, nitr))
  edit_rsi(to = 'R',
           rows = which(tbl$fullname %like% '^Proteus penneri'),
           cols = c(aminopenicillines, czol, cfur, tetracyclines, tige, polymyxines, nitr))
  edit_rsi(to = 'R',
           rows = which(tbl$fullname %like% '^Proteus vulgaris'),
           cols = c(aminopenicillines, czol, cfur, tetracyclines, tige, polymyxines, nitr))
  # Providencia
  edit_rsi(to = 'R',
           rows = which(tbl$fullname %like% '^Providencia rettgeri'),
           cols = c(aminopenicillines, amcl, czol, cfur, tetracyclines, tige, polymyxines, nitr))
  edit_rsi(to = 'R',
           rows = which(tbl$fullname %like% '^Providencia stuartii'),
           cols = c(aminopenicillines, amcl, czol, cfur, tetracyclines, tige, polymyxines, nitr))
  # Raoultella
  edit_rsi(to = 'R',
           rows = which(tbl$fullname %like% '^Raoultella'),
           cols = c(aminopenicillines, tica))
  # Serratia
  edit_rsi(to = 'R',
           rows = which(tbl$fullname %like% '^Serratia marcescens'),
           cols = c(aminopenicillines, amcl, czol, cfox, cfur, tetracyclines[tetracyclines != 'mino'], polymyxines, nitr))
  # Yersinia
  edit_rsi(to = 'R',
           rows = which(tbl$fullname %like% '^Yersinia enterocolitica'),
           cols = c(aminopenicillines, amcl, tica, czol, cfox))
  edit_rsi(to = 'R',
           rows = which(tbl$fullname %like% '^Yersinia pseudotuberculosis'),
           cols = c(poly, coli))
  
  
  # Table 2: Intrinsic resistance in non-fermentative Gram-negative bacteria ----
  if (info == TRUE) {
    cat('...Table 2: Intrinsic resistance in non-fermentative Gram-negative bacteria\n')
  }
  #voortgang <- progress_estimated(8)
  # Intrisiek R voor groep
  edit_rsi(to = 'R',
           rows = which(tbl$genus %in% c('Achromobacter',
                                         'Acinetobacter',
                                         'Alcaligenes',
                                         'Bordatella',
                                         'Burkholderia',
                                         'Elizabethkingia',
                                         'Flavobacterium',
                                         'Ochrobactrum',
                                         'Pseudomonas',
                                         'Stenotrophomonas')),
           cols = c(peni, cfox, cfur, glycopeptiden, fusi, macroliden, linc, streptogramines, rifa, dapt, line))
  # Acinetobacter
  edit_rsi(to = 'R',
           rows = which(tbl$fullname %like% '^Acinetobacter (baumannii|pittii|nosocomialis|calcoaceticus)'),
           cols = c(aminopenicillines, amcl, czol, cfot, cftr, aztr, erta, trim, fosf, tetracyclines[tetracyclines != 'mino']))
  # Achromobacter
  edit_rsi(to = 'R',
           rows = which(tbl$fullname %like% '^Achromobacter (xylosoxydans|xylosoxidans)'),
           cols = c(aminopenicillines, czol, cfot, cftr, erta))
  # Burkholderia
  edit_rsi(to = 'R',
           # onder 'Burkholderia cepacia complex' vallen deze species allemaal: PMID 16217180.
           rows = which(tbl$fullname %like% '^Burkholderia (cepacia|multivorans|cenocepacia|stabilis|vietnamiensis|dolosa|ambifaria|anthina|pyrrocinia|ubonensis)'),
           cols = c(aminopenicillines, amcl, tica, pita, czol, cfot, cftr, aztr, erta, cipr, chlo, aminoglycosiden, trim, fosf, polymyxines))
  # Elizabethkingia
  edit_rsi(to = 'R',
           rows = which(tbl$fullname %like% '^Elizabethkingia meningoseptic(a|um)'),
           cols = c(aminopenicillines, amcl, tica, czol, cfot, cftr, cfta, cfep, aztr, erta, imip, mero, polymyxines))
  # Ochrobactrum
  edit_rsi(to = 'R',
           rows = which(tbl$fullname %like% '^Ochrobactrum anthropi'),
           cols = c(aminopenicillines, amcl, tica, pita, czol, cfot, cftr, cfta, cfep, aztr, erta))
  # Pseudomonas
  edit_rsi(to = 'R',
           rows = which(tbl$fullname %like% '^Pseudomonas aeruginosa'),
           cols = c(aminopenicillines, amcl, czol, cfot, cftr, erta, chlo, kana, neom, trim, trsu, tetracyclines, tige))
  # Stenotrophomonas
  edit_rsi(to = 'R',
           rows = which(tbl$fullname %like% '^Stenotrophomonas maltophilia'),
           cols = c(aminopenicillines, amcl, tica, pita, czol, cfot, cftr, cfta, aztr, erta, imip, mero, aminoglycosiden, trim, fosf, tetr))
  
  
  # Table 3: Intrinsic resistance in other Gram-negative bacteria ----
  if (info == TRUE) {
    cat('...Table 3: Intrinsic resistance in other Gram-negative bacteria\n')
  }
  #voortgang <- progress_estimated(7)
  # Intrisiek R voor groep
  edit_rsi(to = 'R',
           rows = which(tbl$genus %in% c('Haemophilus',
                                         'Moraxella',
                                         'Neisseria',
                                         'Campylobacter')),
           cols = c(glycopeptiden, linc, dapt, line))
  # Haemophilus
  edit_rsi(to = 'R',
           rows = which(tbl$fullname %like% '^Haemophilus influenzae'),
           cols = c(fusi, streptogramines))
  # Moraxella
  edit_rsi(to = 'R',
           rows = which(tbl$fullname %like% '^Moraxella catarrhalis'),
           cols = trim)
  # Neisseria
  edit_rsi(to = 'R',
           rows = which(tbl$genus == 'Neisseria'),
           cols = trim)
  # Campylobacter
  edit_rsi(to = 'R',
           rows = which(tbl$fullname %like% '^Campylobacter fetus'),
           cols = c(fusi, streptogramines, trim, nali))
  edit_rsi(to = 'R',
           rows = which(tbl$fullname %like% '^Campylobacter (jejuni|coli)'),
           cols = c(fusi, streptogramines, trim))
  
  
  # Table 4: Intrinsic resistance in Gram-positive bacteria ----
  if (info == TRUE) {
    cat('...Table 4: Intrinsic resistance in Gram-positive bacteria\n')
  }
  #voortgang <- progress_estimated(14)
  # Intrisiek R voor groep
  edit_rsi(to = 'R',
           rows = which(tbl$gramstain %like% 'Positi(e|)(v|f)'),
           cols = c(aztr, polymyxines, nali))
  # Staphylococcus
  edit_rsi(to = 'R',
           rows = which(tbl$fullname %like% '^Staphylococcus saprophyticus'),
           cols = c(fusi, cfta, fosf, novo))
  edit_rsi(to = 'R',
           rows = which(tbl$fullname %like% '^Staphylococcus (cohnii|xylosus)'),
           cols = c(cfta, novo))
  edit_rsi(to = 'R',
           rows = which(tbl$fullname %like% '^Staphylococcus capitis'),
           cols = c(cfta, fosf))
  edit_rsi(to = 'R',
           rows = which(tbl$fullname %like% '^Staphylococcus (aureus|epidermidis|coagulase negatief|hominis|haemolyticus|intermedius|pseudointermedius)'),
           cols = cfta)
  # Streptococcus
  edit_rsi(to = 'R',
           rows = which(tbl$genus == 'Streptococcus'),
           cols = c(fusi, cfta, aminoglycosiden))
  # Enterococcus
  edit_rsi(to = 'R',
           rows = which(tbl$fullname %like% '^Enterococcus faecalis'),
           cols = c(fusi, cfta, cefalosporines[cefalosporines != cfta], aminoglycosiden, macroliden, clin, qida, trim, trsu))
  edit_rsi(to = 'R',
           rows = which(tbl$fullname %like% '^Enterococcus (gallinarum|casseliflavus)'),
           cols = c(fusi, cfta, cefalosporines[cefalosporines != cfta], aminoglycosiden, macroliden, clin, qida, vanc, trim, trsu))
  edit_rsi(to = 'R',
           rows = which(tbl$fullname %like% '^Enterococcus faecium'),
           cols = c(fusi, cfta, cefalosporines[cefalosporines != cfta], aminoglycosiden, macroliden, trim, trsu))
  # Corynebacterium
  edit_rsi(to = 'R',
           rows = which(tbl$genus == 'Corynebacterium'),
           cols = fosf)
  # Listeria
  edit_rsi(to = 'R',
           rows = which(tbl$fullname %like% '^Listeria monocytogenes'),
           cols = c(cfta, cefalosporines[cefalosporines != cfta]))
  # overig
  edit_rsi(to = 'R',
           rows = which(tbl$genus %in% c('Leuconostoc', 'Pediococcus')),
           cols = c(vanc, teic))
  edit_rsi(to = 'R',
           rows = which(tbl$genus == 'Lactobacillus'),
           cols = c(vanc, teic))
  edit_rsi(to = 'R',
           rows = which(tbl$fullname %like% '^Clostridium (ramosum|innocuum)'),
           cols = vanc)
  
  # Table 8: Interpretive rules for B-lactam agents and Gram-positive cocci ----
  if (info == TRUE) {
    cat('...Table 8: Interpretive rules for B-lactam agents and Gram-positive cocci\n')
  }
  #voortgang <- progress_estimated(2)
  # regel 8.3
  if (!is.na(peni)) {
    edit_rsi(to = 'S',
             rows = which(tbl$fullname %like% '^Streptococcus (pyogenes|agalactiae|dysgalactiae|groep A|groep B|groep C|groep G)'
                          & tbl[, peni] == 'S'),
             cols = c(aminopenicillines, cefalosporines, carbapenems))
  }
  # regel 8.6
  if (!is.na(ampi)) {
    edit_rsi(to = 'R',
             rows = which(tbl$genus == 'Enterococcus'
                          & tbl[, ampi] == 'R'),
             cols = c(ureidopenicillines, carbapenems))
  }
  if (!is.na(amox)) {
    edit_rsi(to = 'R',
             rows = which(tbl$genus == 'Enterococcus'
                          & tbl[, amox] == 'R'),
             cols = c(ureidopenicillines, carbapenems))
  }
  
  # Table 9: Interpretive rules for B-lactam agents and Gram-negative rods ----
  if (info == TRUE) {
    cat('...Table 9: Interpretive rules for B-lactam agents and Gram-negative rods\n')
  }
  #voortgang <- progress_estimated(1)
  # regel 9.3
  if (!is.na(tica) & !is.na(pita)) {
    edit_rsi(to = 'R',
             rows = which(tbl$family == 'Enterobacteriaceae'
                          & tbl[, tica] == 'R'
                          & tbl[, pita] == 'S'),
             cols = pita)
  }
  
  # Table 10: Interpretive rules for B-lactam agents and other Gram-negative bacteria ----
  if (info == TRUE) {
    cat('...Table 10: Interpretive rules for B-lactam agents and other Gram-negative bacteria\n')
  }
  #voortgang <- progress_estimated(1)
  # regel 10.2
  if (!is.na(ampi)) {
    # hiervoor moeten we eerst weten of ze B-lactamase-positief zijn
    # edit_rsi(to = 'R',
    #          rows = which(tbl$fullname %like% '^Haemophilus influenza'
    #                       & tbl[, ampi] == 'R'),
    #          cols = c(ampi, amox, amcl, pita, cfur))
  }
  
  # Table 11: Interpretive rules for macrolides, lincosamides, and streptogramins ----
  if (info == TRUE) {
    cat('...Table 11: Interpretive rules for macrolides, lincosamides, and streptogramins\n')
  }
  # regel 11.1
  if (!is.na(eryt)) {
    if (!is.na(azit)) {
      tbl[, azit] <- tbl[, eryt]
    }
    if (!is.na(clar)) {
      tbl[, clar] <- tbl[, eryt]
    }
  }
  
  # Table 12: Interpretive rules for aminoglycosides ----
  if (info == TRUE) {
    cat('...Table 12: Interpretive rules for aminoglycosides\n')
  }
  #voortgang <- progress_estimated(4)
  # regel 12.2
  if (!is.na(tobr)) {
    edit_rsi(to = 'R',
             rows = which(tbl$genus == 'Staphylococcus'
                          & tbl[, tobr] == 'R'),
             cols = c(kana, amik))
  }
  # regel 12.3
  if (!is.na(gent)) {
    edit_rsi(to = 'R',
             rows = which(tbl$genus == 'Staphylococcus'
                          & tbl[, gent] == 'R'),
             cols = aminoglycosiden)
  }
  # regel 12.8
  if (!is.na(gent) & !is.na(tobr)) {
    edit_rsi(to = 'R',
             rows = which(tbl$family == 'Enterobacteriaceae'
                          & tbl[, gent] == 'I'
                          & tbl[, tobr] == 'S'),
             cols = gent)
  }
  # regel 12.9
  if (!is.na(gent) & !is.na(tobr)) {
    edit_rsi(to = 'R',
             rows = which(tbl$family == 'Enterobacteriaceae'
                          & tbl[, tobr] == 'I'
                          & tbl[, gent] == 'R'),
             cols = tobr)
  }
  
  
  # Table 13: Interpretive rules for quinolones ----
  if (info == TRUE) {
    cat('...Table 13: Interpretive rules for quinolones\n')
  }
  #voortgang <- progress_estimated(4)
  # regel 13.2
  if (!is.na(moxi)) {
    edit_rsi(to = 'R',
             rows = which(tbl$genus == 'Staphylococcus'
                          & tbl[, moxi] == 'R'),
             cols = fluorochinolonen)
  }
  # regel 13.4
  if (!is.na(moxi)) {
    edit_rsi(to = 'R',
             rows = which(tbl$fullname %like% '^Streptococcus pneumoniae'
                          & tbl[, moxi] == 'R'),
             cols = fluorochinolonen)
  }
  # regel 13.5
  if (!is.na(cipr)) {
    edit_rsi(to = 'R',
             rows = which(tbl$family == 'Enterobacteriaceae'
                          & tbl[, cipr] == 'R'),
             cols = fluorochinolonen)
  }
  # regel 13.8
  if (!is.na(cipr)) {
    edit_rsi(to = 'R',
             rows = which(tbl$fullname %like% '^Neisseria gonorrhoeae'
                          & tbl[, cipr] == 'R'),
             cols = fluorochinolonen)
  }
  
  
  # Other ----
  if (info == TRUE) {
    cat('...Other\n')
  }
  #voortgang <- progress_estimated(2)
  if (!is.na(amcl)) {
    edit_rsi(to = 'R',
             rows = which(tbl[, amcl] == 'R'),
             cols = ampi)
  }
  if (!is.na(trsu)) {
    edit_rsi(to = 'R',
             rows = which(tbl[, trsu] == 'R'),
             cols = trim)
  }
  if (!is.na(ampi) & !is.na(amox)) {
    tbl[, amox] <- tbl %>% pull(ampi)
  }
  
  # Toegevoegde kolommen weer verwijderen
  bactlist.ncol <- ncol(AMR::bactlist) - 2
  tbl.ncol <- ncol(tbl)
  tbl <- tbl %>% select(-c((tbl.ncol - bactlist.ncol):tbl.ncol))
  # en eventueel toegevoegde suffix aan bestaande kolommen weer verwijderen
  colnames(tbl) <- gsub("_tempbactlist", "", colnames(tbl))
  
  if (info == TRUE) {
    cat('\nDone.\nExpert rules applied to', total, 'test results.\n')
  }
  
  tbl
}

#' @rdname EUCAST
#' @export
interpretive_reading <- function(...) {
  EUCAST_rules(...)
}

#' Poperties of a microorganism
#'
#' @param bactcode ID of a microorganisme, like \code{"STAAUR} and \code{"ESCCOL}
#' @param property One of the values \code{bactid}, \code{bactsys}, \code{family}, \code{genus}, \code{species}, \code{subspecies}, \code{fullname}, \code{type}, \code{gramstain}, \code{aerobic}
#' @export
#' @importFrom dplyr %>% filter select
#' @seealso \code{\link{bactlist}}
mo_property <- function(bactcode, property = 'fullname') {
  
  mocode <- as.character(bactcode)
  
  for (i in 1:length(mocode)) {
    bug <- mocode[i]
    
    if (!is.na(bug)) {
      result = tryCatch({
        mocode[i] <-
          AMR::bactlist %>%
          filter(bactid == bactcode) %>%
          select(property) %>%
          unlist() %>%
          as.character()
      }, error = function(error_condition) {
        warning('Code ', bug, ' not found in bacteria list.')
      }, finally = {
        if (mocode[i] == bug & !property %in% c('bactid', 'bactsys')) {
          mocode[i] <- NA
        }
      })
    }
    
  }
  mocode
}
