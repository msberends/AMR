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
#' @param col_bactid column name of the bacteria ID in \code{tbl} - values of this column should be present in \code{microorganisms$bactid}, see \code{\link{microorganisms}}
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
#'   EUCAST Expert Rules Version 3.1 (Intrinsic Resistance and Exceptional Phenotypes Tables): \cr
#'   \url{http://www.eucast.org/fileadmin/src/media/PDFs/EUCAST_files/Expert_Rules/Expert_rules_intrinsic_exceptional_V3.1.pdf}
#' @examples
#' a <- data.frame(bactid = c("STAAUR",  # Staphylococcus aureus
#'                            "ENCFAE",  # Enterococcus faecalis
#'                            "ESCCOL",  # Escherichia coli
#'                            "KLEPNE",  # Klebsiella pneumoniae
#'                            "PSEAER"), # Pseudomonas aeruginosa
#'                 vanc = "-",           # Vancomycin
#'                 amox = "-",           # Amoxicillin
#'                 coli = "-",           # Colistin
#'                 cfta = "-",           # Ceftazidime
#'                 cfur = "-",           # Cefuroxime
#'                 stringsAsFactors = FALSE)
#' a
#' 
#' b <- EUCAST_rules(a)
#' b
EUCAST_rules <- function(tbl,
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

  EUCAST_VERSION <- "3.1"
  
  # support using columns as objects; the tidyverse way
  amcl <- quasiquotate(deparse(substitute(amcl)), amcl)
  amik <- quasiquotate(deparse(substitute(amik)), amik)
  amox <- quasiquotate(deparse(substitute(amox)), amox)
  ampi <- quasiquotate(deparse(substitute(ampi)), ampi)
  azit <- quasiquotate(deparse(substitute(azit)), azit)
  aztr <- quasiquotate(deparse(substitute(aztr)), aztr)
  cefa <- quasiquotate(deparse(substitute(cefa)), cefa)
  cfra <- quasiquotate(deparse(substitute(cfra)), cfra)
  cfep <- quasiquotate(deparse(substitute(cfep)), cfep)
  cfot <- quasiquotate(deparse(substitute(cfot)), cfot)
  cfox <- quasiquotate(deparse(substitute(cfox)), cfox)
  cfta <- quasiquotate(deparse(substitute(cfta)), cfta)
  cftr <- quasiquotate(deparse(substitute(cftr)), cftr)
  cfur <- quasiquotate(deparse(substitute(cfur)), cfur)
  chlo <- quasiquotate(deparse(substitute(chlo)), chlo)
  cipr <- quasiquotate(deparse(substitute(cipr)), cipr)
  clar <- quasiquotate(deparse(substitute(clar)), clar)
  clin <- quasiquotate(deparse(substitute(clin)), clin)
  clox <- quasiquotate(deparse(substitute(clox)), clox)
  coli <- quasiquotate(deparse(substitute(coli)), coli)
  czol <- quasiquotate(deparse(substitute(czol)), czol)
  dapt <- quasiquotate(deparse(substitute(dapt)), dapt)
  doxy <- quasiquotate(deparse(substitute(doxy)), doxy)
  erta <- quasiquotate(deparse(substitute(erta)), erta)
  eryt <- quasiquotate(deparse(substitute(eryt)), eryt)
  fosf <- quasiquotate(deparse(substitute(fosf)), fosf)
  fusi <- quasiquotate(deparse(substitute(fusi)), fusi)
  gent <- quasiquotate(deparse(substitute(gent)), gent)
  imip <- quasiquotate(deparse(substitute(imip)), imip)
  kana <- quasiquotate(deparse(substitute(kana)), kana)
  levo <- quasiquotate(deparse(substitute(levo)), levo)
  linc <- quasiquotate(deparse(substitute(linc)), linc)
  line <- quasiquotate(deparse(substitute(line)), line)
  mero <- quasiquotate(deparse(substitute(mero)), mero)
  mino <- quasiquotate(deparse(substitute(mino)), mino)
  moxi <- quasiquotate(deparse(substitute(moxi)), moxi)
  nali <- quasiquotate(deparse(substitute(nali)), nali)
  neom <- quasiquotate(deparse(substitute(neom)), neom)
  neti <- quasiquotate(deparse(substitute(neti)), neti)
  nitr <- quasiquotate(deparse(substitute(nitr)), nitr)
  novo <- quasiquotate(deparse(substitute(novo)), novo)
  norf <- quasiquotate(deparse(substitute(norf)), norf)
  oflo <- quasiquotate(deparse(substitute(oflo)), oflo)
  peni <- quasiquotate(deparse(substitute(peni)), peni)
  pita <- quasiquotate(deparse(substitute(pita)), pita)
  poly <- quasiquotate(deparse(substitute(poly)), poly)
  qida <- quasiquotate(deparse(substitute(qida)), qida)
  rifa <- quasiquotate(deparse(substitute(rifa)), rifa)
  roxi <- quasiquotate(deparse(substitute(roxi)), roxi)
  siso <- quasiquotate(deparse(substitute(siso)), siso)
  teic <- quasiquotate(deparse(substitute(teic)), teic)
  tetr <- quasiquotate(deparse(substitute(tetr)), tetr)
  tica <- quasiquotate(deparse(substitute(tica)), tica)
  tige <- quasiquotate(deparse(substitute(tige)), tige)
  tobr <- quasiquotate(deparse(substitute(tobr)), tobr)
  trim <- quasiquotate(deparse(substitute(trim)), trim)
  trsu <- quasiquotate(deparse(substitute(trsu)), trsu)
  vanc <- quasiquotate(deparse(substitute(vanc)), vanc)

  if (!col_bactid %in% colnames(tbl)) {
    stop('Column ', col_bactid, ' not found.')
  }
  
  # check columns
  col.list <- c(amcl, amik, amox, ampi, azit, aztr, cefa, cfra, cfep, cfot,
                cfox, cfta, cftr, cfur, chlo, cipr, clar, clin, clox, coli,
                czol, dapt, doxy, erta, eryt, fosf, fusi, gent, imip, kana,
                levo, linc, line, mero, mino, moxi, nali, neom, neti, nitr,
                novo, norf, oflo, peni, pita, poly, qida, rifa, roxi, siso,
                teic, tetr, tica, tige, tobr, trim, trsu, vanc)
  col.list <- col.list[!is.na(col.list)]
  col.list.bak <- col.list
  # are they available as upper case then?
  for (i in 1:length(col.list)) {
    if (toupper(col.list[i]) %in% colnames(tbl)) {
      col.list[i] <- toupper(col.list[i])
    } else if (tolower(col.list[i]) %in% colnames(tbl)) {
      col.list[i] <- tolower(col.list[i])
    } else if (!col.list[i] %in% colnames(tbl)) {
      col.list[i] <- NA
    }
  }
  if (!all(col.list %in% colnames(tbl))) {
    if (info == TRUE) {
      cat('\n')
    }
    if (info == TRUE) {
      warning('These columns do not exist and will be ignored: ',
              col.list.bak[!(col.list %in% colnames(tbl))] %>% toString(),
              immediate. = TRUE,
              call. = FALSE)
    }
  }
  
  amcl <- col.list[1]
  amik <- col.list[2]
  amox <- col.list[3]
  ampi <- col.list[4]
  azit <- col.list[5]
  aztr <- col.list[6]
  cefa <- col.list[7]
  cfra <- col.list[8]
  cfep <- col.list[9]
  cfot <- col.list[10]
  cfox <- col.list[11]
  cfta <- col.list[12]
  cftr <- col.list[13]
  cfur <- col.list[14]
  chlo <- col.list[15]
  cipr <- col.list[16]
  clar <- col.list[17]
  clin <- col.list[18]
  clox <- col.list[19]
  coli <- col.list[20]
  czol <- col.list[21]
  dapt <- col.list[22]
  doxy <- col.list[23]
  erta <- col.list[24]
  eryt <- col.list[25]
  fosf <- col.list[26]
  fusi <- col.list[27]
  gent <- col.list[28]
  imip <- col.list[29]
  kana <- col.list[30]
  levo <- col.list[31]
  linc <- col.list[32]
  line <- col.list[33]
  mero <- col.list[34]
  mino <- col.list[35]
  moxi <- col.list[36]
  nali <- col.list[37]
  neom <- col.list[38]
  neti <- col.list[39]
  nitr <- col.list[40]
  novo <- col.list[41]
  norf <- col.list[42]
  oflo <- col.list[43]
  peni <- col.list[44]
  pita <- col.list[45]
  poly <- col.list[46]
  qida <- col.list[47]
  rifa <- col.list[48]
  roxi <- col.list[49]
  siso <- col.list[50]
  teic <- col.list[51]
  tetr <- col.list[52]
  tica <- col.list[53]
  tige <- col.list[54]
  tobr <- col.list[55]
  trim <- col.list[56]
  trsu <- col.list[57]
  vanc <- col.list[58]
  
  total <- 0
  total_rows <- integer(0)
  
  # helper function for editing the table
  edit_rsi <- function(to, rows, cols) {
    cols <- cols[!is.na(cols)]
    if (length(rows) > 0 & length(cols) > 0) {
      tbl[rows, cols] <<- to
      total <<- total + (length(rows) * length(cols))
      total_rows <<- c(total_rows, rows)
    }
  }
  
  # join to microorganisms table
  joinby <- colnames(AMR::microorganisms)[1]
  names(joinby) <- col_bactid
  tbl <- tbl %>% left_join(y = AMR::microorganisms, by = joinby, suffix = c("_tempmicroorganisms", ""))
  
  # antibiotic classes
  aminoglycosides <- c(tobr, gent, kana, neom, neti, siso)
  tetracyclines <- c(doxy, mino, tetr) # since EUCAST v3.1 tige(cycline) is set apart
  polymyxins <- c(poly, coli)
  macrolides <- c(eryt, azit, roxi, clar) # since EUCAST v3.1 clinda is set apart
  glycopeptides <- c(vanc, teic)
  streptogramins <- qida # should officially also be pristinamycin and quinupristin/dalfopristin
  cephalosporins <- c(cfep, cfot, cfox, cfra, cfta, cftr, cfur, czol)
  carbapenems <- c(erta, imip, mero)
  aminopenicillins <- c(ampi, amox)
  ureidopenicillins <- pita # should officially also be azlo and mezlo
  fluoroquinolones <- c(oflo, cipr, norf, levo, moxi)
  
  if (info == TRUE) {
    cat(
      paste0(
        '\nApplying rules to ',
        tbl[!is.na(tbl$genus),] %>% nrow() %>% format(big.mark = ","),
        ' rows according to "EUCAST Expert Rules Version ', EUCAST_VERSION, '"\n')
    )
  }
  
  # Table 1: Intrinsic resistance in Enterobacteriaceae ----
  if (info == TRUE) {
    cat('...Table 1: Intrinsic resistance in Enterobacteriaceae\n')
  }
  # Intrisiek R for this group
  edit_rsi(to = 'R',
           rows = which(tbl$family == 'Enterobacteriaceae'),
           cols = c(peni, glycopeptides, fusi, macrolides, linc, streptogramins, rifa, dapt, line))
  # Citrobacter
  edit_rsi(to = 'R',
           rows = which(tbl$fullname %like% '^Citrobacter (koseri|amalonaticus|sedlakii|farmeri|rodentium)'),
           cols = c(aminopenicillins, tica))
  edit_rsi(to = 'R',
           rows = which(tbl$fullname %like% '^Citrobacter (freundii|braakii|murliniae|werkmanii|youngae)'),
           cols = c(aminopenicillins, amcl, czol, cfox))
  # Enterobacter
  edit_rsi(to = 'R',
           rows = which(tbl$fullname %like% '^Enterobacter cloacae'),
           cols = c(aminopenicillins, amcl, czol, cfox))
  edit_rsi(to = 'R',
           rows = which(tbl$fullname %like% '^Enterobacter aerogenes'),
           cols = c(aminopenicillins, amcl, czol, cfox))
  # Escherichia
  edit_rsi(to = 'R',
           rows = which(tbl$fullname %like% '^Escherichia hermanni'),
           cols = c(aminopenicillins, tica))
  # Hafnia
  edit_rsi(to = 'R',
           rows = which(tbl$fullname %like% '^Hafnia alvei'),
           cols = c(aminopenicillins, amcl, czol, cfox))
  # Klebsiella
  edit_rsi(to = 'R',
           rows = which(tbl$fullname %like% '^Klebsiella'),
           cols = c(aminopenicillins, tica))
  # Morganella / Proteus
  edit_rsi(to = 'R',
           rows = which(tbl$fullname %like% '^Morganella morganii'),
           cols = c(aminopenicillins, amcl, czol, tetracyclines, polymyxins, nitr))
  edit_rsi(to = 'R',
           rows = which(tbl$fullname %like% '^Proteus mirabilis'),
           cols = c(tetracyclines, tige, polymyxins, nitr))
  edit_rsi(to = 'R',
           rows = which(tbl$fullname %like% '^Proteus penneri'),
           cols = c(aminopenicillins, czol, cfur, tetracyclines, tige, polymyxins, nitr))
  edit_rsi(to = 'R',
           rows = which(tbl$fullname %like% '^Proteus vulgaris'),
           cols = c(aminopenicillins, czol, cfur, tetracyclines, tige, polymyxins, nitr))
  # Providencia
  edit_rsi(to = 'R',
           rows = which(tbl$fullname %like% '^Providencia rettgeri'),
           cols = c(aminopenicillins, amcl, czol, cfur, tetracyclines, tige, polymyxins, nitr))
  edit_rsi(to = 'R',
           rows = which(tbl$fullname %like% '^Providencia stuartii'),
           cols = c(aminopenicillins, amcl, czol, cfur, tetracyclines, tige, polymyxins, nitr))
  # Raoultella
  edit_rsi(to = 'R',
           rows = which(tbl$fullname %like% '^Raoultella'),
           cols = c(aminopenicillins, tica))
  # Serratia
  edit_rsi(to = 'R',
           rows = which(tbl$fullname %like% '^Serratia marcescens'),
           cols = c(aminopenicillins, amcl, czol, cfox, cfur, tetracyclines[tetracyclines != 'mino'], polymyxins, nitr))
  # Yersinia
  edit_rsi(to = 'R',
           rows = which(tbl$fullname %like% '^Yersinia enterocolitica'),
           cols = c(aminopenicillins, amcl, tica, czol, cfox))
  edit_rsi(to = 'R',
           rows = which(tbl$fullname %like% '^Yersinia pseudotuberculosis'),
           cols = c(poly, coli))
  
  
  # Table 2: Intrinsic resistance in non-fermentative Gram-negative bacteria ----
  if (info == TRUE) {
    cat('...Table 2: Intrinsic resistance in non-fermentative Gram-negative bacteria\n')
  }
  # Intrisiek R for this group
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
           cols = c(peni, cfox, cfur, glycopeptides, fusi, macrolides, linc, streptogramins, rifa, dapt, line))
  # Acinetobacter
  edit_rsi(to = 'R',
           rows = which(tbl$fullname %like% '^Acinetobacter (baumannii|pittii|nosocomialis|calcoaceticus)'),
           cols = c(aminopenicillins, amcl, czol, cfot, cftr, aztr, erta, trim, fosf, tetracyclines[tetracyclines != 'mino']))
  # Achromobacter
  edit_rsi(to = 'R',
           rows = which(tbl$fullname %like% '^Achromobacter (xylosoxydans|xylosoxidans)'),
           cols = c(aminopenicillins, czol, cfot, cftr, erta))
  # Burkholderia
  edit_rsi(to = 'R',
           # onder 'Burkholderia cepacia complex' vallen deze species allemaal: PMID 16217180.
           rows = which(tbl$fullname %like% '^Burkholderia (cepacia|multivorans|cenocepacia|stabilis|vietnamiensis|dolosa|ambifaria|anthina|pyrrocinia|ubonensis)'),
           cols = c(aminopenicillins, amcl, tica, pita, czol, cfot, cftr, aztr, erta, cipr, chlo, aminoglycosides, trim, fosf, polymyxins))
  # Elizabethkingia
  edit_rsi(to = 'R',
           rows = which(tbl$fullname %like% '^Elizabethkingia meningoseptic(a|um)'),
           cols = c(aminopenicillins, amcl, tica, czol, cfot, cftr, cfta, cfep, aztr, erta, imip, mero, polymyxins))
  # Ochrobactrum
  edit_rsi(to = 'R',
           rows = which(tbl$fullname %like% '^Ochrobactrum anthropi'),
           cols = c(aminopenicillins, amcl, tica, pita, czol, cfot, cftr, cfta, cfep, aztr, erta))
  # Pseudomonas
  edit_rsi(to = 'R',
           rows = which(tbl$fullname %like% '^Pseudomonas aeruginosa'),
           cols = c(aminopenicillins, amcl, czol, cfot, cftr, erta, chlo, kana, neom, trim, trsu, tetracyclines, tige))
  # Stenotrophomonas
  edit_rsi(to = 'R',
           rows = which(tbl$fullname %like% '^Stenotrophomonas maltophilia'),
           cols = c(aminopenicillins, amcl, tica, pita, czol, cfot, cftr, cfta, aztr, erta, imip, mero, aminoglycosides, trim, fosf, tetr))
  
  
  # Table 3: Intrinsic resistance in other Gram-negative bacteria ----
  if (info == TRUE) {
    cat('...Table 3: Intrinsic resistance in other Gram-negative bacteria\n')
  }
  # Intrisiek R for this group
  edit_rsi(to = 'R',
           rows = which(tbl$genus %in% c('Haemophilus',
                                         'Moraxella',
                                         'Neisseria',
                                         'Campylobacter')),
           cols = c(glycopeptides, linc, dapt, line))
  # Haemophilus
  edit_rsi(to = 'R',
           rows = which(tbl$fullname %like% '^Haemophilus influenzae'),
           cols = c(fusi, streptogramins))
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
           cols = c(fusi, streptogramins, trim, nali))
  edit_rsi(to = 'R',
           rows = which(tbl$fullname %like% '^Campylobacter (jejuni|coli)'),
           cols = c(fusi, streptogramins, trim))
  
  
  # Table 4: Intrinsic resistance in Gram-positive bacteria ----
  if (info == TRUE) {
    cat('...Table 4: Intrinsic resistance in Gram-positive bacteria\n')
  }
  # Intrisiek R for this group
  edit_rsi(to = 'R',
           rows = which(tbl$gramstain %like% 'Positi(e|)(v|f)'),
           cols = c(aztr, polymyxins, nali))
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
           cols = c(fusi, cfta, aminoglycosides))
  # Enterococcus
  edit_rsi(to = 'R',
           rows = which(tbl$fullname %like% '^Enterococcus faecalis'),
           cols = c(fusi, cfta, cephalosporins[cephalosporins != cfta], aminoglycosides, macrolides, clin, qida, trim, trsu))
  edit_rsi(to = 'R',
           rows = which(tbl$fullname %like% '^Enterococcus (gallinarum|casseliflavus)'),
           cols = c(fusi, cfta, cephalosporins[cephalosporins != cfta], aminoglycosides, macrolides, clin, qida, vanc, trim, trsu))
  edit_rsi(to = 'R',
           rows = which(tbl$fullname %like% '^Enterococcus faecium'),
           cols = c(fusi, cfta, cephalosporins[cephalosporins != cfta], aminoglycosides, macrolides, trim, trsu))
  # Corynebacterium
  edit_rsi(to = 'R',
           rows = which(tbl$genus == 'Corynebacterium'),
           cols = fosf)
  # Listeria
  edit_rsi(to = 'R',
           rows = which(tbl$fullname %like% '^Listeria monocytogenes'),
           cols = c(cfta, cephalosporins[cephalosporins != cfta]))
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
  # rule 8.3
  if (!is.na(peni)) {
    edit_rsi(to = 'S',
             rows = which(tbl$fullname %like% '^Streptococcus (pyogenes|agalactiae|dysgalactiae|groep A|groep B|groep C|groep G)'
                          & tbl[, peni] == 'S'),
             cols = c(aminopenicillins, cephalosporins, carbapenems))
  }
  # rule 8.6
  if (!is.na(ampi)) {
    edit_rsi(to = 'R',
             rows = which(tbl$genus == 'Enterococcus'
                          & tbl[, ampi] == 'R'),
             cols = c(ureidopenicillins, carbapenems))
  }
  if (!is.na(amox)) {
    edit_rsi(to = 'R',
             rows = which(tbl$genus == 'Enterococcus'
                          & tbl[, amox] == 'R'),
             cols = c(ureidopenicillins, carbapenems))
  }
  
  # Table 9: Interpretive rules for B-lactam agents and Gram-negative rods ----
  if (info == TRUE) {
    cat('...Table 9: Interpretive rules for B-lactam agents and Gram-negative rods\n')
  }
  # rule 9.3
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
  # rule 10.2
  if (!is.na(ampi)) {
    # you should know first if the are B-lactamase positive, so do not run for now
    # edit_rsi(to = 'R',
    #          rows = which(tbl$fullname %like% '^Haemophilus influenza'
    #                       & tbl[, ampi] == 'R'),
    #          cols = c(ampi, amox, amcl, pita, cfur))
  }
  
  # Table 11: Interpretive rules for macrolides, lincosamides, and streptogramins ----
  if (info == TRUE) {
    cat('...Table 11: Interpretive rules for macrolides, lincosamides, and streptogramins\n')
  }
  # rule 11.1
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
  # rule 12.2
  if (!is.na(tobr)) {
    edit_rsi(to = 'R',
             rows = which(tbl$genus == 'Staphylococcus'
                          & tbl[, tobr] == 'R'),
             cols = c(kana, amik))
  }
  # rule 12.3
  if (!is.na(gent)) {
    edit_rsi(to = 'R',
             rows = which(tbl$genus == 'Staphylococcus'
                          & tbl[, gent] == 'R'),
             cols = aminoglycosides)
  }
  # rule 12.8
  if (!is.na(gent) & !is.na(tobr)) {
    edit_rsi(to = 'R',
             rows = which(tbl$family == 'Enterobacteriaceae'
                          & tbl[, gent] == 'I'
                          & tbl[, tobr] == 'S'),
             cols = gent)
  }
  # rule 12.9
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
  # rule 13.2
  if (!is.na(moxi)) {
    edit_rsi(to = 'R',
             rows = which(tbl$genus == 'Staphylococcus'
                          & tbl[, moxi] == 'R'),
             cols = fluoroquinolones)
  }
  # rule 13.4
  if (!is.na(moxi)) {
    edit_rsi(to = 'R',
             rows = which(tbl$fullname %like% '^Streptococcus pneumoniae'
                          & tbl[, moxi] == 'R'),
             cols = fluoroquinolones)
  }
  # rule 13.5
  if (!is.na(cipr)) {
    edit_rsi(to = 'R',
             rows = which(tbl$family == 'Enterobacteriaceae'
                          & tbl[, cipr] == 'R'),
             cols = fluoroquinolones)
  }
  # rule 13.8
  if (!is.na(cipr)) {
    edit_rsi(to = 'R',
             rows = which(tbl$fullname %like% '^Neisseria gonorrhoeae'
                          & tbl[, cipr] == 'R'),
             cols = fluoroquinolones)
  }
  
  
  # Other ----
  if (info == TRUE) {
    cat('...Non-EUCAST: trim = R where trsu = R and ampi = R where amcl = R\n')
  }
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
  
  # Remove added columns again
  microorganisms.ncol <- ncol(AMR::microorganisms) - 2
  tbl.ncol <- ncol(tbl)
  tbl <- tbl %>% select(-c((tbl.ncol - microorganisms.ncol):tbl.ncol))
  # and remove added suffices
  colnames(tbl) <- gsub("_tempmicroorganisms", "", colnames(tbl))
  
  if (info == TRUE) {
    cat('Done.\n\nEUCAST Expert rules applied to',
        total_rows %>% unique() %>% length() %>% format(big.mark = ","),
        'different rows (isolates); edited a total of',
        total %>% format(big.mark = ","), 'test results.\n\n')
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
#' @param bactid ID of a microorganisme, like \code{"STAAUR} and \code{"ESCCOL}
#' @param property One of the values \code{bactid}, \code{bactsys}, \code{family}, \code{genus}, \code{species}, \code{subspecies}, \code{fullname}, \code{type}, \code{gramstain}, \code{aerobic}
#' @export
#' @importFrom dplyr %>% filter select
#' @seealso \code{\link{microorganisms}}
mo_property <- function(bactid, property = 'fullname') {
  
  mocode <- as.character(bactid)
  
  for (i in 1:length(mocode)) {
    bug <- mocode[i]
    
    if (!is.na(bug)) {
      result = tryCatch({
        mocode[i] <-
          AMR::microorganisms %>%
          filter(bactid == bactid) %>%
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
