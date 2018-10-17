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

#' EUCAST rules
#'
#' Apply susceptibility rules as defined by the European Committee on Antimicrobial Susceptibility Testing (EUCAST, \url{http://eucast.org}), see \emph{Source}. This includes (1) expert rules, (2) intrinsic resistance and (3) inferred resistance as defined in their breakpoint tables.
#' @param tbl table with antibiotic columns, like e.g. \code{amox} and \code{amcl}
#' @param col_mo column name of the microbial ID in \code{tbl} - values in this column should be present in \code{microorganisms$mo}, see \code{\link{microorganisms}}
#' @param info print progress
#' @param amcl,amik,amox,ampi,azit,azlo,aztr,cefa,cfep,cfot,cfox,cfra,cfta,cftr,cfur,chlo,cipr,clar,clin,clox,coli,czol,dapt,doxy,erta,eryt,fosf,fusi,gent,imip,kana,levo,linc,line,mero,mezl,mino,moxi,nali,neom,neti,nitr,norf,novo,oflo,oxac,peni,pita,poly,pris,qida,rifa,roxi,siso,teic,tetr,tica,tige,tobr,trim,trsu,vanc column name of an antibiotic. Use \code{NA} to skip a column, like \code{tica = NA}. Non-existing columns will anyway be skipped. See the Antibiotics section for an explanation of the abbreviations.
#' @param col_bactid Deprecated. Use \code{col_mo} instead.
#' @param verbose a logical to indicate whether extensive info should be printed to the console about which rows and columns are effected with their old and new values
#' @param ... parameters that are passed on to \code{EUCAST_rules}
#' @section Antibiotics:
#' Abbrevations of the column containing antibiotics:
#'
#'  \strong{amcl}: amoxicillin+clavulanic acid (\emph{J01CR02}),
#'  \strong{amik}: amikacin (\emph{J01GB06}),
#'  \strong{amox}: amoxicillin (\emph{J01CA04}),
#'  \strong{ampi}: ampicillin (\emph{J01CA01}),
#'  \strong{azit}: azithromycin (\emph{J01FA10}),
#'  \strong{azlo}: azlocillin (\emph{J01CA09}),
#'  \strong{aztr}: aztreonam (\emph{J01DF01}),
#'  \strong{cefa}: cefaloridine (\emph{J01DB02}),
#'  \strong{cfep}: cefepime (\emph{J01DE01}),
#'  \strong{cfot}: cefotaxime (\emph{J01DD01}),
#'  \strong{cfox}: cefoxitin (\emph{J01DC01}),
#'  \strong{cfra}: cefradine (\emph{J01DB09}),
#'  \strong{cfta}: ceftazidime (\emph{J01DD02}),
#'  \strong{cftr}: ceftriaxone (\emph{J01DD04}),
#'  \strong{cfur}: cefuroxime (\emph{J01DC02}),
#'  \strong{chlo}: chloramphenicol (\emph{J01BA01}),
#'  \strong{cipr}: ciprofloxacin (\emph{J01MA02}),
#'  \strong{clar}: clarithromycin (\emph{J01FA09}),
#'  \strong{clin}: clindamycin (\emph{J01FF01}),
#'  \strong{clox}: flucloxacillin (\emph{J01CF05}),
#'  \strong{coli}: colistin (\emph{J01XB01}),
#'  \strong{czol}: cefazolin (\emph{J01DB04}),
#'  \strong{dapt}: daptomycin (\emph{J01XX09}),
#'  \strong{doxy}: doxycycline (\emph{J01AA02}),
#'  \strong{erta}: ertapenem (\emph{J01DH03}),
#'  \strong{eryt}: erythromycin (\emph{J01FA01}),
#'  \strong{fosf}: fosfomycin (\emph{J01XX01}),
#'  \strong{fusi}: fusidic acid (\emph{J01XC01}),
#'  \strong{gent}: gentamicin (\emph{J01GB03}),
#'  \strong{imip}: imipenem (\emph{J01DH51}),
#'  \strong{kana}: kanamycin (\emph{J01GB04}),
#'  \strong{levo}: levofloxacin (\emph{J01MA12}),
#'  \strong{linc}: lincomycin (\emph{J01FF02}),
#'  \strong{line}: linezolid (\emph{J01XX08}),
#'  \strong{mero}: meropenem (\emph{J01DH02}),
#'  \strong{mezl}: mezlocillin (\emph{J01CA10}),
#'  \strong{mino}: minocycline (\emph{J01AA08}),
#'  \strong{moxi}: moxifloxacin (\emph{J01MA14}),
#'  \strong{nali}: nalidixic acid (\emph{J01MB02}),
#'  \strong{neom}: neomycin (\emph{J01GB05}),
#'  \strong{neti}: netilmicin (\emph{J01GB07}),
#'  \strong{nitr}: nitrofurantoin (\emph{J01XE01}),
#'  \strong{norf}: norfloxacin (\emph{J01MA06}),
#'  \strong{novo}: novobiocin (an ATCvet code: \emph{QJ01XX95}),
#'  \strong{oflo}: ofloxacin (\emph{J01MA01}),
#'  \strong{peni}: penicillin (\emph{J01RA01}),
#'  \strong{pita}: piperacillin+tazobactam (\emph{J01CR05}),
#'  \strong{poly}: polymyxin B (\emph{J01XB02}),
#'  \strong{pris}: pristinamycin (\emph{J01FG01}),
#'  \strong{qida}: quinupristin/dalfopristin (\emph{J01FG02}),
#'  \strong{rifa}: rifampicin (\emph{J04AB02}),
#'  \strong{roxi}: roxithromycin (\emph{J01FA06}),
#'  \strong{siso}: sisomicin (\emph{J01GB08}),
#'  \strong{teic}: teicoplanin (\emph{J01XA02}),
#'  \strong{tetr}: tetracycline (\emph{J01AA07}),
#'  \strong{tica}: ticarcillin (\emph{J01CA13}),
#'  \strong{tige}: tigecycline (\emph{J01AA12}),
#'  \strong{tobr}: tobramycin (\emph{J01GB01}),
#'  \strong{trim}: trimethoprim (\emph{J01EA01}),
#'  \strong{trsu}: sulfamethoxazole and trimethoprim (\emph{J01EE01}),
#'  \strong{vanc}: vancomycin (\emph{J01XA01}).
#' @keywords interpretive eucast reading resistance
#' @rdname EUCAST
#' @export
#' @importFrom dplyr %>% select pull
#' @importFrom crayon bold bgGreen bgYellow bgRed black green blue
#' @return table with edited variables of antibiotics.
#' @source
#'   \itemize{
#'     \item{
#'       EUCAST Expert Rules. Version 2.0, 2012. \cr
#'       Leclercq et al. \strong{EUCAST expert rules in antimicrobial susceptibility testing.} \emph{Clin Microbiol Infect.} 2013;19(2):141-60. \cr
#'       \url{https://doi.org/10.1111/j.1469-0691.2011.03703.x}
#'     }
#'     \item{
#'       EUCAST Expert Rules, Intrinsic Resistance and Exceptional Phenotypes Tables. Version 3.1, 2016.  \cr
#'       \url{http://www.eucast.org/fileadmin/src/media/PDFs/EUCAST_files/Expert_Rules/Expert_rules_intrinsic_exceptional_V3.1.pdf}
#'     }
#'     \item{
#'       EUCAST Breakpoint tables for interpretation of MICs and zone diameters. Version 8.1, 2018. \cr
#'       \url{http://www.eucast.org/fileadmin/src/media/PDFs/EUCAST_files/Breakpoint_tables/v_8.1_Breakpoint_Tables.xlsx}
#'     }
#'   }
#' @examples
#' a <- EUCAST_rules(septic_patients)
#'
#' a <- data.frame(mo = c("Staphylococcus aureus",
#'                        "Enterococcus faecalis",
#'                        "Escherichia coli",
#'                        "Klebsiella pneumoniae",
#'                        "Pseudomonas aeruginosa"),
#'                 vanc = "-",       # Vancomycin
#'                 amox = "-",       # Amoxicillin
#'                 coli = "-",       # Colistin
#'                 cfta = "-",       # Ceftazidime
#'                 cfur = "-",       # Cefuroxime
#'                 stringsAsFactors = FALSE)
#' a
#'
#' b <- EUCAST_rules(a) # 11 results will be forced as R due to intrinsic resistance
#' b
EUCAST_rules <- function(tbl,
                         col_mo = 'mo',
                         info = TRUE,
                         amcl = 'amcl',
                         amik = 'amik',
                         amox = 'amox',
                         ampi = 'ampi',
                         azit = 'azit',
                         azlo = 'azlo',
                         aztr = 'aztr',
                         cefa = 'cefa',
                         cfep = 'cfep',
                         cfot = 'cfot',
                         cfox = 'cfox',
                         cfra = 'cfra',
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
                         mezl = 'mezl',
                         mino = 'mino',
                         moxi = 'moxi',
                         nali = 'nali',
                         neom = 'neom',
                         neti = 'neti',
                         nitr = 'nitr',
                         norf = 'norf',
                         novo = 'novo',
                         oflo = 'oflo',
                         oxac = 'oxac',
                         peni = 'peni',
                         pita = 'pita',
                         poly = 'poly',
                         pris = 'pris',
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
                         vanc = 'vanc',
                         col_bactid = 'bactid',
                         verbose = FALSE) {

  EUCAST_VERSION_BREAKPOINTS <- "8.1, 2018"
  EUCAST_VERSION_EXPERT_RULES <- "3.1, 2016"

  if (col_bactid %in% colnames(tbl)) {
    col_mo <- col_bactid
    warning("Use of `col_bactid` is deprecated. Use `col_mo` instead.")
  }
  if (!col_mo %in% colnames(tbl)) {
    stop('Column ', col_mo, ' not found.', call. = FALSE)
  }

  warned <- FALSE
  changed_results <- 0

  txt_error <- function() { cat("", bgRed(black(" ERROR ")), "\n") }
  txt_warning <- function() { if (warned == FALSE) { cat("", bgYellow(black(" WARNING ")), "\n") }; warned <<- TRUE }
  txt_ok <- function() {
    if (warned == FALSE) {
      if (changed_results > 0) {
        if (changed_results == 1) {
          cat(blue(" (1 change)\n"))
        } else {
          cat(blue(paste0(" (", changed_results, " changes)\n")))
        }
      } else {
        cat(green(" (no changes)\n"))
      }
      warned <<- FALSE
    }
  }

  # check columns
  col.list <- c(amcl, amik, amox, ampi, azit, azlo, aztr, cefa, cfra, cfep, cfot,
                cfox, cfta, cftr, cfur, chlo, cipr, clar, clin, clox, coli,
                czol, dapt, doxy, erta, eryt, fosf, fusi, gent, imip, kana,
                levo, linc, line, mero, mezl, mino, moxi, nali, neom, neti, nitr,
                novo, norf, oflo, oxac, peni, pita, poly, pris,  qida, rifa, roxi, siso,
                teic, tetr, tica, tige, tobr, trim, trsu, vanc)
  col.list <- check_available_columns(tbl = tbl, col.list = col.list, info = info)
  amcl <- col.list[amcl]
  amik <- col.list[amik]
  amox <- col.list[amox]
  ampi <- col.list[ampi]
  azit <- col.list[azit]
  azlo <- col.list[azlo]
  aztr <- col.list[aztr]
  cefa <- col.list[cefa]
  cfep <- col.list[cfep]
  cfot <- col.list[cfot]
  cfox <- col.list[cfox]
  cfra <- col.list[cfra]
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
  mezl <- col.list[mezl]
  mino <- col.list[mino]
  moxi <- col.list[moxi]
  nali <- col.list[nali]
  neom <- col.list[neom]
  neti <- col.list[neti]
  nitr <- col.list[nitr]
  norf <- col.list[norf]
  novo <- col.list[novo]
  oflo <- col.list[oflo]
  oxac <- col.list[oxac]
  peni <- col.list[peni]
  pita <- col.list[pita]
  poly <- col.list[poly]
  pris <- col.list[pris]
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

  amount_changed <- 0
  amount_affected_rows <- integer(0)
  verbose_info <- ""

  # helper function for editing the table
  edit_rsi <- function(to, rule, rows, cols) {
    cols <- cols[!is.na(cols)]
    if (length(rows) > 0 & length(cols) > 0) {
      before <- as.character(unlist(as.list(tbl_original[rows, cols])))
      tryCatch(
        # insert into original table
        tbl_original[rows, cols] <<- to,
        warning = function(w) {
          if (w$message %like% 'invalid factor level') {
            warning('Value "', to, '" could not be applied to column(s) `', paste(cols, collapse = '`, `'), '` because this value is not an existing factor level.', call. = FALSE)
          } else {
            warning(w$message, call. = FALSE)
          }
          txt_warning()
        },
        error = function(e) {
          txt_error()
          stop(e, call. = FALSE)
        }
      )
      after <- as.character(unlist(as.list(tbl_original[rows, cols])))
      amount_changed <<- amount_changed + sum(before != after, na.rm = TRUE)
      amount_affected_rows <<- unique(c(amount_affected_rows, rows))
      changed_results <<- changed_results + sum(before != after, na.rm = TRUE) # will be reset at start of every rule

      if (verbose == TRUE) {
        verbose_info <<- paste0(verbose_info,
                                "\n\nRule Type: ", rule[1],
                                "\nRule Set:  ", rule[2],
                                "\nEffect:    Forced '", to, "' (",
                                "actually changed ", sum(before != after, na.rm = TRUE),
                                "): cols '", paste(cols, collapse = "', '"),
                                "' of rows ", paste(rows, collapse = ", "))
      }
    }
  }
  na.rm <- function(col) {
    if (is.na(col)) {
      ""
    } else {
      col
    }
  }

  # join to microorganisms data set
  col_mo_original <- NULL
  if (!tbl %>% pull(col_mo) %>% is.mo()) {
    col_mo_original <- tbl %>% pull(col_mo)
    tbl[, col_mo] <- as.mo(tbl[, col_mo])
  }
  tbl_original <- tbl
  tbl <- tbl %>% left_join_microorganisms(by = col_mo, suffix = c("_oldcols", ""))

  # antibiotic classes
  aminoglycosides <- c(tobr, gent, kana, neom, neti, siso)
  tetracyclines <- c(doxy, mino, tetr) # since EUCAST v3.1 tige(cycline) is set apart
  polymyxins <- c(poly, coli)
  macrolides <- c(eryt, azit, roxi, clar) # since EUCAST v3.1 clinda is set apart
  glycopeptides <- c(vanc, teic)
  streptogramins <- c(qida, pris) # should officially also be quinupristin/dalfopristin
  cephalosporins <- c(cfep, cfot, cfox, cfra, cfta, cftr, cfur, czol)
  carbapenems <- c(erta, imip, mero)
  aminopenicillins <- c(ampi, amox)
  ureidopenicillins <- c(pita, azlo, mezl)
  fluoroquinolones <- c(oflo, cipr, norf, levo, moxi)
  all_betalactam <- c(aminopenicillins, ureidopenicillins, cephalosporins, carbapenems, amcl, oxac, clox, peni)

  if (info == TRUE) {
    cat("Rules by the European Committee on Antimicrobial Susceptibility Testing (EUCAST)\n")
  }

  # since ampicillin ^= amoxicillin, get the first from the latter
  if (!is.na(ampi) & !is.na(amox)) {
    rule_group <- "Ampicillin susceptibility"
    rule <- "Get ampicillin results from amoxicillin where ampicillin is missing"
    edit_rsi(to = 'S',
             rule = c(rule_group, rule),
             rows = which(tbl[, amox] == 'S' & !tbl[, ampi] %in% c("S", "I", "R")),
             cols = ampi)
    edit_rsi(to = 'I',
             rule = c(rule_group, rule),
             rows = which(tbl[, amox] == 'I' & !tbl[, ampi] %in% c("S", "I", "R")),
             cols = ampi)
    edit_rsi(to = 'R',
             rule = c(rule_group, rule),
             rows = which(tbl[, amox] == 'R' & !tbl[, ampi] %in% c("S", "I", "R")),
             cols = ampi)
  }

  # BREAKPOINTS -------------------------------------------------------------

  if (info == TRUE) {
    cat(bold(paste0('\nEUCAST Clinical Breakpoints (v', EUCAST_VERSION_BREAKPOINTS, ')\n')))
  }
  rule_group <- "Breakpoints"

  # Enterobacteriales (Order) ----
  rule <- 'Enterobacteriales (Order)'
  if (info == TRUE) {
    warned <- FALSE
    changed_results <- 0
    cat(rule)
  }

  if (!is.na(ampi)) {
    edit_rsi(to = 'S',
             rule = c(rule_group, rule),
             rows = which(tbl$order == 'Enterobacteriales'
                          & tbl[, ampi] == 'S'),
             cols = amox)
  }
  if (!is.na(ampi)) {
    edit_rsi(to = 'I',
             rule = c(rule_group, rule),
             rows = which(tbl$order == 'Enterobacteriales'
                          & tbl[, ampi] == 'I'),
             cols = amox)
  }
  if (!is.na(ampi)) {
    edit_rsi(to = 'R',
             rule = c(rule_group, rule),
             rows = which(tbl$order == 'Enterobacteriales'
                          & tbl[, ampi] == 'R'),
             cols = amox)
  }
  if (info == TRUE) {
    txt_ok()
  }
  # Staphylococcus ----
  rule <- 'Staphylococcus'
  if (info == TRUE) {
    warned <- FALSE
    changed_results <- 0
    cat(rule)
  }
  if (!is.na(peni) & !is.na(cfox)) {
    edit_rsi(to = 'S',
             rule = c(rule_group, rule),
             rows = which(tbl$genus == "Staphylococcus"
                          & tbl[, peni] == 'S'
                          & tbl[, cfox] == 'S'),
             cols = c(ampi, amox, pita, tica))
    edit_rsi(to = 'S',
             rule = c(rule_group, rule),
             rows = which(tbl$genus == "Staphylococcus"
                          & tbl[, peni] == 'R'
                          & tbl[, cfox] == 'S'),
             cols = c(oxac, clox))
  }
  if (!is.na(cfox)) {
    edit_rsi(to = 'R',
             rule = c(rule_group, rule),
             rows = which(tbl$genus == "Staphylococcus"
                          & tbl[, cfox] == 'R'),
             cols = all_betalactam)
  }
  if (!is.na(ampi)) {
    edit_rsi(to = 'S',
             rule = c(rule_group, rule),
             rows = which(tbl$fullname %like% "^Staphylococcus saprophyticus"
                          & tbl[, ampi] == 'S'),
             cols = c(ampi, amox, amcl, pita))
  }
  if (!is.na(cfox)) {
    # inferred from cefoxitin
    edit_rsi(to = 'S',
             rule = c(rule_group, rule),
             rows = which(tbl$genus == "Staphylococcus"
                          & tbl[, cfox] == 'S'),
             cols = c(carbapenems, cephalosporins[cephalosporins != na.rm(cfta)]))
    edit_rsi(to = 'I',
             rule = c(rule_group, rule),
             rows = which(tbl$genus == "Staphylococcus"
                          & tbl[, cfox] == 'I'),
             cols = c(carbapenems, cephalosporins[cephalosporins != na.rm(cfta)]))
    edit_rsi(to = 'R',
             rule = c(rule_group, rule),
             rows = which(tbl$genus == "Staphylococcus"
                          & tbl[, cfox] == 'R'),
             cols = c(carbapenems, cephalosporins[cephalosporins != na.rm(cfta)]))
  }
  if (!is.na(norf)) {
    edit_rsi(to = 'S',
             rule = c(rule_group, rule),
             rows = which(tbl$genus == "Staphylococcus"
                          & tbl[, norf] == 'S'),
             cols = c(cipr, levo, moxi, oflo))
  }
  if (!is.na(eryt)) {
    edit_rsi(to = 'S',
             rule = c(rule_group, rule),
             rows = which(tbl$genus == "Staphylococcus"
                          & tbl[, eryt] == 'S'),
             cols = c(azit, clar, roxi))
    edit_rsi(to = 'I',
             rule = c(rule_group, rule),
             rows = which(tbl$genus == "Staphylococcus"
                          & tbl[, eryt] == 'I'),
             cols = c(azit, clar, roxi))
    edit_rsi(to = 'R',
             rule = c(rule_group, rule),
             rows = which(tbl$genus == "Staphylococcus"
                          & tbl[, eryt] == 'R'),
             cols = c(azit, clar, roxi))
  }
  if (!is.na(tetr)) {
    edit_rsi(to = 'S',
             rule = c(rule_group, rule),
             rows = which(tbl$genus == "Staphylococcus"
                          & tbl[, tetr] == 'S'),
             cols = c(doxy, mino))
  }
  if (info == TRUE) {
    txt_ok()
  }
  # Enterococcus ----
  rule <- 'Enterococcus'
  if (info == TRUE) {
    warned <- FALSE
    changed_results <- 0
    cat(rule)
  }
  if (!is.na(peni)) {
    edit_rsi(to = 'R',
             rule = c(rule_group, rule),
             rows = which(tbl$fullname %like% "^Enterococcus faecium"
                          & tbl[, peni] == 'R'),
             cols = all_betalactam)
  }
  if (!is.na(ampi)) {
    edit_rsi(to = 'S',
             rule = c(rule_group, rule),
             rows = which(tbl$genus == "Enterococcus"
                          & tbl[, ampi] == 'S'),
             cols = c(amox, amcl, pita))
    edit_rsi(to = 'I',
             rule = c(rule_group, rule),
             rows = which(tbl$genus == "Enterococcus"
                          & tbl[, ampi] == 'I'),
             cols = c(amox, amcl, pita))
    edit_rsi(to = 'R',
             rule = c(rule_group, rule),
             rows = which(tbl$genus == "Enterococcus"
                          & tbl[, ampi] == 'R'),
             cols = c(amox, amcl, pita))
  }
  if (!is.na(norf)) {
    edit_rsi(to = 'S',
             rule = c(rule_group, rule),
             rows = which(tbl$genus == "Enterococcus"
                          & tbl[, norf] == 'S'),
             cols = c(cipr, levo))
    edit_rsi(to = 'I',
             rule = c(rule_group, rule),
             rows = which(tbl$genus == "Enterococcus"
                          & tbl[, norf] == 'I'),
             cols = c(cipr, levo))
    edit_rsi(to = 'R',
             rule = c(rule_group, rule),
             rows = which(tbl$genus == "Enterococcus"
                          & tbl[, norf] == 'R'),
             cols = c(cipr, levo))
  }
  if (info == TRUE) {
    txt_ok()
  }
  # Streptococcus A, B, C, G----
  rule <- 'Streptococcus A, B, C, G'
  if (info == TRUE) {
    warned <- FALSE
    changed_results <- 0
    cat(rule)
  }
  if (!is.na(peni)) {
    edit_rsi(to = 'S',
             rule = c(rule_group, rule),
             rows = which(tbl$fullname %like% "^Streptococcus (pyogenes|agalactiae|dysgalactiae|group A|group B|group C|group G)"
                          & tbl[, peni] == 'S'),
             cols = c(aminopenicillins, ureidopenicillins, cephalosporins, carbapenems, clox, amcl))
    edit_rsi(to = 'I',
             rule = c(rule_group, rule),
             rows = which(tbl$fullname %like% "^Streptococcus (pyogenes|agalactiae|dysgalactiae|group A|group B|group C|group G)"
                          & tbl[, peni] == 'I'),
             cols = c(aminopenicillins, ureidopenicillins, cephalosporins, carbapenems, clox, amcl))
    edit_rsi(to = 'R',
             rule = c(rule_group, rule),
             rows = which(tbl$fullname %like% "^Streptococcus (pyogenes|agalactiae|dysgalactiae|group A|group B|group C|group G)"
                          & tbl[, peni] == 'R'),
             cols = c(aminopenicillins, ureidopenicillins, cephalosporins, carbapenems, clox, amcl))
  }
  if (!is.na(norf)) {
    edit_rsi(to = 'S',
             rule = c(rule_group, rule),
             rows = which(tbl$fullname %like% "^Streptococcus (pyogenes|agalactiae|dysgalactiae|group A|group B|group C|group G)"
                          & tbl[, norf] == 'S'),
             cols = c(levo, moxi))
  }
  if (!is.na(eryt)) {
    edit_rsi(to = 'S',
             rule = c(rule_group, rule),
             rows = which(tbl$fullname %like% "^Streptococcus (pyogenes|agalactiae|dysgalactiae|group A|group B|group C|group G)"
                          & tbl[, eryt] == 'S'),
             cols = c(azit, clar, roxi))
    edit_rsi(to = 'I',
             rule = c(rule_group, rule),
             rows = which(tbl$fullname %like% "^Streptococcus (pyogenes|agalactiae|dysgalactiae|group A|group B|group C|group G)"
                          & tbl[, eryt] == 'I'),
             cols = c(azit, clar, roxi))
    edit_rsi(to = 'R',
             rule = c(rule_group, rule),
             rows = which(tbl$fullname %like% "^Streptococcus (pyogenes|agalactiae|dysgalactiae|group A|group B|group C|group G)"
                          & tbl[, eryt] == 'R'),
             cols = c(azit, clar, roxi))
  }
  if (!is.na(tetr)) {
    edit_rsi(to = 'S',
             rule = c(rule_group, rule),
             rows = which(tbl$fullname %like% "^Streptococcus (pyogenes|agalactiae|dysgalactiae|group A|group B|group C|group G)"
                          & tbl[, tetr] == 'S'),
             cols = c(doxy, mino))
  }
  if (info == TRUE) {
    txt_ok()
  }
  # Streptococcus pneumoniae ----
  rule <- 'Streptococcus pneumoniae'
  if (info == TRUE) {
    warned <- FALSE
    changed_results <- 0
    cat(rule)
  }
  if (!is.na(peni)) {
    edit_rsi(to = 'S',
             rule = c(rule_group, rule),
             rows = which(tbl$fullname %like% "^Streptococcus pneumoniae"
                          & tbl[, peni] == 'S'),
             cols = c(ampi, amox, amcl, pita))
  }
  if (!is.na(ampi)) {
    edit_rsi(to = 'S',
             rule = c(rule_group, rule),
             rows = which(tbl$fullname %like% "^Streptococcus pneumoniae"
                          & tbl[, ampi] == 'S'),
             cols = c(ampi, amox, amcl, pita))
    edit_rsi(to = 'I',
             rule = c(rule_group, rule),
             rows = which(tbl$fullname %like% "^Streptococcus pneumoniae"
                          & tbl[, ampi] == 'I'),
             cols = c(ampi, amox, amcl, pita))
    edit_rsi(to = 'R',
             rule = c(rule_group, rule),
             rows = which(tbl$fullname %like% "^Streptococcus pneumoniae"
                          & tbl[, ampi] == 'R'),
             cols = c(ampi, amox, amcl, pita))
  }
  if (!is.na(norf)) {
    edit_rsi(to = 'S',
             rule = c(rule_group, rule),
             rows = which(tbl$fullname %like% "^Streptococcus pneumoniae"
                          & tbl[, norf] == 'S'),
             cols = c(levo, moxi))
  }
  if (!is.na(eryt)) {
    edit_rsi(to = 'S',
             rule = c(rule_group, rule),
             rows = which(tbl$fullname %like% "^Streptococcus pneumoniae"
                          & tbl[, eryt] == 'S'),
             cols = c(azit, clar, roxi))
    edit_rsi(to = 'I',
             rule = c(rule_group, rule),
             rows = which(tbl$fullname %like% "^Streptococcus pneumoniae"
                          & tbl[, eryt] == 'I'),
             cols = c(azit, clar, roxi))
    edit_rsi(to = 'R',
             rule = c(rule_group, rule),
             rows = which(tbl$fullname %like% "^Streptococcus pneumoniae"
                          & tbl[, eryt] == 'R'),
             cols = c(azit, clar, roxi))
  }
  if (!is.na(tetr)) {
    edit_rsi(to = 'S',
             rule = c(rule_group, rule),
             rows = which(tbl$fullname %like% "^Streptococcus pneumoniae"
                          & tbl[, tetr] == 'S'),
             cols = c(doxy, mino))
  }
  if (info == TRUE) {
    txt_ok()
  }
  # Viridans group streptococci ----
  rule <- 'Viridans group streptococci'
  if (info == TRUE) {
    warned <- FALSE
    changed_results <- 0
    cat(rule)
  }
  viridans_group <- c("anginosus", "australis", "bovis", "constellatus", "cristatus",
                      "equinus", "gallolyticus", "gordonii", "infantarius", "infantis",
                      "intermedius", "mitis", "mutans", "oligofermentans", "oralis",
                      "parasanguinis", "peroris", "pseudopneumoniae", "salivarius",
                      "sanguinis", "sinensis", "sobrinus", "thermophilus", "vestibularis")
  if (!is.na(ampi)) {
    edit_rsi(to = 'S',
             rule = c(rule_group, rule),
             rows = which(tbl$genus == "Streptococcus" & tbl$species %in% viridans_group
                          & tbl[, ampi] == 'S'),
             cols = c(ampi, amox, amcl, pita))
    edit_rsi(to = 'I',
             rule = c(rule_group, rule),
             rows = which(tbl$genus == "Streptococcus" & tbl$species %in% viridans_group
                          & tbl[, ampi] == 'I'),
             cols = c(ampi, amox, amcl, pita))
    edit_rsi(to = 'R',
             rule = c(rule_group, rule),
             rows = which(tbl$genus == "Streptococcus" & tbl$species %in% viridans_group
                          & tbl[, ampi] == 'R'),
             cols = c(ampi, amox, amcl, pita))
  }
  if (info == TRUE) {
    txt_ok()
  }
  # Haemophilus influenzae ----
  rule <- 'Haemophilus influenzae'
  if (info == TRUE) {
    warned <- FALSE
    changed_results <- 0
    cat(rule)
  }
  if (!is.na(ampi)) {
    edit_rsi(to = 'S',
             rule = c(rule_group, rule),
             rows = which(tbl$fullname %like% "^Haemophilus influenzae"
                          & tbl[, ampi] == 'S'),
             cols = c(amox, pita))
    edit_rsi(to = 'I',
             rule = c(rule_group, rule),
             rows = which(tbl$fullname %like% "^Haemophilus influenzae"
                          & tbl[, ampi] == 'I'),
             cols = c(amox, pita))
    edit_rsi(to = 'R',
             rule = c(rule_group, rule),
             rows = which(tbl$fullname %like% "^Haemophilus influenzae"
                          & tbl[, ampi] == 'R'),
             cols = c(amox, pita))
  }
  if (!is.na(amcl)) {
    edit_rsi(to = 'S',
             rule = c(rule_group, rule),
             rows = which(tbl$fullname %like% "^Haemophilus influenzae"
                          & tbl[, amcl] == 'S'),
             cols = pita)
    edit_rsi(to = 'I',
             rule = c(rule_group, rule),
             rows = which(tbl$fullname %like% "^Haemophilus influenzae"
                          & tbl[, amcl] == 'I'),
             cols = pita)
    edit_rsi(to = 'R',
             rule = c(rule_group, rule),
             rows = which(tbl$fullname %like% "^Haemophilus influenzae"
                          & tbl[, amcl] == 'R'),
             cols = pita)
  }
  if (!is.na(nali)) {
    edit_rsi(to = 'S',
             rule = c(rule_group, rule),
             rows = which(tbl$fullname %like% "^Haemophilus influenzae"
                          & tbl[, nali] == 'S'),
             cols = c(cipr, levo, moxi, oflo))
  }
  if (!is.na(tetr)) {
    edit_rsi(to = 'S',
             rule = c(rule_group, rule),
             rows = which(tbl$fullname %like% "^Haemophilus influenzae"
                          & tbl[, tetr] == 'S'),
             cols = c(doxy, mino))
  }
  if (info == TRUE) {
    txt_ok()
  }
  # Moraxella catarrhalis ----
  rule <- 'Moraxella catarrhalis'
  if (info == TRUE) {
    warned <- FALSE
    changed_results <- 0
    cat(rule)
  }
  if (!is.na(amcl)) {
    edit_rsi(to = 'S',
             rule = c(rule_group, rule),
             rows = which(tbl$fullname %like% "^Moraxella catarrhalis"
                          & tbl[, amcl] == 'S'),
             cols = pita)
    edit_rsi(to = 'I',
             rule = c(rule_group, rule),
             rows = which(tbl$fullname %like% "^Moraxella catarrhalis"
                          & tbl[, amcl] == 'I'),
             cols = pita)
    edit_rsi(to = 'R',
             rule = c(rule_group, rule),
             rows = which(tbl$fullname %like% "^Moraxella catarrhalis"
                          & tbl[, amcl] == 'R'),
             cols = pita)
  }
  if (!is.na(nali)) {
    edit_rsi(to = 'S',
             rule = c(rule_group, rule),
             rows = which(tbl$fullname %like% "^Moraxella catarrhalis"
                          & tbl[, nali] == 'S'),
             cols = c(cipr, levo, moxi, oflo))
  }
  if (!is.na(eryt)) {
    edit_rsi(to = 'S',
             rule = c(rule_group, rule),
             rows = which(tbl$fullname %like% "^Moraxella catarrhalis"
                          & tbl[, eryt] == 'S'),
             cols = c(azit, clar, roxi))
    edit_rsi(to = 'I',
             rule = c(rule_group, rule),
             rows = which(tbl$fullname %like% "^Moraxella catarrhalis"
                          & tbl[, eryt] == 'I'),
             cols = c(azit, clar, roxi))
    edit_rsi(to = 'R',
             rule = c(rule_group, rule),
             rows = which(tbl$fullname %like% "^Moraxella catarrhalis"
                          & tbl[, eryt] == 'R'),
             cols = c(azit, clar, roxi))
  }
  if (!is.na(tetr)) {
    edit_rsi(to = 'S',
             rule = c(rule_group, rule),
             rows = which(tbl$fullname %like% "^Moraxella catarrhalis"
                          & tbl[, tetr] == 'S'),
             cols = c(doxy, mino))
  }
  if (info == TRUE) {
    txt_ok()
  }
  # Anaerobic Gram positives ----
  rule <- 'Anaerobic Gram positives'
  if (info == TRUE) {
    warned <- FALSE
    changed_results <- 0
    cat(rule)
  }
  if (!is.na(peni)) {
    edit_rsi(to = 'S',
             rule = c(rule_group, rule),
             rows = which((tbl$genus %in% c("Clostridium", "Actinomyces", "Propionibacterium",
                                            "Bifidobacterium", "Eggerthella", "Eubacterium",
                                            "Lactobacillus ", "Actinomyces")
                           | tbl$fullname %like% "^Propionibacterium acnes")
                          & tbl[, peni] == 'S'),
             cols = c(ampi, amox, pita, tica))
    edit_rsi(to = 'I',
             rule = c(rule_group, rule),
             rows = which((tbl$genus %in% c("Clostridium", "Actinomyces", "Propionibacterium",
                                            "Bifidobacterium", "Eggerthella", "Eubacterium",
                                            "Lactobacillus ", "Actinomyces")
                           | tbl$fullname %like% "^Propionibacterium acnes")
                          & tbl[, peni] == 'I'),
             cols = c(ampi, amox, pita, tica))
    edit_rsi(to = 'R',
             rule = c(rule_group, rule),
             rows = which((tbl$genus %in% c("Clostridium", "Actinomyces", "Propionibacterium",
                                            "Bifidobacterium", "Eggerthella", "Eubacterium",
                                            "Lactobacillus ", "Actinomyces")
                           | tbl$fullname %like% "^Propionibacterium acnes")
                          & tbl[, peni] == 'R'),
             cols = c(ampi, amox, pita, tica))
  }
  if (info == TRUE) {
    txt_ok()
  }
  # Anaerobic Gram negatives ----
  rule <- 'Anaerobic Gram negatives'
  if (info == TRUE) {
    warned <- FALSE
    changed_results <- 0
    cat(rule)
  }
  if (!is.na(peni)) {
    edit_rsi(to = 'S',
             rule = c(rule_group, rule),
             rows = which(tbl$genus %in% c("Bacteroides", "Prevotella", "Porphyromonas",
                                           "Fusobacterium", "Bilophila ", "Mobiluncus")
                          & tbl[, peni] == 'S'),
             cols = c(ampi, amox, pita, tica))
    edit_rsi(to = 'I',
             rule = c(rule_group, rule),
             rows = which(tbl$genus %in% c("Bacteroides", "Prevotella", "Porphyromonas",
                                           "Fusobacterium", "Bilophila ", "Mobiluncus")
                          & tbl[, peni] == 'I'),
             cols = c(ampi, amox, pita, tica))
    edit_rsi(to = 'R',
             rule = c(rule_group, rule),
             rows = which(tbl$genus %in% c("Bacteroides", "Prevotella", "Porphyromonas",
                                           "Fusobacterium", "Bilophila ", "Mobiluncus")
                          & tbl[, peni] == 'R'),
             cols = c(ampi, amox, pita, tica))
  }
  if (info == TRUE) {
    txt_ok()
  }
  # Pasteurella multocida ----
  rule <- 'Pasteurella multocida'
  if (info == TRUE) {
    warned <- FALSE
    changed_results <- 0
    cat(rule)
  }
  if (!is.na(peni)) {
    edit_rsi(to = 'S',
             rule = c(rule_group, rule),
             rows = which(tbl$fullname %like% "^Pasteurella multocida"
                          & tbl[, peni] == 'S'),
             cols = c(ampi, amox))
    edit_rsi(to = 'I',
             rule = c(rule_group, rule),
             rows = which(tbl$fullname %like% "^Pasteurella multocida"
                          & tbl[, peni] == 'I'),
             cols = c(ampi, amox))
    edit_rsi(to = 'R',
             rule = c(rule_group, rule),
             rows = which(tbl$fullname %like% "^Pasteurella multocida"
                          & tbl[, peni] == 'R'),
             cols = c(ampi, amox))
  }
  if (info == TRUE) {
    txt_ok()
  }
  # Campylobacter jejuni ----
  rule <- 'Campylobacter jejuni'
  if (info == TRUE) {
    warned <- FALSE
    changed_results <- 0
    cat(rule)
  }
  if (!is.na(eryt)) {
    edit_rsi(to = 'S',
             rule = c(rule_group, rule),
             rows = which(tbl$fullname %like% "^Campylobacter jejuni"
                          & tbl[, eryt] == 'S'),
             cols = c(azit, clar))
    edit_rsi(to = 'I',
             rule = c(rule_group, rule),
             rows = which(tbl$fullname %like% "^Campylobacter jejuni"
                          & tbl[, eryt] == 'I'),
             cols = c(azit, clar))
    edit_rsi(to = 'R',
             rule = c(rule_group, rule),
             rows = which(tbl$fullname %like% "^Campylobacter jejuni"
                          & tbl[, eryt] == 'R'),
             cols = c(azit, clar))
  }
  if (!is.na(tetr)) {
    edit_rsi(to = 'S',
             rule = c(rule_group, rule),
             rows = which(tbl$fullname %like% "^Campylobacter jejuni"
                          & tbl[, tetr] == 'S'),
             cols = doxy)
    edit_rsi(to = 'I',
             rule = c(rule_group, rule),
             rows = which(tbl$fullname %like% "^Campylobacter jejuni"
                          & tbl[, tetr] == 'I'),
             cols = doxy)
    edit_rsi(to = 'R',
             rule = c(rule_group, rule),
             rows = which(tbl$fullname %like% "^Campylobacter jejuni"
                          & tbl[, tetr] == 'R'),
             cols = doxy)
  }
  if (info == TRUE) {
    txt_ok()
  }
  # Aerococcus sanguinicola/urinae ----
  rule <- 'Aerococcus sanguinicola/urinae'
  if (info == TRUE) {
    warned <- FALSE
    changed_results <- 0
    cat(rule)
  }
  if (!is.na(cipr)) {
    edit_rsi(to = 'S',
             rule = c(rule_group, rule),
             rows = which(tbl$fullname %like% "^Aerococcus (sanguinicola|urinae)"
                          & tbl[, cipr] == 'S'),
             cols = levo)
    edit_rsi(to = 'I',
             rule = c(rule_group, rule),
             rows = which(tbl$fullname %like% "^Aerococcus (sanguinicola|urinae)"
                          & tbl[, cipr] == 'I'),
             cols = levo)
    edit_rsi(to = 'R',
             rule = c(rule_group, rule),
             rows = which(tbl$fullname %like% "^Aerococcus (sanguinicola|urinae)"
                          & tbl[, cipr] == 'R'),
             cols = levo)
  }
  if (info == TRUE) {
    txt_ok()
  }
  # Kingella kingae ----
  rule <- 'Kingella kingae'
  if (info == TRUE) {
    warned <- FALSE
    changed_results <- 0
    cat(rule)
  }
  if (!is.na(peni)) {
    edit_rsi(to = 'S',
             rule = c(rule_group, rule),
             rows = which(tbl$fullname %like% "^Kingella kingae"
                          & tbl[, peni] == 'S'),
             cols = c(ampi, amox))
    edit_rsi(to = 'I',
             rule = c(rule_group, rule),
             rows = which(tbl$fullname %like% "^Kingella kingae"
                          & tbl[, peni] == 'I'),
             cols = c(ampi, amox))
    edit_rsi(to = 'R',
             rule = c(rule_group, rule),
             rows = which(tbl$fullname %like% "^Kingella kingae"
                          & tbl[, peni] == 'R'),
             cols = c(ampi, amox))
  }
  if (!is.na(eryt)) {
    edit_rsi(to = 'S',
             rule = c(rule_group, rule),
             rows = which(tbl$fullname %like% "^Kingella kingae"
                          & tbl[, eryt] == 'S'),
             cols = c(azit, clar))
    edit_rsi(to = 'I',
             rule = c(rule_group, rule),
             rows = which(tbl$fullname %like% "^Kingella kingae"
                          & tbl[, eryt] == 'I'),
             cols = c(azit, clar))
    edit_rsi(to = 'R',
             rule = c(rule_group, rule),
             rows = which(tbl$fullname %like% "^Kingella kingae"
                          & tbl[, eryt] == 'R'),
             cols = c(azit, clar))
  }
  if (!is.na(tetr)) {
    edit_rsi(to = 'S',
             rule = c(rule_group, rule),
             rows = which(tbl$fullname %like% "^Kingella kingae"
                          & tbl[, tetr] == 'S'),
             cols = doxy)
    edit_rsi(to = 'I',
             rule = c(rule_group, rule),
             rows = which(tbl$fullname %like% "^Kingella kingae"
                          & tbl[, tetr] == 'I'),
             cols = doxy)
    edit_rsi(to = 'R',
             rule = c(rule_group, rule),
             rows = which(tbl$fullname %like% "^Kingella kingae"
                          & tbl[, tetr] == 'R'),
             cols = doxy)
  }
  if (info == TRUE) {
    txt_ok()
  }

  # EXPERT RULES AND INTRINSIC RESISTANCE -----------------------------------

  if (info == TRUE) {
    cat(bold(paste0('\nEUCAST Expert Rules, Intrinsic Resistance and Exceptional Phenotypes (v', EUCAST_VERSION_EXPERT_RULES, ')\n')))
  }
  rule_group <- "Expert Rules"

  # Table 1: Intrinsic resistance in Enterobacteriaceae ----
  rule <- 'Table 1:  Intrinsic resistance in Enterobacteriaceae'
  if (info == TRUE) {
    warned <- FALSE
    changed_results <- 0
    cat(rule)
  }
  # Intrinsic R for this group
  edit_rsi(to = 'R',
           rule = c(rule_group, rule),
           rows = which(tbl$family == 'Enterobacteriaceae'),
           cols = c(peni, glycopeptides, fusi, macrolides, linc, streptogramins, rifa, dapt, line))
  # Citrobacter
  edit_rsi(to = 'R',
           rule = c(rule_group, rule),
           rows = which(tbl$fullname %like% '^Citrobacter (koseri|amalonaticus|sedlakii|farmeri|rodentium)'),
           cols = c(aminopenicillins, tica))
  edit_rsi(to = 'R',
           rule = c(rule_group, rule),
           rows = which(tbl$fullname %like% '^Citrobacter (freundii|braakii|murliniae|werkmanii|youngae)'),
           cols = c(aminopenicillins, amcl, czol, cfox))
  # Enterobacter
  edit_rsi(to = 'R',
           rule = c(rule_group, rule),
           rows = which(tbl$fullname %like% '^Enterobacter cloacae'),
           cols = c(aminopenicillins, amcl, czol, cfox))
  edit_rsi(to = 'R',
           rule = c(rule_group, rule),
           rows = which(tbl$fullname %like% '^Enterobacter aerogenes'),
           cols = c(aminopenicillins, amcl, czol, cfox))
  # Escherichia
  edit_rsi(to = 'R',
           rule = c(rule_group, rule),
           rows = which(tbl$fullname %like% '^Escherichia hermanni'),
           cols = c(aminopenicillins, tica))
  # Hafnia
  edit_rsi(to = 'R',
           rule = c(rule_group, rule),
           rows = which(tbl$fullname %like% '^Hafnia alvei'),
           cols = c(aminopenicillins, amcl, czol, cfox))
  # Klebsiella
  edit_rsi(to = 'R',
           rule = c(rule_group, rule),
           rows = which(tbl$fullname %like% '^Klebsiella'),
           cols = c(aminopenicillins, tica))
  # Morganella / Proteus
  edit_rsi(to = 'R',
           rule = c(rule_group, rule),
           rows = which(tbl$fullname %like% '^Morganella morganii'),
           cols = c(aminopenicillins, amcl, czol, tetracyclines, polymyxins, nitr))
  edit_rsi(to = 'R',
           rule = c(rule_group, rule),
           rows = which(tbl$fullname %like% '^Proteus mirabilis'),
           cols = c(tetracyclines, tige, polymyxins, nitr))
  edit_rsi(to = 'R',
           rule = c(rule_group, rule),
           rows = which(tbl$fullname %like% '^Proteus penneri'),
           cols = c(aminopenicillins, czol, cfur, tetracyclines, tige, polymyxins, nitr))
  edit_rsi(to = 'R',
           rule = c(rule_group, rule),
           rows = which(tbl$fullname %like% '^Proteus vulgaris'),
           cols = c(aminopenicillins, czol, cfur, tetracyclines, tige, polymyxins, nitr))
  # Providencia
  edit_rsi(to = 'R',
           rule = c(rule_group, rule),
           rows = which(tbl$fullname %like% '^Providencia rettgeri'),
           cols = c(aminopenicillins, amcl, czol, cfur, tetracyclines, tige, polymyxins, nitr))
  edit_rsi(to = 'R',
           rule = c(rule_group, rule),
           rows = which(tbl$fullname %like% '^Providencia stuartii'),
           cols = c(aminopenicillins, amcl, czol, cfur, tetracyclines, tige, polymyxins, nitr))
  # Raoultella
  edit_rsi(to = 'R',
           rule = c(rule_group, rule),
           rows = which(tbl$fullname %like% '^Raoultella'),
           cols = c(aminopenicillins, tica))
  # Serratia
  edit_rsi(to = 'R',
           rule = c(rule_group, rule),
           rows = which(tbl$fullname %like% '^Serratia marcescens'),
           cols = c(aminopenicillins, amcl, czol, cfox, cfur, tetracyclines[tetracyclines != na.rm(mino)], polymyxins, nitr))
  # Yersinia
  edit_rsi(to = 'R',
           rule = c(rule_group, rule),
           rows = which(tbl$fullname %like% '^Yersinia enterocolitica'),
           cols = c(aminopenicillins, amcl, tica, czol, cfox))
  edit_rsi(to = 'R',
           rule = c(rule_group, rule),
           rows = which(tbl$fullname %like% '^Yersinia pseudotuberculosis'),
           cols = c(poly, coli))
  if (info == TRUE) {
    txt_ok()
  }

  # Table 2: Intrinsic resistance in non-fermentative Gram-negative bacteria ----
  rule <- 'Table 2:  Intrinsic resistance in non-fermentative Gram-negative bacteria'
  if (info == TRUE) {
    warned <- FALSE
    changed_results <- 0
    cat(rule)
  }
  # Intrinsic R for this group
  edit_rsi(to = 'R',
           rule = c(rule_group, rule),
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
           rule = c(rule_group, rule),
           rows = which(tbl$fullname %like% '^Acinetobacter (baumannii|pittii|nosocomialis|calcoaceticus)'),
           cols = c(aminopenicillins, amcl, czol, cfot, cftr, aztr, erta, trim, fosf, tetracyclines[tetracyclines != na.rm(mino)]))
  # Achromobacter
  edit_rsi(to = 'R',
           rule = c(rule_group, rule),
           rows = which(tbl$fullname %like% '^Achromobacter (xylosoxydans|xylosoxidans)'),
           cols = c(aminopenicillins, czol, cfot, cftr, erta))
  # Burkholderia
  edit_rsi(to = 'R',
           rule = c(rule_group, rule),
           # onder 'Burkholderia cepacia complex' vallen deze species allemaal: PMID 16217180.
           rows = which(tbl$fullname %like% '^Burkholderia (cepacia|multivorans|cenocepacia|stabilis|vietnamiensis|dolosa|ambifaria|anthina|pyrrocinia|ubonensis)'),
           cols = c(aminopenicillins, amcl, tica, pita, czol, cfot, cftr, aztr, erta, cipr, chlo, aminoglycosides, trim, fosf, polymyxins))
  # Elizabethkingia
  edit_rsi(to = 'R',
           rule = c(rule_group, rule),
           rows = which(tbl$fullname %like% '^Elizabethkingia meningoseptic(a|um)'),
           cols = c(aminopenicillins, amcl, tica, czol, cfot, cftr, cfta, cfep, aztr, erta, imip, mero, polymyxins))
  # Ochrobactrum
  edit_rsi(to = 'R',
           rule = c(rule_group, rule),
           rows = which(tbl$fullname %like% '^Ochrobactrum anthropi'),
           cols = c(aminopenicillins, amcl, tica, pita, czol, cfot, cftr, cfta, cfep, aztr, erta))
  # Pseudomonas
  edit_rsi(to = 'R',
           rule = c(rule_group, rule),
           rows = which(tbl$fullname %like% '^Pseudomonas aeruginosa'),
           cols = c(aminopenicillins, amcl, czol, cfot, cftr, erta, chlo, kana, neom, trim, trsu, tetracyclines, tige))
  # Stenotrophomonas
  edit_rsi(to = 'R',
           rule = c(rule_group, rule),
           rows = which(tbl$fullname %like% '^Stenotrophomonas maltophilia'),
           cols = c(aminopenicillins, amcl, tica, pita, czol, cfot, cftr, cfta, aztr, erta, imip, mero, aminoglycosides, trim, fosf, tetr))
  if (info == TRUE) {
    txt_ok()
  }

  # Table 3: Intrinsic resistance in other Gram-negative bacteria ----
  rule <- 'Table 3:  Intrinsic resistance in other Gram-negative bacteria'
  if (info == TRUE) {
    warned <- FALSE
    changed_results <- 0
    cat(rule)
  }
  # Intrinsic R for this group
  edit_rsi(to = 'R',
           rule = c(rule_group, rule),
           rows = which(tbl$genus %in% c('Haemophilus',
                                         'Moraxella',
                                         'Neisseria',
                                         'Campylobacter')),
           cols = c(glycopeptides, linc, dapt, line))
  # Haemophilus
  edit_rsi(to = 'R',
           rule = c(rule_group, rule),
           rows = which(tbl$fullname %like% '^Haemophilus influenzae'),
           cols = c(fusi, streptogramins))
  # Moraxella
  edit_rsi(to = 'R',
           rule = c(rule_group, rule),
           rows = which(tbl$fullname %like% '^Moraxella catarrhalis'),
           cols = trim)
  # Neisseria
  edit_rsi(to = 'R',
           rule = c(rule_group, rule),
           rows = which(tbl$genus == 'Neisseria'),
           cols = trim)
  # Campylobacter
  edit_rsi(to = 'R',
           rule = c(rule_group, rule),
           rows = which(tbl$fullname %like% '^Campylobacter fetus'),
           cols = c(fusi, streptogramins, trim, nali))
  edit_rsi(to = 'R',
           rule = c(rule_group, rule),
           rows = which(tbl$fullname %like% '^Campylobacter (jejuni|coli)'),
           cols = c(fusi, streptogramins, trim))
  if (info == TRUE) {
    txt_ok()
  }

  # Table 4: Intrinsic resistance in Gram-positive bacteria ----
  rule <- 'Table 4:  Intrinsic resistance in Gram-positive bacteria'
  if (info == TRUE) {
    warned <- FALSE
    changed_results <- 0
    cat(rule)
  }
  # Intrinsic R for this group
  edit_rsi(to = 'R',
           rule = c(rule_group, rule),
           rows = which(tbl$gramstain == "Gram positive"),
           cols = c(aztr, polymyxins, nali))
  # Staphylococcus
  edit_rsi(to = 'R',
           rule = c(rule_group, rule),
           rows = which(tbl$fullname %like% '^Staphylococcus saprophyticus'),
           cols = c(fusi, cfta, fosf, novo))
  edit_rsi(to = 'R',
           rule = c(rule_group, rule),
           rows = which(tbl$fullname %like% '^Staphylococcus (cohnii|xylosus)'),
           cols = c(cfta, novo))
  edit_rsi(to = 'R',
           rule = c(rule_group, rule),
           rows = which(tbl$fullname %like% '^Staphylococcus capitis'),
           cols = c(cfta, fosf))
  edit_rsi(to = 'R',
           rule = c(rule_group, rule),
           rows = which(tbl$fullname %like% '^Staphylococcus (aureus|epidermidis|coagulase negatief|hominis|haemolyticus|intermedius|pseudointermedius)'),
           cols = cfta)
  # Streptococcus
  # rule 4.5
  edit_rsi(to = 'R',
           rule = c(rule_group, rule),
           rows = which(tbl$genus == 'Streptococcus'),
           cols = c(fusi, aminoglycosides))
  # Enterococcus
  edit_rsi(to = 'R',
           rule = c(rule_group, rule),
           rows = which(tbl$fullname %like% '^Enterococcus faecalis'),
           cols = c(fusi, cfta, cephalosporins[cephalosporins != na.rm(cfta)], aminoglycosides, macrolides, clin, qida, trim, trsu))
  edit_rsi(to = 'R',
           rule = c(rule_group, rule),
           rows = which(tbl$fullname %like% '^Enterococcus (gallinarum|casseliflavus)'),
           cols = c(fusi, cfta, cephalosporins[cephalosporins != na.rm(cfta)], aminoglycosides, macrolides, clin, qida, vanc, trim, trsu))
  edit_rsi(to = 'R',
           rule = c(rule_group, rule),
           rows = which(tbl$fullname %like% '^Enterococcus faecium'),
           cols = c(fusi, cfta, cephalosporins[cephalosporins != na.rm(cfta)], aminoglycosides, macrolides, trim, trsu))
  # Corynebacterium
  edit_rsi(to = 'R',
           rule = c(rule_group, rule),
           rows = which(tbl$genus == 'Corynebacterium'),
           cols = fosf)
  # Listeria
  edit_rsi(to = 'R',
           rule = c(rule_group, rule),
           rows = which(tbl$fullname %like% '^Listeria monocytogenes'),
           cols = c(cfta, cephalosporins[cephalosporins != na.rm(cfta)]))
  # other
  edit_rsi(to = 'R',
           rule = c(rule_group, rule),
           rows = which(tbl$genus %in% c('Leuconostoc', 'Pediococcus')),
           cols = glycopeptides)
  edit_rsi(to = 'R',
           rule = c(rule_group, rule),
           rows = which(tbl$genus == 'Lactobacillus'),
           cols = glycopeptides)
  edit_rsi(to = 'R',
           rule = c(rule_group, rule),
           rows = which(tbl$fullname %like% '^Clostridium (ramosum|innocuum)'),
           cols = vanc)
  if (info == TRUE) {
    txt_ok()
  }

  # Table 8: Interpretive rules for B-lactam agents and Gram-positive cocci ----
  rule <- 'Table 8:  Interpretive rules for B-lactam agents and Gram-positive cocci'
  if (info == TRUE) {
    warned <- FALSE
    changed_results <- 0
    cat(rule)
  }
  # rule 8.3
  if (!is.na(peni)) {
    edit_rsi(to = 'S',
             rule = c(rule_group, rule),
             rows = which(tbl$fullname %like% '^Streptococcus (pyogenes|agalactiae|dysgalactiae|group A|group B|group C|group G)'
                          & tbl[, peni] == 'S'),
             cols = c(aminopenicillins, cephalosporins, carbapenems))
  }
  # rule 8.6
  if (!is.na(ampi)) {
    edit_rsi(to = 'R',
             rule = c(rule_group, rule),
             rows = which(tbl$genus == 'Enterococcus'
                          & tbl[, ampi] == 'R'),
             cols = c(ureidopenicillins, carbapenems))
  }
  if (!is.na(amox)) {
    edit_rsi(to = 'R',
             rule = c(rule_group, rule),
             rows = which(tbl$genus == 'Enterococcus'
                          & tbl[, amox] == 'R'),
             cols = c(ureidopenicillins, carbapenems))
  }
  if (info == TRUE) {
    txt_ok()
  }

  # Table 9: Interpretive rules for B-lactam agents and Gram-negative rods ----
  rule <- 'Table 9:  Interpretive rules for B-lactam agents and Gram-negative rods'
  if (info == TRUE) {
    warned <- FALSE
    changed_results <- 0
    cat(rule)
  }
  # rule 9.3
  if (!is.na(tica) & !is.na(pita)) {
    edit_rsi(to = 'R',
             rule = c(rule_group, rule),
             rows = which(tbl$family == 'Enterobacteriaceae'
                          & tbl[, tica] == 'R'
                          & tbl[, pita] == 'S'),
             cols = pita)
  }
  if (info == TRUE) {
    txt_ok()
  }

  # Table 10: Interpretive rules for B-lactam agents and other Gram-negative bacteria ----
  rule <- 'Table 10: Interpretive rules for B-lactam agents and other Gram-negative bacteria'
  if (info == TRUE) {
    warned <- FALSE
    changed_results <- 0
    cat(rule)
  }
  # rule 10.2
  # if (!is.na(ampi)) {
  # you should know first if the are B-lactamase positive, so do not run for now
  # edit_rsi(to = 'R',
  #          rule = c(rule_group, rule),
  #          rows = which(tbl$fullname %like% '^Haemophilus influenza'
  #                       & tbl[, ampi] == 'R'),
  #          cols = c(ampi, amox, amcl, pita, cfur))
  # }
  if (info == TRUE) {
    txt_ok()
  }

  # Table 11: Interpretive rules for macrolides, lincosamides, and streptogramins ----
  rule <- 'Table 11: Interpretive rules for macrolides, lincosamides, and streptogramins'
  if (info == TRUE) {
    warned <- FALSE
    changed_results <- 0
    cat(rule)
  }
  # rule 11.1
  if (!is.na(eryt)) {
    edit_rsi(to = 'S',
             rule = c(rule_group, rule),
             rows = which(tbl[, eryt] == 'S'),
             cols = c(azit, clar))
    edit_rsi(to = 'I',
             rule = c(rule_group, rule),
             rows = which(tbl[, eryt] == 'I'),
             cols = c(azit, clar))
    edit_rsi(to = 'R',
             rule = c(rule_group, rule),
             rows = which(tbl[, eryt] == 'R'),
             cols = c(azit, clar))
  }
  if (info == TRUE) {
    txt_ok()
  }

  # Table 12: Interpretive rules for aminoglycosides ----
  rule <- 'Table 12: Interpretive rules for aminoglycosides'
  if (info == TRUE) {
    warned <- FALSE
    changed_results <- 0
    cat(rule)
  }
  # rule 12.2
  if (!is.na(tobr)) {
    edit_rsi(to = 'R',
             rule = c(rule_group, rule),
             rows = which(tbl$genus == 'Staphylococcus'
                          & tbl[, tobr] == 'R'),
             cols = c(kana, amik))
  }
  # rule 12.3
  if (!is.na(gent)) {
    edit_rsi(to = 'R',
             rule = c(rule_group, rule),
             rows = which(tbl$genus == 'Staphylococcus'
                          & tbl[, gent] == 'R'),
             cols = aminoglycosides)
  }
  # rule 12.8
  if (!is.na(gent) & !is.na(tobr)) {
    edit_rsi(to = 'R',
             rule = c(rule_group, rule),
             rows = which(tbl$family == 'Enterobacteriaceae'
                          & tbl[, gent] == 'I'
                          & tbl[, tobr] == 'S'),
             cols = gent)
  }
  # rule 12.9
  if (!is.na(gent) & !is.na(tobr)) {
    edit_rsi(to = 'R',
             rule = c(rule_group, rule),
             rows = which(tbl$family == 'Enterobacteriaceae'
                          & tbl[, tobr] == 'I'
                          & tbl[, gent] == 'R'),
             cols = tobr)
  }
  if (info == TRUE) {
    txt_ok()
  }


  # Table 13: Interpretive rules for quinolones ----
  rule <- 'Table 13: Interpretive rules for quinolones'
  if (info == TRUE) {
    warned <- FALSE
    changed_results <- 0
    cat(rule)
  }
  # rule 13.2
  if (!is.na(moxi)) {
    edit_rsi(to = 'R',
             rule = c(rule_group, rule),
             rows = which(tbl$genus == 'Staphylococcus'
                          & tbl[, moxi] == 'R'),
             cols = fluoroquinolones)
  }
  # rule 13.4
  if (!is.na(moxi)) {
    edit_rsi(to = 'R',
             rule = c(rule_group, rule),
             rows = which(tbl$fullname %like% '^Streptococcus pneumoniae'
                          & tbl[, moxi] == 'R'),
             cols = fluoroquinolones)
  }
  # rule 13.5
  if (!is.na(cipr)) {
    edit_rsi(to = 'R',
             rule = c(rule_group, rule),
             rows = which(tbl$family == 'Enterobacteriaceae'
                          & tbl[, cipr] == 'R'),
             cols = fluoroquinolones)
  }
  # rule 13.8
  if (!is.na(cipr)) {
    edit_rsi(to = 'R',
             rule = c(rule_group, rule),
             rows = which(tbl$fullname %like% '^Neisseria gonorrhoeae'
                          & tbl[, cipr] == 'R'),
             cols = fluoroquinolones)
  }
  if (info == TRUE) {
    txt_ok()
  }

  # OTHER RULES -------------------------------------------------------------

  if (info == TRUE) {
    cat(bold('\nOther rules\n'))
  }
  rule_group <- "Other rules"

  rule <- 'Non-EUCAST: ampicillin = R where amoxicillin/clav acid = R'
  if (info == TRUE) {
    warned <- FALSE
    changed_results <- 0
    cat(rule)
  }
  if (!is.na(amcl)) {
    edit_rsi(to = 'R',
             rule = c(rule_group, rule),
             rows = which(tbl[, amcl] == 'R'),
             cols = ampi)
  }
  if (info == TRUE) {
    txt_ok()
  }
  rule <- 'Non-EUCAST: trimethoprim = R where trimethoprim/sulfa = R'
  if (info == TRUE) {
    warned <- FALSE
    changed_results <- 0
    cat(rule)
  }
  if (!is.na(trsu)) {
    edit_rsi(to = 'R',
             rule = c(rule_group, rule),
             rows = which(tbl[, trsu] == 'R'),
             cols = trim)
  }
  if (info == TRUE) {
    txt_ok()
  }
  rule <- 'Non-EUCAST: amoxicillin/clav acid = S where ampicillin = S'
  if (info == TRUE) {
    warned <- FALSE
    changed_results <- 0
    cat(rule)
  }
  if (!is.na(amcl)) {
    edit_rsi(to = 'S',
             rule = c(rule_group, rule),
             rows = which(tbl[, ampi] == 'S'),
             cols = amcl)
  }
  if (info == TRUE) {
    txt_ok()
  }
  rule <- 'Non-EUCAST: trimethoprim/sulfa = S where trimethoprim = S'
  if (info == TRUE) {
    warned <- FALSE
    changed_results <- 0
    cat(rule)
  }
  if (!is.na(trsu)) {
    edit_rsi(to = 'S',
             rule = c(rule_group, rule),
             rows = which(tbl[, trim] == 'S'),
             cols = trsu)
  }
  if (info == TRUE) {
    txt_ok()
  }

  # restore old col_mo values if needed
  if (!is.null(col_mo_original)) {
    tbl_original[, col_mo] <- col_mo_original
  }

  if (info == TRUE) {
    cat(bold('\n=> EUCAST rules affected',
             amount_affected_rows %>% length() %>% format(big.mark = ","),
             'out of', nrow(tbl_original) %>% format(big.mark = ","),
             'rows -- changed',
             amount_changed %>% format(big.mark = ","), 'test results.\n\n'))
  }

  if (verbose_info != "") {
    message("Verbose information:", verbose_info)
  }

  tbl_original
}

#' @rdname EUCAST
#' @export
interpretive_reading <- function(...) {
  EUCAST_rules(...)
}
