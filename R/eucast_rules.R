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

#' EUCAST rules
#'
#' Apply susceptibility rules as defined by the European Committee on Antimicrobial Susceptibility Testing (EUCAST, \url{http://eucast.org}), see \emph{Source}. This includes (1) expert rules, (2) intrinsic resistance and (3) inferred resistance as defined in their breakpoint tables.
#' @param tbl table with antibiotic columns, like e.g. \code{amox} and \code{amcl}
#' @param info print progress
#' @param rules a character vector that specifies which rules should be applied - one or more of \code{c("breakpoints", "expert", "other", "all")}
#' @param verbose a logical to indicate whether extensive info should be returned as a \code{data.frame} with info about which rows and columns are effected
#' @param amcl,amik,amox,ampi,azit,azlo,aztr,cefa,cfep,cfot,cfox,cfra,cfta,cftr,cfur,chlo,cipr,clar,clin,clox,coli,czol,dapt,doxy,erta,eryt,fosf,fusi,gent,imip,kana,levo,linc,line,mero,mezl,mino,moxi,nali,neom,neti,nitr,norf,novo,oflo,oxac,peni,pipe,pita,poly,pris,qida,rifa,roxi,siso,teic,tetr,tica,tige,tobr,trim,trsu,vanc column name of an antibiotic, see Antibiotics
#' @param ... parameters that are passed on to \code{eucast_rules}
#' @inheritParams first_isolate
#' @section Antibiotics:
#' To define antibiotics column names, leave as it is to determine it automatically with \code{\link{guess_ab}} or input a text (case-insensitive) or use \code{NULL} to skip a column (e.g. \code{tica = NULL}). Non-existing columns will anyway be skipped with a warning.
#'
#' Abbrevations of the column containing antibiotics in the form: \strong{abbreviation}: generic name (\emph{ATC code})
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
#'  \strong{pipe}: piperacillin (\emph{J01CA12}),
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
#' @rdname eucast_rules
#' @export
#' @importFrom dplyr %>% select pull mutate_at vars
#' @importFrom crayon bold bgGreen bgYellow bgRed black green blue italic strip_style
#' @return The input of \code{tbl}, possibly with edited values of antibiotics. Or, if \code{verbose = TRUE}, a \code{data.frame} with verbose info.
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
#' @inheritSection AMR Read more on our website!
#' @examples
#' a <- eucast_rules(septic_patients)
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
#'                 peni = "S",       # Benzylpenicillin
#'                 cfox = "S",       # Cefoxitin
#'                 stringsAsFactors = FALSE)
#'
#' a
#' #                       mo vanc amox coli cfta cfur peni cfox
#' # 1  Staphylococcus aureus    -    -    -    -    -    S    S
#' # 2  Enterococcus faecalis    -    -    -    -    -    S    S
#' # 3       Escherichia coli    -    -    -    -    -    S    S
#' # 4  Klebsiella pneumoniae    -    -    -    -    -    S    S
#' # 5 Pseudomonas aeruginosa    -    -    -    -    -    S    S
#'
#' b <- eucast_rules(a, "mo") # 18 results are forced as R or S
#'
#' b
#' #                       mo vanc amox coli cfta cfur peni cfox
#' # 1  Staphylococcus aureus    -    S    R    R    S    S    S
#' # 2  Enterococcus faecalis    -    -    R    R    R    S    R
#' # 3       Escherichia coli    R    -    -    -    -    R    S
#' # 4  Klebsiella pneumoniae    R    R    -    -    -    R    S
#' # 5 Pseudomonas aeruginosa    R    R    -    -    R    R    R
eucast_rules <- function(tbl,
                         col_mo = NULL,
                         info = TRUE,
                         rules = c("breakpoints", "expert", "other", "all"),
                         verbose = FALSE,
                         amcl = guess_ab(),
                         amik = guess_ab(),
                         amox = guess_ab(),
                         ampi = guess_ab(),
                         azit = guess_ab(),
                         azlo = guess_ab(),
                         aztr = guess_ab(),
                         cefa = guess_ab(),
                         cfep = guess_ab(),
                         cfot = guess_ab(),
                         cfox = guess_ab(),
                         cfra = guess_ab(),
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
                         mezl = guess_ab(),
                         mino = guess_ab(),
                         moxi = guess_ab(),
                         nali = guess_ab(),
                         neom = guess_ab(),
                         neti = guess_ab(),
                         nitr = guess_ab(),
                         norf = guess_ab(),
                         novo = guess_ab(),
                         oflo = guess_ab(),
                         oxac = guess_ab(),
                         peni = guess_ab(),
                         pipe = guess_ab(),
                         pita = guess_ab(),
                         poly = guess_ab(),
                         pris = guess_ab(),
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

  EUCAST_VERSION_BREAKPOINTS <- "8.1, 2018"
  EUCAST_VERSION_EXPERT_RULES <- "3.1, 2016"

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

  if (!all(rules %in% c("breakpoints", "expert", "other", "all"))) {
    stop("Parameter `rules` must be one or more of:  'breakpoints', 'expert', 'other', 'all'.")
  }

  if (is.null(col_mo)) {
    stop("Parameter `col_mo` must be set")
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
  if (identical(amcl, as.name("guess_ab"))) { amcl <- guess_ab(tbl, "amcl", verbose = verbose) }
  if (identical(amik, as.name("guess_ab"))) { amik <- guess_ab(tbl, "amik", verbose = verbose) }
  if (identical(amox, as.name("guess_ab"))) { amox <- guess_ab(tbl, "amox", verbose = verbose) }
  if (identical(ampi, as.name("guess_ab"))) { ampi <- guess_ab(tbl, "ampi", verbose = verbose) }
  if (identical(azit, as.name("guess_ab"))) { azit <- guess_ab(tbl, "azit", verbose = verbose) }
  if (identical(azlo, as.name("guess_ab"))) { azlo <- guess_ab(tbl, "azlo", verbose = verbose) }
  if (identical(aztr, as.name("guess_ab"))) { aztr <- guess_ab(tbl, "aztr", verbose = verbose) }
  if (identical(cefa, as.name("guess_ab"))) { cefa <- guess_ab(tbl, "cefa", verbose = verbose) }
  if (identical(cfep, as.name("guess_ab"))) { cfep <- guess_ab(tbl, "cfep", verbose = verbose) }
  if (identical(cfot, as.name("guess_ab"))) { cfot <- guess_ab(tbl, "cfot", verbose = verbose) }
  if (identical(cfox, as.name("guess_ab"))) { cfox <- guess_ab(tbl, "cfox", verbose = verbose) }
  if (identical(cfra, as.name("guess_ab"))) { cfra <- guess_ab(tbl, "cfra", verbose = verbose) }
  if (identical(cfta, as.name("guess_ab"))) { cfta <- guess_ab(tbl, "cfta", verbose = verbose) }
  if (identical(cftr, as.name("guess_ab"))) { cftr <- guess_ab(tbl, "cftr", verbose = verbose) }
  if (identical(cfur, as.name("guess_ab"))) { cfur <- guess_ab(tbl, "cfur", verbose = verbose) }
  if (identical(chlo, as.name("guess_ab"))) { chlo <- guess_ab(tbl, "chlo", verbose = verbose) }
  if (identical(cipr, as.name("guess_ab"))) { cipr <- guess_ab(tbl, "cipr", verbose = verbose) }
  if (identical(clar, as.name("guess_ab"))) { clar <- guess_ab(tbl, "clar", verbose = verbose) }
  if (identical(clin, as.name("guess_ab"))) { clin <- guess_ab(tbl, "clin", verbose = verbose) }
  if (identical(clox, as.name("guess_ab"))) { clox <- guess_ab(tbl, "clox", verbose = verbose) }
  if (identical(coli, as.name("guess_ab"))) { coli <- guess_ab(tbl, "coli", verbose = verbose) }
  if (identical(czol, as.name("guess_ab"))) { czol <- guess_ab(tbl, "czol", verbose = verbose) }
  if (identical(dapt, as.name("guess_ab"))) { dapt <- guess_ab(tbl, "dapt", verbose = verbose) }
  if (identical(doxy, as.name("guess_ab"))) { doxy <- guess_ab(tbl, "doxy", verbose = verbose) }
  if (identical(erta, as.name("guess_ab"))) { erta <- guess_ab(tbl, "erta", verbose = verbose) }
  if (identical(eryt, as.name("guess_ab"))) { eryt <- guess_ab(tbl, "eryt", verbose = verbose) }
  if (identical(fosf, as.name("guess_ab"))) { fosf <- guess_ab(tbl, "fosf", verbose = verbose) }
  if (identical(fusi, as.name("guess_ab"))) { fusi <- guess_ab(tbl, "fusi", verbose = verbose) }
  if (identical(gent, as.name("guess_ab"))) { gent <- guess_ab(tbl, "gent", verbose = verbose) }
  if (identical(imip, as.name("guess_ab"))) { imip <- guess_ab(tbl, "imip", verbose = verbose) }
  if (identical(kana, as.name("guess_ab"))) { kana <- guess_ab(tbl, "kana", verbose = verbose) }
  if (identical(levo, as.name("guess_ab"))) { levo <- guess_ab(tbl, "levo", verbose = verbose) }
  if (identical(linc, as.name("guess_ab"))) { linc <- guess_ab(tbl, "linc", verbose = verbose) }
  if (identical(line, as.name("guess_ab"))) { line <- guess_ab(tbl, "line", verbose = verbose) }
  if (identical(mero, as.name("guess_ab"))) { mero <- guess_ab(tbl, "mero", verbose = verbose) }
  if (identical(mezl, as.name("guess_ab"))) { mezl <- guess_ab(tbl, "mezl", verbose = verbose) }
  if (identical(mino, as.name("guess_ab"))) { mino <- guess_ab(tbl, "mino", verbose = verbose) }
  if (identical(moxi, as.name("guess_ab"))) { moxi <- guess_ab(tbl, "moxi", verbose = verbose) }
  if (identical(nali, as.name("guess_ab"))) { nali <- guess_ab(tbl, "nali", verbose = verbose) }
  if (identical(neom, as.name("guess_ab"))) { neom <- guess_ab(tbl, "neom", verbose = verbose) }
  if (identical(neti, as.name("guess_ab"))) { neti <- guess_ab(tbl, "neti", verbose = verbose) }
  if (identical(nitr, as.name("guess_ab"))) { nitr <- guess_ab(tbl, "nitr", verbose = verbose) }
  if (identical(norf, as.name("guess_ab"))) { norf <- guess_ab(tbl, "norf", verbose = verbose) }
  if (identical(novo, as.name("guess_ab"))) { novo <- guess_ab(tbl, "novo", verbose = verbose) }
  if (identical(oflo, as.name("guess_ab"))) { oflo <- guess_ab(tbl, "oflo", verbose = verbose) }
  if (identical(oxac, as.name("guess_ab"))) { oxac <- guess_ab(tbl, "oxac", verbose = verbose) }
  if (identical(peni, as.name("guess_ab"))) { peni <- guess_ab(tbl, "peni", verbose = verbose) }
  if (identical(pipe, as.name("guess_ab"))) { pipe <- guess_ab(tbl, "pipe", verbose = verbose) }
  if (identical(pita, as.name("guess_ab"))) { pita <- guess_ab(tbl, "pita", verbose = verbose) }
  if (identical(poly, as.name("guess_ab"))) { poly <- guess_ab(tbl, "poly", verbose = verbose) }
  if (identical(pris, as.name("guess_ab"))) { pris <- guess_ab(tbl, "pris", verbose = verbose) }
  if (identical(qida, as.name("guess_ab"))) { qida <- guess_ab(tbl, "qida", verbose = verbose) }
  if (identical(rifa, as.name("guess_ab"))) { rifa <- guess_ab(tbl, "rifa", verbose = verbose) }
  if (identical(roxi, as.name("guess_ab"))) { roxi <- guess_ab(tbl, "roxi", verbose = verbose) }
  if (identical(siso, as.name("guess_ab"))) { siso <- guess_ab(tbl, "siso", verbose = verbose) }
  if (identical(teic, as.name("guess_ab"))) { teic <- guess_ab(tbl, "teic", verbose = verbose) }
  if (identical(tetr, as.name("guess_ab"))) { tetr <- guess_ab(tbl, "tetr", verbose = verbose) }
  if (identical(tica, as.name("guess_ab"))) { tica <- guess_ab(tbl, "tica", verbose = verbose) }
  if (identical(tige, as.name("guess_ab"))) { tige <- guess_ab(tbl, "tige", verbose = verbose) }
  if (identical(tobr, as.name("guess_ab"))) { tobr <- guess_ab(tbl, "tobr", verbose = verbose) }
  if (identical(trim, as.name("guess_ab"))) { trim <- guess_ab(tbl, "trim", verbose = verbose) }
  if (identical(trsu, as.name("guess_ab"))) { trsu <- guess_ab(tbl, "trsu", verbose = verbose) }
  if (identical(vanc, as.name("guess_ab"))) { vanc <- guess_ab(tbl, "vanc", verbose = verbose) }
  col.list <- c(amcl, amik, amox, ampi, azit, azlo, aztr, cefa, cfra, cfep, cfot,
                cfox, cfta, cftr, cfur, chlo, cipr, clar, clin, clox, coli,
                czol, dapt, doxy, erta, eryt, fosf, fusi, gent, imip, kana,
                levo, linc, line, mero, mezl, mino, moxi, nali, neom, neti, nitr,
                novo, norf, oflo, oxac, peni, pipe, pita, poly, pris, qida, rifa,
                roxi, siso, teic, tetr, tica, tige, tobr, trim, trsu, vanc)
  if (length(col.list) < 63) {
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
  pipe <- col.list[pipe]
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

  number_changed <- 0
  number_affected_rows <- integer(0)
  verbose_info <- data.frame(rule_type = character(0),
                             rule_set = character(0),
                             force_to = character(0),
                             found = integer(0),
                             changed = integer(0),
                             target_columns = integer(0),
                             target_rows = integer(0),
                             stringsAsFactors = FALSE)

  # helper function for editing the table
  edit_rsi <- function(to, rule, rows, cols) {
    cols <- unique(cols[!is.na(cols) & !is.null(cols)])
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
      suppressMessages(
        suppressWarnings(
          tbl[rows, cols] <<- to
      ))
      after <- as.character(unlist(as.list(tbl_original[rows, cols])))
      number_changed <<- number_changed + sum(before != after, na.rm = TRUE)
      number_affected_rows <<- unique(c(number_affected_rows, rows))
      changed_results <<- changed_results + sum(before != after, na.rm = TRUE) # will be reset at start of every rule

      if (verbose == TRUE) {
        for (i in 1:length(cols)) {
          # add new row for every affected column
          verbose_new <- data.frame(rule_type = strip_style(rule[1]),
                                    rule_set = strip_style(rule[2]),
                                    force_to = to,
                                    found = length(before),
                                    changed = sum(before != after, na.rm = TRUE),
                                    target_column = cols[i],
                                    stringsAsFactors = FALSE)
          verbose_new$target_rows <- list(unname(rows))
          rownames(verbose_new) <- NULL
          verbose_info <<- rbind(verbose_info, verbose_new)
        }

      }
    }
  }
  na.rm <- function(col) {
    if (is.null(col)) {
      ""
    } else {
      col
    }
  }

  # save original table
  tbl_original <- tbl

  # join to microorganisms data set
  tbl <- tbl %>%
    mutate_at(vars(col_mo), as.mo) %>%
    left_join_microorganisms(by = col_mo, suffix = c("_oldcols", "")) %>%
    as.data.frame(stringsAsFactors = FALSE)

  if (info == TRUE) {
    cat("\nRules by the European Committee on Antimicrobial Susceptibility Testing (EUCAST)\n")
  }

  # since ampicillin ^= amoxicillin, get the first from the latter (not in original EUCAST table)
  if (!is.null(ampi) & !is.null(amox)) {
    if (verbose == TRUE) {
      cat(bgGreen("\n VERBOSE: transforming",
                  length(which(tbl[, amox] == "S" & !tbl[, ampi] %in% c("S", "I", "R"))),
                  "empty ampicillin fields to 'S' based on amoxicillin. "))
      cat(bgGreen("\n VERBOSE: transforming",
                  length(which(tbl[, amox] == "I" & !tbl[, ampi] %in% c("S", "I", "R"))),
                  "empty ampicillin fields to 'I' based on amoxicillin. "))
      cat(bgGreen("\n VERBOSE: transforming",
                  length(which(tbl[, amox] == "R" & !tbl[, ampi] %in% c("S", "I", "R"))),
                  "empty ampicillin fields to 'R' based on amoxicillin. \n"))
    }
    tbl[which(tbl[, amox] == "S" & !tbl[, ampi] %in% c("S", "I", "R")), ampi] <- "S"
    tbl[which(tbl[, amox] == "I" & !tbl[, ampi] %in% c("S", "I", "R")), ampi] <- "I"
    tbl[which(tbl[, amox] == "R" & !tbl[, ampi] %in% c("S", "I", "R")), ampi] <- "R"
  } else if (is.null(ampi) & !is.null(amox)) {
    # ampicillin column is missing, but amoxicillin is available
    message(blue(paste0("NOTE: Using column `", bold(amox), "` as input for ampicillin (J01CA01) since many EUCAST rules depend on it.")))
    ampi <- amox
  }

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
  ureidopenicillins <- c(pipe, pita, azlo, mezl)
  fluoroquinolones <- c(oflo, cipr, norf, levo, moxi)
  all_betalactam <- c(aminopenicillins, ureidopenicillins, cephalosporins, carbapenems, amcl, oxac, clox, peni)

  if (any(c("all", "breakpoints") %in% rules)) {
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

    if (!is.null(ampi)) {
      edit_rsi(to = 'S',
               rule = c(rule_group, rule),
               rows = which(tbl$order == 'Enterobacteriales'
                            & tbl[, ampi] == 'S'),
               cols = amox)
    }
    if (!is.null(ampi)) {
      edit_rsi(to = 'I',
               rule = c(rule_group, rule),
               rows = which(tbl$order == 'Enterobacteriales'
                            & tbl[, ampi] == 'I'),
               cols = amox)
    }
    if (!is.null(ampi)) {
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
    rule <- italic('Staphylococcus')
    if (info == TRUE) {
      warned <- FALSE
      changed_results <- 0
      cat(rule)
    }
    if (!is.null(peni) & !is.null(cfox)) {
      edit_rsi(to = 'S',
               rule = c(rule_group, rule),
               rows = which(tbl$genus == "Staphylococcus"
                            & tbl[, peni] == 'S'
                            & tbl[, cfox] == 'S'),
               cols = c(ampi, amox, pipe, tica))
      edit_rsi(to = 'S',
               rule = c(rule_group, rule),
               rows = which(tbl$genus == "Staphylococcus"
                            & tbl[, peni] == 'R'
                            & tbl[, cfox] == 'S'),
               cols = c(oxac, clox))
    }
    if (!is.null(cfox)) {
      edit_rsi(to = 'R',
               rule = c(rule_group, rule),
               rows = which(tbl$genus == "Staphylococcus"
                            & tbl[, cfox] == 'R'),
               cols = all_betalactam)
    }
    if (!is.null(ampi)) {
      edit_rsi(to = 'S',
               rule = c(rule_group, rule),
               rows = which(tbl$fullname %like% "^Staphylococcus saprophyticus"
                            & tbl[, ampi] == 'S'),
               cols = c(amox, amcl, pipe, pita))
    }
    if (!is.null(cfox)) {
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
    if (!is.null(norf)) {
      edit_rsi(to = 'S',
               rule = c(rule_group, rule),
               rows = which(tbl$genus == "Staphylococcus"
                            & tbl[, norf] == 'S'),
               cols = c(cipr, levo, moxi, oflo))
    }
    if (!is.null(eryt)) {
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
    if (!is.null(tetr)) {
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
    rule <- italic('Enterococcus')
    if (info == TRUE) {
      warned <- FALSE
      changed_results <- 0
      cat(rule)
    }
    if (!is.null(ampi)) { # penicillin group
      edit_rsi(to = 'R',
               rule = c(rule_group, rule),
               rows = which(tbl$fullname %like% "^Enterococcus faecium"
                            & tbl[, ampi] == 'R'),
               cols = all_betalactam)
    }
    if (!is.null(ampi)) {
      edit_rsi(to = 'S',
               rule = c(rule_group, rule),
               rows = which(tbl$genus == "Enterococcus"
                            & tbl[, ampi] == 'S'),
               cols = c(amox, amcl, pipe, pita))
      edit_rsi(to = 'I',
               rule = c(rule_group, rule),
               rows = which(tbl$genus == "Enterococcus"
                            & tbl[, ampi] == 'I'),
               cols = c(amox, amcl, pipe, pita))
      edit_rsi(to = 'R',
               rule = c(rule_group, rule),
               rows = which(tbl$genus == "Enterococcus"
                            & tbl[, ampi] == 'R'),
               cols = c(amox, amcl, pipe, pita))
    }
    if (!is.null(norf)) {
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
    # Streptococcus groups A, B, C, G----
    rule <- paste(italic('Streptococcus'), 'groups A, B, C, G')
    if (info == TRUE) {
      warned <- FALSE
      changed_results <- 0
      cat(rule)
    }
    if (!is.null(peni)) {
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
    if (!is.null(norf)) {
      edit_rsi(to = 'S',
               rule = c(rule_group, rule),
               rows = which(tbl$fullname %like% "^Streptococcus (pyogenes|agalactiae|dysgalactiae|group A|group B|group C|group G)"
                            & tbl[, norf] == 'S'),
               cols = c(levo, moxi))
    }
    if (!is.null(eryt)) {
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
    if (!is.null(tetr)) {
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
    rule <- italic('Streptococcus pneumoniae')
    if (info == TRUE) {
      warned <- FALSE
      changed_results <- 0
      cat(rule)
    }
    if (!is.null(peni)) {
      edit_rsi(to = 'S',
               rule = c(rule_group, rule),
               rows = which(tbl$fullname %like% "^Streptococcus pneumoniae"
                            & tbl[, peni] == 'S'),
               cols = c(ampi, amox, amcl, pipe, pita))
    }
    if (!is.null(ampi)) {
      edit_rsi(to = 'S',
               rule = c(rule_group, rule),
               rows = which(tbl$fullname %like% "^Streptococcus pneumoniae"
                            & tbl[, ampi] == 'S'),
               cols = c(amox, amcl, pipe, pita))
      edit_rsi(to = 'I',
               rule = c(rule_group, rule),
               rows = which(tbl$fullname %like% "^Streptococcus pneumoniae"
                            & tbl[, ampi] == 'I'),
               cols = c(amox, amcl, pipe, pita))
      edit_rsi(to = 'R',
               rule = c(rule_group, rule),
               rows = which(tbl$fullname %like% "^Streptococcus pneumoniae"
                            & tbl[, ampi] == 'R'),
               cols = c(amox, amcl, pipe, pita))
    }
    if (!is.null(norf)) {
      edit_rsi(to = 'S',
               rule = c(rule_group, rule),
               rows = which(tbl$fullname %like% "^Streptococcus pneumoniae"
                            & tbl[, norf] == 'S'),
               cols = c(levo, moxi))
    }
    if (!is.null(eryt)) {
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
    if (!is.null(tetr)) {
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
    if (!is.null(peni)) {
      edit_rsi(to = 'S',
               rule = c(rule_group, rule),
               rows = which(tbl$genus == "Streptococcus" & tbl$species %in% viridans_group
                            & tbl[, peni] == 'S'),
               cols = c(ampi, amox, amcl, pipe, pita))
    }
    if (!is.null(ampi)) {
      edit_rsi(to = 'S',
               rule = c(rule_group, rule),
               rows = which(tbl$genus == "Streptococcus" & tbl$species %in% viridans_group
                            & tbl[, ampi] == 'S'),
               cols = c(amox, amcl, pipe, pita))
      edit_rsi(to = 'I',
               rule = c(rule_group, rule),
               rows = which(tbl$genus == "Streptococcus" & tbl$species %in% viridans_group
                            & tbl[, ampi] == 'I'),
               cols = c(amox, amcl, pipe, pita))
      edit_rsi(to = 'R',
               rule = c(rule_group, rule),
               rows = which(tbl$genus == "Streptococcus" & tbl$species %in% viridans_group
                            & tbl[, ampi] == 'R'),
               cols = c(amox, amcl, pipe, pita))
    }
    if (info == TRUE) {
      txt_ok()
    }
    # Haemophilus influenzae ----
    rule <- italic('Haemophilus influenzae')
    if (info == TRUE) {
      warned <- FALSE
      changed_results <- 0
      cat(rule)
    }
    if (!is.null(ampi)) {
      edit_rsi(to = 'S',
               rule = c(rule_group, rule),
               rows = which(tbl$fullname %like% "^Haemophilus influenzae"
                            & tbl[, ampi] == 'S'),
               cols = c(amox, pipe))
      edit_rsi(to = 'I',
               rule = c(rule_group, rule),
               rows = which(tbl$fullname %like% "^Haemophilus influenzae"
                            & tbl[, ampi] == 'I'),
               cols = c(amox, pipe))
      edit_rsi(to = 'R',
               rule = c(rule_group, rule),
               rows = which(tbl$fullname %like% "^Haemophilus influenzae"
                            & tbl[, ampi] == 'R'),
               cols = c(amox, pipe))
    }
    if (!is.null(peni)) {
      edit_rsi(to = 'S',
               rule = c(rule_group, rule),
               rows = which(tbl$fullname %like% "^Haemophilus influenzae"
                            & tbl[, peni] == 'S'),
               cols = c(ampi, amox, amcl, pipe, pita))
    }
    if (!is.null(amcl)) {
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
    if (!is.null(nali)) {
      edit_rsi(to = 'S',
               rule = c(rule_group, rule),
               rows = which(tbl$fullname %like% "^Haemophilus influenzae"
                            & tbl[, nali] == 'S'),
               cols = c(cipr, levo, moxi, oflo))
    }
    if (!is.null(tetr)) {
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
    rule <- italic('Moraxella catarrhalis')
    if (info == TRUE) {
      warned <- FALSE
      changed_results <- 0
      cat(rule)
    }
    if (!is.null(amcl)) {
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
    if (!is.null(nali)) {
      edit_rsi(to = 'S',
               rule = c(rule_group, rule),
               rows = which(tbl$fullname %like% "^Moraxella catarrhalis"
                            & tbl[, nali] == 'S'),
               cols = c(cipr, levo, moxi, oflo))
    }
    if (!is.null(eryt)) {
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
    if (!is.null(tetr)) {
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
    if (!is.null(peni)) {
      edit_rsi(to = 'S',
               rule = c(rule_group, rule),
               rows = which(tbl$genus %in% c("Clostridium", "Actinomyces", "Propionibacterium",
                                             "Cutibacterium", # new name of Propionibacterium
                                             "Bifidobacterium", "Eggerthella", "Eubacterium",
                                             "Lactobacillus ", "Actinomyces")
                            & tbl[, peni] == 'S'),
               cols = c(ampi, amox, pipe, pita, tica))
      edit_rsi(to = 'I',
               rule = c(rule_group, rule),
               rows = which(tbl$genus %in% c("Clostridium", "Actinomyces", "Propionibacterium",
                                             "Cutibacterium", # new name of Propionibacterium
                                             "Bifidobacterium", "Eggerthella", "Eubacterium",
                                             "Lactobacillus ", "Actinomyces")
                            & tbl[, peni] == 'I'),
               cols = c(ampi, amox, pipe, pita, tica))
      edit_rsi(to = 'R',
               rule = c(rule_group, rule),
               rows = which(tbl$genus %in% c("Clostridium", "Actinomyces", "Propionibacterium",
                                             "Cutibacterium", # new name of Propionibacterium
                                             "Bifidobacterium", "Eggerthella", "Eubacterium",
                                             "Lactobacillus ", "Actinomyces")
                            & tbl[, peni] == 'R'),
               cols = c(ampi, amox, pipe, pita, tica))
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
    if (!is.null(peni)) {
      edit_rsi(to = 'S',
               rule = c(rule_group, rule),
               rows = which(tbl$genus %in% c("Bacteroides", "Prevotella", "Porphyromonas",
                                             "Fusobacterium", "Bilophila ", "Mobiluncus")
                            & tbl[, peni] == 'S'),
               cols = c(ampi, amox, pipe, pita, tica))
      edit_rsi(to = 'I',
               rule = c(rule_group, rule),
               rows = which(tbl$genus %in% c("Bacteroides", "Prevotella", "Porphyromonas",
                                             "Fusobacterium", "Bilophila ", "Mobiluncus")
                            & tbl[, peni] == 'I'),
               cols = c(ampi, amox, pipe, pita, tica))
      edit_rsi(to = 'R',
               rule = c(rule_group, rule),
               rows = which(tbl$genus %in% c("Bacteroides", "Prevotella", "Porphyromonas",
                                             "Fusobacterium", "Bilophila ", "Mobiluncus")
                            & tbl[, peni] == 'R'),
               cols = c(ampi, amox, pipe, pita, tica))
    }
    if (info == TRUE) {
      txt_ok()
    }
    # Pasteurella multocida ----
    rule <- italic('Pasteurella multocida')
    if (info == TRUE) {
      warned <- FALSE
      changed_results <- 0
      cat(rule)
    }
    if (!is.null(peni)) {
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
    # Campylobacter jejuni and coli ----
    rule <- paste(italic('Campylobacter jejuni'), 'and', italic('C. coli'))
    if (info == TRUE) {
      warned <- FALSE
      changed_results <- 0
      cat(rule)
    }
    if (!is.null(eryt)) {
      edit_rsi(to = 'S',
               rule = c(rule_group, rule),
               rows = which(tbl$fullname %like% "^Campylobacter (jejuni|coli)"
                            & tbl[, eryt] == 'S'),
               cols = c(azit, clar))
      edit_rsi(to = 'I',
               rule = c(rule_group, rule),
               rows = which(tbl$fullname %like% "^Campylobacter (jejuni|coli)"
                            & tbl[, eryt] == 'I'),
               cols = c(azit, clar))
      edit_rsi(to = 'R',
               rule = c(rule_group, rule),
               rows = which(tbl$fullname %like% "^Campylobacter (jejuni|coli)"
                            & tbl[, eryt] == 'R'),
               cols = c(azit, clar))
    }
    if (!is.null(tetr)) {
      edit_rsi(to = 'S',
               rule = c(rule_group, rule),
               rows = which(tbl$fullname %like% "^Campylobacter (jejuni|coli)"
                            & tbl[, tetr] == 'S'),
               cols = doxy)
      edit_rsi(to = 'I',
               rule = c(rule_group, rule),
               rows = which(tbl$fullname %like% "^Campylobacter (jejuni|coli)"
                            & tbl[, tetr] == 'I'),
               cols = doxy)
      edit_rsi(to = 'R',
               rule = c(rule_group, rule),
               rows = which(tbl$fullname %like% "^Campylobacter (jejuni|coli)"
                            & tbl[, tetr] == 'R'),
               cols = doxy)
    }
    if (info == TRUE) {
      txt_ok()
    }
    # Aerococcus sanguinicola/urinae ----
    rule <- paste(italic('Aerococcus sanguinicola'), 'and', italic('A. urinae'))
    if (info == TRUE) {
      warned <- FALSE
      changed_results <- 0
      cat(rule)
    }
    if (!is.null(norf)) {
      edit_rsi(to = 'S',
               rule = c(rule_group, rule),
               rows = which(tbl$fullname %like% "^Aerococcus (sanguinicola|urinae)"
                            & tbl[, norf] == 'S'),
               cols = fluoroquinolones)
      edit_rsi(to = 'I',
               rule = c(rule_group, rule),
               rows = which(tbl$fullname %like% "^Aerococcus (sanguinicola|urinae)"
                            & tbl[, norf] == 'I'),
               cols = fluoroquinolones)
      edit_rsi(to = 'R',
               rule = c(rule_group, rule),
               rows = which(tbl$fullname %like% "^Aerococcus (sanguinicola|urinae)"
                            & tbl[, norf] == 'R'),
               cols = fluoroquinolones)
    }
    if (!is.null(cipr)) {
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
    rule <- italic('Kingella kingae')
    if (info == TRUE) {
      warned <- FALSE
      changed_results <- 0
      cat(rule)
    }
    if (!is.null(peni)) {
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
    if (!is.null(eryt)) {
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
    if (!is.null(tetr)) {
      edit_rsi(to = 'S',
               rule = c(rule_group, rule),
               rows = which(tbl$fullname %like% "^Kingella kingae"
                            & tbl[, tetr] == 'S'),
               cols = doxy)
    }
    if (info == TRUE) {
      txt_ok()
    }

  } # end of breakpoints
  if (any(c("all", "expert") %in% rules)) {

    # EXPERT RULES AND INTRINSIC RESISTANCE -----------------------------------

    if (info == TRUE) {
      cat(bold(paste0('\nEUCAST Expert Rules, Intrinsic Resistance and Exceptional Phenotypes (v', EUCAST_VERSION_EXPERT_RULES, ')\n')))
    }
    rule_group <- "Expert Rules"

    # Table 1: Intrinsic resistance in Enterobacteriaceae ----
    rule <- paste('Table 1:  Intrinsic resistance in', italic('Enterobacteriaceae'))
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
             # the 'Burkholderia cepacia complex' are all these species: (PMID 16217180)
             rows = which(tbl$fullname %like% '^Burkholderia (cepacia|multivorans|cenocepacia|stabilis|vietnamiensis|dolosa|ambifaria|anthina|pyrrocinia|ubonensis)'),
             cols = c(aminopenicillins, amcl, tica, pipe, pita, czol, cfot, cftr, aztr, erta, cipr, chlo, aminoglycosides, trim, fosf, polymyxins))
    # Elizabethkingia
    edit_rsi(to = 'R',
             rule = c(rule_group, rule),
             rows = which(tbl$fullname %like% '^Elizabethkingia meningoseptic(a|um)'),
             cols = c(aminopenicillins, amcl, tica, czol, cfot, cftr, cfta, cfep, aztr, erta, imip, mero, polymyxins))
    # Ochrobactrum
    edit_rsi(to = 'R',
             rule = c(rule_group, rule),
             rows = which(tbl$fullname %like% '^Ochrobactrum anthropi'),
             cols = c(aminopenicillins, amcl, tica, pipe, pita, czol, cfot, cftr, cfta, cfep, aztr, erta))
    # Pseudomonas
    edit_rsi(to = 'R',
             rule = c(rule_group, rule),
             rows = which(tbl$fullname %like% '^Pseudomonas aeruginosa'),
             cols = c(aminopenicillins, amcl, czol, cfot, cftr, erta, chlo, kana, neom, trim, trsu, tetracyclines, tige))
    # Stenotrophomonas
    edit_rsi(to = 'R',
             rule = c(rule_group, rule),
             rows = which(tbl$fullname %like% '^Stenotrophomonas maltophilia'),
             cols = c(aminopenicillins, amcl, tica, pipe, pita, czol, cfot, cftr, cfta, aztr, erta, imip, mero, aminoglycosides, trim, fosf, tetr))
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
    if (!is.null(peni)) {
      edit_rsi(to = 'S',
               rule = c(rule_group, rule),
               rows = which(tbl$fullname %like% '^Streptococcus (pyogenes|agalactiae|dysgalactiae|group A|group B|group C|group G)'
                            & tbl[, peni] == 'S'),
               cols = c(aminopenicillins, cephalosporins, carbapenems))
    }
    # rule 8.6
    if (!is.null(ampi)) {
      edit_rsi(to = 'R',
               rule = c(rule_group, rule),
               rows = which(tbl$genus == 'Enterococcus'
                            & tbl[, ampi] == 'R'),
               cols = c(ureidopenicillins, carbapenems))
    }
    if (!is.null(amox)) {
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
    if (!is.null(tica) & !is.null(pipe)) {
      edit_rsi(to = 'R',
               rule = c(rule_group, rule),
               rows = which(tbl$family == 'Enterobacteriaceae'
                            & tbl[, tica] == 'R'
                            & tbl[, pipe] == 'S'),
               cols = pipe)
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
    # if (!is.null(ampi)) {
    # you should know first if the are B-lactamase positive, so do not run for now
    # edit_rsi(to = 'R',
    #          rule = c(rule_group, rule),
    #          rows = which(tbl$fullname %like% '^Haemophilus influenza'
    #                       & tbl[, ampi] == 'R'),
    #          cols = c(ampi, amox, amcl, pipe, pita, cfur))
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
    if (!is.null(eryt)) {
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
    if (!is.null(tobr)) {
      edit_rsi(to = 'R',
               rule = c(rule_group, rule),
               rows = which(tbl$genus == 'Staphylococcus'
                            & tbl[, tobr] == 'R'),
               cols = c(kana, amik))
    }
    # rule 12.3
    if (!is.null(gent)) {
      edit_rsi(to = 'R',
               rule = c(rule_group, rule),
               rows = which(tbl$genus == 'Staphylococcus'
                            & tbl[, gent] == 'R'),
               cols = aminoglycosides)
    }
    # rule 12.8
    if (!is.null(gent) & !is.null(tobr)) {
      edit_rsi(to = 'R',
               rule = c(rule_group, rule),
               rows = which(tbl$family == 'Enterobacteriaceae'
                            & tbl[, gent] == 'I'
                            & tbl[, tobr] == 'S'),
               cols = gent)
    }
    # rule 12.9
    if (!is.null(gent) & !is.null(tobr)) {
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
    if (!is.null(moxi)) {
      edit_rsi(to = 'R',
               rule = c(rule_group, rule),
               rows = which(tbl$genus == 'Staphylococcus'
                            & tbl[, moxi] == 'R'),
               cols = fluoroquinolones)
    }
    # rule 13.4
    if (!is.null(moxi)) {
      edit_rsi(to = 'R',
               rule = c(rule_group, rule),
               rows = which(tbl$fullname %like% '^Streptococcus pneumoniae'
                            & tbl[, moxi] == 'R'),
               cols = fluoroquinolones)
    }
    # rule 13.5
    if (!is.null(cipr)) {
      edit_rsi(to = 'R',
               rule = c(rule_group, rule),
               rows = which(tbl$family == 'Enterobacteriaceae'
                            & tbl[, cipr] == 'R'),
               cols = fluoroquinolones)
    }
    # rule 13.8
    if (!is.null(cipr)) {
      edit_rsi(to = 'R',
               rule = c(rule_group, rule),
               rows = which(tbl$fullname %like% '^Neisseria gonorrhoeae'
                            & tbl[, cipr] == 'R'),
               cols = fluoroquinolones)
    }
    if (info == TRUE) {
      txt_ok()
    }

  } # end of expert rules
  if (any(c("all", "other") %in% rules)) {

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
    if (!is.null(amcl)) {
      edit_rsi(to = 'R',
               rule = c(rule_group, rule),
               rows = which(tbl[, amcl] == 'R'),
               cols = ampi)
    }
    if (info == TRUE) {
      txt_ok()
    }
    rule <- 'Non-EUCAST: piperacillin = R where piperacillin/tazobactam = R'
    if (info == TRUE) {
      warned <- FALSE
      changed_results <- 0
      cat(rule)
    }
    if (!is.null(pita)) {
      edit_rsi(to = 'R',
               rule = c(rule_group, rule),
               rows = which(tbl[, pita] == 'R'),
               cols = pipe)
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
    if (!is.null(trsu)) {
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
    if (!is.null(ampi)) {
      edit_rsi(to = 'S',
               rule = c(rule_group, rule),
               rows = which(tbl[, ampi] == 'S'),
               cols = amcl)
    }
    if (info == TRUE) {
      txt_ok()
    }
    rule <- 'Non-EUCAST: piperacillin/tazobactam = S where piperacillin = S'
    if (info == TRUE) {
      warned <- FALSE
      changed_results <- 0
      cat(rule)
    }
    if (!is.null(pipe)) {
      edit_rsi(to = 'S',
               rule = c(rule_group, rule),
               rows = which(tbl[, pipe] == 'S'),
               cols = pita)
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
    if (!is.null(trim)) {
      edit_rsi(to = 'S',
               rule = c(rule_group, rule),
               rows = which(tbl[, trim] == 'S'),
               cols = trsu)
    }
    if (info == TRUE) {
      txt_ok()
    }

  } # end of other rules

  # restore old col_mo values if needed
  # if (!is.null(col_mo_original)) {
  #   tbl_original[, col_mo] <- col_mo_original
  # }

  if (info == TRUE) {
    if (verbose == TRUE) {
      wouldve <- "would have "
    } else {
      wouldve <- ""
    }
    if (number_changed == 0) {
      colour <- green
    } else {
      colour <- blue
    }
    decimal.mark <- getOption("OutDec")
    big.mark <- ifelse(decimal.mark != ",", ",", ".")
    cat(bold(paste('\n=> EUCAST rules', paste0(wouldve, 'affected'),
             number_affected_rows %>% length() %>% format(big.mark = big.mark, decimal.mark = decimal.mark),
             'out of', nrow(tbl_original) %>% format(big.mark = big.mark, decimal.mark = decimal.mark),
             'rows ->',
             colour(paste0(wouldve, 'changed'),
                    number_changed %>% format(big.mark = big.mark, decimal.mark = decimal.mark), 'test results.\n\n'))))
  }

  if (verbose == TRUE) {
    return(verbose_info)
  }

  tbl_original
}

#' @rdname eucast_rules
#' @export
EUCAST_rules <- function(...) {
  .Deprecated("eucast_rules")
  eucast_rules(...)
}

#' @rdname eucast_rules
#' @export
interpretive_reading <- function(...) {
  .Deprecated("eucast_rules")
  eucast_rules(...)
}
