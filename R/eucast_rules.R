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

# global variables
EUCAST_RULES_FILE_LOCATION <- system.file("eucast/eucast_rules.tsv", package = "AMR")
EUCAST_VERSION_BREAKPOINTS <- "9.0, 2019"
EUCAST_VERSION_EXPERT_RULES <- "3.1, 2016"

#' EUCAST rules
#'
#' Apply susceptibility rules as defined by the European Committee on Antimicrobial Susceptibility Testing (EUCAST, \url{http://eucast.org}), see \emph{Source}. This includes (1) expert rules, (2) intrinsic resistance and (3) inferred resistance as defined in their breakpoint tables.
#' @param x data with antibiotic columns, like e.g. \code{amox} and \code{amcl}
#' @param info print progress
#' @param rules a character vector that specifies which rules should be applied - one or more of \code{c("breakpoints", "expert", "other", "all")}
#' @param verbose a logical to indicate whether extensive info should be returned as a \code{data.frame} with info about which rows and columns are effected. It runs all EUCAST rules, but will not be applied to an output - only an informative \code{data.frame} with changes will be returned as output.
#' @param amcl,amik,amox,ampi,azit,azlo,aztr,cefa,cfep,cfot,cfox,cfra,cfta,cftr,cfur,chlo,cipr,clar,clin,clox,coli,czol,dapt,doxy,erta,eryt,fosf,fusi,gent,imip,kana,levo,linc,line,mero,mezl,mino,moxi,nali,neom,neti,nitr,norf,novo,oflo,oxac,peni,pipe,pita,poly,pris,qida,rifa,roxi,siso,teic,tetr,tica,tige,tobr,trim,trsu,vanc column name of an antibiotic, see Antibiotics
#' @param ... parameters that are passed on to \code{eucast_rules}
#' @inheritParams first_isolate
#' @details
#' \strong{NOTE:} This function does not translate MIC values to RSI values. It only applies (1) inferred susceptibility and resistance based on results of other antibiotics and (2) intrinsic resistance based on taxonomic properties of a microorganism.
#'
#' The file used for applying all EUCAST rules can be retrieved with \code{\link{eucast_rules_file}()}. It returns an easily readable data set containing all rules. The original TSV file (tab separated file) that is being read by this function can be found when running this command: \cr
#' \code{AMR::EUCAST_RULES_FILE_LOCATION} (without brackets).
#'
#' In the source code it is located under \href{https://gitlab.com/msberends/AMR/blob/master/inst/eucast/eucast_rules.tsv}{\code{./inst/eucast/eucast_rules.tsv}}.
#'
#' \strong{Note:} When ampicillin (J01CA01) is not available but amoxicillin (J01CA04) is, the latter will be used for all rules where there is a dependency on ampicillin. These drugs are interchangeable when it comes to expression of antimicrobial resistance.
#' @section Antibiotics:
#' To define antibiotics column names, leave as it is to determine it automatically with \code{\link{guess_ab_col}} or input a text (case-insensitive) or use \code{NULL} to skip a column (e.g. \code{tica = NULL}). Non-existing columns will anyway be skipped with a warning.
#'
#' Abbrevations of the column containing antibiotics in the form: \strong{abbreviation}: generic name (\emph{ATC code})
#'
#'  \strong{amcl}: amoxicillin+clavulanic acid (\href{https://www.whocc.no/atc_ddd_index/?code=J01CR02}{J01CR02}),
#'  \strong{amik}: amikacin (\href{https://www.whocc.no/atc_ddd_index/?code=J01GB06}{J01GB06}),
#'  \strong{amox}: amoxicillin (\href{https://www.whocc.no/atc_ddd_index/?code=J01CA04}{J01CA04}),
#'  \strong{ampi}: ampicillin (\href{https://www.whocc.no/atc_ddd_index/?code=J01CA01}{J01CA01}),
#'  \strong{azit}: azithromycin (\href{https://www.whocc.no/atc_ddd_index/?code=J01FA10}{J01FA10}),
#'  \strong{azlo}: azlocillin (\href{https://www.whocc.no/atc_ddd_index/?code=J01CA09}{J01CA09}),
#'  \strong{aztr}: aztreonam (\href{https://www.whocc.no/atc_ddd_index/?code=J01DF01}{J01DF01}),
#'  \strong{cefa}: cefaloridine (\href{https://www.whocc.no/atc_ddd_index/?code=J01DB02}{J01DB02}),
#'  \strong{cfep}: cefepime (\href{https://www.whocc.no/atc_ddd_index/?code=J01DE01}{J01DE01}),
#'  \strong{cfot}: cefotaxime (\href{https://www.whocc.no/atc_ddd_index/?code=J01DD01}{J01DD01}),
#'  \strong{cfox}: cefoxitin (\href{https://www.whocc.no/atc_ddd_index/?code=J01DC01}{J01DC01}),
#'  \strong{cfra}: cefradine (\href{https://www.whocc.no/atc_ddd_index/?code=J01DB09}{J01DB09}),
#'  \strong{cfta}: ceftazidime (\href{https://www.whocc.no/atc_ddd_index/?code=J01DD02}{J01DD02}),
#'  \strong{cftr}: ceftriaxone (\href{https://www.whocc.no/atc_ddd_index/?code=J01DD04}{J01DD04}),
#'  \strong{cfur}: cefuroxime (\href{https://www.whocc.no/atc_ddd_index/?code=J01DC02}{J01DC02}),
#'  \strong{chlo}: chloramphenicol (\href{https://www.whocc.no/atc_ddd_index/?code=J01BA01}{J01BA01}),
#'  \strong{cipr}: ciprofloxacin (\href{https://www.whocc.no/atc_ddd_index/?code=J01MA02}{J01MA02}),
#'  \strong{clar}: clarithromycin (\href{https://www.whocc.no/atc_ddd_index/?code=J01FA09}{J01FA09}),
#'  \strong{clin}: clindamycin (\href{https://www.whocc.no/atc_ddd_index/?code=J01FF01}{J01FF01}),
#'  \strong{clox}: flucloxacillin (\href{https://www.whocc.no/atc_ddd_index/?code=J01CF05}{J01CF05}),
#'  \strong{coli}: colistin (\href{https://www.whocc.no/atc_ddd_index/?code=J01XB01}{J01XB01}),
#'  \strong{czol}: cefazolin (\href{https://www.whocc.no/atc_ddd_index/?code=J01DB04}{J01DB04}),
#'  \strong{dapt}: daptomycin (\href{https://www.whocc.no/atc_ddd_index/?code=J01XX09}{J01XX09}),
#'  \strong{doxy}: doxycycline (\href{https://www.whocc.no/atc_ddd_index/?code=J01AA02}{J01AA02}),
#'  \strong{erta}: ertapenem (\href{https://www.whocc.no/atc_ddd_index/?code=J01DH03}{J01DH03}),
#'  \strong{eryt}: erythromycin (\href{https://www.whocc.no/atc_ddd_index/?code=J01FA01}{J01FA01}),
#'  \strong{fosf}: fosfomycin (\href{https://www.whocc.no/atc_ddd_index/?code=J01XX01}{J01XX01}),
#'  \strong{fusi}: fusidic acid (\href{https://www.whocc.no/atc_ddd_index/?code=J01XC01}{J01XC01}),
#'  \strong{gent}: gentamicin (\href{https://www.whocc.no/atc_ddd_index/?code=J01GB03}{J01GB03}),
#'  \strong{imip}: imipenem (\href{https://www.whocc.no/atc_ddd_index/?code=J01DH51}{J01DH51}),
#'  \strong{kana}: kanamycin (\href{https://www.whocc.no/atc_ddd_index/?code=J01GB04}{J01GB04}),
#'  \strong{levo}: levofloxacin (\href{https://www.whocc.no/atc_ddd_index/?code=J01MA12}{J01MA12}),
#'  \strong{linc}: lincomycin (\href{https://www.whocc.no/atc_ddd_index/?code=J01FF02}{J01FF02}),
#'  \strong{line}: linezolid (\href{https://www.whocc.no/atc_ddd_index/?code=J01XX08}{J01XX08}),
#'  \strong{mero}: meropenem (\href{https://www.whocc.no/atc_ddd_index/?code=J01DH02}{J01DH02}),
#'  \strong{mezl}: mezlocillin (\href{https://www.whocc.no/atc_ddd_index/?code=J01CA10}{J01CA10}),
#'  \strong{mino}: minocycline (\href{https://www.whocc.no/atc_ddd_index/?code=J01AA08}{J01AA08}),
#'  \strong{moxi}: moxifloxacin (\href{https://www.whocc.no/atc_ddd_index/?code=J01MA14}{J01MA14}),
#'  \strong{nali}: nalidixic acid (\href{https://www.whocc.no/atc_ddd_index/?code=J01MB02}{J01MB02}),
#'  \strong{neom}: neomycin (\href{https://www.whocc.no/atc_ddd_index/?code=J01GB05}{J01GB05}),
#'  \strong{neti}: netilmicin (\href{https://www.whocc.no/atc_ddd_index/?code=J01GB07}{J01GB07}),
#'  \strong{nitr}: nitrofurantoin (\href{https://www.whocc.no/atc_ddd_index/?code=J01XE01}{J01XE01}),
#'  \strong{norf}: norfloxacin (\href{https://www.whocc.no/atc_ddd_index/?code=J01MA06}{J01MA06}),
#'  \strong{novo}: novobiocin (an ATCvet code: \href{https://www.whocc.no/atc_ddd_index/?code=QJ01XX95}{QJ01XX95}),
#'  \strong{oflo}: ofloxacin (\href{https://www.whocc.no/atc_ddd_index/?code=J01MA01}{J01MA01}),
#'  \strong{peni}: (benzyl)penicillin (\href{https://www.whocc.no/atc_ddd_index/?code=J01CE01}{J01CE01}),
#'  \strong{pipe}: piperacillin (\href{https://www.whocc.no/atc_ddd_index/?code=J01CA12}{J01CA12}),
#'  \strong{pita}: piperacillin+tazobactam (\href{https://www.whocc.no/atc_ddd_index/?code=J01CR05}{J01CR05}),
#'  \strong{poly}: polymyxin B (\href{https://www.whocc.no/atc_ddd_index/?code=J01XB02}{J01XB02}),
#'  \strong{pris}: pristinamycin (\href{https://www.whocc.no/atc_ddd_index/?code=J01FG01}{J01FG01}),
#'  \strong{qida}: quinupristin/dalfopristin (\href{https://www.whocc.no/atc_ddd_index/?code=J01FG02}{J01FG02}),
#'  \strong{rifa}: rifampicin (\href{https://www.whocc.no/atc_ddd_index/?code=J04AB02}{J04AB02}),
#'  \strong{roxi}: roxithromycin (\href{https://www.whocc.no/atc_ddd_index/?code=J01FA06}{J01FA06}),
#'  \strong{siso}: sisomicin (\href{https://www.whocc.no/atc_ddd_index/?code=J01GB08}{J01GB08}),
#'  \strong{teic}: teicoplanin (\href{https://www.whocc.no/atc_ddd_index/?code=J01XA02}{J01XA02}),
#'  \strong{tetr}: tetracycline (\href{https://www.whocc.no/atc_ddd_index/?code=J01AA07}{J01AA07}),
#'  \strong{tica}: ticarcillin (\href{https://www.whocc.no/atc_ddd_index/?code=J01CA13}{J01CA13}),
#'  \strong{tige}: tigecycline (\href{https://www.whocc.no/atc_ddd_index/?code=J01AA12}{J01AA12}),
#'  \strong{tobr}: tobramycin (\href{https://www.whocc.no/atc_ddd_index/?code=J01GB01}{J01GB01}),
#'  \strong{trim}: trimethoprim (\href{https://www.whocc.no/atc_ddd_index/?code=J01EA01}{J01EA01}),
#'  \strong{trsu}: sulfamethoxazole and trimethoprim (\href{https://www.whocc.no/atc_ddd_index/?code=J01EE01}{J01EE01}),
#'  \strong{vanc}: vancomycin (\href{https://www.whocc.no/atc_ddd_index/?code=J01XA01}{J01XA01}).
#' @keywords interpretive eucast reading resistance
#' @rdname eucast_rules
#' @export
#' @importFrom dplyr %>% select pull mutate_at vars group_by summarise n
#' @importFrom crayon bold bgGreen bgYellow bgRed black green blue italic strip_style
#' @return The input of \code{tbl_}, possibly with edited values of antibiotics. Or, if \code{verbose = TRUE}, a \code{data.frame} with all original and new values of the affected bug-drug combinations.
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
#'       EUCAST Breakpoint tables for interpretation of MICs and zone diameters. Version 9.0, 2019. \cr
#'       \url{http://www.eucast.org/fileadmin/src/media/PDFs/EUCAST_files/Breakpoint_tables/v_9.0_Breakpoint_Tables.xlsx}
#'     }
#'   }
#'
#'   For editing the reference file (which is available with \code{\link{eucast_rules_file}}), these values can all be used for target antibiotics: aminoglycosides, tetracyclines, polymyxins, macrolides, glycopeptides, streptogramins, cephalosporins, cephalosporins_without_cfta, carbapenems, aminopenicillins, ureidopenicillins, fluoroquinolones, all_betalactams, and all separate four letter codes like amcl. They can be separated by comma: \code{"amcl, fluoroquinolones"}. The mo_property can be any column name from the \code{\link{microorganisms}} data set, or \code{genus_species} or \code{gramstain}. This file contains references to the 'Burkholderia cepacia complex'. The species in this group can be found in: LiPuma JJ, 2015 (PMID 16217180).
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
#'
#' # apply EUCAST rules: 18 results are forced as R or S
#' b <- eucast_rules(a)
#'
#' b
#' #                       mo vanc amox coli cfta cfur peni cfox
#' # 1  Staphylococcus aureus    -    S    R    R    S    S    S
#' # 2  Enterococcus faecalis    -    -    R    R    R    S    R
#' # 3       Escherichia coli    R    -    -    -    -    R    S
#' # 4  Klebsiella pneumoniae    R    R    -    -    -    R    S
#' # 5 Pseudomonas aeruginosa    R    R    -    -    R    R    R
#'
#'
#' # do not apply EUCAST rules, but rather get a a data.frame
#' # with 18 rows, containing all details about the transformations:
#' c <- eucast_rules(a, verbose = TRUE)
eucast_rules <- function(x,
                         col_mo = NULL,
                         info = TRUE,
                         rules = c("breakpoints", "expert", "other", "all"),
                         verbose = FALSE,
                         amcl = guess_ab_col(),
                         amik = guess_ab_col(),
                         amox = guess_ab_col(),
                         ampi = guess_ab_col(),
                         azit = guess_ab_col(),
                         azlo = guess_ab_col(),
                         aztr = guess_ab_col(),
                         cefa = guess_ab_col(),
                         cfep = guess_ab_col(),
                         cfot = guess_ab_col(),
                         cfox = guess_ab_col(),
                         cfra = guess_ab_col(),
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
                         mezl = guess_ab_col(),
                         mino = guess_ab_col(),
                         moxi = guess_ab_col(),
                         nali = guess_ab_col(),
                         neom = guess_ab_col(),
                         neti = guess_ab_col(),
                         nitr = guess_ab_col(),
                         norf = guess_ab_col(),
                         novo = guess_ab_col(),
                         oflo = guess_ab_col(),
                         oxac = guess_ab_col(),
                         peni = guess_ab_col(),
                         pipe = guess_ab_col(),
                         pita = guess_ab_col(),
                         poly = guess_ab_col(),
                         pris = guess_ab_col(),
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
                         ...) {


  # support old `tbl` parameter
  if ("tbl" %in% names(list(...))) {
    x <- list(...)$tbl
  }

  tbl_ <- x

  if (!is.data.frame(tbl_)) {
    stop("`tbl_` must be a data frame.", call. = FALSE)
  }

  # try to find columns based on type
  # -- mo
  if (is.null(col_mo)) {
    col_mo <- search_type_in_df(tbl = tbl_, type = "mo")
  }
  if (is.null(col_mo)) {
    stop("`col_mo` must be set.", call. = FALSE)
  }

  if (!all(rules %in% c("breakpoints", "expert", "other", "all"))) {
    stop("`rules` must be one or more of:  'breakpoints', 'expert', 'other', 'all'.")
  }

  if (is.null(col_mo)) {
    stop("`col_mo` must be set")
  }

  warned <- FALSE

  txt_error <- function() { cat("", bgRed(black(" ERROR ")), "\n") }
  txt_warning <- function() { if (warned == FALSE) { cat("", bgYellow(black(" WARNING ")), "\n") }; warned <<- TRUE }
  txt_ok <- function(no_of_changes) {
    if (warned == FALSE) {
      if (no_of_changes > 0) {
        if (no_of_changes == 1) {
          cat(blue(" (1 new change)\n"))
        } else {
          cat(blue(paste0(" (", no_of_changes, " new changes)\n")))
        }
      } else {
        cat(green(" (no new changes)\n"))
      }
      warned <<- FALSE
    }
  }

  # check columns
  if (identical(amcl, as.name("guess_ab_col"))) { amcl <- guess_ab_col(tbl_, "amcl", verbose = verbose) }
  if (identical(amik, as.name("guess_ab_col"))) { amik <- guess_ab_col(tbl_, "amik", verbose = verbose) }
  if (identical(amox, as.name("guess_ab_col"))) { amox <- guess_ab_col(tbl_, "amox", verbose = verbose) }
  if (identical(ampi, as.name("guess_ab_col"))) { ampi <- guess_ab_col(tbl_, "ampi", verbose = verbose) }
  if (identical(azit, as.name("guess_ab_col"))) { azit <- guess_ab_col(tbl_, "azit", verbose = verbose) }
  if (identical(azlo, as.name("guess_ab_col"))) { azlo <- guess_ab_col(tbl_, "azlo", verbose = verbose) }
  if (identical(aztr, as.name("guess_ab_col"))) { aztr <- guess_ab_col(tbl_, "aztr", verbose = verbose) }
  if (identical(cefa, as.name("guess_ab_col"))) { cefa <- guess_ab_col(tbl_, "cefa", verbose = verbose) }
  if (identical(cfep, as.name("guess_ab_col"))) { cfep <- guess_ab_col(tbl_, "cfep", verbose = verbose) }
  if (identical(cfot, as.name("guess_ab_col"))) { cfot <- guess_ab_col(tbl_, "cfot", verbose = verbose) }
  if (identical(cfox, as.name("guess_ab_col"))) { cfox <- guess_ab_col(tbl_, "cfox", verbose = verbose) }
  if (identical(cfra, as.name("guess_ab_col"))) { cfra <- guess_ab_col(tbl_, "cfra", verbose = verbose) }
  if (identical(cfta, as.name("guess_ab_col"))) { cfta <- guess_ab_col(tbl_, "cfta", verbose = verbose) }
  if (identical(cftr, as.name("guess_ab_col"))) { cftr <- guess_ab_col(tbl_, "cftr", verbose = verbose) }
  if (identical(cfur, as.name("guess_ab_col"))) { cfur <- guess_ab_col(tbl_, "cfur", verbose = verbose) }
  if (identical(chlo, as.name("guess_ab_col"))) { chlo <- guess_ab_col(tbl_, "chlo", verbose = verbose) }
  if (identical(cipr, as.name("guess_ab_col"))) { cipr <- guess_ab_col(tbl_, "cipr", verbose = verbose) }
  if (identical(clar, as.name("guess_ab_col"))) { clar <- guess_ab_col(tbl_, "clar", verbose = verbose) }
  if (identical(clin, as.name("guess_ab_col"))) { clin <- guess_ab_col(tbl_, "clin", verbose = verbose) }
  if (identical(clox, as.name("guess_ab_col"))) { clox <- guess_ab_col(tbl_, "clox", verbose = verbose) }
  if (identical(coli, as.name("guess_ab_col"))) { coli <- guess_ab_col(tbl_, "coli", verbose = verbose) }
  if (identical(czol, as.name("guess_ab_col"))) { czol <- guess_ab_col(tbl_, "czol", verbose = verbose) }
  if (identical(dapt, as.name("guess_ab_col"))) { dapt <- guess_ab_col(tbl_, "dapt", verbose = verbose) }
  if (identical(doxy, as.name("guess_ab_col"))) { doxy <- guess_ab_col(tbl_, "doxy", verbose = verbose) }
  if (identical(erta, as.name("guess_ab_col"))) { erta <- guess_ab_col(tbl_, "erta", verbose = verbose) }
  if (identical(eryt, as.name("guess_ab_col"))) { eryt <- guess_ab_col(tbl_, "eryt", verbose = verbose) }
  if (identical(fosf, as.name("guess_ab_col"))) { fosf <- guess_ab_col(tbl_, "fosf", verbose = verbose) }
  if (identical(fusi, as.name("guess_ab_col"))) { fusi <- guess_ab_col(tbl_, "fusi", verbose = verbose) }
  if (identical(gent, as.name("guess_ab_col"))) { gent <- guess_ab_col(tbl_, "gent", verbose = verbose) }
  if (identical(imip, as.name("guess_ab_col"))) { imip <- guess_ab_col(tbl_, "imip", verbose = verbose) }
  if (identical(kana, as.name("guess_ab_col"))) { kana <- guess_ab_col(tbl_, "kana", verbose = verbose) }
  if (identical(levo, as.name("guess_ab_col"))) { levo <- guess_ab_col(tbl_, "levo", verbose = verbose) }
  if (identical(linc, as.name("guess_ab_col"))) { linc <- guess_ab_col(tbl_, "linc", verbose = verbose) }
  if (identical(line, as.name("guess_ab_col"))) { line <- guess_ab_col(tbl_, "line", verbose = verbose) }
  if (identical(mero, as.name("guess_ab_col"))) { mero <- guess_ab_col(tbl_, "mero", verbose = verbose) }
  if (identical(mezl, as.name("guess_ab_col"))) { mezl <- guess_ab_col(tbl_, "mezl", verbose = verbose) }
  if (identical(mino, as.name("guess_ab_col"))) { mino <- guess_ab_col(tbl_, "mino", verbose = verbose) }
  if (identical(moxi, as.name("guess_ab_col"))) { moxi <- guess_ab_col(tbl_, "moxi", verbose = verbose) }
  if (identical(nali, as.name("guess_ab_col"))) { nali <- guess_ab_col(tbl_, "nali", verbose = verbose) }
  if (identical(neom, as.name("guess_ab_col"))) { neom <- guess_ab_col(tbl_, "neom", verbose = verbose) }
  if (identical(neti, as.name("guess_ab_col"))) { neti <- guess_ab_col(tbl_, "neti", verbose = verbose) }
  if (identical(nitr, as.name("guess_ab_col"))) { nitr <- guess_ab_col(tbl_, "nitr", verbose = verbose) }
  if (identical(norf, as.name("guess_ab_col"))) { norf <- guess_ab_col(tbl_, "norf", verbose = verbose) }
  if (identical(novo, as.name("guess_ab_col"))) { novo <- guess_ab_col(tbl_, "novo", verbose = verbose) }
  if (identical(oflo, as.name("guess_ab_col"))) { oflo <- guess_ab_col(tbl_, "oflo", verbose = verbose) }
  if (identical(oxac, as.name("guess_ab_col"))) { oxac <- guess_ab_col(tbl_, "oxac", verbose = verbose) }
  if (identical(peni, as.name("guess_ab_col"))) { peni <- guess_ab_col(tbl_, "peni", verbose = verbose) }
  if (identical(pipe, as.name("guess_ab_col"))) { pipe <- guess_ab_col(tbl_, "pipe", verbose = verbose) }
  if (identical(pita, as.name("guess_ab_col"))) { pita <- guess_ab_col(tbl_, "pita", verbose = verbose) }
  if (identical(poly, as.name("guess_ab_col"))) { poly <- guess_ab_col(tbl_, "poly", verbose = verbose) }
  if (identical(pris, as.name("guess_ab_col"))) { pris <- guess_ab_col(tbl_, "pris", verbose = verbose) }
  if (identical(qida, as.name("guess_ab_col"))) { qida <- guess_ab_col(tbl_, "qida", verbose = verbose) }
  if (identical(rifa, as.name("guess_ab_col"))) { rifa <- guess_ab_col(tbl_, "rifa", verbose = verbose) }
  if (identical(roxi, as.name("guess_ab_col"))) { roxi <- guess_ab_col(tbl_, "roxi", verbose = verbose) }
  if (identical(siso, as.name("guess_ab_col"))) { siso <- guess_ab_col(tbl_, "siso", verbose = verbose) }
  if (identical(teic, as.name("guess_ab_col"))) { teic <- guess_ab_col(tbl_, "teic", verbose = verbose) }
  if (identical(tetr, as.name("guess_ab_col"))) { tetr <- guess_ab_col(tbl_, "tetr", verbose = verbose) }
  if (identical(tica, as.name("guess_ab_col"))) { tica <- guess_ab_col(tbl_, "tica", verbose = verbose) }
  if (identical(tige, as.name("guess_ab_col"))) { tige <- guess_ab_col(tbl_, "tige", verbose = verbose) }
  if (identical(tobr, as.name("guess_ab_col"))) { tobr <- guess_ab_col(tbl_, "tobr", verbose = verbose) }
  if (identical(trim, as.name("guess_ab_col"))) { trim <- guess_ab_col(tbl_, "trim", verbose = verbose) }
  if (identical(trsu, as.name("guess_ab_col"))) { trsu <- guess_ab_col(tbl_, "trsu", verbose = verbose) }
  if (identical(vanc, as.name("guess_ab_col"))) { vanc <- guess_ab_col(tbl_, "vanc", verbose = verbose) }
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
  col.list <- check_available_columns(tbl = tbl_, col.list = col.list, info = info)
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

  ab_missing <- function(ab) {
    all(ab %in% c(NULL, NA))
  }

  verbose_info <- data.frame(row = integer(0),
                             col = character(0),
                             mo_fullname = character(0),
                             old = character(0),
                             new = character(0),
                             rule = character(0),
                             rule_group = character(0),
                             rule_name = character(0),
                             stringsAsFactors = FALSE)

  # helper function for editing the table
  edit_rsi <- function(to, rule, rows, cols) {
    cols <- unique(cols[!is.na(cols) & !is.null(cols)])
    if (length(rows) > 0 & length(cols) > 0) {
      before_df <- tbl_original
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

      tbl_[rows, cols] <<- tbl_original[rows, cols]

      after <- as.character(unlist(as.list(tbl_original[rows, cols])))

      # before_df might not be a data.frame, but a tibble of data.table instead
      old <- as.data.frame(before_df, stringsAsFactors = FALSE)[rows,]
      no_of_changes_this_run <- 0
      for (i in 1:length(cols)) {
        verbose_new <- data.frame(row = rows,
                                  col = cols[i],
                                  mo_fullname = tbl_[rows, "fullname"],
                                  old = as.character(old[, cols[i]]),
                                  new = as.character(tbl_[rows, cols[i]]),
                                  rule = strip_style(rule[1]),
                                  rule_group = strip_style(rule[2]),
                                  rule_name = strip_style(rule[3]),
                                  stringsAsFactors = FALSE)
        colnames(verbose_new) <- c("row", "col", "mo_fullname", "old", "new", "rule", "rule_group", "rule_name")
        verbose_new <- verbose_new %>% filter(old != new | is.na(old))
        verbose_info <<- rbind(verbose_info, verbose_new)
        no_of_changes_this_run <- no_of_changes_this_run + nrow(verbose_new)
      }
      # return number of (new) changes
      return(no_of_changes_this_run)
    }
    # return number of (new) changes: none.
    return(0)
  }

  # save original table
  tbl_original <- tbl_

  # join to microorganisms data set
  suppressWarnings(
    tbl_ <- tbl_ %>%
      mutate_at(vars(col_mo), as.mo) %>%
      left_join_microorganisms(by = col_mo, suffix = c("_oldcols", "")) %>%
      mutate(gramstain = mo_gramstain(pull(., col_mo), language = "en"),
             genus_species = paste(genus, species)) %>%
      as.data.frame(stringsAsFactors = FALSE)
  )

  if (info == TRUE) {
    cat(paste0(
      "\nRules by the ", bold("European Committee on Antimicrobial Susceptibility Testing (EUCAST)"),
      "\n", blue("http://eucast.org/"), "\n"))
  }

  # since ampicillin ^= amoxicillin, get the first from the latter (not in original EUCAST table)
  if (!ab_missing(ampi) & !ab_missing(amox)) {
    if (verbose == TRUE) {
      cat("\n VERBOSE: transforming",
          length(which(tbl_[, amox] == "S" & !tbl_[, ampi] %in% c("S", "I", "R"))),
          "empty ampicillin fields to 'S' based on amoxicillin. ")
      cat("\n VERBOSE: transforming",
          length(which(tbl_[, amox] == "I" & !tbl_[, ampi] %in% c("S", "I", "R"))),
          "empty ampicillin fields to 'I' based on amoxicillin. ")
      cat("\n VERBOSE: transforming",
          length(which(tbl_[, amox] == "R" & !tbl_[, ampi] %in% c("S", "I", "R"))),
          "empty ampicillin fields to 'R' based on amoxicillin. \n")
    }
    tbl_[which(tbl_[, amox] == "S" & !tbl_[, ampi] %in% c("S", "I", "R")), ampi] <- "S"
    tbl_[which(tbl_[, amox] == "I" & !tbl_[, ampi] %in% c("S", "I", "R")), ampi] <- "I"
    tbl_[which(tbl_[, amox] == "R" & !tbl_[, ampi] %in% c("S", "I", "R")), ampi] <- "R"
  } else if (ab_missing(ampi) & !ab_missing(amox)) {
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
  cephalosporins_without_cfta <- cephalosporins[cephalosporins != ifelse(is.null(cfta), "", cfta)]
  carbapenems <- c(erta, imip, mero)
  aminopenicillins <- c(ampi, amox)
  ureidopenicillins <- c(pipe, pita, azlo, mezl)
  fluoroquinolones <- c(oflo, cipr, norf, levo, moxi)
  all_betalactams <- c(aminopenicillins, ureidopenicillins, cephalosporins, carbapenems, amcl, oxac, clox, peni)

  # Help function to get available antibiotic column names ------------------
  get_antibiotic_columns <- function(x, df) {
    x <- trimws(unlist(strsplit(x, ",", fixed = TRUE)))
    y <- character(0)
    for (i in 1:length(x)) {
      y <- c(y, tryCatch(get(x[i]), error = function(e) ""))
    }
    y[y != "" & y %in% colnames(df)]
  }

  eucast_rules_df <- eucast_rules_file()
  no_of_changes <- 0
  for (i in 1:nrow(eucast_rules_df)) {

    rule_previous <- eucast_rules_df[max(1, i - 1), "reference.rule"]
    rule_current <- eucast_rules_df[i, "reference.rule"]
    rule_next <- eucast_rules_df[min(nrow(eucast_rules_df), i + 1), "reference.rule"]
    rule_group_previous <- eucast_rules_df[max(1, i - 1), "reference.rule_group"]
    rule_group_current <- eucast_rules_df[i, "reference.rule_group"]
    rule_group_next <- eucast_rules_df[min(nrow(eucast_rules_df), i + 1), "reference.rule_group"]
    if (is.na(eucast_rules_df[i, 4])) {
      rule_text <- paste(eucast_rules_df[i, 6], "=", eucast_rules_df[i, 7])
    } else {
      rule_text <- paste("if", eucast_rules_df[i, 4], "=", eucast_rules_df[i, 5],
                         "then", eucast_rules_df[i, 6], "=", eucast_rules_df[i, 7])
    }
    if (i == 1) {
      rule_previous <- ""
      rule_group_previous <- ""
    }
    if (i == nrow(eucast_rules_df)) {
      rule_next <- ""
      rule_group_next <- ""
    }

    # don't apply rules if user doesn't want to apply them
    if (rule_group_current %like% "breakpoint" & !any(c("all", "breakpoints") %in% rules)) {
      next
    }
    if (rule_group_current %like% "expert" & !any(c("all", "expert") %in% rules)) {
      next
    }
    if (rule_group_current %like% "other" & !any(c("all", "other") %in% rules)) {
      next
    }


    if (info == TRUE) {
      # Print rule (group) ------------------------------------------------------
      if (rule_group_current != rule_group_previous) {
        # is new rule group, one of Breakpoints, Expert Rules and Other
        cat(bold(
          case_when(
            rule_group_current %like% "breakpoint" ~
              paste0("\nEUCAST Clinical Breakpoints (v", EUCAST_VERSION_BREAKPOINTS, ")\n"),
            rule_group_current %like% "expert" ~
              paste0("\nEUCAST Expert Rules, Intrinsic Resistance and Exceptional Phenotypes (v", EUCAST_VERSION_EXPERT_RULES, ")\n"),
            TRUE ~
              "\nOther rules\n"
          )
        ))
      }
      # Print rule  -------------------------------------------------------------
      if (rule_current != rule_previous) {
        # is new rule within group, print its name
        if (rule_current %in% c(AMR::microorganisms$family,
                                AMR::microorganisms$fullname)) {
          cat(italic(rule_current))
        } else {
          cat(rule_current)
        }
        warned <- FALSE
      }
    }

    # Get rule from file ------------------------------------------------------
    col_mo_property <- eucast_rules_df[i, 1]
    like_is_one_of <- eucast_rules_df[i, 2]

    # be sure to comprise all coagulase-negative/-positive Staphylococci when they are mentioned
    if (eucast_rules_df[i, 3] %like% "coagulase-") {
      suppressWarnings(
        all_staph <- AMR::microorganisms %>%
          filter(genus == "Staphylococcus") %>%
          mutate(CNS_CPS = mo_fullname(mo, Becker = "all"))
      )
      if (eucast_rules_df[i, 3] %like% "coagulase-") {
        eucast_rules_df[i, 3] <- paste0("^(",
                                        paste0(all_staph %>%
                                                 filter(CNS_CPS %like% "coagulase-negative") %>%
                                                 pull(fullname),
                                               collapse = "|"),
                                        ")$")
      } else {
        eucast_rules_df[i, 3] <- paste0("^(",
                                        paste0(all_staph %>%
                                                 filter(CNS_CPS %like% "coagulase-positive") %>%
                                                 pull(fullname),
                                               collapse = "|"),
                                        ")$")
      }
      like_is_one_of <- "like"
    }

    if (like_is_one_of == "is") {
      mo_value <- paste0("^", eucast_rules_df[i, 3], "$")
    } else if (like_is_one_of == "one_of") {
      # "Clostridium, Actinomyces, ..." -> "^(Clostridium|Actinomyces|...)$"
      mo_value <- paste0("^(",
                         paste(trimws(unlist(strsplit(eucast_rules_df[i, 3], ",", fixed = TRUE))),
                               collapse = "|"),
                         ")$")
    } else if (like_is_one_of == "like") {
      mo_value <- eucast_rules_df[i, 3]
    } else {
      stop("invalid like_is_one_of", call. = FALSE)
    }

    source_antibiotics <- eucast_rules_df[i, 4]
    source_value <- trimws(unlist(strsplit(eucast_rules_df[i, 5], ",", fixed = TRUE)))
    target_antibiotics <- eucast_rules_df[i, 6]
    target_value <- eucast_rules_df[i, 7]

    if (is.na(source_antibiotics)) {
      rows <- tryCatch(which(tbl_[, col_mo_property] %like% mo_value),
                       error = function(e) integer(0))
    } else {
      source_antibiotics <- get_antibiotic_columns(source_antibiotics, tbl_)
      if (length(source_value) == 1 & length(source_antibiotics) > 1) {
        source_value <- rep(source_value, length(source_antibiotics))
      }
      if (length(source_antibiotics) == 0) {
        rows <- integer(0)
      } else if (length(source_antibiotics) == 1) {
        rows <-  tryCatch(which(tbl_[, col_mo_property] %like% mo_value
                                & tbl_[, source_antibiotics[1L]] == source_value[1L]),
                          error = function(e) integer(0))
      } else if (length(source_antibiotics) == 2) {
        rows <-  tryCatch(which(tbl_[, col_mo_property] %like% mo_value
                                & tbl_[, source_antibiotics[1L]] == source_value[1L]
                                & tbl_[, source_antibiotics[2L]] == source_value[2L]),
                          error = function(e) integer(0))
      } else if (length(source_antibiotics) == 3) {
        rows <-  tryCatch(which(tbl_[, col_mo_property] %like% mo_value
                                & tbl_[, source_antibiotics[1L]] == source_value[1L]
                                & tbl_[, source_antibiotics[2L]] == source_value[2L]
                                & tbl_[, source_antibiotics[3L]] == source_value[3L]),
                          error = function(e) integer(0))
      } else {
        stop("only 3 antibiotics supported for source_antibiotics ", call. = FALSE)
      }
    }

    cols <- get_antibiotic_columns(target_antibiotics, tbl_)

    # Apply rule on data ------------------------------------------------------
    # this will return the unique number of changes
    no_of_changes <- no_of_changes + edit_rsi(to = target_value,
                                              rule = c(rule_text, rule_group_current, rule_current),
                                              rows = rows,
                                              cols = cols)

    # Print number of new changes ---------------------------------------------
    if (info == TRUE & rule_next != rule_current) {
      # print only on last one of rules in this group
      txt_ok(no_of_changes = no_of_changes)
      no_of_changes <- 0
    }
  }

  # Print overview ----------------------------------------------------------
  if (info == TRUE) {
    if (verbose == TRUE) {
      wouldve <- "would have "
    } else {
      wouldve <- ""
    }

    verbose_info <- verbose_info %>%
      arrange(row, rule_group, rule_name, col)

    decimal.mark <- getOption("OutDec")
    big.mark <- ifelse(decimal.mark != ",", ",", ".")
    formatnr <- function(x) {
      trimws(format(x, big.mark = big.mark, decimal.mark = decimal.mark))
    }

    cat(paste0("\n", silver(strrep("-", options()$width - 1)), "\n"))
    cat(bold(paste('EUCAST rules', paste0(wouldve, 'affected'),
                   formatnr(n_distinct(verbose_info$row)),
                   'out of', formatnr(nrow(tbl_original)),
                   'rows, making a total of', formatnr(nrow(verbose_info)), 'edits\n')))

    # print added values ----
    if (verbose_info %>% filter(is.na(old)) %>% nrow() == 0) {
      colour <- cat # is function
    } else {
      colour <- blue # is function
    }
    cat(colour(paste0("=> ", wouldve, "added ",
                      bold(formatnr(verbose_info %>%
                                      filter(is.na(old)) %>%
                                      nrow()), "test results"),
                      "\n")))
    if (verbose_info %>% filter(is.na(old)) %>% nrow() > 0) {
      verbose_info %>%
        filter(is.na(old)) %>%
        # sort it well: S < I < R
        mutate(new = as.rsi(new)) %>%
        group_by(new) %>%
        summarise(n = n()) %>%
        mutate(plural = ifelse(n > 1, "s", ""),
               txt = paste0(formatnr(n), " test result", plural, " added as ", new)) %>%
        pull(txt) %>%
        paste("   -", ., collapse = "\n") %>%
        cat()
    }

    # print changed values ----
    if (verbose_info %>% filter(!is.na(old)) %>% nrow() == 0) {
      colour <- cat # is function
    } else {
      colour <- blue # is function
    }
    cat(colour(paste0("\n=> ", wouldve, "changed ",
                      bold(formatnr(verbose_info %>%
                                      filter(!is.na(old)) %>%
                                      nrow()), "test results"),
                      "\n")))
    if (verbose_info %>% filter(!is.na(old)) %>% nrow() > 0) {
      verbose_info %>%
        filter(!is.na(old)) %>%
        # sort it well: S < I < R
        mutate(old = as.rsi(old),
               new = as.rsi(new)) %>%
        group_by(old, new) %>%
        summarise(n = n()) %>%
        mutate(plural = ifelse(n > 1, "s", ""),
               txt = paste0(formatnr(n), " test result", plural, " changed from ", old, " to ", new)) %>%
        pull(txt) %>%
        paste("   -", ., collapse = "\n") %>%
        cat()
      cat("\n")
    }
    cat(paste0(silver(strrep("-", options()$width - 1)), "\n"))

    if (verbose == FALSE & nrow(verbose_info) > 0) {
      cat(paste("\nUse", bold("verbose = TRUE"), "to get a data.frame with all specified edits instead.\n"))
    }
  }

  # Return data set ---------------------------------------------------------
  if (verbose == TRUE) {
    verbose_info
  } else {
    tbl_original
  }
}

#' @rdname eucast_rules
#' @importFrom dplyr %>% arrange
#' @export
eucast_rules_file <- function() {
  utils::read.delim(file = EUCAST_RULES_FILE_LOCATION,
                    sep = "\t",
                    stringsAsFactors = FALSE,
                    header = TRUE,
                    strip.white = TRUE,
                    na = c(NA, "", NULL)) %>%
    arrange(reference.rule_group,
            reference.rule)
}
