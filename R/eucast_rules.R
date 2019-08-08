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
EUCAST_VERSION_BREAKPOINTS <- "9.0, 2019"
EUCAST_VERSION_EXPERT_RULES <- "3.1, 2016"

#' EUCAST rules
#'
#' Apply susceptibility rules as defined by the European Committee on Antimicrobial Susceptibility Testing (EUCAST, \url{http://eucast.org}), see \emph{Source}. This includes (1) expert rules, (2) intrinsic resistance and (3) inferred resistance as defined in their breakpoint tables.
#' @param x data with antibiotic columns, like e.g. \code{AMX} and \code{AMC}
#' @param info print progress
#' @param rules a character vector that specifies which rules should be applied - one or more of \code{c("breakpoints", "expert", "other", "all")}
#' @param verbose a logical to turn Verbose mode on and off (default is off). In Verbose mode, the function does not apply rules to the data, but instead returns a \code{data.frame} with extensive info about which rows and columns would be effected and in which way.
#' @param ... column name of an antibiotic, see section Antibiotics
#' @inheritParams first_isolate
#' @details
#' \strong{Note:} This function does not translate MIC values to RSI values. Use \code{\link{as.rsi}} for that. \cr
#' \strong{Note:} When ampicillin (AMP, J01CA01) is not available but amoxicillin (AMX, J01CA04) is, the latter will be used for all rules where there is a dependency on ampicillin. These drugs are interchangeable when it comes to expression of antimicrobial resistance.
#'
#' The file containing all EUCAST rules is located here: \url{https://gitlab.com/msberends/AMR/blob/master/data-raw/eucast_rules.tsv}.
#'
#' @section Antibiotics:
#' To define antibiotics column names, leave as it is to determine it automatically with \code{\link{guess_ab_col}} or input a text (case-insensitive), or use \code{NULL} to skip a column (e.g. \code{TIC = NULL} to skip ticarcillin). Manually defined but non-existing columns will be skipped with a warning.
#'
#' The following antibiotics are used for the functions \code{\link{eucast_rules}} and \code{\link{mdro}}. These are shown in the format '\strong{antimicrobial ID}: name (\emph{ATC code})', sorted by name:
#'
#'  \strong{AMK}: amikacin (\href{https://www.whocc.no/atc_ddd_index/?code=J01GB06}{J01GB06}),
#'  \strong{AMX}: amoxicillin (\href{https://www.whocc.no/atc_ddd_index/?code=J01CA04}{J01CA04}),
#'  \strong{AMC}: amoxicillin/clavulanic acid (\href{https://www.whocc.no/atc_ddd_index/?code=J01CR02}{J01CR02}),
#'  \strong{AMP}: ampicillin (\href{https://www.whocc.no/atc_ddd_index/?code=J01CA01}{J01CA01}),
#'  \strong{AZM}: azithromycin (\href{https://www.whocc.no/atc_ddd_index/?code=J01FA10}{J01FA10}),
#'  \strong{AZL}: azlocillin (\href{https://www.whocc.no/atc_ddd_index/?code=J01CA09}{J01CA09}),
#'  \strong{ATM}: aztreonam (\href{https://www.whocc.no/atc_ddd_index/?code=J01DF01}{J01DF01}),
#'  \strong{CAP}: capreomycin (\href{https://www.whocc.no/atc_ddd_index/?code=J04AB30}{J04AB30}),
#'  \strong{RID}: cefaloridine (\href{https://www.whocc.no/atc_ddd_index/?code=J01DB02}{J01DB02}),
#'  \strong{CZO}: cefazolin (\href{https://www.whocc.no/atc_ddd_index/?code=J01DB04}{J01DB04}),
#'  \strong{FEP}: cefepime (\href{https://www.whocc.no/atc_ddd_index/?code=J01DE01}{J01DE01}),
#'  \strong{CTX}: cefotaxime (\href{https://www.whocc.no/atc_ddd_index/?code=J01DD01}{J01DD01}),
#'  \strong{FOX}: cefoxitin (\href{https://www.whocc.no/atc_ddd_index/?code=J01DC01}{J01DC01}),
#'  \strong{CED}: cefradine (\href{https://www.whocc.no/atc_ddd_index/?code=J01DB09}{J01DB09}),
#'  \strong{CAZ}: ceftazidime (\href{https://www.whocc.no/atc_ddd_index/?code=J01DD02}{J01DD02}),
#'  \strong{CRO}: ceftriaxone (\href{https://www.whocc.no/atc_ddd_index/?code=J01DD04}{J01DD04}),
#'  \strong{CXM}: cefuroxime (\href{https://www.whocc.no/atc_ddd_index/?code=J01DC02}{J01DC02}),
#'  \strong{CHL}: chloramphenicol (\href{https://www.whocc.no/atc_ddd_index/?code=J01BA01}{J01BA01}),
#'  \strong{CIP}: ciprofloxacin (\href{https://www.whocc.no/atc_ddd_index/?code=J01MA02}{J01MA02}),
#'  \strong{CLR}: clarithromycin (\href{https://www.whocc.no/atc_ddd_index/?code=J01FA09}{J01FA09}),
#'  \strong{CLI}: clindamycin (\href{https://www.whocc.no/atc_ddd_index/?code=J01FF01}{J01FF01}),
#'  \strong{COL}: colistin (\href{https://www.whocc.no/atc_ddd_index/?code=J01XB01}{J01XB01}),
#'  \strong{DAP}: daptomycin (\href{https://www.whocc.no/atc_ddd_index/?code=J01XX09}{J01XX09}),
#'  \strong{DOX}: doxycycline (\href{https://www.whocc.no/atc_ddd_index/?code=J01AA02}{J01AA02}),
#'  \strong{ETP}: ertapenem (\href{https://www.whocc.no/atc_ddd_index/?code=J01DH03}{J01DH03}),
#'  \strong{ERY}: erythromycin (\href{https://www.whocc.no/atc_ddd_index/?code=J01FA01}{J01FA01}),
#'  \strong{ETH}: ethambutol (\href{https://www.whocc.no/atc_ddd_index/?code=J04AK02}{J04AK02}),
#'  \strong{FLC}: flucloxacillin (\href{https://www.whocc.no/atc_ddd_index/?code=J01CF05}{J01CF05}),
#'  \strong{FOS}: fosfomycin (\href{https://www.whocc.no/atc_ddd_index/?code=J01XX01}{J01XX01}),
#'  \strong{FUS}: fusidic acid (\href{https://www.whocc.no/atc_ddd_index/?code=J01XC01}{J01XC01}),
#'  \strong{GAT}: gatifloxacin (\href{https://www.whocc.no/atc_ddd_index/?code=J01MA16}{J01MA16}),
#'  \strong{GEN}: gentamicin (\href{https://www.whocc.no/atc_ddd_index/?code=J01GB03}{J01GB03}),
#'  \strong{IPM}: imipenem (\href{https://www.whocc.no/atc_ddd_index/?code=J01DH51}{J01DH51}),
#'  \strong{INH}: isoniazid (\href{https://www.whocc.no/atc_ddd_index/?code=J04AC01}{J04AC01}),
#'  \strong{KAN}: kanamycin (\href{https://www.whocc.no/atc_ddd_index/?code=J01GB04}{J01GB04}),
#'  \strong{LVX}: levofloxacin (\href{https://www.whocc.no/atc_ddd_index/?code=J01MA12}{J01MA12}),
#'  \strong{LIN}: lincomycin (\href{https://www.whocc.no/atc_ddd_index/?code=J01FF02}{J01FF02}),
#'  \strong{LNZ}: linezolid (\href{https://www.whocc.no/atc_ddd_index/?code=J01XX08}{J01XX08}),
#'  \strong{MEM}: meropenem (\href{https://www.whocc.no/atc_ddd_index/?code=J01DH02}{J01DH02}),
#'  \strong{MTR}: metronidazole (\href{https://www.whocc.no/atc_ddd_index/?code=J01XD01}{J01XD01}),
#'  \strong{MEZ}: mezlocillin (\href{https://www.whocc.no/atc_ddd_index/?code=J01CA10}{J01CA10}),
#'  \strong{MNO}: minocycline (\href{https://www.whocc.no/atc_ddd_index/?code=J01AA08}{J01AA08}),
#'  \strong{MFX}: moxifloxacin (\href{https://www.whocc.no/atc_ddd_index/?code=J01MA14}{J01MA14}),
#'  \strong{NAL}: nalidixic acid (\href{https://www.whocc.no/atc_ddd_index/?code=J01MB02}{J01MB02}),
#'  \strong{NEO}: neomycin (\href{https://www.whocc.no/atc_ddd_index/?code=J01GB05}{J01GB05}),
#'  \strong{NET}: netilmicin (\href{https://www.whocc.no/atc_ddd_index/?code=J01GB07}{J01GB07}),
#'  \strong{NIT}: nitrofurantoin (\href{https://www.whocc.no/atc_ddd_index/?code=J01XE01}{J01XE01}),
#'  \strong{NOR}: norfloxacin (\href{https://www.whocc.no/atc_ddd_index/?code=J01MA06}{J01MA06}),
#'  \strong{NOV}: novobiocin (an ATCvet code: \href{https://www.whocc.no/atc_ddd_index/?code=QJ01XX95}{QJ01XX95}),
#'  \strong{OFX}: ofloxacin (\href{https://www.whocc.no/atc_ddd_index/?code=J01MA01}{J01MA01}),
#'  \strong{OXA}: oxacillin (\href{https://www.whocc.no/atc_ddd_index/?code=J01CF04}{J01CF04}),
#'  \strong{PEN}: penicillin G (\href{https://www.whocc.no/atc_ddd_index/?code=J01CE01}{J01CE01}),
#'  \strong{PIP}: piperacillin (\href{https://www.whocc.no/atc_ddd_index/?code=J01CA12}{J01CA12}),
#'  \strong{TZP}: piperacillin/tazobactam (\href{https://www.whocc.no/atc_ddd_index/?code=J01CR05}{J01CR05}),
#'  \strong{PLB}: polymyxin B (\href{https://www.whocc.no/atc_ddd_index/?code=J01XB02}{J01XB02}),
#'  \strong{PRI}: pristinamycin (\href{https://www.whocc.no/atc_ddd_index/?code=J01FG01}{J01FG01}),
#'  \strong{PZA}: pyrazinamide (\href{https://www.whocc.no/atc_ddd_index/?code=J04AK01}{J04AK01}),
#'  \strong{QDA}: quinupristin/dalfopristin (\href{https://www.whocc.no/atc_ddd_index/?code=J01FG02}{J01FG02}),
#'  \strong{RIB}: rifabutin (\href{https://www.whocc.no/atc_ddd_index/?code=J04AB04}{J04AB04}),
#'  \strong{RIF}: rifampicin (\href{https://www.whocc.no/atc_ddd_index/?code=J04AB02}{J04AB02}),
#'  \strong{RIF}: rifampin (\href{https://www.whocc.no/atc_ddd_index/?code=J04AB02}{J04AB02}),
#'  \strong{RFP}: rifapentine (\href{https://www.whocc.no/atc_ddd_index/?code=J04AB05}{J04AB05}),
#'  \strong{RXT}: roxithromycin (\href{https://www.whocc.no/atc_ddd_index/?code=J01FA06}{J01FA06}),
#'  \strong{SIS}: sisomicin (\href{https://www.whocc.no/atc_ddd_index/?code=J01GB08}{J01GB08}),
#'  \strong{TEC}: teicoplanin (\href{https://www.whocc.no/atc_ddd_index/?code=J01XA02}{J01XA02}),
#'  \strong{TCY}: tetracycline (\href{https://www.whocc.no/atc_ddd_index/?code=J01AA07}{J01AA07}),
#'  \strong{TIC}: ticarcillin (\href{https://www.whocc.no/atc_ddd_index/?code=J01CA13}{J01CA13}),
#'  \strong{TGC}: tigecycline (\href{https://www.whocc.no/atc_ddd_index/?code=J01AA12}{J01AA12}),
#'  \strong{TOB}: tobramycin (\href{https://www.whocc.no/atc_ddd_index/?code=J01GB01}{J01GB01}),
#'  \strong{TMP}: trimethoprim (\href{https://www.whocc.no/atc_ddd_index/?code=J01EA01}{J01EA01}),
#'  \strong{SXT}: trimethoprim/sulfamethoxazole (\href{https://www.whocc.no/atc_ddd_index/?code=J01EE01}{J01EE01}),
#'  \strong{VAN}: vancomycin (\href{https://www.whocc.no/atc_ddd_index/?code=J01XA01}{J01XA01}).
#' @keywords interpretive eucast reading resistance
#' @rdname eucast_rules
#' @export
#' @importFrom dplyr %>% select pull mutate_at vars group_by summarise n
#' @importFrom crayon bold bgGreen bgYellow bgRed black green blue italic strip_style white red
#' @importFrom utils menu
#' @return The input of \code{x}, possibly with edited values of antibiotics. Or, if \code{verbose = TRUE}, a \code{data.frame} with all original and new values of the affected bug-drug combinations.
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
#' @inheritSection AMR Read more on our website!
#' @examples
#' a <- data.frame(mo = c("Staphylococcus aureus",
#'                        "Enterococcus faecalis",
#'                        "Escherichia coli",
#'                        "Klebsiella pneumoniae",
#'                        "Pseudomonas aeruginosa"),
#'                 VAN = "-",       # Vancomycin
#'                 AMX = "-",       # Amoxicillin
#'                 COL = "-",       # Colistin
#'                 CAZ = "-",       # Ceftazidime
#'                 CXM = "-",       # Cefuroxime
#'                 PEN = "S",       # Penicillin G
#'                 FOX = "S",       # Cefoxitin
#'                 stringsAsFactors = FALSE)
#'
#' a
#' #                       mo  VAN  AMX  COL  CAZ  CXM  PEN  FOX
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
#' #                       mo  VAN  AMX  COL  CAZ  CXM  PEN  FOX
#' # 1  Staphylococcus aureus    -    S    R    R    S    S    S
#' # 2  Enterococcus faecalis    -    -    R    R    R    S    R
#' # 3       Escherichia coli    R    -    -    -    -    R    S
#' # 4  Klebsiella pneumoniae    R    R    -    -    -    R    S
#' # 5 Pseudomonas aeruginosa    R    R    -    -    R    R    R
#'
#'
#' # do not apply EUCAST rules, but rather get a data.frame
#' # with 18 rows, containing all details about the transformations:
#' c <- eucast_rules(a, verbose = TRUE)
eucast_rules <- function(x,
                         col_mo = NULL,
                         info = TRUE,
                         rules = c("breakpoints", "expert", "other", "all"),
                         verbose = FALSE,
                         ...) {
  
  if (verbose == TRUE & interactive()) {
    txt <- paste0("WARNING: In Verbose mode, the eucast_rules() function does not apply rules to the data, but instead returns a data set in logbook form: with extensive info about which rows and columns would be effected and in which way.",
                  "\n\nThis may overwrite your existing data if you use e.g.:",
                  "\ndata <- eucast_rules(data, verbose = TRUE)\n\nDo you want to continue?")
    if ("rstudioapi" %in% rownames(installed.packages())) {
      q_continue <- rstudioapi::showQuestion("Using verbose = TRUE with eucast_rules()", txt)
    } else {
      q_continue <- menu(choices = c("OK", "Cancel"), graphics = TRUE, title = txt)
    }
    if (q_continue %in% c(FALSE, 2)) {
      return(invisible())
    }
  }
  
  if (!is.data.frame(x)) {
    stop("`x` must be a data frame.", call. = FALSE)
  }

  # try to find columns based on type
  # -- mo
  if (is.null(col_mo)) {
    col_mo <- search_type_in_df(x = x, type = "mo")
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

  decimal.mark <- getOption("OutDec")
  big.mark <- ifelse(decimal.mark != ",", ",", ".")
  formatnr <- function(x) {
    trimws(format(x, big.mark = big.mark, decimal.mark = decimal.mark))
  }

  warned <- FALSE

  txt_error <- function() { cat("", bgRed(white(" ERROR ")), "\n\n") }
  txt_warning <- function() { if (warned == FALSE) { cat("", bgYellow(black(" WARNING "))) }; warned <<- TRUE }
  txt_ok <- function(no_of_changes) {
    if (warned == FALSE) {
      if (no_of_changes > 0) {
        if (no_of_changes == 1) {
          cat(blue(" (1 value changed)\n"))
        } else {
          cat(blue(paste0(" (", formatnr(no_of_changes), " values changed)\n")))
        }
      } else {
        cat(green(" (no values changed)\n"))
      }
      warned <<- FALSE
    }
  }

  cols_ab <- get_column_abx(x = x,
                            soft_dependencies = c("AMC",
                                                  "AMK",
                                                  "AMX",
                                                  "AMP",
                                                  "AZM",
                                                  "AZL",
                                                  "ATM",
                                                  "RID",
                                                  "FEP",
                                                  "CTX",
                                                  "FOX",
                                                  "CED",
                                                  "CAZ",
                                                  "CRO",
                                                  "CXM",
                                                  "CHL",
                                                  "CIP",
                                                  "CLR",
                                                  "CLI",
                                                  "FLC",
                                                  "COL",
                                                  "CZO",
                                                  "DAP",
                                                  "DOX",
                                                  "ETP",
                                                  "ERY",
                                                  "FOS",
                                                  "FUS",
                                                  "GEN",
                                                  "IPM",
                                                  "KAN",
                                                  "LVX",
                                                  "LIN",
                                                  "LNZ",
                                                  "MEM",
                                                  "MEZ",
                                                  "MNO",
                                                  "MFX",
                                                  "NAL",
                                                  "NEO",
                                                  "NET",
                                                  "NIT",
                                                  "NOR",
                                                  "NOV",
                                                  "OFX",
                                                  "OXA",
                                                  "PEN",
                                                  "PIP",
                                                  "TZP",
                                                  "PLB",
                                                  "PRI",
                                                  "QDA",
                                                  "RIF",
                                                  "RXT",
                                                  "SIS",
                                                  "TEC",
                                                  "TCY",
                                                  "TIC",
                                                  "TGC",
                                                  "TOB",
                                                  "TMP",
                                                  "SXT",
                                                  "VAN"),
                            hard_dependencies = NULL,
                            verbose = verbose,
                            ...)

  AMC <- cols_ab['AMC']
  AMK <- cols_ab['AMK']
  AMP <- cols_ab['AMP']
  AMX <- cols_ab['AMX']
  ATM <- cols_ab['ATM']
  AZL <- cols_ab['AZL']
  AZM <- cols_ab['AZM']
  CAZ <- cols_ab['CAZ']
  CED <- cols_ab['CED']
  CHL <- cols_ab['CHL']
  CIP <- cols_ab['CIP']
  CLI <- cols_ab['CLI']
  CLR <- cols_ab['CLR']
  COL <- cols_ab['COL']
  CRO <- cols_ab['CRO']
  CTX <- cols_ab['CTX']
  CXM <- cols_ab['CXM']
  CZO <- cols_ab['CZO']
  DAP <- cols_ab['DAP']
  DOX <- cols_ab['DOX']
  ERY <- cols_ab['ERY']
  ETP <- cols_ab['ETP']
  FEP <- cols_ab['FEP']
  FLC <- cols_ab['FLC']
  FOS <- cols_ab['FOS']
  FOX <- cols_ab['FOX']
  FUS <- cols_ab['FUS']
  GEN <- cols_ab['GEN']
  IPM <- cols_ab['IPM']
  KAN <- cols_ab['KAN']
  LIN <- cols_ab['LIN']
  LNZ <- cols_ab['LNZ']
  LVX <- cols_ab['LVX']
  MEM <- cols_ab['MEM']
  MEZ <- cols_ab['MEZ']
  MFX <- cols_ab['MFX']
  MNO <- cols_ab['MNO']
  NAL <- cols_ab['NAL']
  NEO <- cols_ab['NEO']
  NET <- cols_ab['NET']
  NIT <- cols_ab['NIT']
  NOR <- cols_ab['NOR']
  NOV <- cols_ab['NOV']
  OFX <- cols_ab['OFX']
  OXA <- cols_ab['OXA']
  PEN <- cols_ab['PEN']
  PIP <- cols_ab['PIP']
  PLB <- cols_ab['PLB']
  PRI <- cols_ab['PRI']
  QDA <- cols_ab['QDA']
  RID <- cols_ab['RID']
  RIF <- cols_ab['RIF']
  RXT <- cols_ab['RXT']
  SIS <- cols_ab['SIS']
  SXT <- cols_ab['SXT']
  TCY <- cols_ab['TCY']
  TEC <- cols_ab['TEC']
  TGC <- cols_ab['TGC']
  TIC <- cols_ab['TIC']
  TMP <- cols_ab['TMP']
  TOB <- cols_ab['TOB']
  TZP <- cols_ab['TZP']
  VAN <- cols_ab['VAN']

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
      before_df <- x_original

      tryCatch(
        # insert into original table
        x_original[rows, cols] <<- to,
        warning = function(w) {
          if (w$message %like% 'invalid factor level') {
            x_original <<- x_original %>% mutate_at(vars(cols), ~factor(x = as.character(.), levels = c(to, levels(.))))
            x <<- x %>% mutate_at(vars(cols), ~factor(x = as.character(.), levels = c(to, levels(.))))
            x_original[rows, cols] <<- to
            warning('Value "', to, '" added to the factor levels of column(s) `', paste(cols, collapse = '`, `'), '` because this value was not an existing factor level.\nA better way is to use as.rsi() on beforehand on antibiotic columns to guarantee the right structure.', call. = FALSE)
            txt_warning()
            warned <<- FALSE
          } else {
            warning(w$message, call. = FALSE)
            txt_warning()
            cat("\n") # txt_warning() does not append a "\n" on itself
          }
        },
        error = function(e) {
          txt_error()
          stop(paste0("In row(s) ", paste(rows[1:min(length(rows), 10)], collapse = ","), 
                      ifelse(length(rows) > 10, "...", ""),
                      ' while writing value "', to, 
                      '" to column(s) `', paste(cols, collapse = "`, `"), 
                      "` (data class: ", paste(class(x_original), collapse = "/"), "):\n", e$message),
               call. = FALSE)
        }
      )
      
      tryCatch(
        x[rows, cols] <<- x_original[rows, cols],
        error = function(e) {
          stop(paste0("Error in row(s) ", paste(rows[1:min(length(rows), 10)], collapse = ","),
                      '... while writing value "', to, 
                      '" to column(s) `', paste(cols, collapse = "`, `"), 
                      "` (data class:", paste(class(x), collapse = "/"), "):\n", e$message), call. = FALSE)
        }
      )
      
      # before_df might not be a data.frame, but a tibble or data.table instead
      old <- as.data.frame(before_df, stringsAsFactors = FALSE)[rows,]
      no_of_changes_this_run <- 0
      for (i in 1:length(cols)) {
        verbose_new <- data.frame(row = rows,
                                  col = cols[i],
                                  mo_fullname = x[rows, "fullname"],
                                  old = as.character(old[, cols[i]]),
                                  new = as.character(x[rows, cols[i]]),
                                  rule = strip_style(rule[1]),
                                  rule_group = strip_style(rule[2]),
                                  rule_name = strip_style(rule[3]),
                                  stringsAsFactors = FALSE)
        colnames(verbose_new) <- c("row", "col", "mo_fullname", "old", "new", "rule", "rule_group", "rule_name")
        verbose_new <- verbose_new %>% filter(old != new | is.na(old))
        # save changes to data set 'verbose_info'
        verbose_info <<- rbind(verbose_info, verbose_new)
        no_of_changes_this_run <- no_of_changes_this_run + nrow(verbose_new)
      }
      # after the applied changes: return number of (new) changes
      return(no_of_changes_this_run)
    }
    # no changes were applied: return number of (new) changes: none.
    return(0)
  }

  # save original table
  x_original <- x

  # join to microorganisms data set
  suppressWarnings(
    x <- x %>%
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
  if (!ab_missing(AMP) & !ab_missing(AMX)) {
    if (verbose == TRUE) {
      cat("\n VERBOSE: transforming",
          length(which(x[, AMX] == "S" & !x[, AMP] %in% c("S", "I", "R"))),
          "empty ampicillin fields to 'S' based on amoxicillin. ")
      cat("\n VERBOSE: transforming",
          length(which(x[, AMX] == "I" & !x[, AMP] %in% c("S", "I", "R"))),
          "empty ampicillin fields to 'I' based on amoxicillin. ")
      cat("\n VERBOSE: transforming",
          length(which(x[, AMX] == "R" & !x[, AMP] %in% c("S", "I", "R"))),
          "empty ampicillin fields to 'R' based on amoxicillin. \n")
    }
    x[which(x[, AMX] == "S" & !x[, AMP] %in% c("S", "I", "R")), AMP] <- "S"
    x[which(x[, AMX] == "I" & !x[, AMP] %in% c("S", "I", "R")), AMP] <- "I"
    x[which(x[, AMX] == "R" & !x[, AMP] %in% c("S", "I", "R")), AMP] <- "R"
  } else if (ab_missing(AMP) & !ab_missing(AMX)) {
    # ampicillin column is missing, but amoxicillin is available
    message(blue(paste0("NOTE: Using column `", bold(AMX), "` as input for ampicillin (J01CA01) since many EUCAST rules depend on it.")))
    AMP <- AMX
  }

  # antibiotic classes
  aminoglycosides <- c(TOB, GEN, KAN, NEO, NET, SIS)
  tetracyclines <- c(DOX, MNO, TCY) # since EUCAST v3.1 tigecycline (TGC) is set apart
  polymyxins <- c(PLB, COL)
  macrolides <- c(ERY, AZM, RXT, CLR) # since EUCAST v3.1 clinda is set apart
  glycopeptides <- c(VAN, TEC)
  streptogramins <- c(QDA, PRI) # should officially also be quinupristin/dalfopristin
  aminopenicillins <- c(AMP, AMX)
  cephalosporins <- c(FEP, CTX, FOX, CED, CAZ, CRO, CXM, CZO)
  cephalosporins_without_CAZ <- cephalosporins[cephalosporins != ifelse(is.null(CAZ), "", CAZ)]
  carbapenems <- c(ETP, IPM, MEM)
  ureidopenicillins <- c(PIP, TZP, AZL, MEZ)
  all_betalactams <- c(aminopenicillins, cephalosporins, carbapenems, ureidopenicillins, AMC, OXA, FLC, PEN)
  fluoroquinolones <- c(OFX, CIP, NOR, LVX, MFX)

  # Help function to get available antibiotic column names ------------------
  get_antibiotic_columns <- function(x, df) {
    x <- trimws(unlist(strsplit(x, ",", fixed = TRUE)))
    y <- character(0)
    for (i in 1:length(x)) {
      if (is.function(get(x[i]))) {
        stop("Column ", x[i], " is also a function. Please create an issue on github.com/msberends/AMR/issues.")
      }
      y <- c(y, tryCatch(get(x[i]), error = function(e) ""))
    }
    y[y != "" & y %in% colnames(df)]
  }
  get_antibiotic_names <- function(x) {
    x %>%
      strsplit(",") %>%
      unlist() %>%
      trimws() %>%
      sapply(function(x) if(x %in% AMR::antibiotics$ab) ab_name(x, language = NULL, tolower = TRUE) else x) %>%
      sort() %>%
      paste(collapse = ", ")
  }

  eucast_rules_df <- eucast_rules_file # internal data file
  no_of_changes <- 0
  for (i in 1:nrow(eucast_rules_df)) {

    rule_previous <- eucast_rules_df[max(1, i - 1), "reference.rule"]
    rule_current <- eucast_rules_df[i, "reference.rule"]
    rule_next <- eucast_rules_df[min(nrow(eucast_rules_df), i + 1), "reference.rule"]
    rule_group_previous <- eucast_rules_df[max(1, i - 1), "reference.rule_group"]
    rule_group_current <- eucast_rules_df[i, "reference.rule_group"]
    rule_group_next <- eucast_rules_df[min(nrow(eucast_rules_df), i + 1), "reference.rule_group"]
    if (is.na(eucast_rules_df[i, 4])) {
      rule_text <- paste0("always report as '", eucast_rules_df[i, 7], "': ", get_antibiotic_names(eucast_rules_df[i, 6]))
    } else {
      rule_text <- paste0("report as '", eucast_rules_df[i, 7], "' when ",
                          get_antibiotic_names(eucast_rules_df[i, 4]), " is '", eucast_rules_df[i, 5], "': ",
                          get_antibiotic_names(eucast_rules_df[i, 6]))
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
      rows <- tryCatch(which(x[, col_mo_property] %like% mo_value),
                       error = function(e) integer(0))
    } else {
      source_antibiotics <- get_antibiotic_columns(source_antibiotics, x)
      if (length(source_value) == 1 & length(source_antibiotics) > 1) {
        source_value <- rep(source_value, length(source_antibiotics))
      }
      if (length(source_antibiotics) == 0) {
        rows <- integer(0)
      } else if (length(source_antibiotics) == 1) {
        rows <-  tryCatch(which(x[, col_mo_property] %like% mo_value
                                & x[, source_antibiotics[1L]] == source_value[1L]),
                          error = function(e) integer(0))
      } else if (length(source_antibiotics) == 2) {
        rows <-  tryCatch(which(x[, col_mo_property] %like% mo_value
                                & x[, source_antibiotics[1L]] == source_value[1L]
                                & x[, source_antibiotics[2L]] == source_value[2L]),
                          error = function(e) integer(0))
      } else if (length(source_antibiotics) == 3) {
        rows <-  tryCatch(which(x[, col_mo_property] %like% mo_value
                                & x[, source_antibiotics[1L]] == source_value[1L]
                                & x[, source_antibiotics[2L]] == source_value[2L]
                                & x[, source_antibiotics[3L]] == source_value[3L]),
                          error = function(e) integer(0))
      } else {
        stop("only 3 antibiotics supported for source_antibiotics ", call. = FALSE)
      }
    }

    cols <- get_antibiotic_columns(target_antibiotics, x)

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

    cat(paste0("\n", silver(strrep("-", options()$width - 1)), "\n"))
    cat(bold(paste('EUCAST rules', paste0(wouldve, 'affected'),
                   formatnr(n_distinct(verbose_info$row)),
                   'out of', formatnr(nrow(x_original)),
                   'rows, making a total of', formatnr(nrow(verbose_info)), 'edits\n')))

    n_added <- verbose_info %>% filter(is.na(old)) %>% nrow()
    n_changed <- verbose_info %>% filter(!is.na(old)) %>% nrow()

    # print added values ----
    if (n_added == 0) {
      colour <- cat # is function
    } else {
      colour <- blue # is function
    }
    cat(colour(paste0("=> ", wouldve, "added ",
                      bold(formatnr(verbose_info %>%
                                      filter(is.na(old)) %>%
                                      nrow()), "test results"),
                      "\n")))
    if (n_added > 0) {
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
    if (n_changed == 0) {
      colour <- cat # is function
    } else {
      colour <- blue # is function
    }
    if (n_added + n_changed > 0) {
      cat("\n")
    }
    cat(colour(paste0("=> ", wouldve, "changed ",
                      bold(formatnr(verbose_info %>%
                                      filter(!is.na(old)) %>%
                                      nrow()), "test results"),
                      "\n")))
    if (n_changed > 0) {
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
      cat(paste("\nUse", bold("eucast_rules(..., verbose = TRUE)"), "(on your original data) to get a data.frame with all specified edits instead.\n\n"))
    } else if (verbose == TRUE) {
      cat(paste0("\nUsed 'Verbose mode' (", bold("verbose = TRUE"), "), which returns a data.frame with all specified edits.\nUse ", bold("verbose = FALSE"), " to apply the rules on your data.\n\n"))
    }
  }

  # Return data set ---------------------------------------------------------
  if (verbose == TRUE) {
    verbose_info
  } else {
    x_original
  }
}

