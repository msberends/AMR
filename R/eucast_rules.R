# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Analysis                              #
#                                                                      #
# SOURCE                                                               #
# https://github.com/msberends/AMR                                     #
#                                                                      #
# LICENCE                                                              #
# (c) 2018-2020 Berends MS, Luz CF et al.                              #
#                                                                      #
# This R package is free software; you can freely use and distribute   #
# it for both personal and commercial purposes under the terms of the  #
# GNU General Public License version 2.0 (GNU GPL-2), as published by  #
# the Free Software Foundation.                                        #
#                                                                      #
# We created this package for both routine data analysis and academic  #
# research and it was publicly released in the hope that it will be    #
# useful, but it comes WITHOUT ANY WARRANTY OR LIABILITY.              #
# Visit our website for more info: https://msberends.github.io/AMR.    #
# ==================================================================== #

# global variables
EUCAST_VERSION_BREAKPOINTS <- "10.0, 2020"
EUCAST_VERSION_EXPERT_RULES <- "3.1, 2016"

#' Apply EUCAST rules
#' 
#' @description
#' Apply susceptibility rules as defined by the European Committee on Antimicrobial Susceptibility Testing (EUCAST, <http://eucast.org>), see *Source*. This includes (1) expert rules and intrinsic resistance and (2) inferred resistance as defined in their breakpoint tables. 
#' 
#' To improve the interpretation of the antibiogram before EUCAST rules are applied, some non-EUCAST rules are applied at default, see Details.
#' @inheritSection lifecycle Maturing lifecycle
#' @param x data with antibiotic columns, like e.g. `AMX` and `AMC`
#' @param info print progress
#' @param rules a character vector that specifies which rules should be applied. Must be one or more of `"breakpoints"`, `"expert"`, `"other"`, `"all"`, and defaults to `c("breakpoints", "expert")`. The default value can be set to another value using e.g. `options(AMR.eucast_rules = "all")`.
#' @param verbose a logical to turn Verbose mode on and off (default is off). In Verbose mode, the function does not apply rules to the data, but instead returns a data set in logbook form with extensive info about which rows and columns would be effected and in which way.
#' @param ... column name of an antibiotic, please see section *Antibiotics* below
#' @inheritParams first_isolate
#' @details
#' **Note:** This function does not translate MIC values to RSI values. Use [as.rsi()] for that. \cr
#' **Note:** When ampicillin (AMP, J01CA01) is not available but amoxicillin (AMX, J01CA04) is, the latter will be used for all rules where there is a dependency on ampicillin. These drugs are interchangeable when it comes to expression of antimicrobial resistance.
#'
#' Before further processing, some non-EUCAST rules can be applied to improve the efficacy of the EUCAST rules. These non-EUCAST rules, that are then applied to all isolates, are:
#' - Inherit amoxicillin (AMX) from ampicillin (AMP), where amoxicillin (AMX) is unavailable;
#' - Inherit ampicillin (AMP) from amoxicillin (AMX), where ampicillin (AMP) is unavailable;
#' - Set amoxicillin (AMX) = R where amoxicillin/clavulanic acid (AMC) = R;
#' - Set piperacillin (PIP) = R where piperacillin/tazobactam (TZP) = R;
#' - Set trimethoprim (TMP) = R where trimethoprim/sulfamethoxazole (SXT) = R;
#' - Set amoxicillin/clavulanic acid (AMC) = S where amoxicillin (AMX) = S;
#' - Set piperacillin/tazobactam (TZP) = S where piperacillin (PIP) = S;
#' - Set trimethoprim/sulfamethoxazole (SXT) = S where trimethoprim (TMP) = S.
#' 
#' These rules are not applied at default, since they are not approved by EUCAST. To use these rules, please use `eucast_rules(..., rules = "all")`, or set the default behaviour of the `[eucast_rules()]` function with `options(AMR.eucast_rules = "all")` (or any other valid input value(s) to the `rules` parameter).
#'
#' The file containing all EUCAST rules is located here: <https://github.com/msberends/AMR/blob/master/data-raw/eucast_rules.tsv>.
#'
#' @section Antibiotics:
#' To define antibiotics column names, leave as it is to determine it automatically with [guess_ab_col()] or input a text (case-insensitive), or use `NULL` to skip a column (e.g. `TIC = NULL` to skip ticarcillin). Manually defined but non-existing columns will be skipped with a warning.
#'
#' The following antibiotics are used for the functions [eucast_rules()] and [mdro()]. These are shown below in the format '**antimicrobial ID**: name ([ATC code](https://www.whocc.no/atc/structure_and_principles/))', sorted by name:
#'
#'  **AMK**: amikacin ([J01GB06](https://www.whocc.no/atc_ddd_index/?code=J01GB06)),
#'  **AMX**: amoxicillin ([J01CA04](https://www.whocc.no/atc_ddd_index/?code=J01CA04)),
#'  **AMC**: amoxicillin/clavulanic acid ([J01CR02](https://www.whocc.no/atc_ddd_index/?code=J01CR02)),
#'  **AMP**: ampicillin ([J01CA01](https://www.whocc.no/atc_ddd_index/?code=J01CA01)),
#'  **SAM**: ampicillin/sulbactam ([J01CR01](https://www.whocc.no/atc_ddd_index/?code=J01CR01)),
#'  **AZM**: azithromycin ([J01FA10](https://www.whocc.no/atc_ddd_index/?code=J01FA10)),
#'  **AZL**: azlocillin ([J01CA09](https://www.whocc.no/atc_ddd_index/?code=J01CA09)),
#'  **ATM**: aztreonam ([J01DF01](https://www.whocc.no/atc_ddd_index/?code=J01DF01)),
#'  **CAP**: capreomycin ([J04AB30](https://www.whocc.no/atc_ddd_index/?code=J04AB30)),
#'  **RID**: cefaloridine ([J01DB02](https://www.whocc.no/atc_ddd_index/?code=J01DB02)),
#'  **CZO**: cefazolin ([J01DB04](https://www.whocc.no/atc_ddd_index/?code=J01DB04)),
#'  **FEP**: cefepime ([J01DE01](https://www.whocc.no/atc_ddd_index/?code=J01DE01)),
#'  **CTX**: cefotaxime ([J01DD01](https://www.whocc.no/atc_ddd_index/?code=J01DD01)),
#'  **CTT**: cefotetan ([J01DC05](https://www.whocc.no/atc_ddd_index/?code=J01DC05)),
#'  **FOX**: cefoxitin ([J01DC01](https://www.whocc.no/atc_ddd_index/?code=J01DC01)),
#'  **CPT**: ceftaroline ([J01DI02](https://www.whocc.no/atc_ddd_index/?code=J01DI02)),
#'  **CAZ**: ceftazidime ([J01DD02](https://www.whocc.no/atc_ddd_index/?code=J01DD02)),
#'  **CRO**: ceftriaxone ([J01DD04](https://www.whocc.no/atc_ddd_index/?code=J01DD04)),
#'  **CXM**: cefuroxime ([J01DC02](https://www.whocc.no/atc_ddd_index/?code=J01DC02)),
#'  **CED**: cephradine ([J01DB09](https://www.whocc.no/atc_ddd_index/?code=J01DB09)),
#'  **CHL**: chloramphenicol ([J01BA01](https://www.whocc.no/atc_ddd_index/?code=J01BA01)),
#'  **CIP**: ciprofloxacin ([J01MA02](https://www.whocc.no/atc_ddd_index/?code=J01MA02)),
#'  **CLR**: clarithromycin ([J01FA09](https://www.whocc.no/atc_ddd_index/?code=J01FA09)),
#'  **CLI**: clindamycin ([J01FF01](https://www.whocc.no/atc_ddd_index/?code=J01FF01)),
#'  **COL**: colistin ([J01XB01](https://www.whocc.no/atc_ddd_index/?code=J01XB01)),
#'  **DAP**: daptomycin ([J01XX09](https://www.whocc.no/atc_ddd_index/?code=J01XX09)),
#'  **DOR**: doripenem ([J01DH04](https://www.whocc.no/atc_ddd_index/?code=J01DH04)),
#'  **DOX**: doxycycline ([J01AA02](https://www.whocc.no/atc_ddd_index/?code=J01AA02)),
#'  **ETP**: ertapenem ([J01DH03](https://www.whocc.no/atc_ddd_index/?code=J01DH03)),
#'  **ERY**: erythromycin ([J01FA01](https://www.whocc.no/atc_ddd_index/?code=J01FA01)),
#'  **ETH**: ethambutol ([J04AK02](https://www.whocc.no/atc_ddd_index/?code=J04AK02)),
#'  **FLC**: flucloxacillin ([J01CF05](https://www.whocc.no/atc_ddd_index/?code=J01CF05)),
#'  **FOS**: fosfomycin ([J01XX01](https://www.whocc.no/atc_ddd_index/?code=J01XX01)),
#'  **FUS**: fusidic acid ([J01XC01](https://www.whocc.no/atc_ddd_index/?code=J01XC01)),
#'  **GAT**: gatifloxacin ([J01MA16](https://www.whocc.no/atc_ddd_index/?code=J01MA16)),
#'  **GEN**: gentamicin ([J01GB03](https://www.whocc.no/atc_ddd_index/?code=J01GB03)),
#'  **GEH**: gentamicin-high (no ATC code),
#'  **IPM**: imipenem ([J01DH51](https://www.whocc.no/atc_ddd_index/?code=J01DH51)),
#'  **INH**: isoniazid ([J04AC01](https://www.whocc.no/atc_ddd_index/?code=J04AC01)),
#'  **KAN**: kanamycin ([J01GB04](https://www.whocc.no/atc_ddd_index/?code=J01GB04)),
#'  **LVX**: levofloxacin ([J01MA12](https://www.whocc.no/atc_ddd_index/?code=J01MA12)),
#'  **LIN**: lincomycin ([J01FF02](https://www.whocc.no/atc_ddd_index/?code=J01FF02)),
#'  **LNZ**: linezolid ([J01XX08](https://www.whocc.no/atc_ddd_index/?code=J01XX08)),
#'  **MEM**: meropenem ([J01DH02](https://www.whocc.no/atc_ddd_index/?code=J01DH02)),
#'  **MTR**: metronidazole ([J01XD01](https://www.whocc.no/atc_ddd_index/?code=J01XD01)),
#'  **MEZ**: mezlocillin ([J01CA10](https://www.whocc.no/atc_ddd_index/?code=J01CA10)),
#'  **MNO**: minocycline ([J01AA08](https://www.whocc.no/atc_ddd_index/?code=J01AA08)),
#'  **MFX**: moxifloxacin ([J01MA14](https://www.whocc.no/atc_ddd_index/?code=J01MA14)),
#'  **NAL**: nalidixic acid ([J01MB02](https://www.whocc.no/atc_ddd_index/?code=J01MB02)),
#'  **NEO**: neomycin ([J01GB05](https://www.whocc.no/atc_ddd_index/?code=J01GB05)),
#'  **NET**: netilmicin ([J01GB07](https://www.whocc.no/atc_ddd_index/?code=J01GB07)),
#'  **NIT**: nitrofurantoin ([J01XE01](https://www.whocc.no/atc_ddd_index/?code=J01XE01)),
#'  **NOR**: norfloxacin ([J01MA06](https://www.whocc.no/atc_ddd_index/?code=J01MA06)),
#'  **NOV**: novobiocin ([QJ01XX95](https://www.whocc.no/atc_ddd_index/?code=QJ01XX95)),
#'  **OFX**: ofloxacin ([J01MA01](https://www.whocc.no/atc_ddd_index/?code=J01MA01)),
#'  **OXA**: oxacillin ([J01CF04](https://www.whocc.no/atc_ddd_index/?code=J01CF04)),
#'  **PEN**: penicillin G ([J01CE01](https://www.whocc.no/atc_ddd_index/?code=J01CE01)),
#'  **PIP**: piperacillin ([J01CA12](https://www.whocc.no/atc_ddd_index/?code=J01CA12)),
#'  **TZP**: piperacillin/tazobactam ([J01CR05](https://www.whocc.no/atc_ddd_index/?code=J01CR05)),
#'  **PLB**: polymyxin B ([J01XB02](https://www.whocc.no/atc_ddd_index/?code=J01XB02)),
#'  **PRI**: pristinamycin ([J01FG01](https://www.whocc.no/atc_ddd_index/?code=J01FG01)),
#'  **PZA**: pyrazinamide ([J04AK01](https://www.whocc.no/atc_ddd_index/?code=J04AK01)),
#'  **QDA**: quinupristin/dalfopristin ([J01FG02](https://www.whocc.no/atc_ddd_index/?code=J01FG02)),
#'  **RIB**: rifabutin ([J04AB04](https://www.whocc.no/atc_ddd_index/?code=J04AB04)),
#'  **RIF**: rifampicin ([J04AB02](https://www.whocc.no/atc_ddd_index/?code=J04AB02)),
#'  **RFP**: rifapentine ([J04AB05](https://www.whocc.no/atc_ddd_index/?code=J04AB05)),
#'  **RXT**: roxithromycin ([J01FA06](https://www.whocc.no/atc_ddd_index/?code=J01FA06)),
#'  **SIS**: sisomicin ([J01GB08](https://www.whocc.no/atc_ddd_index/?code=J01GB08)),
#'  **STH**: streptomycin-high (no ATC code),
#'  **TEC**: teicoplanin ([J01XA02](https://www.whocc.no/atc_ddd_index/?code=J01XA02)),
#'  **TLV**: telavancin ([J01XA03](https://www.whocc.no/atc_ddd_index/?code=J01XA03)),
#'  **TCY**: tetracycline ([J01AA07](https://www.whocc.no/atc_ddd_index/?code=J01AA07)),
#'  **TIC**: ticarcillin ([J01CA13](https://www.whocc.no/atc_ddd_index/?code=J01CA13)),
#'  **TCC**: ticarcillin/clavulanic acid ([J01CR03](https://www.whocc.no/atc_ddd_index/?code=J01CR03)),
#'  **TGC**: tigecycline ([J01AA12](https://www.whocc.no/atc_ddd_index/?code=J01AA12)),
#'  **TOB**: tobramycin ([J01GB01](https://www.whocc.no/atc_ddd_index/?code=J01GB01)),
#'  **TMP**: trimethoprim ([J01EA01](https://www.whocc.no/atc_ddd_index/?code=J01EA01)),
#'  **SXT**: trimethoprim/sulfamethoxazole ([J01EE01](https://www.whocc.no/atc_ddd_index/?code=J01EE01)),
#'  **VAN**: vancomycin ([J01XA01](https://www.whocc.no/atc_ddd_index/?code=J01XA01)).
#' @aliases EUCAST
#' @rdname eucast_rules
#' @export
#' @return The input of `x`, possibly with edited values of antibiotics. Or, if `verbose = TRUE`, a [`data.frame`] with all original and new values of the affected bug-drug combinations.
#' @source
#' - EUCAST Expert Rules. Version 2.0, 2012. \cr
#'   Leclercq et al. **EUCAST expert rules in antimicrobial susceptibility testing.** *Clin Microbiol Infect.* 2013;19(2):141-60. \cr
#'   <https://doi.org/10.1111/j.1469-0691.2011.03703.x>
#' - EUCAST Expert Rules, Intrinsic Resistance and Exceptional Phenotypes Tables. Version 3.1, 2016.  \cr
#'   <http://www.eucast.org/fileadmin/src/media/PDFs/EUCAST_files/Expert_Rules/Expert_rules_intrinsic_exceptional_V3.1.pdf>
#' - EUCAST Breakpoint tables for interpretation of MICs and zone diameters. Version 9.0, 2019. \cr
#'   <http://www.eucast.org/fileadmin/src/media/PDFs/EUCAST_files/Breakpoint_tables/v_9.0_Breakpoint_Tables.xlsx>
#' @inheritSection AMR Reference data publicly available
#' @inheritSection AMR Read more on our website!
#' @examples
#' \donttest{
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
#' }
eucast_rules <- function(x,
                         col_mo = NULL,
                         info = interactive(),
                         rules = getOption("AMR.eucast_rules", default = c("breakpoints", "expert")),
                         verbose = FALSE,
                         ...) {
  
  check_dataset_integrity()
  
  if (verbose == TRUE & interactive()) {
    txt <- paste0("WARNING: In Verbose mode, the eucast_rules() function does not apply rules to the data, but instead returns a data set in logbook form with extensive info about which rows and columns would be effected and in which way.",
                  "\n\nThis may overwrite your existing data if you use e.g.:",
                  "\ndata <- eucast_rules(data, verbose = TRUE)\n\nDo you want to continue?")
    if ("rstudioapi" %in% rownames(utils::installed.packages())) {
      showQuestion <- import_fn("showQuestion", "rstudioapi")
      q_continue <- showQuestion("Using verbose = TRUE with eucast_rules()", txt)
    } else {
      q_continue <- utils::menu(choices = c("OK", "Cancel"), graphics = FALSE, title = txt)
    }
    if (q_continue %in% c(FALSE, 2)) {
      message("Cancelled, returning original data")
      return(x)
    }
  }
  
  stop_ifnot(is.data.frame(x), "`x` must be a data frame")
  
  # try to find columns based on type
  # -- mo
  if (is.null(col_mo)) {
    col_mo <- search_type_in_df(x = x, type = "mo")
  }
  stop_if(is.null(col_mo), "`col_mo` must be set")
  
  stop_ifnot(all(rules %in% c("breakpoints", "expert", "other", "all")),
             '`rules` must be one or more of: "breakpoints", "expert", "other", "all".')
  
  decimal.mark <- getOption("OutDec")
  big.mark <- ifelse(decimal.mark != ",", ",", ".")
  formatnr <- function(x, big = big.mark, dec = decimal.mark) {
    trimws(format(x, big.mark = big, decimal.mark = dec))
  }
  
  warned <- FALSE
  warn_lacking_rsi_class <- FALSE
  
  txt_error <- function() {
    if (info == TRUE) cat("", font_red_bg(font_white(" ERROR ")), "\n\n") 
  }
  txt_warning <- function() {
    if (warned == FALSE) {
      if (info == TRUE) cat("", font_yellow_bg(font_black(" WARNING ")))
    }
    warned <<- TRUE 
  }
  txt_ok <- function(no_added, no_changed) {
    if (warned == FALSE) {
      if (no_added + no_changed == 0) {
        cat(font_subtle(" (no changes)\n"))
      } else {
        # opening
        cat(font_grey(" ("))
        # additions
        if (no_added > 0) {
          if (no_added == 1) {
            cat(font_green("1 value added"))
          } else {
            cat(font_green(formatnr(no_added), "values added"))
          }
        }
        # separator
        if (no_added > 0 & no_changed > 0) {
          cat(font_grey(", "))
        }
        # changes
        if (no_changed > 0) {
          if (no_changed == 1) {
            cat(font_blue("1 value changed"))
          } else {
            cat(font_blue(formatnr(no_changed), "values changed"))
          }
        } 
        # closing
        cat(font_grey(")\n"))
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
  
  AMC <- cols_ab["AMC"]
  AMK <- cols_ab["AMK"]
  AMP <- cols_ab["AMP"]
  AMX <- cols_ab["AMX"]
  ATM <- cols_ab["ATM"]
  AZL <- cols_ab["AZL"]
  AZM <- cols_ab["AZM"]
  CAZ <- cols_ab["CAZ"]
  CED <- cols_ab["CED"]
  CHL <- cols_ab["CHL"]
  CIP <- cols_ab["CIP"]
  CLI <- cols_ab["CLI"]
  CLR <- cols_ab["CLR"]
  COL <- cols_ab["COL"]
  CRO <- cols_ab["CRO"]
  CTX <- cols_ab["CTX"]
  CXM <- cols_ab["CXM"]
  CZO <- cols_ab["CZO"]
  DAP <- cols_ab["DAP"]
  DOX <- cols_ab["DOX"]
  ERY <- cols_ab["ERY"]
  ETP <- cols_ab["ETP"]
  FEP <- cols_ab["FEP"]
  FLC <- cols_ab["FLC"]
  FOS <- cols_ab["FOS"]
  FOX <- cols_ab["FOX"]
  FUS <- cols_ab["FUS"]
  GEN <- cols_ab["GEN"]
  IPM <- cols_ab["IPM"]
  KAN <- cols_ab["KAN"]
  LIN <- cols_ab["LIN"]
  LNZ <- cols_ab["LNZ"]
  LVX <- cols_ab["LVX"]
  MEM <- cols_ab["MEM"]
  MEZ <- cols_ab["MEZ"]
  MFX <- cols_ab["MFX"]
  MNO <- cols_ab["MNO"]
  NAL <- cols_ab["NAL"]
  NEO <- cols_ab["NEO"]
  NET <- cols_ab["NET"]
  NIT <- cols_ab["NIT"]
  NOR <- cols_ab["NOR"]
  NOV <- cols_ab["NOV"]
  OFX <- cols_ab["OFX"]
  OXA <- cols_ab["OXA"]
  PEN <- cols_ab["PEN"]
  PIP <- cols_ab["PIP"]
  PLB <- cols_ab["PLB"]
  PRI <- cols_ab["PRI"]
  QDA <- cols_ab["QDA"]
  RID <- cols_ab["RID"]
  RIF <- cols_ab["RIF"]
  RXT <- cols_ab["RXT"]
  SAM <- cols_ab["SAM"]
  SIS <- cols_ab["SIS"]
  SXT <- cols_ab["SXT"]
  TCY <- cols_ab["TCY"]
  TEC <- cols_ab["TEC"]
  TGC <- cols_ab["TGC"]
  TIC <- cols_ab["TIC"]
  TMP <- cols_ab["TMP"]
  TOB <- cols_ab["TOB"]
  TZP <- cols_ab["TZP"]
  VAN <- cols_ab["VAN"]
  
  ab_missing <- function(ab) {
    all(ab %in% c(NULL, NA))
  }
  
  verbose_info <- data.frame(row = integer(0),
                             col = character(0),
                             mo_fullname = character(0),
                             old = as.rsi(character(0)),
                             new = as.rsi(character(0)),
                             rule = character(0),
                             rule_group = character(0),
                             rule_name = character(0),
                             stringsAsFactors = FALSE)
  
  # helper function for editing the table
  edit_rsi <- function(to, rule, rows, cols) {
    cols <- unique(cols[!is.na(cols) & !is.null(cols)])
    if (length(rows) > 0 & length(cols) > 0) {
      before_df <- x_original
      if (any(!sapply(x[, cols, drop = FALSE], is.rsi), na.rm = TRUE)) {
        warn_lacking_rsi_class <<- TRUE
      }
      tryCatch(
        # insert into original table
        x_original[rows, cols] <<- to,
        warning = function(w) {
          if (w$message %like% "invalid factor level") {
            xyz <- sapply(cols, function(col) {
              x_original[, col] <<- factor(x = as.character(pull(x_original, col)), levels = c(to, levels(pull(x_original, col))))
              x[, col] <<- factor(x = as.character(pull(x, col)), levels = c(to, levels(pull(x, col))))
              invisible()
            })
            x_original[rows, cols] <<- to
            warning('Value "', to, '" added to the factor levels of column(s) `', paste(cols, collapse = "`, `"), "` because this value was not an existing factor level.\nA better way is to use as.rsi() on beforehand on antimicrobial columns to guarantee the right structure.", call. = FALSE)
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
                      " while writing value '", to, 
                      "' to column(s) `", paste(cols, collapse = "`, `"), 
                      "`:\n", e$message),
               call. = FALSE)
        }
      )
      
      tryCatch(
        x[rows, cols] <<- x_original[rows, cols],
        error = function(e) {
          stop(paste0("In row(s) ", paste(rows[1:min(length(rows), 10)], collapse = ","),
                      "... while writing value '", to, 
                      "' to column(s) `", paste(cols, collapse = "`, `"), 
                      "`:\n", e$message), call. = FALSE)
        }
      )
      
      # before_df might not be a data.frame, but a tibble or data.table instead
      old <- as.data.frame(before_df, stringsAsFactors = FALSE)[rows, ]
      track_changes <- list(added = 0,
                            changed = 0)
      for (i in seq_len(length(cols))) {
        verbose_new <- data.frame(row = rows,
                                  col = cols[i],
                                  mo_fullname = x[rows, "fullname"],
                                  old = as.rsi(as.character(old[, cols[i]]), warn = FALSE),
                                  new = as.rsi(as.character(x[rows, cols[i]])),
                                  rule = font_stripstyle(rule[1]),
                                  rule_group = font_stripstyle(rule[2]),
                                  rule_name = font_stripstyle(rule[3]),
                                  stringsAsFactors = FALSE)
        colnames(verbose_new) <- c("row", "col", "mo_fullname", "old", "new", "rule", "rule_group", "rule_name")
        verbose_new <- verbose_new %>% filter(old != new | is.na(old))
        # save changes to data set 'verbose_info'
        verbose_info <<- rbind(verbose_info, verbose_new)
        # count adds and changes
        track_changes$added <- track_changes$added + verbose_new %>% filter(is.na(old)) %>% nrow()
        track_changes$changed <- track_changes$changed + verbose_new %>% filter(!is.na(old)) %>% nrow()
      }
      # after the applied changes: return list with counts of added and changed
      return(track_changes)
    }
    # no changes were applied: return number of (new) changes: none.
    return(list(added = 0,
                changed = 0))
  }
  
  # save original table
  x_original <- x
  x_original_attr <- attributes(x)
  x_original <- as.data.frame(x_original, stringsAsFactors = FALSE) # no tibbles, data.tables, etc.
  
  # join to microorganisms data set
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  x[, col_mo] <- as.mo(x[, col_mo, drop = TRUE])
  x <- x %>%
    left_join_microorganisms(by = col_mo, suffix = c("_oldcols", ""))
  x$gramstain <- mo_gramstain(x[, col_mo, drop = TRUE], language = NULL)
  x$genus_species <- paste(x$genus, x$species)
  
  if (ab_missing(AMP) & !ab_missing(AMX)) {
    # ampicillin column is missing, but amoxicillin is available
    message(font_blue(paste0("NOTE: Using column `", font_bold(AMX), "` as input for ampicillin (J01CA01) since many EUCAST rules depend on it.")))
    AMP <- AMX
  }
  
  # nolint start
  # antibiotic classes
  aminoglycosides <- c(TOB, GEN, KAN, NEO, NET, SIS)
  tetracyclines <- c(DOX, MNO, TCY) # since EUCAST v3.1 tigecycline (TGC) is set apart
  polymyxins <- c(PLB, COL)
  macrolides <- c(ERY, AZM, RXT, CLR) # since EUCAST v3.1 clinda is set apart
  glycopeptides <- c(VAN, TEC)
  streptogramins <- c(QDA, PRI) # should officially also be quinupristin/dalfopristin
  aminopenicillins <- c(AMP, AMX)
  cephalosporins <- c(FEP, CTX, FOX, CED, CAZ, CRO, CXM, CZO)
  cephalosporins_except_CAZ <- cephalosporins[cephalosporins != ifelse(is.null(CAZ), "", CAZ)]
  carbapenems <- c(ETP, IPM, MEM)
  ureidopenicillins <- c(PIP, TZP, AZL, MEZ)
  all_betalactams <- c(aminopenicillins, cephalosporins, carbapenems, ureidopenicillins, AMC, OXA, FLC, PEN)
  fluoroquinolones <- c(OFX, CIP, NOR, LVX, MFX)
  # nolint end
  
  # Help function to get available antibiotic column names ------------------
  get_antibiotic_columns <- function(x, df) {
    x <- trimws(unlist(strsplit(x, ",", fixed = TRUE)))
    y <- character(0)
    for (i in seq_len(length(x))) {
      if (is.function(get(x[i]))) {
        stop("Column ", x[i], " is also a function. Please create an issue on github.com/msberends/AMR/issues.")
      }
      y <- c(y, tryCatch(get(x[i]), error = function(e) ""))
    }
    y[y != "" & y %in% colnames(df)]
  }
  get_antibiotic_names <- function(x) {
    x <- x %>%
      strsplit(",") %>%
      unlist() %>%
      trimws() %>%
      sapply(function(x) if (x %in% antibiotics$ab) ab_name(x, language = NULL, tolower = TRUE) else x) %>%
      sort() %>%
      paste(collapse = ", ")
    x <- gsub("_", " ", x, fixed = TRUE)
    x <- gsub("except CAZ", paste("except", ab_name("CAZ", language = NULL, tolower = TRUE)), x, fixed = TRUE)
    x
  }
  format_antibiotic_names <- function(ab_names, ab_results) {
    ab_names <- trimws(unlist(strsplit(ab_names, ",")))
    ab_results <- trimws(unlist(strsplit(ab_results, ",")))
    if (length(ab_results) == 1) {
      if (length(ab_names) == 1) {
        # like FOX S
        x <- paste(ab_names, "is")
      } else if (length(ab_names) == 2) {
        # like PEN,FOX S
        x <- paste(paste0(ab_names, collapse = " and "), "are both")
      } else {
        # like PEN,FOX,GEN S (although dependency on > 2 ABx does not exist at the moment)
        x <- paste(paste0(ab_names, collapse = " and "), "are all")
      }
      return(paste0(x, " '", ab_results, "'"))
    } else {
      if (length(ab_names) == 2) {
        # like PEN,FOX S,R
        paste0(ab_names[1], " is '", ab_results[1], "' and ", 
               ab_names[2], " is '", ab_results[2], "'")
      } else {
        # like PEN,FOX,GEN S,R,R (although dependency on > 2 ABx does not exist at the moment)
        paste0(ab_names[1], " is '", ab_results[1], "' and ", 
               ab_names[2], " is '", ab_results[2], "' and ", 
               ab_names[3], " is '", ab_results[3], "'")
      }
    }
  }
  
  as.rsi_no_warning <- function(x) suppressWarnings(as.rsi(x))
  no_added <- 0
  no_changed <- 0
  
  # Other rules: enzyme inhibitors ------------------------------------------
  if (any(c("all", "other") %in% rules)) {
    if (info == TRUE) {
      cat(font_bold(paste0("\nRules by this AMR package (",
                           font_red(paste0("v", utils::packageVersion("AMR"), ", ", 
                                           format(utils::packageDate("AMR"), "%Y"))), "), see ?eucast_rules\n")))
    }
    
    ab_enzyme <- subset(antibiotics, name %like% "/")[, c("ab", "name")]
    ab_enzyme$base_name <- gsub("^([a-zA-Z0-9]+).*", "\\1", ab_enzyme$name)
    ab_enzyme$base_ab <- as.ab(ab_enzyme$base_name)
    for (i in seq_len(nrow(ab_enzyme))) {
      if (all(c(ab_enzyme[i, ]$ab, ab_enzyme[i, ]$base_ab) %in% names(cols_ab), na.rm = TRUE)) {
        ab_name_base <- ab_name(cols_ab[ab_enzyme[i, ]$base_ab], language = NULL, tolower = TRUE)
        ab_name_enzyme <- ab_name(cols_ab[ab_enzyme[i, ]$ab], language = NULL, tolower = TRUE)
        
        # Set base to R where base + enzyme inhibitor is R
        rule_current <- paste0("Set ", ab_name_base, " (", cols_ab[ab_enzyme[i, ]$base_ab], ") = R where ",
                               ab_name_enzyme, " (", cols_ab[ab_enzyme[i, ]$ab], ") = R")
        if (info == TRUE) {
          cat(rule_current)
        }
        run_changes <- edit_rsi(to = "R",
                                rule = c(rule_current, "Other rules", ""),
                                rows = which(as.rsi_no_warning(x[, cols_ab[ab_enzyme[i, ]$ab]]) == "R"),
                                cols = cols_ab[ab_enzyme[i, ]$base_ab])
        no_added <- no_added + run_changes$added
        no_changed <- no_changed + run_changes$changed
        # Print number of new changes
        if (info == TRUE) {
          # print only on last one of rules in this group
          txt_ok(no_added = no_added, no_changed = no_changed)
          # and reset counters
          no_added <- 0
          no_changed <- 0
        }
        
        # Set base + enzyme inhibitor to S where base is S
        rule_current <- paste0("Set ", ab_name_enzyme, " (", cols_ab[ab_enzyme[i, ]$ab], ") = S where ",
                               ab_name_base, " (", cols_ab[ab_enzyme[i, ]$base_ab], ") = S")
        if (info == TRUE) {
          cat(rule_current)
        }
        run_changes <- edit_rsi(to = "S",
                                rule = c(rule_current, "Other rules", ""),
                                rows = which(as.rsi_no_warning(x[, cols_ab[ab_enzyme[i, ]$base_ab]]) == "S"),
                                cols = cols_ab[ab_enzyme[i, ]$ab])
        no_added <- no_added + run_changes$added
        no_changed <- no_changed + run_changes$changed
        # Print number of new changes
        if (info == TRUE) {
          # print only on last one of rules in this group
          txt_ok(no_added = no_added, no_changed = no_changed)
          # and reset counters
          no_added <- 0
          no_changed <- 0
        }
      }
    }
    
  } else {
    if (info == TRUE) {
      cat(font_red("\nSkipping inheritance rules defined by this package, such as setting trimethoprim (TMP) = R where trimethoprim/sulfamethoxazole (SXT) = R.\nUse eucast_rules(..., rules = \"all\") to also apply those rules.\n"))
    }
  }
  
  # Official EUCAST rules ---------------------------------------------------
  eucast_notification_shown <- FALSE
  if (!is.null(list(...)$eucast_rules_df)) {
    # this allows: eucast_rules(x, eucast_rules_df = AMR:::eucast_rules_file %>% filter(is.na(have_these_values)))
    eucast_rules_df <- list(...)$eucast_rules_df
  } else {
    # otherwise internal data file, created in data-raw/internals.R
    eucast_rules_df <- eucast_rules_file
  }
  for (i in seq_len(nrow(eucast_rules_df))) {
    
    rule_previous <- eucast_rules_df[max(1, i - 1), "reference.rule"]
    rule_current <- eucast_rules_df[i, "reference.rule"]
    rule_next <- eucast_rules_df[min(nrow(eucast_rules_df), i + 1), "reference.rule"]
    rule_group_previous <- eucast_rules_df[max(1, i - 1), "reference.rule_group"]
    rule_group_current <- eucast_rules_df[i, "reference.rule_group"]
    if (is.na(eucast_rules_df[i, 4])) {
      rule_text <- paste0("always report as '", eucast_rules_df[i, 7], "': ", get_antibiotic_names(eucast_rules_df[i, 6]))
    } else {
      rule_text <- paste0("report as '", eucast_rules_df[i, 7], "' when ",
                          format_antibiotic_names(ab_names = get_antibiotic_names(eucast_rules_df[i, 4]),
                                                  ab_results = eucast_rules_df[i, 5]), ": ",
                          get_antibiotic_names(eucast_rules_df[i, 6]))
    }
    if (i == 1) {
      rule_previous <- ""
      rule_group_previous <- ""
    }
    if (i == nrow(eucast_rules_df)) {
      rule_next <- ""
    }
    
    # don't apply rules if user doesn't want to apply them
    if (rule_group_current %like% "breakpoint" & !any(c("all", "breakpoints") %in% rules)) {
      next
    }
    if (rule_group_current %like% "expert" & !any(c("all", "expert") %in% rules)) {
      next
    }
    
    if (info == TRUE & !rule_group_current %like% "other" & eucast_notification_shown == FALSE) {
      cat(paste0("\n", font_grey(strrep("-", options()$width - 1)),
                 "\nRules by the ", font_bold("European Committee on Antimicrobial Susceptibility Testing (EUCAST)"),
                 "\n", font_blue("http://eucast.org/"), "\n"))
      eucast_notification_shown <- TRUE
    }
    
    if (info == TRUE) {
      # Print rule (group) ------------------------------------------------------
      if (rule_group_current != rule_group_previous) {
        # is new rule group, one of Breakpoints, Expert Rules and Other
        cat(font_bold(
          ifelse(
            rule_group_current %like% "breakpoint",
            paste0("\nEUCAST Clinical Breakpoints (",
                   font_red(paste0("v", EUCAST_VERSION_BREAKPOINTS)), ")\n"),
            ifelse(
              rule_group_current %like% "expert",
              paste0("\nEUCAST Expert Rules, Intrinsic Resistance and Exceptional Phenotypes (", 
                     font_red(paste0("v", EUCAST_VERSION_EXPERT_RULES)), ")\n"),
              ""))))
      }
      # Print rule  -------------------------------------------------------------
      if (rule_current != rule_previous) {
        # is new rule within group, print its name
        if (rule_current %in% c(microorganisms$family,
                                microorganisms$fullname)) {
          cat(font_italic(rule_current))
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
    if (eucast_rules_df[i, 3] %like% "coagulase") {
      all_staph <- microorganisms[which(microorganisms$genus == "Staphylococcus"), ]
      all_staph$CNS_CPS <- suppressWarnings(mo_name(all_staph$mo, Becker = "all", language = NULL))
      if (eucast_rules_df[i, 3] %like% "coagulase") {
        eucast_rules_df[i, 3] <- paste0("^(", paste0(all_staph[which(all_staph$CNS_CPS %like% "negative"),
                                                               "fullname", 
                                                               drop = TRUE],
                                                     collapse = "|"),
                                        ")$")
      } else {
        eucast_rules_df[i, 3] <- paste0("^(", paste0(all_staph[which(all_staph$CNS_CPS %like% "positive"),
                                                               "fullname", 
                                                               drop = TRUE],
                                                     collapse = "|"),
                                        ")$")
      }
      like_is_one_of <- "like"
    }
    
    if (like_is_one_of == "is") {
      # so e.g. 'Enterococcus' will turn into '^Enterococcus$'
      mo_value <- paste0("^", eucast_rules_df[i, 3], "$")
    } else if (like_is_one_of == "one_of") {
      # so 'Clostridium, Actinomyces, ...' will turn into '^(Clostridium|Actinomyces|...)$'
      mo_value <- paste0("^(",
                         paste(trimws(unlist(strsplit(eucast_rules_df[i, 3], ",", fixed = TRUE))),
                               collapse = "|"),
                         ")$")
    } else if (like_is_one_of == "like") {
      mo_value <- eucast_rules_df[i, 3]
    } else {
      stop("invalid value for column 'like.is.one_of'", call. = FALSE)
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
                                & as.rsi_no_warning(x[, source_antibiotics[1L]]) == source_value[1L]),
                          error = function(e) integer(0))
      } else if (length(source_antibiotics) == 2) {
        rows <-  tryCatch(which(x[, col_mo_property] %like% mo_value
                                & as.rsi_no_warning(x[, source_antibiotics[1L]]) == source_value[1L]
                                & as.rsi_no_warning(x[, source_antibiotics[2L]]) == source_value[2L]),
                          error = function(e) integer(0))
      } else if (length(source_antibiotics) == 3) {
        rows <-  tryCatch(which(x[, col_mo_property] %like% mo_value
                                & as.rsi_no_warning(x[, source_antibiotics[1L]]) == source_value[1L]
                                & as.rsi_no_warning(x[, source_antibiotics[2L]]) == source_value[2L]
                                & as.rsi_no_warning(x[, source_antibiotics[3L]]) == source_value[3L]),
                          error = function(e) integer(0))
      } else {
        stop("only 3 antibiotics supported for source_antibiotics ", call. = FALSE)
      }
    }
    
    cols <- get_antibiotic_columns(target_antibiotics, x)
    
    # Apply rule on data ------------------------------------------------------
    # this will return the unique number of changes
    run_changes <- edit_rsi(to = target_value,
                            rule = c(rule_text, rule_group_current, rule_current),
                            rows = rows,
                            cols = cols)
    no_added <- no_added + run_changes$added
    no_changed <- no_changed + run_changes$changed
    
    # Print number of new changes ---------------------------------------------
    if (info == TRUE & rule_next != rule_current) {
      # print only on last one of rules in this group
      txt_ok(no_added = no_added, no_changed = no_changed)
      # and reset counters
      no_added <- 0
      no_changed <- 0
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
    
    cat(paste0("\n", font_grey(strrep("-", options()$width - 1)), "\n"))
    cat(font_bold(paste("The rules", paste0(wouldve, "affected"),
                        formatnr(n_distinct(verbose_info$row)),
                        "out of", formatnr(nrow(x_original)),
                        "rows, making a total of", formatnr(nrow(verbose_info)), "edits\n")))
    
    n_added <- verbose_info %>% filter(is.na(old)) %>% nrow()
    n_changed <- verbose_info %>% filter(!is.na(old)) %>% nrow()
    
    # print added values ----
    if (n_added == 0) {
      colour <- cat # is function
    } else {
      colour <- font_green # is function
    }
    cat(colour(paste0("=> ", wouldve, "added ",
                      font_bold(formatnr(verbose_info %>%
                                           filter(is.na(old)) %>%
                                           nrow()), "test results"),
                      "\n")))
    if (n_added > 0) {
      added_summary <- verbose_info %>%
        filter(is.na(old)) %>%
        group_by(new) %>%
        summarise(n = n())
      cat(paste("   -", 
                paste0(formatnr(added_summary$n), " test result", ifelse(added_summary$n > 1, "s", ""), 
                       " added as ", added_summary$new), collapse = "\n"))
    }
    
    # print changed values ----
    if (n_changed == 0) {
      colour <- cat # is function
    } else {
      colour <- font_blue # is function
    }
    if (n_added + n_changed > 0) {
      cat("\n")
    }
    cat(colour(paste0("=> ", wouldve, "changed ",
                      font_bold(formatnr(verbose_info %>%
                                           filter(!is.na(old)) %>%
                                           nrow()), "test results"),
                      "\n")))
    if (n_changed > 0) {
      changed_summary <- verbose_info %>%
        filter(!is.na(old)) %>%
        group_by(old, new) %>%
        summarise(n = n())
      cat(paste("   -", 
                paste0(formatnr(changed_summary$n), " test result", ifelse(changed_summary$n > 1, "s", ""), " changed from ", 
                       changed_summary$old, " to ", changed_summary$new), collapse = "\n"))
      cat("\n")
    }
    cat(paste0(font_grey(strrep("-", options()$width - 1)), "\n"))
    
    if (verbose == FALSE & nrow(verbose_info) > 0) {
      cat(paste("\nUse", font_bold("eucast_rules(..., verbose = TRUE)"), "(on your original data) to get a data.frame with all specified edits instead.\n\n"))
    } else if (verbose == TRUE) {
      cat(paste0("\nUsed 'Verbose mode' (", font_bold("verbose = TRUE"), "), which returns a data.frame with all specified edits.\nUse ", font_bold("verbose = FALSE"), " to apply the rules on your data.\n\n"))
    }
  }
  
  if (isTRUE(warn_lacking_rsi_class)) {
    warning("Not all columns with antimicrobial results are of class <rsi>.\n",
            "Transform eligible columns to class <rsi> on beforehand: your_data %>% mutate_if(is.rsi.eligible, as.rsi)",
            call. = FALSE)
  }
  
  # Return data set ---------------------------------------------------------
  if (verbose == TRUE) {
    rownames(verbose_info) <- NULL
    verbose_info
  } else {
    # reset original attributes
    attributes(x_original) <- x_original_attr
    x_original
  }
}
