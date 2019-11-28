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
#' Determine which isolates are multidrug-resistant organisms (MDRO) according to international and national guidelines.
#' @param guideline a specific guideline to follow. When left empty, the publication by Magiorakos *et al.* (2012, Clinical Microbiology and Infection) will be followed, please see *Details*.
#' @param info a logical to indicate whether progress should be printed to the console
#' @inheritParams eucast_rules
#' @param pct_required_classes minimal required percentage of antimicrobial classes that must be available per isolate, rounded down. For example, with the default guideline, 17 antimicrobial classes must be available for *S. aureus*. Setting this `pct_required_classes` argument to `0.5` (default) means that for every *S. aureus* isolate at least 8 different classes must be available. Any lower number of available classes will return `NA` for that isolate.
#' @param combine_SI a logical to indicate whether all values of S and I must be merged into one, so resistance is only considered when isolates are R, not I. As this is the default behaviour of the [mdro()] function, it follows the redefinition by EUCAST about the interpretion of I (increased exposure) in 2019, see section 'Interpretation of S, I and R' below. When using `combine_SI = FALSE`, resistance is considered when isolates are R or I.
#' @param verbose a logical to turn Verbose mode on and off (default is off). In Verbose mode, the function does not return the MDRO results, but instead returns a data set in logbook form with extensive info about which isolates would be MDRO-positive, or why they are not.
#' @inheritSection eucast_rules Antibiotics
#' @details 
#' For the `pct_required_classes` argument, values above 1 will be divided by 100. This is to support both fractions (`0.75` or `3/4`) and percentages (`75`).
#' 
#' Currently supported guidelines are (case-insensitive):
#' 
#' - `guideline = "CMI2012"`\cr
#'   Magiorakos AP, Srinivasan A *et al.* "Multidrug-resistant, extensively drug-resistant and pandrug-resistant bacteria: an international expert proposal for interim standard definitions for acquired resistance." Clinical Microbiology and Infection (2012) ([link](https://www.clinicalmicrobiologyandinfection.com/article/S1198-743X(14)61632-3/fulltext))
#' - `guideline = "EUCAST"`\cr
#'   The European international guideline - EUCAST Expert Rules Version 3.1 "Intrinsic Resistance and Exceptional Phenotypes Tables" ([link](http://www.eucast.org/fileadmin/src/media/PDFs/EUCAST_files/Expert_Rules/Expert_rules_intrinsic_exceptional_V3.1.pdf))
#' - `guideline = "TB"`\cr
#'   The international guideline for multi-drug resistant tuberculosis - World Health Organization "Companion handbook to the WHO guidelines for the programmatic management of drug-resistant tuberculosis" ([link](https://www.who.int/tb/publications/pmdt_companionhandbook/en/))
#' - `guideline = "MRGN"`\cr
#'   The German national guideline - Mueller et al. (2015) Antimicrobial Resistance and Infection Control 4:7. DOI: 10.1186/s13756-015-0047-6
#' - `guideline = "BRMO"`\cr
#'   The Dutch national guideline - Rijksinstituut voor Volksgezondheid en Milieu "WIP-richtlijn BRMO (Bijzonder Resistente Micro-Organismen) [ZKH]" ([link](https://www.rivm.nl/Documenten_en_publicaties/Professioneel_Praktisch/Richtlijnen/Infectieziekten/WIP_Richtlijnen/WIP_Richtlijnen/Ziekenhuizen/WIP_richtlijn_BRMO_Bijzonder_Resistente_Micro_Organismen_ZKH))
#'
#' Please suggest your own (country-specific) guidelines by letting us know: <https://gitlab.com/msberends/AMR/issues/new>.
#' 
#' **Note:** Every test that involves the Enterobacteriaceae family, will internally be performed using its newly named order Enterobacterales, since the Enterobacteriaceae family has been taxonomically reclassified by Adeolu *et al.* in 2016. Before that, Enterobacteriaceae was the only family under the Enterobacteriales (with an i) order. All species under the old Enterobacteriaceae family are still under the new Enterobacterales (without an i) order, but divided into multiple families. The way tests are performed now by this [mdro()] function makes sure that results from before 2016 and after 2016 are identical.
#' @inheritSection as.rsi Interpretation of S, I and R
#' @return
#' - CMI 2012 paper - function [mdr_cmi2012()] or [mdro()]:\cr
#'   Ordered [`factor`] with levels `Negative` < `Multi-drug-resistant (MDR)` < `Extensively drug-resistant (XDR)` < `Pandrug-resistant (PDR)`
#' - TB guideline - function [mdr_tb()] or [`mdro(..., guideline = "TB")`][mdro()]:\cr
#'   Ordered [`factor`] with levels `Negative` < `Mono-resistant` < `Poly-resistant` < `Multi-drug-resistant` < `Extensively drug-resistant`
#' - German guideline - function [mrgn()] or [`mdro(..., guideline = "MRGN")`][mdro()]:\cr
#'   Ordered [`factor`] with levels `Negative` < `3MRGN` < `4MRGN`
#' - Everything else:\cr
#'   Ordered [`factor`] with levels `Negative` < `Positive, unconfirmed` < `Positive`. The value `"Positive, unconfirmed"` means that, according to the guideline, it is not entirely sure if the isolate is multi-drug resistant and this should be confirmed with additional (e.g. molecular) tests
#' @rdname mdro
#' @aliases MDR XDR PDR BRMO 3MRGN 4MRGN
#' @importFrom dplyr %>% filter_at vars all_vars pull mutate_at
#' @importFrom crayon blue bold italic red
#' @importFrom cleaner percentage
#' @export
#' @inheritSection AMR Read more on our website!
#' @source
#' Please see *Details* for the list of publications used for this function.
#' @examples
#' \donttest{
#' library(dplyr)
#' 
#' example_isolates %>%
#'   mdro() %>%
#'   freq()
#'   
#' example_isolates %>%
#'   mutate(EUCAST = eucast_exceptional_phenotypes(.),
#'          BRMO = brmo(.),
#'          MRGN = mrgn(.))
#' }
mdro <- function(x,
                 guideline = "CMI2012",
                 col_mo = NULL,
                 info = TRUE,
                 pct_required_classes = 0.5,
                 combine_SI = TRUE,
                 verbose = FALSE,
                 ...) {
  
  if (verbose == TRUE & interactive()) {
    txt <- paste0("WARNING: In Verbose mode, the mdro() function does not return the MDRO results, but instead returns a data set in logbook form with extensive info about which isolates would be MDRO-positive, or why they are not.",
                  "\n\nThis may overwrite your existing data if you use e.g.:",
                  "\ndata <- mdro(data, verbose = TRUE)\n\nDo you want to continue?")
    if ("rstudioapi" %in% rownames(utils::installed.packages())) {
      q_continue <- rstudioapi::showQuestion("Using verbose = TRUE with mdro()", txt)
    } else {
      q_continue <- menu(choices = c("OK", "Cancel"), graphics = TRUE, title = txt)
    }
    if (q_continue %in% c(FALSE, 2)) {
      message("Cancelled, returning original data")
      return(x)
    }
  }
  
  if (!is.data.frame(x)) {
    stop("`x` must be a data frame.", call. = FALSE)
  }
  if (!is.numeric(pct_required_classes)) {
    stop("`pct_required_classes` must be numeric.", call. = FALSE)
  }
  if (pct_required_classes > 1) {
    # allow pct_required_classes = 75 -> pct_required_classes = 0.75
    pct_required_classes <- pct_required_classes / 100
  }

  if (!is.null(list(...)$country)) {
    warning("Using `country` is deprecated, use `guideline` instead. Please see ?mdro.", call. = FALSE)
    guideline <- list(...)$country
  }
  if (length(guideline) > 1) {
    stop("`guideline` must be a length one character string.", call. = FALSE)
  }
  
  if (is.null(guideline)) {
    # default to the paper by Magiorakos et al. (2012)
    guideline <- "cmi2012"
  }
  if (tolower(guideline) == "nl") {
    guideline <- "BRMO"
  }
  if (tolower(guideline) == "de") {
    guideline <- "MRGN"
  }
  if (!tolower(guideline) %in% c("brmo", "mrgn", "eucast", "tb", "cmi2012")) {
    stop("invalid guideline: ", guideline, call. = FALSE)
  }
  guideline <- list(code = tolower(guideline))
  
  # try to find columns based on type
  # -- mo
  if (is.null(col_mo)) {
    col_mo <- search_type_in_df(x = x, type = "mo")
  }
  if (is.null(col_mo) & guideline$code == "tb") {
    message(blue("NOTE: No column found as input for `col_mo`,",
                 bold("assuming all records contain", italic("Mycobacterium tuberculosis.\n"))))
    x$mo <- AMR::as.mo("Mycobacterium tuberculosis")
    col_mo <- "mo"
  }
  if (is.null(col_mo)) {
    stop("`col_mo` must be set.", call. = FALSE)
  }
  
  if (guideline$code == "cmi2012") {
    guideline$name <- "Multidrug-resistant, extensively drug-resistant and pandrug-resistant bacteria: an international expert proposal for interim standard definitions for acquired resistance."
    guideline$author <- "Magiorakos AP, Srinivasan A, Carey RB, ..., Vatopoulos A, Weber JT, Monnet DL"
    guideline$version <- "N/A"
    guideline$source <- "Clinical Microbiology and Infection 18:3, 2012. DOI: 10.1111/j.1469-0691.2011.03570.x"
    
  } else if (guideline$code == "eucast") {
    guideline$name <- "EUCAST Expert Rules, \"Intrinsic Resistance and Exceptional Phenotypes Tables\""
    guideline$author <- "EUCAST (European Committee on Antimicrobial Susceptibility Testing)"
    guideline$version <- "3.1"
    guideline$source <- "http://www.eucast.org/fileadmin/src/media/PDFs/EUCAST_files/Expert_Rules/Expert_rules_intrinsic_exceptional_V3.1.pdf"
    
  } else if (guideline$code == "tb") {
    guideline$name <- "Companion handbook to the WHO guidelines for the programmatic management of drug-resistant tuberculosis"
    guideline$author <- "WHO (World Health Organization)"
    guideline$version <- "WHO/HTM/TB/2014.11"
    guideline$source <- "https://www.who.int/tb/publications/pmdt_companionhandbook/en/"
    
    # support per country:
  } else if (guideline$code == "mrgn") {
    guideline$name <- "Cross-border comparison of the Dutch and German guidelines on multidrug-resistant Gram-negative microorganisms"
    guideline$author <- "M\u00fcller J, Voss A, K\u00f6ck R, ..., Kern WV, Wendt C, Friedrich AW"
    guideline$version <- "N/A"
    guideline$source <- "Antimicrobial Resistance and Infection Control 4:7, 2015. DOI: 10.1186/s13756-015-0047-6"
    
  } else if (guideline$code == "brmo") {
    guideline$name <- "WIP-Richtlijn Bijzonder Resistente Micro-organismen (BRMO)"
    guideline$author <- "RIVM (Rijksinstituut voor de Volksgezondheid)"
    guideline$version <- "Revision as of December 2017"
    guideline$source <- "https://www.rivm.nl/Documenten_en_publicaties/Professioneel_Praktisch/Richtlijnen/Infectieziekten/WIP_Richtlijnen/WIP_Richtlijnen/Ziekenhuizen/WIP_richtlijn_BRMO_Bijzonder_Resistente_Micro_Organismen_ZKH"
  } else {
    stop("This guideline is currently unsupported: ", guideline$code, call. = FALSE)
  }
  
  if (guideline$code == "cmi2012") {
    cols_ab <- get_column_abx(x = x,
                              soft_dependencies = c(
                                # table 1 (S aureus):
                                "GEN",
                                "RIF",
                                "CPT",
                                "OXA",
                                "CIP",
                                "MFX",
                                "SXT",
                                "FUS",
                                "VAN",
                                "TEC",
                                "TLV",
                                "TGC",
                                "CLI",
                                "DAP",
                                "ERY",
                                "LNZ",
                                "CHL",
                                "FOS",
                                "QDA",
                                "TCY",
                                "DOX",
                                "MNO",
                                # table 2 (Enterococcus)
                                "GEH",
                                "STH",
                                "IPM",
                                "MEM",
                                "DOR",
                                "CIP",
                                "LVX",
                                "MFX",
                                "VAN",
                                "TEC",
                                "TGC",
                                "DAP",
                                "LNZ",
                                "AMP",
                                "QDA",
                                "DOX",
                                "MNO",
                                # table 3 (Enterobacteriaceae)
                                "GEN",
                                "TOB",
                                "AMK",
                                "NET",
                                "CPT",
                                "TCC",
                                "TZP",
                                "ETP",
                                "IPM",
                                "MEM",
                                "DOR",
                                "CZO",
                                "CXM",
                                "CTX",
                                "CAZ",
                                "FEP",
                                "FOX",
                                "CTT",
                                "CIP",
                                "SXT",
                                "TGC",
                                "ATM",
                                "AMP",
                                "AMC",
                                "SAM",
                                "CHL",
                                "FOS",
                                "COL",
                                "TCY",
                                "DOX",
                                "MNO",
                                # table 4 (Pseudomonas)
                                "GEN",
                                "TOB",
                                "AMK",
                                "NET",
                                "IPM",
                                "MEM",
                                "DOR",
                                "CAZ",
                                "FEP",
                                "CIP",
                                "LVX",
                                "TCC",
                                "TZP",
                                "ATM",
                                "FOS",
                                "COL",
                                "PLB",
                                # table 5 (Acinetobacter)
                                "GEN",
                                "TOB",
                                "AMK",
                                "NET",
                                "IPM",
                                "MEM",
                                "DOR",
                                "CIP",
                                "LVX",
                                "TZP",
                                "TCC",
                                "CTX",
                                "CRO",
                                "CAZ",
                                "FEP",
                                "SXT",
                                "SAM",
                                "COL",
                                "PLB",
                                "TCY",
                                "DOX",
                                "MNO"),
                              verbose = verbose, ...)
  } else if (guideline$code == "tb") {
    cols_ab <- get_column_abx(x = x,
                              soft_dependencies = c("CAP",
                                                    "ETH",
                                                    "GAT",
                                                    "INH",
                                                    "PZA",
                                                    "RIF",
                                                    "RIB",
                                                    "RFP"),
                              verbose = verbose, ...)
  } else if (guideline$code == "mrgn") {
    cols_ab <- get_column_abx(x = x,
                              soft_dependencies = c("PIP",
                                                    "CTX",
                                                    "CAZ",
                                                    "IPM",
                                                    "MEM",
                                                    "CIP"),
                              verbose = verbose, ...)
  } else {
    cols_ab <- get_column_abx(x = x, verbose = verbose, ...)
  }
  
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
  CPT <- cols_ab["CPT"]
  CRO <- cols_ab["CRO"]
  CTT <- cols_ab["CTT"]
  CTX <- cols_ab["CTX"]
  CXM <- cols_ab["CXM"]
  CZO <- cols_ab["CZO"]
  DAP <- cols_ab["DAP"]
  DOR <- cols_ab["DOR"]
  DOX <- cols_ab["DOX"]
  ERY <- cols_ab["ERY"]
  ETP <- cols_ab["ETP"]
  FEP <- cols_ab["FEP"]
  FLC <- cols_ab["FLC"]
  FOS <- cols_ab["FOS"]
  FOX <- cols_ab["FOX"]
  FUS <- cols_ab["FUS"]
  GEH <- cols_ab["GEH"]
  GEN <- cols_ab["GEN"]
  IPM <- cols_ab["IPM"]
  KAN <- cols_ab["KAN"]
  LIN <- cols_ab["LIN"]
  LNZ <- cols_ab["LNZ"]
  LVX <- cols_ab["LVX"]
  MEM <- cols_ab["MEM"]
  MEZ <- cols_ab["MEZ"]
  MTR <- cols_ab["MTR"]
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
  STH <- cols_ab["STH"]
  SXT <- cols_ab["SXT"]
  TCC <- cols_ab["TCC"]
  TCY <- cols_ab["TCY"]
  TEC <- cols_ab["TEC"]
  TGC <- cols_ab["TGC"]
  TIC <- cols_ab["TIC"]
  TLV <- cols_ab["TLV"]
  TMP <- cols_ab["TMP"]
  TOB <- cols_ab["TOB"]
  TZP <- cols_ab["TZP"]
  VAN <- cols_ab["VAN"]
  # additional for TB
  CAP <- cols_ab["CAP"]
  ETH <- cols_ab["ETH"]
  GAT <- cols_ab["GAT"]
  INH <- cols_ab["INH"]
  PZA <- cols_ab["PZA"]
  RIF <- cols_ab["RIF"]
  RIB <- cols_ab["RIB"]
  RFP <- cols_ab["RFP"]
  abx_tb <- c(CAP, ETH, GAT, INH, PZA, RIF, RIB, RFP)
  abx_tb <- abx_tb[!is.na(abx_tb)]
  if (guideline$code == "tb" & length(abx_tb) == 0) {
    stop("No antimycobacterials found in data set.", call. = FALSE)
  }

  if (combine_SI == TRUE) {
    search_result <- "R"
  } else {
    search_result <- c("R", "I")
  }
  
  if (info == TRUE) {
    if (combine_SI == TRUE) {
      cat(red("\nOnly results with 'R' are considered as resistance. Use `combine_SI = FALSE` to also consider 'I' as resistance.\n"))
    } else {
      cat(red("\nResults with 'R' or 'I' are considered as resistance. Use `combine_SI = TRUE` to only consider 'R' as resistance.\n"))
    }
    cat("\nDetermining multidrug-resistant organisms (MDRO), according to:\n",
        bold("Guideline: "), italic(guideline$name), "\n",
        bold("Version:   "), guideline$version, "\n",
        bold("Author:    "),    guideline$author, "\n",
        bold("Source:    "), guideline$source, "\n",
        "\n", sep = "")
  }
  
  ab_missing <- function(ab) {
    isTRUE(ab %in% c(NULL, NA)) | length(ab) == 0
  }
  ab_NA <- function(x) {
    x[!is.na(x)]
  }
  
  verbose_df <- NULL
  
  # antibiotic classes
  aminoglycosides <- c(TOB, GEN)
  cephalosporins <- c(FEP, CTX, FOX, CED, CAZ, CRO, CXM, CZO)
  cephalosporins_3rd <- c(CTX, CRO, CAZ)
  carbapenems <- c(ETP, IPM, MEM)
  fluoroquinolones <- c(OFX, CIP, LVX, MFX)
  
  # helper function for editing the table
  trans_tbl <- function(to, rows, cols, any_all) {
    cols <- cols[!ab_missing(cols)]
    cols <- cols[!is.na(cols)]
    if (length(rows) > 0 & length(cols) > 0) {
      x <<- x %>% mutate_at(vars(cols), as.rsi)
      x[rows, "columns_nonsusceptible"] <<- sapply(rows, 
                                                   function(row, group_vct = cols) {
                                                     cols_nonsus <- sapply(x[row, group_vct, drop = FALSE], 
                                                                           function(y) y %in% search_result)
                                                     paste(sort(c(unlist(strsplit(x[row, "columns_nonsusceptible", drop = TRUE], ", ")),
                                                                  names(cols_nonsus)[cols_nonsus])), 
                                                           collapse = ", ")
                                                   })
      
      if (any_all == "any") {
        search_function <- dplyr::any_vars
      } else if (any_all == "all") {
        search_function <- dplyr::all_vars
      }
      row_filter <- x %>%
        filter_at(vars(cols), search_function(. %in% search_result)) %>%
        pull("row_number")
      rows <- rows[rows %in% row_filter]
      x[rows, "MDRO"] <<- to
      x[rows, "reason"] <<- paste0(any_all, " of the required antibiotics ", ifelse(any_all == "any", "is", "are"), " R")
    }
  }
  trans_tbl2 <- function(txt, rows, lst) {
    if (info == TRUE) {
      message(blue(txt, "..."), appendLF = FALSE)
    }
    if (length(rows) > 0) {
      # function specific for the CMI paper of 2012 (Magiorakos et al.)
      lst_vector <- unlist(lst)[!is.na(unlist(lst))]
      x <<- x %>% mutate_at(vars(lst_vector), as.rsi)
      x[rows, "classes_in_guideline"] <<- length(lst)
      x[rows, "classes_available"] <<- sapply(rows, 
                                              function(row, group_tbl = lst) {
                                                sum(sapply(group_tbl, function(group) any(unlist(x[row, group[!is.na(group)], drop = TRUE]) %in% c("S", "I", "R"))))
                                              })
      
      if (verbose == TRUE) {
        x[rows, "columns_nonsusceptible"] <<- sapply(rows, 
                                                      function(row, group_vct = lst_vector) {
                                                        cols_nonsus <- sapply(x[row, group_vct, drop = FALSE], function(y) y %in% search_result)
                                                        paste(sort(names(cols_nonsus)[cols_nonsus]), collapse = ", ")
                                                      })
      }
      x[rows, "classes_affected"] <<- sapply(rows, 
                                            function(row, group_tbl = lst) {
                                              sum(sapply(group_tbl, 
                                                         function(group) {
                                                           any(unlist(x[row, group[!is.na(group)], drop = TRUE]) %in% search_result, na.rm = TRUE)
                                                         }),
                                                  na.rm = TRUE) 
                                            })
      # for PDR; all agents are R (or I if combine_SI = FALSE)
      x[filter_at(x[rows, ],
                  vars(lst_vector),
                  all_vars(. %in% search_result))$row_number, "classes_affected"] <<- 999
    }
    
    if (info == TRUE) {
      message(blue(" OK"))
    }
  }
  
  x <- x %>%
    mutate_at(vars(col_mo), as.mo) %>%
    # join to microorganisms data set
    left_join_microorganisms(by = col_mo) %>%
    # add unavailable to where genus is available
    mutate(MDRO = ifelse(!is.na(genus), 1, NA_integer_),
           row_number = seq_len(nrow(.)),
           reason = paste0("not covered by ", toupper(guideline$code), " guideline"),
           columns_nonsusceptible = "") %>% 
    # transform to data.frame so subsetting is possible with x[y, z] (might not be the case with tibble/data.table/...)
    as.data.frame(stringsAsFactors = FALSE)
  
  if (guideline$code == "cmi2012") {
    # CMI, 2012 ---------------------------------------------------------------
    # Non-susceptible = R and I
    # (see header 'Approaches to Creating Definitions for MDR, XDR and PDR' in paper)
    
    # take amoxicillin if ampicillin is unavailable
    if (is.na(AMP) & !is.na(AMX)) {
      if (verbose == TRUE) {
        message(blue("NOTE: Filling ampicillin (AMP) results with amoxicillin (AMX) results"))
      }
      AMP <- AMX
    }
    # take ceftriaxone if cefotaxime is unavailable and vice versa
    if (is.na(CRO) & !is.na(CTX)) {
      if (verbose == TRUE) {
        message(blue("NOTE: Filling ceftriaxone (CRO) results with cefotaxime (CTX) results"))
      }
      CRO <- CTX
    }
    if (is.na(CTX) & !is.na(CRO)) {
      if (verbose == TRUE) {
        message(blue("NOTE: Filling cefotaxime (CTX) results with ceftriaxone (CRO) results"))
      }
      CTX <- CRO
    }
    
    # intrinsic resistant must not be considered for the determination of MDR,
    # so let's just remove them, meticulously following the paper
    x[which(x$genus == "Enterococcus" & x$species == "faecium"), ab_NA(IPM)] <- NA
    x[which(x$genus == "Enterococcus" & x$species == "faecalis"), ab_NA(QDA)] <- NA
    x[which((x$genus == "Providencia" & x$species == "rettgeri")
            | (x$genus == "Providencia" & x$species == "stuartii")), ab_NA(c(GEN, TOB, NET))] <- NA
    x[which(x$genus == "Escherichia" & x$species == "hermannii"), ab_NA(c(TCC, TZP))] <- NA
    x[which((x$genus == "Citrobacter" & x$species == "freundii")
            | (x$genus == "Enterobacter" & x$species == "aerogenes")
            | (x$genus == "Enterobacter" & x$species == "cloacae")
            | (x$genus == "Hafnia" & x$species == "alvei")
            | (x$genus == "Morganella" & x$species == "morganii")
            | (x$genus == "Proteus" & x$species == "penneri")
            | (x$genus == "Proteus" & x$species == "vulgaris")
            | (x$genus == "Serratia" & x$species == "marcescens")), ab_NA(CZO)] <- NA
    x[which((x$genus == "Morganella" & x$species == "morganii")
            | (x$genus == "Proteus" & x$species == "penneri")
            | (x$genus == "Proteus" & x$species == "vulgaris")
            | (x$genus == "Serratia" & x$species == "marcescens")), ab_NA(CXM)] <- NA
    x[which((x$genus == "Morganella" & x$species == "morganii")
            | (x$genus == "Proteus" & x$species == "mirabilis")
            | (x$genus == "Proteus" & x$species == "penneri")
            | (x$genus == "Proteus" & x$species == "vulgaris")
            | (x$genus == "Providencia" & x$species == "rettgeri")
            | (x$genus == "Providencia" & x$species == "stuartii")), ab_NA(TGC)] <- NA
    x[which((x$genus == "Citrobacter" & x$species == "koseri")
            | (x$genus == "Citrobacter" & x$species == "freundii")
            | (x$genus == "Enterobacter" & x$species == "aerogenes")
            | (x$genus == "Enterobacter" & x$species == "cloacae")
            | (x$genus == "Escherichia" & x$species == "hermannii")
            | (x$genus == "Hafnia" & x$species == "alvei")
            | (x$genus == "Klebsiella")
            | (x$genus == "Morganella" & x$species == "morganii")
            | (x$genus == "Proteus" & x$species == "penneri")
            | (x$genus == "Proteus" & x$species == "vulgaris")
            | (x$genus == "Providencia" & x$species == "rettgeri")
            | (x$genus == "Providencia" & x$species == "stuartii")
            | (x$genus == "Serratia" & x$species == "marcescens")), ab_NA(AMP)] <- NA
    x[which((x$genus == "Citrobacter" & x$species == "freundii")
            | (x$genus == "Enterobacter" & x$species == "aerogenes")
            | (x$genus == "Enterobacter" & x$species == "cloacae")
            | (x$genus == "Hafnia" & x$species == "alvei")
            | (x$genus == "Morganella" & x$species == "morganii")
            | (x$genus == "Providencia" & x$species == "rettgeri")
            | (x$genus == "Providencia" & x$species == "stuartii")
            | (x$genus == "Serratia" & x$species == "marcescens")), ab_NA(AMC)] <- NA
    x[which((x$genus == "Citrobacter" & x$species == "freundii")
            | (x$genus == "Citrobacter" & x$species == "koseri")
            | (x$genus == "Enterobacter" & x$species == "aerogenes")
            | (x$genus == "Enterobacter" & x$species == "cloacae")
            | (x$genus == "Hafnia" & x$species == "alvei")
            | (x$genus == "Providencia" & x$species == "rettgeri")
            | (x$genus == "Serratia" & x$species == "marcescens")), ab_NA(SAM)] <- NA
    x[which((x$genus == "Morganella" & x$species == "morganii")
            | (x$genus == "Proteus" & x$species == "mirabilis")
            | (x$genus == "Proteus" & x$species == "penneri")
            | (x$genus == "Proteus" & x$species == "vulgaris")
            | (x$genus == "Providencia" & x$species == "rettgeri")
            | (x$genus == "Providencia" & x$species == "stuartii")
            | (x$genus == "Serratia" & x$species == "marcescens")), ab_NA(COL)] <- NA
    x[which((x$genus == "Morganella" & x$species == "morganii")
            | (x$genus == "Proteus" & x$species == "mirabilis")
            | (x$genus == "Proteus" & x$species == "penneri")
            | (x$genus == "Proteus" & x$species == "vulgaris")
            | (x$genus == "Providencia" & x$species == "rettgeri")
            | (x$genus == "Providencia" & x$species == "stuartii")), ab_NA(TCY)] <- NA
    x[which((x$genus == "Morganella" & x$species == "morganii")
            | (x$genus == "Proteus" & x$species == "penneri")
            | (x$genus == "Proteus" & x$species == "vulgaris")
            | (x$genus == "Providencia" & x$species == "rettgeri")
            | (x$genus == "Providencia" & x$species == "stuartii")), ab_NA(c(DOX, MNO))] <- NA
    
    x$classes_in_guideline <- NA_integer_
    x$classes_available <- NA_integer_
    x$classes_affected <- NA_integer_
    
    # now add the MDR levels to the data
    trans_tbl(2,
              which(x$genus == "Staphylococcus" & x$species == "aureus"),
              c(OXA, FOX),
              "any")
    trans_tbl2(paste("Table 1 -", italic("Staphylococcus aureus")),
               which(x$genus == "Staphylococcus" & x$species == "aureus"),
               list(GEN,
                    RIF,
                    CPT,
                    c(OXA, FOX),
                    c(CIP, MFX),
                    SXT,
                    FUS,
                    c(VAN, TEC, TLV),
                    TGC,
                    CLI,
                    DAP,
                    ERY,
                    LNZ,
                    CHL,
                    FOS,
                    QDA,
                    c(TCY, DOX, MNO)))
    trans_tbl2(paste("Table 2 -", italic("Enterococcus"), "spp."),
               which(x$genus == "Enterococcus"),
               list(GEH,
                    STH,
                    c(IPM, MEM, DOR),
                    c(CIP, LVX, MFX),
                    c(VAN, TEC),
                    TGC,
                    DAP,
                    LNZ,
                    AMP,
                    QDA,
                    c(DOX, MNO)))
    trans_tbl2(paste0("Table 3 - ", italic("Enterobacteriaceae")),
               # this new order was previously 'Enterobacteriales' and contained only the family 'Enterobacteriaceae':
               which(x$order == "Enterobacterales"),
               list(c(GEN, TOB, AMK, NET),
                    CPT,
                    c(TCC, TZP),
                    c(ETP, IPM, MEM, DOR),
                    CZO,
                    CXM,
                    c(CTX, CAZ, FEP),
                    c(FOX, CTT),
                    CIP,
                    SXT,
                    TGC,
                    ATM,
                    AMP,
                    c(AMC, SAM),
                    CHL,
                    FOS,
                    COL,
                    c(TCY, DOX, MNO)))
    trans_tbl2(paste("Table 4 -", italic("Pseudomonas aeruginosa")),
               which(x$genus == "Pseudomonas" & x$species == "aeruginosa"),
               list(c(GEN, TOB, AMK, NET),
                    c(IPM, MEM, DOR),
                    c(CAZ, FEP),
                    c(CIP, LVX),
                    c(TCC, TZP),
                    ATM,
                    FOS,
                    c(COL, PLB)))
    trans_tbl2(paste("Table 5 -", italic("Acinetobacter"), "spp."),
               which(x$genus == "Acinetobacter"),
               list(c(GEN, TOB, AMK, NET),
                    c(IPM, MEM, DOR),
                    c(CIP, LVX),
                    c(TZP, TCC),
                    c(CTX, CRO, CAZ, FEP),
                    SXT,
                    SAM,
                    c(COL, PLB),
                    c(TCY, DOX, MNO)))
    
    # now set MDROs: 
    # MDR (=2): >=3 classes affected
    x[which(x$classes_affected >= 3), "MDRO"] <- 2
    if (verbose == TRUE) {
      x[which(x$classes_affected >= 3), "reason"] <- paste0("at least 3 classes contain R or I: ", x$classes_affected[which(x$classes_affected >= 3)],
                                                            " out of ", x$classes_available[which(x$classes_affected >= 3)], " available classes")
    }
    
    # XDR (=3): all but <=2 classes affected
    x[which((x$classes_in_guideline - x$classes_affected) <= 2), "MDRO"] <- 3
    if (verbose == TRUE) {
      x[which(x$MDRO == 3), "reason"] <- paste0("less than 3 classes remain susceptible (", x$classes_in_guideline[which((x$classes_in_guideline - x$classes_affected) <= 2)] - x$classes_affected[which(x$MDRO == 3)],
                                                                                       " out of ", x$classes_in_guideline[which(x$MDRO == 3)], " classes)")
    }
    
    # PDR (=4): all agents are R 
    x[which(x$classes_affected == 999 & x$classes_in_guideline == x$classes_available), "MDRO"] <- 4
    if (verbose == TRUE) {
      x[which(x$MDRO == 4), "reason"] <- paste("all antibiotics in all", x$classes_in_guideline[which(x$MDRO == 4)], "classes were tested R or I")
    }
    
    # not enough classes available
    x[which(x$MDRO %in% c(1, 3) & x$classes_available < base::floor(x$classes_in_guideline * pct_required_classes)), "MDRO"] <- -1
    if (verbose == TRUE) {
      x[which(x$MDRO == -1), "reason"] <- paste0("not enough classes available: ", x$classes_available[which(x$MDRO == -1)], 
                                                 " of required ", (base::floor(x$classes_in_guideline * pct_required_classes))[which(x$MDRO == -1)], 
                                                 " (~", percentage(pct_required_classes), " of ", x$classes_in_guideline[which(x$MDRO == -1)], ")")
    }
    
    # add antibiotic names of resistant ones to verbose output
    
  }
  
  if (guideline$code == "eucast") {
    # EUCAST ------------------------------------------------------------------
    # Table 5
    trans_tbl(3,
              which(x$order == "Enterobacterales"
                    | x$fullname %like% "^Pseudomonas aeruginosa"
                    | x$genus == "Acinetobacter"),
              COL,
              "all")
    trans_tbl(3,
              which(x$fullname %like% "^Salmonella Typhi"),
              c(carbapenems, fluoroquinolones),
              "any")
    trans_tbl(3,
              which(x$fullname %like% "^Haemophilus influenzae"),
              c(cephalosporins_3rd, carbapenems, fluoroquinolones),
              "any")
    trans_tbl(3,
              which(x$fullname %like% "^Moraxella catarrhalis"),
              c(cephalosporins_3rd, fluoroquinolones),
              "any")
    trans_tbl(3,
              which(x$fullname %like% "^Neisseria meningitidis"),
              c(cephalosporins_3rd, fluoroquinolones),
              "any")
    trans_tbl(3,
              which(x$fullname %like% "^Neisseria gonorrhoeae"),
              AZM,
              "any")
    # Table 6
    trans_tbl(3,
              which(x$fullname %like% "^(Coagulase-negative|Staphylococcus (aureus|epidermidis|hominis|haemolyticus|intermedius|pseudointermedius))"),
              c(VAN, TEC, DAP, LNZ, QDA, TGC),
              "any")
    trans_tbl(3,
              which(x$genus == "Corynebacterium"),
              c(VAN, TEC, DAP, LNZ, QDA, TGC),
              "any")
    trans_tbl(3,
              which(x$fullname %like% "^Streptococcus pneumoniae"),
              c(carbapenems, VAN, TEC, DAP, LNZ, QDA, TGC, RIF),
              "any")
    trans_tbl(3, # Sr. groups A/B/C/G
              which(x$fullname %like% "^Streptococcus (group (A|B|C|G)|pyogenes|agalactiae|equisimilis|equi|zooepidemicus|dysgalactiae|anginosus)"),
              c(PEN, cephalosporins, VAN, TEC, DAP, LNZ, QDA, TGC),
              "any")
    trans_tbl(3,
              which(x$genus == "Enterococcus"),
              c(DAP, LNZ, TGC, TEC),
              "any")
    trans_tbl(3,
              which(x$fullname %like% "^Enterococcus faecalis"),
              c(AMP, AMX),
              "any")
    # Table 7
    trans_tbl(3,
              which(x$genus == "Bacteroides"),
              MTR,
              "any")
    trans_tbl(3,
              which(x$fullname %like% "^Clostridium difficile"),
              c(MTR, VAN),
              "any")
  }
  
  if (guideline$code == "mrgn") {
    # Germany -----------------------------------------------------------------
    CTX_or_CAZ <- CTX %or% CAZ
    IPM_or_MEM <- IPM %or% MEM
    x$missing <- NA_character_
    if (is.na(PIP)) PIP <- "missing"
    if (is.na(CTX_or_CAZ)) CTX_or_CAZ <- "missing"
    if (is.na(IPM_or_MEM)) IPM_or_MEM <- "missing"
    if (is.na(IPM)) IPM <- "missing"
    if (is.na(MEM)) MEM <- "missing"
    if (is.na(CIP)) CIP <- "missing"
    
    # Table 1
    x[which((x$order == "Enterobacterales" |  # following in fact the old Enterobacteriaceae classification
               x$fullname %like% "^Acinetobacter baumannii") &
              x[, PIP] == "R" &
              x[, CTX_or_CAZ] == "R" &
              x[, IPM_or_MEM] == "S" &
              x[, CIP] == "R"),
      "MDRO"] <- 2 # 2 = 3MRGN
    
    x[which((x$order == "Enterobacterales" |  # following in fact the old Enterobacteriaceae classification
               x$fullname %like% "^Acinetobacter baumannii") &
              x[, PIP] == "R" &
              x[, CTX_or_CAZ] == "R" &
              x[, IPM_or_MEM] == "R" &
              x[, CIP] == "R"),
      "MDRO"] <- 3 # 3 = 4MRGN, overwrites 3MRGN if applicable
    
    x[which((x$order == "Enterobacterales" |  # following in fact the old Enterobacteriaceae classification
               x$fullname %like% "^Acinetobacter baumannii") &
              x[, IPM] == "R" | x[, MEM] == "R"),
      "MDRO"] <- 3 # 3 = 4MRGN, always when imipenem or meropenem is R
    
    x[which(x$fullname %like% "^Pseudomonas aeruginosa" &
              (x[, PIP] == "S") +
              (x[, CTX_or_CAZ] == "S") +
              (x[, IPM_or_MEM] == "S") +
              (x[, CIP] == "S") == 1),
      "MDRO"] <- 2 # 2 = 3MRGN, if only 1 group is S
    
    x[which((x$fullname %like% "^Pseudomonas aeruginosa") &
              x[, PIP] == "R" &
              x[, CTX_or_CAZ] == "R" &
              x[, IPM_or_MEM] == "R" &
              x[, CIP] == "R"),
      "MDRO"] <- 3 # 3 = 4MRGN
  }
  
  if (guideline$code == "brmo") {
    # Netherlands -------------------------------------------------------------
    aminoglycosides <- aminoglycosides[!is.na(aminoglycosides)]
    fluoroquinolones <- fluoroquinolones[!is.na(fluoroquinolones)]
    carbapenems <- carbapenems[!is.na(carbapenems)]
    amino <- AMX %or% AMP
    third <- CAZ %or% CTX
    ESBLs <- c(amino, third)
    ESBLs <- ESBLs[!is.na(ESBLs)]
    if (length(ESBLs) != 2) {
      ESBLs <- character(0)
    }
    
    # Table 1
    trans_tbl(3,
              which(x$order == "Enterobacterales"), # following in fact the old Enterobacteriaceae classification
              c(aminoglycosides, fluoroquinolones),
              "all")
    
    trans_tbl(2,
              which(x$order == "Enterobacterales"), # following in fact the old Enterobacteriaceae classification
              carbapenems,
              "any")
    
    trans_tbl(2,
              which(x$order == "Enterobacterales"), # following in fact the old Enterobacteriaceae classification
              ESBLs,
              "all")
    
    # Table 2
    trans_tbl(2,
              which(x$genus == "Acinetobacter"),
              c(carbapenems),
              "any")
    trans_tbl(3,
              which(x$genus == "Acinetobacter"),
              c(aminoglycosides, fluoroquinolones),
              "all")
    
    trans_tbl(3,
              which(x$fullname %like% "^Stenotrophomonas maltophilia"),
              SXT,
              "all")
    
    if (!ab_missing(MEM) & !ab_missing(IPM)
        & !ab_missing(GEN) & !ab_missing(TOB)
        & !ab_missing(CIP)
        & !ab_missing(CAZ)
        & !ab_missing(TZP)) {
      x$psae <- 0
      x[which(x[, MEM] == "R" | x[, IPM] == "R"), "psae"] <- 1 + x[which(x[, MEM] == "R" | x[, IPM] == "R"), "psae"]
      x[which(x[, GEN] == "R" & x[, TOB] == "R"), "psae"] <- 1 + x[which(x[, GEN] == "R" & x[, TOB] == "R"), "psae"]
      x[which(x[, CIP] == "R"), "psae"] <- 1 + x[which(x[, CIP] == "R"), "psae"]
      x[which(x[, CAZ] == "R"), "psae"] <- 1 + x[which(x[, CAZ] == "R"), "psae"]
      x[which(x[, TZP] == "R"), "psae"] <- 1 + x[which(x[, TZP] == "R"), "psae"]
    } else {
      x$psae <- 0
    }
    x[which(
      x$fullname %like% "Pseudomonas aeruginosa"
      & x$psae >= 3), "MDRO"] <- 3
    
    # Table 3
    trans_tbl(3,
              which(x$fullname %like% "Streptococcus pneumoniae"),
              PEN,
              "all")
    trans_tbl(3,
              which(x$fullname %like% "Streptococcus pneumoniae"),
              VAN,
              "all")
    trans_tbl(3,
              which(x$fullname %like% "Enterococcus faecium"),
              c(PEN, VAN),
              "all")
  }
  
  prepare_drug <- function(ab) {
    # returns vector values of drug
    # if `ab` is a column name, looks up the values in `x`
    if (length(ab) == 1 & is.character(ab)) {
      if (ab %in% colnames(x)) {
        ab <- as.data.frame(x)[, ab]
      }
    }
    ab <- as.character(as.rsi(ab))
    ab[is.na(ab)] <- ""
    ab
  }
  drug_is_R <- function(ab) {
    # returns logical vector
    ab <- prepare_drug(ab)
    if (length(ab) == 1) {
      rep(ab, NROW(x)) == "R"
    } else {
      ab == "R"
    }
  }
  drug_is_not_R <- function(ab) {
    # returns logical vector
    ab <- prepare_drug(ab)
    if (length(ab) == 1) {
      rep(ab, NROW(x)) != "R"
    } else {
      ab != "R"
    }
  }
  
  if (guideline$code == "tb") {
    # Tuberculosis ------------------------------------------------------------
    x <- x %>%
      mutate(mono_count = 0,
             mono_count = ifelse(drug_is_R(INH), mono_count + 1, mono_count),
             mono_count = ifelse(drug_is_R(RIF), mono_count + 1, mono_count),
             mono_count = ifelse(drug_is_R(ETH), mono_count + 1, mono_count),
             mono_count = ifelse(drug_is_R(PZA), mono_count + 1, mono_count),
             mono_count = ifelse(drug_is_R(RIB), mono_count + 1, mono_count),
             mono_count = ifelse(drug_is_R(RFP), mono_count + 1, mono_count),
             # from here on logicals
             mono = mono_count > 0,
             poly = ifelse(mono_count > 1 & drug_is_not_R(RIF) & drug_is_not_R(INH),
                           TRUE, FALSE),
             mdr = ifelse(drug_is_R(RIF) & drug_is_R(INH),
                          TRUE, FALSE),
             xdr = ifelse(drug_is_R(LVX) | drug_is_R(MFX) | drug_is_R(GAT),
                          TRUE, FALSE),
             second = ifelse(drug_is_R(CAP) | drug_is_R(KAN) | drug_is_R(AMK),
                             TRUE, FALSE),
             xdr = ifelse(mdr & xdr & second, TRUE, FALSE)) %>%
      mutate(MDRO = case_when(xdr ~ 5,
                              mdr ~ 4,
                              poly ~ 3,
                              mono ~ 2,
                              TRUE ~ 1),
             # keep all real TB, make other species NA
             MDRO = ifelse(x$fullname == "Mycobacterium tuberculosis", MDRO, NA_real_))
  }
  
  if (info == TRUE) {
    cat(bold(paste0("=> Found ", sum(x$MDRO %in% c(2:5), na.rm = TRUE), " MDROs out of ", sum(!is.na(x$MDRO)), 
                    " tested isolates (", percentage(sum(x$MDRO %in% c(2:5), na.rm = TRUE) / sum(!is.na(x$MDRO))), ")\n")))
  }
  
  # some more info on negative results
  if (verbose == TRUE) {
    if (guideline$code == "cmi2012") {
        x[which(x$MDRO == 1 & !is.na(x$classes_affected)), "reason"] <- paste0(x$classes_affected[which(x$MDRO == 1 & !is.na(x$classes_affected))], " of ", x$classes_available[which(x$MDRO == 1 & !is.na(x$classes_affected))], " available classes contain R or I (3 required for MDR)")
    } else {
      x[which(x$MDRO == 1), "reason"] <- "too few antibiotics are R"
    }
  }
  
  # Results ----
  if (guideline$code == "cmi2012") {
    if (any(x$MDRO == -1, na.rm = TRUE)) {
      warning("NA introduced for isolates where the available percentage of antimicrobial classes was below ",
              percentage(pct_required_classes), " (set with `pct_required_classes`)")
      # set these -1s to NA
      x[which(x$MDRO == -1), "MDRO"] <- NA_integer_
    }
    x$MDRO <- factor(x = x$MDRO,
                     levels = 1:4,
                     labels = c("Negative", "Multi-drug-resistant (MDR)", 
                                "Extensively drug-resistant (XDR)", "Pandrug-resistant (PDR)"),
                     ordered = TRUE)
  } else if (guideline$code == "tb") {
    x$MDRO <- factor(x = x$MDRO,
                     levels = 1:5,
                     labels = c("Negative", "Mono-resistant", "Poly-resistant", 
                                "Multi-drug-resistant", "Extensively drug-resistant"),
                     ordered = TRUE)
  } else if (guideline$code == "mrgn") {
    x$MDRO <- factor(x = x$MDRO,
                     levels = 1:3,
                     labels = c("Negative", "3MRGN", "4MRGN"),
                     ordered = TRUE)
  } else {
    x$MDRO <- factor(x = x$MDRO,
                     levels = 1:3,
                     labels = c("Negative", "Positive, unconfirmed", "Positive"),
                     ordered = TRUE)
  }
  
  if (verbose == TRUE) {
    x[, c("row_number",
          col_mo,
          "MDRO",
          "reason",
          "columns_nonsusceptible")]
    #x
  } else {
    x$MDRO
  }
  
}

#' @rdname mdro
#' @export
brmo <- function(x, guideline = "BRMO", ...) {
  mdro(x, guideline = "BRMO", ...)
}

#' @rdname mdro
#' @export
mrgn <- function(x, guideline = "MRGN", ...) {
  mdro(x = x, guideline = "MRGN", ...)
}

#' @rdname mdro
#' @export
mdr_tb <- function(x, guideline = "TB", ...) {
  mdro(x = x, guideline = "TB", ...)
}

#' @rdname mdro
#' @export
mdr_cmi2012 <- function(x, guideline = "CMI2012", ...) {
  mdro(x = x, guideline = "CMI2012", ...)
}


#' @rdname mdro
#' @export
eucast_exceptional_phenotypes <- function(x, guideline = "EUCAST", ...) {
  mdro(x = x, guideline = "EUCAST", ...)
}
