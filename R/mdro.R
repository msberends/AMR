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
#' Determine which isolates are multidrug-resistant organisms (MDRO) according to (country-specific) guidelines.
#' @param x table with antibiotic columns, like e.g. \code{AMX} and \code{AMC}
#' @param guideline a specific guideline to mention, see Details. EUCAST guidelines will be used when left empty, see Details.
#' @param info print progress
#' @inheritParams eucast_rules
#' @param verbose print additional info: missing antibiotic columns per parameter
#' @inheritSection eucast_rules Antibiotics
#' @details Currently supported guidelines are:
#' \itemize{
#'   \item{\code{guideline = "EUCAST"}: EUCAST Expert Rules Version 3.1 "Intrinsic Resistance and Exceptional Phenotypes Tables" (\href{http://www.eucast.org/fileadmin/src/media/PDFs/EUCAST_files/Expert_Rules/Expert_rules_intrinsic_exceptional_V3.1.pdf}{link})}
#'   \item{\code{guideline = "TB"}: World Health Organization "Companion handbook to the WHO guidelines for the programmatic management of drug-resistant tuberculosis" (\href{https://www.who.int/tb/publications/pmdt_companionhandbook/en/}{link})}
#'   \item{\code{guideline = "MRGN"}: (work in progress)}
#'   \item{\code{guideline = "BRMO"}: Rijksinstituut voor Volksgezondheid en Milieu "WIP-richtlijn BRMO (Bijzonder Resistente Micro-Organismen) [ZKH]" (\href{https://www.rivm.nl/Documenten_en_publicaties/Professioneel_Praktisch/Richtlijnen/Infectieziekten/WIP_Richtlijnen/WIP_Richtlijnen/Ziekenhuizen/WIP_richtlijn_BRMO_Bijzonder_Resistente_Micro_Organismen_ZKH}{link})}
#' }
#'
#' Please suggest your own (country-specific) guidelines by letting us know: \url{https://gitlab.com/msberends/AMR/issues/new}.
#' @return For TB (\code{mdr_tb()}): Ordered factor with levels \code{Negative < Mono-resistance < Poly-resistance < Multidrug resistance < Extensive drug resistance}.
#'
#' For everything else: Ordered factor with levels \code{Negative < Positive, unconfirmed < Positive}. The value \code{"Positive, unconfirmed"} means that, according to the guideline, it is not entirely sure if the isolate is multi-drug resistant and this should be confirmed with additional (e.g. molecular) tests.
#' @rdname mdro
#' @importFrom dplyr %>%
#' @importFrom crayon red blue bold
#' @export
#' @inheritSection AMR Read more on our website!
#' @source
#' EUCAST Expert Rules Version 3.1 "Intrinsic Resistance and Exceptional Phenotypes Tables" (\href{http://www.eucast.org/fileadmin/src/media/PDFs/EUCAST_files/Expert_Rules/Expert_rules_intrinsic_exceptional_V3.1.pdf}{link})
#'
#' World Health Organization "Companion handbook to the WHO guidelines for the programmatic management of drug-resistant tuberculosis" (\href{https://www.who.int/tb/publications/pmdt_companionhandbook/en/}{link})
#'
#' Rijksinstituut voor Volksgezondheid en Milieu "WIP-richtlijn BRMO (Bijzonder Resistente Micro-Organismen) [ZKH]" (\href{https://www.rivm.nl/Documenten_en_publicaties/Professioneel_Praktisch/Richtlijnen/Infectieziekten/WIP_Richtlijnen/WIP_Richtlijnen/Ziekenhuizen/WIP_richtlijn_BRMO_Bijzonder_Resistente_Micro_Organismen_ZKH}{link})
#' @examples
#' library(dplyr)
#'
#' septic_patients %>%
#'   mutate(EUCAST = mdro(.),
#'          BRMO = brmo(.))
mdro <- function(x,
                 guideline = NULL,
                 col_mo = NULL,
                 info = TRUE,
                 verbose = FALSE,
                 ...) {

  if (!is.data.frame(x)) {
    stop("`x` must be a data frame.", call. = FALSE)
  }


  if (!is.null(list(...)$country)) {
    warning("Using `country` is deprecated, use `guideline` instead. Please see ?mdro.", call. = FALSE)
    guideline <- list(...)$country
  }
  if (length(guideline) > 1) {
    stop("`guideline` must be a length one character string.", call. = FALSE)
  }

  if (is.null(guideline)) {
    guideline <- "eucast"
  }
  if (tolower(guideline) == "nl") {
    guideline <- "BRMO"
  }
  if (tolower(guideline) == "de") {
    guideline <- "MRGN"
  }
  if (!tolower(guideline) %in% c("brmo", "mrgn", "eucast", "tb")) {
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

  if (guideline$code == "eucast") {
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
    guideline$name <- "Germany"
    guideline$name <- ""
    guideline$version <- ""
    guideline$source <- ""
  } else if (guideline$code == "brmo") {
    guideline$name <- "WIP-Richtlijn Bijzonder Resistente Micro-organismen (BRMO)"
    guideline$author <- "RIVM (Rijksinstituut voor de Volksgezondheid)"
    guideline$version <- "Revision as of December 2017"
    guideline$source <- "https://www.rivm.nl/Documenten_en_publicaties/Professioneel_Praktisch/Richtlijnen/Infectieziekten/WIP_Richtlijnen/WIP_Richtlijnen/Ziekenhuizen/WIP_richtlijn_BRMO_Bijzonder_Resistente_Micro_Organismen_ZKH"
  } else {
    stop("This guideline is currently unsupported: ", guideline$code, call. = FALSE)
  }

  if (info == TRUE) {
    cat("Determining multidrug-resistant organisms (MDRO), according to:\n",
        "Guideline: ", red(guideline$name), "\n",
        "Version:   ", red(guideline$version), "\n",
        "Author:    ",    red(guideline$author), "\n",
        "Source:    ", blue(guideline$source), "\n",
        "\n", sep = "")
  }

  if (guideline$code == "tb") {
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
  PEN <- cols_ab["PEN"]
  PIP <- cols_ab["PIP"]
  PLB <- cols_ab["PLB"]
  PRI <- cols_ab["PRI"]
  QDA <- cols_ab["QDA"]
  RID <- cols_ab["RID"]
  RIF <- cols_ab["RIF"]
  RXT <- cols_ab["RXT"]
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

  ab_missing <- function(ab) {
    isTRUE(ab %in% c(NULL, NA)) | length(ab) == 0
  }

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
      if (any_all == "any") {
        row_filter <- which(x[, cols] == "R")
      } else if (any_all == "all") {
        row_filter <- x %>%
          mutate(index = 1:nrow(.)) %>%
          filter_at(vars(cols), all_vars(. == "R")) %>%
          pull((index))
      }
      rows <- rows[rows %in% row_filter]
      x[rows, "MDRO"] <<- to
    }
  }

  x <- x %>%
    mutate_at(vars(col_mo), as.mo) %>%
    # join to microorganisms data set
    left_join_microorganisms(by = col_mo) %>%
    # add unconfirmed to where genus is available
    mutate(MDRO = ifelse(!is.na(genus), 1, NA_integer_))

  if (guideline$code == "eucast") {
    # EUCAST ------------------------------------------------------------------
    # Table 5
    trans_tbl(3,
              which(x$family == "Enterobacteriaceae"
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
              which(x$fullname %like% "^Staphylococcus (aureus|epidermidis|coagulase negatief|hominis|haemolyticus|intermedius|pseudointermedius)"),
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
              which(x$fullname %like% "^Streptococcus (pyogenes|agalactiae|equisimilis|equi|zooepidemicus|dysgalactiae|anginosus)"),
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
    stop("We are still working on German guidelines in this beta version.", call. = FALSE)
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
              which(x$family == "Enterobacteriaceae"),
              c(aminoglycosides, fluoroquinolones),
              "all")

    trans_tbl(2,
              which(x$family == "Enterobacteriaceae"),
              carbapenems,
              "any")

    trans_tbl(2,
              which(x$family == "Enterobacteriaceae"),
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
        & !ab_missing(TZP) ) {
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
      & x$psae >= 3
    ), "MDRO"] <- 3

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
      mutate(mdr_tb = case_when(xdr ~ 5,
                                mdr ~ 4,
                                poly ~ 3,
                                mono ~ 2,
                                TRUE ~ 1),
             # keep all real TB, make other species NA
             mdr_tb = ifelse(x$fullname == "Mycobacterium tuberculosis", mdr_tb, NA_real_))
  }

  # return results
  if (guideline$code == "tb") {
    factor(x =  x$mdr_tb,
           levels = 1:5,
           labels = c("Negative", "Mono-resistance", "Poly-resistance", "Multidrug resistance", "Extensive drug resistance"),
           ordered = TRUE)
  } else {
    factor(x =  x$MDRO,
           levels = 1:3,
           labels = c("Negative", "Positive, unconfirmed", "Positive"),
           ordered = TRUE)
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
eucast_exceptional_phenotypes <- function(x, guideline = "EUCAST", ...) {
  mdro(x = x, guideline = "EUCAST", ...)
}
