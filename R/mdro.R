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
#' Determine which isolates are multidrug-resistant organisms (MDRO) according to country-specific guidelines.
#' @param x table with antibiotic columns, like e.g. \code{AMX} and \code{AMC}
#' @param country country code to determine guidelines. EUCAST rules will be used when left empty, see Details. Should be or a code from the \href{https://en.wikipedia.org/wiki/ISO_3166-1_alpha-2#Officially_assigned_code_elements}{list of ISO 3166-1 alpha-2 country codes}. Case-insensitive. Currently supported are \code{de} (Germany) and \code{nl} (the Netherlands).
#' @param info print progress
#' @inheritParams eucast_rules
#' @param verbose print additional info: missing antibiotic columns per parameter
#' @inheritSection eucast_rules Antibiotics
#' @details When \code{country} will be left blank, guidelines will be taken from EUCAST Expert Rules Version 3.1 "Intrinsic Resistance and Exceptional Phenotypes Tables" (\href{http://www.eucast.org/fileadmin/src/media/PDFs/EUCAST_files/Expert_Rules/Expert_rules_intrinsic_exceptional_V3.1.pdf}{link}).
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
mdro <- function(x,
                 country = NULL,
                 col_mo = NULL,
                 info = TRUE,
                 verbose = FALSE,
                 ...) {

  tbl_ <- x

  if (!is.data.frame(tbl_)) {
    stop("`x` must be a data frame.", call. = FALSE)
  }

  # try to find columns based on type
  # -- mo
  if (is.null(col_mo)) {
    col_mo <- search_type_in_df(tbl = tbl_, type = "mo")
  }
  if (is.null(col_mo)) {
    stop("`col_mo` must be set.", call. = FALSE)
  }

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
    guideline$version <- 'Revision as of December 2017'
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


  cols_ab <- get_column_abx(tbl = x,
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
  MTR <- cols_ab['MTR']
  MFX <- cols_ab['MFX']
  MNO <- cols_ab['MNO']
  NAL <- cols_ab['NAL']
  NEO <- cols_ab['NEO']
  NET <- cols_ab['NET']
  NIT <- cols_ab['NIT']
  NOR <- cols_ab['NOR']
  NOV <- cols_ab['NOV']
  OFX <- cols_ab['OFX']
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
        col_filter <- which(tbl_[, cols] == 'R')
      } else if (any_all == "all") {
        col_filter <- tbl_ %>%
          mutate(index = 1:nrow(.)) %>%
          filter_at(vars(cols), all_vars(. == "R")) %>%
          pull((index))
      }
      rows <- rows[rows %in% col_filter]
      tbl_[rows, 'MDRO'] <<- to
    }
  }

  tbl_ <- tbl_ %>%
    mutate_at(vars(col_mo), as.mo) %>%
    # join to microorganisms data set
    left_join_microorganisms(by = col_mo) %>%
    # add unconfirmed to where genus is available
    mutate(MDRO = ifelse(!is.na(genus), 1, NA_integer_))

  if (guideline$country$code == 'eucast') {
    # EUCAST ------------------------------------------------------------------
    # Table 5
    trans_tbl(3,
              which(tbl_$family == 'Enterobacteriaceae'
                    | tbl_$fullname %like% '^Pseudomonas aeruginosa'
                    | tbl_$genus == 'Acinetobacter'),
              COL,
              "all")
    trans_tbl(3,
              which(tbl_$fullname %like% '^Salmonella Typhi'),
              c(carbapenems, fluoroquinolones),
              "any")
    trans_tbl(3,
              which(tbl_$fullname %like% '^Haemophilus influenzae'),
              c(cephalosporins_3rd, carbapenems, fluoroquinolones),
              "any")
    trans_tbl(3,
              which(tbl_$fullname %like% '^Moraxella catarrhalis'),
              c(cephalosporins_3rd, fluoroquinolones),
              "any")
    trans_tbl(3,
              which(tbl_$fullname %like% '^Neisseria meningitidis'),
              c(cephalosporins_3rd, fluoroquinolones),
              "any")
    trans_tbl(3,
              which(tbl_$fullname %like% '^Neisseria gonorrhoeae'),
              AZM,
              "any")
    # Table 6
    trans_tbl(3,
              which(tbl_$fullname %like% '^Staphylococcus (aureus|epidermidis|coagulase negatief|hominis|haemolyticus|intermedius|pseudointermedius)'),
              c(VAN, TEC, DAP, LNZ, QDA, TGC),
              "any")
    trans_tbl(3,
              which(tbl_$genus == 'Corynebacterium'),
              c(VAN, TEC, DAP, LNZ, QDA, TGC),
              "any")
    trans_tbl(3,
              which(tbl_$fullname %like% '^Streptococcus pneumoniae'),
              c(carbapenems, VAN, TEC, DAP, LNZ, QDA, TGC, RIF),
              "any")
    trans_tbl(3, # Sr. groups A/B/C/G
              which(tbl_$fullname %like% '^Streptococcus (pyogenes|agalactiae|equisimilis|equi|zooepidemicus|dysgalactiae|anginosus)'),
              c(PEN, cephalosporins, VAN, TEC, DAP, LNZ, QDA, TGC),
              "any")
    trans_tbl(3,
              which(tbl_$genus == 'Enterococcus'),
              c(DAP, LNZ, TGC, TEC),
              "any")
    trans_tbl(3,
              which(tbl_$fullname %like% '^Enterococcus faecalis'),
              c(AMP, AMX),
              "any")
    # Table 7
    trans_tbl(3,
              which(tbl_$genus == 'Bacteroides'),
               MTR,
              "any")
    trans_tbl(3,
              which(tbl_$fullname %like% '^Clostridium difficile'),
              c( MTR, VAN),
              "any")
  }

  if (guideline$country$code == 'de') {
    # Germany -----------------------------------------------------------------
    stop("We are still working on German guidelines in this beta version.", call. = FALSE)
  }

  if (guideline$country$code == 'nl') {
    # Netherlands -------------------------------------------------------------
    aminoglycosides <- aminoglycosides[!ab_missing(aminoglycosides)]
    fluoroquinolones <- fluoroquinolones[!ab_missing(fluoroquinolones)]
    carbapenems <- carbapenems[!ab_missing(carbapenems)]

    # Table 1
    trans_tbl(3,
              which(tbl_$family == 'Enterobacteriaceae'),
              c(aminoglycosides, fluoroquinolones),
              "all")

    trans_tbl(2,
              which(tbl_$family == 'Enterobacteriaceae'),
              c(carbapenems),
              "any")

    # Table 2
    trans_tbl(2,
              which(tbl_$genus == 'Acinetobacter'),
              c(carbapenems),
              "any")
    trans_tbl(3,
              which(tbl_$genus == 'Acinetobacter'),
              c(aminoglycosides, fluoroquinolones),
              "all")

    trans_tbl(3,
              which(tbl_$fullname %like% '^Stenotrophomonas maltophilia'),
              SXT,
              "all")

    if (!ab_missing(MEM) & !ab_missing(IPM)
        & !ab_missing(GEN) & !ab_missing(TOB)
        & !ab_missing(CIP)
        & !ab_missing(CAZ)
        & !ab_missing(TZP) ) {
      tbl_$psae <- 0
      tbl_[which(tbl_[, MEM] == "R" | tbl_[, IPM] == "R"), "psae"] <- 1 + tbl_[which(tbl_[, MEM] == "R" | tbl_[, IPM] == "R"), "psae"]
      tbl_[which(tbl_[, GEN] == "R" & tbl_[, TOB] == "R"), "psae"] <- 1 + tbl_[which(tbl_[, GEN] == "R" & tbl_[, TOB] == "R"), "psae"]
      tbl_[which(tbl_[, CIP] == "R"), "psae"] <- 1 + tbl_[which(tbl_[, CIP] == "R"), "psae"]
      tbl_[which(tbl_[, CAZ] == "R"), "psae"] <- 1 + tbl_[which(tbl_[, CAZ] == "R"), "psae"]
      tbl_[which(tbl_[, TZP] == "R"), "psae"] <- 1 + tbl_[which(tbl_[, TZP] == "R"), "psae"]
    } else {
      tbl_$psae <- 0
    }
    tbl_[which(
      tbl_$fullname %like% 'Pseudomonas aeruginosa'
      & tbl_$psae >= 3
    ), 'MDRO'] <- 3

    # Table 3
    trans_tbl(3,
              which(tbl_$fullname %like% 'Streptococcus pneumoniae'),
              PEN,
              "all")
    trans_tbl(3,
              which(tbl_$fullname %like% 'Streptococcus pneumoniae'),
              VAN,
              "all")
    trans_tbl(3,
              which(tbl_$fullname %like% 'Enterococcus faecium'),
              c(PEN, VAN),
              "all")
  }

  factor(x =  tbl_$MDRO,
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
mrgn <- function(x, country = "de", ...) {
  mdro(x = x, country = "de", ...)
}

#' @rdname mdro
#' @export
eucast_exceptional_phenotypes <- function(x, country = "EUCAST", ...) {
  mdro(x = x, country = "EUCAST", ...)
}

is_ESBL <- function(x, col_mo = NULL, ...) {
  col_mo <- get_column_mo(tbl = x, col_mo = col_mo)
  cols_ab <- get_column_abx(tbl = x,
                            soft_dependencies = c("AMX", "AMP"),
                            hard_dependencies = c("CAZ"),
                            ...)

  if (!any(c("AMX", "AMP") %in% names(cols_ab))) {
    # both ampicillin and amoxicillin are missing
    generate_warning_abs_missing(c("AMX", "AMP"), any = TRUE)
    return(rep(NA, nrow(x)))
  }

  ESBLs <- rep(NA, nrow(x))

  # first make all eligible cases FALSE
  ESBLs[which(mo_family(x[, col_mo]) == "Enterobacteriaceae"
              & x[, get_ab_col(cols_ab, "AMX")] %in% c("R", "I", "S")
              & x[, get_ab_col(cols_ab, "AMX")] %in% c("R", "I", "S")
              & x[, get_ab_col(cols_ab, "AMX")] %in% c("R", "I", "S")
              )] <- FALSE
  # now make the positives cases TRUE
  ESBLs[which(!is.na(ESBLs)
              & x[, get_ab_col(cols_ab, "AMX")]  == "R"
              & x[, get_ab_col(cols_ab, "CAZ")] == "R")] <- TRUE
  ESBLs

}

is_3MRGN <- function(x, ...) {

}

is_4MRGN <- function(x, ...) {

}

get_column_mo <- function(tbl, col_mo = NULL) {
  # throws a blue note about which column will be used if guessed
  if (is.null(col_mo)) {
    col_mo <- search_type_in_df(tbl = tbl, type = "mo")
  }
  if (is.null(col_mo)) {
    stop("`col_mo` must be set.", call. = FALSE)
  }
  col_mo
}

