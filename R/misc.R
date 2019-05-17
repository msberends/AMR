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

# No export, no Rd
addin_insert_in <- function() {
  rstudioapi::insertText(" %in% ")
}

# No export, no Rd
addin_insert_like <- function() {
  rstudioapi::insertText(" %like% ")
}

# No export, no Rd
# works exactly like round(), but rounds `round(44.55, 1)` as 44.6 instead of 44.5
# and adds decimal zeroes until `digits` is reached when force_zero = TRUE
round2 <- function(x, digits = 0, force_zero = TRUE) {
  # https://stackoverflow.com/a/12688836/4575331
  val <- (trunc((abs(x) * 10 ^ digits) + 0.5) / 10 ^ digits) * sign(x)
  if (digits > 0 & force_zero == TRUE) {
    val[val != as.integer(val)] <- paste0(val[val != as.integer(val)],
                                          strrep("0", max(0, digits - nchar(gsub(".*[.](.*)$", "\\1", val[val != as.integer(val)])))))
  }
  val
}

# Coefficient of variation (CV)
cv <- function(x, na.rm = TRUE) {
  stats::sd(x, na.rm = na.rm) / base::abs(base::mean(x, na.rm = na.rm))
}

# Coefficient of dispersion, or coefficient of quartile variation (CQV).
# (Bonett et al., 2006: Confidence interval for a coefficient of quartile variation).
cqv <- function(x, na.rm = TRUE) {
  fives <- stats::fivenum(x, na.rm = na.rm)
  (fives[4] - fives[2]) / (fives[4] + fives[2])
}

# show bytes as kB/MB/GB
# size_humanreadable(123456) # 121 kB
# size_humanreadable(12345678) # 11.8 MB
size_humanreadable <- function(bytes, decimals = 1) {
  bytes <- bytes %>% as.double()
  # Adapted from:
  # http://jeffreysambells.com/2012/10/25/human-readable-filesize-php
  size <- c('B','kB','MB','GB','TB','PB','EB','ZB','YB')
  factor <- floor((nchar(bytes) - 1) / 3)
  # added slight improvement; no decimals for B and kB:
  decimals <- rep(decimals, length(bytes))
  decimals[size[factor + 1] %in% c('B', 'kB')] <- 0

  out <- paste(sprintf(paste0("%.", decimals, "f"), bytes / (1024 ^ factor)), size[factor + 1])
  out
}

percent_scales <- scales::percent
# No export, no Rd
# based on scales::percent
percent <- function(x, round = 1, force_zero = FALSE, decimal.mark = getOption("OutDec"), ...) {
  x <- percent_scales(x = as.double(x),
                      accuracy = 1 / 10 ^ round,
                      decimal.mark = decimal.mark,
                      ...)
  if (force_zero == FALSE) {
    x <- gsub("([.]%|%%)", "%", paste0(gsub("0+%$", "", x), "%"))
  }
  x
}

#' @importFrom crayon blue bold red
#' @importFrom dplyr %>% pull
search_type_in_df <- function(tbl, type) {
  # try to find columns based on type
  found <- NULL

  colnames(tbl) <- trimws(colnames(tbl))

  # -- mo
  if (type == "mo") {
    if ("mo" %in% lapply(tbl, class)) {
      found <- colnames(tbl)[lapply(tbl, class) == "mo"][1]
    } else if (any(colnames(tbl) %like% "^(mo|microorganism|organism|bacteria)")) {
      found <- colnames(tbl)[colnames(tbl) %like% "^(mo|microorganism|organism|bacteria)"][1]
    } else if (any(colnames(tbl) %like% "species")) {
      found <- colnames(tbl)[colnames(tbl) %like% "species"][1]
    }

  }
  # -- key antibiotics
  if (type == "keyantibiotics") {
    if (any(colnames(tbl) %like% "^key.*(ab|antibiotics)")) {
      found <- colnames(tbl)[colnames(tbl) %like% "^key.*(ab|antibiotics)"][1]
    }
  }
  # -- date
  if (type == "date") {
    if (any(colnames(tbl) %like% "^(specimen date|specimen_date|spec_date)")) {
      # WHONET support
      found <- colnames(tbl)[colnames(tbl) %like% "^(specimen date|specimen_date|spec_date)"][1]
      if (!any(class(tbl %>% pull(found)) %in% c("Date", "POSIXct"))) {
        stop(red(paste0("ERROR: Found column `", bold(found), "` to be used as input for `col_", type,
                        "`, but this column contains no valid dates. Transform its values to valid dates first.")),
             call. = FALSE)
      }
    } else {
      for (i in 1:ncol(tbl)) {
        if (any(class(tbl %>% pull(i)) %in% c("Date", "POSIXct"))) {
          found <- colnames(tbl)[i]
          break
        }
      }
    }
  }
  # -- patient id
  if (type == "patient_id") {
    if (any(colnames(tbl) %like% "^(identification |patient|patid)")) {
      found <- colnames(tbl)[colnames(tbl) %like% "^(identification |patient|patid)"][1]
    }
  }
  # -- specimen
  if (type == "specimen") {
    if (any(colnames(tbl) %like% "(specimen type|spec_type)")) {
      found <- colnames(tbl)[colnames(tbl) %like% "(specimen type|spec_type)"][1]
    } else if (any(colnames(tbl) %like% "^(specimen)")) {
      found <- colnames(tbl)[colnames(tbl) %like% "^(specimen)"][1]
    }
  }

  if (!is.null(found)) {
    msg <- paste0("NOTE: Using column `", bold(found), "` as input for `col_", type, "`.")
    if (type %in% c("keyantibiotics", "specimen")) {
      msg <- paste(msg, "Use", bold(paste0("col_", type), "= FALSE"), "to prevent this.")
    }
    message(blue(msg))
  }
  found
}

get_ab_col <- function(columns, ab) {
  columns[names(columns) == ab]
}

get_column_abx <- function(tbl,
                           soft_dependencies = NULL,
                           hard_dependencies = NULL,
                           verbose = FALSE,
                           AMC = guess_ab_col(),
                           AMK = guess_ab_col(),
                           AMX = guess_ab_col(),
                           AMP = guess_ab_col(),
                           AZM = guess_ab_col(),
                           AZL = guess_ab_col(),
                           ATM = guess_ab_col(),
                           RID = guess_ab_col(),
                           FEP = guess_ab_col(),
                           CTX = guess_ab_col(),
                           FOX = guess_ab_col(),
                           CED = guess_ab_col(),
                           CAZ = guess_ab_col(),
                           CRO = guess_ab_col(),
                           CXM = guess_ab_col(),
                           CHL = guess_ab_col(),
                           CIP = guess_ab_col(),
                           CLR = guess_ab_col(),
                           CLI = guess_ab_col(),
                           FLC = guess_ab_col(),
                           COL = guess_ab_col(),
                           CZO = guess_ab_col(),
                           DAP = guess_ab_col(),
                           DOX = guess_ab_col(),
                           ETP = guess_ab_col(),
                           ERY = guess_ab_col(),
                           FOS = guess_ab_col(),
                           FUS = guess_ab_col(),
                           GEN = guess_ab_col(),
                           IPM = guess_ab_col(),
                           KAN = guess_ab_col(),
                           LVX = guess_ab_col(),
                           LIN = guess_ab_col(),
                           LNZ = guess_ab_col(),
                           MEM = guess_ab_col(),
                           MTR = guess_ab_col(),
                           MEZ = guess_ab_col(),
                           MNO = guess_ab_col(),
                           MFX = guess_ab_col(),
                           NAL = guess_ab_col(),
                           NEO = guess_ab_col(),
                           NET = guess_ab_col(),
                           NIT = guess_ab_col(),
                           NOR = guess_ab_col(),
                           NOV = guess_ab_col(),
                           OFX = guess_ab_col(),
                           OXA = guess_ab_col(),
                           PEN = guess_ab_col(),
                           PIP = guess_ab_col(),
                           TZP = guess_ab_col(),
                           PLB = guess_ab_col(),
                           PRI = guess_ab_col(),
                           QDA = guess_ab_col(),
                           RIF = guess_ab_col(),
                           RXT = guess_ab_col(),
                           SIS = guess_ab_col(),
                           TEC = guess_ab_col(),
                           TCY = guess_ab_col(),
                           TIC = guess_ab_col(),
                           TGC = guess_ab_col(),
                           TOB = guess_ab_col(),
                           TMP = guess_ab_col(),
                           SXT = guess_ab_col(),
                           VAN = guess_ab_col()) {
  # check columns
  if (identical(AMC, as.name("guess_ab_col"))) AMC <- guess_ab_col(tbl, "AMC", verbose = verbose)
  if (identical(AMK, as.name("guess_ab_col"))) AMK <- guess_ab_col(tbl, "AMK", verbose = verbose)
  if (identical(AMX, as.name("guess_ab_col"))) AMX <- guess_ab_col(tbl, "AMX", verbose = verbose)
  if (identical(AMP, as.name("guess_ab_col"))) AMP <- guess_ab_col(tbl, "AMP", verbose = verbose)
  if (identical(AZM, as.name("guess_ab_col"))) AZM <- guess_ab_col(tbl, "AZM", verbose = verbose)
  if (identical(AZL, as.name("guess_ab_col"))) AZL <- guess_ab_col(tbl, "AZL", verbose = verbose)
  if (identical(ATM, as.name("guess_ab_col"))) ATM <- guess_ab_col(tbl, "ATM", verbose = verbose)
  if (identical(RID, as.name("guess_ab_col"))) RID <- guess_ab_col(tbl, "RID", verbose = verbose)
  if (identical(FEP, as.name("guess_ab_col"))) FEP <- guess_ab_col(tbl, "FEP", verbose = verbose)
  if (identical(CTX, as.name("guess_ab_col"))) CTX <- guess_ab_col(tbl, "CTX", verbose = verbose)
  if (identical(FOX, as.name("guess_ab_col"))) FOX <- guess_ab_col(tbl, "FOX", verbose = verbose)
  if (identical(CED, as.name("guess_ab_col"))) CED <- guess_ab_col(tbl, "CED", verbose = verbose)
  if (identical(CAZ, as.name("guess_ab_col"))) CAZ <- guess_ab_col(tbl, "CAZ", verbose = verbose)
  if (identical(CRO, as.name("guess_ab_col"))) CRO <- guess_ab_col(tbl, "CRO", verbose = verbose)
  if (identical(CXM, as.name("guess_ab_col"))) CXM <- guess_ab_col(tbl, "CXM", verbose = verbose)
  if (identical(CHL, as.name("guess_ab_col"))) CHL <- guess_ab_col(tbl, "CHL", verbose = verbose)
  if (identical(CIP, as.name("guess_ab_col"))) CIP <- guess_ab_col(tbl, "CIP", verbose = verbose)
  if (identical(CLR, as.name("guess_ab_col"))) CLR <- guess_ab_col(tbl, "CLR", verbose = verbose)
  if (identical(CLI, as.name("guess_ab_col"))) CLI <- guess_ab_col(tbl, "CLI", verbose = verbose)
  if (identical(FLC, as.name("guess_ab_col"))) FLC <- guess_ab_col(tbl, "FLC", verbose = verbose)
  if (identical(COL, as.name("guess_ab_col"))) COL <- guess_ab_col(tbl, "COL", verbose = verbose)
  if (identical(CZO, as.name("guess_ab_col"))) CZO <- guess_ab_col(tbl, "CZO", verbose = verbose)
  if (identical(DAP, as.name("guess_ab_col"))) DAP <- guess_ab_col(tbl, "DAP", verbose = verbose)
  if (identical(DOX, as.name("guess_ab_col"))) DOX <- guess_ab_col(tbl, "DOX", verbose = verbose)
  if (identical(ETP, as.name("guess_ab_col"))) ETP <- guess_ab_col(tbl, "ETP", verbose = verbose)
  if (identical(ERY, as.name("guess_ab_col"))) ERY <- guess_ab_col(tbl, "ERY", verbose = verbose)
  if (identical(FOS, as.name("guess_ab_col"))) FOS <- guess_ab_col(tbl, "FOS", verbose = verbose)
  if (identical(FUS, as.name("guess_ab_col"))) FUS <- guess_ab_col(tbl, "FUS", verbose = verbose)
  if (identical(GEN, as.name("guess_ab_col"))) GEN <- guess_ab_col(tbl, "GEN", verbose = verbose)
  if (identical(IPM, as.name("guess_ab_col"))) IPM <- guess_ab_col(tbl, "IPM", verbose = verbose)
  if (identical(KAN, as.name("guess_ab_col"))) KAN <- guess_ab_col(tbl, "KAN", verbose = verbose)
  if (identical(LVX, as.name("guess_ab_col"))) LVX <- guess_ab_col(tbl, "LVX", verbose = verbose)
  if (identical(LIN, as.name("guess_ab_col"))) LIN <- guess_ab_col(tbl, "LIN", verbose = verbose)
  if (identical(LNZ, as.name("guess_ab_col"))) LNZ <- guess_ab_col(tbl, "LNZ", verbose = verbose)
  if (identical(MEM, as.name("guess_ab_col"))) MEM <- guess_ab_col(tbl, "MEM", verbose = verbose)
  if (identical(MTR, as.name("guess_ab_col"))) MTR <- guess_ab_col(tbl, "MTR", verbose = verbose)
  if (identical(MEZ, as.name("guess_ab_col"))) MEZ <- guess_ab_col(tbl, "MEZ", verbose = verbose)
  if (identical(MNO, as.name("guess_ab_col"))) MNO <- guess_ab_col(tbl, "MNO", verbose = verbose)
  if (identical(MFX, as.name("guess_ab_col"))) MFX <- guess_ab_col(tbl, "MFX", verbose = verbose)
  if (identical(NAL, as.name("guess_ab_col"))) NAL <- guess_ab_col(tbl, "NAL", verbose = verbose)
  if (identical(NEO, as.name("guess_ab_col"))) NEO <- guess_ab_col(tbl, "NEO", verbose = verbose)
  if (identical(NET, as.name("guess_ab_col"))) NET <- guess_ab_col(tbl, "NET", verbose = verbose)
  if (identical(NIT, as.name("guess_ab_col"))) NIT <- guess_ab_col(tbl, "NIT", verbose = verbose)
  if (identical(NOR, as.name("guess_ab_col"))) NOR <- guess_ab_col(tbl, "NOR", verbose = verbose)
  if (identical(NOV, as.name("guess_ab_col"))) NOV <- guess_ab_col(tbl, "NOV", verbose = verbose)
  if (identical(OFX, as.name("guess_ab_col"))) OFX <- guess_ab_col(tbl, "OFX", verbose = verbose)
  if (identical(OXA, as.name("guess_ab_col"))) OXA <- guess_ab_col(tbl, "OXA", verbose = verbose)
  if (identical(PEN, as.name("guess_ab_col"))) PEN <- guess_ab_col(tbl, "PEN", verbose = verbose)
  if (identical(PIP, as.name("guess_ab_col"))) PIP <- guess_ab_col(tbl, "PIP", verbose = verbose)
  if (identical(TZP, as.name("guess_ab_col"))) TZP <- guess_ab_col(tbl, "TZP", verbose = verbose)
  if (identical(PLB, as.name("guess_ab_col"))) PLB <- guess_ab_col(tbl, "PLB", verbose = verbose)
  if (identical(PRI, as.name("guess_ab_col"))) PRI <- guess_ab_col(tbl, "PRI", verbose = verbose)
  if (identical(QDA, as.name("guess_ab_col"))) QDA <- guess_ab_col(tbl, "QDA", verbose = verbose)
  if (identical(RIF, as.name("guess_ab_col"))) RIF <- guess_ab_col(tbl, "RIF", verbose = verbose)
  if (identical(RXT, as.name("guess_ab_col"))) RXT <- guess_ab_col(tbl, "RXT", verbose = verbose)
  if (identical(SIS, as.name("guess_ab_col"))) SIS <- guess_ab_col(tbl, "SIS", verbose = verbose)
  if (identical(TEC, as.name("guess_ab_col"))) TEC <- guess_ab_col(tbl, "TEC", verbose = verbose)
  if (identical(TCY, as.name("guess_ab_col"))) TCY <- guess_ab_col(tbl, "TCY", verbose = verbose)
  if (identical(TIC, as.name("guess_ab_col"))) TIC <- guess_ab_col(tbl, "TIC", verbose = verbose)
  if (identical(TGC, as.name("guess_ab_col"))) TGC <- guess_ab_col(tbl, "TGC", verbose = verbose)
  if (identical(TOB, as.name("guess_ab_col"))) TOB <- guess_ab_col(tbl, "TOB", verbose = verbose)
  if (identical(TMP, as.name("guess_ab_col"))) TMP <- guess_ab_col(tbl, "TMP", verbose = verbose)
  if (identical(SXT, as.name("guess_ab_col"))) SXT <- guess_ab_col(tbl, "SXT", verbose = verbose)
  if (identical(VAN, as.name("guess_ab_col"))) VAN <- guess_ab_col(tbl, "VAN", verbose = verbose)
  columns_available <- c(AMC = AMC, AMK = AMK, AMX = AMX, AMP = AMP, AZM = AZM,
                         AZL = AZL, ATM = ATM, RID = RID, FEP = FEP, CTX = CTX,
                         FOX = FOX, CED = CED, CAZ = CAZ, CRO = CRO, CXM = CXM,
                         CHL = CHL, CIP = CIP, CLR = CLR, CLI = CLI, FLC = FLC,
                         COL = COL, CZO = CZO, DAP = DAP, DOX = DOX, ETP = ETP,
                         ERY = ERY, FOS = FOS, FUS = FUS, GEN = GEN, IPM = IPM,
                         KAN = KAN, LVX = LVX, LIN = LIN, LNZ = LNZ, MEM = MEM,
                         MTR = MTR, MEZ = MEZ, MNO = MNO, MFX = MFX, NAL = NAL,
                         NEO = NEO, NET = NET, NIT = NIT, NOR = NOR, NOV = NOV,
                         OFX = OFX, OXA = OXA, PEN = PEN, PIP = PIP, TZP = TZP,
                         PLB = PLB, PRI = PRI, QDA = QDA, RIF = RIF, RXT = RXT,
                         SIS = SIS, TEC = TEC, TCY = TCY, TIC = TIC, TGC = TGC,
                         TOB = TOB, TMP = TMP, SXT = SXT, VAN = VAN)

  if (!is.null(hard_dependencies)) {
    if (!all(hard_dependencies %in% names(columns_available[!is.na(columns_available)]))) {
      # missing a hard dependency will return NA and consequently the data will not be analysed
      missing <- hard_dependencies[!hard_dependencies %in% names(columns_available[!is.na(columns_available)])]
      generate_warning_abs_missing(missing, any = FALSE)
      return(NA)
    }
  }
  if (!is.null(soft_dependencies)) {
    if (!all(soft_dependencies %in% names(columns_available[!is.na(columns_available)]))) {
      # missing a soft dependency may lower the reliability
      missing <- soft_dependencies[!soft_dependencies %in% names(columns_available[!is.na(columns_available)])]
      missing <- paste0("`", missing, "` (", ab_name(missing, tolower = TRUE), ")")
      warning('Reliability might be improved if these antimicrobial results would be available too: ', paste(missing, collapse = ", "),
              immediate. = TRUE,
              call. = FALSE)
    }
  }
  #deps <- c(soft_dependencies, hard_dependencies)
  #if (length(deps) > 0) {
  #  columns_available[names(columns_available) %in% deps]
  #} else {
    columns_available
  #}
}

generate_warning_abs_missing <- function(missing, any = FALSE) {
  missing <- paste0("`", missing, "` (", ab_name(missing, tolower = TRUE), ")")
  if (any == TRUE) {
    any_txt <- c(" any of", "is")
  } else {
    any_txt <- c("", "are")
  }
  warning(paste0("Introducing NAs since", any_txt[1], " these antimicrobials ", any_txt[2], " required: ",
                 paste(missing, collapse = ", ")),
          immediate. = TRUE,
          call. = FALSE)
}


stopifnot_installed_package <- function(package) {
  if (!package %in% base::rownames(utils::installed.packages())) {
    stop("this function requires the ", package, " package.", call. = FALSE)
  }
}

# translate strings based on inst/translations.tsv
#' @importFrom dplyr %>% filter
t <- function(from, language = get_locale()) {
  # if (getOption("AMR_locale", "en") != language) {
  #   language <- getOption("AMR_locale", "en")
  # }

  if (is.null(language)) {
    return(from)
  }
  if (language %in% c("en", "")) {
    return(from)
  }

  df_trans <- utils::read.table(file = system.file("translations.tsv", package = "AMR"),
                                sep = "\t",
                                stringsAsFactors = FALSE,
                                header = TRUE,
                                blank.lines.skip = TRUE,
                                fill = TRUE,
                                strip.white = TRUE,
                                encoding = "UTF-8",
                                fileEncoding = "UTF-8",
                                na.strings = c(NA, "", NULL))

  if (!language %in% df_trans$lang) {
    stop("Unsupported language: '", language, "' - use one of: ",
         paste0("'", sort(unique(df_trans$lang)), "'", collapse = ", "),
         call. = FALSE)
  }

  df_trans <- df_trans %>% filter(lang == language)

  # default case sensitive if value if 'ignore.case' is missing:
  df_trans$ignore.case[is.na(df_trans$ignore.case)] <- FALSE
  # default not using regular expressions (fixed = TRUE) if 'fixed' is missing:
  df_trans$fixed[is.na(df_trans$fixed)] <- TRUE

  # check if text to look for is in one of the patterns
  any_form_in_patterns <- tryCatch(any(from %like% paste0("(", paste(df_trans$pattern, collapse = "|"), ")")),
                                   error = function(e) {
                                     warning("Translation not possible. Please open an issue on GitLab (https://gitlab.com/msberends/AMR/issues) or GitHub (https://github.com/msberends/AMR/issues).", call. = FALSE)
                                     return(FALSE)
                                   })
  if (NROW(df_trans) == 0 | !any_form_in_patterns) {
    return(from)
  }

  for (i in 1:nrow(df_trans)) {
    from <- gsub(x = from,
                 pattern = df_trans$pattern[i],
                 replacement = df_trans$replacement[i],
                 fixed = df_trans$fixed[i],
                 ignore.case = df_trans$ignore.case[i])
  }

  # force UTF-8 for diacritics
  base::enc2utf8(from)

}
