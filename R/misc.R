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

get_column_abx <- function(x,
                           soft_dependencies = NULL,
                           hard_dependencies = NULL,
                           verbose = FALSE,
                           ...) {

  # determine from given data set
  df_trans <- data.frame(colnames = colnames(x),
                         abcode = suppressWarnings(as.ab(colnames(x))))
  df_trans <- df_trans[!is.na(df_trans$abcode),]
  x <- as.character(df_trans$colnames)
  names(x) <- df_trans$abcode

  # add from self-defined dots (...):
  # get_column_abx(septic_patients %>% rename(thisone = AMX), amox = "thisone")
  dots <- list(...)
  if (length(dots) > 0) {
    dots <- unlist(dots)
    newnames <- suppressWarnings(as.ab(names(dots)))
    if (any(is.na(newnames))) {
      warning("Invalid antibiotic reference(s): ", toString(names(dots)[is.na(newnames)]),
              call. = FALSE, immediate. = TRUE)
    }
    names(dots) <- newnames
    dots <- dots[!is.na(names(dots))]
    # merge, but overwrite automatically determined ones by 'dots'
    x <- c(x[!x %in% dots & !names(x) %in% names(dots)], dots)
  }

  # sort on name
  x <- x[sort(names(x))]
  duplies <- x[base::duplicated(x)]

  if (verbose == TRUE) {
    for (i in 1:length(x)) {
      if (x[i] %in% duplies) {
        message(red(paste0("NOTE: Using column `", bold(x[i]), "` as input for ", names(x)[i],
                           " (", ab_name(names(x)[i], language = "en", tolower = TRUE), ") [DUPLICATED USE].")))
      } else {
        message(blue(paste0("NOTE: Using column `", bold(x[i]), "` as input for ", names(x)[i],
                            " (", ab_name(names(x)[i], language = "en", tolower = TRUE), ").")))
      }
    }
  }

  if (n_distinct(x) != length(x)) {
    msg_txt <- paste("Column(s)", paste0("'", duplies, "'", collapse = "'"), "used for more than one antibiotic.")
    if (verbose == FALSE) {
      msg_txt <- paste(msg_txt, "Use verbose = TRUE to see which antibiotics are used by which columns.")
    }
    stop(msg_txt, call. = FALSE)
  }

  if (!is.null(hard_dependencies)) {
    if (!all(hard_dependencies %in% names(x))) {
      # missing a hard dependency will return NA and consequently the data will not be analysed
      missing <- hard_dependencies[!hard_dependencies %in% names(x)]
      generate_warning_abs_missing(missing, any = FALSE)
      return(NA)
    }
  }
  if (!is.null(soft_dependencies)) {
    if (!all(soft_dependencies %in% names(x))) {
      # missing a soft dependency may lower the reliability
      missing <- soft_dependencies[!soft_dependencies %in% names(x)]
      missing <- paste0("`", missing, "` (", ab_name(missing, tolower = TRUE), ")")
      warning('Reliability might be improved if these antimicrobial results would be available too: ', paste(missing, collapse = ", "),
              immediate. = TRUE,
              call. = FALSE)
    }
  }
  x
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

"%or%" <- function(x, y) {
  ifelse(!is.na(x), x, ifelse(!is.na(y), y, NA))
}
