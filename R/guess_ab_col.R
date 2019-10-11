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

#' Guess antibiotic column
#'
#' This tries to find a column name in a data set based on information from the \code{\link{antibiotics}} data set. Also supports WHONET abbreviations.
#' @param x a \code{data.frame}
#' @param search_string a text to search \code{x} for, will be checked with \code{\link{as.ab}} if this value is not a column in \code{x}
#' @param verbose a logical to indicate whether additional info should be printed
#' @details You can look for an antibiotic (trade) name or abbreviation and it will search \code{x} and the \code{\link{antibiotics}} data set for any column containing a name or code of that antibiotic. \strong{Longer columns names take precendence over shorter column names.}
#' @importFrom dplyr %>% select filter_all any_vars
#' @importFrom crayon blue
#' @return A column name of \code{x}, or \code{NULL} when no result is found.
#' @export
#' @inheritSection AMR Read more on our website!
#' @examples
#' df <- data.frame(amox = "S",
#'                  tetr = "R")
#'
#' guess_ab_col(df, "amoxicillin")
#' # [1] "amox"
#' guess_ab_col(df, "J01AA07") # ATC code of tetracycline
#' # [1] "tetr"
#'
#' guess_ab_col(df, "J01AA07", verbose = TRUE)
#' # Note: Using column `tetr` as input for "J01AA07".
#' # [1] "tetr"
#'
#' # WHONET codes
#' df <- data.frame(AMP_ND10 = "R",
#'                  AMC_ED20 = "S")
#' guess_ab_col(df, "ampicillin")
#' # [1] "AMP_ND10"
#' guess_ab_col(df, "J01CR02")
#' # [1] "AMC_ED20"
#' guess_ab_col(df, as.ab("augmentin"))
#' # [1] "AMC_ED20"
#'
#' # Longer names take precendence:
#' df <- data.frame(AMP_ED2 = "S",
#'                  AMP_ED20 = "S")
#' guess_ab_col(df, "ampicillin")
#' # [1] "AMP_ED20"
guess_ab_col <- function(x = NULL, search_string = NULL, verbose = FALSE) {
  if (is.null(x) & is.null(search_string)) {
    return(as.name("guess_ab_col"))
  }
  if (!is.data.frame(x)) {
    stop("`x` must be a data.frame")
  }

  if (length(search_string) > 1) {
    warning("argument 'search_string' has length > 1 and only the first element will be used")
    search_string <- search_string[1]
  }
  search_string <- as.character(search_string)

  if (search_string %in% colnames(x)) {
    ab_result <- search_string
  } else {
    search_string.ab <- suppressWarnings(as.ab(search_string))
    if (search_string.ab %in% colnames(x)) {
      ab_result <- colnames(x)[colnames(x) == search_string.ab][1L]

    } else if (any(tolower(colnames(x)) %in% tolower(unlist(ab_property(search_string.ab, "abbreviations"))))) {
      ab_result <- colnames(x)[tolower(colnames(x)) %in% tolower(unlist(ab_property(search_string.ab, "abbreviations")))][1L]

    } else {
      # sort colnames on length - longest first
      cols <- colnames(x[, x %>% colnames() %>% nchar() %>% order() %>% rev()])
      df_trans <- data.frame(cols = cols,
                             abs = suppressWarnings(as.ab(cols)),
                             stringsAsFactors = FALSE)
      ab_result <- df_trans[which(df_trans$abs == search_string.ab), "cols"]
      ab_result <- ab_result[!is.na(ab_result)][1L]
    }
  }

  if (length(ab_result) == 0) {
    if (verbose == TRUE) {
      message(paste0("No column found as input for `", search_string,
                     "` (", ab_name(search_string, language = "en", tolower = TRUE), ")."))
    }
    return(NULL)
  } else {
    if (verbose == TRUE) {
      message(blue(paste0("NOTE: Using column `", bold(ab_result), "` as input for `", search_string,
                          "` (", ab_name(search_string, language = "en", tolower = TRUE), ").")))
    }
    return(ab_result)
  }
}


#' @importFrom crayon blue bold
#' @importFrom dplyr %>% mutate arrange pull
get_column_abx <- function(x,
                           soft_dependencies = NULL,
                           hard_dependencies = NULL,
                           verbose = FALSE,
                           ...) {

  message(blue("NOTE: Auto-guessing columns suitable for analysis..."))
  
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  x_bak <- x
  # only check columns that are a valid AB code, ATC code, name, abbreviation or synonym,
  # or already have the rsi class (as.rsi) 
  # and that have no more than 50% invalid values
  vectr_antibiotics <- unique(toupper(unlist(AMR::antibiotics[, c("ab", "atc", "name", "abbreviations", "synonyms")])))
  vectr_antibiotics <- vectr_antibiotics[!is.na(vectr_antibiotics) & nchar(vectr_antibiotics) >= 3]
  x_columns <- sapply(colnames(x), function(col, df = x_bak) {
    if (toupper(col) %in% vectr_antibiotics | 
        is.rsi(as.data.frame(df)[, col]) |
        is.rsi.eligible(as.data.frame(df)[, col], threshold = 0.5)) {
      return(col)
    } else {
      return(NA_character_)
    }
  })
  x_columns <- x_columns[!is.na(x_columns)]
  x <- x[, x_columns, drop = FALSE] # without drop = TRUE, x will become a vector when x_columns is length 1

  df_trans <- data.frame(colnames = colnames(x),
                         abcode = suppressWarnings(as.ab(colnames(x))))
  df_trans <- df_trans[!is.na(df_trans$abcode), ]
  x <- as.character(df_trans$colnames)
  names(x) <- df_trans$abcode

  # add from self-defined dots (...):
  # such as get_column_abx(example_isolates %>% rename(thisone = AMX), amox = "thisone")
  dots <- list(...)
  if (length(dots) > 0) {
    newnames <- suppressWarnings(as.ab(names(dots)))
    if (any(is.na(newnames))) {
      warning("Invalid antibiotic reference(s): ", toString(names(dots)[is.na(newnames)]),
              call. = FALSE, immediate. = TRUE)
    }
    # turn all NULLs to NAs
    dots <- unlist(lapply(dots, function(x) if (is.null(x)) NA else x))
    names(dots) <- newnames
    dots <- dots[!is.na(names(dots))]
    # merge, but overwrite automatically determined ones by 'dots'
    x <- c(x[!x %in% dots & !names(x) %in% names(dots)], dots)
    # delete NAs, this will make e.g. eucast_rules(... TMP = NULL) work to prevent TMP from being used
    x <- x[!is.na(x)]
  }

  # sort on name
  x <- x[order(names(x), x)]
  duplicates <- x[base::duplicated(x)]
  x <- x[!names(x) %in% names(duplicates)]
  
  if (verbose == TRUE) {
    for (i in seq_len(length(x))) {
      message(blue(paste0("NOTE: Using column `", bold(x[i]), "` as input for `", names(x)[i],
                          "` (", ab_name(names(x)[i], tolower = TRUE), ").")))
    }
  } else if (length(duplicates) > 0) {
    for (i in seq_len(length(duplicates))) {
      warning(red(paste0("Using column `", bold(duplicates[i]), "` as input for `", names(x[which(x == duplicates[i])]), 
                         "` (", ab_name(names(x[names(which(x == duplicates))[i]]), tolower = TRUE), 
                         "), although it was matched for multiple antibiotics or columns.")), call. = FALSE)
    }
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
      missing_txt <- data.frame(missing = missing,
                                missing_names = AMR::ab_name(missing, tolower = TRUE),
                                stringsAsFactors = FALSE) %>%
        mutate(txt = paste0(bold(missing), " (", missing_names, ")")) %>%
        arrange(missing_names) %>%
        pull(txt)
      message(blue("NOTE: Reliability might be improved if these antimicrobial results would be available too:",
                   paste(missing_txt, collapse = ", ")))
    }
  }
  x
}

generate_warning_abs_missing <- function(missing, any = FALSE) {
  missing <- paste0(missing, " (", ab_name(missing, tolower = TRUE), ")")
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
