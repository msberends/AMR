# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Data Analysis for R                   #
#                                                                      #
# SOURCE                                                               #
# https://github.com/msberends/AMR                                     #
#                                                                      #
# LICENCE                                                              #
# (c) 2018-2022 Berends MS, Luz CF et al.                              #
# Developed at the University of Groningen, the Netherlands, in        #
# collaboration with non-profit organisations Certe Medical            #
# Diagnostics & Advice, and University Medical Center Groningen.       # 
#                                                                      #
# This R package is free software; you can freely use and distribute   #
# it for both personal and commercial purposes under the terms of the  #
# GNU General Public License version 2.0 (GNU GPL-2), as published by  #
# the Free Software Foundation.                                        #
# We created this package for both routine data analysis and academic  #
# research and it was publicly released in the hope that it will be    #
# useful, but it comes WITHOUT ANY WARRANTY OR LIABILITY.              #
#                                                                      #
# Visit our website for the full manual and a complete tutorial about  #
# how to conduct AMR data analysis: https://msberends.github.io/AMR/   #
# ==================================================================== #

#' Guess Antibiotic Column
#'
#' This tries to find a column name in a data set based on information from the [antibiotics] data set. Also supports WHONET abbreviations.
#' @inheritSection lifecycle Stable Lifecycle
#' @param x a [data.frame]
#' @param search_string a text to search `x` for, will be checked with [as.ab()] if this value is not a column in `x`
#' @param verbose a [logical] to indicate whether additional info should be printed
#' @param only_rsi_columns a [logical] to indicate whether only antibiotic columns must be detected that were transformed to class `<rsi>` (see [as.rsi()]) on beforehand (defaults to `FALSE`)
#' @details You can look for an antibiotic (trade) name or abbreviation and it will search `x` and the [antibiotics] data set for any column containing a name or code of that antibiotic. **Longer columns names take precedence over shorter column names.**
#' @return A column name of `x`, or `NULL` when no result is found.
#' @export
#' @inheritSection AMR Read more on Our Website!
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
#' # NOTE: Using column 'tetr' as input for J01AA07 (tetracycline).
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
guess_ab_col <- function(x = NULL, search_string = NULL, verbose = FALSE, only_rsi_columns = FALSE) {
  meet_criteria(x, allow_class = "data.frame", allow_NULL = TRUE)
  meet_criteria(search_string, allow_class = "character", has_length = 1, allow_NULL = TRUE)
  meet_criteria(verbose, allow_class = "logical", has_length = 1)
  meet_criteria(only_rsi_columns, allow_class = "logical", has_length = 1)
  
  if (is.null(x) & is.null(search_string)) {
    return(as.name("guess_ab_col"))
  } else {
    meet_criteria(search_string, allow_class = "character", has_length = 1, allow_NULL = FALSE)
  }
  
  all_found <- get_column_abx(x, info = verbose, only_rsi_columns = only_rsi_columns,
                              verbose = verbose, fn = "guess_ab_col")
  search_string.ab <- suppressWarnings(as.ab(search_string))
  ab_result <- unname(all_found[names(all_found) == search_string.ab])
  
  if (length(ab_result) == 0) {
    if (verbose == TRUE) {
      message_("No column found as input for ", search_string,
               " (", ab_name(search_string, language = NULL, tolower = TRUE), ").",
               add_fn = font_black,
               as_note = FALSE)
    }
    return(NULL)
  } else {
    if (verbose == TRUE) {
      message_("Using column '", font_bold(ab_result), "' as input for ", search_string,
               " (", ab_name(search_string, language = NULL, tolower = TRUE), ").")
    }
    return(ab_result)
  }
}

get_column_abx <- function(x,
                           ...,
                           soft_dependencies = NULL,
                           hard_dependencies = NULL,
                           verbose = FALSE,
                           info = TRUE,
                           only_rsi_columns = FALSE,
                           sort = TRUE,
                           reuse_previous_result = TRUE,
                           fn = NULL) {
  # check if retrieved before, then get it from package environment
  if (isTRUE(reuse_previous_result) && identical(unique_call_id(entire_session = FALSE,
                                                                match_fn = fn),
                                                 pkg_env$get_column_abx.call)) {
    # so within the same call, within the same environment, we got here again.
    # but we could've come from another function within the same call, so now only check the columns that changed
    
    # first remove the columns that are not existing anymore
    previous <- pkg_env$get_column_abx.out
    current <- previous[previous %in% colnames(x)]
    
    # then compare columns in current call with columns in original call
    new_cols <- colnames(x)[!colnames(x) %in% pkg_env$get_column_abx.checked_cols]
    if (length(new_cols) > 0) {
      # these columns did not exist in the last call, so add them
      new_cols_rsi <- get_column_abx(x[, new_cols, drop = FALSE], reuse_previous_result = FALSE, info = FALSE, sort = FALSE)
      current <- c(current, new_cols_rsi)
      # order according to columns in current call
      current <- current[match(colnames(x)[colnames(x) %in% current], current)]
    }
    
    # update pkg environment to improve speed on next run
    pkg_env$get_column_abx.out <- current
    pkg_env$get_column_abx.checked_cols <- colnames(x)

    # and return right values
    return(pkg_env$get_column_abx.out)
  }
  
  meet_criteria(x, allow_class = "data.frame")
  meet_criteria(soft_dependencies, allow_class = "character", allow_NULL = TRUE)
  meet_criteria(hard_dependencies, allow_class = "character", allow_NULL = TRUE)
  meet_criteria(verbose, allow_class = "logical", has_length = 1)
  meet_criteria(info, allow_class = "logical", has_length = 1)
  meet_criteria(only_rsi_columns, allow_class = "logical", has_length = 1)
  meet_criteria(sort, allow_class = "logical", has_length = 1)
  
  if (info == TRUE) {
    message_("Auto-guessing columns suitable for analysis", appendLF = FALSE, as_note = FALSE)
  }
  
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  x.bak <- x
  if (only_rsi_columns == TRUE) {
    x <- x[, which(is.rsi(x)), drop = FALSE]
  }

  if (NROW(x) > 10000) {
    # only test maximum of 10,000 values per column
    if (info == TRUE) {
      message_(" (using only ", font_bold("the first 10,000 rows"), ")...",
               appendLF = FALSE, 
               as_note = FALSE)
    }
    x <- x[1:10000, , drop = FALSE]
  } else if (info == TRUE) {
    message_("...", appendLF = FALSE, as_note = FALSE)
  }

  # only check columns that are a valid AB code, ATC code, name, abbreviation or synonym,
  # or already have the <rsi> class (as.rsi) 
  # and that they have no more than 50% invalid values
  vectr_antibiotics <- unlist(AB_lookup$generalised_all)
  vectr_antibiotics <- vectr_antibiotics[!is.na(vectr_antibiotics) & nchar(vectr_antibiotics) >= 3]
  x_columns <- vapply(FUN.VALUE = character(1), 
                      colnames(x),
                      function(col, df = x) {
                        if (generalise_antibiotic_name(col) %in% vectr_antibiotics || 
                            is.rsi(x[, col, drop = TRUE]) ||
                            is.rsi.eligible(x[, col, drop = TRUE], threshold = 0.5)
                        ) {
                          return(col)
                        } else {
                          return(NA_character_)
                        }
                      }, USE.NAMES = FALSE)
  
  x_columns <- x_columns[!is.na(x_columns)]
  x <- x[, x_columns, drop = FALSE] # without drop = FALSE, x will become a vector when x_columns is length 1
  df_trans <- data.frame(colnames = colnames(x),
                         abcode = suppressWarnings(as.ab(colnames(x), info = FALSE)),
                         stringsAsFactors = FALSE)
  df_trans <- df_trans[!is.na(df_trans$abcode), , drop = FALSE]
  out <- as.character(df_trans$colnames)
  names(out) <- df_trans$abcode
  
  # add from self-defined dots (...):
  # such as get_column_abx(example_isolates %>% rename(thisone = AMX), amox = "thisone")
  all_okay <- TRUE
  dots <- list(...)
  # remove data.frames, since this is also used running `eucast_rules(eucast_rules_df = df)`
  dots <- dots[!vapply(FUN.VALUE = logical(1), dots, is.data.frame)]
  if (length(dots) > 0) {
    newnames <- suppressWarnings(as.ab(names(dots), info = FALSE))
    if (any(is.na(newnames))) {
      if (info == TRUE) {
        message_(" WARNING", add_fn = list(font_yellow, font_bold), as_note = FALSE)
      }
      warning_("Invalid antibiotic reference(s): ", vector_and(names(dots)[is.na(newnames)], quotes = FALSE),
               call = FALSE,
               immediate = TRUE)
      all_okay <- FALSE
    }
    unexisting_cols <- which(!vapply(FUN.VALUE = logical(1), dots, function(col) all(col %in% x_columns)))
    if (length(unexisting_cols) > 0) {
      if (info == TRUE) {
        message_(" ERROR", add_fn = list(font_red, font_bold), as_note = FALSE)
      }
      stop_("Column(s) not found: ", vector_and(unlist(dots[[unexisting_cols]]), quotes = FALSE),
            call = FALSE)
      all_okay <- FALSE
    }
    # turn all NULLs to NAs
    dots <- unlist(lapply(dots, function(dot) if (is.null(dot)) NA else dot))
    names(dots) <- newnames
    dots <- dots[!is.na(names(dots))]
    # merge, but overwrite automatically determined ones by 'dots'
    out <- c(out[!out %in% dots & !names(out) %in% names(dots)], dots)
    # delete NAs, this will make e.g. eucast_rules(... TMP = NULL) work to prevent TMP from being used
    out <- out[!is.na(out)]
  }
  
  if (length(out) == 0) {
    if (info == TRUE & all_okay == TRUE) {
      message_("No columns found.")
    }
    pkg_env$get_column_abx.call <- unique_call_id(entire_session = FALSE, match_fn = fn)
    pkg_env$get_column_abx.checked_cols <- colnames(x.bak)
    pkg_env$get_column_abx.out <- out
    return(out)
  }
  
  # sort on name
  if (sort == TRUE) {
    out <- out[order(names(out), out)]
  }
  # only keep the first hits, no duplicates
  duplicates <- c(out[duplicated(names(out))], out[duplicated(unname(out))])
  if (length(duplicates) > 0) {
    all_okay <- FALSE
  }
  
  if (info == TRUE) {
    if (all_okay == TRUE) {
      message_(" OK.", add_fn = list(font_green, font_bold), as_note = FALSE)
    } else {
      message_(" WARNING.", add_fn = list(font_yellow, font_bold), as_note = FALSE)
    }
    for (i in seq_len(length(out))) {
      if (verbose == TRUE & !names(out[i]) %in% names(duplicates)) {
        message_("Using column '", font_bold(out[i]), "' as input for ", names(out)[i],
                 " (", ab_name(names(out)[i], tolower = TRUE, language = NULL), ").")
      }
      if (names(out[i]) %in% names(duplicates)) {
        already_set_as <- out[unname(out) == unname(out[i])][1L]
        warning_(paste0("Column '", font_bold(out[i]), "' will not be used for ", 
                        names(out)[i], " (", ab_name(names(out)[i], tolower = TRUE, language = NULL), ")",
                        ", as it is already set for ", 
                        names(already_set_as), " (", ab_name(names(already_set_as), tolower = TRUE, language = NULL), ")"),
                 add_fn = font_red,
                 immediate = verbose)
      }
    }
  }
  
  out <- out[!duplicated(names(out))]
  out <- out[!duplicated(unname(out))]
  if (sort == TRUE) {
    out <- out[order(names(out), out)]
  }
  
  if (!is.null(hard_dependencies)) {
    hard_dependencies <- unique(hard_dependencies)
    if (!all(hard_dependencies %in% names(out))) {
      # missing a hard dependency will return NA and consequently the data will not be analysed
      missing <- hard_dependencies[!hard_dependencies %in% names(out)]
      generate_warning_abs_missing(missing, any = FALSE)
      return(NA)
    }
  }
  if (!is.null(soft_dependencies)) {
    soft_dependencies <- unique(soft_dependencies)
    if (info == TRUE & !all(soft_dependencies %in% names(out))) {
      # missing a soft dependency may lower the reliability
      missing <- soft_dependencies[!soft_dependencies %in% names(out)]
      missing_msg <- vector_and(paste0(ab_name(missing, tolower = TRUE, language = NULL), 
                                       " (", font_bold(missing, collapse = NULL), ")"), 
                                quotes = FALSE)
      message_("Reliability would be improved if these antimicrobial results would be available too: ",
               missing_msg)
    }
  }
  
  pkg_env$get_column_abx.call <- unique_call_id(entire_session = FALSE, match_fn = fn)
  pkg_env$get_column_abx.checked_cols <- colnames(x.bak)
  pkg_env$get_column_abx.out <- out
  out
}

get_ab_from_namespace <- function(x, cols_ab) {
  # cols_ab comes from get_column_abx()
  
  x <- trimws(unique(toupper(unlist(strsplit(x, ",")))))
  x_new <- character()
  for (val in x) {
    if (paste0("AB_", val) %in% ls(envir = asNamespace("AMR"))) {
      # antibiotic group names, as defined in data-raw/_internals.R, such as `AB_CARBAPENEMS`
      val <- eval(parse(text = paste0("AB_", val)), envir = asNamespace("AMR"))
    } else if (val %in% AB_lookup$ab) {
      # separate drugs, such as `AMX`
      val <- as.ab(val)
    } else {
      stop_("unknown antimicrobial agent (group): ", val, call = FALSE)
    }
    x_new <- c(x_new, val)
  }
  x_new <- unique(x_new)
  out <- cols_ab[match(x_new, names(cols_ab))]
  out[!is.na(out)]
}

generate_warning_abs_missing <- function(missing, any = FALSE) {
  missing <- paste0(missing, " (", ab_name(missing, tolower = TRUE, language = NULL), ")")
  if (any == TRUE) {
    any_txt <- c(" any of", "is")
  } else {
    any_txt <- c("", "are")
  }
  warning_(paste0("Introducing NAs since", any_txt[1], " these antimicrobials ", any_txt[2], " required: ",
                  vector_and(missing, quotes = FALSE)),
           immediate = TRUE)
}
