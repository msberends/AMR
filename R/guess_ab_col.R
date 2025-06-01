# ==================================================================== #
# TITLE:                                                               #
# AMR: An R Package for Working with Antimicrobial Resistance Data     #
#                                                                      #
# SOURCE CODE:                                                         #
# https://github.com/msberends/AMR                                     #
#                                                                      #
# PLEASE CITE THIS SOFTWARE AS:                                        #
# Berends MS, Luz CF, Friedrich AW, et al. (2022).                     #
# AMR: An R Package for Working with Antimicrobial Resistance Data.    #
# Journal of Statistical Software, 104(3), 1-31.                       #
# https://doi.org/10.18637/jss.v104.i03                                #
#                                                                      #
# Developed at the University of Groningen and the University Medical  #
# Center Groningen in The Netherlands, in collaboration with many      #
# colleagues from around the world, see our website.                   #
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
# how to conduct AMR data analysis: https://amr-for-r.org              #
# ==================================================================== #

#' Guess Antibiotic Column
#'
#' This tries to find a column name in a data set based on information from the [antimicrobials] data set. Also supports WHONET abbreviations.
#' @param x A [data.frame].
#' @param search_string A text to search `x` for, will be checked with [as.ab()] if this value is not a column in `x`.
#' @param verbose A [logical] to indicate whether additional info should be printed.
#' @param only_sir_columns A [logical] to indicate whether only antimicrobial columns must be included that were transformed to class [sir][as.sir()] on beforehand. Defaults to `FALSE` if no columns of `x` have a class [sir][as.sir()].
#' @details You can look for an antibiotic (trade) name or abbreviation and it will search `x` and the [antimicrobials] data set for any column containing a name or code of that antibiotic.
#' @return A column name of `x`, or `NULL` when no result is found.
#' @export
#' @examples
#' df <- data.frame(
#'   amox = "S",
#'   tetr = "R"
#' )
#'
#' guess_ab_col(df, "amoxicillin")
#' guess_ab_col(df, "J01AA07") # ATC code of tetracycline
#'
#' guess_ab_col(df, "J01AA07", verbose = TRUE)
#'
#' # WHONET codes
#' df <- data.frame(
#'   AMP_ND10 = "R",
#'   AMC_ED20 = "S"
#' )
#' guess_ab_col(df, "ampicillin")
#' guess_ab_col(df, "J01CR02")
#' guess_ab_col(df, "augmentin")
guess_ab_col <- function(x = NULL, search_string = NULL, verbose = FALSE, only_sir_columns = FALSE) {
  meet_criteria(x, allow_class = "data.frame", allow_NULL = TRUE)
  meet_criteria(search_string, allow_class = "character", has_length = 1, allow_NULL = TRUE)
  meet_criteria(verbose, allow_class = "logical", has_length = 1)
  meet_criteria(only_sir_columns, allow_class = "logical", has_length = 1)

  if (is.null(x) && is.null(search_string)) {
    return(as.name("guess_ab_col"))
  } else {
    meet_criteria(search_string, allow_class = "character", has_length = 1, allow_NULL = FALSE)
  }

  all_found <- get_column_abx(x,
    info = verbose, only_sir_columns = only_sir_columns,
    verbose = verbose, fn = "guess_ab_col"
  )
  search_string.ab <- suppressWarnings(as.ab(search_string))
  ab_result <- unname(all_found[names(all_found) == search_string.ab])

  if (length(ab_result) == 0) {
    if (isTRUE(verbose)) {
      message_("No column found as input for ", search_string,
        " (", ab_name(search_string, language = NULL, tolower = TRUE), ").",
        add_fn = font_black,
        as_note = FALSE
      )
    }
    return(NULL)
  } else {
    if (isTRUE(verbose)) {
      message_(
        "Using column '", font_bold(ab_result), "' as input for ", search_string,
        " (", ab_name(search_string, language = NULL, tolower = TRUE), ")."
      )
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
                           only_sir_columns = FALSE,
                           sort = TRUE,
                           reuse_previous_result = TRUE,
                           fn = NULL,
                           return_all = FALSE) {
  # check if retrieved before, then get it from package environment
  if (isTRUE(reuse_previous_result) && identical(
    unique_call_id(
      entire_session = FALSE,
      match_fn = fn
    ),
    AMR_env$get_column_abx.call
  )) {
    # so within the same call, within the same environment, we got here again.
    # but we could've come from another function within the same call, so now only check the columns that changed

    # first remove the columns that are not existing anymore
    previous <- AMR_env$get_column_abx.out
    current <- previous[previous %in% colnames(x)]

    # then compare columns in current call with columns in original call
    new_cols <- colnames(x)[!colnames(x) %in% AMR_env$get_column_abx.checked_cols]
    if (length(new_cols) > 0) {
      # these columns did not exist in the last call, so add them
      new_cols_sir <- get_column_abx(x[, new_cols, drop = FALSE], reuse_previous_result = FALSE, info = FALSE, sort = FALSE)
      current <- c(current, new_cols_sir)
      # order according to columns in current call
      current <- current[match(colnames(x)[colnames(x) %in% current], current)]
    }

    # update pkg environment to improve speed on next run
    AMR_env$get_column_abx.out <- current
    AMR_env$get_column_abx.checked_cols <- colnames(x)

    # and return right values
    return(AMR_env$get_column_abx.out)
  }

  meet_criteria(x, allow_class = "data.frame")
  meet_criteria(soft_dependencies, allow_class = "character", allow_NULL = TRUE)
  meet_criteria(hard_dependencies, allow_class = "character", allow_NULL = TRUE)
  meet_criteria(verbose, allow_class = "logical", has_length = 1)
  meet_criteria(info, allow_class = "logical", has_length = 1)
  meet_criteria(only_sir_columns, allow_class = "logical", has_length = 1)
  meet_criteria(sort, allow_class = "logical", has_length = 1)

  if (isTRUE(info)) {
    message_("Auto-guessing columns suitable for analysis", appendLF = FALSE, as_note = FALSE)
  }

  x <- as.data.frame(x, stringsAsFactors = FALSE)
  x.bak <- x
  if (only_sir_columns == TRUE) {
    x <- x[, which(is.sir(x)), drop = FALSE]
  }

  if (NROW(x) > 10000) {
    # only test maximum of 10,000 values per column
    if (isTRUE(info)) {
      message_(" (using only ", font_bold("the first 10,000 rows"), ")...",
        appendLF = FALSE,
        as_note = FALSE
      )
    }
    x <- x[1:10000, , drop = FALSE]
  } else if (isTRUE(info)) {
    message_("...", appendLF = FALSE, as_note = FALSE)
  }

  # only check columns that are a valid AB code, ATC code, name, abbreviation or synonym,
  # or already have the 'sir' class (as.sir)
  # and that they have no more than 50% invalid values
  vectr_antibiotics <- unlist(AMR_env$AB_lookup$generalised_all)
  vectr_antibiotics <- vectr_antibiotics[!is.na(vectr_antibiotics) & nchar(vectr_antibiotics) >= 3]
  x_columns <- vapply(
    FUN.VALUE = character(1),
    colnames(x),
    function(col, df = x) {
      if (generalise_antibiotic_name(col) %in% vectr_antibiotics ||
        is.sir(x[, col, drop = TRUE]) ||
        is_sir_eligible(x[, col, drop = TRUE], threshold = 0.5)
      ) {
        return(col)
      } else {
        return(NA_character_)
      }
    }, USE.NAMES = FALSE
  )

  x_columns <- x_columns[!is.na(x_columns)]
  x <- x[, x_columns, drop = FALSE] # without drop = FALSE, x will become a vector when x_columns is length 1
  df_trans <- data.frame(
    colnames = colnames(x),
    abcode = suppressWarnings(as.ab(colnames(x), info = FALSE)),
    stringsAsFactors = FALSE
  )
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
    if (anyNA(newnames)) {
      if (isTRUE(info)) {
        message_(paste0(font_yellow(font_bold(" WARNING: ")), "some columns returned `NA` for `as.ab()`"), as_note = FALSE)
      }
      warning_("Invalid antibiotic reference(s): ", vector_and(names(dots)[is.na(newnames)], quotes = FALSE),
        call = FALSE,
        immediate = TRUE
      )
      all_okay <- FALSE
    }
    unexisting_cols <- which(!vapply(FUN.VALUE = logical(1), dots, function(col) all(col %in% x_columns)))
    if (length(unexisting_cols) > 0) {
      if (isTRUE(info)) {
        message_(" ERROR", add_fn = list(font_red, font_bold), as_note = FALSE)
      }
      stop_("Column(s) not found: ", vector_and(unlist(dots[[unexisting_cols]]), quotes = FALSE),
        call = FALSE
      )
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
    if (isTRUE(info) && all_okay == TRUE) {
      message_("No columns found.")
    }
    AMR_env$get_column_abx.call <- unique_call_id(entire_session = FALSE, match_fn = fn)
    AMR_env$get_column_abx.checked_cols <- colnames(x.bak)
    AMR_env$get_column_abx.out <- out
    return(out)
  }

  # sort on name
  if (sort == TRUE) {
    out <- out[order(names(out), out)]
  }

  dups <- FALSE

  if (return_all == FALSE) {
    dups <- names(out)[names(out) %in% names(out)[duplicated(names(out))]]
    # only keep the first hits, no duplicates
    duplicates <- c(out[duplicated(names(out))], out[duplicated(unname(out))])
    if (length(duplicates) > 0) {
      all_okay <- FALSE
    }

    if (isTRUE(info)) {
      if (all_okay == TRUE) {
        message_(" OK.", add_fn = list(font_green, font_bold), as_note = FALSE)
      } else if (!isFALSE(dups)) {
        message_(paste0(font_yellow(font_bold(" WARNING: ")), "some results from `as.ab()` are duplicated: ", vector_and(dups, quotes = "`")), as_note = FALSE)
      } else {
        message_(" WARNING.", add_fn = list(font_yellow, font_bold), as_note = FALSE)
      }

      for (i in seq_len(length(out))) {
        if (isTRUE(verbose) && !out[i] %in% duplicates) {
          message_(
            "Using column '", font_bold(out[i]), "' as input for ", names(out)[i],
            " (", ab_name(names(out)[i], tolower = TRUE, language = NULL), ")."
          )
        }
        if (out[i] %in% duplicates) {
          already_set_as <- out[which(out == out[i])[1L]]
          if (names(out)[i] != already_set_as) {
            message_(
              paste0(
                "Column '", font_bold(out[i]), "' will not be used for ",
                names(out)[i], " (", suppressMessages(ab_name(names(out)[i], tolower = TRUE, language = NULL, fast_mode = TRUE)), ")",
                ", as this antimicrobial has already been set."
              ),
              add_fn = font_red
            )
          }
        }
      }
    }

    out <- out[!duplicated(names(out))]
    out <- out[!duplicated(unname(out))]
    if (sort == TRUE) {
      out <- out[order(names(out), out)]
    }
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
    if (isTRUE(info) && !all(soft_dependencies %in% names(out))) {
      # missing a soft dependency may lower the reliability
      missing <- soft_dependencies[!soft_dependencies %in% names(out)]
      missing_msg <- vector_and(
        paste0(
          ab_name(missing, tolower = TRUE, language = NULL),
          " (", font_bold(missing, collapse = NULL), ")"
        ),
        quotes = FALSE
      )
      message_(
        "Reliability would be improved if these antimicrobial results would be available too: ",
        missing_msg
      )
    }
  }

  AMR_env$get_column_abx.call <- unique_call_id(entire_session = FALSE, match_fn = fn)
  AMR_env$get_column_abx.checked_cols <- colnames(x.bak)
  AMR_env$get_column_abx.out <- out
  out
}

get_ab_from_namespace <- function(x, cols_ab) {
  # cols_ab comes from get_column_abx()

  x <- trimws2(unique(toupper(unlist(strsplit(x, ",", fixed = TRUE)))))
  x_new <- character()
  for (val in x) {
    if (paste0("AB_", val) %in% ls(envir = asNamespace("AMR"))) {
      # antibiotic group names, as defined in data-raw/_pre_commit_checks.R, such as `AB_CARBAPENEMS`
      val <- eval(parse(text = paste0("AB_", val)), envir = asNamespace("AMR"))
    } else if (val %in% AMR_env$AB_lookup$ab) {
      # separate drugs, such as `AMX`
      val <- as.ab(val)
    } else {
      stop_("unknown antimicrobial drug (group): ", val, call = FALSE)
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
  warning_(
    paste0(
      "Introducing NAs since", any_txt[1], " these antimicrobials ", any_txt[2], " required: ",
      vector_and(missing, quotes = FALSE)
    ),
    immediate = TRUE
  )
}
