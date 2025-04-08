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
# how to conduct AMR data analysis: https://amr-for-r.org/             #
# ==================================================================== #

dots2vars <- function(...) {
  # this function is to give more informative output about
  # variable names in count_* and proportion_* functions
  dots <- substitute(list(...))
  dots <- as.character(dots)[2:length(dots)]
  paste0(dots[dots != "."], collapse = "+")
}

sir_calc <- function(...,
                     ab_result,
                     minimum = 0,
                     as_percent = FALSE,
                     only_all_tested = FALSE,
                     only_count = FALSE) {
  meet_criteria(ab_result, allow_class = c("character", "numeric", "integer"), has_length = c(1:5))
  meet_criteria(minimum, allow_class = c("numeric", "integer"), has_length = 1, is_positive_or_zero = TRUE, is_finite = TRUE)
  meet_criteria(as_percent, allow_class = "logical", has_length = 1)
  meet_criteria(only_all_tested, allow_class = "logical", has_length = 1)
  meet_criteria(only_count, allow_class = "logical", has_length = 1)

  data_vars <- dots2vars(...)

  dots_df <- switch(1,
    ...
  )
  if (is.data.frame(dots_df)) {
    # make sure to remove all other classes like tibbles, data.tables, etc
    dots_df <- as.data.frame(dots_df, stringsAsFactors = FALSE)
  }

  dots <- eval(substitute(alist(...)))
  stop_if(length(dots) == 0, "no variables selected", call = -2)

  stop_if("also_single_tested" %in% names(dots),
    "`also_single_tested` was replaced by `only_all_tested`.\n",
    "Please read Details in the help page (`?proportion`) as this may have a considerable impact on your analysis.",
    call = -2
  )
  ndots <- length(dots)

  if (is.data.frame(dots_df)) {
    # data.frame passed with other columns, like: example_isolates %pm>% proportion_S(AMC, GEN)

    dots <- as.character(dots)
    # remove first element, it's the data.frame
    if (length(dots) == 1) {
      dots <- character(0)
    } else {
      dots <- dots[2:length(dots)]
    }
    if (length(dots) == 0 || all(dots == "df")) {
      # for complete data.frames, like example_isolates %pm>% select(AMC, GEN) %pm>% proportion_S()
      # and the old sir function, which has "df" as name of the first argument
      x <- dots_df
    } else {
      # get dots that are in column names already, and the ones that will be once evaluated using dots_df or global env
      # this is to support susceptibility(example_isolates, AMC, any_of(some_vector_with_AB_names))
      dots <- c(
        dots[dots %in% colnames(dots_df)],
        eval(parse(text = dots[!dots %in% colnames(dots_df)]), envir = dots_df, enclos = globalenv())
      )
      dots_not_exist <- dots[!dots %in% colnames(dots_df)]
      stop_if(length(dots_not_exist) > 0, "column(s) not found: ", vector_and(dots_not_exist, quotes = TRUE), call = -2)
      x <- dots_df[, dots, drop = FALSE]
    }
  } else if (ndots == 1) {
    # only 1 variable passed (can also be data.frame), like: proportion_S(example_isolates$AMC) and example_isolates$AMC %pm>% proportion_S()
    x <- dots_df
  } else {
    # multiple variables passed without pipe, like: proportion_S(example_isolates$AMC, example_isolates$GEN)
    x <- NULL
    try(x <- as.data.frame(dots, stringsAsFactors = FALSE), silent = TRUE)
    if (is.null(x)) {
      # support for example_isolates %pm>% group_by(ward) %pm>% summarise(amox = susceptibility(GEN, AMX))
      x <- as.data.frame(list(...), stringsAsFactors = FALSE)
    }
  }

  if (is.null(x)) {
    warning_("argument is NULL (check if columns exist): returning NA")
    if (as_percent == TRUE) {
      return(NA_character_)
    } else {
      return(NA_real_)
    }
  }

  print_warning <- FALSE

  ab_result <- as.sir(ab_result)

  if (is.data.frame(x)) {
    sir_integrity_check <- character(0)
    for (i in seq_len(ncol(x))) {
      # check integrity of columns: force 'sir' class
      if (!is.sir(x[, i, drop = TRUE])) {
        sir_integrity_check <- c(sir_integrity_check, as.character(x[, i, drop = TRUE]))
        x[, i] <- suppressWarnings(as.sir(x[, i, drop = TRUE])) # warning will be given later
        print_warning <- TRUE
      }
    }
    if (length(sir_integrity_check) > 0) {
      # this will give a warning for invalid results, of all input columns (so only 1 warning)
      sir_integrity_check <- as.sir(sir_integrity_check)
    }

    x_transposed <- as.list(as.data.frame(t(x), stringsAsFactors = FALSE))
    if (isTRUE(only_all_tested)) {
      # no NAs in any column
      y <- apply(
        X = as.data.frame(lapply(x, as.double), stringsAsFactors = FALSE),
        MARGIN = 1,
        FUN = min
      )
      if ("SDD" %in% ab_result && "SDD" %in% y && message_not_thrown_before("sir_calc", only_count, ab_result, entire_session = TRUE)) {
        message_("Note that `", ifelse(only_count, "count", "proportion"), "_", ifelse("S" %in% ab_result, "S", ""), "I", ifelse("R" %in% ab_result, "R", ""), "()` will also include dose-dependent susceptibility, 'SDD'. This note will be shown once for this session.", as_note = FALSE)
      }
      numerator <- sum(!is.na(y) & y %in% as.double(ab_result), na.rm = TRUE)
      denominator <- sum(vapply(FUN.VALUE = logical(1), x_transposed, function(y) !(anyNA(y))))
    } else {
      # may contain NAs in any column
      other_values <- setdiff(c(NA, levels(ab_result)), ab_result)
      if ("SDD" %in% ab_result && "SDD" %in% unlist(x_transposed) && message_not_thrown_before("sir_calc", only_count, ab_result, entire_session = TRUE)) {
        message_("Note that `", ifelse(only_count, "count", "proportion"), "_", ifelse("S" %in% ab_result, "S", ""), "I", ifelse("R" %in% ab_result, "R", ""), "()` will also include dose-dependent susceptibility, 'SDD'. This note will be shown once for this session.", as_note = FALSE)
      }
      numerator <- sum(vapply(FUN.VALUE = logical(1), x_transposed, function(y) any(y %in% ab_result, na.rm = TRUE)))
      denominator <- sum(vapply(FUN.VALUE = logical(1), x_transposed, function(y) !(all(y %in% other_values) & anyNA(y))))
    }
  } else {
    # x is not a data.frame
    if (!is.sir(x)) {
      x <- as.sir(x)
      print_warning <- TRUE
    }
    if ("SDD" %in% ab_result && "SDD" %in% x && message_not_thrown_before("sir_calc", only_count, ab_result, entire_session = TRUE)) {
      message_("Note that `", ifelse(only_count, "count", "proportion"), "_", ifelse("S" %in% ab_result, "S", ""), "I", ifelse("R" %in% ab_result, "R", ""), "()` will also include dose-dependent susceptibility, 'SDD'. This note will be shown once for this session.", as_note = FALSE)
    }
    numerator <- sum(x %in% ab_result, na.rm = TRUE)
    denominator <- sum(x %in% levels(ab_result), na.rm = TRUE)
  }

  if (print_warning == TRUE) {
    if (message_not_thrown_before("sir_calc")) {
      warning_("Increase speed by transforming to class 'sir' on beforehand:\n",
        "  your_data %>% mutate_if(is_sir_eligible, as.sir)",
        call = FALSE
      )
    }
  }

  if (only_count == TRUE) {
    return(numerator)
  }

  if (denominator < minimum) {
    if (data_vars != "") {
      data_vars <- paste(" for", data_vars)
      # also add group name if used in dplyr::group_by()
      cur_group <- import_fn("cur_group", "dplyr", error_on_fail = FALSE)
      if (!is.null(cur_group)) {
        group_df <- tryCatch(cur_group(), error = function(e) data.frame())
        if (NCOL(group_df) > 0) {
          # transform factors to characters
          group <- vapply(FUN.VALUE = character(1), group_df, function(x) {
            if (is.numeric(x)) {
              format(x)
            } else if (is.logical(x)) {
              as.character(x)
            } else {
              paste0('"', x, '"')
            }
          })
          data_vars <- paste0(data_vars, " in group: ", paste0(names(group), " = ", group, collapse = ", "))
        }
      }
    }
    warning_("Introducing NA: ",
      ifelse(denominator == 0, "no", paste("only", denominator)),
      " results available",
      data_vars,
      " (`minimum` = ", minimum, ").",
      call = FALSE
    )
    fraction <- NA_real_
  } else {
    fraction <- numerator / denominator
    fraction[is.nan(fraction)] <- NA_real_
  }

  if (as_percent == TRUE) {
    percentage(fraction, digits = 1)
  } else {
    fraction
  }
}

sir_calc_df <- function(type, # "proportion", "count" or "both"
                        data,
                        translate_ab = "name",
                        language = get_AMR_locale(),
                        minimum = 30,
                        as_percent = FALSE,
                        combine_SI = TRUE,
                        confidence_level = 0.95) {
  meet_criteria(type, is_in = c("proportion", "count", "both"), has_length = 1)
  meet_criteria(data, allow_class = "data.frame")
  data <- ascertain_sir_classes(data, "data")
  meet_criteria(translate_ab, allow_class = c("character", "logical"), has_length = 1, allow_NA = TRUE)
  language <- validate_language(language)
  meet_criteria(minimum, allow_class = c("numeric", "integer"), has_length = 1, is_positive_or_zero = TRUE, is_finite = TRUE)
  meet_criteria(as_percent, allow_class = "logical", has_length = 1)
  meet_criteria(combine_SI, allow_class = "logical", has_length = 1)
  meet_criteria(confidence_level, allow_class = "numeric", has_length = 1)

  translate_ab <- get_translate_ab(translate_ab)

  data.bak <- data
  # select only groups and antimicrobials
  if (is_null_or_grouped_tbl(data)) {
    data_has_groups <- TRUE
    groups <- get_group_names(data)
    data <- data[, c(groups, colnames(data)[vapply(FUN.VALUE = logical(1), data, is.sir)]), drop = FALSE]
  } else {
    data_has_groups <- FALSE
    data <- data[, colnames(data)[vapply(FUN.VALUE = logical(1), data, is.sir)], drop = FALSE]
  }

  data <- as.data.frame(data, stringsAsFactors = FALSE)
  if (isTRUE(combine_SI)) {
    for (i in seq_len(ncol(data))) {
      if (is.sir(data[, i, drop = TRUE])) {
        data[, i] <- as.character(data[, i, drop = TRUE])
        if ("SDD" %in% data[, i, drop = TRUE] && message_not_thrown_before("sir_calc_df", combine_SI, entire_session = TRUE)) {
          message_("Note that `sir_calc_df()` will also count dose-dependent susceptibility, 'SDD', as 'SI' when `combine_SI = TRUE`. This note will be shown once for this session.", as_note = FALSE)
        }
        data[, i] <- gsub("(I|S|SDD)", "SI", data[, i, drop = TRUE])
      }
    }
  }

  sum_it <- function(.data) {
    out <- data.frame(
      antibiotic = character(0),
      interpretation = character(0),
      value = double(0),
      ci_min = double(0),
      ci_max = double(0),
      isolates = integer(0),
      stringsAsFactors = FALSE
    )
    if (data_has_groups) {
      group_values <- unique(.data[, which(colnames(.data) %in% groups), drop = FALSE])
      rownames(group_values) <- NULL
      .data <- .data[, which(!colnames(.data) %in% groups), drop = FALSE]
    }
    for (i in seq_len(ncol(.data))) {
      values <- .data[, i, drop = TRUE]
      if (isTRUE(combine_SI)) {
        values <- factor(values, levels = c("SI", "R", "NI"), ordered = TRUE)
      } else {
        values <- factor(values, levels = c("S", "SDD", "I", "R", "NI"), ordered = TRUE)
      }
      col_results <- as.data.frame(as.matrix(table(values)), stringsAsFactors = FALSE)
      col_results$interpretation <- rownames(col_results)
      col_results$isolates <- col_results[, 1, drop = TRUE]
      if (NROW(col_results) > 0 && sum(col_results$isolates, na.rm = TRUE) > 0) {
        if (sum(col_results$isolates, na.rm = TRUE) >= minimum) {
          col_results$value <- col_results$isolates / sum(col_results$isolates, na.rm = TRUE)
          ci <- lapply(
            col_results$isolates,
            function(x) {
              stats::binom.test(
                x = x,
                n = sum(col_results$isolates, na.rm = TRUE),
                conf.level = confidence_level
              )$conf.int
            }
          )
          col_results$ci_min <- vapply(FUN.VALUE = double(1), ci, `[`, 1)
          col_results$ci_max <- vapply(FUN.VALUE = double(1), ci, `[`, 2)
        } else {
          col_results$value <- rep(NA_real_, NROW(col_results))
          # confidence intervals also to NA
          col_results$ci_min <- col_results$value
          col_results$ci_max <- col_results$value
        }
        out_new <- data.frame(
          antibiotic = ifelse(isFALSE(translate_ab),
            colnames(.data)[i],
            ab_property(colnames(.data)[i], property = translate_ab, language = language)
          ),
          interpretation = col_results$interpretation,
          value = col_results$value,
          ci_min = col_results$ci_min,
          ci_max = col_results$ci_max,
          isolates = col_results$isolates,
          stringsAsFactors = FALSE
        )
        if (data_has_groups) {
          if (nrow(group_values) < nrow(out_new)) {
            # repeat group_values for the number of rows in out_new
            repeated <- rep(seq_len(nrow(group_values)),
              each = nrow(out_new) / nrow(group_values)
            )
            group_values <- group_values[repeated, , drop = FALSE]
          }
          out_new <- cbind(group_values, out_new)
        }
        out <- rbind_AMR(out, out_new)
      }
    }
    out
  }

  # based on pm_apply_grouped_function
  apply_group <- function(.data, fn, groups, drop = FALSE, ...) {
    grouped <- pm_split_into_groups(.data, groups, drop)
    res <- do.call(rbind_AMR, unname(lapply(grouped, fn, ...)))
    if (any(groups %in% colnames(res))) {
      class(res) <- c("grouped_data", class(res))
      res <- pm_set_groups(res, groups[groups %in% colnames(res)])
    }
    res
  }

  if (data_has_groups) {
    out <- apply_group(data, "sum_it", groups)
  } else {
    out <- sum_it(data)
  }

  # apply factors for right sorting in interpretation
  if (isTRUE(combine_SI)) {
    out$interpretation <- factor(out$interpretation, levels = c("SI", "R"), ordered = TRUE)
  } else {
    # don't use as.sir() here, as it would add the class 'sir' and we would like
    # the same data structure as output, regardless of input
    if (out$value[out$interpretation == "SDD"] > 0) {
      out$interpretation <- factor(out$interpretation, levels = c("S", "SDD", "I", "R"), ordered = TRUE)
    } else {
      out$interpretation <- factor(out$interpretation, levels = c("S", "I", "R"), ordered = TRUE)
    }
  }

  out <- out[!is.na(out$interpretation), , drop = FALSE]

  if (data_has_groups) {
    # ordering by the groups and two more: "antibiotic" and "interpretation"
    out <- pm_ungroup(out[do.call("order", out[, seq_len(length(groups) + 2), drop = FALSE]), , drop = FALSE])
  } else {
    out <- out[order(out$antibiotic, out$interpretation), , drop = FALSE]
  }

  if (type == "proportion") {
    # remove number of isolates
    out <- subset(out, select = -c(isolates))
  } else if (type == "count") {
    # set value to be number of isolates
    out$value <- out$isolates
    # remove redundant columns
    out <- subset(out, select = -c(ci_min, ci_max, isolates))
  }

  as_original_data_class(out, class(data.bak), extra_class = "sir_df") # will remove tibble groups
}
