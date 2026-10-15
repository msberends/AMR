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
# Center Groningen in the Netherlands, in collaboration with many      #
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

dots2vars <- function(...) {
  # this function is to give more informative output about
  # variable names in count_* and proportion_* functions
  dots <- substitute(list(...))
  dots <- as.character(dots)[2:length(dots)]
  paste0(dots[dots != "."], collapse = "+")
}

validate_sir_guideline <- function(guideline, fn, allow_none = FALSE) {
  # returns "EUCAST" or "CLSI" to consider 'I' susceptible or resistant, or NULL to not merge 'I' at all, see the section 'Interpretation of SIR' in ?as.sir
  if (is.null(guideline) || isFALSE(guideline)) {
    stop_ifnot(
      allow_none,
      "{.arg guideline} must be {.val EUCAST} or {.val CLSI} in {.help [{.fun ", fn, "}](AMR::", fn, ")}, since a single proportion or count requires the 'I' category to be considered either susceptible (EUCAST) or resistant (CLSI). See the section 'Interpretation of SIR' in {.help [{.fun as.sir}](AMR::as.sir)}.",
      call = FALSE
    )
    return(NULL)
  }
  stop_ifnot(
    is.character(guideline) && length(guideline) == 1 && !is.na(guideline) && (guideline %like% "EUCAST" || guideline %like% "CLSI"),
    "{.arg guideline} must be {.val EUCAST} or {.val CLSI}", ifelse(allow_none, ", or {.code NULL} to keep S, I, and R separate,", ""), " in {.help [{.fun ", fn, "}](AMR::", fn, ")}. See the section 'Interpretation of SIR' in {.help [{.fun as.sir}](AMR::as.sir)}.",
    call = FALSE
  )
  if (guideline %like% "EUCAST") "EUCAST" else "CLSI"
}

message_guideline_assumed <- function(fn, txt) {
  # informs once per session that EUCAST is assumed, only to be called when the `guideline` argument is missing
  # - not when the AMR_guideline option is set, since that is a deliberate choice
  # - not when knitting R Markdown or Quarto documents, so it will not end up in reports
  if (is.null(getOption("AMR_guideline")) && !isTRUE(getOption("knitr.in.progress")) && message_not_thrown_before(fn, "guideline_assumed", entire_session = TRUE)) {
    message_(
      "{.help [{.fun ", fn, "}](AMR::", fn, ")} assumes the EUCAST guideline and thus considers the 'I' category susceptible", txt,
      " Set the {.arg guideline} argument, or set the {.code AMR_guideline} option to change the default (see {.topic [AMR-options](AMR::AMR-options)}). This message is shown once per session."
    )
  }
}

get_aggregation_guideline <- function(guideline, guideline_missing, combine_SI, combine_SI_missing, fn) {
  # for functions that return results per category: returns "EUCAST" ("SI" vs. "R"), "CLSI" ("S" vs. "IR"), or NULL (S, I, and R separately)
  if (!combine_SI_missing) {
    # deprecated argument, if set it must be followed literally, regardless of the guideline
    meet_criteria(combine_SI, allow_class = "logical", has_length = 1, .call_depth = 1)
    deprecation_warning("combine_SI", "guideline", fn = fn, is_argument = TRUE, extra_msg = "As long as `combine_SI` is set, it takes precedence over `guideline`. Use `guideline = \"EUCAST\"` to merge 'I' with S, `guideline = \"CLSI\"` to merge 'I' with R, or `guideline = NULL` to keep S, I, and R separate.")
    if (isTRUE(combine_SI)) {
      return("EUCAST")
    } else {
      return(NULL)
    }
  }
  guideline <- validate_sir_guideline(guideline, fn = fn, allow_none = TRUE)
  if (isTRUE(guideline_missing) && !is.null(guideline)) {
    message_guideline_assumed(fn, txt = paste0(" (S+SDD+I vs. R). Use {.code guideline = \"CLSI\"} for S+SDD vs. I+R, or {.code guideline = NULL} to keep S, I, and R separate."))
  }
  guideline
}

sir_calc <- function(...,
                     ab_result,
                     minimum = 0,
                     as_percent = FALSE,
                     only_all_tested = FALSE,
                     only_count = FALSE,
                     guideline = getOption("AMR_guideline", "EUCAST")) {
  meet_criteria(ab_result, allow_class = c("character", "sir"), has_length = seq_along(VALID_SIR_LEVELS), is_in = VALID_SIR_LEVELS)
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
  denominator_vals <- levels(ab_result)
  denominator_vals <- denominator_vals[denominator_vals != "NI"]

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
    if (ncol(x) == 1) {
      x <- x[, 1, drop = TRUE]
    } else if (isTRUE(only_all_tested)) {
      # combination therapy of isolates tested for all antimicrobials, see ?sir_combine
      x <- sir_combine_engine(as.list(x), S_values = get_susceptible_values(guideline), only_all_tested = TRUE)
    }
  } else if (!is.sir(x)) {
    x <- as.sir(x)
    print_warning <- TRUE
  }

  if ("SDD" %in% ab_result && "SDD" %in% unlist(lapply(as.list(x), as.character)) && message_not_thrown_before("sir_calc", only_count, ab_result, entire_session = TRUE)) {
    message_("Note that {.fun ", ifelse(only_count, "count", "proportion"), "_", ifelse("S" %in% ab_result, "S", ""), "I", ifelse("R" %in% ab_result, "R", ""), "} will also include dose-dependent susceptibility, {.val SDD}. This note will be shown once for this session.", as_note = FALSE)
  }
  if (is.data.frame(x)) {
    # multiple antimicrobials with only_all_tested = FALSE: an isolate counts if any antimicrobial has a result in `ab_result`, and is considered tested if it has such a result or was tested for all antimicrobials, see ?sir_combine
    # this equals sir_combine(method = "most_susceptible") for susceptibility and sir_combine(method = "most_resistant") for resistance
    # 'NI' is a missing result
    numerator_vals <- intersect(as.character(ab_result), denominator_vals)
    is_counted <- Reduce(`|`, lapply(x, function(col) as.character(col) %in% numerator_vals))
    is_complete <- Reduce(`&`, lapply(x, function(col) as.character(col) %in% denominator_vals))
    numerator <- sum(is_counted)
    denominator <- sum(is_counted | is_complete)
  } else {
    numerator <- sum(x %in% ab_result, na.rm = TRUE)
    denominator <- sum(x %in% denominator_vals, na.rm = TRUE)
  }

  if (print_warning == TRUE) {
    if (message_not_thrown_before("sir_calc")) {
      warning_("Increase speed by transforming to class {.cls sir} on beforehand:\n",
        highlight_code("  your_data %>% mutate_if(is_sir_eligible, as.sir)"),
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
      " (whilst {.arg minimum = ", minimum, "}).",
      call = FALSE
    )
    fraction <- NA_real_
  } else {
    fraction <- numerator / denominator
    fraction[is.nan(fraction)] <- NA_real_
  }

  if (as_percent == TRUE) {
    trimws(percentage(fraction, digits = 1))
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
                        confidence_level = 0.95,
                        guideline = "EUCAST") { # must be "EUCAST", "CLSI", or NULL, see get_aggregation_guideline()
  meet_criteria(type, is_in = c("proportion", "count", "both"), has_length = 1)
  meet_criteria(data, allow_class = "data.frame")
  data <- ascertain_sir_classes(data, "data")
  meet_criteria(translate_ab, allow_class = c("character", "logical"), has_length = 1, allow_NA = TRUE)
  language <- validate_language(language)
  meet_criteria(minimum, allow_class = c("numeric", "integer"), has_length = 1, is_positive_or_zero = TRUE, is_finite = TRUE)
  meet_criteria(as_percent, allow_class = "logical", has_length = 1)
  meet_criteria(confidence_level, allow_class = "numeric", has_length = 1)
  meet_criteria(guideline, allow_class = "character", has_length = 1, is_in = c("EUCAST", "CLSI"), allow_NULL = TRUE)
  # with EUCAST, 'I' is merged with S into 'SI', with CLSI, 'I' is merged with R into 'IR', and with NULL, 'I' is not merged
  combine_SI <- !is.null(guideline)
  combined_levels <- if (identical(guideline, "CLSI")) c("S", "IR") else c("SI", "R")

  translate_ab <- get_translate_ab(translate_ab)

  data.bak <- data
  # select only groups and antibiotics
  if (is_null_or_grouped_tbl(data)) {
    data_has_groups <- TRUE
    groups <- get_group_names(data)
    data <- data[, c(groups, colnames(data)[vapply(FUN.VALUE = logical(1), data, is.sir)]), drop = FALSE]
  } else {
    data_has_groups <- FALSE
    data <- data[, colnames(data)[vapply(FUN.VALUE = logical(1), data, is.sir)], drop = FALSE]
  }

  data <- as.data.frame(data, stringsAsFactors = FALSE)

  for (i in seq_len(ncol(data))) {
    # transform SIR columns
    if (is.sir(data[, i, drop = TRUE])) {
      data[, i] <- as.character(as.sir(data[, i, drop = TRUE]))
      data[which(data[, i, drop = TRUE] %in% c("S", "SDD", "WT")), i] <- "S"
      data[which(data[, i, drop = TRUE] %in% c("R", "NWT", "NS")), i] <- "R"
      if (isTRUE(combine_SI) && identical(guideline, "CLSI")) {
        data[which(data[, i, drop = TRUE] %in% c("I", "R")), i] <- "IR"
      } else if (isTRUE(combine_SI)) {
        data[which(data[, i, drop = TRUE] %in% c("I", "S")), i] <- "SI"
      }
      data[which(!data[, i, drop = TRUE] %in% c("S", "SI", "I", "IR", "R")), i] <- NA_character_
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
        values <- factor(values, levels = c(combined_levels, "NI"), ordered = TRUE)
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
    out$interpretation <- factor(out$interpretation, levels = combined_levels, ordered = TRUE)
  } else {
    # don't use as.sir() here, as it would add the class 'sir' and we would like
    # the same data structure as output, regardless of input
    if (any(out$value[out$interpretation == "SDD"] > 0, na.rm = TRUE)) {
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
