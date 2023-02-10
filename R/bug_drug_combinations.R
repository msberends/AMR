# ==================================================================== #
# TITLE                                                                #
# AMR: An R Package for Working with Antimicrobial Resistance Data     #
#                                                                      #
# SOURCE                                                               #
# https://github.com/msberends/AMR                                     #
#                                                                      #
# CITE AS                                                              #
# Berends MS, Luz CF, Friedrich AW, Sinha BNM, Albers CJ, Glasner C    #
# (2022). AMR: An R Package for Working with Antimicrobial Resistance  #
# Data. Journal of Statistical Software, 104(3), 1-31.                 #
# doi:10.18637/jss.v104.i03                                            #
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
# how to conduct AMR data analysis: https://msberends.github.io/AMR/   #
# ==================================================================== #

#' Determine Bug-Drug Combinations
#'
#' Determine antimicrobial resistance (AMR) of all bug-drug combinations in your data set where at least 30 (default) isolates are available per species. Use [format()] on the result to prettify it to a publishable/printable format, see *Examples*.
#' @inheritParams eucast_rules
#' @param combine_SI a [logical] to indicate whether values S and I should be summed, so resistance will be based on only R, defaults to `TRUE`
#' @param add_ab_group a [logical] to indicate where the group of the antimicrobials must be included as a first column
#' @param remove_intrinsic_resistant [logical] to indicate that rows and columns with 100% resistance for all tested antimicrobials must be removed from the table
#' @param FUN the function to call on the `mo` column to transform the microorganism codes, defaults to [mo_shortname()]
#' @param translate_ab a [character] of length 1 containing column names of the [antibiotics] data set
#' @param ... arguments passed on to `FUN`
#' @inheritParams sir_df
#' @inheritParams base::formatC
#' @details The function [format()] calculates the resistance per bug-drug combination. Use `combine_SI = TRUE` (default) to test R vs. S+I and `combine_SI = FALSE` to test R+I vs. S.
#' @export
#' @rdname bug_drug_combinations
#' @return The function [bug_drug_combinations()] returns a [data.frame] with columns "mo", "ab", "S", "I", "R" and "total".
#' @examples
#' # example_isolates is a data set available in the AMR package.
#' # run ?example_isolates for more info.
#' example_isolates
#' 
#' \donttest{
#' x <- bug_drug_combinations(example_isolates)
#' head(x)
#' format(x, translate_ab = "name (atc)")
#'
#' # Use FUN to change to transformation of microorganism codes
#' bug_drug_combinations(example_isolates,
#'   FUN = mo_gramstain
#' )
#'
#' bug_drug_combinations(example_isolates,
#'   FUN = function(x) {
#'     ifelse(x == as.mo("Escherichia coli"),
#'       "E. coli",
#'       "Others"
#'     )
#'   }
#' )
#' }
bug_drug_combinations <- function(x,
                                  col_mo = NULL,
                                  FUN = mo_shortname,
                                  ...) {
  meet_criteria(x, allow_class = "data.frame", contains_column_class = "sir")
  meet_criteria(col_mo, allow_class = "character", is_in = colnames(x), has_length = 1, allow_NULL = TRUE)
  meet_criteria(FUN, allow_class = "function", has_length = 1)

  # try to find columns based on type
  # -- mo
  if (is.null(col_mo)) {
    col_mo <- search_type_in_df(x = x, type = "mo")
    stop_if(is.null(col_mo), "`col_mo` must be set")
  } else {
    stop_ifnot(col_mo %in% colnames(x), "column '", col_mo, "' (`col_mo`) not found")
  }

  x.bak <- x
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  x[, col_mo] <- FUN(x[, col_mo, drop = TRUE], ...)

  unique_mo <- sort(unique(x[, col_mo, drop = TRUE]))

  # select only groups and antibiotics
  if (is_null_or_grouped_tbl(x.bak)) {
    data_has_groups <- TRUE
    groups <- get_group_names(x.bak)
    x <- x[, c(groups, col_mo, colnames(x)[vapply(FUN.VALUE = logical(1), x, is.sir)]), drop = FALSE]
  } else {
    data_has_groups <- FALSE
    x <- x[, c(col_mo, names(which(vapply(FUN.VALUE = logical(1), x, is.sir)))), drop = FALSE]
  }

  run_it <- function(x) {
    out <- data.frame(
      mo = character(0),
      ab = character(0),
      S = integer(0),
      I = integer(0),
      R = integer(0),
      total = integer(0),
      stringsAsFactors = FALSE
    )
    if (data_has_groups) {
      group_values <- unique(x[, which(colnames(x) %in% groups), drop = FALSE])
      rownames(group_values) <- NULL
      x <- x[, which(!colnames(x) %in% groups), drop = FALSE]
    }

    for (i in seq_len(length(unique_mo))) {
      # filter on MO group and only select SIR columns
      x_mo_filter <- x[which(x[, col_mo, drop = TRUE] == unique_mo[i]), names(which(vapply(FUN.VALUE = logical(1), x, is.sir))), drop = FALSE]
      # turn and merge everything
      pivot <- lapply(x_mo_filter, function(x) {
        m <- as.matrix(table(x))
        data.frame(S = m["S", ], I = m["I", ], R = m["R", ], stringsAsFactors = FALSE)
      })
      merged <- do.call(rbind, pivot)
      out_group <- data.frame(
        mo = rep(unique_mo[i], NROW(merged)),
        ab = rownames(merged),
        S = merged$S,
        I = merged$I,
        R = merged$R,
        total = merged$S + merged$I + merged$R,
        stringsAsFactors = FALSE
      )
      if (data_has_groups) {
        if (nrow(group_values) < nrow(out_group)) {
          # repeat group_values for the number of rows in out_group
          repeated <- rep(seq_len(nrow(group_values)),
            each = nrow(out_group) / nrow(group_values)
          )
          group_values <- group_values[repeated, , drop = FALSE]
        }
        out_group <- cbind(group_values, out_group)
      }
      out <- rbind(out, out_group, stringsAsFactors = FALSE)
    }
    out
  }
  # based on pm_apply_grouped_function
  apply_group <- function(.data, fn, groups, drop = FALSE, ...) {
    grouped <- pm_split_into_groups(.data, groups, drop)
    res <- do.call(rbind, unname(lapply(grouped, fn, ...)))
    if (any(groups %in% colnames(res))) {
      class(res) <- c("grouped_data", class(res))
      res <- pm_set_groups(res, groups[groups %in% colnames(res)])
    }
    res
  }

  if (data_has_groups) {
    out <- apply_group(x, "run_it", groups)
  } else {
    out <- run_it(x)
  }
  rownames(out) <- NULL
  out <- out %>% pm_arrange(mo, ab)
  out <- as_original_data_class(out, class(x.bak)) # will remove tibble groups
  structure(out, class = c("bug_drug_combinations", ifelse(data_has_groups, "grouped", character(0)), class(out)))
}

#' @method format bug_drug_combinations
#' @export
#' @rdname bug_drug_combinations
format.bug_drug_combinations <- function(x,
                                         translate_ab = "name (ab, atc)",
                                         language = get_AMR_locale(),
                                         minimum = 30,
                                         combine_SI = TRUE,
                                         add_ab_group = TRUE,
                                         remove_intrinsic_resistant = FALSE,
                                         decimal.mark = getOption("OutDec"),
                                         big.mark = ifelse(decimal.mark == ",", ".", ","),
                                         ...) {
  meet_criteria(x, allow_class = "data.frame")
  meet_criteria(translate_ab, allow_class = c("character", "logical"), has_length = 1, allow_NA = TRUE)
  language <- validate_language(language)
  meet_criteria(minimum, allow_class = c("numeric", "integer"), has_length = 1, is_positive_or_zero = TRUE, is_finite = TRUE)
  meet_criteria(combine_SI, allow_class = "logical", has_length = 1)
  meet_criteria(add_ab_group, allow_class = "logical", has_length = 1)
  meet_criteria(remove_intrinsic_resistant, allow_class = "logical", has_length = 1)
  meet_criteria(decimal.mark, allow_class = "character", has_length = 1)
  meet_criteria(big.mark, allow_class = "character", has_length = 1)

  x.bak <- x
  if (inherits(x, "grouped")) {
    # bug_drug_combinations() has been run on groups, so de-group here
    warning_("in `format()`: formatting the output of `bug_drug_combinations()` does not support grouped variables, they were ignored")
    x <- as.data.frame(x, stringsAsFactors = FALSE)
    idx <- split(seq_len(nrow(x)), paste0(x$mo, "%%", x$ab))
    x <- data.frame(
      mo = gsub("(.*)%%(.*)", "\\1", names(idx)),
      ab = gsub("(.*)%%(.*)", "\\2", names(idx)),
      S = vapply(FUN.VALUE = double(1), idx, function(i) sum(x$S[i], na.rm = TRUE)),
      I = vapply(FUN.VALUE = double(1), idx, function(i) sum(x$I[i], na.rm = TRUE)),
      R = vapply(FUN.VALUE = double(1), idx, function(i) sum(x$R[i], na.rm = TRUE)),
      total = vapply(FUN.VALUE = double(1), idx, function(i) {
        sum(x$S[i], na.rm = TRUE) +
          sum(x$I[i], na.rm = TRUE) +
          sum(x$R[i], na.rm = TRUE)
      }),
      stringsAsFactors = FALSE
    )
  }

  x <- as.data.frame(x, stringsAsFactors = FALSE)
  x <- subset(x, total >= minimum)

  if (remove_intrinsic_resistant == TRUE) {
    x <- subset(x, R != total)
  }
  if (combine_SI == TRUE) {
    x$isolates <- x$R
  } else {
    x$isolates <- x$R + x$I
  }

  give_ab_name <- function(ab, format, language) {
    format <- tolower(format)
    ab_txt <- rep(format, length(ab))
    for (i in seq_len(length(ab_txt))) {
      ab_txt[i] <- gsub("ab", as.character(as.ab(ab[i])), ab_txt[i], fixed = TRUE)
      ab_txt[i] <- gsub("cid", ab_cid(ab[i]), ab_txt[i], fixed = TRUE)
      ab_txt[i] <- gsub("group", ab_group(ab[i], language = language), ab_txt[i], fixed = TRUE)
      ab_txt[i] <- gsub("atc_group1", ab_atc_group1(ab[i], language = language), ab_txt[i], fixed = TRUE)
      ab_txt[i] <- gsub("atc_group2", ab_atc_group2(ab[i], language = language), ab_txt[i], fixed = TRUE)
      ab_txt[i] <- gsub("atc", ab_atc(ab[i], only_first = TRUE), ab_txt[i], fixed = TRUE)
      ab_txt[i] <- gsub("name", ab_name(ab[i], language = language), ab_txt[i], fixed = TRUE)
      ab_txt[i]
    }
    ab_txt
  }

  remove_NAs <- function(.data) {
    cols <- colnames(.data)
    .data <- as.data.frame(lapply(.data, function(x) ifelse(is.na(x), "", x)),
      stringsAsFactors = FALSE
    )
    colnames(.data) <- cols
    .data
  }

  create_var <- function(.data, ...) {
    dots <- list(...)
    for (i in seq_len(length(dots))) {
      .data[, names(dots)[i]] <- dots[[i]]
    }
    .data
  }

  y <- x %pm>%
    create_var(
      ab = as.ab(x$ab),
      ab_txt = give_ab_name(ab = x$ab, format = translate_ab, language = language)
    ) %pm>%
    pm_group_by(ab, ab_txt, mo) %pm>%
    pm_summarise(
      isolates = sum(isolates, na.rm = TRUE),
      total = sum(total, na.rm = TRUE)
    ) %pm>%
    pm_ungroup()

  y <- y %pm>%
    create_var(txt = paste0(
      percentage(y$isolates / y$total, decimal.mark = decimal.mark, big.mark = big.mark),
      " (", trimws(format(y$isolates, big.mark = big.mark)), "/",
      trimws(format(y$total, big.mark = big.mark)), ")"
    )) %pm>%
    pm_select(ab, ab_txt, mo, txt) %pm>%
    pm_arrange(mo)

  # replace tidyr::pivot_wider() from here
  for (i in unique(y$mo)) {
    mo_group <- y[which(y$mo == i), c("ab", "txt"), drop = FALSE]
    colnames(mo_group) <- c("ab", i)
    rownames(mo_group) <- NULL
    y <- y %pm>%
      pm_left_join(mo_group, by = "ab")
  }
  y <- y %pm>%
    pm_distinct(ab, .keep_all = TRUE) %pm>%
    pm_select(-mo, -txt) %pm>%
    # replace tidyr::pivot_wider() until here
    remove_NAs()

  select_ab_vars <- function(.data) {
    .data[, c("ab_group", "ab_txt", colnames(.data)[!colnames(.data) %in% c("ab_group", "ab_txt", "ab")]), drop = FALSE]
  }

  y <- y %pm>%
    create_var(ab_group = ab_group(y$ab, language = language)) %pm>%
    select_ab_vars() %pm>%
    pm_arrange(ab_group, ab_txt)
  y <- y %pm>%
    create_var(ab_group = ifelse(y$ab_group != pm_lag(y$ab_group) | is.na(pm_lag(y$ab_group)), y$ab_group, ""))

  if (add_ab_group == FALSE) {
    y <- y %pm>%
      pm_select(-ab_group) %pm>%
      pm_rename("Drug" = ab_txt)
    colnames(y)[1] <- translate_into_language(colnames(y)[1], language, only_unknown = FALSE)
  } else {
    y <- y %pm>%
      pm_rename(
        "Group" = ab_group,
        "Drug" = ab_txt
      )
  }

  if (!is.null(language)) {
    colnames(y) <- translate_into_language(colnames(y), language, only_unknown = FALSE)
  }

  if (remove_intrinsic_resistant == TRUE) {
    y <- y[, !vapply(FUN.VALUE = logical(1), y, function(col) all(col %like% "100", na.rm = TRUE) & !anyNA(col)), drop = FALSE]
  }

  rownames(y) <- NULL
  as_original_data_class(y, class(x.bak)) # will remove tibble groups
}

#' @method print bug_drug_combinations
#' @export
print.bug_drug_combinations <- function(x, ...) {
  x_class <- class(x)
  print(
    set_clean_class(x,
      new_class = x_class[!x_class %in% c("bug_drug_combinations", "grouped")]
    ),
    ...
  )
  message_("Use 'format()' on this result to get a publishable/printable format.",
    ifelse(inherits(x, "grouped"), " Note: The grouping variable(s) will be ignored.", ""),
    as_note = FALSE
  )
}
