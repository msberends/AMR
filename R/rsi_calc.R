# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Analysis                              #
#                                                                      #
# SOURCE                                                               #
# https://gitlab.com/msberends/AMR                                     #
#                                                                      #
# LICENCE                                                              #
# (c) 2018-2020 Berends MS, Luz CF et al.                              #
#                                                                      #
# This R package is free software; you can freely use and distribute   #
# it for both personal and commercial purposes under the terms of the  #
# GNU General Public License version 2.0 (GNU GPL-2), as published by  #
# the Free Software Foundation.                                        #
#                                                                      #
# We created this package for both routine data analysis and academic  #
# research and it was publicly released in the hope that it will be    #
# useful, but it comes WITHOUT ANY WARRANTY OR LIABILITY.              #
# Visit our website for more info: https://msberends.gitlab.io/AMR.    #
# ==================================================================== #

dots2vars <- function(...) {
  # this function is to give more informative output about 
  # variable names in count_* and proportion_* functions
  dots <- substitute(list(...))
  paste(as.character(dots)[2:length(dots)], collapse = ", ")
}

rsi_calc <- function(...,
                     ab_result,
                     minimum = 0,
                     as_percent = FALSE,
                     only_all_tested = FALSE,
                     only_count = FALSE) {
  
  data_vars <- dots2vars(...)
  
  if (!is.numeric(minimum)) {
    stop("`minimum` must be numeric", call. = FALSE)
  }
  if (!is.logical(as_percent)) {
    stop("`as_percent` must be logical", call. = FALSE)
  }
  if (!is.logical(only_all_tested)) {
    stop("`only_all_tested` must be logical", call. = FALSE)
  }
  
  dots_df <- ...elt(1) # it needs this evaluation
  dots <- base::eval(base::substitute(base::alist(...)))
  if ("also_single_tested" %in% names(dots)) {
    stop("`also_single_tested` was replaced by `only_all_tested`. Please read Details in the help page (`?proportion`) as this may have a considerable impact on your analysis.", call. = FALSE)
  }
  ndots <- length(dots)
  
  if ("data.frame" %in% class(dots_df)) {
    # data.frame passed with other columns, like: example_isolates %>% proportion_S(amcl, gent)
    dots <- as.character(dots)
    dots <- dots[dots != "."]
    if (length(dots) == 0 | all(dots == "df")) {
      # for complete data.frames, like example_isolates %>% select(amcl, gent) %>% proportion_S()
      # and the old rsi function, which has "df" as name of the first parameter
      x <- dots_df
    } else {
      x <- dots_df[, dots[dots %in% colnames(dots_df)]]
    }
  } else if (ndots == 1) {
    # only 1 variable passed (can also be data.frame), like: proportion_S(example_isolates$amcl) and example_isolates$amcl %>% proportion_S()
    x <- dots_df
  } else {
    # multiple variables passed without pipe, like: proportion_S(example_isolates$amcl, example_isolates$gent)
    x <- NULL
    try(x <- as.data.frame(dots), silent = TRUE)
    if (is.null(x)) {
      # support for example_isolates %>% group_by(hospital_id) %>% summarise(amox = susceptibility(GEN, AMX))
      x <- as.data.frame(list(...))
    }
  }
  
  if (is.null(x)) {
    warning("argument is NULL (check if columns exist): returning NA", call. = FALSE)
    return(NA)
  }
  
  print_warning <- FALSE
  
  ab_result <- as.rsi(ab_result)
  
  if (is.data.frame(x)) {
    rsi_integrity_check <- character(0)
    for (i in seq_len(ncol(x))) {
      # check integrity of columns: force rsi class
      if (!is.rsi(x %>% pull(i))) {
        rsi_integrity_check <- c(rsi_integrity_check, x %>% pull(i) %>% as.character())
        x[, i] <- suppressWarnings(x %>% pull(i) %>% as.rsi()) # warning will be given later
        print_warning <- TRUE
      }
    }
    if (length(rsi_integrity_check) > 0) {
      # this will give a warning for invalid results, of all input columns (so only 1 warning)
      rsi_integrity_check <- as.rsi(rsi_integrity_check)
    }

    if (only_all_tested == TRUE) {
      # THE NUMBER OF ISOLATES WHERE *ALL* ABx ARE S/I/R
      x <- apply(X = as.data.frame(lapply(x, as.integer), stringsAsFactors = FALSE),
                 MARGIN = 1,
                 FUN = base::min)
      numerator <- sum(as.integer(x) %in% as.integer(ab_result), na.rm = TRUE)
      denominator <- length(x) - sum(is.na(x))
      
    } else {
      # THE NUMBER OF ISOLATES WHERE *ANY* ABx IS S/I/R
      other_values <- base::setdiff(c(NA, levels(ab_result)), ab_result)
      other_values_filter <- base::apply(x, 1, function(y) {
        base::all(y %in% other_values) & base::any(is.na(y))
      })
      numerator <- sum(as.logical(by(x, seq_len(nrow(x)), function(row) any(unlist(row) %in% ab_result, na.rm = TRUE))))
      denominator <- nrow(x[!other_values_filter, ])
    }
  } else {
    # x is not a data.frame
    if (!is.rsi(x)) {
      x <- as.rsi(x)
      print_warning <- TRUE
    }
    numerator <- sum(x %in% ab_result, na.rm = TRUE)
    denominator <- sum(x %in% levels(ab_result), na.rm = TRUE)
  }
  
  if (print_warning == TRUE) {
    warning("Increase speed by transforming to class `rsi` on beforehand: df %>% mutate_if(is.rsi.eligible, as.rsi)",
            call. = FALSE)
  }
  
  if (only_count == TRUE) {
    return(numerator)
  }
  
  if (denominator < minimum) {
    if (data_vars != "") {
      data_vars <- paste(" for", data_vars)
    }
    warning("Introducing NA: only ", denominator, " results available", data_vars, " (`minimum` was set to ", minimum, ").", call. = FALSE)
    fraction <- NA
  } else {
    fraction <- numerator / denominator
  }
  
  if (as_percent == TRUE) {
    percentage(fraction, digits = 1)
  } else {
    fraction
  }
}

rsi_calc_df <- function(type, # "proportion", "count" or "both"
                        data,
                        translate_ab = "name",
                        language = get_locale(),
                        minimum = 30,
                        as_percent = FALSE,
                        combine_SI = TRUE,
                        combine_IR = FALSE,
                        combine_SI_missing = FALSE) {
  
  check_dataset_integrity()
  
  if (!"data.frame" %in% class(data)) {
    stop(paste0("`", type, "_df` must be called on a data.frame"), call. = FALSE)
  }
  
  if (isTRUE(combine_IR) & isTRUE(combine_SI_missing)) {
    combine_SI <- FALSE
  }
  if (isTRUE(combine_SI) & isTRUE(combine_IR)) {
    stop("either `combine_SI` or `combine_IR` can be TRUE, not both", call. = FALSE)
  }
  
  if (!any(sapply(data, is.rsi), na.rm = TRUE)) {
    stop("No columns with class 'rsi' found. See ?as.rsi.", call. = FALSE)
  }
  
  if (as.character(translate_ab) %in% c("TRUE", "official")) {
    translate_ab <- "name"
  }

  # select only groups and antibiotics
  if (has_groups(data)) {
    data_has_groups <- TRUE
    groups <- setdiff(names(get_groups(data)), ".rows") # get_groups is from poorman.R
    data <- data[, c(groups, colnames(data)[sapply(data, is.rsi)]), drop = FALSE]
  } else {
    data_has_groups <- FALSE
    data <- data[, colnames(data)[sapply(data, is.rsi)], drop = FALSE]
  }
  
  data <- as.data.frame(data, stringsAsFactors = FALSE)
  if (isTRUE(combine_SI) | isTRUE(combine_IR)) {
    for (i in seq_len(ncol(data))) {
      if (is.rsi(data[, i, drop = TRUE])) {
        data[, i] <- as.character(data[, i, drop = TRUE])
        if (isTRUE(combine_SI)) {
          data[, i] <- gsub("(I|S)", "SI", data[, i, drop = TRUE])
        } else if (isTRUE(combine_IR)) {
          data[, i] <- gsub("(I|R)", "IR", data[, i, drop = TRUE])
        }
      }
    }
  }
  
  sum_it <- function(.data) {
    out <- data.frame(antibiotic = character(0),
                      interpretation = character(0),
                      value = double(0),
                      isolates <- integer(0),
                      stringsAsFactors = FALSE)
    if (data_has_groups) {
      group_values <- unique(.data[, which(colnames(.data) %in% groups), drop = FALSE])
      rownames(group_values) <- NULL
      .data <- .data[, which(!colnames(.data) %in% groups), drop = FALSE]
    }
    for (i in seq_len(ncol(.data))) {
      col_results <- as.data.frame(as.matrix(table(.data[, i, drop = TRUE])))
      col_results$interpretation <- rownames(col_results)
      col_results$isolates <- col_results[, 1, drop = TRUE]
      if (nrow(col_results) > 0) {
        if (sum(col_results$isolates, na.rm = TRUE) >= minimum) {
          col_results$value <- col_results$isolates / sum(col_results$isolates, na.rm = TRUE)
        } else {
          col_results$value <- rep(NA_real_, NROW(col_results))
        }
        out_new <- data.frame(antibiotic = ifelse(isFALSE(translate_ab),
                                                  colnames(.data)[i],
                                                  ab_property(colnames(.data)[i], property = translate_ab, language = language)),
                              interpretation = col_results$interpretation,
                              value = col_results$value,
                              isolates = col_results$isolates,
                              stringsAsFactors = FALSE)
        if (data_has_groups) {
          out_new <- cbind(group_values, out_new)
        }
        out <- rbind(out, out_new)
      }
    }
    out
  }
  
  # support dplyr groups
  apply_group <- function(.data, fn, groups, ...) {
    grouped <- split(x = .data, f = lapply(groups, function(x, .data) as.factor(.data[, x]), .data))
    res <- do.call(rbind, unname(lapply(grouped, fn, ...)))
    if (any(groups %in% colnames(res))) {
      class(res) <- c("grouped_data", class(res))
      attr(res, "groups") <- groups[groups %in% colnames(res)]
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
  } else if (isTRUE(combine_IR)) {
    out$interpretation <- factor(out$interpretation, levels = c("S", "IR"), ordered = TRUE)
  } else {
    out$interpretation <- as.rsi(out$interpretation)
  }
  
  if (data_has_groups) {
    # ordering by the groups and two more: "antibiotic" and "interpretation"
    out <- out[do.call("order", out[, seq_len(length(groups) + 2)]), ]
  } else {
    out <- out[order(out$antibiotic, out$interpretation), ]
  }
  
  if (type == "proportion") {
    out <- subset(out, select = -c(isolates))
  } else if (type == "count") {
    out$value <- out$isolates
    out <- subset(out, select = -c(isolates))
  } 
  
  rownames(out) <- NULL
  out
}
