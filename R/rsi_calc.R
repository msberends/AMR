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

#' @importFrom rlang enquos as_label
dots2vars <- function(...) {
  # this function is to give more informative output about 
  # variable names in count_* and portion_* functions
  paste(
    unlist(
      lapply(enquos(...),
             function(x) {
               l <- as_label(x)
               if (l != ".") {
                 l
               } else {
                 character(0)
               }
             })
    ),
    collapse = ", ")
}

#' @importFrom dplyr %>% pull all_vars any_vars filter_all funs mutate_all
# @importFrom clean percentage
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
    stop("`also_single_tested` was replaced by `only_all_tested`. Please read Details in the help page (`?portion`) as this may have a considerable impact on your analysis.", call. = FALSE)
  }
  ndots <- length(dots)

 if ("data.frame" %in% class(dots_df)) {
   # data.frame passed with other columns, like: example_isolates %>% portion_S(amcl, gent)
   dots <- as.character(dots)
   dots <- dots[dots != "."]
    if (length(dots) == 0 | all(dots == "df")) {
      # for complete data.frames, like example_isolates %>% select(amcl, gent) %>% portion_S()
      # and the old rsi function, that has "df" as name of the first parameter
      x <- dots_df
    } else {
      x <- dots_df[, dots]
    }
  } else if (ndots == 1) {
    # only 1 variable passed (can also be data.frame), like: portion_S(example_isolates$amcl) and example_isolates$amcl %>% portion_S()
    x <- dots_df
  } else {
    # multiple variables passed without pipe, like: portion_S(example_isolates$amcl, example_isolates$gent)
    x <- NULL
    try(x <- as.data.frame(dots), silent = TRUE)
    if (is.null(x)) {
      # support for: with(example_isolates, portion_S(amcl, gent))
      x <- as.data.frame(rlang::list2(...))
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
      x <- apply(X = x %>% mutate_all(as.integer),
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
      numerator <- x %>% filter_all(any_vars(. %in% ab_result)) %>% nrow()
      denominator <- x %>% filter(!other_values_filter) %>% nrow()
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

#' @importFrom dplyr %>% summarise_if mutate select everything bind_rows
#' @importFrom tidyr gather
rsi_calc_df <- function(type, # "portion" or "count"
                        data,
                        translate_ab = "name",
                        language = get_locale(),
                        minimum = 30,
                        as_percent = FALSE,
                        combine_SI = TRUE,
                        combine_IR = FALSE,
                        combine_SI_missing = FALSE) {

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

  get_summaryfunction <- function(int, type) {
    # look for portion_S, count_S, etc:
    int_fn <- get(paste0(type, "_", int), envir = asNamespace("AMR"))

    if (type == "portion") {
      summ <- summarise_if(.tbl = data,
                           .predicate = is.rsi,
                           .funs = int_fn,
                           minimum = minimum,
                           as_percent = as_percent)
    } else if (type == "count") {
      summ <- summarise_if(.tbl = data,
                           .predicate = is.rsi,
                           .funs = int_fn)
    }
    summ %>%
      mutate(interpretation = int) %>%
      select(interpretation, everything())
  }

  resS <- get_summaryfunction("S", type)
  resI <- get_summaryfunction("I", type)
  resR <- get_summaryfunction("R", type)
  resSI <- get_summaryfunction("SI", type)
  resIR <- get_summaryfunction("IR", type)
  data.groups <- group_vars(data)

  if (isFALSE(combine_SI) & isFALSE(combine_IR)) {
    res <- bind_rows(resS, resI, resR) %>%
      mutate(interpretation = factor(interpretation,
                                     levels = c("S", "I", "R"),
                                     ordered = TRUE))

  } else if (isTRUE(combine_IR)) {
    res <- bind_rows(resS, resIR) %>%
      mutate(interpretation = factor(interpretation,
                                     levels = c("S", "IR"),
                                     ordered = TRUE))

  } else if (isTRUE(combine_SI)) {
    res <- bind_rows(resSI, resR) %>%
      mutate(interpretation = factor(interpretation,
                                     levels = c("SI", "R"),
                                     ordered = TRUE))
  }

  res <- res %>%
    gather(antibiotic, value, -interpretation, -data.groups) %>%
    select(antibiotic, everything())

  if (!translate_ab == FALSE) {
    res <- res %>% mutate(antibiotic = AMR::ab_property(antibiotic, property = translate_ab, language = language))
  }

  res
}
