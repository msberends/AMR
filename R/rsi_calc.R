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

#' @importFrom dplyr %>% pull all_vars any_vars filter_all funs mutate_all
rsi_calc <- function(...,
                     type,
                     include_I,
                     minimum,
                     as_percent,
                     also_single_tested,
                     only_count) {

  if (!is.logical(include_I)) {
    stop('`include_I` must be logical', call. = FALSE)
  }
  if (!is.numeric(minimum)) {
    stop('`minimum` must be numeric', call. = FALSE)
  }
  if (!is.logical(as_percent)) {
    stop('`as_percent` must be logical', call. = FALSE)
  }
  if (!is.logical(also_single_tested)) {
    stop('`also_single_tested` must be logical', call. = FALSE)
  }

  dots_df <- ...elt(1) # it needs this evaluation
  dots <- base::eval(base::substitute(base::alist(...)))
  ndots <- length(dots)

 if ("data.frame" %in% class(dots_df)) {
   # data.frame passed with other columns, like:
   #   septic_patients %>% portion_S(amcl, gent)
   dots <- as.character(dots)
   dots <- dots[dots != "."]
    if (length(dots) == 0 | all(dots == "df")) {
      # for complete data.frames, like septic_patients %>% select(amcl, gent) %>% portion_S()
      # and the old rsi function, that has "df" as name of the first parameter
      x <- dots_df
    } else {
      x <- dots_df[, dots]
    }
  } else if (ndots == 1) {
    # only 1 variable passed (can also be data.frame), like:
    #   portion_S(septic_patients$amcl)
    #   septic_patients$amcl %>% portion_S()
    x <- dots_df
  } else {
    # multiple variables passed without pipe, like:
    #   portion_S(septic_patients$amcl, septic_patients$gent)
    x <- NULL
    try(x <- as.data.frame(dots), silent = TRUE)
    if (is.null(x)) {
      # support for: with(septic_patients, portion_S(amcl, gent))
      x <- as.data.frame(rlang::list2(...))
    }
  }

  if (is.null(x)) {
    warning("argument is NULL (check if columns exist): returning NA", call. = FALSE)
    return(NA)
  }

  print_warning <- FALSE

  type_trans <- as.integer(as.rsi(type))
  type_others <- base::setdiff(1:3, type_trans)

  if (is.data.frame(x)) {
    rsi_integrity_check <- character(0)
    for (i in 1:ncol(x)) {
      # check integrity of columns: force rsi class
      if (!is.rsi(x %>% pull(i))) {
        rsi_integrity_check <- c(rsi_integrity_check, x %>% pull(i) %>% as.character())
        x[, i] <- suppressWarnings(as.rsi(x[, i])) # warning will be given later
        print_warning <- TRUE
      }
      x[, i] <- x %>% pull(i) %>% as.integer()
    }
    if (length(rsi_integrity_check) > 0) {
      # this will give a warning for invalid results, of all input columns (so only 1 warning)
      rsi_integrity_check <- as.rsi(rsi_integrity_check)
    }

    if (include_I == TRUE) {
      x <- x %>% mutate_all(funs(ifelse(. == 2, type_trans, .)))
    }

    if (also_single_tested == TRUE) {
      # THE CHANCE THAT AT LEAST ONE RESULT IS type
      found <- x %>% filter_all(any_vars(. == type_trans)) %>% nrow()
      # THE CHANCE THAT AT LEAST ONE RESULT IS type OR ALL ARE TESTED
      total <- found + x %>% filter_all(all_vars(. %in% type_others)) %>% nrow()
    } else {
      x <- apply(X = x,
                 MARGIN = 1,
                 FUN = min)
      found <- sum(as.integer(x) == type_trans, na.rm = TRUE)
      total <- length(x) - sum(is.na(x))
    }
  } else {
    if (!is.rsi(x)) {
      x <- as.rsi(x)
      print_warning <- TRUE
    }
    x <- as.integer(x)
    if (include_I == TRUE) {
      x[x == 2] <- type_trans
    }
    found <- sum(x == type_trans, na.rm = TRUE)
    total <- length(x) - sum(is.na(x))
  }

  if (print_warning == TRUE) {
    warning("Increase speed by transforming to class `rsi` on beforehand: df %>% mutate_if(is.rsi.eligible, as.rsi)",
            call. = FALSE)
  }

  if (only_count == TRUE) {
    return(found)
  }

  if (total < minimum) {
    warning("Introducing NA: only ", total, " results available (minimum set to ", minimum, ").", call. = FALSE)
    result <- NA
  } else {
    result <- found / total
  }

  if (as_percent == TRUE) {
    percent(result, force_zero = TRUE)
  } else {
    result
  }
}

#' @importFrom dplyr %>% summarise_if mutate select everything bind_rows
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
    stop("either `combine_SI` or `combine_IR` can be TRUE", call. = FALSE)
  }

  if (data %>% select_if(is.rsi) %>% ncol() == 0) {
    stop("No columns with class 'rsi' found. See ?as.rsi.", call. = FALSE)
  }

  if (as.character(translate_ab) %in% c("TRUE", "official")) {
    translate_ab <- "name"
  }

  get_summaryfunction <- function(int) {
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

  resS <- get_summaryfunction("S")
  resI <- get_summaryfunction("I")
  resR <- get_summaryfunction("R")
  resSI <- get_summaryfunction("SI")
  resIR <- get_summaryfunction("IR")
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
    tidyr::gather(antibiotic, value, -interpretation, -data.groups) %>%
    select(antibiotic, everything())

  if (!translate_ab == FALSE) {
    res <- res %>% mutate(antibiotic = ab_property(antibiotic, property = translate_ab, language = language))
  }

  res
}
