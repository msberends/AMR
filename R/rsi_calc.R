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
# Visit our website for more info: https://msberends.gitab.io/AMR.     #
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
