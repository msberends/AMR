# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Analysis                              #
#                                                                      #
# AUTHORS                                                              #
# Berends MS (m.s.berends@umcg.nl), Luz CF (c.f.luz@umcg.nl)           #
#                                                                      #
# LICENCE                                                              #
# This program is free software; you can redistribute it and/or modify #
# it under the terms of the GNU General Public License version 2.0,    #
# as published by the Free Software Foundation.                        #
#                                                                      #
# This program is distributed in the hope that it will be useful,      #
# but WITHOUT ANY WARRANTY; without even the implied warranty of       #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        #
# GNU General Public License for more details.                         #
# ==================================================================== #

#' @importFrom dplyr %>% pull
rsi_calc <- function(...,
                     type,
                     include_I,
                     minimum,
                     as_percent,
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
  # check integrity of columns: force rsi class
  if (is.data.frame(x)) {
    for (i in 1:ncol(x)) {
      if (!is.rsi(x %>% pull(i))) {
        x[, i] <- as.rsi(x[, i])
        print_warning <- TRUE
      }
      x[, i] <- x %>% pull(i) %>% as.integer()
    }
    x <- apply(X = x,
               MARGIN = 1,
               FUN = min)
  } else {
    if (!is.rsi(x)) {
      x <- as.rsi(x)
      print_warning <- TRUE
    }
  }

  if (print_warning == TRUE) {
    warning("Increase speed by transforming to class `rsi` on beforehand: df %>% mutate_if(is.rsi.eligible, as.rsi)",
            call. = FALSE)
  }

  if (type == "S") {
    found <- sum(as.integer(x) <= 1 + include_I, na.rm = TRUE)
  } else if (type == "I") {
    found <- sum(as.integer(x) == 2, na.rm = TRUE)
  } else if (type == "R") {
    found <- sum(as.integer(x) >= 3 - include_I, na.rm = TRUE)
  } else {
    stop("invalid type")
  }

  if (only_count == TRUE) {
    return(found)
  }

  total <- length(x) - sum(is.na(x))
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
