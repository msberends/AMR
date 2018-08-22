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

#' @importFrom dplyr %>% bind_cols pull
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

  dots_length <- ...length()
  dots <- ...elt(1) # it needs this evaluation
  dots <- rlang::exprs(...) # or this will be a list without actual values

  if ("data.frame" %in% class(dots[[1]]) & dots_length > 1) {
    # data.frame passed with other columns, like:
    #   septic_patients %>% portion_S(amcl, gent)
    df <- dots[[1]]
    dots_df <- data.frame(col1 = df[,1])
    for (i in 2:dots_length) {
      dots_col <- as.character(dots[[i]])
      if (!dots_col %in% colnames(df)) {
        stop("variable not found: ", dots_col)
      }
      dots_df <- dots_df %>% bind_cols(data.frame(df %>% pull(dots_col)))
    }
    x <- dots_df[, -1]
  } else if (dots_length == 1) {
    # only 1 variable passed (count also be data.frame), like:
    #   portion_S(septic_patients$amcl)
    #   septic_patients$amcl %>% portion_S()
    x <- dots[[1]]
  } else {
    # multiple variables passed without pipe, like:
    #   portion_S(septic_patients$amcl, septic_patients$gent)
    #   with(septic_patients, portion_S(amcl, gent))
    x <- as.data.frame(rlang::list2(...))
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
    return(NA)
  }

  if (as_percent == TRUE) {
    percent(found / total, force_zero = TRUE)
  } else {
    found / total
  }
}
