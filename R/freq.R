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

#' Frequency table
#'
#' Create a frequency table of a vector of data, a single column or a maximum of 9 columns of a data frame. Supports markdown for reports.
#' @param x data
#' @param sort.count Sort on count. Use \code{FALSE} to sort alphabetically on item.
#' @param nmax number of row to print. Use \code{nmax = 0} or \code{nmax = NA} to print all rows.
#' @param na.rm a logical value indicating whether NA values should be removed from the frequency table. The header will always print the amount of\code{NA}s.
#' @param markdown print table in markdown format (this forces \code{nmax = NA})
#' @param toConsole Print table to the console. Use \code{FALSE} to assign the table to an object.
#' @param digits how many significant digits are to be used for numeric values (not for the items themselves, that depends on \code{\link{getOption}("digits")})
#' @param sep a character string to separate the terms when selecting multiple columns
#' @details For numeric values, the next values will be calculated and shown into the header:
#' \itemize{
#'   \item{Mean, using \code{\link[base]{mean}}}
#'   \item{Standard deviation, using \code{\link[stats]{sd}}}
#'   \item{Five numbers of Tukey (min, Q1, median, Q3, max), using \code{\link[stats]{fivenum}}}
#'   \item{Outliers (count and list), using \code{\link{boxplot.stats}}}
#'   \item{Coefficient of variation (CV), the standard deviation divided by the mean}
#'   \item{Coefficient of quartile variation (CQV, sometimes called coefficient of dispersion), calculated as \code{(Q3 - Q1) / (Q3 + Q1)} using \code{\link{quantile}} with \code{type = 6} as quantile algorithm to comply with SPSS standards}
#' }
#' @importFrom stats fivenum sd quantile
#' @importFrom grDevices boxplot.stats
#' @importFrom dplyr %>% select pull n_distinct group_by arrange desc mutate summarise
#' @keywords summary summarise frequency freq
#' @rdname freq
#' @export
#' @examples
#' library(dplyr)
#'
#' freq(septic_patients$hospital_id)
#'
#' septic_patients %>%
#'   filter(hospital_id == "A") %>%
#'   select(bactid) %>%
#'   freq()
#'
#' # select multiple columns; they will be pasted together
#' septic_patients %>%
#'   left_join_microorganisms %>%
#'   filter(hospital_id == "A") %>%
#'   select(genus, species) %>%
#'   freq()
#'
#' # save frequency table to an object
#' years <- septic_patients %>%
#'   mutate(year = format(date, "%Y")) %>%
#'   select(year) %>%
#'   freq(toConsole = FALSE)
freq <- function(x,
                 sort.count = TRUE,
                 nmax = 15,
                 na.rm = TRUE,
                 markdown = FALSE,
                 toConsole = TRUE,
                 digits = 2,
                 sep = " ") {

  mult.columns <- 0

  if (NROW(x) == 0) {
    cat('\nNo observations.\n')
    return(invisible())
  }

  if (!is.null(ncol(x))) {
    if (ncol(x) == 1 & any(class(x) == 'data.frame')) {
      x <- x %>% pull(1)
    } else if (ncol(x) < 10) {

      mult.columns <- ncol(x)

      colnames(x) <- LETTERS[1:ncol(x)]
      if (ncol(x) == 2) {
        x$total <- paste(x$A %>% as.character(),
                          x$B %>% as.character(),
                          sep = sep)
      } else if (ncol(x) == 3) {
        x$total <- paste(x$A %>% as.character(),
                          x$B %>% as.character(),
                          x$C %>% as.character(),
                          sep = sep)
      } else if (ncol(x) == 4) {
        x$total <- paste(x$A %>% as.character(),
                          x$B %>% as.character(),
                          x$C %>% as.character(),
                          x$D %>% as.character(),
                          sep = sep)
      } else if (ncol(x) == 5) {
        x$total <- paste(x$A %>% as.character(),
                          x$B %>% as.character(),
                          x$C %>% as.character(),
                          x$D %>% as.character(),
                          x$E %>% as.character(),
                          sep = sep)
      } else if (ncol(x) == 6) {
        x$total <- paste(x$A %>% as.character(),
                          x$B %>% as.character(),
                          x$C %>% as.character(),
                          x$D %>% as.character(),
                          x$E %>% as.character(),
                          x$F %>% as.character(),
                          sep = sep)
      } else if (ncol(x) == 7) {
        x$total <- paste(x$A %>% as.character(),
                          x$B %>% as.character(),
                          x$C %>% as.character(),
                          x$D %>% as.character(),
                          x$E %>% as.character(),
                          x$F %>% as.character(),
                          x$G %>% as.character(),
                          sep = sep)
      } else if (ncol(x) == 8) {
        x$total <- paste(x$A %>% as.character(),
                          x$B %>% as.character(),
                          x$C %>% as.character(),
                          x$D %>% as.character(),
                          x$E %>% as.character(),
                          x$F %>% as.character(),
                          x$G %>% as.character(),
                          x$H %>% as.character(),
                          sep = sep)
      } else if (ncol(x) == 9) {
        x$total <- paste(x$A %>% as.character(),
                          x$B %>% as.character(),
                          x$C %>% as.character(),
                          x$D %>% as.character(),
                          x$E %>% as.character(),
                          x$F %>% as.character(),
                          x$G %>% as.character(),
                          x$H %>% as.character(),
                          x$I %>% as.character(),
                          sep = sep)
      }

      x <- x$total

    } else {
      stop('A maximum of 9 columns can be analysed at the same time.', call. = FALSE)
    }
  }
  if (markdown == TRUE & toConsole == FALSE) {
    warning('`toConsole = FALSE` will be ignored when `markdown = TRUE`.')
  }

  if (mult.columns > 1) {
    NAs <- x[is.na(x) | x == trimws(strrep2('NA ', mult.columns))]
  } else {
    NAs <- x[is.na(x)]
  }
  if (na.rm == TRUE) {
    x <- x[!x %in% NAs]
  }

  if (missing(sort.count) & any(class(x) %in% c('double', 'integer', 'numeric', 'raw', 'single', 'factor'))) {
    # sort on item/level at default when x is numeric or a factor and sort.count is not set
    sort.count <- FALSE
  }

  header <- character(0)

  markdown_line <- ''
  if (markdown == TRUE) {
    markdown_line <- '\n'
  }
  x_align <- 'l'

  if (mult.columns > 0) {
    header <- header %>% paste0(markdown_line, 'Columns:   ', mult.columns)
  } else {
    header <- header %>% paste0(markdown_line, 'Class:     ', class(x) %>% rev() %>% paste(collapse = " > "))
  }

  if (is.list(x) | is.matrix(x) | is.environment(x) | is.function(x)) {
    cat(header, "\n")
    stop('`freq()` does not support lists, matrices, environments or functions.', call. = FALSE)
  }

  header <- header %>% paste0(markdown_line, '\nLength:    ', (NAs %>% length() + x %>% length()) %>% format(),
                              ' (of which NA: ', NAs %>% length() %>% format(),
                              ' = ', (NAs %>% length() / (NAs %>% length() + x %>% length())) %>% percent(force_zero = TRUE), ')')
  header <- header %>% paste0(markdown_line, '\nUnique:    ', x %>% n_distinct() %>% format())

  header.numbers.done <- FALSE
  if (any(class(x) %in% c('double', 'integer', 'numeric', 'raw', 'single'))) {
    # right align number
    x_align <- 'r'
    header <- header %>% paste0('\n')
    header <- header %>% paste(markdown_line, '\nMean:     ', x %>% base::mean(na.rm = TRUE) %>% format(digits = digits))
    header <- header %>% paste0(markdown_line, '\nStd. dev.: ', x %>% stats::sd(na.rm = TRUE) %>% format(digits = digits),
                                ' (CV: ', x %>% cv(na.rm = TRUE) %>% format(digits = digits), ')')
    header <- header %>% paste0(markdown_line, '\nFive-Num:  ', x %>% stats::fivenum(na.rm = TRUE) %>% format(digits = digits) %>% trimws() %>% paste(collapse = '  |  '),
                                ' (CQV: ', x %>% cqv(na.rm = TRUE) %>% format(digits = digits), ')')
    outlier_length <- length(boxplot.stats(x)$out)
    header <- header %>% paste0(markdown_line, '\nOutliers:  ', outlier_length)
    if (outlier_length > 0) {
      header <- header %>% paste0(' (unique: ', boxplot.stats(x)$out %>% unique() %>% length(), ')')
    }
  }

  formatdates <- "%e %B %Y" # = d mmmm yyyy
  if (any(class(x) == 'hms')) {
    x <- x %>% as.POSIXlt()
    formatdates <- "%H:%M:%S"
  }
  if (any(class(x) %in% c('Date', 'POSIXct', 'POSIXlt'))) {
    header <- header %>% paste0('\n')
    mindatum <- x %>% min()
    maxdatum <- x %>% max()
    header <- header %>% paste0(markdown_line, '\nOldest:    ', mindatum %>% format(formatdates) %>% trimws())
    header <- header %>% paste0(markdown_line, '\nNewest:    ', maxdatum %>% format(formatdates) %>% trimws(),
                                ' (+', difftime(maxdatum, mindatum, units = 'auto') %>% as.double() %>% format(), ')')
  }
  if (any(class(x) == 'POSIXlt')) {
    x <- x %>% format(formatdates)
  }

  if (toConsole == TRUE) {
    cat(header)
  }

  if (all(is.na(x))) {
    cat('\n\nNo observations.\n')
    return(invisible())
  }
  if (n_distinct(x) == length(x)) {
    warning('All observations are unique.', call. = FALSE)
  }

  if (nmax == 0 | is.na(nmax)) {
    nmax <- length(x)
  }
  nmax.1 <- min(length(x), nmax + 1)

  # create table with counts and percentages
  if (any(class(x) == 'factor')) {
    df <- tibble::tibble(Item = x,
                         Fctlvl = x %>% as.integer()) %>%
      group_by(Item, Fctlvl)
    column_names <- c('Item', 'Count', 'Percent', 'Cum. Count', 'Cum. Percent', '(Factor Level)')
    column_align <- c('l', 'r', 'r', 'r', 'r', 'r')
  } else {
    df <- tibble::tibble(Item = x) %>%
      group_by(Item)
    column_names <- c('Item', 'Count', 'Percent', 'Cum. Count', 'Cum. Percent')
    column_align <- c(x_align, 'r', 'r', 'r', 'r')
  }
  df <- df %>%
    summarise(Count = n(),
              Percent = (n() / length(x)) %>% percent(force_zero = TRUE))

  if (df$Item %>% paste(collapse = ',') %like% '\033') {
    df <- df %>%
      mutate(Item = Item %>%
               # remove escape char
               # see https://en.wikipedia.org/wiki/Escape_character#ASCII_escape_character
               gsub('\033', ' ', ., fixed = TRUE))
  }

  # sort according to setting
  if (sort.count == TRUE) {
    df <- df %>% arrange(desc(Count))
  } else {
    if (any(class(x) == 'factor')) {
      df <- df %>% arrange(Fctlvl)
    } else {
      df <- df %>% arrange(Item)
    }
  }

  # add cumulative values
  df$Cum <- cumsum(df$Count)
  df$CumTot <- (df$Cum / sum(df$Count, na.rm = TRUE)) %>% percent(force_zero = TRUE)
  df$Cum <- df$Cum %>% format()

  if (any(class(x) == 'factor')) {
    # put factor last
    df <- df %>% select(Item, Count, Percent, Cum, CumTot, Fctlvl)
  }

  if (markdown == TRUE) {
    tblformat <- 'markdown'
  } else {
    tblformat <- 'pandoc'
  }

  if (toConsole == FALSE) {
    # assign to object
    df[, 3] <- df[, 2] / sum(df[, 2], na.rm = TRUE)
    df[, 4] <- cumsum(df[, 2])
    df[, 5] <- df[, 4] / sum(df[, 2], na.rm = TRUE)
    return(df)

  } else {

    # save old NA setting for kable
    opt.old <- options()$knitr.kable.NA
    options(knitr.kable.NA = "<NA>")

    Count.rest <- sum(df[nmax.1:nrow(df), 'Count'], na.rm = TRUE)
    if (any(class(x) %in% c('double', 'integer', 'numeric', 'raw', 'single'))) {
      df <- df %>% mutate(Item = format(Item))
    }
    df <- df %>% mutate(Count = format(Count))

    if (nrow(df) > nmax.1 & markdown == FALSE) {
      df2 <- df[1:nmax,]
      print(
        knitr::kable(df2,
                     format = tblformat,
                     col.names = column_names,
                     align = column_align,
                     padding = 1)
      )
      cat('... and ',
          format(nrow(df) - nmax),
          ' more ',
          paste0('(n = ',
                 format(Count.rest),
                 '; ',
                 (Count.rest / length(x)) %>% percent(force_zero = TRUE),
                 ')'),
          '. Use `nmax` to show more rows.\n', sep = '')

    } else {
      print(
        knitr::kable(df,
                     format = tblformat,
                     col.names = column_names,
                     align = column_align,
                     padding = 1)
      )
    }
    cat('\n')

    # reset old kable setting
    options(knitr.kable.NA = opt.old)
    return(invisible())
  }
}

#' @rdname freq
#' @export
frequency_tbl <- freq
