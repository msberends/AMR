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
#' Create a frequency table of a vector with items or a data frame. Supports quasiquotation and markdown for reports. \code{top_freq} can be used to get the top/bottom \emph{n} items of a frequency table, with counts as names.
#' @param x vector of any class or a \code{\link{data.frame}}, \code{\link{tibble}} or \code{\link{table}}
#' @param ... up to nine different columns of \code{x} when \code{x} is a \code{data.frame} or \code{tibble}, to calculate frequencies from - see Examples
#' @param sort.count sort on count, i.e. frequencies. This will be \code{TRUE} at default for everything except for factors.
#' @param nmax number of row to print. The default, \code{15}, uses \code{\link{getOption}("max.print.freq")}. Use \code{nmax = 0}, \code{nmax = Inf}, \code{nmax = NULL} or \code{nmax = NA} to print all rows.
#' @param na.rm a logical value indicating whether \code{NA} values should be removed from the frequency table. The header will always print the amount of \code{NA}s.
#' @param row.names a logical value indicating whether row indices should be printed as \code{1:nrow(x)}
#' @param markdown print table in markdown format (this forces \code{nmax = NA})
#' @param digits how many significant digits are to be used for numeric values in the header (not for the items themselves, that depends on \code{\link{getOption}("digits")})
#' @param sep a character string to separate the terms when selecting multiple columns
#' @param f a frequency table
#' @param n number of top \emph{n} items to return, use -n for the bottom \emph{n} items. It will include more than \code{n} rows if there are ties.
#' @details Frequency tables (or frequency distributions) are summaries of the distribution of values in a sample. With the `freq` function, you can create univariate frequency tables. Multiple variables will be pasted into one variable, so it forces a univariate distribution. This package also has a vignette available to explain the use of this function further, run \code{browseVignettes("AMR")} to read it.
#'
#' For numeric values of any class, these additional values will all be calculated with \code{na.rm = TRUE} and shown into the header:
#' \itemize{
#'   \item{Mean, using \code{\link[base]{mean}}}
#'   \item{Standard Deviation, using \code{\link[stats]{sd}}}
#'   \item{Coefficient of Variation (CV), the standard deviation divided by the mean}
#'   \item{Mean Absolute Deviation (MAD), using \code{\link[stats]{mad}}}
#'   \item{Tukey Five-Number Summaries (minimum, Q1, median, Q3, maximum), using \code{\link[stats]{fivenum}}}
#'   \item{Interquartile Range (IQR) calculated as \code{Q3 - Q1} using the Tukey Five-Number Summaries, i.e. \strong{not} using the \code{\link[stats]{quantile}} function}
#'   \item{Coefficient of Quartile Variation (CQV, sometimes called coefficient of dispersion), calculated as \code{(Q3 - Q1) / (Q3 + Q1)} using the Tukey Five-Number Summaries}
#'   \item{Outliers (total count and unique count), using \code{\link[grDevices]{boxplot.stats}}}
#' }
#'
#' For dates and times of any class, these additional values will be calculated with \code{na.rm = TRUE} and shown into the header:
#' \itemize{
#'   \item{Oldest, using \code{\link{min}}}
#'   \item{Newest, using \code{\link{max}}, with difference between newest and oldest}
#'   \item{Median, using \code{\link[stats]{median}}, with percentage since oldest}
#' }
#'
#'
#' The function \code{top_freq} uses \code{\link[dplyr]{top_n}} internally and will include more than \code{n} rows if there are ties.
#' @importFrom stats fivenum sd mad
#' @importFrom grDevices boxplot.stats
#' @importFrom dplyr %>% select pull n_distinct group_by arrange desc mutate summarise n_distinct tibble
#' @importFrom utils browseVignettes installed.packages
#' @keywords summary summarise frequency freq
#' @rdname freq
#' @name freq
#' @return A \code{data.frame} with an additional class \code{"frequency_tbl"}
#' @export
#' @examples
#' library(dplyr)
#'
#' # this all gives the same result:
#' freq(septic_patients$hospital_id)
#' freq(septic_patients[, "hospital_id"])
#' septic_patients$hospital_id %>% freq()
#' septic_patients[, "hospital_id"] %>% freq()
#' septic_patients %>% freq("hospital_id")
#' septic_patients %>% freq(hospital_id)  #<- easiest to remember when you're used to tidyverse
#'
#' # you could also use `select` or `pull` to get your variables
#' septic_patients %>%
#'   filter(hospital_id == "A") %>%
#'   select(bactid) %>%
#'   freq()
#'
#' # multiple selected variables will be pasted together
#' septic_patients %>%
#'   left_join_microorganisms %>%
#'   filter(hospital_id == "A") %>%
#'   freq(genus, species)
#'
#' # get top 10 bugs of hospital A as a vector
#' septic_patients %>%
#'   filter(hospital_id == "A") %>%
#'   freq(bactid) %>%
#'   top_freq(10)
#'
#' # save frequency table to an object
#' years <- septic_patients %>%
#'   mutate(year = format(date, "%Y")) %>%
#'   freq(year)
#'
#' # show only the top 5
#' years %>% print(nmax = 5)
#'
#' # save to an object with formatted percentages
#' years <- format(years)
#'
#' # print a histogram of numeric values
#' septic_patients %>%
#'   freq(age) %>%
#'   hist()  # prettier: ggplot(septic_patients, aes(age)) + geom_histogram()
#'
#' # or print all points to a regular plot
#' septic_patients %>%
#'   freq(age) %>%
#'   plot()
#'
#' # transform to a data.frame or tibble
#' septic_patients %>%
#'   freq(age) %>%
#'   as.data.frame()
#'
#' # or transform (back) to a vector
#' septic_patients %>%
#'   freq(age) %>%
#'   as.vector()
#'
#' identical(septic_patients %>%
#'             freq(age) %>%
#'             as.vector() %>%
#'             sort(),
#'           sort(septic_patients$age)) # TRUE
#'
#' # it also supports `table` objects:
#' table(septic_patients$sex,
#'       septic_patients$age) %>%
#'   freq(sep = " **sep** ")
#'
#' \dontrun{
#' # send frequency table to clipboard (e.g. for pasting in Excel)
#' septic_patients %>%
#'   freq(age) %>%
#'   format() %>%       # this will format the percentages
#'   clipboard_export()
#' }
frequency_tbl <- function(x,
                          ...,
                          sort.count = TRUE,
                          nmax = getOption("max.print.freq"),
                          na.rm = TRUE,
                          row.names = TRUE,
                          markdown = FALSE,
                          digits = 2,
                          sep = " ") {

  mult.columns <- 0

  x.name <- NULL
  cols <- NULL
  if (any(class(x) == 'list')) {
    cols <- names(x)
    x <- as.data.frame(x, stringsAsFactors = FALSE)
    x.name <- "a list"
  } else if (any(class(x) == 'matrix')) {
    x <- as.data.frame(x, stringsAsFactors = FALSE)
    x.name <- "a matrix"
    cols <- colnames(x)
    if (all(cols %like% 'V[0-9]')) {
      cols <- NULL
    }
  }

  if (any(class(x) == 'data.frame')) {
    if (is.null(x.name)) {
      x.name <- deparse(substitute(x))
    }
    if (x.name == ".") {
      x.name <- NULL
    }
    dots <- base::eval(base::substitute(base::alist(...)))
    ndots <- length(dots)

    if (ndots < 10) {
      cols <- as.character(dots)
      if (!all(cols %in% colnames(x))) {
        stop("one or more columns not found: `", paste(cols, collapse = "`, `"), '`', call. = FALSE)
      }
      if (length(cols) > 0) {
        x <- x[, cols]
      }
    } else if (ndots >= 10) {
      stop('A maximum of 9 columns can be analysed at the same time.', call. = FALSE)
    } else {
      cols <- NULL
    }
  } else if (any(class(x) == 'table')) {
    if (!"tidyr" %in% rownames(installed.packages())) {
      stop('transformation from `table` to frequency table requires the tidyr package.', call. = FALSE)
    }
    x <- x %>%
      as.data.frame(stringsAsFactors = FALSE) %>%
      # paste first two columns
      tidyr::unite(col = "Pasted", 1:2, sep = sep, remove = TRUE)
    x <- rep(x %>% pull(Pasted), x %>% pull(Freq))
    x.name <- "a `table` object"
    cols <- NULL
    mult.columns <- 2
  } else {
    x.name <- NULL
    cols <- NULL
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

  if (mult.columns > 1) {
    NAs <- x[is.na(x) | x == trimws(strrep('NA ', mult.columns))]
  } else {
    NAs <- x[is.na(x)]
  }

  if (na.rm == TRUE) {
    x_class <- class(x)
    x <- x[!x %in% NAs]
    class(x) <- x_class
  }

  if (missing(sort.count) & 'factor' %in% class(x)) {
    # sort on factor level at default when x is a factor and sort.count is not set
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

  header <- header %>% paste0(markdown_line, '\nLength:    ', (NAs %>% length() + x %>% length()) %>% format(),
                              ' (of which NA: ', NAs %>% length() %>% format(),
                              ' = ', (NAs %>% length() / (NAs %>% length() + x %>% length())) %>% percent(force_zero = TRUE, round = digits) %>% sub('NaN', '0', ., fixed = TRUE), ')')
  header <- header %>% paste0(markdown_line, '\nUnique:    ', x %>% n_distinct() %>% format())

  if (NROW(x) > 0 & any(class(x) %in% c('double', 'integer', 'numeric', 'raw', 'single'))) {
    # right align number
    Tukey_five <- stats::fivenum(x, na.rm = TRUE)
    x_align <- 'r'
    header <- header %>% paste0('\n')
    header <- header %>% paste(markdown_line, '\nMean:     ', x %>% base::mean(na.rm = TRUE) %>% format(digits = digits))
    header <- header %>% paste0(markdown_line, '\nStd. dev.: ', x %>% stats::sd(na.rm = TRUE) %>% format(digits = digits),
                                ' (CV: ', x %>% cv(na.rm = TRUE) %>% format(digits = digits),
                                ', MAD: ', x %>% stats::mad(na.rm = TRUE) %>% format(digits = digits), ')')
    header <- header %>% paste0(markdown_line, '\nFive-Num:  ', Tukey_five %>% format(digits = digits) %>% trimws() %>% paste(collapse = ' | '),
                                ' (IQR: ', (Tukey_five[4] - Tukey_five[2]) %>% format(digits = digits),
                                ', CQV: ', x %>% cqv(na.rm = TRUE) %>% format(digits = digits), ')')
    outlier_length <- length(boxplot.stats(x)$out)
    header <- header %>% paste0(markdown_line, '\nOutliers:  ', outlier_length)
    if (outlier_length > 0) {
      header <- header %>% paste0(' (unique: ', boxplot.stats(x)$out %>% n_distinct(), ')')
    }
  }
  if (NROW(x) > 0 & any(class(x) == "rsi")) {
    header <- header %>% paste0('\n')
    cnt_S <- sum(x == "S")
    cnt_I <- sum(x == "I")
    cnt_R <- sum(x == "R")
    header <- header %>% paste(markdown_line, '\n%IR:      ',
                               ((cnt_I + cnt_R) / sum(!is.na(x))) %>% percent(force_zero = TRUE, round = digits))
    header <- header %>% paste0(markdown_line, '\nRatio SIR: 1.0 : ',
                                (cnt_I / cnt_S) %>% format(digits = 1, nsmall = 1), " : ",
                                (cnt_R / cnt_S) %>% format(digits = 1, nsmall = 1))
  }

  formatdates <- "%e %B %Y" # = d mmmm yyyy
  if (any(class(x) == 'hms')) {
    x <- x %>% as.POSIXlt()
    formatdates <- "%H:%M:%S"
  }
  if (NROW(x) > 0 & any(class(x) %in% c('Date', 'POSIXct', 'POSIXlt'))) {
    header <- header %>% paste0('\n')
    mindate <- x %>% min(na.rm = TRUE)
    maxdate <- x %>% max(na.rm = TRUE)
    maxdate_days <- difftime(maxdate, mindate, units = 'auto') %>% as.double()
    mediandate <- x %>% median(na.rm = TRUE)
    median_days <- difftime(mediandate, mindate, units = 'auto') %>% as.double()

    if (formatdates == "%H:%M:%S") {
      # hms
      header <- header %>% paste0(markdown_line, '\nEarliest:  ', mindate %>% format(formatdates) %>% trimws())
      header <- header %>% paste0(markdown_line, '\nLatest:    ', maxdate %>% format(formatdates) %>% trimws(),
                                  ' (+', difftime(maxdate, mindate, units = 'mins') %>% as.double() %>% format(digits = digits), ' min.)')
    } else {
      # other date formats
      header <- header %>% paste0(markdown_line, '\nOldest:    ', mindate %>% format(formatdates) %>% trimws())
      header <- header %>% paste0(markdown_line, '\nNewest:    ', maxdate %>% format(formatdates) %>% trimws(),
                                  ' (+', difftime(maxdate, mindate, units = 'auto') %>% as.double() %>% format(digits = digits), ')')
    }
    header <- header %>% paste0(markdown_line, '\nMedian:    ', mediandate %>% format(formatdates) %>% trimws(),
                                ' (~', percent(median_days / maxdate_days, round = 0), ')')
  }
  if (any(class(x) == 'POSIXlt')) {
    x <- x %>% format(formatdates)
  }

  nmax.set <- !missing(nmax)
  if (!nmax.set & is.null(nmax) & is.null(base::getOption("max.print.freq", default = NULL))) {
    # default for max print setting
    nmax <- 15
  } else if (is.null(nmax)) {
    nmax <- length(x)
  }

  if (nmax %in% c(0, Inf, NA, NULL)) {
    nmax <- length(x)
  }

  # create table with counts and percentages
  column_names <- c('Item', 'Count', 'Percent', 'Cum. Count', 'Cum. Percent', '(Factor Level)')
  column_names_df <- c('item', 'count', 'percent', 'cum_count', 'cum_percent', 'factor_level')

  if (any(class(x) == 'factor')) {
    df <- tibble(item = x,
                         fctlvl = x %>% as.integer()) %>%
      group_by(item, fctlvl)
    column_align <- c('l', 'r', 'r', 'r', 'r', 'r')
  } else {
    df <- tibble(item = x) %>%
      group_by(item)
    # strip factor lvl from col names
    column_names <- column_names[1:length(column_names) - 1]
    column_names_df <- column_names_df[1:length(column_names_df) - 1]
    column_align <- c(x_align, 'r', 'r', 'r', 'r')
  }
  df <- df %>% summarise(count = n())

  if (df$item %>% paste(collapse = ',') %like% '\033') {
    # remove escape char
    # see https://en.wikipedia.org/wiki/Escape_character#ASCII_escape_character
    df <- df %>% mutate(item = item %>% gsub('\033', ' ', ., fixed = TRUE))
  }

  # sort according to setting
  if (sort.count == TRUE) {
    df <- df %>% arrange(desc(count), item)
  } else {
    if (any(class(x) == 'factor')) {
      df <- df %>% arrange(fctlvl, item)
    } else {
      df <- df %>% arrange(item)
    }
  }

  df <- as.data.frame(df, stringsAsFactors = FALSE)

  df$percent <- df$count / base::sum(df$count, na.rm = TRUE)
  df$cum_count <- base::cumsum(df$count)
  df$cum_percent <- df$cum_count / base::sum(df$count, na.rm = TRUE)

  if (any(class(x) == 'factor')) {
    # put factor last
    df <- df %>% select(item, count, percent, cum_count, cum_percent, fctlvl)
  }

  colnames(df) <- column_names_df

  class(df) <- c('frequency_tbl', class(df))
  attr(df, 'package') <- 'AMR'

  if (markdown == TRUE) {
    tbl_format <- 'markdown'
  } else {
    tbl_format <- 'pandoc'
  }

  attr(df, 'opt') <- list(data = x.name,
                          vars = cols,
                          header = header,
                          row_names = row.names,
                          column_names = column_names,
                          column_align = column_align,
                          tbl_format = tbl_format,
                          nmax = nmax,
                          nmax.set = nmax.set)

  df
}

#' @rdname freq
#' @export
freq <- frequency_tbl

#' @rdname freq
#' @export
#' @importFrom dplyr top_n pull
top_freq <- function(f, n) {
  if (!'frequency_tbl' %in% class(f)) {
    stop('top_freq can only be applied to frequency tables', call. = FALSE)
  }
  if (!is.numeric(n) | length(n) != 1L) {
    stop('For top_freq, `nmax` must be a number of length 1', call. = FALSE)
  }
  top <- f %>% top_n(n, count)
  vect <- top %>% pull(item)
  names(vect) <- top %>% pull(count)
  if (length(vect) > abs(n)) {
    message("top_freq: selecting ", length(vect), " items instead of ", abs(n), ", because of ties")
  }
  vect
}

#' @rdname freq
#' @exportMethod print.frequency_tbl
#' @importFrom knitr kable
#' @importFrom dplyr n_distinct
#' @export
print.frequency_tbl <- function(x, nmax = getOption("max.print.freq", default = 15), ...) {

  opt <- attr(x, 'opt')

  if (length(opt$vars) == 0) {
    opt$vars <- NULL
  }

  if (!is.null(opt$data) & !is.null(opt$vars)) {
    title <- paste0("of `", paste0(opt$vars, collapse = "` and `"), "` from ", opt$data)
  } else if (!is.null(opt$data) & is.null(opt$vars)) {
    title <- paste("of", opt$data)
  } else if (is.null(opt$data) & !is.null(opt$vars)) {
    title <- paste0("of `", paste0(opt$vars, collapse = "` and `"), "`")
  } else {
    title <- ""
  }

  if (!missing(nmax)) {
    opt$nmax <- nmax
    opt$nmax.set <- TRUE
  }
  dots <- list(...)
  if ("markdown" %in% names(dots)) {
    if (dots$markdown == TRUE) {
      opt$tbl_format <- "markdown"
    } else {
      opt$tbl_format <- "pandoc"
    }
  }

  cat("Frequency table", title, "\n")

  if (!is.null(opt$header)) {
    cat(opt$header)
  }

  if (NROW(x) == 0) {
    cat('\n\nNo observations.\n')
    return(invisible())
  }

  if (all(x$count == 1)) {
    warning('All observations are unique.', call. = FALSE)
  }

  # save old NA setting for kable
  opt.old <- options()$knitr.kable.NA
  options(knitr.kable.NA = "<NA>")

  if (nrow(x) > opt$nmax & opt$tbl_format != "markdown") {

    x.rows <- nrow(x)
    x.unprinted <- base::sum(x[(opt$nmax + 1):nrow(x), 'count'], na.rm = TRUE)
    x.printed <- base::sum(x$count) - x.unprinted

    if (opt$nmax.set == TRUE) {
      nmax <- opt$nmax
    } else {
      nmax <- getOption("max.print.freq", default = 15)
    }

    x <- x[1:nmax,]

    if (opt$nmax.set == TRUE) {
      footer <- paste('[ reached `nmax = ', opt$nmax, '`', sep = '')
    } else {
      footer <- '[ reached getOption("max.print.freq")'
    }
    footer <- paste(footer,
                    ' -- omitted ',
                    format(x.rows - opt$nmax),
                    ' entries, n = ',
                    format(x.unprinted),
                    ' (',
                    (x.unprinted / (x.unprinted + x.printed)) %>% percent(force_zero = TRUE),
                    ') ]\n', sep = '')
  } else {
    footer <- NULL
  }

  if (any(class(x$item) %in% c('double', 'integer', 'numeric', 'raw', 'single'))) {
    x$item <- format(x$item)
  }
  x$count <- format(x$count)
  x$percent <- percent(x$percent, force_zero = TRUE)
  x$cum_count <- format(x$cum_count)
  x$cum_percent <- percent(x$cum_percent, force_zero = TRUE)

  print(
    knitr::kable(x,
                 format = opt$tbl_format,
                 row.names = opt$row_names,
                 col.names = opt$column_names,
                 align = opt$column_align,
                 padding = 1)
  )

  if (!is.null(footer)) {
    cat(footer)
  }

  cat('\n')

  # reset old kable setting
  options(knitr.kable.NA = opt.old)
  return(invisible())

}

#' @noRd
#' @exportMethod as.data.frame.frequency_tbl
#' @export
as.data.frame.frequency_tbl <- function(x, ...) {
  attr(x, 'package') <- NULL
  attr(x, 'opt') <- NULL
  as.data.frame.data.frame(x, ...)
}

#' @noRd
#' @exportMethod as_tibble.frequency_tbl
#' @export
#' @importFrom dplyr as_tibble
as_tibble.frequency_tbl <- function(x, validate = TRUE, ..., rownames = NA) {
  attr(x, 'package') <- NULL
  attr(x, 'opt') <- NULL
  as_tibble(x = as.data.frame(x), validate = validate, ..., rownames = rownames)
}

#' @noRd
#' @exportMethod hist.frequency_tbl
#' @export
#' @importFrom graphics hist
hist.frequency_tbl <- function(x, ...) {
  opt <- attr(x, 'opt')
  if (!is.null(opt$vars)) {
    title <- opt$vars
  } else {
    title <- ""
  }
  hist(as.vector(x), main = paste("Histogram of", title), xlab = title, ...)
}

#' @noRd
#' @exportMethod plot.frequency_tbl
#' @export
plot.frequency_tbl <- function(x, y, ...) {
  opt <- attr(x, 'opt')
  if (!is.null(opt$vars)) {
    title <- opt$vars
  } else {
    title <- ""
  }
  plot(x = x$item, y = x$count, ylab = "Count", xlab = title, ...)
}

#' @noRd
#' @exportMethod as.vector.frequency_tbl
#' @export
as.vector.frequency_tbl <- function(x, mode = "any") {
  as.vector(rep(x$item, x$count), mode = mode)
}

#' @noRd
#' @exportMethod format.frequency_tbl
#' @export
format.frequency_tbl <- function(x, digits = 1, ...) {
  opt <- attr(x, 'opt')
  if (opt$nmax.set == TRUE) {
    nmax <- opt$nmax
  } else {
    nmax <- getOption("max.print.freq", default = 15)
  }

  x <- x[1:nmax,]
  x$percent <- percent(x$percent, round = digits, force_zero = TRUE)
  x$cum_percent <- percent(x$cum_percent, round = digits, force_zero = TRUE)
  base::format.data.frame(x, ...)
}
