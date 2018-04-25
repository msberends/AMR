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

#' Class 'rsi'
#'
#' This transforms a vector to a new class \code{rsi}, which is an ordered factor with levels \code{S < I < R}. Invalid antimicrobial interpretations will be translated as \code{NA} with a warning.
#' @rdname as.rsi
#' @param x vector
#' @return Ordered factor with new class \code{rsi} and new attributes \code{package} and \code{package.version}
#' @export
#' @importFrom dplyr %>%
#' @importFrom utils packageDescription
#' @examples
#' rsi_data <- as.rsi(c(rep("S", 474), rep("I", 36), rep("R", 370)))
#' rsi_data <- as.rsi(c(rep("S", 474), rep("I", 36), rep("R", 370), "A", "B", "C"))
#' is.rsi(rsi_data)
#'
#' plot(rsi_data)    # for percentages
#' barplot(rsi_data) # for frequencies
as.rsi <- function(x) {
  if (is.rsi(x)) {
    x
  } else {

    x <- x %>% unlist()
    x.bak <- x

    na_before <- x[is.na(x) | x == ''] %>% length()
    # remove all spaces
    x <- gsub(' {2,55}', '', x)
    # disallow more than 3 characters
    x[nchar(x) > 3] <- NA
    # remove all invalid characters
    x <- gsub('[^RSI]+', '', x %>% toupper())
    # needed for UMCG in cases of "S;S" but also "S;I"; the latter will be NA:
    x <- gsub('^S+$', 'S', x)
    x <- gsub('^I+$', 'I', x)
    x <- gsub('^R+$', 'R', x)
    x[!x %in% c('S', 'I', 'R')] <- NA
    na_after <- x[is.na(x) | x == ''] %>% length()

    if (na_before != na_after) {
      list_missing <- x.bak[is.na(x) & !is.na(x.bak) & x.bak != ''] %>%
        unique() %>%
        sort()
      list_missing <- paste0('"', list_missing , '"', collapse = ", ")
      warning(na_after - na_before, ' results truncated (',
              round(((na_after - na_before) / length(x)) * 100),
              '%) that were invalid antimicrobial interpretations: ',
              list_missing, call. = FALSE)
    }

    x <- x %>% toupper() %>% factor(levels = c("S", "I", "R"), ordered = TRUE)
    class(x) <- c('rsi', 'ordered', 'factor')
    attr(x, 'package') <- 'AMR'
    attr(x, 'package.version') <- packageDescription('AMR')$Version
    x
  }
}

#' @rdname as.rsi
#' @export
#' @importFrom dplyr %>%
is.rsi <- function(x) {
  class(x) %>% identical(c('rsi', 'ordered', 'factor'))
}

#' @exportMethod print.rsi
#' @export
#' @importFrom dplyr %>%
#' @noRd
print.rsi <- function(x, ...) {
  n_total <- x %>% length()
  x <- x[!is.na(x)]
  n <- x %>% length()
  S <- x[x == 'S'] %>% length()
  I <- x[x == 'I'] %>% length()
  R <- x[x == 'R'] %>% length()
  IR <- x[x %in% c('I', 'R')] %>% length()
  cat("Class 'rsi'\n")
  cat(n, " results (missing: ", n_total - n, ' = ', percent((n_total - n) / n, force_zero = TRUE), ')\n', sep = "")
  cat('\n')
  cat('Sum of S:   ', S, ' (', percent(S / n, force_zero = TRUE), ')\n', sep = "")
  cat('Sum of IR:  ', IR, ' (', percent(IR / n, force_zero = TRUE), ')\n', sep = "")
  cat('- Sum of R: ', R, ' (', percent(R / n, force_zero = TRUE), ')\n', sep = "")
  cat('- Sum of I: ', I, ' (', percent(I / n, force_zero = TRUE), ')\n', sep = "")
}

#' @exportMethod summary.rsi
#' @export
#' @importFrom dplyr %>%
#' @noRd
summary.rsi <- function(object, ...) {
  x <- object
  n_total <- x %>% length()
  x <- x[!is.na(x)]
  n <- x %>% length()
  S <- x[x == 'S'] %>% length()
  I <- x[x == 'I'] %>% length()
  R <- x[x == 'R'] %>% length()
  IR <- x[x %in% c('I', 'R')] %>% length()
  lst <- c('rsi', n_total - n, S, IR, R, I)
  names(lst) <- c("Mode", "<NA>", "Sum S", "Sum IR", "Sum R", "Sum I")
  lst
}

#' @exportMethod plot.rsi
#' @export
#' @importFrom dplyr %>% group_by summarise filter mutate if_else n_distinct
#' @importFrom graphics plot text
#' @noRd
plot.rsi <- function(x, ...) {
  x_name <- deparse(substitute(x))

  data <- data.frame(x = x,
                     y = 1,
                     stringsAsFactors = TRUE) %>%
    group_by(x) %>%
    summarise(n = sum(y)) %>%
    filter(!is.na(x)) %>%
    mutate(s = round((n / sum(n)) * 100, 1))
  data$x <- factor(data$x, levels = c('S', 'I', 'R'), ordered = TRUE)

  ymax <- if_else(max(data$s) > 95, 105, 100)

  plot(x = data$x,
       y = data$s,
       lwd = 2,
       col = c('green', 'orange', 'red'),
       ylim = c(0, ymax),
       ylab = 'Percentage',
       xlab = 'Antimicrobial Interpretation',
       main = paste('Susceptibilty Analysis of', x_name),
       axes = FALSE,
       ...)
  # x axis
  axis(side = 1, at = 1:n_distinct(data$x), labels = levels(data$x), lwd = 0)
  # y axis, 0-100%
  axis(side = 2, at = seq(0, 100, 5))

  text(x = data$x,
       y = data$s + 4,
       labels = paste0(data$s, '% (n = ', data$n, ')'))
}


#' @exportMethod barplot.rsi
#' @export
#' @importFrom dplyr %>% group_by summarise filter mutate if_else n_distinct
#' @importFrom graphics barplot axis
#' @noRd
barplot.rsi <- function(height, ...) {
  x <- height
  x_name <- deparse(substitute(height))

  data <- data.frame(rsi = x, cnt = 1) %>%
    group_by(rsi) %>%
    summarise(cnt = sum(cnt)) %>%
    droplevels()

  barplot(table(x),
          col = c('green3', 'orange2', 'red3'),
          xlab = 'Antimicrobial Interpretation',
          main = paste('Susceptibilty Analysis of', x_name),
          ylab = 'Frequency',
          axes = FALSE,
          ...)
  # y axis, 0-100%
  axis(side = 2, at = seq(0, max(data$cnt) + max(data$cnt) * 1.1, by = 25))
}

#' Class 'mic'
#'
#' This transforms a vector to a new class\code{mic}, which is an ordered factor with valid MIC values as levels. Invalid MIC values will be translated as \code{NA} with a warning.
#' @rdname as.mic
#' @param x vector
#' @param na.rm a logical indicating whether missing values should be removed
#' @return Ordered factor with new class \code{mic} and new attributes \code{package} and \code{package.version}
#' @export
#' @importFrom dplyr %>%
#' @importFrom utils packageDescription
#' @examples
#' mic_data <- as.mic(c(">=32", "1.0", "1", "1.00", 8, "<=0.128", "8", "16", "16"))
#' is.mic(mic_data)
#'
#' plot(mic_data)
#' barplot(mic_data)
as.mic <- function(x, na.rm = FALSE) {
  if (is.mic(x)) {
    x
  } else {
    x <- x %>% unlist()
    if (na.rm == TRUE) {
      x <- x[!is.na(x)]
    }
    x.bak <- x

    # comma to dot
    x <- gsub(',', '.', x, fixed = TRUE)
    # starting dots must start with 0
    x <- gsub('^[.]', '0.', x)
    # <=0.2560.512 should be 0.512
    x <- gsub('.*[.].*[.]', '0.', x)
    # remove ending .0
    x <- gsub('[.]0$', '', x)
    # remove all after last digit
    x <- gsub('[^0-9]$', '', x)
    # remove last zeroes
    x <- gsub('[.]?0+$', '', x)

    lvls <- c("<0.002", "<=0.002", "0.002", ">=0.002", ">0.002",
              "<0.003", "<=0.003", "0.003", ">=0.003", ">0.003",
              "<0.004", "<=0.004", "0.004", ">=0.004", ">0.004",
              "<0.006", "<=0.006", "0.006", ">=0.006", ">0.006",
              "<0.008", "<=0.008", "0.008", ">=0.008", ">0.008",
              "<0.012", "<=0.012", "0.012", ">=0.012", ">0.012",
              "<0.016", "<=0.016", "0.016", ">=0.016", ">0.016",
              "<0.023", "<=0.023", "0.023", ">=0.023", ">0.023",
              "<0.025", "<=0.025", "0.025", ">=0.025", ">0.025",
              "<0.03", "<=0.03", "0.03", ">=0.03", ">0.03",
              "<0.032", "<=0.032", "0.032", ">=0.032", ">0.032",
              "<0.047", "<=0.047", "0.047", ">=0.047", ">0.047",
              "<0.05", "<=0.05", "0.05", ">=0.05", ">0.05",
              "<0.06", "<=0.06", "0.06", ">=0.06", ">0.06",
              "<0.0625", "<=0.0625", "0.0625", ">=0.0625", ">0.0625",
              "<0.063", "<=0.063", "0.063", ">=0.063", ">0.063",
              "<0.064", "<=0.064", "0.064", ">=0.064", ">0.064",
              "<0.09", "<=0.09", "0.09", ">=0.09", ">0.09",
              "<0.094", "<=0.094", "0.094", ">=0.094", ">0.094",
              "<0.12", "<=0.12", "0.12", ">=0.12", ">0.12",
              "<0.125", "<=0.125", "0.125", ">=0.125", ">0.125",
              "<0.128", "<=0.128", "0.128", ">=0.128", ">0.128",
              "<0.16", "<=0.16", "0.16", ">=0.16", ">0.16",
              "<0.19", "<=0.19", "0.19", ">=0.19", ">0.19",
              "<0.25", "<=0.25", "0.25", ">=0.25", ">0.25",
              "<0.256", "<=0.256", "0.256", ">=0.256", ">0.256",
              "<0.32", "<=0.32", "0.32", ">=0.32", ">0.32",
              "<0.38", "<=0.38", "0.38", ">=0.38", ">0.38",
              "<0.5", "<=0.5", "0.5", ">=0.5", ">0.5",
              "<0.512", "<=0.512", "0.512", ">=0.512", ">0.512",
              "<0.64", "<=0.64", "0.64", ">=0.64", ">0.64",
              "<0.75", "<=0.75", "0.75", ">=0.75", ">0.75",
              "<1", "<=1", "1", ">=1", ">1",
              "<1.5", "<=1.5", "1.5", ">=1.5", ">1.5",
              "<2", "<=2", "2", ">=2", ">2",
              "<3", "<=3", "3", ">=3", ">3",
              "<4", "<=4", "4", ">=4", ">4",
              "<6", "<=6", "6", ">=6", ">6",
              "<8", "<=8", "8", ">=8", ">8",
              "<10", "<=10", "10", ">=10", ">10",
              "<12", "<=12", "12", ">=12", ">12",
              "<16", "<=16", "16", ">=16", ">16",
              "<20", "<=20", "20", ">=20", ">20",
              "<24", "<=24", "24", ">=24", ">24",
              "<32", "<=32", "32", ">=32", ">32",
              "<40", "<=40", "40", ">=40", ">40",
              "<48", "<=48", "48", ">=48", ">48",
              "<64", "<=64", "64", ">=64", ">64",
              "<80", "<=80", "80", ">=80", ">80",
              "<96", "<=96", "96", ">=96", ">96",
              "<128", "<=128", "128", ">=128", ">128",
              "<160", "<=160", "160", ">=160", ">160",
              "<256", "<=256", "256", ">=256", ">256",
              "<320", "<=320", "320", ">=320", ">320",
              "<512", "<=512", "512", ">=512", ">512",
              "<1024", "<=1024", "1024", ">=1024", ">1024")
    x <- x %>% as.character()

    na_before <- x[is.na(x) | x == ''] %>% length()
    x[!x %in% lvls] <- NA
    na_after <- x[is.na(x) | x == ''] %>% length()

    if (na_before != na_after) {
      list_missing <- x.bak[is.na(x) & !is.na(x.bak) & x.bak != ''] %>%
        unique() %>%
        sort()
      list_missing <- paste0('"', list_missing , '"', collapse = ", ")
      warning(na_after - na_before, ' results truncated (',
              round(((na_after - na_before) / length(x)) * 100),
              '%) that were invalid MICs: ',
              list_missing, call. = FALSE)
    }

    x <- factor(x = x,
                levels = lvls,
                ordered = TRUE)
    class(x) <- c('mic', 'ordered', 'factor')
    attr(x, 'package') <- 'AMR'
    attr(x, 'package.version') <- packageDescription('AMR')$Version
    x
  }
}

#' @rdname as.mic
#' @export
#' @importFrom dplyr %>%
is.mic <- function(x) {
  class(x) %>% identical(c('mic', 'ordered', 'factor'))
}

#' @exportMethod as.double.mic
#' @export
#' @noRd
as.double.mic <- function(x, ...) {
  as.double(gsub('(<|=|>)+', '', as.character(x)))
}

#' @exportMethod as.integer.mic
#' @export
#' @noRd
as.integer.mic <- function(x, ...) {
  as.integer(gsub('(<|=|>)+', '', as.character(x)))
}

#' @exportMethod as.numeric.mic
#' @export
#' @noRd
as.numeric.mic <- function(x, ...) {
  as.numeric(gsub('(<|=|>)+', '', as.character(x)))
}

#' @exportMethod print.mic
#' @export
#' @importFrom dplyr %>% tibble group_by summarise pull
#' @noRd
print.mic <- function(x, ...) {
  n_total <- x %>% length()
  x <- x[!is.na(x)]
  n <- x %>% length()
  cat("Class 'mic': ", n, " isolates\n", sep = '')
  cat('\n')
  cat('<NA> ', n_total - n, '\n')
  cat('\n')
  tbl <- tibble(x = x, y = 1) %>% group_by(x) %>% summarise(y = sum(y))
  cnt <- tbl %>% pull(y)
  names(cnt) <- tbl %>% pull(x)
  print(cnt)
}

#' @exportMethod summary.mic
#' @export
#' @importFrom dplyr %>%
#' @noRd
summary.mic <- function(object, ...) {
  x <- object
  n_total <- x %>% length()
  x <- x[!is.na(x)]
  n <- x %>% length()
  lst <- c('mic',
           n_total - n,
           sort(x)[1] %>% as.character(),
           sort(x)[n] %>% as.character())
  names(lst) <- c("Mode", "<NA>", "Min.", "Max.")
  lst
}

#' @exportMethod plot.mic
#' @export
#' @importFrom dplyr %>% group_by summarise
#' @importFrom graphics plot text
#' @noRd
plot.mic <- function(x, ...) {
  x_name <- deparse(substitute(x))
  create_barplot_mic(x, x_name, ...)
}

#' @exportMethod barplot.mic
#' @export
#' @importFrom dplyr %>% group_by summarise
#' @importFrom graphics barplot axis
#' @noRd
barplot.mic <- function(height, ...) {
  x_name <- deparse(substitute(height))
  create_barplot_mic(height, x_name, ...)
}

#' @importFrom graphics barplot axis
create_barplot_mic <- function(x, x_name, ...) {
  data <- data.frame(mic = x, cnt = 1) %>%
    group_by(mic) %>%
    summarise(cnt = sum(cnt)) %>%
    droplevels()
  barplot(table(droplevels(x)),
          ylab = 'Frequency',
          xlab = 'MIC value',
          main = paste('MIC values of', x_name),
          axes = FALSE,
          ...)
  axis(2, seq(0, max(data$cnt)))
}
