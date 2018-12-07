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

#' Class 'mic'
#'
#' This transforms a vector to a new class \code{mic}, which is an ordered factor with valid MIC values as levels. Invalid MIC values will be translated as \code{NA} with a warning.
#' @rdname as.mic
#' @param x vector
#' @param na.rm a logical indicating whether missing values should be removed
#' @return Ordered factor with new class \code{mic}
#' @keywords mic
#' @export
#' @importFrom dplyr %>%
#' @seealso \code{\link{as.rsi}}
#' @examples
#' mic_data <- as.mic(c(">=32", "1.0", "1", "1.00", 8, "<=0.128", "8", "16", "16"))
#' is.mic(mic_data)
#'
#' # this can also coerce combined MIC/RSI values:
#' as.mic("<=0.002; S") # will return <=0.002
#'
#' plot(mic_data)
#' barplot(mic_data)
#' freq(mic_data)
as.mic <- function(x, na.rm = FALSE) {
  if (is.mic(x)) {
    x
  } else {
    x <- x %>% unlist()
    if (na.rm == TRUE) {
      x <- x[!is.na(x)]
    }
    x.bak <- x

    # comma to period
    x <- gsub(',', '.', x, fixed = TRUE)
    # remove space between operator and number ("<= 0.002" -> "<=0.002")
    x <- gsub('(<|=|>) +', '\\1', x)
    # starting dots must start with 0
    x <- gsub('^[.]+', '0.', x)
    # <=0.2560.512 should be 0.512
    x <- gsub('.*[.].*[.]', '0.', x)
    # remove ending .0
    x <- gsub('[.]+0$', '', x)
    # remove all after last digit
    x <- gsub('[^0-9]+$', '', x)
    # remove last zeroes
    x <- gsub('([.].?)0+$', '\\1', x)
    # force to be character
    x <- as.character(x)

    # previously unempty values now empty - should return a warning later on
    x[x.bak != "" & x == ""] <- "invalid"

    # these are alllowed MIC values and will become factor levels
    lvls <- c("<0.002", "<=0.002", "0.002", ">=0.002", ">0.002",
              "<0.003", "<=0.003", "0.003", ">=0.003", ">0.003",
              "<0.004", "<=0.004", "0.004", ">=0.004", ">0.004",
              "<0.006", "<=0.006", "0.006", ">=0.006", ">0.006",
              "<0.008", "<=0.008", "0.008", ">=0.008", ">0.008",
              "<0.012", "<=0.012", "0.012", ">=0.012", ">0.012",
              "<0.0125", "<=0.0125", "0.0125", ">=0.0125", ">0.0125",
              "<0.016", "<=0.016", "0.016", ">=0.016", ">0.016",
              "<0.023", "<=0.023", "0.023", ">=0.023", ">0.023",
              "<0.025", "<=0.025", "0.025", ">=0.025", ">0.025",
              "<0.03", "<=0.03", "0.03", ">=0.03", ">0.03",
              "<0.032", "<=0.032", "0.032", ">=0.032", ">0.032",
              "<0.047", "<=0.047", "0.047", ">=0.047", ">0.047",
              "<0.05", "<=0.05", "0.05", ">=0.05", ">0.05",
              "<0.054", "<=0.054", "0.054", ">=0.054", ">0.054",
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
              "<0.23", "<=0.23", "0.23", ">=0.23", ">0.23",
              "<0.25", "<=0.25", "0.25", ">=0.25", ">0.25",
              "<0.256", "<=0.256", "0.256", ">=0.256", ">0.256",
              "<0.28", "<=0.28", "0.28", ">=0.28", ">0.28",
              "<0.3", "<=0.3", "0.3", ">=0.3", ">0.3",
              "<0.32", "<=0.32", "0.32", ">=0.32", ">0.32",
              "<0.36", "<=0.36", "0.36", ">=0.36", ">0.36",
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
              "<5", "<=5", "5", ">=5", ">5",
              "<6", "<=6", "6", ">=6", ">6",
              "<7", "<=7", "7", ">=7", ">7",
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

    x <- factor(x, levels = lvls, ordered = TRUE)
    class(x) <- c('mic', 'ordered', 'factor')
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
  cat("Class 'mic'\n")
  print(as.character(x), quote = FALSE)
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
  c(
    "Class" = 'mic',
    "<NA>" = n_total - n,
    "Min." = sort(x)[1] %>% as.character(),
    "Max." = sort(x)[n] %>% as.character()
  )
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
#' @importFrom graphics barplot axis
#' @noRd
barplot.mic <- function(height, ...) {
  x_name <- deparse(substitute(height))
  create_barplot_mic(height, x_name, ...)
}

#' @importFrom graphics barplot axis
#' @importFrom dplyr %>% group_by summarise
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
