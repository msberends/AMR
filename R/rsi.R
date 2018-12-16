# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Analysis                              #
#                                                                      #
# AUTHORS                                                              #
# Berends MS (m.s.berends@umcg.nl), Luz CF (c.f.luz@umcg.nl)           #
#                                                                      #
# LICENCE                                                              #
# This package is free software; you can redistribute it and/or modify #
# it under the terms of the GNU General Public License version 2.0,    #
# as published by the Free Software Foundation.                        #
#                                                                      #
# This R package is distributed in the hope that it will be useful,    #
# but WITHOUT ANY WARRANTY; without even the implied warranty of       #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        #
# GNU General Public License version 2.0 for more details.             #
# ==================================================================== #

#' Class 'rsi'
#'
#' This transforms a vector to a new class \code{rsi}, which is an ordered factor with levels \code{S < I < R}. Invalid antimicrobial interpretations will be translated as \code{NA} with a warning.
#' @rdname as.rsi
#' @param x vector
#' @details The function \code{is.rsi.eligible} returns \code{TRUE} when a columns contains only valid antimicrobial interpretations (S and/or I and/or R), and \code{FALSE} otherwise.
#' @return Ordered factor with new class \code{rsi}
#' @keywords rsi
#' @export
#' @importFrom dplyr %>%
#' @seealso \code{\link{as.mic}}
#' @examples
#' rsi_data <- as.rsi(c(rep("S", 474), rep("I", 36), rep("R", 370)))
#' rsi_data <- as.rsi(c(rep("S", 474), rep("I", 36), rep("R", 370), "A", "B", "C"))
#' is.rsi(rsi_data)
#'
#' # this can also coerce combined MIC/RSI values:
#' as.rsi("<= 0.002; S") # will return S
#'
#' plot(rsi_data)    # for percentages
#' barplot(rsi_data) # for frequencies
#' freq(rsi_data)    # frequency table with informative header
#'
#' # using dplyr's mutate
#' library(dplyr)
#' septic_patients %>%
#'   mutate_at(vars(peni:rifa), as.rsi)
#'
#' # fastest way to transform all columns with already valid AB results to class `rsi`:
#' septic_patients %>%
#'   mutate_if(is.rsi.eligible,
#'             as.rsi)
as.rsi <- function(x) {
  if (is.rsi(x)) {
    x
  } else if (identical(levels(x), c("S", "I", "R"))) {
    structure(x, class = c('rsi', 'ordered', 'factor'))
  } else {

    x <- x %>% unlist()
    x.bak <- x

    na_before <- x[is.na(x) | x == ''] %>% length()
    # remove all spaces
    x <- gsub(' +', '', x)
    # remove all MIC-like values: numbers, operators and periods
    x <- gsub('[0-9.,;:<=>]+', '', x)
    # disallow more than 3 characters
    x[nchar(x) > 3] <- NA
    # set to capitals
    x <- toupper(x)
    # remove all invalid characters
    x <- gsub('[^RSI]+', '', x)
    # in cases of "S;S" keep S, but in case of "S;I" make it NA
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

    x <- factor(x, levels = c("S", "I", "R"), ordered = TRUE)
    class(x) <- c('rsi', 'ordered', 'factor')
    x
  }
}

#' @rdname as.rsi
#' @export
#' @importFrom dplyr %>%
is.rsi <- function(x) {
  class(x) %>% identical(c('rsi', 'ordered', 'factor'))
}

#' @rdname as.rsi
#' @export
#' @importFrom dplyr %>%
is.rsi.eligible <- function(x) {
  if (is.logical(x)
      | is.numeric(x)
      | is.mo(x)
      | identical(class(x), "Date")
      | is.rsi(x)) {
    # no transformation needed
    FALSE
  } else {
    # check all but a-z
    y <- unique(gsub("[^RSIrsi]+", "", unique(x)))
    !all(y %in% c("", NA_character_)) &
      all(y %in% c("R", "I", "S", "", NA_character_)) &
      max(nchar(as.character(x)), na.rm = TRUE) < 8
  }
}

#' @exportMethod print.rsi
#' @export
#' @importFrom dplyr %>%
#' @noRd
print.rsi <- function(x, ...) {
  cat("Class 'rsi'\n")
  print(as.character(x), quote = FALSE)
}

#' @exportMethod summary.rsi
#' @export
#' @noRd
summary.rsi <- function(object, ...) {
  x <- object
  c(
    "Class" = 'rsi',
    "<NA>" = sum(is.na(x)),
    "Sum S" = sum(x == "S", na.rm = TRUE),
    "Sum IR" = sum(x %in% c("I", "R"), na.rm = TRUE),
    "-Sum R" = sum(x == "R", na.rm = TRUE),
    "-Sum I" = sum(x == "I", na.rm = TRUE)
  )
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
       main = paste('Susceptibility Analysis of', x_name),
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
          main = paste('Susceptibility Analysis of', x_name),
          ylab = 'Frequency',
          axes = FALSE,
          ...)
  # y axis, 0-100%
  axis(side = 2, at = seq(0, max(data$cnt) + max(data$cnt) * 1.1, by = 25))
}
