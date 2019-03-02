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

#' Class 'rsi'
#'
#' This transforms a vector to a new class \code{rsi}, which is an ordered factor with levels \code{S < I < R}. Invalid antimicrobial interpretations will be translated as \code{NA} with a warning.
#' @rdname as.rsi
#' @param x vector
#' @param threshold maximum fraction of \code{x} that is allowed to fail transformation, see Examples
#' @details The function \code{is.rsi.eligible} returns \code{TRUE} when a columns contains only valid antimicrobial interpretations (S and/or I and/or R), and \code{FALSE} otherwise.
#' @return Ordered factor with new class \code{rsi}
#' @keywords rsi
#' @export
#' @importFrom dplyr %>%
#' @seealso \code{\link{as.mic}}
#' @inheritSection AMR Read more on our website!
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
#'
#' # fastest way to transform all columns with already valid AB results to class `rsi`:
#' septic_patients %>%
#'   mutate_if(is.rsi.eligible,
#'             as.rsi)
#'
#' # default threshold of `is.rsi.eligible` is 5%.
#' is.rsi.eligible(WHONET$`First name`) # fails, >80% is invalid
#' is.rsi.eligible(WHONET$`First name`, threshold = 0.9) # succeeds
as.rsi <- function(x) {
  if (is.rsi(x)) {
    x
  } else if (identical(levels(x), c("S", "I", "R"))) {
    structure(x, class = c('rsi', 'ordered', 'factor'))
  } else {
    if (mic_like(x) > 0.5) {
      warning("`as.rsi` is intended to clean antimicrobial interpretations - not to interpret MIC values.", call. = FALSE)
    }

    x <- x %>% unlist()
    x.bak <- x

    na_before <- x[is.na(x) | x == ''] %>% length()
    # remove all spaces
    x <- gsub(' +', '', x)
    # remove all MIC-like values: numbers, operators and periods
    x <- gsub('[0-9.,;:<=>]+', '', x)
    # remove everything between brackets, and 'high' and 'low'
    x <- gsub("([(].*[)])", "", x)
    x <- gsub("(high|low)", "", x, ignore.case = TRUE)
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

mic_like <- function(x) {
  mic <- x %>%
    gsub("[^0-9.,]+", "", .) %>%
    unique()
  mic_valid <- suppressWarnings(as.mic(mic))
  sum(!is.na(mic_valid)) / length(mic)
}

#' @rdname as.rsi
#' @export
is.rsi <- function(x) {
  identical(class(x),
            c('rsi', 'ordered', 'factor'))
}

#' @rdname as.rsi
#' @export
is.rsi.eligible <- function(x, threshold = 0.05) {
  if (NCOL(x) > 1) {
    stop('`x` must be a one-dimensional vector.')
  }
  if (any(c("logical",
            "numeric",
            "integer",
            "mo",
            "Date",
            "POSIXct",
            "rsi",
            "raw",
            "hms")
          %in% class(x))) {
    # no transformation needed
    FALSE
  } else {
    x <- x[!is.na(x) & !is.null(x) & !identical(x, "")]
    if (length(x) == 0) {
      return(FALSE)
    }
    checked <- suppressWarnings(as.rsi(x))
    outcome <- sum(is.na(checked)) / length(x)
    outcome <= threshold
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

#' @exportMethod droplevels.rsi
#' @export
#' @noRd
droplevels.rsi <- function(x, exclude = if(anyNA(levels(x))) NULL else NA, ...) {
  x <- droplevels.factor(x, exclude = exclude, ...)
  class(x) <- c('rsi', 'ordered', 'factor')
  x
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

  suppressWarnings(
    data <- data.frame(x = x,
                       y = 1,
                       stringsAsFactors = TRUE) %>%
      group_by(x) %>%
      summarise(n = sum(y)) %>%
      filter(!is.na(x)) %>%
      mutate(s = round((n / sum(n)) * 100, 1))
  )

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

  suppressWarnings(
    data <- data.frame(rsi = x, cnt = 1) %>%
      group_by(rsi) %>%
      summarise(cnt = sum(cnt)) %>%
      droplevels()
  )

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
