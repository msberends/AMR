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

#' Age in years of individuals
#'
#' Calculates age in years based on a reference date, which is the sytem time at default.
#' @param x date(s), will be coerced with \code{\link{as.POSIXlt}}
#' @param reference reference date(s) (defaults to today), will be coerced with \code{\link{as.POSIXlt}}
#' @return Integer (no decimals)
#' @seealso \code{\link{age_groups}} to splits age into groups
#' @importFrom dplyr if_else
#' @inheritSection AMR Read more on our website!
#' @export
age <- function(x, reference = Sys.Date()) {
  if (length(x) != length(reference)) {
    if (length(reference) == 1) {
      reference <- rep(reference, length(x))
    } else {
      stop("`x` and `reference` must be of same length, or `reference` must be of length 1.")
    }
  }
  x <- base::as.POSIXlt(x)
  reference <- base::as.POSIXlt(reference)
  if (any(reference < x)) {
    stop("`reference` cannot be lower (older) than `x`.")
  }
  years_gap <- reference$year - x$year
  # from https://stackoverflow.com/a/25450756/4575331
  ages <- if_else(reference$mon < x$mon | (reference$mon == x$mon & reference$mday < x$mday),
         as.integer(years_gap - 1),
         as.integer(years_gap))
  if (any(ages > 120)) {
    warning("Some ages are >120.")
  }
  ages
}

#' Split ages into age groups
#'
#' Split ages into age groups defined by the \code{split} parameter. This allows for easier demographic (antimicrobial resistance) analysis.
#' @param x age, e.g. calculated with \code{\link{age}}
#' @param split_at values to split \code{x} at, defaults to age groups 0-11, 12-24, 26-54, 55-74 and 75+. See Details.
#' @details To split ages, the input can be:
#' \itemize{
#'   \item{A numeric vector. A vector of \code{c(10, 20)} will split on 0-9, 10-19 and 20+. A value of only \code{50} will split on 0-49 and 50+.
#'         The default is to split on young children (0-11), youth (12-24), young adults (26-54), middle-aged adults (55-74) and elderly (75+).}
#'   \item{A character:}
#'     \itemize{
#'       \item{\code{"children"}, equivalent of: \code{c(0, 1, 2, 4, 6, 13, 18)}. This will split on 0, 1, 2-3, 4-5, 6-12, 13-17 and 18+.}
#'       \item{\code{"elderly"} or \code{"seniors"}, equivalent of: \code{c(65, 75, 85, 95)}. This will split on 0-64, 65-74, 75-84, 85-94 and 95+.}
#'       \item{\code{"fives"}, equivalent of: \code{1:20 * 5}. This will split on 0-4, 5-9, 10-14, 15-19 and so forth.}
#'       \item{\code{"tens"}, equivalent of: \code{1:10 * 10}. This will split on 0-9, 10-19, 20-29 and so forth.}
#'     }
#' }
#' @keywords age_group age
#' @return Ordered \code{\link{factor}}
#' @seealso \code{\link{age}} to determine ages based on one or more reference dates
#' @export
#' @inheritSection AMR Read more on our website!
#' @examples
#' ages <- c(3, 8, 16, 54, 31, 76, 101, 43, 21)
#'
#' # split into 0-49 and 50+
#' age_groups(ages, 50)
#'
#' # split into 0-19, 20-49 and 50+
#' age_groups(ages, c(20, 50))
#'
#' # split into groups of ten years
#' age_groups(ages, 1:10 * 10)
#' age_groups(ages, split_at = "tens")
#'
#' # split into groups of five years
#' age_groups(ages, 1:20 * 5)
#' age_groups(ages, split_at = "fives")
#'
#' # split specifically for children
#' age_groups(ages, "children")
#' # same:
#' age_groups(ages, c(1, 2, 4, 6, 13, 17))
#'
#' # resistance of ciprofloxacine per age group
#' library(dplyr)
#' septic_patients %>%
#'   mutate(first_isolate = first_isolate(.)) %>%
#'   filter(first_isolate == TRUE,
#'          mo == as.mo("E. coli")) %>%
#'   group_by(age_group = age_groups(age)) %>%
#'   select(age_group,
#'          cipr) %>%
#'   ggplot_rsi(x = "age_group")
age_groups <- function(x, split_at = c(12, 25, 55, 75)) {
  if (is.character(split_at)) {
    split_at <- split_at[1L]
    if (split_at %like% "^child") {
      split_at <- c(0, 1, 2, 4, 6, 13, 18)
    } else if (split_at %like% "^(elder|senior)") {
      split_at <-  c(65, 75, 85, 95)
    } else if (split_at %like% "^five") {
      split_at <- 1:20 * 5
    } else if (split_at %like% "^ten") {
      split_at <- 1:10 * 10
    }
  }
  split_at <- as.integer(split_at)
  if (!is.numeric(x) | !is.numeric(split_at)) {
    stop("`x` and `split_at` must both be numeric.")
  }
  split_at <- sort(unique(split_at))
  if (!split_at[1] == 0) {
    split_at <- c(0, split_at)
  }
  if (length(split_at) == 1) {
    # only 0 available
    stop("invalid value for `split_at`.")
  }

  # turn input values to 'split_at' indices
  y <- x
  labs <- split_at
  for (i in 1:length(split_at)) {
    y[x >= split_at[i]] <- i
    # create labels
      # when age group consists of only one age
      labs[i - 1] <- paste0(unique(c(split_at[i - 1], split_at[i] - 1)), collapse = "-")
  }

  # last category
  labs[length(labs)] <- paste0(split_at[length(split_at)], "+")

  factor(labs[y], levels = labs, ordered = TRUE)
}
