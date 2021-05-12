# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Data Analysis for R                   #
#                                                                      #
# SOURCE                                                               #
# https://github.com/msberends/AMR                                     #
#                                                                      #
# LICENCE                                                              #
# (c) 2018-2021 Berends MS, Luz CF et al.                              #
# Developed at the University of Groningen, the Netherlands, in        #
# collaboration with non-profit organisations Certe Medical            #
# Diagnostics & Advice, and University Medical Center Groningen.       # 
#                                                                      #
# This R package is free software; you can freely use and distribute   #
# it for both personal and commercial purposes under the terms of the  #
# GNU General Public License version 2.0 (GNU GPL-2), as published by  #
# the Free Software Foundation.                                        #
# We created this package for both routine data analysis and academic  #
# research and it was publicly released in the hope that it will be    #
# useful, but it comes WITHOUT ANY WARRANTY OR LIABILITY.              #
#                                                                      #
# Visit our website for the full manual and a complete tutorial about  #
# how to conduct AMR data analysis: https://msberends.github.io/AMR/   #
# ==================================================================== #

#' Skewness of the Sample
#'
#' @description Skewness is a measure of the asymmetry of the probability distribution of a real-valued random variable about its mean.
#'
#' When negative ('left-skewed'): the left tail is longer; the mass of the distribution is concentrated on the right of a histogram. When positive ('right-skewed'): the right tail is longer; the mass of the distribution is concentrated on the left of a histogram. A normal distribution has a skewness of 0.
#' @inheritSection lifecycle Stable Lifecycle
#' @param x a vector of values, a [matrix] or a [data.frame]
#' @param na.rm a [logical] value indicating whether `NA` values should be stripped before the computation proceeds
#' @seealso [kurtosis()]
#' @rdname skewness
#' @inheritSection AMR Read more on Our Website!
#' @export
skewness <- function(x, na.rm = FALSE) {
  meet_criteria(na.rm, allow_class = "logical", has_length = 1)
  UseMethod("skewness")
}

#' @method skewness default
#' @rdname skewness
#' @export
skewness.default <- function(x, na.rm = FALSE) {
  meet_criteria(na.rm, allow_class = "logical", has_length = 1)
  x <- as.vector(x)
  if (na.rm == TRUE) {
    x <- x[!is.na(x)]
  }
  n <- length(x)
  (sum((x - mean(x))^3) / n) / (sum((x - mean(x)) ^ 2) / n) ^ (3 / 2)
}

#' @method skewness matrix
#' @rdname skewness
#' @export
skewness.matrix <- function(x, na.rm = FALSE) {
  meet_criteria(na.rm, allow_class = "logical", has_length = 1)
  apply(x, 2, skewness.default, na.rm = na.rm)
}

#' @method skewness data.frame
#' @rdname skewness
#' @export
skewness.data.frame <- function(x, na.rm = FALSE) {
  meet_criteria(na.rm, allow_class = "logical", has_length = 1)
  vapply(FUN.VALUE = double(1), x, skewness.default, na.rm = na.rm)
}
