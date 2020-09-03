# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Analysis                              #
#                                                                      #
# SOURCE                                                               #
# https://github.com/msberends/AMR                                     #
#                                                                      #
# LICENCE                                                              #
# (c) 2018-2020 Berends MS, Luz CF et al.                              #
#                                                                      #
# This R package is free software; you can freely use and distribute   #
# it for both personal and commercial purposes under the terms of the  #
# GNU General Public License version 2.0 (GNU GPL-2), as published by  #
# the Free Software Foundation.                                        #
#                                                                      #
# We created this package for both routine data analysis and academic  #
# research and it was publicly released in the hope that it will be    #
# useful, but it comes WITHOUT ANY WARRANTY OR LIABILITY.              #
# Visit our website for more info: https://msberends.github.io/AMR.    #
# ==================================================================== #

#' Skewness of the sample
#'
#' @description Skewness is a measure of the asymmetry of the probability distribution of a real-valued random variable about its mean.
#'
#' When negative: the left tail is longer; the mass of the distribution is concentrated on the right of the figure. When positive: the right tail is longer; the mass of the distribution is concentrated on the left of the figure.
#' @inheritSection lifecycle Questioning lifecycle
#' @param x a vector of values, a [`matrix`] or a [`data.frame`]
#' @param na.rm a logical value indicating whether `NA` values should be stripped before the computation proceeds.
#' @seealso [kurtosis()]
#' @rdname skewness
#' @inheritSection AMR Read more on our website!
#' @export
skewness <- function(x, na.rm = FALSE) {
  UseMethod("skewness")
}

#' @method skewness default
#' @rdname skewness
#' @export
skewness.default <- function(x, na.rm = FALSE) {
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
  apply(x, 2, skewness.default, na.rm = na.rm)
}

#' @method skewness data.frame
#' @rdname skewness
#' @export
skewness.data.frame <- function(x, na.rm = FALSE) {
  sapply(x, skewness.default, na.rm = na.rm)
}
