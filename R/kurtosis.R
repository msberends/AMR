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

#' Kurtosis of the sample
#'
#' @description Kurtosis is a measure of the "tailedness" of the probability distribution of a real-valued random variable.
#' @inheritSection lifecycle Questioning lifecycle
#' @param x a vector of values, a [`matrix`] or a [`data.frame`]
#' @param na.rm a logical value indicating whether `NA` values should be stripped before the computation proceeds.
#' @seealso [skewness()]
#' @rdname kurtosis
#' @inheritSection AMR Read more on our website!
#' @export
kurtosis <- function(x, na.rm = FALSE) {
  UseMethod("kurtosis")
}

#' @method kurtosis default
#' @rdname kurtosis
#' @export
kurtosis.default <- function(x, na.rm = FALSE) {
  x <- as.vector(x)
  if (na.rm == TRUE) {
    x <- x[!is.na(x)]
  }
  n <- length(x)
  n * base::sum((x - base::mean(x, na.rm = na.rm))^4, na.rm = na.rm) /
    (base::sum((x - base::mean(x, na.rm = na.rm))^2, na.rm = na.rm)^2)
}

#' @method kurtosis matrix
#' @rdname kurtosis
#' @export
kurtosis.matrix <- function(x, na.rm = FALSE) {
  base::apply(x, 2, kurtosis.default, na.rm = na.rm)
}

#' @method kurtosis data.frame
#' @rdname kurtosis
#' @export
kurtosis.data.frame <- function(x, na.rm = FALSE) {
  base::sapply(x, kurtosis.default, na.rm = na.rm)
}
