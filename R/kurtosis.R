# ==================================================================== #
# TITLE:                                                               #
# AMR: An R Package for Working with Antimicrobial Resistance Data     #
#                                                                      #
# SOURCE CODE:                                                         #
# https://github.com/msberends/AMR                                     #
#                                                                      #
# PLEASE CITE THIS SOFTWARE AS:                                        #
# Berends MS, Luz CF, Friedrich AW, et al. (2022).                     #
# AMR: An R Package for Working with Antimicrobial Resistance Data.    #
# Journal of Statistical Software, 104(3), 1-31.                       #
# https://doi.org/10.18637/jss.v104.i03                                #
#                                                                      #
# Developed at the University of Groningen and the University Medical  #
# Center Groningen in The Netherlands, in collaboration with many      #
# colleagues from around the world, see our website.                   #
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

#' Kurtosis of the Sample
#'
#' @description Kurtosis is a measure of the "tailedness" of the probability distribution of a real-valued random variable. A normal distribution has a kurtosis of 3 and a excess kurtosis of 0.
#' @param x a vector of values, a [matrix] or a [data.frame]
#' @param na.rm a [logical] to indicate whether `NA` values should be stripped before the computation proceeds
#' @param excess a [logical] to indicate whether the *excess kurtosis* should be returned, defined as the kurtosis minus 3.
#' @seealso [skewness()]
#' @rdname kurtosis
#' @export
#' @examples
#' kurtosis(rnorm(10000))
#' kurtosis(rnorm(10000), excess = TRUE)
kurtosis <- function(x, na.rm = FALSE, excess = FALSE) {
  meet_criteria(na.rm, allow_class = "logical", has_length = 1)
  meet_criteria(excess, allow_class = "logical", has_length = 1)
  UseMethod("kurtosis")
}

#' @method kurtosis default
#' @rdname kurtosis
#' @export
kurtosis.default <- function(x, na.rm = FALSE, excess = FALSE) {
  meet_criteria(na.rm, allow_class = "logical", has_length = 1)
  meet_criteria(excess, allow_class = "logical", has_length = 1)
  x <- as.vector(x)
  if (isTRUE(na.rm)) {
    x <- x[!is.na(x)]
  }
  n <- length(x)
  k <- n * sum((x - mean(x, na.rm = na.rm))^4, na.rm = na.rm) /
    (sum((x - mean(x, na.rm = na.rm))^2, na.rm = na.rm)^2)
  k - ifelse(excess, 3, 0)
}

#' @method kurtosis matrix
#' @rdname kurtosis
#' @export
kurtosis.matrix <- function(x, na.rm = FALSE, excess = FALSE) {
  meet_criteria(na.rm, allow_class = "logical", has_length = 1)
  meet_criteria(excess, allow_class = "logical", has_length = 1)
  apply(x, 2, kurtosis.default, na.rm = na.rm, excess = excess)
}

#' @method kurtosis data.frame
#' @rdname kurtosis
#' @export
kurtosis.data.frame <- function(x, na.rm = FALSE, excess = FALSE) {
  meet_criteria(na.rm, allow_class = "logical", has_length = 1)
  meet_criteria(excess, allow_class = "logical", has_length = 1)
  vapply(FUN.VALUE = double(1), x, kurtosis.default, na.rm = na.rm, excess = excess)
}
