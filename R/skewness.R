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

#' Skewness of the Sample
#'
#' @description Skewness is a measure of the asymmetry of the probability distribution of a real-valued random variable about its mean.
#'
#' When negative ('left-skewed'): the left tail is longer; the mass of the distribution is concentrated on the right of a histogram. When positive ('right-skewed'): the right tail is longer; the mass of the distribution is concentrated on the left of a histogram. A normal distribution has a skewness of 0.
#' @param x A vector of values, a [matrix] or a [data.frame]
#' @param na.rm A [logical] value indicating whether `NA` values should be stripped before the computation proceeds
#' @seealso [kurtosis()]
#' @rdname skewness
#' @export
#' @examples
#' skewness(runif(1000))
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
  if (isTRUE(na.rm)) {
    x <- x[!is.na(x)]
  }
  n <- length(x)
  (sum((x - mean(x))^3) / n) / (sum((x - mean(x))^2) / n)^(3 / 2)
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
