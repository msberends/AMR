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

#' Kurtosis of the sample
#'
#' @description Kurtosis is a measure of the "tailedness" of the probability distribution of a real-valued random variable.
#'
#' @param x a vector of values, a \code{matrix} or a \code{data frame}
#' @param na.rm a logical value indicating whether \code{NA} values should be stripped before the computation proceeds.
#' @exportMethod kurtosis
#' @seealso \code{\link{skewness}}
#' @rdname kurtosis
#' @inheritSection AMR Read more on our website!
#' @export
kurtosis <- function(x, na.rm = FALSE) {
  UseMethod("kurtosis")
}

#' @exportMethod kurtosis.default
#' @rdname kurtosis
#' @export
kurtosis.default <- function (x, na.rm = FALSE) {
  x <- as.vector(x)
  if (na.rm == TRUE) {
    x <- x[!is.na(x)]
  }
  n <- length(x)
  n * base::sum((x - base::mean(x, na.rm = na.rm))^4, na.rm = na.rm) /
    (base::sum((x - base::mean(x, na.rm = na.rm))^2, na.rm = na.rm)^2)
}

#' @exportMethod kurtosis.matrix
#' @rdname kurtosis
#' @export
kurtosis.matrix <- function (x, na.rm = FALSE) {
  base::apply(x, 2, kurtosis.default, na.rm = na.rm)
}

#' @exportMethod kurtosis.data.frame
#' @rdname kurtosis
#' @export
kurtosis.data.frame <- function (x, na.rm = FALSE) {
  base::sapply(x, kurtosis.default, na.rm = na.rm)
}
