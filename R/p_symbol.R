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

#' Symbol of a p-value
#'
#' Return the symbol related to the p-value: 0 '`***`' 0.001 '`**`' 0.01 '`*`' 0.05 '`.`' 0.1 ' ' 1. Values above `p = 1` will return `NA`.
#' @inheritSection lifecycle Questioning lifecycle
#' @param p p value
#' @param emptychar text to show when `p > 0.1`
#' @return Text
#' @inheritSection AMR Read more on our website!
#' @export
p_symbol <- function(p, emptychar = " ") {

  p <- as.double(p)
  s <- rep(NA_character_, length(p))

  s[p <= 1] <- emptychar
  s[p <= 0.100] <- "."
  s[p <= 0.050] <- "*"
  s[p <= 0.010] <- "**"
  s[p <= 0.001] <- "***"
  
  s
}
