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
# Visit our website for more info: https://msberends.gitlab.io/AMR.    #
# ==================================================================== #

#' Symbol of a p-value
#'
#' Return the symbol related to the p-value: 0 '`***`' 0.001 '`**`' 0.01 '`*`' 0.05 '`.`' 0.1 ' ' 1. Values above `p = 1` will return `NA`.
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
