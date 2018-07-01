# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Analysis                              #
#                                                                      #
# AUTHORS                                                              #
# Berends MS (m.s.berends@umcg.nl), Luz CF (c.f.luz@umcg.nl)           #
#                                                                      #
# LICENCE                                                              #
# This program is free software; you can redistribute it and/or modify #
# it under the terms of the GNU General Public License version 2.0,    #
# as published by the Free Software Foundation.                        #
#                                                                      #
# This program is distributed in the hope that it will be useful,      #
# but WITHOUT ANY WARRANTY; without even the implied warranty of       #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        #
# GNU General Public License for more details.                         #
# ==================================================================== #

#' Symbol of a p value
#'
#' Return the symbol related to the p value: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#' @param p p value
#' @param emptychar text to show when \code{p > 0.1}
#' @return Text
#' @export
p.symbol <- function(p, emptychar = " ") {
  instelling.oud <- options()$scipen
  options(scipen = 999)
  s <- ''
  s[1:length(p)] <- ''
  for (i in 1:length(p)) {
    if (is.na(p[i])) {
      s[i] <- NA
      next
    }
    if (p[i] > 1) {
      s[i] <- NA
      next
    } else {
      p_test <- p[i]
    }

    if (p_test > 0.1) {
      s[i] <- emptychar
    } else if (p_test > 0.05) {
      s[i] <- '.'
    } else if (p_test > 0.01) {
      s[i] <- '*'
    } else if (p_test > 0.001) {
      s[i] <- '**'
    } else if (p_test >= 0) {
      s[i] <- '***'
    }
  }
  options(scipen = instelling.oud)
  s
}
