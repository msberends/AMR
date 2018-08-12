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

#' Count cases with antimicrobial results
#'
#' This counts all cases where antimicrobial interpretations are available. Its use is equal to \code{\link{n_distinct}}.
#' @param ab1,ab2 vector of antibiotic interpretations, they will be transformed internally with \code{\link{as.rsi}} if needed
#' @export
#' @seealso The \code{\link{portion}} functions to calculate resistance and susceptibility.
#' @examples
#' library(dplyr)
#'
#' septic_patients %>%
#'   group_by(hospital_id) %>%
#'   summarise(cipro_p = portion_S(cipr, as_percent = TRUE),
#'             cipro_n = n_rsi(cipr),
#'             genta_p = portion_S(gent, as_percent = TRUE),
#'             genta_n = n_rsi(gent),
#'             combination_p = portion_S(cipr, gent, as_percent = TRUE),
#'             combination_n = n_rsi(cipr, gent))
n_rsi <- function(ab1, ab2 = NULL) {
  if (NCOL(ab1) > 1) {
    stop('`ab1` must be a vector of antimicrobial interpretations', call. = FALSE)
  }
  if (!is.rsi(ab1)) {
    ab1 <- as.rsi(ab1)
  }
  if (!is.null(ab2)) {
    if (NCOL(ab2) > 1) {
      stop('`ab2` must be a vector of antimicrobial interpretations', call. = FALSE)
    }
    if (!is.rsi(ab2)) {
      ab2 <- as.rsi(ab2)
    }
    sum(!is.na(ab1) & !is.na(ab2))
  } else {
    sum(!is.na(ab1))
  }
}
