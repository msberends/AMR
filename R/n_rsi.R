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
#' This counts all cases where antimicrobial interpretations are available. The way it can be used is equal to \code{\link{n_distinct}}. Its function is equal to \code{count_S(...) + count_IR(...)}.
#' @inheritParams portion
#' @export
#' @seealso \code{\link[AMR]{count}_*} to count resistant and susceptibile isolates per interpretation type.\cr
#' \code{\link{portion}_*} to calculate microbial resistance and susceptibility.
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
n_rsi <- function(...) {
  # only print warnings once, if needed
  count_S(...) + suppressWarnings(count_IR(...))
}
