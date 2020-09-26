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

#' Calculate the matching score for microorganisms
#' 
#' This helper function is used by [as.mo()] to determine the most probable match of taxonomic records, based on user input. 
#' @param x Any user input value(s)
#' @param n A full taxonomic name, that exists in [`microorganisms$fullname`][microorganisms]
#' @section Matching score for microorganisms:
#' With ambiguous user input in [as.mo()] and all the [`mo_*`][mo_property()] functions, the returned results are chosen based on their matching score using [mo_matching_score()]. This matching score \eqn{m}, ranging from 0 to 100%, is calculated as:
#' 
#' \deqn{m_{(x, n)} = \frac{l_{n} - 0.5 \cdot \min \begin{cases}l_{n} \\ \operatorname{lev}(x, n)\end{cases}}{l_{n} \cdot p_{n} \cdot k_{n}}}{m(x, n) = ( l_n * min(l_n, lev(x, n) ) ) / ( l_n * p_n * k_n )}
#' 
#' where:
#' 
#' * \eqn{x} is the user input;
#' * \eqn{n} is a taxonomic name (genus, species and subspecies) as found in [`microorganisms$fullname`][microorganisms];
#' * \eqn{l_{n}}{l_n} is the length of \eqn{n};
#' * \eqn{\operatorname{lev}}{lev} is the [Levenshtein distance function](https://en.wikipedia.org/wiki/Levenshtein_distance);
#' * \eqn{p_{n}}{p_n} is the human pathogenic prevalence of \eqn{n}, categorised into group \eqn{1}, \eqn{2} and \eqn{3} (see *Details* in `?as.mo`), meaning that \eqn{p = \{1, 2 , 3\}}{p = {1, 2, 3}};
#' * \eqn{k_{n}}{k_n} is the kingdom index of \eqn{n}, set as follows: Bacteria = \eqn{1}, Fungi = \eqn{2}, Protozoa = \eqn{3}, Archaea = \eqn{4}, and all others = \eqn{5}, meaning that \eqn{k = \{1, 2 , 3, 4, 5\}}{k = {1, 2, 3, 4, 5}}.
#' 
#' This means that the user input `x = "E. coli"` gets for *Escherichia coli* a matching score of `r percentage(mo_matching_score("E. coli", "Escherichia coli"), 1)` and for *Entamoeba coli* a matching score of `r percentage(mo_matching_score("E. coli", "Entamoeba coli"), 1)`.
#' 
#' All matches are sorted descending on their matching score and for all user input values, the top match will be returned.
#' @export
#' @examples 
#' as.mo("E. coli")
#' mo_uncertainties()
#' 
#' mo_matching_score("E. coli", "Escherichia coli")
mo_matching_score <- function(x, n) {
  # n is always a taxonomically valid full name
  levenshtein <- double(length = length(x))
  if (length(n) == 1) {
    n <- rep(n, length(x))
  }
  if (length(x) == 1) {
    x <- rep(x, length(n))
  }
  for (i in seq_len(length(x))) {
    # determine Levenshtein distance, but maximise to nchar of n
    levenshtein[i] <- min(as.double(utils::adist(x[i], n[i], ignore.case = FALSE)),
                          nchar(n[i]))
  }
  
  # F = length of fullname
  var_F <- nchar(n)
  # L = modified Levenshtein distance
  var_L <- levenshtein
  # P = Prevalence (1 to 3)
  var_P <- MO_lookup[match(n, MO_lookup$fullname), "prevalence", drop = TRUE]
  # K = kingdom index (Bacteria = 1, Fungi = 2, Protozoa = 3, Archaea = 4, others = 5)
  var_K <- MO_lookup[match(n, MO_lookup$fullname), "kingdom_index", drop = TRUE]
  
  # matching score:
  (var_F - 0.5 * var_L) / (var_F * var_P * var_K)
}
