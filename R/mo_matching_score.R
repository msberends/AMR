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
#' @param fullname A full taxonomic name, that exists in [`microorganisms$fullname`][microorganisms]
#' @param uncertainty The level of uncertainty set in [as.mo()], see `allow_uncertain` in that function (here, it defaults to 1, but is automatically determined in [as.mo()] based on the number of transformations needed to get to a result)
#' @details The matching score is based on four parameters:
#' 
#' 1. A human pathogenic prevalence \eqn{P}, that is categorised into group 1, 2 and 3 (see [as.mo()]);
#' 2. A kingdom index \eqn{K} is set as follows: Bacteria = 1, Fungi = 2, Protozoa = 3, Archaea = 4, and all others = 5;
#' 3. The level of uncertainty \eqn{U} that is needed to get to a result (1 to 3, see [as.mo()]);
#' 4. The [Levenshtein distance](https://en.wikipedia.org/wiki/Levenshtein_distance) \eqn{L} is the distance between the user input and all taxonomic full names, with the text length of the user input being the maximum distance. A modified version of the Levenshtein distance \eqn{L'} based on the text length of the full name \eqn{F} is calculated as:
#'  
#' \deqn{L' = F - \frac{0.5L}{F}}{L' = (F - 0.5L) / F}
#'   
#' The final matching score \eqn{M} is calculated as:
#' \deqn{M = L' \times \frac{1}{P K U} = \frac{F - 0.5L}{F P K U}}{M = L' * (1 / (P * K * U)) = (F - 0.5L) / (F * P * K * U)}
#' 
#' @export
#' @examples 
#' as.mo("E. coli")
#' mo_uncertainties()
mo_matching_score <- function(x, fullname, uncertainty = 1) {
  # fullname is always a taxonomically valid full name
  levenshtein <- double(length = length(x))
  if (length(fullname) == 1) {
    fullname <- rep(fullname, length(x))
  }
  if (length(x) == 1) {
    x <- rep(x, length(fullname))
  }
  for (i in seq_len(length(x))) {
    # determine Levenshtein distance, but maximise to nchar of fullname
    levenshtein[i] <- min(as.double(utils::adist(x[i], fullname[i], ignore.case = FALSE)),
                          nchar(fullname[i]))
  }
  
  # F = length of fullname
  var_F <- nchar(fullname)
  # L = modified Levenshtein distance
  var_L <- levenshtein
  # P = Prevalence (1 to 3)
  var_P <- MO_lookup[match(fullname, MO_lookup$fullname), "prevalence", drop = TRUE]
  # K = kingdom index (Bacteria = 1, Fungi = 2, Protozoa = 3, Archaea = 4, others = 5)
  var_K <- MO_lookup[match(fullname, MO_lookup$fullname), "kingdom_index", drop = TRUE]
  # U = uncertainty level (1 to 3), as per as.mo()
  var_U <- uncertainty
  
  # matching score:
  (var_F - 0.5 * var_L) / (var_F * var_P * var_K * var_U)
}
