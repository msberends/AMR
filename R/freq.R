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

#' @importFrom clean freq
#' @export
clean::freq

#' @exportMethod freq.mo
#' @importFrom dplyr n_distinct
#' @importFrom clean freq.default
#' @export
#' @noRd
freq.mo <- function(x, ...) {
  x_noNA <- as.mo(x[!is.na(x)]) # as.mo() to get the newest mo codes
  grams <- mo_gramstain(x_noNA, language = NULL)
  freq.default(x = x, ...,
               .add_header = list(`Gram-negative` = paste0(format(sum(grams == "Gram-negative", na.rm = TRUE),
                                                                  big.mark = ",",
                                                                  decimal.mark = "."),
                                                           " (", percent(sum(grams == "Gram-negative", na.rm = TRUE) / length(grams), force_zero = TRUE, round = 2),
                                                           ")"),
                                  `Gram-positive` = paste0(format(sum(grams == "Gram-positive", na.rm = TRUE),
                                                                  big.mark = ",",
                                                                  decimal.mark = "."),
                                                           " (", percent(sum(grams == "Gram-positive", na.rm = TRUE) / length(grams), force_zero = TRUE, round = 2),
                                                           ")"),
                                  `Unique genera` = n_distinct(mo_genus(x_noNA, language = NULL)),
                                  `Unique species` = n_distinct(paste(mo_genus(x_noNA, language = NULL),
                                                             mo_species(x_noNA, language = NULL)))))
}

#' @exportMethod freq.rsi
#' @importFrom clean freq.default
#' @export
#' @noRd
freq.rsi <- function(x, ...) {
  x_name <- deparse(substitute(x))
  x_name <- gsub(".*[$]", "", x_name)
  ab <- suppressMessages(suppressWarnings(AMR::as.ab(x_name)))
  if (!is.na(ab)) {
    freq.default(x = x, ...,
                 .add_header = list(Drug = paste0(ab_name(ab), " (", ab, ", ", ab_atc(ab), ")"),
                                    group = ab_group(ab),
                                    `%SI` = AMR::portion_SI(x, minimum = 0, as_percent = TRUE)))
  } else {
    freq.default(x = x, ...,
                 .add_header = list(`%SI` = AMR::portion_SI(x, minimum = 0, as_percent = TRUE)))
  }
}
