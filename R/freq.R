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

freq_def <- clean:::freq.default
#' @exportMethod freq.mo
#' @importFrom dplyr n_distinct
#' @importFrom clean freq
#' @export
#' @noRd
freq.mo <- function(x, ...) {
  # replace with freq.default() if next `clean` version is published on CRAN
  freq_def(x = x, ...,
           .add_header = list(families = n_distinct(mo_family(x, language = NULL)),
                              genera = n_distinct(mo_genus(x, language = NULL)),
                              species = n_distinct(paste(mo_genus(x, language = NULL),
                                                         mo_species(x, language = NULL)))))
}

#' @exportMethod freq.rsi
#' @importFrom clean freq
#' @export
#' @noRd
freq.rsi <- function(x, ...) {
  x_name <- deparse(substitute(x))
  x_name <- gsub(".*[$]", "", x_name)
  ab <- suppressMessages(suppressWarnings(AMR::as.ab(x_name)))
  if (!is.na(ab)) {
    freq_def(x = x, ...,
             .add_header = list(Drug = paste0(ab_name(ab), " (", ab, ", ", ab_atc(ab), ")"),
                                group = ab_group(ab),
                                `%SI` = portion_SI(x, minimum = 0, as_percent = TRUE)))
  } else {
    freq_def(x = x, ...,
             .add_header = list(`%SI` = portion_SI(x, minimum = 0, as_percent = TRUE)))
  }
}
