# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Analysis                              #
#                                                                      #
# SOURCE                                                               #
# https://gitlab.com/msberends/AMR                                     #
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
# Visit our website for more info: https://msberends.gitlab.io/AMR.    #
# ==================================================================== #

if ("cleaner" %in% rownames(utils::installed.packages())) {
  freq <- get("freq", envir = asNamespace("cleaner"))
  freq.default <- get("freq.default", envir = asNamespace("cleaner"))
} else {
  freq <- ""
  freq.default <- ""
}

#' @method freq mo
#' @export
#' @noRd
freq.mo <- function(x, ...) {
  x_noNA <- as.mo(x[!is.na(x)]) # as.mo() to get the newest mo codes
  grams <- mo_gramstain(x_noNA, language = NULL)
  digits <- list(...)$digits
  if (is.null(digits)) {
    digits <- 2
  }
  freq.default(x = x, ...,
               .add_header = list(`Gram-negative` = paste0(format(sum(grams == "Gram-negative", na.rm = TRUE),
                                                                  big.mark = ",",
                                                                  decimal.mark = "."),
                                                           " (", percentage(sum(grams == "Gram-negative", na.rm = TRUE) / length(grams), digits = digits),
                                                           ")"),
                                  `Gram-positive` = paste0(format(sum(grams == "Gram-positive", na.rm = TRUE),
                                                                  big.mark = ",",
                                                                  decimal.mark = "."),
                                                           " (", percentage(sum(grams == "Gram-positive", na.rm = TRUE) / length(grams), digits = digits),
                                                           ")"),
                                  `No of genera` = n_distinct(mo_genus(x_noNA, language = NULL)),
                                  `No of species` = n_distinct(paste(mo_genus(x_noNA, language = NULL),
                                                                     mo_species(x_noNA, language = NULL)))))
}

#' @method freq rsi
#' @export
#' @noRd
freq.rsi <- function(x, ...) {
  x_name <- deparse(substitute(x))
  x_name <- gsub(".*[$]", "", x_name)
  ab <- suppressMessages(suppressWarnings(as.ab(x_name)))
  if (!is.na(ab)) {
    freq.default(x = x, ...,
                 .add_header = list(Drug = paste0(ab_name(ab), " (", ab, ", ", ab_atc(ab), ")"),
                                    group = ab_group(ab),
                                    `%SI` = susceptibility(x, minimum = 0, as_percent = TRUE)))
  } else {
    freq.default(x = x, ...,
                 .add_header = list(`%SI` = susceptibility(x, minimum = 0, as_percent = TRUE)))
  }
}
