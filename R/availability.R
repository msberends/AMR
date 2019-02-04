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
# Visit our website for more info: https://msberends.gitab.io/AMR.     #
# ==================================================================== #

#' Check availability of columns
#'
#' Easy check for availability of columns in a data set. This makes it easy to get an idea of which antibiotic combination can be used for calculation with e.g. \code{\link{portion_IR}}.
#' @param tbl a \code{data.frame} or \code{list}
#' @return \code{data.frame} with column names of \code{tbl} as row names and columns: \code{percent_IR}, \code{count}, \code{percent}, \code{visual_availability}.
#' @export
#' @examples
#' availability(septic_patients)
#'
#' library(dplyr)
#' septic_patients %>% availability()
#'
#' septic_patients %>%
#'   select_if(is.rsi) %>%
#'   availability()
#'
#' septic_patients %>%
#'   filter(mo == as.mo("E. coli")) %>%
#'   select_if(is.rsi) %>%
#'   availability()
availability <- function(tbl) {
  x <- base::sapply(tbl, function(x) { 1 - base::sum(base::is.na(x)) / base::length(x) })
  n <- base::sapply(tbl, function(x) base::length(x[!base::is.na(x)]))
  IR <- base::sapply(tbl, function(x) base::ifelse(is.rsi(x), base::round(portion_IR(x, minimum = 0) * 100, 1), "NaN"))
  IR <- paste0(IR, "%")
  IR <- gsub("NaN%", "", IR)
  max_chars <- 50
  x_chars <- strrep("#", round(x, digits = 2) / (1 / max_chars))
  x_chars_empty <- strrep("-", max_chars - nchar(x_chars))
  # x_abnames <- character(length(x))
  # for (i in 1:length(x)) {
  #   if (tbl %>% pull(i) %>% is.rsi()) {
  #     x_abnames[i] <- atc_name(colnames(tbl)[i])
  #   }
  # }
  data.frame(percent_IR = IR,
             count = n,
             percent = paste0(round(x * 100, 1), "%"),
             visual_availabilty = paste0("|", x_chars, x_chars_empty, "|"))
}
