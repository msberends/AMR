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

#' Check availability of columns
#'
#' Easy check for availability of columns in a data set. This makes it easy to get an idea of which antibiotic combination can be used for calculation with e.g. \code{\link{portion_IR}}.
#' @param tbl a \code{data.frame} or \code{list}
#' @param width number of characters to present the visual availability, defaults to filling the width of the console
#' @return \code{data.frame} with column names of \code{tbl} as row names
#' @inheritSection AMR Read more on our website!
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
availability <- function(tbl, width = NULL) {
  x <- base::sapply(tbl, function(x) { 1 - base::sum(base::is.na(x)) / base::length(x) })
  n <- base::sapply(tbl, function(x) base::length(x[!base::is.na(x)]))
  IR <- base::sapply(tbl, function(x) base::ifelse(is.rsi(x), portion_IR(x, minimum = 0), NA))
  IR_print <- character(length(IR))
  IR_print[!is.na(IR)] <- percent(IR[!is.na(IR)], round = 1, force_zero = TRUE)
  IR_print[is.na(IR)] <- ""

  if (is.null(width)) {
    width <- options()$width -
      (max(nchar(colnames(tbl))) +
         # count col
         8 +
         # available % column
         10 +
         # resistant % column
         10 +
         # extra margin
         5)
    width <- width / 2
  }

  if (length(IR[is.na(IR)]) == ncol(tbl)) {
    width <- width * 2 + 10
  }

  x_chars_IR <- strrep("#", round(width * IR, digits = 2))
  x_chars_S <- strrep("-", width - nchar(x_chars_IR))
  vis_resistance <- paste0("|", x_chars_IR, x_chars_S, "|")
  vis_resistance[is.na(IR)] <- ""

  x_chars <- strrep("#", round(x, digits = 2) / (1 / width))
  x_chars_empty <- strrep("-", width - nchar(x_chars))

  df <- data.frame(count = n,
                   available = percent(x, round = 1, force_zero = TRUE),
                   visual_availabilty = paste0("|", x_chars, x_chars_empty, "|"),
                   resistant = IR_print,
                   visual_resistance = vis_resistance)
  if (length(IR[is.na(IR)]) == ncol(tbl)) {
    df[,1:3]
  } else {
    df
  }
}
