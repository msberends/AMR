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

#' Calculate resistance of isolates
#'
#' This function is deprecated. Use the \code{\link{portion}} functions instead.
#' @inheritParams portion
#' @param ab1,ab2  vector (or column) with antibiotic interpretations. It will be transformed internally with \code{\link{as.rsi}} if needed.
#' @param interpretation antimicrobial interpretation to check for
#' @param ... deprecated parameters to support usage on older versions
#' @importFrom dplyr tibble case_when
#' @export
rsi <- function(ab1,
                ab2 = NULL,
                interpretation = "IR",
                minimum = 30,
                as_percent = FALSE,
                ...) {

  if (all(is.null(ab2))) {
    df <- tibble(ab1 = ab1)
  } else {
    df <- tibble(ab1 = ab1,
                 ab2 = ab2)
  }

  result <- case_when(
    interpretation == "S"             ~ portion_S(df, minimum = minimum, as_percent = FALSE),
    interpretation %in% c("SI", "IS") ~ portion_SI(df, minimum = minimum, as_percent = FALSE),
    interpretation == "I"             ~ portion_I(df, minimum = minimum, as_percent = FALSE),
    interpretation %in% c("RI", "IR") ~ portion_IR(df, minimum = minimum, as_percent = FALSE),
    interpretation == "R"             ~ portion_R(df, minimum = minimum, as_percent = FALSE),
    TRUE ~ -1
  )
  if (result == -1) {
    stop("invalid interpretation")
  }

  .Deprecated(new = paste0("portion_", interpretation))

  if (as_percent == TRUE) {
    percent(result, force_zero = TRUE)
  } else {
    result
  }
}
