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

#' Guess antibiotic column
#'
#' This tries to find a column name in a data set based on information from the \code{\link{antibiotics}} data set. Also supports WHONET abbreviations. You can look for an antibiotic (trade) name or abbreviation and it will search the \code{data.frame} for any column containing a name or ATC code of that antibiotic.
#' @param tbl a \code{data.frame}
#' @param col a character to look for
#' @param verbose a logical to indicate whether additional info should be printed
#' @importFrom dplyr %>% select filter_all any_vars
#' @importFrom crayon blue
#' @export
#' @inheritSection AMR Read more on our website!
#' @examples
#' df <- data.frame(amox = "S",
#'                  tetr = "R")
#'
#' guess_ab_col(df, "amoxicillin")
#' # [1] "amox"
#' guess_ab_col(df, "J01AA07") # ATC code of tetracycline
#' # [1] "tetr"
#'
#' guess_ab_col(df, "J01AA07", verbose = TRUE)
#' # using column `tetr` for col "J01AA07"
#' # [1] "tetr"
#'
#' # WHONET codes
#' df <- data.frame(AMP_ND10 = "R",
#'                  AMC_ED20 = "S")
#' guess_ab_col(df, "ampicillin")
#' # [1] "AMP_ND10"
#' guess_ab_col(df, "J01CR02")
#' # [1] "AMC_ED20"
#' guess_ab_col(df, as.ab("augmentin"))
#' # [1] "AMC_ED20"
guess_ab_col <- function(tbl = NULL, col = NULL, verbose = FALSE) {
  if (is.null(tbl) & is.null(col)) {
    return(as.name("guess_ab_col"))
  }

  if (length(col) > 1) {
    warning("argument 'col' has length > 1 and only the first element will be used")
    col <- col[1]
  }
  col <- as.character(col)
  if (!is.data.frame(tbl)) {
    stop("`tbl` must be a data.frame")
  }

  if (col %in% colnames(tbl)) {
    ab_result <- col
  } else {
    # sort colnames on length - longest first
    cols <- colnames(tbl[, tbl %>% colnames() %>% nchar() %>% order() %>% rev()])
    df_trans <- data.frame(cols = cols,
                           abs = suppressWarnings(as.ab(cols)),
                           stringsAsFactors = FALSE)
    ab_result <- df_trans[which(df_trans$abs == as.ab(col)), "cols"]
    ab_result <- ab_result[!is.na(ab_result)][1L]
  }

  if (length(ab_result) == 0) {
    if (verbose == TRUE) {
      message('No column found as input for `', col, '`.')
    }
    return(NULL)
  } else {
    if (verbose == TRUE) {
      message(blue(paste0("NOTE: Using column `", bold(ab_result), "` as input for `", col, "`.")))
    }
    return(ab_result)
  }
}
