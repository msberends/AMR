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

#' Guess antibiotic column
#'
#' This tries to find a column name in a data set based on information from the \code{\link{antibiotics}} data set.
#' @param tbl a \code{data.frame}
#' @param col a character to look for
#' @param verbose a logical to indicate whether additional info should be printed
#' @importFrom dplyr %>% select filter_all any_vars
#' @export
#' @inheritSection AMR Read more on our website!
# @examples
#
guess_ab <- function(tbl = NULL, col = NULL, verbose = FALSE) {
  if (is.null(tbl) & is.null(col)) {
    return(as.name("guess_ab"))
  }
  #stop("This function should not be called directly.")
  if (length(col) > 1) {
    warning("argument 'col' has length > 1 and only the first element will be used")
    col <- col[1]
  }
  if (!is.data.frame(tbl)) {
    stop("`tbl` must be a data.frame")
  }
  tbl_names <- colnames(tbl)
  ab_result <- antibiotics %>%
    select(atc:trade_name) %>%
    filter_all(any_vars(tolower(.) == tolower(col)))
  if (nrow(ab_result) > 1) {
    # get most likely one
    if (col %in% ab_result$atc) {
      ab_result <- ab_result %>% filter(atc == col)
    } else if (col %in% ab_result$certe) {
      ab_result <- ab_result %>% filter(certe == col)
    } else if (col %in% ab_result$umcg) {
      ab_result <- ab_result %>% filter(umcg == col)
    } else if (col %in% ab_result$umcg) {
      ab_result <- ab_result %>% filter(official == col)
    } else {
      ab_result <- ab_result[1,]
    }
  }
  tbl_result <- tbl_names[tbl_names %in% ab_result]
  if (length(tbl_result) > 1) {
    tbl_result <- tbl_result[1]
    warning('using column `', tbl_result, '` for col "', col, '"', call. = FALSE)
  } else if (length(tbl_result) == 0) {
    if (verbose == TRUE) {
      message('no result found for col "', col, '"')
    }
    return(NULL)
  } else if (verbose == TRUE) {
    message('using column `', tbl_result, '` for col "', col, '"')
  }
  tbl_result
}
