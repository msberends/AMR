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
#' This tries to find a column name in a data set based on information from the \code{\link{antibiotics}} data set. You can look for an antibiotic (trade) of abbreviation and it will search the data for any column containing a name or ATC code of that antibiotic.
#' @param tbl a \code{data.frame}
#' @param col a character to look for
#' @param verbose a logical to indicate whether additional info should be printed
#' @importFrom dplyr %>% select filter_all any_vars
#' @export
#' @inheritSection AMR Read more on our website!
#' @examples
#' df <- data.frame(amox = "S",
#'                  tetr = "R")
#'
#' guess_ab_col(df, "amoxicillin")
#' # [1] "amox"
#' guess_ab_col(df, "J01AA07") # ATC code of Tetracycline
#' # [1] "tetr"
#'
#' guess_ab_col(df, "J01AA07", verbose = TRUE)
#' # using column `tetr` for col "J01AA07"
#' # [1] "tetr"
guess_ab_col <- function(tbl = NULL, col = NULL, verbose = FALSE) {
  if (is.null(tbl) & is.null(col)) {
    return(as.name("guess_ab_col"))
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
  if (col %in% tbl_names) {
    return(col)
  }
  ab_result <- antibiotics %>%
    select(atc:trade_name) %>%
    filter_all(any_vars(tolower(.) == tolower(col))) %>%
    filter_all(any_vars(. %in% tbl_names))

  if (nrow(ab_result) == 0 & nchar(col) > 4) {
    # use like when col >= 5 characters
    ab_result <- antibiotics %>%
      select(atc:trade_name) %>%
      filter_all(any_vars(tolower(.) %like% tolower(col))) %>%
      filter_all(any_vars(. %in% tbl_names))
  }

  if (nrow(ab_result) > 1) {
    # looking more and more for reliable hit
    ab_result_1 <- ab_result %>% filter(tolower(atc) == tolower(col))
    if (nrow(ab_result_1) == 0) {
      ab_result_1 <- ab_result %>% filter(tolower(certe) == tolower(col))
    }
    if (nrow(ab_result_1) == 0) {
      ab_result_1 <- ab_result %>% filter(tolower(umcg) == tolower(col))
    }
    if (nrow(ab_result_1) == 0) {
      ab_result_1 <- ab_result %>% filter(tolower(official) == tolower(col))
    }
    if (nrow(ab_result_1) == 0) {
      ab_result_1 <- ab_result[1, ]
    }
    ab_result <- ab_result_1
  }

  if (length(ab_result) == 0) {
    if (verbose == TRUE) {
      message('no result found for col "', col, '"')
    }
    return(NULL)
  } else {
    result <- tbl_names[tbl_names %in% ab_result]
    if (length(result) == 0) {
      if (verbose == TRUE) {
        message('no result found for col "', col, '"')
      }
      return(NULL)
    }
    if (verbose == TRUE) {
      message('using column `', result, '` for col "', col, '"')
    }
    return(result)
  }
}
