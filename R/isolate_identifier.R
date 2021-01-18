# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Analysis for R                        #
#                                                                      #
# SOURCE                                                               #
# https://github.com/msberends/AMR                                     #
#                                                                      #
# LICENCE                                                              #
# (c) 2018-2021 Berends MS, Luz CF et al.                              #
# Developed at the University of Groningen, the Netherlands, in        #
# collaboration with non-profit organisations Certe Medical            #
# Diagnostics & Advice, and University Medical Center Groningen.       # 
#                                                                      #
# This R package is free software; you can freely use and distribute   #
# it for both personal and commercial purposes under the terms of the  #
# GNU General Public License version 2.0 (GNU GPL-2), as published by  #
# the Free Software Foundation.                                        #
# We created this package for both routine data analysis and academic  #
# research and it was publicly released in the hope that it will be    #
# useful, but it comes WITHOUT ANY WARRANTY OR LIABILITY.              #
#                                                                      #
# Visit our website for the full manual and a complete tutorial about  #
# how to conduct AMR analysis: https://msberends.github.io/AMR/        #
# ==================================================================== #

#' Create Identifier of an Isolate
#' 
#' This function will paste the microorganism code with all antimicrobial results into one string for each row in a data set. This is useful to compare isolates, e.g. between institutions or regions, when there is no genotyping available.
#' @inheritSection lifecycle Maturing Lifecycle
#' @inheritParams eucast_rules
#' @param cols_ab a character vector of column names of `x`, or (a combination with) an [antibiotic selector function]([ab_class()]), such as [carbapenems()] and [aminoglycosides()]
#' @export
#' @inheritSection AMR Read more on Our Website!
#' @examples 
#' # automatic selection of microorganism and antibiotics (i.e., all <rsi> columns, see ?as.rsi)
#' x <- isolate_identifier(example_isolates)
#' 
#' # ignore microorganism codes, only use antimicrobial results
#' x <- isolate_identifier(example_isolates, col_mo = FALSE, cols_ab = c("AMX", "TZP", "GEN", "TOB"))
#' 
#' # select antibiotics from certain antibiotic classes
#' x <- isolate_identifier(example_isolates, cols_ab = c(carbapenems(), aminoglycosides()))
isolate_identifier <- function(x, col_mo = NULL, cols_ab = NULL) {
  if (is.null(col_mo)) {
    col_mo <- search_type_in_df(x, "mo")
  }
  if (isFALSE(col_mo)) {
    # is FALSE then ignore mo column
    x$col_mo <- ""
    col_mo <- "col_mo"
  } else if (!is.null(col_mo)) {
    x[, col_mo] <- paste0(as.mo(x[, col_mo, drop = TRUE]), "|")
  }
  
  cols_ab <- deparse(substitute(cols_ab)) # support ab class selectors: isolate_identifier(x, cols_ab = carbapenems())
  if (identical(cols_ab, "NULL")) {
    cols_ab <- colnames(x)[vapply(FUN.VALUE = logical(1), x, is.rsi)]
  } else {
    cols_ab <- tryCatch(colnames(x[, eval(parse(text = cols_ab), envir = parent.frame())]),
                        # tryCatch adds 4 calls, so total is -5
                        error = function(e) stop_(e$message, call = -5))
  }
  if (length(cols_ab) == 0) {
    warning_("no columns with antimicrobial agents found", call = TRUE)
  }
  
  out <- x[, c(col_mo, cols_ab), drop = FALSE]
  out <- do.call(paste, c(out, sep = ""))
  out <- gsub("NA", ".", out, fixed = TRUE)
  set_clean_class(out, new_class = c("isolate_identifier", "character"))
}

#' @method print isolate_identifier
#' @export
#' @noRd
print.isolate_identifier <- function(x, ...) {
  print(as.character(x), ...)
}

#' @method [ isolate_identifier
#' @export
#' @noRd
"[.isolate_identifier" <- function(x, ...) {
  y <- NextMethod()
  attributes(y) <- attributes(x)
  y
}
#' @method [[ isolate_identifier
#' @export
#' @noRd
"[[.isolate_identifier" <- function(x, ...) {
  y <- NextMethod()
  attributes(y) <- attributes(x)
  y
}
#' @method [<- isolate_identifier
#' @export
#' @noRd
"[<-.isolate_identifier" <- function(i, j, ..., value) {
  y <- NextMethod()
  attributes(y) <- attributes(i)
  y
}
#' @method [[<- isolate_identifier
#' @export
#' @noRd
"[[<-.isolate_identifier" <- function(i, j, ..., value) {
  y <- NextMethod()
  attributes(y) <- attributes(i)
  y
}
#' @method c isolate_identifier
#' @export
#' @noRd
c.isolate_identifier <- function(x, ...) {
  y <- NextMethod()
  attributes(y) <- attributes(x)
  y
}

#' @method unique isolate_identifier
#' @export
#' @noRd
unique.isolate_identifier <- function(x, incomparables = FALSE, ...) {
  y <- NextMethod()
  attributes(y) <- attributes(x)
  y
}
