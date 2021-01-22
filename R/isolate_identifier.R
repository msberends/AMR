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
#' @inheritSection lifecycle Experimental Lifecycle
#' @inheritParams eucast_rules
#' @param cols_ab a character vector of column names of `x`, or (a combination with) an [antibiotic selector function]([ab_class()]), such as [carbapenems()] and [aminoglycosides()]
#' @rdname isolate_identifier
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
    if (is.null(col_mo)) {
      # no column found, then ignore the argument
      col_mo <- FALSE
    }
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
  
  # cope with empty values
  if (length(cols_ab) == 0 && all(x[, col_mo, drop = TRUE] == "", na.rm = TRUE)) {
    warning_("in isolate_identifier(): no column with microorganisms and no columns with antimicrobial agents found", call = FALSE)
  } else if (length(cols_ab) == 0) {
    warning_("in isolate_identifier(): no columns with antimicrobial agents found", call = FALSE)
  }
  
  out <- x[, c(col_mo, cols_ab), drop = FALSE]
  out <- do.call(paste, c(out, sep = ""))
  out <- gsub("NA", ".", out, fixed = TRUE)
  out <- set_clean_class(out, new_class = c("isolate_identifier", "character"))
  attr(out, "ab") <- cols_ab
  out
}

#' @method all.equal isolate_identifier
#' @inheritParams base::all.equal
#' @param ignore_empty_results a logical to indicate whether empty results must be ignored, so that only values R, S and I will be compared
#' @rdname isolate_identifier
#' @export
all.equal.isolate_identifier <- function(target, current, ignore_empty_results = TRUE, ...) {
  meet_criteria(target, allow_class = "isolate_identifier")
  meet_criteria(current, allow_class = "isolate_identifier")
  meet_criteria(ignore_empty_results, allow_class = "logical", has_length = 1)
  
  if (isTRUE(all.equal.character(target, current))) {
    return(TRUE)
  }
  # vectorise over both target and current
  if (length(target) > 1 && length(current) == 1) {
    current <- rep(current, length(target))
  } else if (length(current) > 1 && length(target) == 1) {
    target <- rep(target, length(current))
  }
  stop_if(length(target) != length(current),
          "length of `target` and `current` must be the same, or one must be 1")
  
  get_vector <- function(x) {
    if (grepl("|", x, fixed = TRUE)) {
      mo <- gsub("(.*)\\|.*", "\\1", x)
    } else {
      mo <- NULL
    }
    if (grepl("|", x, fixed = TRUE)) {
      ab <- gsub(".*\\|(.*)", "\\1", x)
    } else {
      ab <- x
    }
    ab <- strsplit(ab, "")[[1L]]
    if (is.null(mo)) {
      out <- as.character(ab)
      names(out) <- attributes(x)$ab
    } else {
      out <- as.character(c(mo, ab))
      names(out) <- c("mo", attributes(x)$ab)
    }
    out
  }
  
  # run it
  for (i in seq_len(length(target))) {
    if (i == 1) {
      df <- data.frame(object = paste0(c("target[", "current["), i, "]"))
    }
    trgt <- get_vector(target[i])
    crnt <- get_vector(current[i])
    if (ignore_empty_results == TRUE) {
      diff <- names(trgt[trgt != crnt & trgt != "." & crnt != "."])  
    } else {
      diff <- names(trgt[trgt != crnt])
    }
    
  }
  
  stop("THIS FUNCTION IS WORK IN PROGRESS AND NOT AVAILABLE IN THIS BETA VERSION")
  
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
