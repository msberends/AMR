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

#' Pattern Matching
#'
#' Convenient wrapper around [grep()] to match a pattern: `x %like% pattern`. It always returns a [`logical`] vector and is always case-insensitive (use `x %like_case% pattern` for case-sensitive matching). Also, `pattern` can be as long as `x` to compare items of each index in both vectors, or they both can have the same length to iterate over all cases.
#' @inheritSection lifecycle Stable lifecycle
#' @param x a character vector where matches are sought, or an object which can be coerced by [as.character()] to a character vector.
#' @param pattern a character string containing a regular expression (or [`character`] string for `fixed = TRUE`) to be matched in the given character vector. Coerced by [as.character()] to a character string if possible.  If a [`character`] vector of length 2 or more is supplied, the first element is used with a warning.
#' @param ignore.case if `FALSE`, the pattern matching is *case sensitive* and if `TRUE`, case is ignored during matching.
#' @return A [`logical`] vector
#' @name like
#' @rdname like
#' @export
#' @details
#' The `%like%` function:
#' * Is case insensitive (use `%like_case%` for case-sensitive matching)
#' * Supports multiple patterns
#' * Checks if `pattern` is a regular expression and sets `fixed = TRUE` if not, to greatly improve speed
#' * Tries again with `perl = TRUE` if regex fails
#' 
#' Using RStudio? This function can also be inserted from the Addins menu and can have its own Keyboard Shortcut like `Ctrl+Shift+L` or `Cmd+Shift+L` (see `Tools` > `Modify Keyboard Shortcuts...`).
#' @source Idea from the [`like` function from the `data.table` package](https://github.com/Rdatatable/data.table/blob/master/R/like.R)
#' @seealso [base::grep()]
#' @inheritSection AMR Read more on our website!
#' @examples
#' # simple test
#' a <- "This is a test"
#' b <- "TEST"
#' a %like% b
#' #> TRUE
#' b %like% a
#' #> FALSE
#'
#' # also supports multiple patterns, length must be equal to x
#' a <- c("Test case", "Something different", "Yet another thing")
#' b <- c(     "case",           "diff",      "yet")
#' a %like% b
#' #> TRUE TRUE TRUE
#'
#' # get isolates whose name start with 'Ent' or 'ent'
#' \dontrun{
#' library(dplyr)
#' example_isolates %>%
#'   filter(mo_name(mo) %like% "^ent") %>% 
#'   freq(mo)
#' }
like <- function(x, pattern, ignore.case = TRUE) {
  # set to fixed if no regex found
  fixed <- all(!grepl("[$.^*?+}{|)(]", pattern))
  if (ignore.case == TRUE) {
    # set here, otherwise if fixed = TRUE, this warning will be thrown: argument 'ignore.case = TRUE' will be ignored
    x <- tolower(x)
    pattern <- tolower(pattern)
  }
  
  if (length(pattern) > 1) {
    if (length(x) != length(pattern)) {
      if (length(x) == 1) {
        x <- rep(x, length(pattern))
      }
      # return TRUE for every 'x' that matches any 'pattern', FALSE otherwise
      res <- sapply(pattern, function(pttrn) base::grepl(pttrn, x, ignore.case = FALSE, fixed = fixed))
      res2 <- as.logical(rowSums(res))
      # get only first item of every hit in pattern
      res2[duplicated(res)] <- FALSE
      res2[rowSums(res) == 0] <- NA
      return(res2)
    } else {
      # x and pattern are of same length, so items with each other
      res <- vector(length = length(pattern))
      for (i in seq_len(length(res))) {
        if (is.factor(x[i])) {
          res[i] <- as.integer(x[i]) %in% base::grep(pattern[i], levels(x[i]), ignore.case = FALSE, fixed = fixed)
        } else {
          res[i] <- base::grepl(pattern[i], x[i], ignore.case = FALSE, fixed = fixed)
        }
      }
      return(res)
    }
  }

  # the regular way how grepl works; just one pattern against one or more x
  if (is.factor(x)) {
    as.integer(x) %in% base::grep(pattern, levels(x), ignore.case = FALSE, fixed = fixed)
  } else {
    tryCatch(base::grepl(pattern, x, ignore.case = FALSE, fixed = fixed),
             error = function(e) ifelse(grepl("Invalid regexp", e$message),
                                        # try with perl = TRUE:
                                        return(base::grepl(pattern = pattern, x = x,
                                                           ignore.case = FALSE, 
                                                           fixed = fixed,
                                                           perl = TRUE)),
                                        # stop otherwise
                                        stop(e$message)))
  }
}

#' @rdname like
#' @export
"%like%" <- function(x, pattern) {
  like(x, pattern, ignore.case = TRUE)
}

#' @rdname like
#' @export
"%like_case%" <- function(x, pattern) {
  like(x, pattern, ignore.case = FALSE)
}
