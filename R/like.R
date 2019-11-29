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

#' Pattern Matching
#'
#' Convenient wrapper around [base::grep()] to match a pattern: `a %like% b`. It always returns a [`logical`] vector and is always case-insensitive (use `a %like_case% b` for case-sensitive matching). Also, `pattern` (*b*) can be as long as `x` (*a*) to compare items of each index in both vectors, or can both have the same length to iterate over all cases.
#' @param x a character vector where matches are sought, or an object which can be coerced by [as.character()] to a character vector.
#' @param pattern a character string containing a regular expression (or [`character`] string for `fixed = TRUE`) to be matched in the given character vector. Coerced by [as.character()] to a character string if possible.  If a [`character`] vector of length 2 or more is supplied, the first element is used with a warning.
#' @param ignore.case if `FALSE`, the pattern matching is *case sensitive* and if `TRUE`, case is ignored during matching.
#' @return A [`logical`] vector
#' @name like
#' @rdname like
#' @export
#' @details Using RStudio? This function can also be inserted from the Addins menu and can have its own Keyboard Shortcut like `Ctrl+Shift+L` or `Cmd+Shift+L` (see `Tools` > `Modify Keyboard Shortcuts...`).
#' @source Idea from the [`like` function from the `data.table` package](https://github.com/Rdatatable/data.table/blob/master/R/like.R), but made it case insensitive at default and let it support multiple patterns. Also, if the regex fails the first time, it tries again with `perl = TRUE`.
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
#' b <- c("case", "diff", "yet")
#' a %like% b
#' #> TRUE TRUE TRUE
#'
#' # get frequencies of bacteria whose name start with 'Ent' or 'ent'
#' library(dplyr)
#' example_isolates %>%
#'   filter(mo_genus(mo) %like% '^ent') %>%
#'   freq(mo_fullname(mo))
like <- function(x, pattern, ignore.case = TRUE) {
  if (length(pattern) > 1) {
    if (length(x) != length(pattern)) {
      if (length(x) == 1) {
        x <- rep(x, length(pattern))
      }
      # return TRUE for every 'x' that matches any 'pattern', FALSE otherwise
      res <- sapply(pattern, function(pttrn) base::grepl(pttrn, x, ignore.case = ignore.case))
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
          res[i] <- as.integer(x[i]) %in% base::grep(pattern[i], levels(x[i]), ignore.case = ignore.case)
        } else {
          res[i] <- base::grepl(pattern[i], x[i], ignore.case = ignore.case)
        }
      }
      return(res)
    }
  }

  # the regular way how grepl works; just one pattern against one or more x
  if (is.factor(x)) {
    as.integer(x) %in% base::grep(pattern, levels(x), ignore.case = ignore.case)
  } else {
    tryCatch(base::grepl(pattern, x, ignore.case = ignore.case),
             error = function(e) ifelse(grepl("Invalid regexp", e$message),
                                        # try with perl = TRUE:
                                        return(base::grepl(pattern = pattern, x = x,
                                                                 ignore.case = ignore.case, perl = TRUE)),
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
