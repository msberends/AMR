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
#' Convenient wrapper around \code{\link[base]{grep}} to match a pattern: \code{a \%like\% b}. It always returns a \code{logical} vector and is always case-insensitive. Also, \code{pattern} (\code{b}) can be as long as \code{x} (\code{a}) to compare items of each index in both vectors.
#' @inheritParams base::grepl
#' @return A \code{logical} vector
#' @name like
#' @rdname like
#' @export
#' @details Using RStudio? This function can also be inserted from the Addins menu and can have its own Keyboard Shortcut like Ctrl+Shift+L or Cmd+Shift+L (see Tools > Modify Keyboard Shortcuts...).
#' @source Idea from the \href{https://github.com/Rdatatable/data.table/blob/master/R/like.R}{\code{like} function from the \code{data.table} package}, but made it case insensitive at default and let it support multiple patterns.
#' @seealso \code{\link[base]{grep}}
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
#' library(clean)
#' septic_patients %>%
#'   left_join_microorganisms() %>%
#'   filter(genus %like% '^ent') %>%
#'   freq(genus, species)
like <- function(x, pattern) {
  if (length(pattern) > 1) {
    if (length(x) != length(pattern)) {
      pattern <- pattern[1]
      warning('only the first element of argument `pattern` used for `%like%`', call. = TRUE)
    } else {
      # x and pattern are of same length, so items with each other
      res <- vector(length = length(pattern))
      for (i in 1:length(res)) {
        if (is.factor(x[i])) {
          res[i] <- as.integer(x[i]) %in% base::grep(pattern[i], levels(x[i]), ignore.case = TRUE)
        } else {
          res[i] <- base::grepl(pattern[i], x[i], ignore.case = TRUE)
        }
      }
      return(res)
    }
  }

  # the regular way how grepl works; just one pattern against one or more x
  if (is.factor(x)) {
    as.integer(x) %in% base::grep(pattern, levels(x), ignore.case = TRUE)
  } else {
    tryCatch(base::grepl(pattern, x, ignore.case = TRUE),
             error = function(e) ifelse(test = grepl("Invalid regexp", e$message),
                                        # try with perl = TRUE:
                                        yes = return(base::grepl(pattern, x, ignore.case = TRUE, perl = TRUE)),
                                        no = stop(e$message)))
  }
}

#' @rdname like
#' @export
"%like%" <- like
