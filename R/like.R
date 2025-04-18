# ==================================================================== #
# TITLE:                                                               #
# AMR: An R Package for Working with Antimicrobial Resistance Data     #
#                                                                      #
# SOURCE CODE:                                                         #
# https://github.com/msberends/AMR                                     #
#                                                                      #
# PLEASE CITE THIS SOFTWARE AS:                                        #
# Berends MS, Luz CF, Friedrich AW, et al. (2022).                     #
# AMR: An R Package for Working with Antimicrobial Resistance Data.    #
# Journal of Statistical Software, 104(3), 1-31.                       #
# https://doi.org/10.18637/jss.v104.i03                                #
#                                                                      #
# Developed at the University of Groningen and the University Medical  #
# Center Groningen in The Netherlands, in collaboration with many      #
# colleagues from around the world, see our website.                   #
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
# how to conduct AMR data analysis: https://amr-for-r.org              #
# ==================================================================== #

#' Vectorised Pattern Matching with Keyboard Shortcut
#'
#' Convenient wrapper around [grepl()] to match a pattern: `x %like% pattern`. It always returns a [`logical`] vector and is always case-insensitive (use `x %like_case% pattern` for case-sensitive matching). Also, `pattern` can be as long as `x` to compare items of each index in both vectors, or they both can have the same length to iterate over all cases.
#' @param x A [character] vector where matches are sought, or an object which can be coerced by [as.character()] to a [character] vector.
#' @param pattern A [character] vector containing regular expressions (or a [character] string for `fixed = TRUE`) to be matched in the given [character] vector. Coerced by [as.character()] to a [character] string if possible.
#' @param ignore.case If `FALSE`, the pattern matching is *case sensitive* and if `TRUE`, case is ignored during matching.
#' @return A [logical] vector
#' @name like
#' @rdname like
#' @export
#' @details
#' These [like()] and `%like%`/`%unlike%` functions:
#' * Are case-insensitive (use `%like_case%`/`%unlike_case%` for case-sensitive matching)
#' * Support multiple patterns
#' * Check if `pattern` is a valid regular expression and sets `fixed = TRUE` if not, to greatly improve speed (vectorised over `pattern`)
#' * Always use compatibility with Perl unless `fixed = TRUE`, to greatly improve speed
#'
#' Using RStudio? The `%like%`/`%unlike%` functions can also be directly inserted in your code from the Addins menu and can have its own keyboard shortcut like `Shift+Ctrl+L` or `Shift+Cmd+L` (see menu `Tools` > `Modify Keyboard Shortcuts...`). If you keep pressing your shortcut, the inserted text will be iterated over `%like%` -> `%unlike%` -> `%like_case%` -> `%unlike_case%`.
#' @source Idea from the [`like` function from the `data.table` package](https://github.com/Rdatatable/data.table/blob/ec1259af1bf13fc0c96a1d3f9e84d55d8106a9a4/R/like.R), although altered as explained in *Details*.
#' @seealso [grepl()]

#' @examples
#' # data.table has a more limited version of %like%, so unload it:
#' try(detach("package:data.table", unload = TRUE), silent = TRUE)
#'
#' a <- "This is a test"
#' b <- "TEST"
#' a %like% b
#' b %like% a
#'
#' # also supports multiple patterns
#' a <- c("Test case", "Something different", "Yet another thing")
#' b <- c("case", "diff", "yet")
#' a %like% b
#' a %unlike% b
#'
#' a[1] %like% b
#' a %like% b[1]
#'
#' \donttest{
#' # get isolates whose name start with 'Entero' (case-insensitive)
#' example_isolates[which(mo_name() %like% "^entero"), ]
#'
#' if (require("dplyr")) {
#'   example_isolates %>%
#'     filter(mo_name() %like% "^ent")
#' }
#' }
like <- function(x, pattern, ignore.case = TRUE) {
  meet_criteria(x, allow_NA = TRUE)
  meet_criteria(pattern, allow_NA = FALSE)
  meet_criteria(ignore.case, allow_class = "logical", has_length = 1)

  if (all(is.na(x))) {
    return(rep(FALSE, length(x)))
  }

  # set to fixed if no valid regex (vectorised)
  fixed <- !is_valid_regex(pattern)

  if (ignore.case == TRUE) {
    # set here, otherwise if fixed = TRUE, this warning will be thrown: argument `ignore.case = TRUE` will be ignored
    x <- tolower(x)
    pattern <- tolower(pattern)
  }

  if (is.factor(x)) {
    x <- as.character(x)
  }

  if (length(pattern) == 1) {
    grepl(pattern, x, ignore.case = FALSE, fixed = fixed, perl = !fixed)
  } else {
    if (length(x) == 1) {
      x <- rep(x, length(pattern))
    } else if (length(pattern) != length(x)) {
      stop_(
        "arguments `x` and `pattern` must be of same length, or either one must be 1 ",
        "(`x` has length ", length(x), " and `pattern` has length ", length(pattern), ")"
      )
    }
    unlist(
      Map(
        f = grepl,
        x = x,
        pattern = pattern,
        fixed = fixed,
        perl = !fixed,
        MoreArgs = list(ignore.case = FALSE),
        USE.NAMES = FALSE
      )
    )
  }
}

#' @rdname like
#' @export
"%like%" <- function(x, pattern) {
  like(x, pattern, ignore.case = TRUE)
}

#' @rdname like
#' @export
"%unlike%" <- function(x, pattern) {
  !like(x, pattern, ignore.case = TRUE)
}

#' @rdname like
#' @export
"%like_case%" <- function(x, pattern) {
  like(x, pattern, ignore.case = FALSE)
}

#' @rdname like
#' @export
"%unlike_case%" <- function(x, pattern) {
  !like(x, pattern, ignore.case = FALSE)
}
