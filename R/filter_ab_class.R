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

#' Filter isolates on result in antimicrobial class
#'
#' Filter isolates on results in specific antimicrobial classes. This makes it easy to filter on isolates that were tested for e.g. any aminoglycoside.
#' @inheritSection lifecycle Stable lifecycle
#' @param x a data set
#' @param ab_class an antimicrobial class, like `"carbapenems"`, as can be found in [`antibiotics$group`][antibiotics]
#' @param result an antibiotic result: S, I or R (or a combination of more of them)
#' @param scope the scope to check which variables to check, can be `"any"` (default) or `"all"`
#' @param ... parameters passed on to `filter_at` from the `dplyr` package
#' @details The `group` column in [antibiotics] data set will be searched for `ab_class` (case-insensitive). If no results are found, the `atc_group1` and `atc_group2` columns will be searched. Next, `x` will be checked for column names with a value in any abbreviations, codes or official names found in the [antibiotics] data set.
#' @rdname filter_ab_class
#' @export
#' @examples
#' \dontrun{
#' library(dplyr)
#'
#' # filter on isolates that have any result for any aminoglycoside
#' example_isolates %>% filter_aminoglycosides()
#'
#' # this is essentially the same as (but without determination of column names):
#' example_isolates %>%
#'   filter_at(.vars = vars(c("GEN", "TOB", "AMK", "KAN")),
#'             .vars_predicate = any_vars(. %in% c("S", "I", "R")))
#'
#'
#' # filter on isolates that show resistance to ANY aminoglycoside
#' example_isolates %>% filter_aminoglycosides("R")
#'
#' # filter on isolates that show resistance to ALL aminoglycosides
#' example_isolates %>% filter_aminoglycosides("R", "all")
#'
#' # filter on isolates that show resistance to
#' # any aminoglycoside and any fluoroquinolone
#' example_isolates %>%
#'   filter_aminoglycosides("R") %>%
#'   filter_fluoroquinolones("R")
#'
#' # filter on isolates that show resistance to
#' # all aminoglycosides and all fluoroquinolones
#' example_isolates %>%
#'   filter_aminoglycosides("R", "all") %>%
#'   filter_fluoroquinolones("R", "all")
#' }
filter_ab_class <- function(x,
                            ab_class,
                            result = NULL,
                            scope = "any",
                            ...) {
  
  check_dataset_integrity()
  
  if (!is.data.frame(x)) {
    stop("`x` must be a data frame.", call. = FALSE)
  }
  
  # save to return later
  x_class <- class(x)
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  
  scope <- scope[1L]
  if (is.null(result)) {
    result <- c("S", "I", "R")
  }
  # make result = "SI" work too:
  result <- unlist(strsplit(result, ""))
  
  if (!all(result %in% c("S", "I", "R"))) {
    stop("`result` must be one or more of: S, I, R", call. = FALSE)
  }
  if (!all(scope %in% c("any", "all"))) {
    stop("`scope` must be one of: any, all", call. = FALSE)
  }
  
  # get only columns with class ab, mic or disk - those are AMR results
  vars_df <- colnames(x)[sapply(x, function(y) is.rsi(y) | is.mic(y) | is.disk(y))]
  vars_df_ab <- suppressWarnings(as.ab(vars_df))
  # get the columns with a group names in the chosen ab class
  vars_df <- vars_df[which(ab_group(vars_df_ab) %like% ab_class | 
                             ab_atc_group1(vars_df_ab) %like% ab_class |
                             ab_atc_group2(vars_df_ab) %like% ab_class)]
  ab_group <- find_ab_group(ab_class)
  
  if (length(vars_df) > 0) {
    if (length(result) == 1) {
      operator <- " is "
    } else {
      operator <- " is one of "
    }
    if (scope == "any") {
      scope_txt <- " or "
      scope_fn <- any
    } else {
      scope_txt <- " and "
      scope_fn <- all
      if (length(vars_df) > 1) {
        operator <- gsub("is", "are", operator)
      }
    }
    if (length(vars_df) > 1) {
      scope <- paste(scope, "of columns ")
    } else {
      scope <- "column "
    }
    message(font_blue(paste0("Filtering on ", ab_group, ": ", scope,
                             paste0(font_bold(paste0("`", vars_df, "`"), collapse = NULL), collapse = scope_txt), operator, toString(result))))
    filtered <<- as.logical(by(x, seq_len(nrow(x)),
                              function(row) scope_fn(unlist(row[, vars_df]) %in% result, na.rm = TRUE)))
    x <- x[which(filtered), , drop = FALSE]
  } else {
    message(font_blue(paste0("NOTE: no antimicrobial agents of class ", ab_group, 
                             " (such as ", find_ab_names(ab_group), 
                             ") found, data left unchanged.")))
  }
  class(x) <- x_class
  x
}

#' @rdname filter_ab_class
#' @export
filter_aminoglycosides <- function(x,
                                   result = NULL,
                                   scope = "any",
                                   ...) {
  filter_ab_class(x = x,
                  ab_class = "aminoglycoside",
                  result = result,
                  scope = scope,
                  ...)
}

#' @rdname filter_ab_class
#' @export
filter_carbapenems <- function(x,
                               result = NULL,
                               scope = "any",
                               ...) {
  filter_ab_class(x = x,
                  ab_class = "carbapenem",
                  result = result,
                  scope = scope,
                  ...)
}

#' @rdname filter_ab_class
#' @export
filter_cephalosporins <- function(x,
                                  result = NULL,
                                  scope = "any",
                                  ...) {
  filter_ab_class(x = x,
                  ab_class = "cephalosporin",
                  result = result,
                  scope = scope,
                  ...)
}

#' @rdname filter_ab_class
#' @export
filter_1st_cephalosporins <- function(x,
                                      result = NULL,
                                      scope = "any",
                                      ...) {
  filter_ab_class(x = x,
                  ab_class = "cephalosporins (1st gen.)",
                  result = result,
                  scope = scope,
                  ...)
}

#' @rdname filter_ab_class
#' @export
filter_2nd_cephalosporins <- function(x,
                                      result = NULL,
                                      scope = "any",
                                      ...) {
  filter_ab_class(x = x,
                  ab_class = "cephalosporins (2nd gen.)",
                  result = result,
                  scope = scope,
                  ...)
}

#' @rdname filter_ab_class
#' @export
filter_3rd_cephalosporins <- function(x,
                                      result = NULL,
                                      scope = "any",
                                      ...) {
  filter_ab_class(x = x,
                  ab_class = "cephalosporins (3rd gen.)",
                  result = result,
                  scope = scope,
                  ...)
}

#' @rdname filter_ab_class
#' @export
filter_4th_cephalosporins <- function(x,
                                      result = NULL,
                                      scope = "any",
                                      ...) {
  filter_ab_class(x = x,
                  ab_class = "cephalosporins (4th gen.)",
                  result = result,
                  scope = scope,
                  ...)
}

#' @rdname filter_ab_class
#' @export
filter_5th_cephalosporins <- function(x,
                                      result = NULL,
                                      scope = "any",
                                      ...) {
  filter_ab_class(x = x,
                  ab_class = "cephalosporins (5th gen.)",
                  result = result,
                  scope = scope,
                  ...)
}

#' @rdname filter_ab_class
#' @export
filter_fluoroquinolones <- function(x,
                                    result = NULL,
                                    scope = "any",
                                    ...) {
  filter_ab_class(x = x,
                  ab_class = "fluoroquinolone",
                  result = result,
                  scope = scope,
                  ...)
}

#' @rdname filter_ab_class
#' @export
filter_glycopeptides <- function(x,
                                 result = NULL,
                                 scope = "any",
                                 ...) {
  filter_ab_class(x = x,
                  ab_class = "glycopeptide",
                  result = result,
                  scope = scope,
                  ...)
}

#' @rdname filter_ab_class
#' @export
filter_macrolides <- function(x,
                              result = NULL,
                              scope = "any",
                              ...) {
  filter_ab_class(x = x,
                  ab_class = "macrolide",
                  result = result,
                  scope = scope,
                  ...)
}

#' @rdname filter_ab_class
#' @export
filter_tetracyclines <- function(x,
                                 result = NULL,
                                 scope = "any",
                                 ...) {
  filter_ab_class(x = x,
                  ab_class = "tetracycline",
                  result = result,
                  scope = scope,
                  ...)
}

find_ab_group <- function(ab_class) {
  ifelse(ab_class %in% c("aminoglycoside",
                         "carbapenem",
                         "cephalosporin",
                         "fluoroquinolone",
                         "glycopeptide",
                         "macrolide",
                         "tetracycline"),
         paste0(ab_class, "s"),
         antibiotics %>%
           subset(group %like% ab_class | 
                    atc_group1 %like% ab_class | 
                    atc_group2 %like% ab_class) %>%
           pull(group) %>%
           unique() %>%
           tolower() %>%
           paste(collapse = "/")
  )
}

find_ab_names <- function(ab_group) {
  drugs <- antibiotics[which(antibiotics$group %like% ab_group), "name"]
  paste0(ab_name(sample(drugs, size = min(4, length(drugs)), replace = FALSE),
                 tolower = TRUE, language = NULL), 
         collapse = ", ")
}
