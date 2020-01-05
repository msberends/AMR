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

#' Filter isolates on result in antibiotic class
#'
#' Filter isolates on results in specific antibiotic variables based on their class (ATC groups). This makes it easy to get a list of isolates that were tested for e.g. any aminoglycoside.
#' @inheritSection lifecycle Stable lifecycle
#' @param x a data set
#' @param ab_class an antimicrobial class, like `"carbapenems"`, as can be found in [`AMR::antibiotics$group`][antibiotics]
#' @param result an antibiotic result: S, I or R (or a combination of more of them)
#' @param scope the scope to check which variables to check, can be `"any"` (default) or `"all"`
#' @param ... parameters passed on to `filter_at` from the `dplyr` package
#' @details The `group` column in [antibiotics] data set will be searched for `ab_class` (case-insensitive). If no results are found, the `atc_group1` and `atc_group2` columns will be searched. Next, `x` will be checked for column names with a value in any abbreviations, codes or official names found in the [antibiotics] data set.
#' @rdname filter_ab_class
#' @importFrom dplyr filter_at %>% select vars any_vars all_vars
#' @importFrom crayon bold blue
#' @export
#' @examples
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
filter_ab_class <- function(x,
                            ab_class,
                            result = NULL,
                            scope = "any",
                            ...) {
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

  vars_df <- colnames(x)[tolower(colnames(x)) %in% tolower(ab_class_vars(ab_class))]
  ab_group <- find_ab_group(ab_class)

  if (length(vars_df) > 0) {
    if (length(result) == 1) {
      operator <- " is "
    } else {
      operator <- " is one of "
    }
    if (scope == "any") {
      scope_txt <- " or "
      scope_fn <- any_vars
    } else {
      scope_txt <- " and "
      scope_fn <- all_vars
      if (length(vars_df) > 1) {
        operator <- gsub("is", "are", operator)
      }
    }
    if (length(vars_df) > 1) {
      scope <- paste(scope, "of columns ")
    } else {
      scope <- "column "
    }
    message(blue(paste0("Filtering on ", ab_group, ": ", scope,
                        paste(bold(paste0("`", vars_df, "`")), collapse = scope_txt), operator, toString(result))))
    x %>%
      filter_at(vars(vars_df),
                scope_fn(. %in% result),
                ...)
  } else {
    warning(paste0("no antibiotics of class ", ab_group, " found, leaving data unchanged"), call. = FALSE)
    x
  }
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

#' @importFrom dplyr %>% filter_at vars any_vars select
ab_class_vars <- function(ab_class) {
  ab_class <- gsub("[^a-z0-9]+", ".*", ab_class)
  ab_vars <- AMR::antibiotics %>%
    filter(group %like% ab_class) %>% 
    select(ab:name, abbreviations, synonyms) %>%
    unlist() %>%
    as.matrix() %>%
    as.character() %>%
    paste(collapse = "|") %>%
    strsplit("|", fixed = TRUE) %>%
    unlist() %>%
    unique()
  ab_vars <- ab_vars[!ab_vars %in% c(NA, "", "NA") & nchar(ab_vars) > 2]
  if (length(ab_vars) == 0) {
    # try again, searching atc_group1 and atc_group2 columns
    ab_vars <- AMR::antibiotics %>%
      filter_at(vars(c("atc_group1", "atc_group2")), any_vars(. %like% ab_class)) %>% 
      select(ab:name, abbreviations, synonyms) %>%
      unlist() %>%
      as.matrix() %>%
      as.character() %>%
      paste(collapse = "|") %>%
      strsplit("|", fixed = TRUE) %>%
      unlist() %>%
      unique()
    ab_vars <- ab_vars[!ab_vars %in% c(NA, "", "NA") & nchar(ab_vars) > 2]
  }
  ab_vars
}

#' @importFrom dplyr %>% filter pull
find_ab_group <- function(ab_class) {
  ifelse(ab_class %in% c("aminoglycoside",
                         "carbapenem",
                         "cephalosporin",
                         "fluoroquinolone",
                         "glycopeptide",
                         "macrolide",
                         "tetracycline"),
         paste0(ab_class, "s"),
         AMR::antibiotics %>%
           filter(ab %in% ab_class_vars(ab_class)) %>%
           pull(group) %>%
           unique() %>%
           tolower() %>%
           paste(collapse = "/")
  )
}
