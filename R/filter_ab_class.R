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

#' Filter isolates on result in antibiotic class
#'
#' Filter isolates on results in specific antibiotic variables based on their class (ATC groups). This makes it easy to get a list of isolates that were tested for e.g. any aminoglycoside.
#' @param tbl a data set
#' @param ab_class an antimicrobial class, like \code{"carbapenems"}. More specifically, this should be a text that can be found in a 4th level ATC group (chemical subgroup) or a 5th level ATC group (chemical substance), please see \href{https://www.whocc.no/atc/structure_and_principles/}{this explanation on the WHOCC website}.
#' @param result an antibiotic result: S, I or R (or a combination of more of them)
#' @param scope the scope to check which variables to check, can be \code{"any"} (default) or \code{"all"}
#' @param ... parameters passed on to \code{\link[dplyr]{filter_at}}
#' @details The \code{\link{antibiotics}} data set will be searched for \code{ab_class} in the columns \code{atc_group1} and \code{atc_group2} (case-insensitive). Next, \code{tbl} will be checked for column names with a value in any abbreviations, codes or official names found in the \code{antibiotics} data set.
#' @rdname filter_ab_class
#' @keywords filter fillter_class
#' @importFrom dplyr filter_at %>% select vars any_vars all_vars
#' @importFrom crayon bold blue
#' @export
#' @examples
#' library(dplyr)
#'
#' # filter on isolates that have any result for any aminoglycoside
#' septic_patients %>% filter_aminoglycosides()
#'
#' # this is essentially the same as:
#' septic_patients %>%
#'   filter_at(.vars = vars(c("gent", "tobr", "amik", "kana")),
#'             .vars_predicate = any_vars(. %in% c("S", "I", "R")))
#'
#'
#' # filter on isolates that show resistance to ANY aminoglycoside
#' septic_patients %>% filter_aminoglycosides("R")
#'
#' # filter on isolates that show resistance to ALL aminoglycosides
#' septic_patients %>% filter_aminoglycosides("R", "all")
#'
#' # filter on isolates that show resistance to
#' # any aminoglycoside and any fluoroquinolone
#' septic_patients %>%
#'   filter_aminoglycosides("R") %>%
#'   filter_fluoroquinolones("R")
#'
#' # filter on isolates that show resistance to
#' # all aminoglycosides and all fluoroquinolones
#' septic_patients %>%
#'   filter_aminoglycosides("R", "all") %>%
#'   filter_fluoroquinolones("R", "all")
filter_ab_class <- function(tbl,
                            ab_class,
                            result = NULL,
                            scope = "any",
                            ...) {
  scope <- scope[1L]
  if (is.null(result)) {
    result <- c("S", "I", "R")
  }
  # make result = "IR" work too:
  result <- unlist(strsplit(result, ""))

  if (!all(result %in% c("S", "I", "R"))) {
    stop("`result` must be one or more of: S, I, R", call. = FALSE)
  }
  if (!all(scope %in% c("any", "all"))) {
    stop("`scope` must be one of: any, all", call. = FALSE)
  }

  vars_df <- colnames(tbl)[tolower(colnames(tbl)) %in% tolower(ab_class_vars(ab_class))]
  atc_groups <- ab_class_atcgroups(ab_class)

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
      scope <- paste(scope, "of ")
    } else {
      scope <- ""
    }
    message(blue(paste0("Filtering on ", atc_groups, ": ", scope,
                        paste(bold(paste0("`", vars_df, "`")), collapse = scope_txt), operator, toString(result))))
    tbl %>%
      filter_at(vars(vars_df),
                scope_fn(. %in% result),
                ...)
  } else {
    warning(paste0("no antibiotics of class ", atc_groups, " found, leaving data unchanged"), call. = FALSE)
    tbl
  }
}

#' @rdname filter_ab_class
#' @export
filter_aminoglycosides <- function(tbl,
                                   result = NULL,
                                   scope = "any",
                                   ...) {
  filter_ab_class(tbl = tbl,
                  ab_class = "aminoglycoside",
                  result = result,
                  scope = scope,
                  ...)
}

#' @rdname filter_ab_class
#' @export
filter_carbapenems <- function(tbl,
                               result = NULL,
                               scope = "any",
                               ...) {
  filter_ab_class(tbl = tbl,
                  ab_class = "carbapenem",
                  result = result,
                  scope = scope,
                  ...)
}

#' @rdname filter_ab_class
#' @export
filter_cephalosporins <- function(tbl,
                                  result = NULL,
                                  scope = "any",
                                  ...) {
  filter_ab_class(tbl = tbl,
                  ab_class = "cephalosporin",
                  result = result,
                  scope = scope,
                  ...)
}

#' @rdname filter_ab_class
#' @export
filter_1st_cephalosporins <- function(tbl,
                                      result = NULL,
                                      scope = "any",
                                      ...) {
  filter_ab_class(tbl = tbl,
                  ab_class = "first-generation cephalosporin",
                  result = result,
                  scope = scope,
                  ...)
}

#' @rdname filter_ab_class
#' @export
filter_2nd_cephalosporins <- function(tbl,
                                      result = NULL,
                                      scope = "any",
                                      ...) {
  filter_ab_class(tbl = tbl,
                  ab_class = "second-generation cephalosporin",
                  result = result,
                  scope = scope,
                  ...)
}

#' @rdname filter_ab_class
#' @export
filter_3rd_cephalosporins <- function(tbl,
                                      result = NULL,
                                      scope = "any",
                                      ...) {
  filter_ab_class(tbl = tbl,
                  ab_class = "third-generation cephalosporin",
                  result = result,
                  scope = scope,
                  ...)
}

#' @rdname filter_ab_class
#' @export
filter_4th_cephalosporins <- function(tbl,
                                      result = NULL,
                                      scope = "any",
                                      ...) {
  filter_ab_class(tbl = tbl,
                  ab_class = "fourth-generation cephalosporin",
                  result = result,
                  scope = scope,
                  ...)
}

#' @rdname filter_ab_class
#' @export
filter_fluoroquinolones <- function(tbl,
                                    result = NULL,
                                    scope = "any",
                                    ...) {
  filter_ab_class(tbl = tbl,
                  ab_class = "fluoroquinolone",
                  result = result,
                  scope = scope,
                  ...)
}

#' @rdname filter_ab_class
#' @export
filter_glycopeptides <- function(tbl,
                                 result = NULL,
                                 scope = "any",
                                 ...) {
  filter_ab_class(tbl = tbl,
                  ab_class = "glycopeptide",
                  result = result,
                  scope = scope,
                  ...)
}

#' @rdname filter_ab_class
#' @export
filter_macrolides <- function(tbl,
                              result = NULL,
                              scope = "any",
                              ...) {
  filter_ab_class(tbl = tbl,
                  ab_class = "macrolide",
                  result = result,
                  scope = scope,
                  ...)
}

#' @rdname filter_ab_class
#' @export
filter_tetracyclines <- function(tbl,
                                 result = NULL,
                                 scope = "any",
                                 ...) {
  filter_ab_class(tbl = tbl,
                  ab_class = "tetracycline",
                  result = result,
                  scope = scope,
                  ...)
}

#' @importFrom dplyr %>% filter_at vars any_vars select
ab_class_vars <- function(ab_class) {
  ab_vars <- AMR::antibiotics %>%
    filter_at(vars(c("atc_group1", "atc_group2")), any_vars(. %like% ab_class)) %>%
    select(atc:trade_name) %>%
    as.matrix() %>%
    as.character() %>%
    paste(collapse = "|") %>%
    strsplit("|", fixed = TRUE) %>%
    unlist() %>%
    unique()
  ab_vars[!is.na(ab_vars)]
}

#' @importFrom dplyr %>% filter pull
ab_class_atcgroups <- function(ab_class) {
  ifelse(ab_class %in% c("aminoglycoside",
                         "carbapenem",
                         "cephalosporin",
                         "first-generation cephalosporin",
                         "second-generation cephalosporin",
                         "third-generation cephalosporin",
                         "fourth-generation cephalosporin",
                         "fluoroquinolone",
                         "glycopeptide",
                         "macrolide",
                         "tetracycline"),
         paste0(ab_class, "s"),
         AMR::antibiotics %>%
           filter(atc %in% ab_class_vars(ab_class)) %>%
           pull("atc_group2") %>%
           unique() %>%
           tolower() %>%
           paste(collapse = "/")
  )
}
