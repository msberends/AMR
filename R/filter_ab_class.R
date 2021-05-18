# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Data Analysis for R                   #
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
# how to conduct AMR data analysis: https://msberends.github.io/AMR/   #
# ==================================================================== #

#' Filter Isolates on Result in Antimicrobial Class
#'
#' Filter isolates on results in specific antimicrobial classes. This makes it easy to filter on isolates that were tested for e.g. any aminoglycoside, or to filter on carbapenem-resistant isolates without the need to specify the drugs.
#' @inheritSection lifecycle Stable Lifecycle
#' @param x a data set
#' @param ab_class an antimicrobial class, like `"carbapenems"`. The columns `group`, `atc_group1` and `atc_group2` of the [antibiotics] data set will be searched (case-insensitive) for this value.
#' @param result an antibiotic result: S, I or R (or a combination of more of them)
#' @param scope the scope to check which variables to check, can be `"any"` (default) or `"all"`
#' @param only_rsi_columns a [logical] to indicate whether only columns must be included that were transformed to class `<rsi>` (see [as.rsi()]) on beforehand (defaults to `FALSE`)
#' @param ... arguments passed on to [filter_ab_class()]
#' @details All columns of `x` will be searched for known antibiotic names, abbreviations, brand names and codes (ATC, EARS-Net, WHO, etc.). This means that a filter function like e.g. [filter_aminoglycosides()] will include column names like 'gen', 'genta', 'J01GB03', 'tobra', 'Tobracin', etc.
#' 
#' The group of betalactams consists of all carbapenems, cephalosporins and penicillins.
#' @rdname filter_ab_class
#' @seealso [antibiotic_class_selectors()] for the `select()` equivalent.
#' @export
#' @examples
#' x <- filter_carbapenems(example_isolates)
#' \donttest{
#' # base R filter options (requires R >= 3.2)
#' example_isolates[filter_carbapenems(), ]
#' example_isolates[which(filter_carbapenems() & mo_is_gram_negative()), ]
#' 
#' if (require("dplyr")) {
#'
#'   # filter on isolates that have any result for any aminoglycoside
#'   example_isolates %>% filter_aminoglycosides()
#'   example_isolates %>% filter_ab_class("aminoglycoside")
#' 
#'   # this is essentially the same as (but without determination of column names):
#'   example_isolates %>%
#'     filter_at(.vars = vars(c("GEN", "TOB", "AMK", "KAN")),
#'               .vars_predicate = any_vars(. %in% c("S", "I", "R")))
#'
#'
#'   # filter on isolates that show resistance to ANY aminoglycoside
#'   example_isolates %>% filter_aminoglycosides("R", "any")
#'  
#'   # filter on isolates that show resistance to ALL aminoglycosides
#'   example_isolates %>% filter_aminoglycosides("R", "all")
#'  
#'   # filter on isolates that show resistance to
#'   # any aminoglycoside and any fluoroquinolone
#'   example_isolates %>%
#'     filter_aminoglycosides("R") %>%
#'     filter_fluoroquinolones("R")
#'  
#'   # filter on isolates that show resistance to
#'   # all aminoglycosides and all fluoroquinolones
#'   example_isolates %>%
#'     filter_aminoglycosides("R", "all") %>%
#'     filter_fluoroquinolones("R", "all")
#'   
#'   # with dplyr 1.0.0 and higher (that adds 'across()'), this is all equal:
#'   # (though the row names on the first are more correct)
#'   example_isolates %>% filter_carbapenems("R", "all")
#'   example_isolates %>% filter(across(carbapenems(), ~. == "R"))
#'   example_isolates %>% filter(across(carbapenems(), function(x) x == "R"))
#'   example_isolates %>% filter(filter_carbapenems("R", "all"))
#' }
#' }
filter_ab_class <- function(x,
                            ab_class,
                            result = NULL,
                            scope = "any",
                            only_rsi_columns = FALSE,
                            ...) {
  .call_depth <- list(...)$`.call_depth`
  if (is.null(.call_depth)) {
    .call_depth <- 0
  }
  .fn <- list(...)$`.fn`
  if (is.null(.fn)) {
    .fn <- "filter_ab_class"
  }
  
  return_only_row_indices <- FALSE
  if (missing(x) || is_null_or_grouped_tbl(x)) {
    # when `x` is left blank, auto determine it (get_current_data() also contains dplyr::cur_data_all())
    # is also fix for using a grouped df as input (a dot as first argument)
    x <- get_current_data(arg_name = "x", call = -2 - .call_depth)
    return_only_row_indices <- TRUE
  }
  meet_criteria(x, allow_class = "data.frame", .call_depth = .call_depth)
  meet_criteria(ab_class, allow_class = "character", has_length = 1, .call_depth = .call_depth)
  if (!is.null(result)) {
    result <- toupper(result)
  }
  meet_criteria(result, allow_class = "character", has_length = c(1, 2, 3), is_in = c("S", "I", "R"), allow_NULL = TRUE, .call_depth = .call_depth)
  meet_criteria(scope, allow_class = "character", has_length = 1, is_in = c("all", "any"), .call_depth = .call_depth)
  meet_criteria(only_rsi_columns, allow_class = "logical", has_length = 1, .call_depth = .call_depth)
  
  check_dataset_integrity()
  
  # save to return later
  x.bak <- x
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  
  if (is.null(result)) {
    result <- c("S", "I", "R")
  }
  # make result = "SI" works too:
  result <- unlist(strsplit(result, ""))
  
  # get all columns in data with names that resemble antibiotics
  ab_in_data <- get_column_abx(x, info = FALSE, only_rsi_columns = only_rsi_columns, sort = FALSE)
  
  # improve speed here so it will only run once when e.g. in one select call
  if (!identical(pkg_env$filter_ab_selector, unique_call_id())) {
    ab_in_data <- get_column_abx(x, info = FALSE, only_rsi_columns = only_rsi_columns, sort = FALSE)
    pkg_env$filter_ab_selector <- unique_call_id()
    pkg_env$filter_ab_selector_cols <- ab_in_data
  } else {
    ab_in_data <- pkg_env$filter_ab_selector_cols
  }
  
  if (length(ab_in_data) == 0) {
    message_("No columns with antibiotic test results found (see ?as.rsi), data left unchanged.")
    return(x.bak)
  }
  # get reference data
  ab_class.bak <- ab_class
  ab_class <- gsub("[^a-zA-Z|0-9]+", ".*", ab_class)
  ab_class <- gsub("(ph|f)", "(ph|f)", ab_class)
  ab_class <- gsub("(t|th)", "(t|th)", ab_class)
  ab_reference <- subset(antibiotics,
                         group %like% ab_class | 
                           atc_group1 %like% ab_class | 
                           atc_group2 %like% ab_class)
  if (nrow(ab_reference) == 0) {
    message_("Unknown antimicrobial class '", ab_class.bak, "', data left unchanged.")
    return(x.bak)
  }
  ab_group <- find_ab_group(ab_class.bak)
  # get the columns with a group names in the chosen ab class
  agents <- ab_in_data[names(ab_in_data) %in% ab_reference$ab]
  if (length(agents) == 0) {
    message_("No antimicrobial agents of class '", ab_group, 
             "' found (such as ", find_ab_names(ab_class, 2), 
             ")",
             ifelse(only_rsi_columns == TRUE, " with class <rsi>,", ","),
             " data left unchanged.")
    return(x.bak)
  }
  
  if (scope == "any") {
    scope_txt <- " or "
    scope_fn <- any
  } else {
    scope_txt <- " and "
    scope_fn <- all
  }
  if (length(agents) > 1) {
    operator <- " are"
    scope <- paste("values in", scope, "of columns ")
  } else {
    operator <- " is"
    scope <- "value in column "
  }
  if (length(result) > 1) {
    operator <- paste(operator, "either")
  }
  
  # sort columns on official name
  agents <- agents[order(ab_name(names(agents), language = NULL))]
  
  agents_formatted <- paste0("'", font_bold(agents, collapse = NULL), "'")
  agents_names <- ab_name(names(agents), tolower = TRUE, language = NULL)
  need_name <- tolower(gsub("[^a-zA-Z]", "", agents)) != tolower(gsub("[^a-zA-Z]", "", agents_names))
  agents_formatted[need_name] <- paste0(agents_formatted[need_name],
                                        " (", agents_names[need_name], ")")
  
  message_("Applying `", .fn, "()`: ", scope,
           vector_or(agents_formatted, quotes = FALSE, last_sep = scope_txt),
           operator, " ", vector_or(result, quotes = TRUE))
  x_transposed <- as.list(as.data.frame(t(x[, agents, drop = FALSE]), stringsAsFactors = FALSE))
  filtered <- vapply(FUN.VALUE = logical(1), x_transposed, function(y) scope_fn(y %in% result, na.rm = TRUE))
  
  if (return_only_row_indices == TRUE) {
    filtered
  } else {
    # this returns the original data with the filtering, also preserving attributes (such as dplyr groups)
    x.bak[which(filtered), , drop = FALSE]
  }
}

#' @rdname filter_ab_class
#' @export
filter_aminoglycosides <- function(x,
                                   result = NULL,
                                   scope = "any",
                                   only_rsi_columns = FALSE,
                                   ...) {
  filter_ab_class(x = x,
                  ab_class = "aminoglycoside",
                  result = result,
                  scope = scope,
                  only_rsi_columns = only_rsi_columns,
                  .call_depth = 1,
                  .fn = "filter_aminoglycosides",
                  ...)
}

#' @rdname filter_ab_class
#' @export
filter_betalactams <- function(x,
                               result = NULL,
                               scope = "any",
                               only_rsi_columns = FALSE,
                               ...) {
  filter_ab_class(x = x,
                  ab_class = "carbapenem|cephalosporin|penicillin",
                  result = result,
                  scope = scope,
                  only_rsi_columns = only_rsi_columns,
                  .call_depth = 1,
                  .fn = "filter_betalactams",
                  ...)
}
#' @rdname filter_ab_class
#' @export
filter_carbapenems <- function(x,
                               result = NULL,
                               scope = "any",
                               only_rsi_columns = FALSE,
                               ...) {
  filter_ab_class(x = x,
                  ab_class = "carbapenem",
                  result = result,
                  scope = scope,
                  only_rsi_columns = only_rsi_columns,
                  .call_depth = 1,
                  .fn = "filter_carbapenems",
                  ...)
}

#' @rdname filter_ab_class
#' @export
filter_cephalosporins <- function(x,
                                  result = NULL,
                                  scope = "any",
                                  only_rsi_columns = FALSE,
                                  ...) {
  filter_ab_class(x = x,
                  ab_class = "cephalosporin",
                  result = result,
                  scope = scope,
                  only_rsi_columns = only_rsi_columns,
                  .call_depth = 1,
                  .fn = "filter_cephalosporins",
                  ...)
}

#' @rdname filter_ab_class
#' @export
filter_1st_cephalosporins <- function(x,
                                      result = NULL,
                                      scope = "any",
                                      only_rsi_columns = FALSE,
                                      ...) {
  filter_ab_class(x = x,
                  ab_class = "cephalosporins (1st gen.)",
                  result = result,
                  scope = scope,
                  only_rsi_columns = only_rsi_columns,
                  .call_depth = 1,
                  .fn = "filter_1st_cephalosporins",
                  ...)
}

#' @rdname filter_ab_class
#' @export
filter_2nd_cephalosporins <- function(x,
                                      result = NULL,
                                      scope = "any",
                                      only_rsi_columns = FALSE,
                                      ...) {
  filter_ab_class(x = x,
                  ab_class = "cephalosporins (2nd gen.)",
                  result = result,
                  scope = scope,
                  only_rsi_columns = only_rsi_columns,
                  .call_depth = 1,
                  .fn = "filter_2nd_cephalosporins",
                  ...)
}

#' @rdname filter_ab_class
#' @export
filter_3rd_cephalosporins <- function(x,
                                      result = NULL,
                                      scope = "any",
                                      only_rsi_columns = FALSE,
                                      ...) {
  filter_ab_class(x = x,
                  ab_class = "cephalosporins (3rd gen.)",
                  result = result,
                  scope = scope,
                  only_rsi_columns = only_rsi_columns,
                  .call_depth = 1,
                  .fn = "filter_3rd_cephalosporins",
                  ...)
}

#' @rdname filter_ab_class
#' @export
filter_4th_cephalosporins <- function(x,
                                      result = NULL,
                                      scope = "any",
                                      only_rsi_columns = FALSE,
                                      ...) {
  filter_ab_class(x = x,
                  ab_class = "cephalosporins (4th gen.)",
                  result = result,
                  scope = scope,
                  only_rsi_columns = only_rsi_columns,
                  .call_depth = 1,
                  .fn = "filter_4th_cephalosporins",
                  ...)
}

#' @rdname filter_ab_class
#' @export
filter_5th_cephalosporins <- function(x,
                                      result = NULL,
                                      scope = "any",
                                      only_rsi_columns = FALSE,
                                      ...) {
  filter_ab_class(x = x,
                  ab_class = "cephalosporins (5th gen.)",
                  result = result,
                  scope = scope,
                  only_rsi_columns = only_rsi_columns,
                  .call_depth = 1,
                  .fn = "filter_5th_cephalosporins",
                  ...)
}

#' @rdname filter_ab_class
#' @export
filter_fluoroquinolones <- function(x,
                                    result = NULL,
                                    scope = "any",
                                    only_rsi_columns = FALSE,
                                    ...) {
  filter_ab_class(x = x,
                  ab_class = "fluoroquinolone",
                  result = result,
                  scope = scope,
                  only_rsi_columns = only_rsi_columns,
                  .call_depth = 1,
                  .fn = "filter_fluoroquinolones",
                  ...)
}

#' @rdname filter_ab_class
#' @export
filter_glycopeptides <- function(x,
                                 result = NULL,
                                 scope = "any",
                                 only_rsi_columns = FALSE,
                                 ...) {
  filter_ab_class(x = x,
                  ab_class = "glycopeptide",
                  result = result,
                  scope = scope,
                  only_rsi_columns = only_rsi_columns,
                  .call_depth = 1,
                  .fn = "filter_glycopeptides",
                  ...)
}

#' @rdname filter_ab_class
#' @export
filter_macrolides <- function(x,
                              result = NULL,
                              scope = "any",
                              only_rsi_columns = FALSE,
                              ...) {
  filter_ab_class(x = x,
                  ab_class = "macrolide",
                  result = result,
                  scope = scope,
                  only_rsi_columns = only_rsi_columns,
                  .call_depth = 1,
                  .fn = "filter_macrolides",
                  ...)
}

#' @rdname filter_ab_class
#' @export
filter_oxazolidinones <- function(x,
                                  result = NULL,
                                  scope = "any",
                                  only_rsi_columns = FALSE,
                                  ...) {
  filter_ab_class(x = x,
                  ab_class = "oxazolidinone",
                  result = result,
                  scope = scope,
                  only_rsi_columns = only_rsi_columns,
                  .call_depth = 1,
                  .fn = "filter_oxazolidinones",
                  ...)
}

#' @rdname filter_ab_class
#' @export
filter_penicillins <- function(x,
                               result = NULL,
                               scope = "any",
                               only_rsi_columns = FALSE,
                               ...) {
  filter_ab_class(x = x,
                  ab_class = "penicillin",
                  result = result,
                  scope = scope,
                  only_rsi_columns = only_rsi_columns,
                  .call_depth = 1,
                  .fn = "filter_penicillins",
                  ...)
}

#' @rdname filter_ab_class
#' @export
filter_tetracyclines <- function(x,
                                 result = NULL,
                                 scope = "any",
                                 only_rsi_columns = FALSE,
                                 ...) {
  filter_ab_class(x = x,
                  ab_class = "tetracycline",
                  result = result,
                  scope = scope,
                  only_rsi_columns = only_rsi_columns,
                  .call_depth = 1,
                  .fn = "filter_tetracyclines",
                  ...)
}

find_ab_group <- function(ab_class) {
  ab_class[ab_class == "carbapenem|cephalosporin|penicillin"] <- "betalactam"
  ab_class <- gsub("[^a-zA-Z0-9]", ".*", ab_class)
  ifelse(ab_class %in% c("aminoglycoside",
                         "betalactam",
                         "carbapenem",
                         "cephalosporin",
                         "fluoroquinolone",
                         "glycopeptide",
                         "macrolide",
                         "oxazolidinone",
                         "tetracycline"),
         paste0(ab_class, "s"),
         antibiotics %pm>%
           subset(group %like% ab_class | 
                    atc_group1 %like% ab_class | 
                    atc_group2 %like% ab_class) %pm>%
           pm_pull(group) %pm>%
           unique() %pm>%
           tolower() %pm>%
           sort() %pm>% 
           paste(collapse = "/")
  )
}

find_ab_names <- function(ab_group, n = 3) {
  ab_group <- gsub("[^a-zA-Z|0-9]", ".*", ab_group)
  
  # try popular first, they have DDDs
  drugs <- antibiotics[which((!is.na(antibiotics$iv_ddd) | !is.na(antibiotics$oral_ddd)) &
                               antibiotics$name %unlike% " " &
                               antibiotics$group %like% ab_group &
                               antibiotics$ab %unlike% "[0-9]$"), ]$name
  if (length(drugs) < n) {
    # now try it all
    drugs <- antibiotics[which((antibiotics$group %like% ab_group |
                                  antibiotics$atc_group1 %like% ab_group |
                                  antibiotics$atc_group2 %like% ab_group) &
                                 antibiotics$ab %unlike% "[0-9]$"), ]$name
  }
  vector_or(ab_name(sample(drugs, size = min(n, length(drugs)), replace = FALSE),
                    tolower = TRUE,
                    language = NULL),
            quotes = FALSE)
}
