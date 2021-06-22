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

#' Deprecated Functions
#'
#' These functions are so-called '[Deprecated]'. **They will be removed in a future release.** Using the functions will give a warning with the name of the function it has been replaced by (if there is one).
#' @details All antibiotic class selectors (such as [carbapenems()], [aminoglycosides()]) can now be used for filtering as well, making all their accompanying `filter_*()` functions redundant (such as [filter_carbapenems()], [filter_aminoglycosides()]).
#' @inheritSection lifecycle Retired Lifecycle
#' @inheritSection AMR Read more on Our Website!
#' @keywords internal
#' @name AMR-deprecated
#' @export
p_symbol <- function(p, emptychar = " ") {
  .Deprecated(package = "AMR", new = "cleaner::p_symbol")
  
  p <- as.double(p)
  s <- rep(NA_character_, length(p))
  
  s[p <= 1] <- emptychar
  s[p <= 0.100] <- "."
  s[p <= 0.050] <- "*"
  s[p <= 0.010] <- "**"
  s[p <= 0.001] <- "***"
  
  s
}

#' @name AMR-deprecated
#' @export
filter_first_weighted_isolate <- function(x = NULL,
                                          col_date = NULL,
                                          col_patient_id = NULL,
                                          col_mo = NULL,
                                          ...) {
  
  .Deprecated(old = "filter_first_weighted_isolate()",
              new = "filter_first_isolate()",
              package = "AMR")
  
  if (is_null_or_grouped_tbl(x)) {
    # when `x` is left blank, auto determine it (get_current_data() also contains dplyr::cur_data_all())
    # is also fix for using a grouped df as input (a dot as first argument)
    x <- tryCatch(get_current_data(arg_name = "x", call = -2), error = function(e) x)
  }
  meet_criteria(x, allow_class = "data.frame") # also checks dimensions to be >0
  meet_criteria(col_date, allow_class = "character", has_length = 1, allow_NULL = TRUE, is_in = colnames(x))
  meet_criteria(col_patient_id, allow_class = "character", has_length = 1, allow_NULL = TRUE, is_in = colnames(x))
  meet_criteria(col_mo, allow_class = "character", has_length = 1, allow_NULL = TRUE, is_in = colnames(x))
  
  filter_first_isolate(x = x, col_date = col_date, col_patient_id = col_patient_id, col_mo = col_mo, ...)
}

#' @name AMR-deprecated
#' @export
key_antibiotics <- function(x = NULL,
                            col_mo = NULL,
                            universal_1 = guess_ab_col(x, "amoxicillin"),
                            universal_2 = guess_ab_col(x, "amoxicillin/clavulanic acid"),
                            universal_3 = guess_ab_col(x, "cefuroxime"),
                            universal_4 = guess_ab_col(x, "piperacillin/tazobactam"),
                            universal_5 = guess_ab_col(x, "ciprofloxacin"),
                            universal_6 = guess_ab_col(x, "trimethoprim/sulfamethoxazole"),
                            GramPos_1 = guess_ab_col(x, "vancomycin"),
                            GramPos_2 = guess_ab_col(x, "teicoplanin"),
                            GramPos_3 = guess_ab_col(x, "tetracycline"),
                            GramPos_4 = guess_ab_col(x, "erythromycin"),
                            GramPos_5 = guess_ab_col(x, "oxacillin"),
                            GramPos_6 = guess_ab_col(x, "rifampin"),
                            GramNeg_1 = guess_ab_col(x, "gentamicin"),
                            GramNeg_2 = guess_ab_col(x, "tobramycin"),
                            GramNeg_3 = guess_ab_col(x, "colistin"),
                            GramNeg_4 = guess_ab_col(x, "cefotaxime"),
                            GramNeg_5 = guess_ab_col(x, "ceftazidime"),
                            GramNeg_6 = guess_ab_col(x, "meropenem"),
                            warnings = TRUE,
                            ...) {
  
  .Deprecated(old = "key_antibiotics()",
              new = "key_antimicrobials()",
              package = "AMR")
  
  if (is_null_or_grouped_tbl(x)) {
    # when `x` is left blank, auto determine it (get_current_data() also contains dplyr::cur_data_all())
    # is also fix for using a grouped df as input (a dot as first argument)
    x <- tryCatch(get_current_data(arg_name = "x", call = -2), error = function(e) x)
  }

  key_antimicrobials(x = x, 
                     col_mo = col_mo, 
                     universal = c(universal_1, universal_2, universal_3, universal_4, universal_5, universal_6),
                     gram_negative = c(GramNeg_1, GramNeg_2, GramNeg_3, GramNeg_4, GramNeg_5, GramNeg_6),
                     gram_positive = c(GramPos_1, GramPos_2, GramPos_3, GramPos_4, GramPos_5, GramPos_6),
                     antifungal = NULL,
                     only_rsi_columns = FALSE,
                     ...)
}

#' @name AMR-deprecated
#' @export
key_antibiotics_equal <- function(y,
                                  z,
                                  type = "keyantimicrobials",
                                  ignore_I = TRUE,
                                  points_threshold = 2,
                                  info = FALSE,
                                  na.rm = TRUE,
                                  ...) {
  
  .Deprecated(old = "key_antibiotics_equal()",
              new = "antimicrobials_equal()",
              package = "AMR")
  
  antimicrobials_equal(y = y,
                       z = z,
                       type = type,
                       ignore_I = ignore_I,
                       points_threshold = points_threshold,
                       info = info)
}


#' @name AMR-deprecated
#' @export
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
  .x_name <- list(...)$`.x_name`
  if (is.null(.x_name)) {
    .x_name <- deparse(substitute(x))
  }
  .fn <- list(...)$`.fn`
  if (is.null(.fn)) {
    .fn <- "filter_ab_class"
  }
  .fn_old <- .fn
  # new way: using the ab selectors
  .fn <- gsub("filter_", "", .fn, fixed = TRUE)
  .fn <- gsub("^([1-5][a-z]+)_cephalosporins", "cephalosporins_\\1", .fn)
  
  if (missing(x) || is_null_or_grouped_tbl(x)) {
    # when `x` is left blank, auto determine it (get_current_data() also contains dplyr::cur_data_all())
    # is also fix for using a grouped df as input (a dot as first argument)
    x <- get_current_data(arg_name = "x", call = -2 - .call_depth)
    .x_name <- "your_data"
  }
  meet_criteria(x, allow_class = "data.frame", .call_depth = .call_depth)
  meet_criteria(ab_class, allow_class = "character", has_length = 1, .call_depth = .call_depth)
  if (!is.null(result)) {
    # make result = "SI" works too:
    result <- toupper(unlist(strsplit(result, "")))
  }
  meet_criteria(result, allow_class = "character", has_length = c(1, 2, 3), is_in = c("S", "I", "R"), allow_NULL = TRUE, .call_depth = .call_depth)
  meet_criteria(scope, allow_class = "character", has_length = 1, is_in = c("all", "any"), .call_depth = .call_depth)
  meet_criteria(only_rsi_columns, allow_class = "logical", has_length = 1, .call_depth = .call_depth)
  
  if (is.null(result)) {
    result <- c("S", "I", "R")
  }
  
  # get e.g. carbapenems() from filter_carbapenems()
  fn <- get(.fn, envir = asNamespace("AMR"))
  if (scope == "any") {
    scope_fn <- any
  } else {
    scope_fn <- all
  }
  
  # be nice here, be VERY extensive about how the AB selectors have taken over this function
  deprecated_fn <- paste0(.fn, "(", ifelse(.fn == "ab_class", paste0("\"", ab_class, "\""), ""), ")",
                          ifelse(length(result) > 1,
                                 paste0(", c(", paste0("\"", result, "\"", collapse = ", "), ")"),
                                 ifelse(is.null(result),
                                        "",
                                        paste0(" == \"", result, "\""))))
  if (.x_name == ".") {
    .x_name <- "your_data"
  }
  warning_(paste0("`", .fn_old, "()` is deprecated. Use the antibiotic selector `", .fn, "()` instead.\n",
                  "In dplyr:\n",
                  "  - ", .x_name, " %>% filter(", scope, "(", deprecated_fn, "))\n",
                  ifelse(length(result) > 1, 
                         paste0("  - ", .x_name, " %>% filter(", scope, "(", 
                                .fn, "(", ifelse(.fn == "ab_class", paste0("\"", ab_class, "\""), ""), ") == \"R\"))\n"),
                         ""),
                  "In base R:\n",
                  "  - ", .x_name, "[", scope, "(", deprecated_fn, "), ]\n",
                  ifelse(length(result) > 1, 
                         paste0("  - ", .x_name, "[", scope, "(", 
                                .fn, "(", ifelse(.fn == "ab_class", paste0("\"", ab_class, "\""), ""), ") == \"R\"), ]\n"),
                         ""),
                  "  - subset(", .x_name, ", ", scope, "(", deprecated_fn, "))",
                  ifelse(length(result) > 1, 
                         paste0("\n  - subset(", .x_name, ", ", scope, "(", 
                                .fn, "(", ifelse(.fn == "ab_class", paste0("\"", ab_class, "\""), ""), ") == \"R\"))"),
                         "")),
           call = FALSE)
  
  if (.fn == "ab_class") {
    subset(x, scope_fn(fn(ab_class = ab_class), result))
  } else {
    subset(x, scope_fn(fn(), result))
  }
}

#' @name AMR-deprecated
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
                  .x_name = deparse(substitute(x)),
                  ...)
}

#' @name AMR-deprecated
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
                  .x_name = deparse(substitute(x)),
                  ...)
}
#' @name AMR-deprecated
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
                  .x_name = deparse(substitute(x)),
                  ...)
}

#' @name AMR-deprecated
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
                  .x_name = deparse(substitute(x)),
                  ...)
}

#' @name AMR-deprecated
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
                  .x_name = deparse(substitute(x)),
                  ...)
}

#' @name AMR-deprecated
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
                  .x_name = deparse(substitute(x)),
                  ...)
}

#' @name AMR-deprecated
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
                  .x_name = deparse(substitute(x)),
                  ...)
}

#' @name AMR-deprecated
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
                  .x_name = deparse(substitute(x)),
                  ...)
}

#' @name AMR-deprecated
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
                  .x_name = deparse(substitute(x)),
                  ...)
}

#' @name AMR-deprecated
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
                  .x_name = deparse(substitute(x)),
                  ...)
}

#' @name AMR-deprecated
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
                  .x_name = deparse(substitute(x)),
                  ...)
}

#' @name AMR-deprecated
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
                  .x_name = deparse(substitute(x)),
                  ...)
}

#' @name AMR-deprecated
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
                  .x_name = deparse(substitute(x)),
                  ...)
}

#' @name AMR-deprecated
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
                  .x_name = deparse(substitute(x)),
                  ...)
}

#' @name AMR-deprecated
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
                  .x_name = deparse(substitute(x)),
                  ...)
}
