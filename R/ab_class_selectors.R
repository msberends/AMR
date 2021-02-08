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

#' Antibiotic Class Selectors
#' 
#' These functions help to select the columns of antibiotics that are of a specific antibiotic class, without the need to define the columns or antibiotic abbreviations. \strong{\Sexpr{ifelse(as.double(R.Version()$major) + (as.double(R.Version()$minor) / 10) < 3.2, paste0("NOTE: THESE FUNCTIONS DO NOT WORK ON YOUR CURRENT R VERSION. These functions require R version 3.2 or later - you have ", R.version.string, "."), "")}}
#' @inheritSection lifecycle Stable Lifecycle
#' @param only_rsi_columns a logical to indicate whether only columns of class [`<rsi>`]([rsi]) must be selected (defaults to `FALSE`)
#' @inheritParams filter_ab_class 
#' @details \strong{\Sexpr{ifelse(as.double(R.Version()$major) + (as.double(R.Version()$minor) / 10) < 3.2, paste0("NOTE: THESE FUNCTIONS DO NOT WORK ON YOUR CURRENT R VERSION. These functions require R version 3.2 or later - you have ", R.version.string, "."), "")}}
#' 
#' All columns will be searched for known antibiotic names, abbreviations, brand names and codes (ATC, EARS-Net, WHO, etc.) in the [antibiotics] data set. This means that a selector like e.g. [aminoglycosides()] will pick up column names like 'gen', 'genta', 'J01GB03', 'tobra', 'Tobracin', etc.
#' @rdname antibiotic_class_selectors
#' @seealso [filter_ab_class()] for the `filter()` equivalent.
#' @name antibiotic_class_selectors
#' @export
#' @inheritSection AMR Reference Data Publicly Available
#' @inheritSection AMR Read more on Our Website!
#' @examples 
#' # `example_isolates` is a data set available in the AMR package.
#' # See ?example_isolates.
#' 
#' # this will select columns 'IPM' (imipenem) and 'MEM' (meropenem):
#' example_isolates[, c(carbapenems())]
#' # this will select columns 'mo', 'AMK', 'GEN', 'KAN' and 'TOB':
#' example_isolates[, c("mo", aminoglycosides())]
#' 
#' if (require("dplyr")) {
#' 
#'   # this will select columns 'IPM' (imipenem) and 'MEM' (meropenem):
#'   example_isolates %>% 
#'     select(carbapenems())
#'     
#'   # this will select columns 'mo', 'AMK', 'GEN', 'KAN' and 'TOB':
#'   example_isolates %>% 
#'     select(mo, aminoglycosides())
#'     
#'   # this will select columns 'mo' and all antimycobacterial drugs ('RIF'):
#'   example_isolates %>% 
#'     select(mo, ab_class("mycobact"))
#'     
#'     
#'   # get bug/drug combinations for only macrolides in Gram-positives:
#'   example_isolates %>% 
#'     filter(mo_is_gram_positive()) %>% 
#'     select(mo, macrolides()) %>% 
#'     bug_drug_combinations() %>%
#'     format()
#'     
#'     
#'   data.frame(some_column = "some_value",
#'              J01CA01 = "S") %>%   # ATC code of ampicillin
#'     select(penicillins())         # only the 'J01CA01' column will be selected
#'     
#'     
#'   # with dplyr 1.0.0 and higher (that adds 'across()'), this is equal:
#'   # (though the row names on the first are more correct)
#'   example_isolates %>% filter_carbapenems("R", "all")
#'   example_isolates %>% filter(across(carbapenems(), ~. == "R"))
#' }
ab_class <- function(ab_class, 
                     only_rsi_columns = FALSE) {
  ab_selector(ab_class, function_name = "ab_class", only_rsi_columns = only_rsi_columns)
}

#' @rdname antibiotic_class_selectors
#' @export
aminoglycosides <- function(only_rsi_columns = FALSE) {
  ab_selector("aminoglycoside", function_name = "aminoglycosides", only_rsi_columns = only_rsi_columns)
}

#' @rdname antibiotic_class_selectors
#' @export
carbapenems <- function(only_rsi_columns = FALSE) {
  ab_selector("carbapenem", function_name = "carbapenems", only_rsi_columns = only_rsi_columns)
}

#' @rdname antibiotic_class_selectors
#' @export
cephalosporins <- function(only_rsi_columns = FALSE) {
  ab_selector("cephalosporin", function_name = "cephalosporins", only_rsi_columns = only_rsi_columns)
}

#' @rdname antibiotic_class_selectors
#' @export
cephalosporins_1st <- function(only_rsi_columns = FALSE) {
  ab_selector("cephalosporins.*1", function_name = "cephalosporins_1st", only_rsi_columns = only_rsi_columns)
}

#' @rdname antibiotic_class_selectors
#' @export
cephalosporins_2nd <- function(only_rsi_columns = FALSE) {
  ab_selector("cephalosporins.*2", function_name = "cephalosporins_2nd", only_rsi_columns = only_rsi_columns)
}

#' @rdname antibiotic_class_selectors
#' @export
cephalosporins_3rd <- function(only_rsi_columns = FALSE) {
  ab_selector("cephalosporins.*3", function_name = "cephalosporins_3rd", only_rsi_columns = only_rsi_columns)
}

#' @rdname antibiotic_class_selectors
#' @export
cephalosporins_4th <- function(only_rsi_columns = FALSE) {
  ab_selector("cephalosporins.*4", function_name = "cephalosporins_4th", only_rsi_columns = only_rsi_columns)
}

#' @rdname antibiotic_class_selectors
#' @export
cephalosporins_5th <- function(only_rsi_columns = FALSE) {
  ab_selector("cephalosporins.*5", function_name = "cephalosporins_5th", only_rsi_columns = only_rsi_columns)
}

#' @rdname antibiotic_class_selectors
#' @export
fluoroquinolones <- function(only_rsi_columns = FALSE) {
  ab_selector("fluoroquinolone", function_name = "fluoroquinolones", only_rsi_columns = only_rsi_columns)
}

#' @rdname antibiotic_class_selectors
#' @export
glycopeptides <- function(only_rsi_columns = FALSE) {
  ab_selector("glycopeptide", function_name = "glycopeptides", only_rsi_columns = only_rsi_columns)
}

#' @rdname antibiotic_class_selectors
#' @export
macrolides <- function(only_rsi_columns = FALSE) {
  ab_selector("macrolide", function_name = "macrolides", only_rsi_columns = only_rsi_columns)
}

#' @rdname antibiotic_class_selectors
#' @export
oxazolidinones <- function(only_rsi_columns = FALSE) {
  ab_selector("oxazolidinone", function_name = "oxazolidinones", only_rsi_columns = only_rsi_columns)
}

#' @rdname antibiotic_class_selectors
#' @export
penicillins <- function(only_rsi_columns = FALSE) {
  ab_selector("penicillin", function_name = "penicillins", only_rsi_columns = only_rsi_columns)
}

#' @rdname antibiotic_class_selectors
#' @export
tetracyclines <- function(only_rsi_columns = FALSE) {
  ab_selector("tetracycline", function_name = "tetracyclines", only_rsi_columns = only_rsi_columns)
}

ab_selector <- function(ab_class,
                        function_name,
                        only_rsi_columns) {
  meet_criteria(ab_class, allow_class = "character", has_length = 1, .call_depth = 1)
  meet_criteria(function_name, allow_class = "character", has_length = 1, .call_depth = 1)
  meet_criteria(only_rsi_columns, allow_class = "logical", has_length = 1, .call_depth = 1)
  
  if (as.double(R.Version()$major) + (as.double(R.Version()$minor) / 10) < 3.2) {
    warning_("antibiotic class selectors such as ", function_name, 
             "() require R version 3.2 or later - you have ", R.version.string,
             call = FALSE)
    return(NULL)
  }
  
  vars_df <- get_current_data(arg_name = NA, call = -3)

  # improve speed here so it will only run once when e.g. in one select call
  if (!identical(pkg_env$ab_selector, unique_call_id())) {
    ab_in_data <- get_column_abx(vars_df, info = FALSE, only_rsi_columns = only_rsi_columns)
    pkg_env$ab_selector <- unique_call_id()
    pkg_env$ab_selector_cols <- ab_in_data
  } else {
    ab_in_data <- pkg_env$ab_selector_cols
  }
  
  if (length(ab_in_data) == 0) {
    message_("No antimicrobial agents found.")
    return(NULL)
  }
  
  ab_reference <- subset(antibiotics,
                         group %like% ab_class | 
                           atc_group1 %like% ab_class | 
                           atc_group2 %like% ab_class)
  ab_group <- find_ab_group(ab_class)
  if (ab_group == "") {
    ab_group <- paste0("'", ab_class, "'")
    examples <- ""
  } else {
    examples <- paste0(" (such as ", find_ab_names(ab_class, 2), ")")
  }
  # get the columns with a group names in the chosen ab class
  agents <- ab_in_data[names(ab_in_data) %in% ab_reference$ab]
  if (message_not_thrown_before(function_name)) {
    if (length(agents) == 0) {
      message_("No antimicrobial agents of class ", ab_group, " found", examples, ".")
    } else {
      agents_formatted <- paste0("'", font_bold(agents, collapse = NULL), "'")
      agents_names <- ab_name(names(agents), tolower = TRUE, language = NULL)
      need_name <- tolower(gsub("[^a-zA-Z]", "", agents)) != tolower(gsub("[^a-zA-Z]", "", agents_names))
      agents_formatted[need_name] <- paste0(agents_formatted[need_name],
                                            " (", agents_names[need_name], ")")
      message_("Selecting ", ab_group, ": ",
               ifelse(length(agents) == 1, "column ", "columns "),
               vector_and(agents_formatted, quotes = FALSE),
               as_note = FALSE,
               extra_indent = 6)
    }
    remember_thrown_message(function_name)
  }
  unname(agents)
}
