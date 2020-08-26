# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Analysis                              #
#                                                                      #
# SOURCE                                                               #
# https://github.com/msberends/AMR                                     #
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
# Visit our website for more info: https://msberends.github.io/AMR.    #
# ==================================================================== #

#' Antibiotic class selectors
#' 
#' Use these selection helpers inside any function that allows [Tidyverse selection helpers](https://tidyselect.r-lib.org/reference/language.html), like `dplyr::select()` or `tidyr::pivot_longer()`. They help to select the columns of antibiotics that are of a specific antibiotic class, without the need to define the columns or antibiotic abbreviations.
#' @inheritParams filter_ab_class 
#' @details All columns will be searched for known antibiotic names, abbreviations, brand names and codes (ATC, EARS-Net, WHO, etc.) in the [antibiotics] data set. This means that a selector like e.g. [aminoglycosides()] will pick up column names like 'gen', 'genta', 'J01GB03', 'tobra', 'Tobracin', etc.
#' 
#' **N.B. These functions only work if the `tidyselect` package is installed**, that comes with the `dplyr` package. An error will be thrown if the `tidyselect` package is not installed, or if the functions are used outside a function that allows Tidyverse selections like `select()` or `pivot_longer()`.
#' @rdname antibiotic_class_selectors
#' @seealso [filter_ab_class()] for the `filter()` equivalent.
#' @name antibiotic_class_selectors
#' @export
#' @inheritSection AMR Reference data publicly available
#' @inheritSection AMR Read more on our website!
#' @examples 
#' \dontrun{
#'   library(dplyr)
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
#'     filter(mo_gramstain(mo) %like% "pos") %>% 
#'     select(mo, macrolides()) %>% 
#'     bug_drug_combinations() %>%
#'     format()
#'     
#'     
#'   data.frame(irrelevant = "value",
#'              J01CA01 = "S") %>%   # ATC code of ampicillin
#'     select(penicillins())         # the 'J01CA01' column will be selected
#'
#' }
ab_class <- function(ab_class) {
  ab_selector(ab_class, function_name = "ab_class")
}

#' @rdname antibiotic_class_selectors
#' @export
aminoglycosides <- function() {
  ab_selector("aminoglycoside", function_name = "aminoglycosides")
}

#' @rdname antibiotic_class_selectors
#' @export
carbapenems <- function() {
  ab_selector("carbapenem", function_name = "carbapenems")
}

#' @rdname antibiotic_class_selectors
#' @export
cephalosporins <- function() {
  ab_selector("cephalosporin", function_name = "cephalosporins")
}

#' @rdname antibiotic_class_selectors
#' @export
cephalosporins_1st <- function() {
  ab_selector("cephalosporins.*1", function_name = "cephalosporins_1st")
}

#' @rdname antibiotic_class_selectors
#' @export
cephalosporins_2nd <- function() {
  ab_selector("cephalosporins.*2", function_name = "cephalosporins_2nd")
}

#' @rdname antibiotic_class_selectors
#' @export
cephalosporins_3rd <- function() {
  ab_selector("cephalosporins.*3", function_name = "cephalosporins_3rd")
}

#' @rdname antibiotic_class_selectors
#' @export
cephalosporins_4th <- function() {
  ab_selector("cephalosporins.*4", function_name = "cephalosporins_4th")
}

#' @rdname antibiotic_class_selectors
#' @export
cephalosporins_5th <- function() {
  ab_selector("cephalosporins.*5", function_name = "cephalosporins_5th")
}

#' @rdname antibiotic_class_selectors
#' @export
fluoroquinolones <- function() {
  ab_selector("fluoroquinolone", function_name = "fluoroquinolones")
}

#' @rdname antibiotic_class_selectors
#' @export
glycopeptides <- function() {
  ab_selector("glycopeptide", function_name = "glycopeptides")
}

#' @rdname antibiotic_class_selectors
#' @export
macrolides <- function() {
  ab_selector("macrolide", function_name = "macrolides")
}

#' @rdname antibiotic_class_selectors
#' @export
penicillins <- function() {
  ab_selector("penicillin", function_name = "penicillins")
}

#' @rdname antibiotic_class_selectors
#' @export
tetracyclines <- function() {
  ab_selector("tetracycline", function_name = "tetracyclines")
}

ab_selector <- function(ab_class, function_name) {
  peek_vars_tidyselect <- import_fn("peek_vars", "tidyselect")
  vars_vct <- peek_vars_tidyselect(fn = function_name)
  vars_df <- data.frame(as.list(vars_vct))[0, , drop = FALSE]
  colnames(vars_df) <- vars_vct
  ab_in_data <- suppressMessages(get_column_abx(vars_df))
  
  if (length(ab_in_data) == 0) {
    message(font_blue("NOTE: no antimicrobial agents found."))
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
  if (length(agents) == 0) {
    message(font_blue(paste0("NOTE: No antimicrobial agents of class ", ab_group, 
                             " found", examples, ".")))
  } else {
    message(font_blue(paste0("Selecting ", ab_group, ": ",
                             paste(paste0("`", font_bold(agents, collapse = NULL),
                                          "` (", ab_name(names(agents), tolower = TRUE, language = NULL), ")"),
                                   collapse = ", "))))
  }
  unname(agents)
}
