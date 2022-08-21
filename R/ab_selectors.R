# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Data Analysis for R                   #
#                                                                      #
# SOURCE                                                               #
# https://github.com/msberends/AMR                                     #
#                                                                      #
# LICENCE                                                              #
# (c) 2018-2022 Berends MS, Luz CF et al.                              #
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

#' Antibiotic Selectors
#' 
#' These functions allow for filtering rows and selecting columns based on antibiotic test results that are of a specific antibiotic class or group, without the need to define the columns or antibiotic abbreviations. In short, if you have a column name that resembles an antimicrobial agent, it will be picked up by any of these functions that matches its pharmaceutical class: "cefazolin", "CZO" and "J01DB04" will all be picked up by [cephalosporins()].
#' @param ab_class an antimicrobial class or a part of it, such as `"carba"` and `"carbapenems"`. The columns `group`, `atc_group1` and `atc_group2` of the [antibiotics] data set will be searched (case-insensitive) for this value.
#' @param filter an [expression] to be evaluated in the [antibiotics] data set, such as `name %like% "trim"`
#' @param only_rsi_columns a [logical] to indicate whether only columns of class `<rsi>` must be selected (defaults to `FALSE`), see [as.rsi()]
#' @param only_treatable a [logical] to indicate whether agents that are only for laboratory tests should be excluded (defaults to `TRUE`), such as gentamicin-high (`GEH`) and imipenem/EDTA (`IPE`)
#' @param ... ignored, only in place to allow future extensions
#' @details
#' These functions can be used in data set calls for selecting columns and filtering rows. They are heavily inspired by the [Tidyverse selection helpers][tidyselect::language] such as [`everything()`][tidyselect::everything()], but also work in base \R and not only in `dplyr` verbs. Nonetheless, they are very convenient to use with `dplyr` functions such as [`select()`][dplyr::select()], [`filter()`][dplyr::filter()] and [`summarise()`][dplyr::summarise()], see *Examples*.
#' 
#' All columns in the data in which these functions are called will be searched for known antibiotic names, abbreviations, brand names, and codes (ATC, EARS-Net, WHO, etc.) according to the [antibiotics] data set. This means that a selector such as [aminoglycosides()] will pick up column names like 'gen', 'genta', 'J01GB03', 'tobra', 'Tobracin', etc. 
#' 
#' The [ab_class()] function can be used to filter/select on a manually defined antibiotic class. It searches for results in the [antibiotics] data set within the columns `group`, `atc_group1` and `atc_group2`.
#' @section Full list of supported (antibiotic) classes:
#' 
#' `r paste0(" * ", na.omit(sapply(DEFINED_AB_GROUPS, function(ab) ifelse(tolower(gsub("^AB_", "", ab)) %in% ls(envir = asNamespace("AMR")), paste0("[", tolower(gsub("^AB_", "", ab)), "()] can select: \\cr ", vector_and(paste0(ab_name(eval(parse(text = ab), envir = asNamespace("AMR")), language = NULL, tolower = TRUE), " (", eval(parse(text = ab), envir = asNamespace("AMR")), ")"), quotes = FALSE, sort = TRUE)), character(0)), USE.NAMES = FALSE)), "\n", collapse = "")`
#' @rdname antibiotic_class_selectors
#' @name antibiotic_class_selectors
#' @return (internally) a [character] vector of column names, with additional class `"ab_selector"`
#' @export
#' @inheritSection AMR Reference Data Publicly Available

#' @examples 
#' # `example_isolates` is a data set available in the AMR package.
#' # See ?example_isolates.
#' df <- example_isolates[ , c("hospital_id", "mo",
#'                             "AMP", "AMC", "TZP", "CXM", "CRO", "GEN",
#'                             "TOB", "COL", "IPM", "MEM", "TEC", "VAN")]
#' 
#' # base R ------------------------------------------------------------------
#' 
#' # select columns 'IPM' (imipenem) and 'MEM' (meropenem)
#' df[, carbapenems()]
#' 
#' # select columns 'mo', 'AMK', 'GEN', 'KAN' and 'TOB'
#' df[, c("mo", aminoglycosides())]
#' 
#' # select only antibiotic columns with DDDs for oral treatment
#' df[, administrable_per_os()]
#' 
#' # filter using any() or all()
#' df[any(carbapenems() == "R"), ]
#' subset(df, any(carbapenems() == "R"))
#' 
#' # filter on any or all results in the carbapenem columns (i.e., IPM, MEM):
#' df[any(carbapenems()), ]
#' df[all(carbapenems()), ]
#' 
#' # filter with multiple antibiotic selectors using c()
#' df[all(c(carbapenems(), aminoglycosides()) == "R"), ]
#' 
#' # filter + select in one go: get penicillins in carbapenems-resistant strains
#' df[any(carbapenems() == "R"), penicillins()]
#' 
#' # You can combine selectors with '&' to be more specific. For example,
#' # penicillins() would select benzylpenicillin ('peni G') and
#' # administrable_per_os() would select erythromycin. Yet, when combined these
#' # drugs are both omitted since benzylpenicillin is not administrable per os
#' # and erythromycin is not a penicillin:
#' df[, penicillins() & administrable_per_os()]
#' 
#' # ab_selector() applies a filter in the `antibiotics` data set and is thus very
#' # flexible. For instance, to select antibiotic columns with an oral DDD of at
#' # least 1 gram:
#' df[, ab_selector(oral_ddd > 1 & oral_units == "g")]
#' 
#' # dplyr -------------------------------------------------------------------
#' \donttest{
#' if (require("dplyr")) {
#' 
#'   # get AMR for all aminoglycosides e.g., per hospital:
#'   df %>%
#'     group_by(hospital_id) %>% 
#'     summarise(across(aminoglycosides(), resistance))
#'     
#'   # You can combine selectors with '&' to be more specific:
#'   df %>%
#'     select(penicillins() & administrable_per_os())
#'   
#'   # get AMR for only drugs that matter - no intrinsic resistance:
#'   df %>%
#'     filter(mo_genus() %in% c("Escherichia", "Klebsiella")) %>% 
#'     group_by(hospital_id) %>% 
#'     summarise(across(not_intrinsic_resistant(), resistance))
#'     
#'   # get susceptibility for antibiotics whose name contains "trim":
#'   df %>%
#'     filter(first_isolate()) %>% 
#'     group_by(hospital_id) %>% 
#'     summarise(across(ab_selector(name %like% "trim"), susceptibility))
#' 
#'   # this will select columns 'IPM' (imipenem) and 'MEM' (meropenem):
#'   df %>% 
#'     select(carbapenems())
#'     
#'   # this will select columns 'mo', 'AMK', 'GEN', 'KAN' and 'TOB':
#'   df %>% 
#'     select(mo, aminoglycosides())
#'     
#'  # any() and all() work in dplyr's filter() too:
#'  df %>% 
#'     filter(any(aminoglycosides() == "R"),
#'            all(cephalosporins_2nd() == "R"))
#'     
#'  # also works with c():
#'  df %>% 
#'     filter(any(c(carbapenems(), aminoglycosides()) == "R"))
#'     
#'  # not setting any/all will automatically apply all():
#'  df %>% 
#'     filter(aminoglycosides() == "R")
#'     
#'   # this will select columns 'mo' and all antimycobacterial drugs ('RIF'):
#'   df %>% 
#'     select(mo, ab_class("mycobact"))
#'     
#'   # get bug/drug combinations for only glycopeptides in Gram-positives:
#'   df %>% 
#'     filter(mo_is_gram_positive()) %>% 
#'     select(mo, glycopeptides()) %>% 
#'     bug_drug_combinations() %>%
#'     format()
#'     
#'   data.frame(some_column = "some_value",
#'              J01CA01 = "S") %>%   # ATC code of ampicillin
#'     select(penicillins())         # only the 'J01CA01' column will be selected
#'     
#'     
#'   # with recent versions of dplyr this is all equal:
#'   x <- df[carbapenems() == "R", ]
#'   y <- df %>% filter(carbapenems() == "R")
#'   z <- df %>% filter(if_all(carbapenems(), ~.x == "R"))
#'   identical(x, y)
#'   identical(y, z)
#' }
#' }
ab_class <- function(ab_class, 
                     only_rsi_columns = FALSE,
                     only_treatable = TRUE,
                     ...) {
  meet_criteria(ab_class, allow_class = "character", has_length = 1, allow_NULL = TRUE)
  meet_criteria(only_rsi_columns, allow_class = "logical", has_length = 1)
  meet_criteria(only_treatable, allow_class = "logical", has_length = 1)
  ab_select_exec(NULL, only_rsi_columns = only_rsi_columns, ab_class_args = ab_class, only_treatable = only_treatable)
}

#' @rdname antibiotic_class_selectors
#' @details The [ab_selector()] function can be used to internally filter the [antibiotics] data set on any results, see *Examples*. It allows for filtering on a (part of) a certain name, and/or a group name or even a minimum of DDDs for oral treatment. This function yields the highest flexibility, but is also the least user-friendly, since it requires a hard-coded filter to set.
#' @export
ab_selector <- function(filter, 
                        only_rsi_columns = FALSE,
                        only_treatable = TRUE,
                        ...) {
  meet_criteria(only_rsi_columns, allow_class = "logical", has_length = 1)
  meet_criteria(only_treatable, allow_class = "logical", has_length = 1)
  
  # get_current_data() has to run each time, for cases where e.g., filter() and select() are used in same call
  # but it only takes a couple of milliseconds
  vars_df <- get_current_data(arg_name = NA, call = -2)
  # to improve speed, get_column_abx() will only run once when e.g. in a select or group call
  ab_in_data <- get_column_abx(vars_df, info = FALSE, only_rsi_columns = only_rsi_columns,
                               sort = FALSE, fn = "ab_selector")
  call <- substitute(filter)
  agents <- tryCatch(AMR::antibiotics[which(eval(call, envir = AMR::antibiotics)), "ab", drop = TRUE],
                     error = function(e) stop_(e$message, call = -5))
  agents <- ab_in_data[ab_in_data %in% agents]
  message_agent_names(function_name = "ab_selector", 
                      agents = agents,
                      ab_group = NULL,
                      examples = "",
                      call = call)
  structure(unname(agents),
            class = c("ab_selector", "character"))
}

#' @rdname antibiotic_class_selectors
#' @export
aminoglycosides <- function(only_rsi_columns = FALSE, only_treatable = TRUE, ...) {
  meet_criteria(only_rsi_columns, allow_class = "logical", has_length = 1)
  meet_criteria(only_treatable, allow_class = "logical", has_length = 1)
  ab_select_exec("aminoglycosides", only_rsi_columns = only_rsi_columns, only_treatable = only_treatable)
}

#' @rdname antibiotic_class_selectors
#' @export
aminopenicillins <- function(only_rsi_columns = FALSE, ...) {
  meet_criteria(only_rsi_columns, allow_class = "logical", has_length = 1)
  ab_select_exec("aminopenicillins", only_rsi_columns = only_rsi_columns)
}

#' @rdname antibiotic_class_selectors
#' @export
antifungals <- function(only_rsi_columns = FALSE, ...) {
  meet_criteria(only_rsi_columns, allow_class = "logical", has_length = 1)
  ab_select_exec("antifungals", only_rsi_columns = only_rsi_columns)
}

#' @rdname antibiotic_class_selectors
#' @export
antimycobacterials <- function(only_rsi_columns = FALSE, ...) {
  meet_criteria(only_rsi_columns, allow_class = "logical", has_length = 1)
  ab_select_exec("antimycobacterials", only_rsi_columns = only_rsi_columns)
}

#' @rdname antibiotic_class_selectors
#' @export
betalactams <- function(only_rsi_columns = FALSE, only_treatable = TRUE, ...) {
  meet_criteria(only_rsi_columns, allow_class = "logical", has_length = 1)
  meet_criteria(only_treatable, allow_class = "logical", has_length = 1)
  ab_select_exec("betalactams", only_rsi_columns = only_rsi_columns, only_treatable = only_treatable)
}

#' @rdname antibiotic_class_selectors
#' @export
carbapenems <- function(only_rsi_columns = FALSE, only_treatable = TRUE, ...) {
  meet_criteria(only_rsi_columns, allow_class = "logical", has_length = 1)
  meet_criteria(only_treatable, allow_class = "logical", has_length = 1)
  ab_select_exec("carbapenems", only_rsi_columns = only_rsi_columns, only_treatable = only_treatable)
}

#' @rdname antibiotic_class_selectors
#' @export
cephalosporins <- function(only_rsi_columns = FALSE, ...) {
  meet_criteria(only_rsi_columns, allow_class = "logical", has_length = 1)
  ab_select_exec("cephalosporins", only_rsi_columns = only_rsi_columns)
}

#' @rdname antibiotic_class_selectors
#' @export
cephalosporins_1st <- function(only_rsi_columns = FALSE, ...) {
  meet_criteria(only_rsi_columns, allow_class = "logical", has_length = 1)
  ab_select_exec("cephalosporins_1st", only_rsi_columns = only_rsi_columns)
}

#' @rdname antibiotic_class_selectors
#' @export
cephalosporins_2nd <- function(only_rsi_columns = FALSE, ...) {
  meet_criteria(only_rsi_columns, allow_class = "logical", has_length = 1)
  ab_select_exec("cephalosporins_2nd", only_rsi_columns = only_rsi_columns)
}

#' @rdname antibiotic_class_selectors
#' @export
cephalosporins_3rd <- function(only_rsi_columns = FALSE, ...) {
  meet_criteria(only_rsi_columns, allow_class = "logical", has_length = 1)
  ab_select_exec("cephalosporins_3rd", only_rsi_columns = only_rsi_columns)
}

#' @rdname antibiotic_class_selectors
#' @export
cephalosporins_4th <- function(only_rsi_columns = FALSE, ...) {
  meet_criteria(only_rsi_columns, allow_class = "logical", has_length = 1)
  ab_select_exec("cephalosporins_4th", only_rsi_columns = only_rsi_columns)
}

#' @rdname antibiotic_class_selectors
#' @export
cephalosporins_5th <- function(only_rsi_columns = FALSE, ...) {
  meet_criteria(only_rsi_columns, allow_class = "logical", has_length = 1)
  ab_select_exec("cephalosporins_5th", only_rsi_columns = only_rsi_columns)
}

#' @rdname antibiotic_class_selectors
#' @export
fluoroquinolones <- function(only_rsi_columns = FALSE, ...) {
  meet_criteria(only_rsi_columns, allow_class = "logical", has_length = 1)
  ab_select_exec("fluoroquinolones", only_rsi_columns = only_rsi_columns)
}

#' @rdname antibiotic_class_selectors
#' @export
glycopeptides <- function(only_rsi_columns = FALSE, ...) {
  meet_criteria(only_rsi_columns, allow_class = "logical", has_length = 1)
  ab_select_exec("glycopeptides", only_rsi_columns = only_rsi_columns)
}

#' @rdname antibiotic_class_selectors
#' @export
lincosamides <- function(only_rsi_columns = FALSE, ...) {
  meet_criteria(only_rsi_columns, allow_class = "logical", has_length = 1)
  ab_select_exec("lincosamides", only_rsi_columns = only_rsi_columns)
}

#' @rdname antibiotic_class_selectors
#' @export
lipoglycopeptides <- function(only_rsi_columns = FALSE, ...) {
  meet_criteria(only_rsi_columns, allow_class = "logical", has_length = 1)
  ab_select_exec("lipoglycopeptides", only_rsi_columns = only_rsi_columns)
}

#' @rdname antibiotic_class_selectors
#' @export
macrolides <- function(only_rsi_columns = FALSE, ...) {
  meet_criteria(only_rsi_columns, allow_class = "logical", has_length = 1)
  ab_select_exec("macrolides", only_rsi_columns = only_rsi_columns)
}

#' @rdname antibiotic_class_selectors
#' @export
oxazolidinones <- function(only_rsi_columns = FALSE, ...) {
  meet_criteria(only_rsi_columns, allow_class = "logical", has_length = 1)
  ab_select_exec("oxazolidinones", only_rsi_columns = only_rsi_columns)
}

#' @rdname antibiotic_class_selectors
#' @export
penicillins <- function(only_rsi_columns = FALSE, ...) {
  meet_criteria(only_rsi_columns, allow_class = "logical", has_length = 1)
  ab_select_exec("penicillins", only_rsi_columns = only_rsi_columns)
}

#' @rdname antibiotic_class_selectors
#' @export
polymyxins <- function(only_rsi_columns = FALSE, only_treatable = TRUE, ...) {
  meet_criteria(only_rsi_columns, allow_class = "logical", has_length = 1)
  meet_criteria(only_treatable, allow_class = "logical", has_length = 1)
  ab_select_exec("polymyxins", only_rsi_columns = only_rsi_columns, only_treatable = only_treatable)
}

#' @rdname antibiotic_class_selectors
#' @export
streptogramins <- function(only_rsi_columns = FALSE, ...) {
  meet_criteria(only_rsi_columns, allow_class = "logical", has_length = 1)
  ab_select_exec("streptogramins", only_rsi_columns = only_rsi_columns)
}

#' @rdname antibiotic_class_selectors
#' @export
quinolones <- function(only_rsi_columns = FALSE, ...) {
  meet_criteria(only_rsi_columns, allow_class = "logical", has_length = 1)
  ab_select_exec("quinolones", only_rsi_columns = only_rsi_columns)
}

#' @rdname antibiotic_class_selectors
#' @export
tetracyclines <- function(only_rsi_columns = FALSE, ...) {
  meet_criteria(only_rsi_columns, allow_class = "logical", has_length = 1)
  ab_select_exec("tetracyclines", only_rsi_columns = only_rsi_columns)
}

#' @rdname antibiotic_class_selectors
#' @export
trimethoprims <- function(only_rsi_columns = FALSE, ...) {
  meet_criteria(only_rsi_columns, allow_class = "logical", has_length = 1)
  ab_select_exec("trimethoprims", only_rsi_columns = only_rsi_columns)
}

#' @rdname antibiotic_class_selectors
#' @export
ureidopenicillins <- function(only_rsi_columns = FALSE, ...) {
  meet_criteria(only_rsi_columns, allow_class = "logical", has_length = 1)
  ab_select_exec("ureidopenicillins", only_rsi_columns = only_rsi_columns)
}

#' @rdname antibiotic_class_selectors
#' @details The [administrable_per_os()] and [administrable_iv()] functions also rely on the [antibiotics] data set - antibiotic columns will be matched where a DDD (defined daily dose) for resp. oral and IV treatment is available in the [antibiotics] data set.
#' @export
administrable_per_os <- function(only_rsi_columns = FALSE, ...) {
  meet_criteria(only_rsi_columns, allow_class = "logical", has_length = 1)
  # get_current_data() has to run each time, for cases where e.g., filter() and select() are used in same call
  # but it only takes a couple of milliseconds
  vars_df <- get_current_data(arg_name = NA, call = -2)
  # to improve speed, get_column_abx() will only run once when e.g. in a select or group call
  ab_in_data <- get_column_abx(vars_df, info = FALSE, only_rsi_columns = only_rsi_columns,
                               sort = FALSE, fn = "administrable_per_os")
  agents_all <- antibiotics[which(!is.na(antibiotics$oral_ddd)), "ab", drop = TRUE]
  agents <- antibiotics[which(antibiotics$ab %in% ab_in_data & !is.na(antibiotics$oral_ddd)), "ab", drop = TRUE]
  agents <- ab_in_data[ab_in_data %in% agents]
  message_agent_names(function_name = "administrable_per_os", 
                      agents = agents,
                      ab_group = "administrable_per_os",
                      examples = paste0(" (such as ", 
                                        vector_or(ab_name(sample(agents_all,
                                                                 size = min(5, length(agents_all)),
                                                                 replace = FALSE),
                                                          tolower = TRUE,
                                                          language = NULL),
                                                  quotes = FALSE), 
                                        ")"))
  structure(unname(agents),
            class = c("ab_selector", "character"))
}

#' @rdname antibiotic_class_selectors
#' @export
administrable_iv <- function(only_rsi_columns = FALSE, ...) {
  meet_criteria(only_rsi_columns, allow_class = "logical", has_length = 1)
  # get_current_data() has to run each time, for cases where e.g., filter() and select() are used in same call
  # but it only takes a couple of milliseconds
  vars_df <- get_current_data(arg_name = NA, call = -2)
  # to improve speed, get_column_abx() will only run once when e.g. in a select or group call
  ab_in_data <- get_column_abx(vars_df, info = FALSE, only_rsi_columns = only_rsi_columns,
                               sort = FALSE, fn = "administrable_iv")
  agents_all <- antibiotics[which(!is.na(antibiotics$iv_ddd)), "ab", drop = TRUE]
  agents <- antibiotics[which(antibiotics$ab %in% ab_in_data & !is.na(antibiotics$iv_ddd)), "ab", drop = TRUE]
  agents <- ab_in_data[ab_in_data %in% agents]
  message_agent_names(function_name = "administrable_iv", 
                      agents = agents,
                      ab_group = "administrable_iv",
                      examples = "")
  structure(unname(agents),
            class = c("ab_selector", "character"))
}

#' @rdname antibiotic_class_selectors
#' @inheritParams eucast_rules
#' @details The [not_intrinsic_resistant()] function can be used to only select antibiotic columns that pose no intrinsic resistance for the microorganisms in the data set. For example, if a data set contains only microorganism codes or names of *E. coli* and *K. pneumoniae* and contains a column "vancomycin", this column will be removed (or rather, unselected) using this function. It currently applies `r format_eucast_version_nr(names(EUCAST_VERSION_EXPERT_RULES[length(EUCAST_VERSION_EXPERT_RULES)]))` to determine intrinsic resistance, using the [eucast_rules()] function internally. Because of this determination, this function is quite slow in terms of performance.
#' @export
not_intrinsic_resistant <- function(only_rsi_columns = FALSE, col_mo = NULL, version_expertrules = 3.3, ...) {
  meet_criteria(only_rsi_columns, allow_class = "logical", has_length = 1)
  # get_current_data() has to run each time, for cases where e.g., filter() and select() are used in same call
  # but it only takes a couple of milliseconds
  vars_df <- get_current_data(arg_name = NA, call = -2)
  # to improve speed, get_column_abx() will only run once when e.g. in a select or group call
  ab_in_data <- get_column_abx(vars_df, info = FALSE, only_rsi_columns = only_rsi_columns,
                               sort = FALSE, fn = "not_intrinsic_resistant")
  # intrinsic vars
  vars_df_R <- tryCatch(sapply(eucast_rules(vars_df,
                                            col_mo = col_mo,
                                            version_expertrules = version_expertrules,
                                            rules = "expert",
                                            info = FALSE),
                               function(col) tryCatch(!any(is.na(col)) && all(col == "R"),
                                                      error = function(e) FALSE)),
                        error = function(e) stop_("in not_intrinsic_resistant(): ", e$message, call = FALSE))
  
  agents <- ab_in_data[ab_in_data %in% names(vars_df_R[which(vars_df_R)])]
  if (length(agents) > 0 &&
      message_not_thrown_before("not_intrinsic_resistant", sort(agents))) {
    agents_formatted <- paste0("'", font_bold(agents, collapse = NULL), "'")
    agents_names <- ab_name(names(agents), tolower = TRUE, language = NULL)
    need_name <- generalise_antibiotic_name(agents) != generalise_antibiotic_name(agents_names)
    agents_formatted[need_name] <- paste0(agents_formatted[need_name], " (", agents_names[need_name], ")")
    message_("For `not_intrinsic_resistant()` removing ",
             ifelse(length(agents) == 1, "column ", "columns "),
             vector_and(agents_formatted, quotes = FALSE, sort = FALSE))
  }
  
  vars_df_R <- names(vars_df_R)[which(!vars_df_R)]
  # find columns that are abx, but also intrinsic R
  out <- unname(intersect(ab_in_data, vars_df_R))
  structure(out,
            class = c("ab_selector", "character"))
}

ab_select_exec <- function(function_name,
                           only_rsi_columns = FALSE,
                           only_treatable = FALSE,
                           ab_class_args = NULL) {
  # get_current_data() has to run each time, for cases where e.g., filter() and select() are used in same call
  # but it only takes a couple of milliseconds
  vars_df <- get_current_data(arg_name = NA, call = -3)
  # to improve speed, get_column_abx() will only run once when e.g. in a select or group call
  ab_in_data <- get_column_abx(vars_df, info = FALSE, only_rsi_columns = only_rsi_columns,
                               sort = FALSE, fn = function_name)
  # untreatable drugs
  if (only_treatable == TRUE) {
    untreatable <- antibiotics[which(antibiotics$name %like% "-high|EDTA|polysorbate|macromethod|screening|/nacubactam"), "ab", drop = TRUE]
    if (any(untreatable %in% names(ab_in_data))) {
      if (message_not_thrown_before(function_name, "ab_class", "untreatable", entire_session = TRUE)) {
        warning_("in `", function_name, "()`: some agents were ignored since they cannot be used for treating patients: ",
                 vector_and(ab_name(names(ab_in_data)[names(ab_in_data) %in% untreatable],
                                    language = NULL,
                                    tolower = TRUE),
                            quotes = FALSE,
                            sort = TRUE), ". They can be included using `", function_name, "(only_treatable = FALSE)`. ",
                 "This warning will be shown once per session.")
      }
      ab_in_data <- ab_in_data[!names(ab_in_data) %in% untreatable]
    }
  }
  
  if (length(ab_in_data) == 0) {
    message_("No antimicrobial agents found in the data.")
    return(NULL)
  }
  
  if (is.null(ab_class_args)) {
    # their upper case equivalent are vectors with class <ab>, created in data-raw/_internals.R
    # carbapenems() gets its codes from AMR:::AB_CARBAPENEMS
    abx <- get(paste0("AB_", toupper(function_name)), envir = asNamespace("AMR"))  
    ab_group <- function_name
    examples <- paste0(" (such as ", vector_or(ab_name(sample(abx, size = min(2, length(abx)), replace = FALSE),
                                                       tolower = TRUE,
                                                       language = NULL),
                                               quotes = FALSE), ")")
  } else {
    # this for the 'manual' ab_class() function
    abx <- subset(AB_lookup,
                  group %like% ab_class_args |
                    atc_group1 %like% ab_class_args |
                    atc_group2 %like% ab_class_args)$ab
    ab_group <- find_ab_group(ab_class_args)
    function_name <- "ab_class"
    examples <- paste0(" (such as ", find_ab_names(ab_class_args, 2), ")")
  }
  
  # get the columns with a group names in the chosen ab class
  agents <- ab_in_data[names(ab_in_data) %in% abx]
  
  message_agent_names(function_name = function_name, 
                      agents = agents,
                      ab_group = ab_group,
                      examples = examples,
                      ab_class_args = ab_class_args)
  
  structure(unname(agents),
            class = c("ab_selector", "character"))
}

#' @method c ab_selector
#' @export
#' @noRd
c.ab_selector <- function(...) {
  structure(unlist(lapply(list(...), as.character)),
            class = c("ab_selector", "character"))
}

all_any_ab_selector <- function(type, ..., na.rm = TRUE) {
  cols_ab <- c(...)
  result <- cols_ab[toupper(cols_ab) %in% c("R", "S", "I")]
  if (length(result) == 0) {
    message_("Filtering ", type, " of columns ", vector_and(font_bold(cols_ab, collapse = NULL), quotes = "'"), ' to contain value "R", "S" or "I"')
    result <- c("R", "S", "I")
  }
  cols_ab <- cols_ab[!cols_ab %in% result]
  df <- get_current_data(arg_name = NA, call = -3)
  
  if (type == "all") {
    scope_fn <- all
  } else {
    scope_fn <- any
  }
  
  x_transposed <- as.list(as.data.frame(t(df[, cols_ab, drop = FALSE]), stringsAsFactors = FALSE))
  vapply(FUN.VALUE = logical(1),
         X = x_transposed,
         FUN = function(y) scope_fn(y %in% result, na.rm = na.rm),
         USE.NAMES = FALSE)
}

#' @method all ab_selector
#' @export
#' @noRd
all.ab_selector <- function(..., na.rm = FALSE) {
  all_any_ab_selector("all", ..., na.rm = na.rm)
}

#' @method any ab_selector
#' @export
#' @noRd
any.ab_selector <- function(..., na.rm = FALSE) {
  all_any_ab_selector("any", ..., na.rm = na.rm)
}


#' @method all ab_selector_any_all
#' @export
#' @noRd
all.ab_selector_any_all <- function(..., na.rm = FALSE) {
  # this is all() on a logical vector from `==.ab_selector` or `!=.ab_selector`
  # e.g., example_isolates %>% filter(all(carbapenems() == "R"))
  # so just return the vector as is, only correcting for na.rm
  out <- unclass(c(...))
  if (na.rm == TRUE) {
    out <- out[!is.na(out)]
  }
  out
}

#' @method any ab_selector_any_all
#' @export
#' @noRd
any.ab_selector_any_all <- function(..., na.rm = FALSE) {
  # this is any() on a logical vector from `==.ab_selector` or `!=.ab_selector`
  # e.g., example_isolates %>% filter(any(carbapenems() == "R"))
  # so just return the vector as is, only correcting for na.rm
  out <- unclass(c(...))
  if (na.rm == TRUE) {
    out <- out[!is.na(out)]
  }
  out
}

#' @method == ab_selector
#' @export
#' @noRd
`==.ab_selector` <- function(e1, e2) {
  calls <- as.character(match.call())
  fn_name <- calls[2]
  fn_name <- gsub("^(c\\()(.*)(\\))$", "\\2", fn_name)
  if (is_any(fn_name)) {
    type <- "any"
  } else if (is_all(fn_name)) {
    type <- "all"
  } else {
    type <- "all"
    if (length(e1) > 1) {
      message_("Assuming a filter on ", type, " ", length(e1), " ", gsub("[\\(\\)]", "", fn_name),
               ". Wrap around `all()` or `any()` to prevent this note.")
    }
  }
  structure(all_any_ab_selector(type = type, e1, e2),
            class = c("ab_selector_any_all", "logical"))
}

#' @method != ab_selector
#' @export
#' @noRd
`!=.ab_selector` <- function(e1, e2) {
  calls <- as.character(match.call())
  fn_name <- calls[2]
  fn_name <- gsub("^(c\\()(.*)(\\))$", "\\2", fn_name)
  if (is_any(fn_name)) {
    type <- "any"
  } else if (is_all(fn_name)) {
    type <- "all"
  } else {
    type <- "all"
    if (length(e1) > 1) {
      message_("Assuming a filter on ", type, " ", length(e1), " ", gsub("[\\(\\)]", "", fn_name),
               ". Wrap around `all()` or `any()` to prevent this note.")
    }
  }
  # this is `!=`, so turn around the values
  rsi <- c("R", "S", "I")
  e2 <- rsi[rsi != e2]
  structure(all_any_ab_selector(type = type, e1, e2),
            class = c("ab_selector_any_all", "logical"))
}

#' @method & ab_selector
#' @export
#' @noRd
`&.ab_selector` <- function(e1, e2) {
  # this is only required for base R, since tidyselect has already implemented this
  # e.g., for: example_isolates[, penicillins() & administrable_per_os()]
  structure(intersect(unclass(e1), unclass(e2)),
            class = c("ab_selector", "character"))
}
#' @method | ab_selector
#' @export
#' @noRd
`|.ab_selector` <- function(e1, e2) {
  # this is only required for base R, since tidyselect has already implemented this
  # e.g., for: example_isolates[, penicillins() | administrable_per_os()]
  structure(union(unclass(e1), unclass(e2)),
            class = c("ab_selector", "character"))
}

is_any <- function(el1) {
  syscalls <- paste0(trimws(deparse(sys.calls())), collapse = " ")
  el1 <- gsub("(.*),.*", "\\1", el1)
  syscalls %like% paste0("[^_a-zA-Z0-9]any\\(", "(c\\()?", el1)
}
is_all <- function(el1) {
  syscalls <- paste0(trimws(deparse(sys.calls())), collapse = " ")
  el1 <- gsub("(.*),.*", "\\1", el1)
  syscalls %like% paste0("[^_a-zA-Z0-9]all\\(", "(c\\()?", el1)
}

find_ab_group <- function(ab_class_args) {
  ab_class_args <- gsub("[^a-zA-Z0-9]", ".*", ab_class_args)
  AB_lookup %pm>%
    subset(group %like% ab_class_args | 
             atc_group1 %like% ab_class_args | 
             atc_group2 %like% ab_class_args) %pm>%
    pm_pull(group) %pm>%
    unique() %pm>%
    tolower() %pm>%
    sort() %pm>% 
    paste(collapse = "/")
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
  if (length(drugs) == 0) {
    return("??")
  }
  vector_or(ab_name(sample(drugs, size = min(n, length(drugs)), replace = FALSE),
                    tolower = TRUE,
                    language = NULL),
            quotes = FALSE)
}

message_agent_names <- function(function_name, agents, ab_group = NULL, examples = "", ab_class_args = NULL, call = NULL) {
  if (message_not_thrown_before(function_name, sort(agents))) {
    if (length(agents) == 0) {
      if (is.null(ab_group)) {
        message_("For `", function_name, "()` no antimicrobial agents found", examples, ".")
      } else if (ab_group == "administrable_per_os") {
        message_("No orally administrable agents found", examples, ".")
      } else if (ab_group == "administrable_iv") {
        message_("No IV administrable agents found", examples, ".")
      } else {
        message_("No antimicrobial agents of class '", ab_group, "' found", examples, ".")
      }
    } else {
      agents_formatted <- paste0("'", font_bold(agents, collapse = NULL), "'")
      agents_names <- ab_name(names(agents), tolower = TRUE, language = NULL)
      need_name <- generalise_antibiotic_name(agents) != generalise_antibiotic_name(agents_names)
      agents_formatted[need_name] <- paste0(agents_formatted[need_name], " (", agents_names[need_name], ")")
      message_("For `", function_name, "(",
               ifelse(function_name == "ab_class", 
                      paste0("\"", ab_class_args, "\""),
                      ifelse(!is.null(call),
                             paste0(deparse(call), collapse = " "),
                             "")),
               ")` using ",
               ifelse(length(agents) == 1, "column ", "columns "),
               vector_and(agents_formatted, quotes = FALSE, sort = FALSE))
    }
  }
}
