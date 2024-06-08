# ==================================================================== #
# TITLE:                                                               #
# AMR: An R Package for Working with Antimicrobial Resistance Data     #
#                                                                      #
# SOURCE CODE:                                                         #
# https://github.com/msberends/AMR                                     #
#                                                                      #
# PLEASE CITE THIS SOFTWARE AS:                                        #
# Berends MS, Luz CF, Friedrich AW, Sinha BNM, Albers CJ, Glasner C    #
# (2022). AMR: An R Package for Working with Antimicrobial Resistance  #
# Data. Journal of Statistical Software, 104(3), 1-31.                 #
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
# how to conduct AMR data analysis: https://msberends.github.io/AMR/   #
# ==================================================================== #

#' Antibiotic Selectors
#'
#' @description These functions allow for filtering rows and selecting columns based on antibiotic test results that are of a specific antibiotic class or group (according to the [antibiotics] data set), without the need to define the columns or antibiotic abbreviations.
#'
#' In short, if you have a column name that resembles an antimicrobial drug, it will be picked up by any of these functions that matches its pharmaceutical class: "cefazolin", "kefzol", "CZO" and "J01DB04" will all be picked up by [cephalosporins()].
#' @param ab_class an antimicrobial class or a part of it, such as `"carba"` and `"carbapenems"`. The columns `group`, `atc_group1` and `atc_group2` of the [antibiotics] data set will be searched (case-insensitive) for this value.
#' @param filter an [expression] to be evaluated in the [antibiotics] data set, such as `name %like% "trim"`
#' @param only_sir_columns a [logical] to indicate whether only columns of class `sir` must be selected (default is `FALSE`), see [as.sir()]
#' @param only_treatable a [logical] to indicate whether antimicrobial drugs should be excluded that are only for laboratory tests (default is `TRUE`), such as gentamicin-high (`GEH`) and imipenem/EDTA (`IPE`)
#' @param ... ignored, only in place to allow future extensions
#' @details
#' These functions can be used in data set calls for selecting columns and filtering rows. They work with base \R, the Tidyverse, and `data.table`. They are heavily inspired by the [Tidyverse selection helpers][tidyselect::language] such as [`everything()`][tidyselect::everything()], but are not limited to `dplyr` verbs. Nonetheless, they are very convenient to use with `dplyr` functions such as [`select()`][dplyr::select()], [`filter()`][dplyr::filter()] and [`summarise()`][dplyr::summarise()], see *Examples*.
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
#' example_isolates
#'
#'
#' # Examples sections below are split into 'base R', 'dplyr', and 'data.table':
#'
#'
#' # base R ------------------------------------------------------------------
#'
#' # select columns 'IPM' (imipenem) and 'MEM' (meropenem)
#' example_isolates[, carbapenems()]
#'
#' # select columns 'mo', 'AMK', 'GEN', 'KAN' and 'TOB'
#' example_isolates[, c("mo", aminoglycosides())]
#'
#' # select only antibiotic columns with DDDs for oral treatment
#' example_isolates[, administrable_per_os()]
#'
#' # filter using any() or all()
#' example_isolates[any(carbapenems() == "R"), ]
#' subset(example_isolates, any(carbapenems() == "R"))
#'
#' # filter on any or all results in the carbapenem columns (i.e., IPM, MEM):
#' example_isolates[any(carbapenems()), ]
#' example_isolates[all(carbapenems()), ]
#'
#' # filter with multiple antibiotic selectors using c()
#' example_isolates[all(c(carbapenems(), aminoglycosides()) == "R"), ]
#'
#' # filter + select in one go: get penicillins in carbapenem-resistant strains
#' example_isolates[any(carbapenems() == "R"), penicillins()]
#'
#' # You can combine selectors with '&' to be more specific. For example,
#' # penicillins() would select benzylpenicillin ('peni G') and
#' # administrable_per_os() would select erythromycin. Yet, when combined these
#' # drugs are both omitted since benzylpenicillin is not administrable per os
#' # and erythromycin is not a penicillin:
#' example_isolates[, penicillins() & administrable_per_os()]
#'
#' # ab_selector() applies a filter in the `antibiotics` data set and is thus
#' # very flexible. For instance, to select antibiotic columns with an oral DDD
#' # of at least 1 gram:
#' example_isolates[, ab_selector(oral_ddd > 1 & oral_units == "g")]
#'
#' \donttest{
#' # dplyr -------------------------------------------------------------------
#'
#' if (require("dplyr")) {
#'   tibble(kefzol = random_sir(5)) %>%
#'     select(cephalosporins())
#' }
#'
#' if (require("dplyr")) {
#'   # get AMR for all aminoglycosides e.g., per ward:
#'   example_isolates %>%
#'     group_by(ward) %>%
#'     summarise(across(aminoglycosides(), resistance))
#' }
#' if (require("dplyr")) {
#'   # You can combine selectors with '&' to be more specific:
#'   example_isolates %>%
#'     select(penicillins() & administrable_per_os())
#' }
#' if (require("dplyr")) {
#'   # get AMR for only drugs that matter - no intrinsic resistance:
#'   example_isolates %>%
#'     filter(mo_genus() %in% c("Escherichia", "Klebsiella")) %>%
#'     group_by(ward) %>%
#'     summarise(across(not_intrinsic_resistant(), resistance))
#' }
#' if (require("dplyr")) {
#'   # get susceptibility for antibiotics whose name contains "trim":
#'   example_isolates %>%
#'     filter(first_isolate()) %>%
#'     group_by(ward) %>%
#'     summarise(across(ab_selector(name %like% "trim"), susceptibility))
#' }
#' if (require("dplyr")) {
#'   # this will select columns 'IPM' (imipenem) and 'MEM' (meropenem):
#'   example_isolates %>%
#'     select(carbapenems())
#' }
#' if (require("dplyr")) {
#'   # this will select columns 'mo', 'AMK', 'GEN', 'KAN' and 'TOB':
#'   example_isolates %>%
#'     select(mo, aminoglycosides())
#' }
#' if (require("dplyr")) {
#'   # any() and all() work in dplyr's filter() too:
#'   example_isolates %>%
#'     filter(
#'       any(aminoglycosides() == "R"),
#'       all(cephalosporins_2nd() == "R")
#'     )
#' }
#' if (require("dplyr")) {
#'   # also works with c():
#'   example_isolates %>%
#'     filter(any(c(carbapenems(), aminoglycosides()) == "R"))
#' }
#' if (require("dplyr")) {
#'   # not setting any/all will automatically apply all():
#'   example_isolates %>%
#'     filter(aminoglycosides() == "R")
#' }
#' if (require("dplyr")) {
#'   # this will select columns 'mo' and all antimycobacterial drugs ('RIF'):
#'   example_isolates %>%
#'     select(mo, ab_class("mycobact"))
#' }
#' if (require("dplyr")) {
#'   # get bug/drug combinations for only glycopeptides in Gram-positives:
#'   example_isolates %>%
#'     filter(mo_is_gram_positive()) %>%
#'     select(mo, glycopeptides()) %>%
#'     bug_drug_combinations() %>%
#'     format()
#' }
#' if (require("dplyr")) {
#'   data.frame(
#'     some_column = "some_value",
#'     J01CA01 = "S"
#'   ) %>% # ATC code of ampicillin
#'     select(penicillins()) # only the 'J01CA01' column will be selected
#' }
#' if (require("dplyr")) {
#'   # with recent versions of dplyr, this is all equal:
#'   x <- example_isolates[carbapenems() == "R", ]
#'   y <- example_isolates %>% filter(carbapenems() == "R")
#'   z <- example_isolates %>% filter(if_all(carbapenems(), ~ .x == "R"))
#'   identical(x, y) && identical(y, z)
#' }
#'
#'
#' # data.table --------------------------------------------------------------
#'
#' # data.table is supported as well, just use it in the same way as with
#' # base R, but add `with = FALSE` if using a single AB selector.
#'
#' if (require("data.table")) {
#'   dt <- as.data.table(example_isolates)
#'
#'   # this does not work, it returns column *names*
#'   dt[, carbapenems()]
#' }
#' if (require("data.table")) {
#'   # so `with = FALSE` is required
#'   dt[, carbapenems(), with = FALSE]
#' }
#'
#' # for multiple selections or AB selectors, `with = FALSE` is not needed:
#' if (require("data.table")) {
#'   dt[, c("mo", aminoglycosides())]
#' }
#' if (require("data.table")) {
#'   dt[, c(carbapenems(), aminoglycosides())]
#' }
#'
#' # row filters are also supported:
#' if (require("data.table")) {
#'   dt[any(carbapenems() == "S"), ]
#' }
#' if (require("data.table")) {
#'   dt[any(carbapenems() == "S"), penicillins(), with = FALSE]
#' }
#' }
ab_class <- function(ab_class,
                     only_sir_columns = FALSE,
                     only_treatable = TRUE,
                     ...) {
  meet_criteria(ab_class, allow_class = "character", has_length = 1, allow_NULL = TRUE)
  meet_criteria(only_sir_columns, allow_class = "logical", has_length = 1)
  meet_criteria(only_treatable, allow_class = "logical", has_length = 1)
  ab_select_exec(NULL, only_sir_columns = only_sir_columns, ab_class_args = ab_class, only_treatable = only_treatable)
}

#' @rdname antibiotic_class_selectors
#' @details The [ab_selector()] function can be used to internally filter the [antibiotics] data set on any results, see *Examples*. It allows for filtering on a (part of) a certain name, and/or a group name or even a minimum of DDDs for oral treatment. This function yields the highest flexibility, but is also the least user-friendly, since it requires a hard-coded filter to set.
#' @export
ab_selector <- function(filter,
                        only_sir_columns = FALSE,
                        only_treatable = TRUE,
                        ...) {
  meet_criteria(only_sir_columns, allow_class = "logical", has_length = 1)
  meet_criteria(only_treatable, allow_class = "logical", has_length = 1)

  # get_current_data() has to run each time, for cases where e.g., filter() and select() are used in same call
  # but it only takes a couple of milliseconds
  vars_df <- get_current_data(arg_name = NA, call = -2)
  # to improve speed, get_column_abx() will only run once when e.g. in a select or group call
  ab_in_data <- get_column_abx(vars_df,
    info = FALSE, only_sir_columns = only_sir_columns,
    sort = FALSE, fn = "ab_selector"
  )
  call <- substitute(filter)
  agents <- tryCatch(AMR_env$AB_lookup[which(eval(call, envir = AMR_env$AB_lookup)), "ab", drop = TRUE],
    error = function(e) stop_(e$message, call = -5)
  )
  agents <- ab_in_data[ab_in_data %in% agents]
  message_agent_names(
    function_name = "ab_selector",
    agents = agents,
    ab_group = NULL,
    examples = "",
    call = call
  )
  structure(unname(agents),
    class = c("ab_selector", "character")
  )
}

#' @rdname antibiotic_class_selectors
#' @export
aminoglycosides <- function(only_sir_columns = FALSE, only_treatable = TRUE, ...) {
  meet_criteria(only_sir_columns, allow_class = "logical", has_length = 1)
  meet_criteria(only_treatable, allow_class = "logical", has_length = 1)
  ab_select_exec("aminoglycosides", only_sir_columns = only_sir_columns, only_treatable = only_treatable)
}

#' @rdname antibiotic_class_selectors
#' @export
aminopenicillins <- function(only_sir_columns = FALSE, ...) {
  meet_criteria(only_sir_columns, allow_class = "logical", has_length = 1)
  ab_select_exec("aminopenicillins", only_sir_columns = only_sir_columns)
}

#' @rdname antibiotic_class_selectors
#' @export
antifungals <- function(only_sir_columns = FALSE, ...) {
  meet_criteria(only_sir_columns, allow_class = "logical", has_length = 1)
  ab_select_exec("antifungals", only_sir_columns = only_sir_columns)
}

#' @rdname antibiotic_class_selectors
#' @export
antimycobacterials <- function(only_sir_columns = FALSE, ...) {
  meet_criteria(only_sir_columns, allow_class = "logical", has_length = 1)
  ab_select_exec("antimycobacterials", only_sir_columns = only_sir_columns)
}

#' @rdname antibiotic_class_selectors
#' @export
betalactams <- function(only_sir_columns = FALSE, only_treatable = TRUE, ...) {
  meet_criteria(only_sir_columns, allow_class = "logical", has_length = 1)
  meet_criteria(only_treatable, allow_class = "logical", has_length = 1)
  ab_select_exec("betalactams", only_sir_columns = only_sir_columns, only_treatable = only_treatable)
}

#' @rdname antibiotic_class_selectors
#' @export
carbapenems <- function(only_sir_columns = FALSE, only_treatable = TRUE, ...) {
  meet_criteria(only_sir_columns, allow_class = "logical", has_length = 1)
  meet_criteria(only_treatable, allow_class = "logical", has_length = 1)
  ab_select_exec("carbapenems", only_sir_columns = only_sir_columns, only_treatable = only_treatable)
}

#' @rdname antibiotic_class_selectors
#' @export
cephalosporins <- function(only_sir_columns = FALSE, ...) {
  meet_criteria(only_sir_columns, allow_class = "logical", has_length = 1)
  ab_select_exec("cephalosporins", only_sir_columns = only_sir_columns)
}

#' @rdname antibiotic_class_selectors
#' @export
cephalosporins_1st <- function(only_sir_columns = FALSE, ...) {
  meet_criteria(only_sir_columns, allow_class = "logical", has_length = 1)
  ab_select_exec("cephalosporins_1st", only_sir_columns = only_sir_columns)
}

#' @rdname antibiotic_class_selectors
#' @export
cephalosporins_2nd <- function(only_sir_columns = FALSE, ...) {
  meet_criteria(only_sir_columns, allow_class = "logical", has_length = 1)
  ab_select_exec("cephalosporins_2nd", only_sir_columns = only_sir_columns)
}

#' @rdname antibiotic_class_selectors
#' @export
cephalosporins_3rd <- function(only_sir_columns = FALSE, ...) {
  meet_criteria(only_sir_columns, allow_class = "logical", has_length = 1)
  ab_select_exec("cephalosporins_3rd", only_sir_columns = only_sir_columns)
}

#' @rdname antibiotic_class_selectors
#' @export
cephalosporins_4th <- function(only_sir_columns = FALSE, ...) {
  meet_criteria(only_sir_columns, allow_class = "logical", has_length = 1)
  ab_select_exec("cephalosporins_4th", only_sir_columns = only_sir_columns)
}

#' @rdname antibiotic_class_selectors
#' @export
cephalosporins_5th <- function(only_sir_columns = FALSE, ...) {
  meet_criteria(only_sir_columns, allow_class = "logical", has_length = 1)
  ab_select_exec("cephalosporins_5th", only_sir_columns = only_sir_columns)
}

#' @rdname antibiotic_class_selectors
#' @export
fluoroquinolones <- function(only_sir_columns = FALSE, ...) {
  meet_criteria(only_sir_columns, allow_class = "logical", has_length = 1)
  ab_select_exec("fluoroquinolones", only_sir_columns = only_sir_columns)
}

#' @rdname antibiotic_class_selectors
#' @export
glycopeptides <- function(only_sir_columns = FALSE, ...) {
  meet_criteria(only_sir_columns, allow_class = "logical", has_length = 1)
  ab_select_exec("glycopeptides", only_sir_columns = only_sir_columns)
}

#' @rdname antibiotic_class_selectors
#' @export
lincosamides <- function(only_sir_columns = FALSE, only_treatable = TRUE, ...) {
  meet_criteria(only_sir_columns, allow_class = "logical", has_length = 1)
  meet_criteria(only_treatable, allow_class = "logical", has_length = 1)
  ab_select_exec("lincosamides", only_sir_columns = only_sir_columns, only_treatable = only_treatable)
}

#' @rdname antibiotic_class_selectors
#' @export
lipoglycopeptides <- function(only_sir_columns = FALSE, ...) {
  meet_criteria(only_sir_columns, allow_class = "logical", has_length = 1)
  ab_select_exec("lipoglycopeptides", only_sir_columns = only_sir_columns)
}

#' @rdname antibiotic_class_selectors
#' @export
macrolides <- function(only_sir_columns = FALSE, ...) {
  meet_criteria(only_sir_columns, allow_class = "logical", has_length = 1)
  ab_select_exec("macrolides", only_sir_columns = only_sir_columns)
}

#' @rdname antibiotic_class_selectors
#' @export
nitrofurans <- function(only_sir_columns = FALSE, ...) {
  meet_criteria(only_sir_columns, allow_class = "logical", has_length = 1)
  ab_select_exec("nitrofurans", only_sir_columns = only_sir_columns)
}

#' @rdname antibiotic_class_selectors
#' @export
oxazolidinones <- function(only_sir_columns = FALSE, ...) {
  meet_criteria(only_sir_columns, allow_class = "logical", has_length = 1)
  ab_select_exec("oxazolidinones", only_sir_columns = only_sir_columns)
}

#' @rdname antibiotic_class_selectors
#' @export
penicillins <- function(only_sir_columns = FALSE, ...) {
  meet_criteria(only_sir_columns, allow_class = "logical", has_length = 1)
  ab_select_exec("penicillins", only_sir_columns = only_sir_columns)
}

#' @rdname antibiotic_class_selectors
#' @export
polymyxins <- function(only_sir_columns = FALSE, only_treatable = TRUE, ...) {
  meet_criteria(only_sir_columns, allow_class = "logical", has_length = 1)
  meet_criteria(only_treatable, allow_class = "logical", has_length = 1)
  ab_select_exec("polymyxins", only_sir_columns = only_sir_columns, only_treatable = only_treatable)
}

#' @rdname antibiotic_class_selectors
#' @export
quinolones <- function(only_sir_columns = FALSE, ...) {
  meet_criteria(only_sir_columns, allow_class = "logical", has_length = 1)
  ab_select_exec("quinolones", only_sir_columns = only_sir_columns)
}

#' @rdname antibiotic_class_selectors
#' @export
rifamycins <- function(only_sir_columns = FALSE, ...) {
  meet_criteria(only_sir_columns, allow_class = "logical", has_length = 1)
  ab_select_exec("rifamycins", only_sir_columns = only_sir_columns)
}

#' @rdname antibiotic_class_selectors
#' @export
streptogramins <- function(only_sir_columns = FALSE, ...) {
  meet_criteria(only_sir_columns, allow_class = "logical", has_length = 1)
  ab_select_exec("streptogramins", only_sir_columns = only_sir_columns)
}

#' @rdname antibiotic_class_selectors
#' @export
tetracyclines <- function(only_sir_columns = FALSE, ...) {
  meet_criteria(only_sir_columns, allow_class = "logical", has_length = 1)
  ab_select_exec("tetracyclines", only_sir_columns = only_sir_columns)
}

#' @rdname antibiotic_class_selectors
#' @export
trimethoprims <- function(only_sir_columns = FALSE, ...) {
  meet_criteria(only_sir_columns, allow_class = "logical", has_length = 1)
  ab_select_exec("trimethoprims", only_sir_columns = only_sir_columns)
}

#' @rdname antibiotic_class_selectors
#' @export
ureidopenicillins <- function(only_sir_columns = FALSE, ...) {
  meet_criteria(only_sir_columns, allow_class = "logical", has_length = 1)
  ab_select_exec("ureidopenicillins", only_sir_columns = only_sir_columns)
}

#' @rdname antibiotic_class_selectors
#' @details The [administrable_per_os()] and [administrable_iv()] functions also rely on the [antibiotics] data set - antibiotic columns will be matched where a DDD (defined daily dose) for resp. oral and IV treatment is available in the [antibiotics] data set.
#' @export
administrable_per_os <- function(only_sir_columns = FALSE, ...) {
  meet_criteria(only_sir_columns, allow_class = "logical", has_length = 1)
  # get_current_data() has to run each time, for cases where e.g., filter() and select() are used in same call
  # but it only takes a couple of milliseconds
  vars_df <- get_current_data(arg_name = NA, call = -2)
  # to improve speed, get_column_abx() will only run once when e.g. in a select or group call
  ab_in_data <- get_column_abx(vars_df,
    info = FALSE, only_sir_columns = only_sir_columns,
    sort = FALSE, fn = "administrable_per_os"
  )
  agents_all <- AMR_env$AB_lookup[which(!is.na(AMR_env$AB_lookup$oral_ddd)), "ab", drop = TRUE]
  agents <- AMR_env$AB_lookup[which(AMR_env$AB_lookup$ab %in% ab_in_data & !is.na(AMR_env$AB_lookup$oral_ddd)), "ab", drop = TRUE]
  agents <- ab_in_data[ab_in_data %in% agents]
  message_agent_names(
    function_name = "administrable_per_os",
    agents = agents,
    ab_group = "administrable_per_os",
    examples = paste0(
      " (such as ",
      vector_or(
        ab_name(
          sample(agents_all,
            size = min(5, length(agents_all)),
            replace = FALSE
          ),
          tolower = TRUE,
          language = NULL
        ),
        quotes = FALSE
      ),
      ")"
    )
  )
  structure(unname(agents),
    class = c("ab_selector", "character")
  )
}

#' @rdname antibiotic_class_selectors
#' @export
administrable_iv <- function(only_sir_columns = FALSE, ...) {
  meet_criteria(only_sir_columns, allow_class = "logical", has_length = 1)
  # get_current_data() has to run each time, for cases where e.g., filter() and select() are used in same call
  # but it only takes a couple of milliseconds
  vars_df <- get_current_data(arg_name = NA, call = -2)
  # to improve speed, get_column_abx() will only run once when e.g. in a select or group call
  ab_in_data <- get_column_abx(vars_df,
    info = FALSE, only_sir_columns = only_sir_columns,
    sort = FALSE, fn = "administrable_iv"
  )
  agents_all <- AMR_env$AB_lookup[which(!is.na(AMR_env$AB_lookup$iv_ddd)), "ab", drop = TRUE]
  agents <- AMR_env$AB_lookup[which(AMR_env$AB_lookup$ab %in% ab_in_data & !is.na(AMR_env$AB_lookup$iv_ddd)), "ab", drop = TRUE]
  agents <- ab_in_data[ab_in_data %in% agents]
  message_agent_names(
    function_name = "administrable_iv",
    agents = agents,
    ab_group = "administrable_iv",
    examples = ""
  )
  structure(unname(agents),
    class = c("ab_selector", "character")
  )
}

#' @rdname antibiotic_class_selectors
#' @inheritParams eucast_rules
#' @details The [not_intrinsic_resistant()] function can be used to only select antibiotic columns that pose no intrinsic resistance for the microorganisms in the data set. For example, if a data set contains only microorganism codes or names of *E. coli* and *K. pneumoniae* and contains a column "vancomycin", this column will be removed (or rather, unselected) using this function. It currently applies `r format_eucast_version_nr(names(EUCAST_VERSION_EXPERT_RULES[1]))` to determine intrinsic resistance, using the [eucast_rules()] function internally. Because of this determination, this function is quite slow in terms of performance.
#' @export
not_intrinsic_resistant <- function(only_sir_columns = FALSE, col_mo = NULL, version_expertrules = 3.3, ...) {
  meet_criteria(only_sir_columns, allow_class = "logical", has_length = 1)
  # get_current_data() has to run each time, for cases where e.g., filter() and select() are used in same call
  # but it only takes a couple of milliseconds
  vars_df <- get_current_data(arg_name = NA, call = -2)
  # to improve speed, get_column_abx() will only run once when e.g. in a select or group call
  ab_in_data <- get_column_abx(vars_df,
    info = FALSE, only_sir_columns = only_sir_columns,
    sort = FALSE, fn = "not_intrinsic_resistant"
  )
  # intrinsic vars
  vars_df_R <- tryCatch(
    sapply(
      eucast_rules(vars_df,
        col_mo = col_mo,
        version_expertrules = version_expertrules,
        rules = "expert",
        info = FALSE
      ),
      function(col) {
        tryCatch(!any(is.na(col)) && all(col == "R"),
          error = function(e) FALSE
        )
      }
    ),
    error = function(e) stop_("in not_intrinsic_resistant(): ", e$message, call = FALSE)
  )
  
  agents <- ab_in_data[ab_in_data %in% names(vars_df_R[which(vars_df_R)])]
  if (length(agents) > 0 &&
      message_not_thrown_before("not_intrinsic_resistant", sort(agents))) {
    agents_formatted <- paste0("'", font_bold(agents, collapse = NULL), "'")
    agents_names <- ab_name(names(agents), tolower = TRUE, language = NULL)
    need_name <- generalise_antibiotic_name(agents) != generalise_antibiotic_name(agents_names)
    agents_formatted[need_name] <- paste0(agents_formatted[need_name], " (", agents_names[need_name], ")")
    message_(
      "For `not_intrinsic_resistant()` removing ",
      ifelse(length(agents) == 1, "column ", "columns "),
      vector_and(agents_formatted, quotes = FALSE, sort = FALSE)
    )
  }
  
  vars_df_R <- names(vars_df_R)[which(!vars_df_R)]
  # find columns that are abx, but also intrinsic R
  out <- unname(intersect(ab_in_data, vars_df_R))
  structure(out,
    class = c("ab_selector", "character")
  )
}

ab_select_exec <- function(function_name,
                           only_sir_columns = FALSE,
                           only_treatable = FALSE,
                           ab_class_args = NULL) {
  # get_current_data() has to run each time, for cases where e.g., filter() and select() are used in same call
  # but it only takes a couple of milliseconds
  vars_df <- get_current_data(arg_name = NA, call = -3)
  # to improve speed, get_column_abx() will only run once when e.g. in a select or group call
  ab_in_data <- get_column_abx(vars_df,
    info = FALSE, only_sir_columns = only_sir_columns,
    sort = FALSE, fn = function_name
  )

  # untreatable drugs
  if (only_treatable == TRUE) {
    untreatable <- AMR_env$AB_lookup[which(AMR_env$AB_lookup$name %like% "-high|EDTA|polysorbate|macromethod|screening|nacubactam"), "ab", drop = TRUE]
    if (any(untreatable %in% names(ab_in_data))) {
      if (message_not_thrown_before(function_name, "ab_class", "untreatable")) {
        warning_(
          "in `", function_name, "()`: some drugs were ignored since they cannot be used for treating patients: ",
          vector_and(
            ab_name(names(ab_in_data)[names(ab_in_data) %in% untreatable],
              language = NULL,
              tolower = TRUE
            ),
            quotes = FALSE,
            sort = TRUE
          ), ". They can be included using `", function_name, "(only_treatable = FALSE)`."
        )
      }
      ab_in_data <- ab_in_data[!names(ab_in_data) %in% untreatable]
    }
  }

  if (length(ab_in_data) == 0) {
    message_("No antimicrobial drugs found in the data.")
    return(NULL)
  }

  if (is.null(ab_class_args) || isTRUE(function_name %in% c("antifungals", "antimycobacterials"))) {
    ab_group <- NULL
    if (isTRUE(function_name == "antifungals")) {
      abx <- AMR_env$AB_lookup$ab[which(AMR_env$AB_lookup$group == "Antifungals")]
    } else if (isTRUE(function_name == "antimycobacterials")) {
      abx <- AMR_env$AB_lookup$ab[which(AMR_env$AB_lookup$group == "Antimycobacterials")]
    } else {
      # their upper case equivalent are vectors with class 'ab', created in data-raw/_pre_commit_checks.R
      # carbapenems() gets its codes from AMR:::AB_CARBAPENEMS
      abx <- get(paste0("AB_", toupper(function_name)), envir = asNamespace("AMR"))
      # manually added codes from add_custom_antimicrobials() must also be supported
      if (length(AMR_env$custom_ab_codes) > 0) {
        custom_ab <- AMR_env$AB_lookup[which(AMR_env$AB_lookup$ab %in% AMR_env$custom_ab_codes), ]
        check_string <- paste0(custom_ab$group, custom_ab$atc_group1, custom_ab$atc_group2)
        if (function_name == "betalactams") {
          find_group <- "beta-lactams"
        } else if (function_name %like% "cephalosporins_") {
          find_group <- gsub("_(.*)$", paste0(" (\\1 gen.)"), function_name)
        } else {
          find_group <- function_name
        }
        abx <- c(abx, custom_ab$ab[which(check_string %like% find_group)])
      }
      ab_group <- function_name
    }
    examples <- paste0(" (such as ", vector_or(
      ab_name(sample(abx, size = min(2, length(abx)), replace = FALSE),
        tolower = TRUE,
        language = NULL
      ),
      quotes = FALSE
    ), ")")
  } else {
    # this for the 'manual' ab_class() function
    abx <- subset(
      AMR_env$AB_lookup,
      group %like% ab_class_args |
        atc_group1 %like% ab_class_args |
        atc_group2 %like% ab_class_args
    )$ab
    ab_group <- find_ab_group(ab_class_args)
    function_name <- "ab_class"
    examples <- paste0(" (such as ", find_ab_names(ab_class_args, 2), ")")
  }

  # get the columns with a group names in the chosen ab class
  agents <- ab_in_data[names(ab_in_data) %in% abx]

  message_agent_names(
    function_name = function_name,
    agents = agents,
    ab_group = ab_group,
    examples = examples,
    ab_class_args = ab_class_args
  )

  structure(unname(agents),
    class = c("ab_selector", "character")
  )
}

#' @method c ab_selector
#' @export
#' @noRd
c.ab_selector <- function(...) {
  structure(unlist(lapply(list(...), as.character)),
    class = c("ab_selector", "character")
  )
}

all_any_ab_selector <- function(type, ..., na.rm = TRUE) {
  cols_ab <- c(...)
  result <- cols_ab[toupper(cols_ab) %in% c("S", "SDD", "I", "R", "N")]
  if (length(result) == 0) {
    message_("Filtering ", type, " of columns ", vector_and(font_bold(cols_ab, collapse = NULL), quotes = "'"), ' to contain value "S", "I" or "R"')
    result <- c("S", "SDD", "I", "R", "N")
  }
  cols_ab <- cols_ab[!cols_ab %in% result]
  df <- get_current_data(arg_name = NA, call = -3)

  if (type == "all") {
    scope_fn <- all
  } else {
    scope_fn <- any
  }

  x_transposed <- as.list(as.data.frame(t(df[, cols_ab, drop = FALSE]), stringsAsFactors = FALSE))
  vapply(
    FUN.VALUE = logical(1),
    X = x_transposed,
    FUN = function(y) scope_fn(y %in% result, na.rm = na.rm),
    USE.NAMES = FALSE
  )
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
  if (isTRUE(na.rm)) {
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
  if (isTRUE(na.rm)) {
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
      message_(
        "Assuming a filter on ", type, " ", length(e1), " ", gsub("[\\(\\)]", "", fn_name),
        ". Wrap around `all()` or `any()` to prevent this note."
      )
    }
  }
  structure(all_any_ab_selector(type = type, e1, e2),
    class = c("ab_selector_any_all", "logical")
  )
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
      message_(
        "Assuming a filter on ", type, " ", length(e1), " ", gsub("[\\(\\)]", "", fn_name),
        ". Wrap around `all()` or `any()` to prevent this note."
      )
    }
  }
  # this is `!=`, so turn around the values
  sir <- c("S", "SDD", "I", "R", "N")
  e2 <- sir[sir != e2]
  structure(all_any_ab_selector(type = type, e1, e2),
    class = c("ab_selector_any_all", "logical")
  )
}

#' @method & ab_selector
#' @export
#' @noRd
`&.ab_selector` <- function(e1, e2) {
  # this is only required for base R, since tidyselect has already implemented this
  # e.g., for: example_isolates[, penicillins() & administrable_per_os()]
  structure(intersect(unclass(e1), unclass(e2)),
    class = c("ab_selector", "character")
  )
}
#' @method | ab_selector
#' @export
#' @noRd
`|.ab_selector` <- function(e1, e2) {
  # this is only required for base R, since tidyselect has already implemented this
  # e.g., for: example_isolates[, penicillins() | administrable_per_os()]
  structure(union(unclass(e1), unclass(e2)),
    class = c("ab_selector", "character")
  )
}

is_any <- function(el1) {
  syscalls <- paste0(trimws2(deparse(sys.calls())), collapse = " ")
  el1 <- gsub("(.*),.*", "\\1", el1)
  syscalls %like% paste0("[^_a-zA-Z0-9]any\\(", "(c\\()?", el1)
}
is_all <- function(el1) {
  syscalls <- paste0(trimws2(deparse(sys.calls())), collapse = " ")
  el1 <- gsub("(.*),.*", "\\1", el1)
  syscalls %like% paste0("[^_a-zA-Z0-9]all\\(", "(c\\()?", el1)
}

find_ab_group <- function(ab_class_args) {
  ab_class_args <- gsub("[^a-zA-Z0-9]", ".*", ab_class_args)
  AMR_env$AB_lookup %pm>%
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
  drugs <- AMR_env$AB_lookup[which((!is.na(AMR_env$AB_lookup$iv_ddd) | !is.na(AMR_env$AB_lookup$oral_ddd)) &
    AMR_env$AB_lookup$name %unlike% " " &
    AMR_env$AB_lookup$group %like% ab_group &
    AMR_env$AB_lookup$ab %unlike% "[0-9]$"), ]$name
  if (length(drugs) < n) {
    # now try it all
    drugs <- AMR_env$AB_lookup[which((AMR_env$AB_lookup$group %like% ab_group |
      AMR_env$AB_lookup$atc_group1 %like% ab_group |
      AMR_env$AB_lookup$atc_group2 %like% ab_group) &
      AMR_env$AB_lookup$ab %unlike% "[0-9]$"), ]$name
  }
  if (length(drugs) == 0) {
    return("??")
  }
  vector_or(
    ab_name(sample(drugs, size = min(n, length(drugs)), replace = FALSE),
      tolower = TRUE,
      language = NULL
    ),
    quotes = FALSE
  )
}

message_agent_names <- function(function_name, agents, ab_group = NULL, examples = "", ab_class_args = NULL, call = NULL) {
  if (message_not_thrown_before(function_name, sort(agents))) {
    if (length(agents) == 0) {
      if (is.null(ab_group)) {
        message_("For `", function_name, "()` no antimicrobial drugs found", examples, ".")
      } else if (ab_group == "administrable_per_os") {
        message_("No orally administrable drugs found", examples, ".")
      } else if (ab_group == "administrable_iv") {
        message_("No IV administrable drugs found", examples, ".")
      } else {
        message_("No antimicrobial drugs of class '", ab_group, "' found", examples, ".")
      }
    } else {
      agents_formatted <- paste0("'", font_bold(agents, collapse = NULL), "'")
      agents_names <- ab_name(names(agents), tolower = TRUE, language = NULL)
      need_name <- generalise_antibiotic_name(agents) != generalise_antibiotic_name(agents_names)
      agents_formatted[need_name] <- paste0(agents_formatted[need_name], " (", agents_names[need_name], ")")
      message_(
        "For `", function_name, "(",
        ifelse(function_name == "ab_class",
          paste0("\"", ab_class_args, "\""),
          ifelse(!is.null(call),
            paste0(deparse(call), collapse = " "),
            ""
          )
        ),
        ")` using ",
        ifelse(length(agents) == 1, "column ", "columns "),
        vector_and(agents_formatted, quotes = FALSE, sort = FALSE)
      )
    }
  }
}
