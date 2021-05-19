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
#' These functions help to filter and select columns with antibiotic test results that are of a specific antibiotic class, without the need to define the columns or antibiotic abbreviations. \strong{\Sexpr{ifelse(as.double(R.Version()$major) + (as.double(R.Version()$minor) / 10) < 3.2, paste0("NOTE: THESE FUNCTIONS DO NOT WORK ON YOUR CURRENT R VERSION. These functions require R version 3.2 or later - you have ", R.version.string, "."), "")}}
#' @inheritSection lifecycle Stable Lifecycle
#' @param ab_class an antimicrobial class, like `"carbapenems"`. The columns `group`, `atc_group1` and `atc_group2` of the [antibiotics] data set will be searched (case-insensitive) for this value.
#' @param only_rsi_columns a [logical] to indicate whether only columns of class `<rsi>` must be selected (defaults to `FALSE`), see [as.rsi()]
#' @details \strong{\Sexpr{ifelse(as.double(R.Version()$major) + (as.double(R.Version()$minor) / 10) < 3.2, paste0("NOTE: THESE FUNCTIONS DO NOT WORK ON YOUR CURRENT R VERSION. These functions require R version 3.2 or later - you have ", R.version.string, "."), "")}}
#' 
#' 
#' These functions can be used in data set calls for selecting columns and filtering rows, see *Examples*. They support base R, but work more convenient in dplyr functions such as [`select()`][dplyr::select()], [`filter()`][dplyr::filter()] and [`summarise()`][dplyr::summarise()].
#' 
#' All columns in the data in which these functions are called will be searched for known antibiotic names, abbreviations, brand names, and codes (ATC, EARS-Net, WHO, etc.) in the [antibiotics] data set. This means that a selector like e.g. [aminoglycosides()] will pick up column names like 'gen', 'genta', 'J01GB03', 'tobra', 'Tobracin', etc.
#' 
#' The group of betalactams consists of all carbapenems, cephalosporins and penicillins.
#' @rdname antibiotic_class_selectors
#' @name antibiotic_class_selectors
#' @export
#' @inheritSection AMR Reference Data Publicly Available
#' @inheritSection AMR Read more on Our Website!
#' @examples 
#' # `example_isolates` is a data set available in the AMR package.
#' # See ?example_isolates.
#' 
#' # Base R ------------------------------------------------------------------
#' 
#' # select columns 'IPM' (imipenem) and 'MEM' (meropenem)
#' example_isolates[, carbapenems()]
#' 
#' # select columns 'mo', 'AMK', 'GEN', 'KAN' and 'TOB'
#' example_isolates[, c("mo", aminoglycosides())]
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
#' # filter + select in one go: get penicillins in carbapenems-resistant strains
#' example_isolates[any(carbapenems() == "R"), penicillins()]
#' 
#' 
#' # dplyr -------------------------------------------------------------------
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
#'  # any() and all() work in dplyr's filter() too:
#'  example_isolates %>% 
#'     filter(any(aminoglycosides() == "R"),
#'            all(cephalosporins_2nd() == "R"))
#'     
#'  # also works with c():
#'  example_isolates %>% 
#'     filter(any(c(carbapenems(), aminoglycosides()) == "R"))
#'     
#'  # not setting any/all will automatically apply all():
#'  example_isolates %>% 
#'     filter(aminoglycosides() == "R")
#'  #> i Assuming a filter on all 4 aminoglycosides.
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
#'   # with dplyr 1.0.0 and higher (that adds 'across()'), this is all equal:
#'   # (though the row names on the first are more correct)
#'   example_isolates[carbapenems() == "R", ]
#'   example_isolates %>% filter(carbapenems() == "R")
#'   example_isolates %>% filter(across(carbapenems(), ~.x == "R"))
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
betalactams <- function(only_rsi_columns = FALSE) {
  ab_selector("carbapenem|cephalosporin|penicillin", function_name = "betalactams", only_rsi_columns = only_rsi_columns)
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
  
  if (current_R_older_than(3.2)) {
    warning_("antibiotic class selectors such as ", function_name, 
             "() require R version 3.2 or later - you have ", R.version.string,
             call = FALSE)
    return(NULL)
  }
  
  vars_df <- get_current_data(arg_name = NA, call = -3)

  # improve speed here so it will only run once when e.g. in one select call
  if (!identical(pkg_env$ab_selector, unique_call_id())) {
    ab_in_data <- get_column_abx(vars_df, info = FALSE, only_rsi_columns = only_rsi_columns, sort = FALSE)
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
      message_("For `", function_name, "(", ifelse(function_name == "ab_class", paste0("\"", ab_class, "\""), ""), ")` using ",
               ifelse(length(agents) == 1, "column: ", "columns: "),
               vector_and(agents_formatted, quotes = FALSE))
    }
    remember_thrown_message(function_name)
  }
  
  if (!is.null(attributes(vars_df)$type) &&
      attributes(vars_df)$type %in% c("dplyr_cur_data_all", "base_R") &&
      !any(as.character(sys.calls()) %like% paste0("(across|if_any|if_all)\\((c\\()?[a-z(), ]*", function_name))) {
    structure(unname(agents),
              class = c("ab_selector", "character"))
  } else {
    # don't return with "ab_selector" class if method is a dplyr selector,
    # dplyr::select() will complain:
    # > Subscript has the wrong type `ab_selector`.
    # > It must be numeric or character.
    unname(agents)
  }
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
  # this is all() for 
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
  # keep only the ... in c(...)
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
  # keep only the ... in c(...)
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

is_any <- function(el1) {
  syscall <- paste0(trimws(deparse(sys.calls()[[1]])), collapse = " ")
  el1 <- gsub("(.*),.*", "\\1", el1)
  syscall %like% paste0("[^_a-zA-Z0-9]any\\(", "(c\\()?", el1)
}
is_all <- function(el1) {
  syscall <- paste0(trimws(deparse(sys.calls()[[1]])), collapse = " ")
  el1 <- gsub("(.*),.*", "\\1", el1)
  syscall %like% paste0("[^_a-zA-Z0-9]all\\(", "(c\\()?", el1)
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
