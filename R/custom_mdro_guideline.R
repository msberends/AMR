# ==================================================================== #
# TITLE:                                                               #
# AMR: An R Package for Working with Antimicrobial Resistance Data     #
#                                                                      #
# SOURCE CODE:                                                         #
# https://github.com/msberends/AMR                                     #
#                                                                      #
# PLEASE CITE THIS SOFTWARE AS:                                        #
# Berends MS, Luz CF, Friedrich AW, et al. (2022).                     #
# AMR: An R Package for Working with Antimicrobial Resistance Data.    #
# Journal of Statistical Software, 104(3), 1-31.                       #
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
# how to conduct AMR data analysis: https://amr-for-r.org              #
# ==================================================================== #

#' Define Custom MDRO Guideline
#'
#' Define custom a MDRO guideline for your organisation or specific analysis and use the output of this function in [mdro()].
#' @param ... Guideline rules in [formula][base::tilde] notation, see below for instructions, and in *Examples*.
#' @inheritParams mdro
#' @details
#' Using a custom MDRO guideline is of importance if you have custom rules to determine MDROs in your hospital, e.g., rules that are dependent on ward, state of contact isolation or other variables in your data.
#' @section How it works:
#'
#' ### Basics
#'
#' If you are familiar with the [`case_when()`][dplyr::case_when()] function of the `dplyr` package, you will recognise the input method to set your own rules. Rules must be set using what \R considers to be the 'formula notation'. The rule itself is written *before* the tilde (`~`) and the consequence of the rule is written *after* the tilde:
#'
#' ```r
#' custom <- custom_mdro_guideline(CIP == "R" & age > 60 ~ "Elderly Type A",
#'                                 ERY == "R" & age > 60 ~ "Elderly Type B")
#' ```
#'
#' If a row/an isolate matches the first rule, the value after the first `~` (in this case *'Elderly Type A'*) will be set as MDRO value. Otherwise, the second rule will be tried and so on. The number of rules is unlimited.
#'
#' You can print the rules set in the console for an overview. Colours will help reading it if your console supports colours.
#'
#' ```r
#' custom
#' #> A set of custom MDRO rules:
#' #>   1. If CIP is R and age is higher than 60 then: Elderly Type A
#' #>   2. If ERY is R and age is higher than 60 then: Elderly Type B
#' #>   3. Otherwise: Negative
#'
#' #> Unmatched rows will return NA.
#' #> Results will be of class 'factor', with ordered levels: Negative < Elderly Type A < Elderly Type B
#' ```
#'
#' The outcome of the function can be used for the `guideline` argument in the [mdro()] function:
#'
#' ```r
#' x <- mdro(example_isolates, guideline = custom)
#' #> Determining MDROs based on custom rules, resulting in factor levels: Negative < Elderly Type A < Elderly Type B.
#' #> - Custom MDRO rule 1: CIP == "R" & age > 60 (198 rows matched)
#' #> - Custom MDRO rule 2: ERY == "R" & age > 60 (732 rows matched)
#' #> => Found 930 custom defined MDROs out of 2000 isolates (46.5%)
#'
#' table(x)
#' #> x
#' #>       Negative  Elderly Type A  Elderly Type B
#' #>           1070             198             732
#' ```
#'
#' Rules can also be combined with other custom rules by using [c()]:
#'
#' ```r
#' x <- mdro(example_isolates,
#'           guideline = c(custom,
#'                         custom_mdro_guideline(ERY == "R" & age > 50 ~ "Elderly Type C")))
#' #> Determining MDROs based on custom rules, resulting in factor levels: Negative < Elderly Type A < Elderly Type B < Elderly Type C.
#' #> - Custom MDRO rule 1: CIP == "R" & age > 60 (198 rows matched)
#' #> - Custom MDRO rule 2: ERY == "R" & age > 60 (732 rows matched)
#' #> - Custom MDRO rule 3: ERY == "R" & age > 50 (109 rows matched)
#' #> => Found 1039 custom defined MDROs out of 2000 isolates (52.0%)
#'
#' table(x)
#' #> x
#' #>       Negative  Elderly Type A  Elderly Type B  Elderly Type C
#' #>            961             198             732             109
#' ```
#'
#' ### Sharing rules among multiple users
#'
#' The rules set (the `custom` object in this case) could be exported to a shared file location using [saveRDS()] if you collaborate with multiple users. The custom rules set could then be imported using [readRDS()].
#'
#' ### Usage of multiple antimicrobials and antimicrobial group names
#'
#' You can define antimicrobial groups instead of single antimicrobials for the rule itself, which is the part *before* the tilde (~). Use [any()] or [all()] to specify the scope of the antimicrobial group:
#'
#' ```r
#' custom_mdro_guideline(
#'   AMX == "R"                       ~ "My MDRO #1",
#'   any(cephalosporins_2nd() == "R") ~ "My MDRO #2",
#'   all(glycopeptides() == "R")      ~ "My MDRO #3"
#' )
#' ```
#'
#' These `r length(DEFINED_AB_GROUPS)` antimicrobial groups are allowed in the rules (case-insensitive) and can be used in any combination:
#'
#' `r paste0("  * ", sapply(DEFINED_AB_GROUPS, function(x) paste0(tolower(gsub("^AB_", "", x)), "\\cr(", vector_and(ab_name(eval(parse(text = x), envir = asNamespace("AMR")), language = NULL, tolower = TRUE), quotes = FALSE), ")"), USE.NAMES = FALSE), "\n", collapse = "")`
#' @returns A [list] containing the custom rules
#' @rdname custom_mdro_guideline
#' @export
#' @examples
#' x <- custom_mdro_guideline(
#'   CIP == "R" & age > 60 ~ "Elderly Type A",
#'   ERY == "R" & age > 60 ~ "Elderly Type B"
#' )
#' x
#'
#' # run the custom rule set (verbose = TRUE will return a logbook instead of the data set):
#' out <- mdro(example_isolates, guideline = x)
#' table(out)
#'
#' out <- mdro(example_isolates, guideline = x, verbose = TRUE)
#' head(out)
#'
#' # you can create custom guidelines using selectors (see ?antimicrobial_selectors)
#' my_guideline <- custom_mdro_guideline(
#'   AMX == "R" ~ "Custom MDRO 1",
#'   all(cephalosporins_2nd() == "R") ~ "Custom MDRO 2"
#' )
#' my_guideline
#'
#' out <- mdro(example_isolates, guideline = my_guideline)
#' table(out)
custom_mdro_guideline <- function(..., as_factor = TRUE) {
  meet_criteria(as_factor, allow_class = "logical", has_length = 1)

  dots <- tryCatch(list(...),
    error = function(e) "error"
  )
  stop_if(
    identical(dots, "error"),
    "rules must be a valid formula inputs (e.g., using '~'), see `?mdro`"
  )
  n_dots <- length(dots)
  stop_if(n_dots == 0, "no custom rules were set. Please read the documentation using `?mdro`.")
  out <- vector("list", n_dots)
  for (i in seq_len(n_dots)) {
    stop_ifnot(
      inherits(dots[[i]], "formula"),
      "rule ", i, " must be a valid formula input (e.g., using '~'), see `?mdro`"
    )

    # Query
    qry <- dots[[i]][[2]]
    if (inherits(qry, "call")) {
      qry <- as.expression(qry)
    }
    qry <- as.character(qry)
    # these will prevent vectorisation, so replace them:
    qry <- gsub("&&", "&", qry, fixed = TRUE)
    qry <- gsub("||", "|", qry, fixed = TRUE)
    # support filter()-like writing: custom_mdro_guideline('CIP == "R", AMX == "S"' ~ "result 1")
    qry <- gsub(" *, *", " & ", qry)
    # format nicely, setting spaces around operators
    qry <- gsub(" *([&|+-/*^><==]+) *", " \\1 ", qry)
    qry <- gsub("'", "\"", qry, fixed = TRUE)
    qry <- as.expression(qry)
    out[[i]]$query <- qry

    # Value
    val <- tryCatch(eval(dots[[i]][[3]]), error = function(e) NULL)
    stop_if(is.null(val), "rule ", i, " must return a valid value, it now returns an error: ", tryCatch(eval(dots[[i]][[3]]), error = function(e) e$message))
    stop_if(length(val) > 1, "rule ", i, " must return a value of length 1, not ", length(val))
    out[[i]]$value <- as.character(val)
  }

  names(out) <- paste0("rule", seq_len(n_dots))
  out <- set_clean_class(out, new_class = c("custom_mdro_guideline", "list"))
  attr(out, "values") <- unname(c("Negative", vapply(FUN.VALUE = character(1), unclass(out), function(x) x$value)))
  attr(out, "as_factor") <- as_factor
  out
}

#' @method c custom_mdro_guideline
#' @rdname custom_mdro_guideline
#' @export
c.custom_mdro_guideline <- function(x, ..., as_factor = NULL) {
  if (length(list(...)) == 0) {
    return(x)
  }
  if (!is.null(as_factor)) {
    meet_criteria(as_factor, allow_class = "logical", has_length = 1)
  } else {
    as_factor <- attributes(x)$as_factor
  }
  for (g in list(...)) {
    stop_ifnot(inherits(g, "custom_mdro_guideline"),
      "for combining custom MDRO guidelines, all rules must be created with `custom_mdro_guideline()`",
      call = FALSE
    )
    vals <- attributes(x)$values
    if (!all(attributes(g)$values %in% vals)) {
      vals <- unname(unique(c(vals, attributes(g)$values)))
    }
    attributes(g) <- NULL
    x <- c(unclass(x), unclass(g))
    attr(x, "values") <- vals
  }
  names(x) <- paste0("rule", seq_len(length(x)))
  x <- set_clean_class(x, new_class = c("custom_mdro_guideline", "list"))
  attr(x, "values") <- vals
  attr(x, "as_factor") <- as_factor
  x
}

#' @method as.list custom_mdro_guideline
#' @noRd
#' @export
as.list.custom_mdro_guideline <- function(x, ...) {
  c(x, ...)
}

#' @method print custom_mdro_guideline
#' @rdname custom_mdro_guideline
#' @export
print.custom_mdro_guideline <- function(x, ...) {
  cat("A set of custom MDRO rules:\n")
  for (i in seq_len(length(x))) {
    rule <- x[[i]]
    rule$query <- format_custom_query_rule(rule$query)
    cat("  ", i, ". ", font_bold("If "), font_blue(rule$query), font_bold(" then: "), font_red(rule$value), "\n", sep = "")
  }
  cat("  ", i + 1, ". ", font_bold("Otherwise: "), font_red(paste0("Negative")), "\n", sep = "")
  cat("\nUnmatched rows will return ", font_red("NA"), ".\n", sep = "")
  if (isTRUE(attributes(x)$as_factor)) {
    cat("Results will be of class 'factor', with ordered levels: ", paste0(attributes(x)$values, collapse = " < "), "\n", sep = "")
  } else {
    cat("Results will be of class 'character'.\n")
  }
}

run_custom_mdro_guideline <- function(df, guideline, info) {
  n_dots <- length(guideline)
  stop_if(n_dots == 0, "no custom guidelines set", call = -2)
  out <- character(length = NROW(df))
  reasons <- character(length = NROW(df))
  for (i in seq_len(n_dots)) {
    qry <- tryCatch(eval(parse(text = guideline[[i]]$query), envir = df, enclos = parent.frame()),
      error = function(e) {
        AMR_env$err_msg <- e$message
        return("error")
      }
    )
    if (identical(qry, "error")) {
      warning_("in `custom_mdro_guideline()`: rule ", i,
        " (`", as.character(guideline[[i]]$query), "`) was ignored because of this error message: ",
        AMR_env$err_msg,
        call = FALSE,
        add_fn = font_red
      )
      next
    }
    stop_ifnot(is.logical(qry), "in custom_mdro_guideline(): rule ", i, " (`", guideline[[i]]$query,
      "`) must return `TRUE` or `FALSE`, not ",
      format_class(class(qry), plural = FALSE),
      call = FALSE
    )

    new_mdros <- which(qry == TRUE & out == "")

    if (isTRUE(info)) {
      cat(word_wrap(
        "- Custom MDRO rule ", i, ": `", as.character(guideline[[i]]$query),
        "` (", length(new_mdros), " rows matched)"
      ), "\n", sep = "")
    }
    val <- guideline[[i]]$value
    out[new_mdros] <- val
    reasons[new_mdros] <- paste0(
      "matched rule ",
      gsub("rule", "", names(guideline)[i], fixed = TRUE), ": ", as.character(guideline[[i]]$query)
    )
  }
  out[out == ""] <- "Negative"
  reasons[out == "Negative"] <- "no rules matched"

  if (isTRUE(attributes(guideline)$as_factor)) {
    out <- factor(out, levels = attributes(guideline)$values, ordered = TRUE)
  }

  all_nonsusceptible_columns <- as.data.frame(t(df[, is.sir(df), drop = FALSE] == "R"))
  all_nonsusceptible_columns <- vapply(
    FUN.VALUE = character(1),
    all_nonsusceptible_columns,
    function(x) paste0(rownames(all_nonsusceptible_columns)[which(x)], collapse = ", ")
  )
  all_nonsusceptible_columns[is.na(out)] <- NA_character_

  data.frame(
    row_number = seq_len(NROW(df)),
    MDRO = out,
    reason = reasons,
    all_nonsusceptible_columns = all_nonsusceptible_columns,
    stringsAsFactors = FALSE
  )
}
