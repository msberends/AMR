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

#' Define Custom EUCAST Rules
#'
#' Define custom EUCAST rules for your organisation or specific analysis and use the output of this function in [eucast_rules()].
#' @param ... Rules in [formula][base::tilde] notation, see below for instructions, and in *Examples*.
#' @details
#' Some organisations have their own adoption of EUCAST rules. This function can be used to define custom EUCAST rules to be used in the [eucast_rules()] function.
#'
#' ### Basics
#'
#' If you are familiar with the [`case_when()`][dplyr::case_when()] function of the `dplyr` package, you will recognise the input method to set your own rules. Rules must be set using what \R considers to be the 'formula notation'. The rule itself is written *before* the tilde (`~`) and the consequence of the rule is written *after* the tilde:
#'
#' ```r
#' x <- custom_eucast_rules(TZP == "S" ~ aminopenicillins == "S",
#'                          TZP == "R" ~ aminopenicillins == "R")
#' ```
#'
#' These are two custom EUCAST rules: if TZP (piperacillin/tazobactam) is "S", all aminopenicillins (ampicillin and amoxicillin) must be made "S", and if TZP is "R", aminopenicillins must be made "R". These rules can also be printed to the console, so it is immediately clear how they work:
#'
#' ```r
#' x
#' #> A set of custom EUCAST rules:
#' #>
#' #>   1. If TZP is "S" then set to  S :
#' #>      amoxicillin (AMX), ampicillin (AMP)
#' #>
#' #>   2. If TZP is "R" then set to  R :
#' #>      amoxicillin (AMX), ampicillin (AMP)
#' ```
#'
#' The rules (the part *before* the tilde, in above example `TZP == "S"` and `TZP == "R"`) must be evaluable in your data set: it should be able to run as a filter in your data set without errors. This means for the above example that the column `TZP` must exist. We will create a sample data set and test the rules set:
#'
#' ```r
#' df <- data.frame(mo = c("Escherichia coli", "Klebsiella pneumoniae"),
#'                  TZP = as.sir("R"),
#'                  ampi = as.sir("S"),
#'                  cipro = as.sir("S"))
#' df
#' #>                      mo TZP ampi cipro
#' #> 1      Escherichia coli   R    S     S
#' #> 2 Klebsiella pneumoniae   R    S     S
#'
#' eucast_rules(df,
#'              rules = "custom",
#'              custom_rules = x,
#'              info = FALSE,
#'              overwrite = TRUE)
#' #>                      mo TZP ampi cipro
#' #> 1      Escherichia coli   R    R     S
#' #> 2 Klebsiella pneumoniae   R    R     S
#' ```
#'
#' ### Using taxonomic properties in rules
#'
#' There is one exception in columns used for the rules: all column names of the [microorganisms] data set can also be used, but do not have to exist in the data set. These column names are: `r vector_and(colnames(microorganisms), sort = FALSE)`. Thus, this next example will work as well, despite the fact that the `df` data set does not contain a column `genus`:
#'
#' ```r
#' y <- custom_eucast_rules(
#'   TZP == "S" & genus == "Klebsiella" ~ aminopenicillins == "S",
#'   TZP == "R" & genus == "Klebsiella" ~ aminopenicillins == "R"
#' )
#'
#' eucast_rules(df,
#'              rules = "custom",
#'              custom_rules = y,
#'              info = FALSE,
#'              overwrite = TRUE)
#' #>                      mo TZP ampi cipro
#' #> 1      Escherichia coli   R    S     S
#' #> 2 Klebsiella pneumoniae   R    R     S
#' ```
#'
#' ### Sharing rules among multiple users
#'
#' The rules set (the `y` object in this case) could be exported to a shared file location using [saveRDS()] if you collaborate with multiple users. The custom rules set could then be imported using [readRDS()].
#'
#' ### Usage of multiple antimicrobials and antimicrobial group names
#'
#' You can define antimicrobial groups instead of single antimicrobials for the rule consequence, which is the part *after* the tilde (~). In the examples above, the antimicrobial group `aminopenicillins` includes both ampicillin and amoxicillin.
#'
#' Rules can also be applied to multiple antimicrobials and antimicrobial groups simultaneously. Use the `c()` function to combine multiple antimicrobials. For instance, the following example sets all aminopenicillins and ureidopenicillins to "R" if column TZP (piperacillin/tazobactam) is "R":
#'
#' ```r
#' x <- custom_eucast_rules(TZP == "R" ~ c(aminopenicillins, ureidopenicillins) == "R")
#' x
#' #> A set of custom EUCAST rules:
#' #>
#' #>   1. If TZP is "R" then set to "R":
#' #>      amoxicillin (AMX), ampicillin (AMP), azlocillin (AZL), mezlocillin (MEZ), piperacillin (PIP), piperacillin/tazobactam (TZP)
#' ```
#'
#' These `r length(DEFINED_AB_GROUPS)` antimicrobial groups are allowed in the rules (case-insensitive) and can be used in any combination:
#'
#' `r paste0("  * ", sapply(DEFINED_AB_GROUPS, function(x) paste0(tolower(gsub("^AB_", "", x)), "\\cr(", vector_and(ab_name(eval(parse(text = x), envir = asNamespace("AMR")), language = NULL, tolower = TRUE), quotes = FALSE), ")"), USE.NAMES = FALSE), "\n", collapse = "")`
#' @returns A [list] containing the custom rules
#' @export
#' @examples
#' x <- custom_eucast_rules(
#'   AMC == "R" & genus == "Klebsiella" ~ aminopenicillins == "R",
#'   AMC == "I" & genus == "Klebsiella" ~ aminopenicillins == "I"
#' )
#' x
#'
#' # run the custom rule set (verbose = TRUE will return a logbook instead of the data set):
#' eucast_rules(example_isolates,
#'   rules = "custom",
#'   custom_rules = x,
#'   info = FALSE,
#'   overwrite = TRUE,
#'   verbose = TRUE
#' )
#'
#' # combine rule sets
#' x2 <- c(
#'   x,
#'   custom_eucast_rules(TZP == "R" ~ carbapenems == "R")
#' )
#' x2
custom_eucast_rules <- function(...) {
  dots <- tryCatch(list(...),
    error = function(e) "error"
  )
  stop_if(
    identical(dots, "error"),
    "rules must be a valid formula inputs (e.g., using '~'), see `?custom_eucast_rules`"
  )
  n_dots <- length(dots)
  stop_if(n_dots == 0, "no custom rules were set. Please read the documentation using `?custom_eucast_rules`.")
  out <- vector("list", n_dots)
  for (i in seq_len(n_dots)) {
    stop_ifnot(
      inherits(dots[[i]], "formula"),
      "rule ", i, " must be a valid formula input (e.g., using '~'), see `?custom_eucast_rules`"
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
    # format nicely, setting spaces around operators
    qry <- gsub(" *([&|+-/*^><==]+) *", " \\1 ", qry)
    qry <- gsub(" ?, ?", ", ", qry)
    qry <- gsub("'", "\"", qry, fixed = TRUE)
    out[[i]]$query <- as.expression(qry)

    # Resulting rule
    result <- dots[[i]][[3]]
    stop_ifnot(
      deparse(result) %like% "==",
      "the result of rule ", i, " (the part after the `~`) must contain `==`, such as in `... ~ ampicillin == \"R\"`, see `?custom_eucast_rules`"
    )
    result_group <- as.character(result)[[2]]
    result_group <- as.character(str2lang(result_group))
    result_group <- result_group[result_group != "c"]
    result_group_agents <- character(0)
    for (j in seq_len(length(result_group))) {
      if (paste0("AB_", toupper(result_group[j]), "S") %in% DEFINED_AB_GROUPS) {
        # support for e.g. 'aminopenicillin' if user meant 'aminopenicillins'
        result_group[j] <- paste0(result_group[j], "s")
      }
      if (paste0("AB_", toupper(result_group[j])) %in% DEFINED_AB_GROUPS) {
        result_group_agents <- c(
          result_group_agents,
          eval(parse(text = paste0("AB_", toupper(result_group[j]))), envir = asNamespace("AMR"))
        )
      } else {
        out_group <- tryCatch(
          suppressWarnings(as.ab(result_group[j],
            fast_mode = TRUE,
            flag_multiple_results = FALSE
          )),
          error = function(e) NA_character_
        )
        if (!all(is.na(out_group))) {
          result_group_agents <- c(result_group_agents, out_group)
        }
      }
    }
    result_group_agents <- result_group_agents[!is.na(result_group_agents)]

    stop_if(
      length(result_group_agents) == 0,
      "this result of rule ", i, " could not be translated to a single antimicrobial drug/group: \"",
      as.character(result)[[2]], "\".\n\nThe input can be a name or code of an antimicrobial drug, or be one of: ",
      vector_or(tolower(gsub("AB_", "", DEFINED_AB_GROUPS)), quotes = FALSE), "."
    )
    result_value <- as.character(result)[[3]]
    result_value[result_value == "NA"] <- NA
    stop_ifnot(
      result_value %in% c("S", "SDD", "I", "R", "NI", NA),
      "the resulting value of rule ", i, " must be either \"S\", \"SDD\", \"I\", \"R\", \"NI\" or NA"
    )
    result_value <- as.sir(result_value)

    out[[i]]$result_group <- result_group_agents
    out[[i]]$result_value <- result_value
  }

  names(out) <- paste0("rule", seq_len(n_dots))
  set_clean_class(out, new_class = c("custom_eucast_rules", "list"))
}

#' @method c custom_eucast_rules
#' @noRd
#' @export
c.custom_eucast_rules <- function(x, ...) {
  if (length(list(...)) == 0) {
    return(x)
  }
  out <- unclass(x)
  for (e in list(...)) {
    out <- c(out, unclass(e))
  }
  names(out) <- paste0("rule", seq_len(length(out)))
  set_clean_class(out, new_class = c("custom_eucast_rules", "list"))
}

#' @method as.list custom_eucast_rules
#' @noRd
#' @export
as.list.custom_eucast_rules <- function(x, ...) {
  c(x, ...)
}

#' @method print custom_eucast_rules
#' @export
#' @noRd
print.custom_eucast_rules <- function(x, ...) {
  cat("A set of custom EUCAST rules:\n")
  for (i in seq_len(length(x))) {
    rule <- x[[i]]
    rule$query <- format_custom_query_rule(rule$query)
    if (is.na(rule$result_value)) {
      val <- font_red("<NA>")
    } else if (rule$result_value == "R") {
      val <- font_rose_bg(" R ")
    } else if (rule$result_value == "S") {
      val <- font_green_bg(" S ")
    } else {
      val <- font_orange_bg(" I ")
    }
    agents <- paste0(
      font_blue(ab_name(rule$result_group, language = NULL, tolower = TRUE),
        collapse = NULL
      ),
      " (", rule$result_group, ")"
    )
    agents <- sort(agents)
    rule_if <- word_wrap(
      paste0(
        i, ". ", font_bold("If "), font_blue(rule$query), font_bold(" then "),
        "set to {result}:"
      ),
      extra_indent = 5
    )
    rule_if <- gsub("{result}", val, rule_if, fixed = TRUE)
    rule_then <- paste0("     ", word_wrap(paste0(agents, collapse = ", "), extra_indent = 5))
    cat("\n  ", rule_if, "\n", rule_then, "\n", sep = "")
  }
}
