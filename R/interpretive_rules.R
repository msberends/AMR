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
# Center Groningen in the Netherlands, in collaboration with many      #
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

# ====================================================== #
# || Change the EUCAST version numbers in R/globals.R || #
# ====================================================== #

format_eucast_version_nr <- function(version, markdown = TRUE) {
  # for documentation - adds title, version number, year and url in markdown language
  lst <- c(EUCAST_VERSION_BREAKPOINTS, EUCAST_VERSION_EXPECTED_PHENOTYPES, EUCAST_VERSION_EXPERT_RULES)
  version <- format(unique(version), nsmall = 1)
  txt <- character(0)
  for (i in seq_len(length(version))) {
    v <- version[i]
    if (markdown == TRUE) {
      txt <- c(txt, paste0(
        "[", lst[[v]]$title, " ", lst[[v]]$version_txt, "](", lst[[v]]$url, ")",
        " (", lst[[v]]$year, ")"
      ))
    } else {
      txt <- c(txt, paste0(
        lst[[version]]$title, " ", lst[[v]]$version_txt,
        " (", lst[[v]]$year, ")"
      ))
    }
  }
  vector_and(txt, quotes = FALSE)
}

#' Apply Interpretive Rules
#'
#' @description
#' Apply rules from clinical breakpoints notes and expected resistant phenotypes as defined by e.g. the European Committee on Antimicrobial Susceptibility Testing (EUCAST, <https://www.eucast.org>), see *Source*. Use [eucast_dosage()] to get a [data.frame] with advised dosages of a certain bug-drug combination, which is based on the [dosage] data set.
#'
#' To improve the interpretation of the antibiogram before CLSI/EUCAST interpretive rules are applied, some AMR-specific rules can be applied at default, see *Details*.
#' @param x A data set with antimicrobials columns, such as `amox`, `AMX` and `AMC`.
#' @param info A [logical] to indicate whether progress should be printed to the console - the default is only print while in interactive sessions.
#' @param guideline A guideline name, either "EUCAST" (default) or "CLSI". This can be set with the package option [`AMR_guideline`][AMR-options].
#' @param rules A [character] vector that specifies which rules should be applied. Must be one or more of `"breakpoints"`, `"expected_phenotypes"`, `"expert"`, `"other"`, `"custom"`, `"all"`, and defaults to `c("breakpoints", "expected_phenotypes")`. The default value can be set to another value using the package option [`AMR_interpretive_rules`][AMR-options]: `options(AMR_interpretive_rules = "all")`. If using `"custom"`, be sure to fill in argument `custom_rules` too. Custom rules can be created with [custom_interpretive_rules()].
#' @param verbose A [logical] to turn Verbose mode on and off (default is off). In Verbose mode, the function does not apply rules to the data, but instead returns a data set in logbook form with extensive info about which rows and columns would be effected and in which way. Using Verbose mode takes a lot more time.
#' @param version_breakpoints The version number to use for the EUCAST Clinical Breakpoints guideline. Can be `r vector_or(names(EUCAST_VERSION_BREAKPOINTS), documentation = TRUE, reverse = TRUE)`.
#' @param version_expected_phenotypes The version number to use for the EUCAST Expected Phenotypes. Can be `r vector_or(names(EUCAST_VERSION_EXPECTED_PHENOTYPES), documentation = TRUE, reverse = TRUE)`.
#' @param version_expertrules The version number to use for the EUCAST Expert Rules and Intrinsic Resistance guideline. Can be `r vector_or(names(EUCAST_VERSION_EXPERT_RULES), documentation = TRUE, reverse = TRUE)`.
#' @param ampc_cephalosporin_resistance (only applies when `rules` contains `"expert"` or `"all"`) a [character] value that should be applied to cefotaxime, ceftriaxone and ceftazidime for AmpC de-repressed cephalosporin-resistant mutants - the default is `NA`. Currently only works when `version_expertrules` is `3.2` and higher; these versions of '*EUCAST Expert Rules on Enterobacterales*' state that results of cefotaxime, ceftriaxone and ceftazidime should be reported with a note, or results should be suppressed (emptied) for these three drugs. A value of `NA` (the default) for this argument will remove results for these three drugs, while e.g. a value of `"R"` will make the results for these drugs resistant. Use `NULL` or `FALSE` to not alter results for these three drugs of AmpC de-repressed cephalosporin-resistant mutants. Using `TRUE` is equal to using `"R"`. \cr For *EUCAST Expert Rules* v3.2, this rule applies to: `r vector_and(gsub("[^a-zA-Z ]+", "", unlist(strsplit(INTERPRETIVE_RULES_DF[which(INTERPRETIVE_RULES_DF$reference.version %in% c(3.2, 3.3) & INTERPRETIVE_RULES_DF$reference.rule %like% "ampc"), "this_value"][1], "|", fixed = TRUE))), quotes = "*")`.
#' @param ... Column names of antimicrobials. To automatically detect antimicrobial column names, do not provide any named arguments; [guess_ab_col()] will then be used for detection. To manually specify a column, provide its name (case-insensitive) as an argument, e.g. `AMX = "amoxicillin"`. To skip a specific antimicrobial, set it to `NULL`, e.g. `TIC = NULL` to exclude ticarcillin. If a manually defined column does not exist in the data, it will be skipped with a warning.
#' @param ab Any (vector of) text that can be coerced to a valid antimicrobial drug code with [as.ab()].
#' @param administration Route of administration, either `r vector_or(dosage$administration, documentation = TRUE)`.
#' @param only_sir_columns A [logical] to indicate whether only antimicrobial columns must be included that were transformed to class [sir][as.sir()] on beforehand. Defaults to `FALSE` if no columns of `x` have a class [sir][as.sir()].
#' @param custom_rules Custom rules to apply, created with [custom_interpretive_rules()].
#' @param overwrite A [logical] indicating whether to overwrite existing SIR values (default: `FALSE`). When `FALSE`, only non-SIR values are modified (i.e., any value that is not already S, I or R). To ensure compliance with EUCAST guidelines, **this should remain** `FALSE`, as EUCAST notes often state that an organism "should be tested for susceptibility to individual agents or be reported resistant".
#' @param add_if_missing A [logical] indicating whether rules should also be applied to missing (`NA`) values (default: `TRUE`). When `FALSE`, rules are only applied to cells that already contain an SIR value; cells with `NA` are left untouched. This is particularly useful when using `overwrite = TRUE` with custom rules and you want to update reported results without imputing values for untested drugs.
#' @inheritParams first_isolate
#' @details
#' **Note:** This function does not translate MIC or disk values to SIR values. Use [as.sir()] for that. \cr
#' **Note:** When ampicillin (AMP, J01CA01) is not available but amoxicillin (AMX, J01CA04) is, the latter will be used for all rules where there is a dependency on ampicillin. These drugs are interchangeable when it comes to expression of antimicrobial resistance. \cr
#'
#' The file containing all interpretive rules is located here: <https://github.com/msberends/AMR/blob/main/data-raw/interpretive_rules.tsv>.  **Note:** Old taxonomic names are replaced with the current taxonomy where applicable. For example, *Ochrobactrum anthropi* was renamed to *Brucella anthropi* in 2020; older interpretive rules may not yet contain this new taxonomic name. The `AMR` package contains the full microbial taxonomy updated until `r documentation_date(max(TAXONOMY_VERSION$GBIF$accessed_date, TAXONOMY_VERSION$LPSN$accessed_date))`, see [microorganisms].
#'
#' ### Custom Rules
#'
#' Custom rules can be created using [custom_interpretive_rules()], e.g.:
#'
#' ```r
#' x <- custom_interpretive_rules(AMC == "R" & genus == "Klebsiella" ~ aminopenicillins == "R",
#'                                AMC == "I" & genus == "Klebsiella" ~ aminopenicillins == "I")
#'
#' interpretive_rules(example_isolates, rules = "custom", custom_rules = x)
#' ```
#'
#' ### 'Other' Rules
#'
#' Before further processing, two non-EUCAST rules about drug combinations can be applied to improve the efficacy of the EUCAST rules, and the reliability of your data (analysis). These rules are:
#'
#' 1. A drug **with** enzyme inhibitor will be set to S if the same drug **without** enzyme inhibitor is S
#' 2. A drug **without** enzyme inhibitor will be set to R if the same drug **with** enzyme inhibitor is R
#'
#' Important examples include amoxicillin and amoxicillin/clavulanic acid, and trimethoprim and trimethoprim/sulfamethoxazole. Needless to say, for these rules to work, both drugs must be available in the data set.
#'
#' Since these rules are not officially approved by EUCAST, they are not applied at default. To use these rules, include `"other"` to the `rules` argument, or use `interpretive_rules(..., rules = "all")`. You can also set the package option [`AMR_interpretive_rules`][AMR-options], i.e. run `options(AMR_interpretive_rules = "all")`.
#' @aliases EUCAST
#' @rdname interpretive_rules
#' @export
#' @return The input of `x`, possibly with edited values of antimicrobials. Or, if `verbose = TRUE`, a [data.frame] with all original and new values of the affected bug-drug combinations.
#' @references
#' - EUCAST Expert Rules. Version 2.0, 2012.\cr
#'   Leclercq et al. **EUCAST expert rules in antimicrobial susceptibility testing.** *Clin Microbiol Infect.* 2013;19(2):141-60; \doi{https://doi.org/10.1111/j.1469-0691.2011.03703.x}
#' - EUCAST Expected Phenotypes. [(link)](https://www.eucast.org/bacteria/important-additional-information/expected-phenotypes/)
#' - EUCAST Breakpoint tables for interpretation of MICs and zone diameters. [(link)](https://www.eucast.org/bacteria/clinical-breakpoints-and-interpretation/clinical-breakpoint-tables/)
#' @inheritSection AMR Download Our Reference Data
#' @examples
#' \donttest{
#' a <- data.frame(
#'   mo = c(
#'     "Staphylococcus aureus",
#'     "Enterococcus faecalis",
#'     "Escherichia coli",
#'     "Klebsiella pneumoniae",
#'     "Pseudomonas aeruginosa"
#'   ),
#'   VAN = "-", # Vancomycin
#'   AMX = "-", # Amoxicillin
#'   COL = "-", # Colistin
#'   CAZ = "-", # Ceftazidime
#'   CXM = "-", # Cefuroxime
#'   PEN = "S", # Benzylpenicillin
#'   FOX = "S", # Cefoxitin
#'   stringsAsFactors = FALSE
#' )
#'
#' head(a)
#'
#'
#' # apply EUCAST rules: some results wil be changed
#' b <- eucast_rules(a, overwrite = TRUE)
#'
#' head(b)
#'
#'
#' # do not apply EUCAST rules, but rather get a data.frame
#' # containing all details about the transformations:
#' c <- eucast_rules(a, overwrite = TRUE, verbose = TRUE)
#' head(c)
#' }
#'
#' # Dosage guidelines:
#'
#' eucast_dosage(c("tobra", "genta", "cipro"), "iv")
#'
#' eucast_dosage(c("tobra", "genta", "cipro"), "iv", version_breakpoints = 10)
interpretive_rules <- function(x,
                               col_mo = NULL,
                               info = interactive(),
                               rules = getOption("AMR_interpretive_rules", default = c("breakpoints", "expected_phenotypes")),
                               guideline = getOption("AMR_guideline", "EUCAST"),
                               verbose = FALSE,
                               version_breakpoints = 16.0,
                               version_expected_phenotypes = 1.2,
                               version_expertrules = 3.3,
                               ampc_cephalosporin_resistance = NA,
                               only_sir_columns = any(is.sir(x)),
                               custom_rules = NULL,
                               overwrite = FALSE,
                               add_if_missing = TRUE,
                               ...) {
  meet_criteria(x, allow_class = "data.frame")
  meet_criteria(col_mo, allow_class = "character", has_length = 1, is_in = colnames(x), allow_NULL = TRUE)
  if (guideline %like% "EUCAST") {
    guideline <- "EUCAST"
  } else if (guideline %like% "CLSI") {
    guideline <- "CLSI"
  }
  meet_criteria(guideline, allow_class = "character", has_length = 1, is_in = c("EUCAST", "CLSI"))
  meet_criteria(info, allow_class = "logical", has_length = 1)
  meet_criteria(rules, allow_class = "character", has_length = c(1, 2, 3, 4, 5, 6), is_in = c("breakpoints", "expected_phenotypes", "expert", "other", "all", "custom"))
  meet_criteria(verbose, allow_class = "logical", has_length = 1)
  meet_criteria(version_breakpoints, allow_class = c("numeric", "integer"), has_length = 1, is_in = as.double(names(EUCAST_VERSION_BREAKPOINTS)))
  meet_criteria(version_expected_phenotypes, allow_class = c("numeric", "integer"), has_length = 1, is_in = as.double(names(EUCAST_VERSION_EXPECTED_PHENOTYPES)))
  meet_criteria(version_expertrules, allow_class = c("numeric", "integer"), has_length = 1, is_in = as.double(names(EUCAST_VERSION_EXPERT_RULES)))
  meet_criteria(ampc_cephalosporin_resistance, allow_class = c("logical", "character", "sir"), has_length = 1, allow_NA = TRUE, allow_NULL = TRUE)
  meet_criteria(only_sir_columns, allow_class = "logical", has_length = 1)
  meet_criteria(custom_rules, allow_class = c("custom_interpretive_rules", "custom_eucast_rules"), allow_NULL = TRUE)
  meet_criteria(overwrite, allow_class = "logical", has_length = 1)
  meet_criteria(add_if_missing, allow_class = "logical", has_length = 1)

  stop_if(
    !overwrite && !add_if_missing,
    "Either set {.arg overwrite} or {.arg add_if_missing} to {.code TRUE}, or both."
  )

  stop_if(
    !is.na(ampc_cephalosporin_resistance) && !any(c("expert", "all") %in% rules),
    "For the {.arg ampc_cephalosporin_resistance} argument to work, the {.arg rules} argument must contain {.code \"expert\"} or {.code \"all\"}."
  )

  add_MO_lookup_to_AMR_env()

  if ("custom" %in% rules && is.null(custom_rules)) {
    warning_("in {.help [{.fun interpretive_rules}](AMR::interpretive_rules)}: no custom rules were set with the {.arg custom_rules} argument",
      immediate = TRUE
    )
    rules <- rules[rules != "custom"]
    if (length(rules) == 0) {
      if (isTRUE(info)) {
        message_("No other rules were set, returning original data", as_note = FALSE)
      }
      return(x)
    }
  }

  x_deparsed <- deparse(substitute(x))
  if (length(x_deparsed) > 1 || any(x_deparsed %unlike% "[a-z]+")) {
    x_deparsed <- "your_data"
  }

  breakpoints_info <- EUCAST_VERSION_BREAKPOINTS[[which(as.double(names(EUCAST_VERSION_BREAKPOINTS)) == version_breakpoints)]]
  expected_phenotypes_info <- EUCAST_VERSION_EXPECTED_PHENOTYPES[[which(as.double(names(EUCAST_VERSION_EXPECTED_PHENOTYPES)) == version_expected_phenotypes)]]
  expertrules_info <- EUCAST_VERSION_EXPERT_RULES[[which(as.double(names(EUCAST_VERSION_EXPERT_RULES)) == version_expertrules)]]

  if (interactive() && isTRUE(verbose) && isTRUE(info)) {
    txt <- paste0(
      "WARNING: In Verbose mode, the interpretive_rules() function does not apply rules to the data, but instead returns a data set in logbook form with comprehensive info about which rows and columns would be effected and in which way.",
      "\n\nThis may overwrite your existing data if you use e.g.:",
      "\ndata <- interpretive_rules(data, verbose = TRUE)\n\nDo you want to continue?"
    )
    showQuestion <- import_fn("showQuestion", "rstudioapi", error_on_fail = FALSE)
    if (!is.null(showQuestion)) {
      q_continue <- showQuestion("Using verbose = TRUE with interpretive_rules()", txt)
    } else {
      q_continue <- utils::menu(choices = c("OK", "Cancel"), graphics = FALSE, title = txt)
    }
    if (q_continue %in% c(FALSE, 2)) {
      message_("Cancelled, returning original data", as_note = FALSE)
      return(x)
    }
  }

  # try to find columns based on type
  # -- mo
  if (is.null(col_mo)) {
    col_mo <- search_type_in_df(x = x, type = "mo", info = info)
    stop_if(is.null(col_mo), "{.arg col_mo} must be set")
  }

  decimal.mark <- getOption("OutDec")
  big.mark <- ifelse(decimal.mark != ",", ",", " ")
  formatnr <- function(x, big = big.mark, dec = decimal.mark) {
    trimws(format(x, big.mark = big, decimal.mark = dec))
  }

  warned <- FALSE
  warn_lacking_sir_class <- character(0)
  overwrite_skipped <- FALSE
  txt_ok <- function(n_added, n_changed, warned = FALSE) {
    if (warned == FALSE) {
      if (n_added + n_changed == 0) {
        cat(" ", font_subtle("(no changes)"), "\n", sep = "")
      } else {
        status <- character(0)
        if (n_added > 0) {
          status <- c(status, font_bold(font_green(formatnr(n_added), ifelse(n_added == 1, "value added", "values added"))))
        }
        if (n_changed > 0) {
          status <- c(status, font_bold(font_blue(formatnr(n_changed), ifelse(n_changed == 1, "value changed", "values changed"))))
        }
        cat(" ", font_subtle("("), paste0(status, collapse = font_subtle(", ")), font_subtle(")"), "\n", sep = "")
      }
      warned <<- FALSE
    }
  }

  cols_ab <- get_column_abx(
    x = x,
    soft_dependencies = c(
      "AMC",
      "AMP",
      "AMX",
      "CIP",
      "ERY",
      "FOX-S",
      "GEN",
      "MFX",
      "NAL-S",
      "NOR-S",
      "OXA-S",
      "PEN-S",
      "PIP",
      "TCY",
      "TIC",
      "TOB"
    ),
    hard_dependencies = NULL,
    verbose = verbose,
    info = info,
    only_sir_columns = only_sir_columns,
    fn = "interpretive_rules",
    ...
  )

  if (!"AMP" %in% names(cols_ab) && "AMX" %in% names(cols_ab)) {
    # ampicillin column is missing, but amoxicillin is available
    if (isTRUE(info)) {
      message_("Using column {.field ", font_bold(cols_ab[names(cols_ab) == "AMX"]), "} as input for ampicillin since many EUCAST rules depend on it.")
    }
    cols_ab <- c(cols_ab, c(AMP = unname(cols_ab[names(cols_ab) == "AMX"])))
  }

  # data preparation ----
  if (isTRUE(info) && NROW(x) > 10000) {
    message_("Preparing data...", appendLF = FALSE, as_note = FALSE)
  }

  # Some helper functions ---------------------------------------------------
  get_antibiotic_names <- function(x) {
    x <- x %pm>%
      strsplit(",") %pm>%
      unlist() %pm>%
      trimws2() %pm>%
      vapply(FUN.VALUE = character(1), function(x) if (x %in% AMR::antimicrobials$ab) ab_name(x, language = NULL, tolower = TRUE, fast_mode = TRUE) else x) %pm>%
      sort() %pm>%
      paste(collapse = ", ")
    x <- gsub("_", " ", x, fixed = TRUE)
    x <- gsub("except CAZ", paste("except", ab_name("CAZ", language = NULL, tolower = TRUE)), x, fixed = TRUE)
    x <- gsub("except TGC", paste("except", ab_name("TGC", language = NULL, tolower = TRUE)), x, fixed = TRUE)
    x <- gsub("cephalosporins (1st|2nd|3rd|4th|5th)", "cephalosporins (\\1 gen.)", x)
    x
  }
  format_antibiotic_names <- function(ab_names, ab_results) {
    ab_names <- trimws2(unlist(strsplit(ab_names, ",")))
    ab_results <- trimws2(unlist(strsplit(ab_results, ",")))
    if (length(ab_results) == 1) {
      if (length(ab_names) == 1) {
        # like FOX S
        x <- paste(ab_names, "is")
      } else if (length(ab_names) == 2) {
        # like PEN,FOX S
        x <- paste(paste0(ab_names, collapse = " and "), "are both")
      } else {
        # like PEN,FOX,GEN S
        x <- paste(paste0(ab_names, collapse = " and "), "are all")
      }
      return(paste0(x, " '", ab_results, "'"))
    } else {
      if (length(ab_names) == 2) {
        # like PEN,FOX S,R
        paste0(
          ab_names[1], " is '", ab_results[1], "' and ",
          ab_names[2], " is '", ab_results[2], "'"
        )
      } else {
        # like PEN,FOX,GEN S,R,R
        paste0(
          ab_names[1], " is '", ab_results[1], "' and ",
          ab_names[2], " is '", ab_results[2], "' and ",
          ab_names[3], " is '", ab_results[3], "'"
        )
      }
    }
  }
  like_unique <- function(vals, pattern) {
    # %like% on the unique values only, as organism columns contain many duplicates
    vals_unique <- unique(vals)
    (vals_unique %like% pattern)[match(vals, vals_unique)]
  }
  sir_values <- function(col) {
    # SIR values of a column of the current data as character, to compare without the overhead of Ops.sir()
    vals <- x[[col]]
    if (!is.sir(vals)) {
      vals <- suppressWarnings(as.sir(vals))
    }
    as.character(vals)
  }
  expand_groups <- function(entry) {
    parts <- trimws(strsplit(entry, ",")[[1]])
    group_names <- tolower(AMR::microorganisms.groups$mo_group_name)
    mo_names <- AMR::microorganisms.groups$mo_name
    group_names_lc <- tolower(group_names)
    result <- unlist(lapply(parts, function(part) {
      match_idx <- which(group_names_lc == tolower(part))
      if (length(match_idx) > 0) {
        mo_names[match_idx]
      } else {
        part
      }
    }))
    # only the ones with genus or genus/species, not subspecies (as genus_species will be matched)
    spaces <- vapply(FUN.VALUE = integer(1), strsplit(result, " "), length)
    result <- result[spaces < 3]
    return(paste0(unique(result), collapse = ", "))
  }
  # the accepted taxa of literal organism names in a rule, e.g. Thomasclavelia ramosa for 'Clostridium ramosum'
  add_MO_lookup_to_AMR_env()
  fullnames_lower <- tolower(AMR_env$MO_lookup$fullname)
  names_to_accepted_mo <- function(mo_names) {
    mo <- AMR_env$MO_lookup$mo[which(fullnames_lower %in% tolower(mo_names))]
    synonym_mo_to_accepted_mo(mo, fill_in_accepted = TRUE)
  }

  # Preparing the data ------------------------------------------------------

  verbose_info <- data.frame(
    rowid = character(0),
    col = character(0),
    mo_fullname = character(0),
    old = as.sir(character(0)),
    new = as.sir(character(0)),
    rule = character(0),
    rule_group = character(0),
    rule_name = character(0),
    rule_source = character(0),
    stringsAsFactors = FALSE
  )

  old_cols <- colnames(x)
  old_attributes <- attributes(x)
  x <- as.data.frame(x, stringsAsFactors = FALSE) # no tibbles, data.tables, etc.
  rownames(x) <- NULL # will later be restored with old_attributes
  # create unique row IDs - combination of the MO and all ABx columns (so they will only run once per unique combination)
  # (with a separator, so that e.g. "S" + "IR" and "SI" + "R" cannot get the same ID)
  x$`.rowid` <- do.call(paste, c(
    lapply(c(col_mo, cols_ab), function(col) {
      vals <- as.character(x[[col]])
      vals[is.na(vals)] <- "."
      vals
    }),
    sep = "\r"
  ))

  # save original table, with the new .rowid column
  x.bak <- x
  # keep only unique rows for MO and ABx
  x <- x %pm>%
    pm_arrange(`.rowid`) %pm>%
    # big speed gain! only analyse unique rows:
    pm_distinct(`.rowid`, .keep_all = TRUE) %pm>%
    as.data.frame(stringsAsFactors = FALSE)
  x[, col_mo] <- as.mo(as.character(x[, col_mo, drop = TRUE]), info = FALSE)
  # rename col_mo to prevent interference with joined columns
  colnames(x)[colnames(x) == col_mo] <- ".col_mo"
  col_mo <- ".col_mo"
  # join to microorganisms data set
  x <- left_join_microorganisms(x, by = col_mo, suffix = c("_oldcols", ""))
  x$gramstain <- mo_gramstain(x[, col_mo, drop = TRUE], language = NULL, info = FALSE)
  x$genus_species <- trimws(paste(x$genus, x$species))
  # the accepted name of each organism, so that a rule can name it by any of its synonyms (see 'Get rule from file')
  x$`.mo_accepted` <- synonym_mo_to_accepted_mo(as.character(x[, col_mo, drop = TRUE]), fill_in_accepted = TRUE)
  if (isTRUE(info) && NROW(x.bak) > 10000) {
    message_("OK.", as_note = FALSE)
  }

  n_added <- 0
  n_changed <- 0

  rule_current <- ""
  rule_group_current <- ""
  rule_group_previous <- ""
  rule_next <- ""
  rule_previous <- ""
  rule_text <- ""

  # >>> Apply Other rules: enzyme inhibitors <<< ------------------------------------------
  if (any(c("all", "other") %in% rules)) {
    if (isTRUE(info)) {
      cat_rules_header(
        title = paste0(
          "Rules by the ", font_bold(paste0("AMR package v", utils::packageDescription("AMR")$Version)),
          " (", format(as.Date(utils::packageDescription("AMR")$Date), format = "%Y"), ")"
        ),
        subtitle = format_inline_("See {.help [{.fun interpretive_rules}](AMR::interpretive_rules)}")
      )
    }
    ab_enzyme <- subset(AMR::antimicrobials, name %like% "/")[, c("ab", "name"), drop = FALSE]
    colnames(ab_enzyme) <- c("enzyme_ab", "enzyme_name")
    ab_enzyme$base_name <- gsub("^([a-zA-Z0-9]+).*", "\\1", ab_enzyme$enzyme_name)
    ab_enzyme$base_ab <- AMR::antimicrobials[match(ab_enzyme$base_name, AMR::antimicrobials$name), "ab", drop = TRUE]
    ab_enzyme <- subset(ab_enzyme, !is.na(base_ab))
    # make ampicillin and amoxicillin interchangable
    ampi <- subset(ab_enzyme, base_ab == "AMX")
    ampi$base_ab <- "AMP"
    ampi$base_name <- ab_name("AMP", language = NULL)
    amox <- subset(ab_enzyme, base_ab == "AMP")
    amox$base_ab <- "AMX"
    amox$base_name <- ab_name("AMX", language = NULL)
    # merge and sort
    ab_enzyme <- rbind_AMR(ab_enzyme, ampi, amox)
    ab_enzyme <- ab_enzyme[order(ab_enzyme$enzyme_name), , drop = FALSE]

    for (i in seq_len(nrow(ab_enzyme))) {
      # check if both base and base + enzyme inhibitor are part of the data set
      if (all(c(ab_enzyme$base_ab[i], ab_enzyme$enzyme_ab[i]) %in% names(cols_ab), na.rm = TRUE)) {
        col_base <- unname(cols_ab[ab_enzyme$base_ab[i]])
        col_enzyme <- unname(cols_ab[ab_enzyme$enzyme_ab[i]])

        ## Set base to R where base + enzyme inhibitor is R ----
        rule_current <- paste0(
          ab_enzyme$base_name[i], " (", col_base, ") = R if ",
          tolower(ab_enzyme$enzyme_name[i]), " (", col_enzyme, ") = R"
        )
        if (isTRUE(info)) {
          cat_rule_name(paste0(
            ab_enzyme$base_name[i], " (", format_field(col_base), ") = ", format_sir_badge("R", guideline = guideline), " if ",
            tolower(ab_enzyme$enzyme_name[i]), " (", format_field(col_enzyme), ") = ", format_sir_badge("R", guideline = guideline)
          ))
        }
        run_changes <- edit_sir(
          x = x,
          to = "R",
          rule = c(
            rule_current, "Other rules", "",
            paste0("Non-EUCAST: AMR package v", utils::packageDescription("AMR")$Version)
          ),
          rows = which(sir_values(col_enzyme) == "R"),
          cols = col_base,
          last_verbose_info = verbose_info,
          original_data = x.bak,
          warned = warned,
          info = info,
          verbose = verbose,
          overwrite = overwrite,
          add_if_missing = add_if_missing
        )
        n_added <- n_added + run_changes$added
        n_changed <- n_changed + run_changes$changed
        verbose_info <- run_changes$verbose_info
        x <- run_changes$output
        warn_lacking_sir_class <- c(warn_lacking_sir_class, run_changes$sir_warn)
        overwrite_skipped <- overwrite_skipped || isTRUE(run_changes$overwrite_skipped)
        # Print number of new changes
        if (isTRUE(info)) {
          # print only on last one of rules in this group
          txt_ok(n_added = n_added, n_changed = n_changed, warned = warned)
          # and reset counters
          n_added <- 0
          n_changed <- 0
        }

        ## Set base + enzyme inhibitor to S where base is S ----
        rule_current <- paste0(
          ab_enzyme$enzyme_name[i], " (", col_enzyme, ") = S if ",
          tolower(ab_enzyme$base_name[i]), " (", col_base, ") = S"
        )
        if (isTRUE(info)) {
          cat_rule_name(paste0(
            ab_enzyme$enzyme_name[i], " (", format_field(col_enzyme), ") = ", format_sir_badge("S", guideline = guideline), " if ",
            tolower(ab_enzyme$base_name[i]), " (", format_field(col_base), ") = ", format_sir_badge("S", guideline = guideline)
          ))
        }
        run_changes <- edit_sir(
          x = x,
          to = "S",
          rule = c(
            rule_current, "Other rules", "",
            paste0("Non-EUCAST: AMR package v", utils::packageDescription("AMR")$Version)
          ),
          rows = which(sir_values(col_base) == "S"),
          cols = col_enzyme,
          last_verbose_info = verbose_info,
          original_data = x.bak,
          warned = warned,
          info = info,
          verbose = verbose,
          overwrite = overwrite,
          add_if_missing = add_if_missing
        )
        n_added <- n_added + run_changes$added
        n_changed <- n_changed + run_changes$changed
        verbose_info <- run_changes$verbose_info
        x <- run_changes$output
        warn_lacking_sir_class <- c(warn_lacking_sir_class, run_changes$sir_warn)
        overwrite_skipped <- overwrite_skipped || isTRUE(run_changes$overwrite_skipped)
        # Print number of new changes
        if (isTRUE(info)) {
          # print only on last one of rules in this group
          txt_ok(n_added = n_added, n_changed = n_changed, warned = warned)
          # and reset counters
          n_added <- 0
          n_changed <- 0
        }
      }
    }
  } else {
    if (isTRUE(info)) {
      cat("\n")
      message_("Skipping inhibitor-inheritance rules defined by this AMR package: setting S to drug+inhibitor where drug is S, and setting R to drug where drug+inhibitor is R. Add {.val other} or {.val all} to the {.arg rules} argument to apply those rules.")
    }
  }

  if (!any(c("all", "custom") %in% rules) && !is.null(custom_rules)) {
    if (isTRUE(info)) {
      message_("Skipping custom interpretive rules, since the {.arg rules} argument does not contain {.code \"custom\"}.")
    }
    custom_rules <- NULL
  }

  # >>> Apply Official interpretive rules <<< ---------------------------------------------------
  eucast_notification_shown <- FALSE
  if (!is.null(list(...)$interpretive_rules_df)) {
    # this allows: interpretive_rules(x, interpretive_rules_df = AMR:::INTERPRETIVE_RULES_DF |> filter(is.na(have_these_values)))
    interpretive_rules_df_total <- list(...)$interpretive_rules_df
  } else if (!is.null(list(...)$eucast_rules_df)) {
    # deprecated parameter name kept for backward compatibility
    interpretive_rules_df_total <- list(...)$eucast_rules_df
    warning("Used interpretive_rules(x, eucast_rules_df = ...) - Do use newer argument interpretive_rules_df now.")
  } else {
    # internal data file, created in data-raw/_pre_commit_checks.R
    interpretive_rules_df_total <- INTERPRETIVE_RULES_DF
  }

  ## filter on guideline provider and user-set guideline versions ----
  interpretive_rules_df <- data.frame()
  if (any(c("all", "breakpoints") %in% rules)) {
    interpretive_rules_df <- interpretive_rules_df %pm>%
      rbind_AMR(interpretive_rules_df_total %pm>%
        subset(rule.provider == guideline & reference.rule_group %like% "breakpoint" & reference.version == version_breakpoints))
  }
  if (any(c("all", "expected_phenotypes") %in% rules)) {
    interpretive_rules_df <- interpretive_rules_df %pm>%
      rbind_AMR(interpretive_rules_df_total %pm>%
        subset(rule.provider == guideline & reference.rule_group %like% "expected" & reference.version == version_expected_phenotypes))
  }
  if (any(c("all", "expert") %in% rules)) {
    interpretive_rules_df <- interpretive_rules_df %pm>%
      rbind_AMR(interpretive_rules_df_total %pm>%
        subset(rule.provider == guideline & reference.rule_group %like% "expert" & reference.version == version_expertrules))
  }
  ## filter out AmpC de-repressed cephalosporin-resistant mutants ----
  # no need to filter on version number here - the rules contain these version number, so are inherently filtered
  # cefotaxime, ceftriaxone, ceftazidime
  if (is.null(ampc_cephalosporin_resistance) || isFALSE(ampc_cephalosporin_resistance)) {
    interpretive_rules_df <- subset(
      interpretive_rules_df,
      reference.rule %unlike% "ampc"
    )
  } else {
    if (isTRUE(ampc_cephalosporin_resistance)) {
      ampc_cephalosporin_resistance <- "R"
    }
    # old R versions (at least 3.0) give an error when assigning to no rows of a data.frame without rows
    ampc_rows <- if (is.null(interpretive_rules_df$reference.rule)) integer(0) else which(interpretive_rules_df$reference.rule %like% "ampc")
    if (length(ampc_rows) > 0) {
      interpretive_rules_df[ampc_rows, "to_value"] <- as.character(ampc_cephalosporin_resistance)
    }
  }

  # sometimes, the screenings are missing but the names are actually available
  # we only hints on remaining rows in `interpretive_rules_df`
  screening_abx <- as.character(AMR::antimicrobials$ab[which(AMR::antimicrobials$ab %like% "-S$")])
  screening_abx <- screening_abx[screening_abx %in% unique(unlist(strsplit(interpretive_rules_df_total$and_these_antibiotics[!is.na(interpretive_rules_df_total$and_these_antibiotics)], ", *")))]
  for (ab_s in screening_abx) {
    ab <- gsub("-S$", "", ab_s)
    if (ab %in% names(cols_ab) && !ab_s %in% names(cols_ab)) {
      if (isTRUE(info)) {
        message_(
          "Using column {.field ", font_bold(cols_ab[names(cols_ab) == ab]),
          "} as ", ab_name(ab_s, language = NULL, tolower = TRUE),
          " since a column {.code ", ab_s, "} is missing but required for the chosen rules"
        )
      }
      cols_ab <- c(cols_ab, stats::setNames(unname(cols_ab[names(cols_ab) == ab]), ab_s))
    }
  }

  ## Go over all rules and apply them ----
  mo_match_cache <- new.env(parent = emptyenv())
  # the type of each rule group, determined once instead of for every rule
  rule_groups <- as.character(interpretive_rules_df$reference.rule_group)
  rule_group_types <- rep("", length(rule_groups))
  rule_group_types[rule_groups %like% "expert"] <- "expert"
  rule_group_types[rule_groups %like% "expected"] <- "expected"
  rule_group_types[rule_groups %like% "breakpoint"] <- "breakpoint"
  rule_group_is_other <- rule_groups %like% "other"
  for (i in seq_len(nrow(interpretive_rules_df))) {
    rule_previous <- interpretive_rules_df[max(1, i - 1), "reference.rule", drop = TRUE]
    rule_current <- interpretive_rules_df[i, "reference.rule", drop = TRUE]
    rule_next <- interpretive_rules_df[min(nrow(interpretive_rules_df), i + 1), "reference.rule", drop = TRUE]
    rule_group_previous <- interpretive_rules_df[max(1, i - 1), "reference.rule_group", drop = TRUE]
    rule_group_current <- interpretive_rules_df[i, "reference.rule_group", drop = TRUE]
    # don't apply rules if user doesn't want to apply them
    if (rule_group_types[i] == "breakpoint" && !any(c("all", "breakpoints") %in% rules)) {
      next
    }
    if (rule_group_types[i] == "expected" && !any(c("all", "expected_phenotypes") %in% rules)) {
      next
    }
    if (rule_group_types[i] == "expert" && !any(c("all", "expert") %in% rules)) {
      next
    }

    if (isFALSE(info) || isFALSE(verbose)) {
      rule_text <- ""
    } else {
      if (is.na(interpretive_rules_df[i, "and_these_antibiotics", drop = TRUE])) {
        rule_text <- paste0("always report as '", interpretive_rules_df[i, "to_value", drop = TRUE], "': ", get_antibiotic_names(interpretive_rules_df[i, "then_change_these_antibiotics", drop = TRUE]))
      } else {
        rule_text <- paste0(
          "report as '", interpretive_rules_df[i, "to_value", drop = TRUE], "' when ",
          format_antibiotic_names(
            ab_names = get_antibiotic_names(interpretive_rules_df[i, "and_these_antibiotics", drop = TRUE]),
            ab_results = interpretive_rules_df[i, "have_these_values", drop = TRUE]
          ), ": ",
          get_antibiotic_names(interpretive_rules_df[i, "then_change_these_antibiotics", drop = TRUE])
        )
      }
    }

    if (i == 1) {
      rule_previous <- ""
      rule_group_previous <- ""
    }
    if (i == nrow(interpretive_rules_df)) {
      rule_next <- ""
    }

    if (isTRUE(info)) {
      ## Print EUCAST intro ------------------------------------------------------
      if (!rule_group_is_other[i] && eucast_notification_shown == FALSE) {
        cat_rules_header(
          title = paste0("Rules by ", font_bold("EUCAST")),
          subtitle = format_inline_("European Committee on Antimicrobial Susceptibility Testing, {.url https://eucast.org/}")
        )
        eucast_notification_shown <- TRUE
      }

      ## Print rule (group) ------------------------------------------------------
      if (rule_group_current != rule_group_previous) {
        # is new rule group, one of Breakpoints, Expected Phenotypes and Expert Rules
        if (rule_group_types[i] == "breakpoint") {
          group_info <- breakpoints_info
        } else if (rule_group_types[i] == "expected") {
          group_info <- expected_phenotypes_info
        } else if (rule_group_types[i] == "expert") {
          group_info <- expertrules_info
        } else {
          group_info <- NULL
        }
        if (!is.null(group_info)) {
          # the first group directly follows the EUCAST header, which already ends with an empty line
          cat(if (i > 1) "\n", wrap_rules_text(paste0(font_bold(group_info$title), " ", font_red(paste0("(", group_info$version_txt, ", ", group_info$year, ")")))), "\n", sep = "")
        }
      }
      ## Print rule  -------------------------------------------------------------
      if (rule_current != rule_previous) {
        # is new rule within group, print its name
        cat_rule_name(italicise_taxonomy(rule_current, type = "ansi"))
        warned <- FALSE
      }
    }

    ## Get rule from file ------------------------------------------------------
    if_mo_property <- trimws(interpretive_rules_df[i, "if_mo_property", drop = TRUE])
    like_is_one_of <- trimws(interpretive_rules_df[i, "like.is.one_of", drop = TRUE])
    mo_value <- trimws(interpretive_rules_df[i, "this_value", drop = TRUE])
    source_antibiotics <- interpretive_rules_df[i, "and_these_antibiotics", drop = TRUE]
    source_value <- trimws(unlist(strsplit(interpretive_rules_df[i, "have_these_values", drop = TRUE], ",", fixed = TRUE)))
    target_antibiotics <- interpretive_rules_df[i, "then_change_these_antibiotics", drop = TRUE]
    target_value <- interpretive_rules_df[i, "to_value", drop = TRUE]

    # the organism columns do not change while applying rules, so rules for the same organisms can reuse the matches
    mo_key <- paste(if_mo_property, like_is_one_of, mo_value, sep = "\r")
    mo_match <- mo_match_cache[[mo_key]]
    if (is.null(mo_match)) {
      # if amo_value contains a group name, expand that name with all species in it
      if (any(trimws(strsplit(mo_value, ",")[[1]]) %in% AMR::microorganisms.groups$mo_group_name, na.rm = TRUE)) {
        like_is_one_of <- "one_of"
        mo_value <- expand_groups(mo_value)
      }

      # Organisms match by the names in the data. Species-level rules with literal names also match by the accepted
      # taxon, so that the rule and the data may use different names of it (e.g. 'Clostridium ramosum' in a rule, while
      # its accepted name is Thomasclavelia ramosa). Not for genera, as a former genus also contains species that now
      # belong to other genera (e.g. Pseudomonas maltophilia), and not for regular expressions, as these would then also
      # match synonyms they were written to exclude (e.g. Klebsiella mobilis for '^Klebsiella(?! aerogenes)').
      accepted_mo_of_rule <- NULL
      if (like_is_one_of %in% c("is", "one_of") && if_mo_property %in% c("genus_species", "fullname")) {
        accepted_mo_of_rule <- names_to_accepted_mo(trimws(strsplit(mo_value, ",", fixed = TRUE)[[1]]))
      }

      if (like_is_one_of == "is") {
        # so e.g. 'Enterococcus' will turn into '^Enterococcus$'
        mo_value <- paste0("^", mo_value, "$")
      } else if (like_is_one_of == "one_of") {
        # so 'Clostridium, Actinomyces, ...' will turn into '^(Clostridium|Actinomyces|...)$'
        mo_value <- paste0(
          "^(",
          paste(trimws(unlist(strsplit(mo_value, ",", fixed = TRUE))),
            collapse = "|"
          ),
          ")$"
        )
      } else if (like_is_one_of != "like") {
        stop("invalid value for column {.field like.is.one_of}", call. = FALSE)
      }

      mo_match <- tryCatch(
        like_unique(x[, if_mo_property, drop = TRUE], mo_value) | like_unique(x$fullname, mo_value),
        error = function(e) rep(FALSE, nrow(x))
      )
      if (!is.null(accepted_mo_of_rule)) {
        mo_match <- mo_match | x$`.mo_accepted` %in% accepted_mo_of_rule
      }
      mo_match_cache[[mo_key]] <- mo_match
    }

    if (is.na(source_antibiotics)) {
      rows <- which(mo_match)
    } else {
      source_antibiotics <- get_ab_from_namespace(source_antibiotics, cols_ab)
      if (length(source_value) == 1 && length(source_antibiotics) > 1) {
        source_value <- rep(source_value, length(source_antibiotics))
      }
      if (length(source_antibiotics) == 0) {
        rows <- integer(0)
      } else if (length(source_antibiotics) == 1) {
        rows <- tryCatch(
          which(mo_match &
            sir_values(source_antibiotics[1L]) == source_value[1L]),
          error = function(e) integer(0)
        )
      } else if (length(source_antibiotics) == 2) {
        rows <- tryCatch(
          which(mo_match &
            sir_values(source_antibiotics[1L]) == source_value[1L] &
            sir_values(source_antibiotics[2L]) == source_value[2L]),
          error = function(e) integer(0)
        )
      } else {
        stop_("only 2 antimicrobials supported for source_antibiotics")
      }
    }

    cols <- get_ab_from_namespace(target_antibiotics, cols_ab)

    ## Apply rule on data ------------------------------------------------------
    # this will return the unique number of changes
    run_changes <- edit_sir(
      x = x,
      to = target_value,
      rule = c(
        rule_text, rule_group_current, rule_current,
        ifelse(rule_group_types[i] == "breakpoint",
          paste0(breakpoints_info$title, " ", breakpoints_info$version_txt, ", ", breakpoints_info$year),
          ifelse(rule_group_types[i] == "expected",
            paste0(expected_phenotypes_info$title, " ", expected_phenotypes_info$version_txt, ", ", expected_phenotypes_info$year),
            paste0(expertrules_info$title, " ", expertrules_info$version_txt, ", ", expertrules_info$year)
          )
        )
      ),
      rows = rows,
      cols = cols,
      last_verbose_info = verbose_info,
      original_data = x.bak,
      warned = warned,
      info = info,
      verbose = verbose,
      overwrite = overwrite,
      add_if_missing = add_if_missing
    )
    n_added <- n_added + run_changes$added
    n_changed <- n_changed + run_changes$changed
    verbose_info <- run_changes$verbose_info
    x <- run_changes$output
    warn_lacking_sir_class <- c(warn_lacking_sir_class, run_changes$sir_warn)
    overwrite_skipped <- overwrite_skipped || isTRUE(run_changes$overwrite_skipped)
    ## Print number of new changes ---------------------------------------------
    if (isTRUE(info) && rule_next != rule_current) {
      # print only on last one of rules in this group
      txt_ok(n_added = n_added, n_changed = n_changed, warned = warned)
      # and reset counters
      n_added <- 0
      n_changed <- 0
    }
  } # end of going over all rules

  # >>> Apply custom rules <<< ----
  if (!is.null(custom_rules)) {
    if (isTRUE(info)) {
      cat_rules_header(title = paste0(font_bold("Custom interpretive rules"), ", set by user"))
    }
    for (i in seq_len(length(custom_rules))) {
      rule <- custom_rules[[i]]
      rows <- tryCatch(which(eval(parse(text = rule$query), envir = x)),
        error = function(e) stop_(conditionMessage(e), " (check available data and compare with the custom rules set)", call = FALSE)
      )
      cols <- as.character(rule$result_group)
      cols <- c(
        cols[cols %in% colnames(x)], # direct column names
        unname(cols_ab[names(cols_ab) %in% cols])
      ) # based on previous cols_ab finding
      cols <- unique(cols)
      target_value <- as.character(rule$result_value)
      rule_text <- paste0(
        "report as '", target_value, "' when ",
        format_custom_query_rule(rule$query, colours = FALSE), ": ",
        get_antibiotic_names(cols)
      )
      if (isTRUE(info)) {
        # formatted like print.custom_interpretive_rules()
        cat_rule_name(paste0(
          font_bold("If "), font_blue(format_custom_query_rule(rule$query, guideline = guideline)),
          font_bold(" then "), "set to ", format_sir_badge(target_value, guideline = guideline), ": ",
          italicise_taxonomy(get_antibiotic_names(cols), type = "ansi")
        ))
        warned <- FALSE
      }
      run_changes <- edit_sir(
        x = x,
        to = target_value,
        rule = c(
          rule_text,
          "Custom interpretive rules",
          paste0("Custom interpretive rule ", i),
          paste0(
            "Object '", deparse(substitute(custom_rules)),
            "' consisting of ", length(custom_rules), " custom rules"
          )
        ),
        rows = rows,
        cols = cols,
        last_verbose_info = verbose_info,
        original_data = x.bak,
        warned = warned,
        info = info,
        verbose = verbose,
        overwrite = overwrite,
        add_if_missing = add_if_missing
      )
      n_added <- n_added + run_changes$added
      n_changed <- n_changed + run_changes$changed
      verbose_info <- run_changes$verbose_info
      x <- run_changes$output
      warn_lacking_sir_class <- c(warn_lacking_sir_class, run_changes$sir_warn)
      overwrite_skipped <- overwrite_skipped || isTRUE(run_changes$overwrite_skipped)
      ## Print number of new changes ---------------------------------------------
      if (isTRUE(info)) {
        # every custom rule is its own group
        txt_ok(n_added = n_added, n_changed = n_changed, warned = warned)
        # and reset counters
        n_added <- 0
        n_changed <- 0
      }
    }
  }

  # Print overview ----------------------------------------------------------
  if (isTRUE(info) || isTRUE(verbose)) {
    verbose_info <- x.bak %pm>%
      pm_mutate(row = pm_row_number()) %pm>%
      pm_select(`.rowid`, row) %pm>%
      pm_right_join(verbose_info,
        by = c(".rowid" = "rowid")
      ) %pm>%
      pm_select(-`.rowid`) %pm>%
      pm_select(row, pm_everything()) %pm>%
      pm_filter(!is.na(new) | is.na(new) & !is.na(old)) %pm>%
      pm_arrange(row, rule_group, rule_name, col)
    if (isFALSE(overwrite)) {
      verbose_info <- verbose_info %pm>%
        pm_filter(!old %in% levels(NA_sir_))
    }
    rownames(verbose_info) <- NULL
  }

  if (isTRUE(info)) {
    if (isTRUE(verbose)) {
      wouldve <- "would have "
    } else {
      wouldve <- ""
    }

    cat_rules_header(title = font_bold("Summary"))
    cat(wrap_rules_text(paste0(
      "The rules ", wouldve, "affected ",
      font_bold(formatnr(pm_n_distinct(verbose_info$row)), "out of", formatnr(nrow(x.bak)), "rows"),
      ", making a total of ",
      font_bold(formatnr(nrow(verbose_info)), ifelse(nrow(verbose_info) == 1, "edit", "edits"))
    )), "\n", sep = "")

    added_rows <- verbose_info[which(is.na(verbose_info$old)), , drop = FALSE]
    changed_rows <- verbose_info[which(!is.na(verbose_info$old)), , drop = FALSE]
    total_n_added <- nrow(added_rows)
    total_n_changed <- nrow(changed_rows)
    arrow <- if (has_cli_rlang()) cli::symbol$arrow_right else "=>"
    test_results <- function(n) paste(formatnr(n), ifelse(n == 1, "test result", "test results"))
    count_old_new <- function(df) {
      # count combinations of old and new values, keeping NA as a value (pm_count() does not support NA)
      old <- as.character(df$old)
      new <- as.character(df$new)
      key <- paste(ifelse(is.na(old), "\001", old), ifelse(is.na(new), "\001", new), sep = "\002")
      tbl <- table(key)
      parts <- strsplit(names(tbl), "\002", fixed = TRUE)
      out <- data.frame(
        old = vapply(FUN.VALUE = character(1), parts, function(p) p[1]),
        new = vapply(FUN.VALUE = character(1), parts, function(p) p[2]),
        n = as.integer(tbl),
        stringsAsFactors = FALSE
      )
      out$old[out$old == "\001"] <- NA_character_
      out$new[out$new == "\001"] <- NA_character_
      # sort on the order of SIR levels, with NA and other values last
      out[order(match(out$old, VALID_SIR_LEVELS), match(out$new, VALID_SIR_LEVELS), out$old, out$new), , drop = FALSE]
    }

    # print added values
    cat(
      if (total_n_added == 0) arrow else font_green(arrow),
      " ", wouldve, "added ",
      if (total_n_added == 0) font_bold(test_results(0)) else font_bold(font_green(test_results(total_n_added))),
      "\n",
      sep = ""
    )
    if (total_n_added > 0) {
      added_summary <- count_old_new(added_rows)
      cat(
        paste0(
          "  ", AMR_env$bullet_icon, " ",
          formatnr(added_summary$n), " test result", ifelse(added_summary$n > 1, "s", ""),
          " ", wouldve, "added as ", format_sir_badge(added_summary$new, guideline = guideline),
          "\n"
        ),
        sep = ""
      )
    }

    # print changed values
    cat(
      if (total_n_changed == 0) arrow else font_blue(arrow),
      " ", wouldve, "changed ",
      if (total_n_changed == 0) font_bold(test_results(0)) else font_bold(font_blue(test_results(total_n_changed))),
      "\n",
      sep = ""
    )
    if (total_n_changed > 0) {
      changed_summary <- count_old_new(changed_rows)
      cat(
        paste0(
          "  ", AMR_env$bullet_icon, " ",
          formatnr(changed_summary$n), " test result", ifelse(changed_summary$n > 1, "s", ""),
          " ", wouldve, "changed from ", format_sir_badge(changed_summary$old, guideline = guideline),
          " to ", format_sir_badge(changed_summary$new, guideline = guideline),
          "\n"
        ),
        sep = ""
      )
    }

    if (isFALSE(verbose) && total_n_added + total_n_changed > 0) {
      cat("\n", wrap_rules_text(format_inline_(
        "Use {.code interpretive_rules(..., verbose = TRUE)} (on your original data) to get a data.frame with all specified edits instead."
      )), "\n\n", sep = "")
    } else if (isTRUE(verbose)) {
      cat("\n", wrap_rules_text(format_inline_(
        "Used 'Verbose mode' ({.code verbose = TRUE}), which returns a data.frame with all specified edits. Use {.code verbose = FALSE} to apply the rules on your data."
      )), "\n\n", sep = "")
    }
  }

  if (isTRUE(overwrite_skipped) && message_not_thrown_before("edit_sir.warning_overwrite")) {
    warning_("in {.help [{.fun interpretive_rules}](AMR::interpretive_rules)}: some columns had SIR values which were not overwritten, since {.code overwrite = FALSE}.")
  }

  if (length(warn_lacking_sir_class) > 0) {
    warn_lacking_sir_class <- unique(warn_lacking_sir_class)
    # take order from original data set
    warn_lacking_sir_class <- warn_lacking_sir_class[order(colnames(x.bak))]
    warn_lacking_sir_class <- warn_lacking_sir_class[!is.na(warn_lacking_sir_class)]
    warning_(
      "in {.help [{.fun interpretive_rules}](AMR::interpretive_rules)}: not all columns with antimicrobial results are of class {.cls sir}. Transform them on beforehand, e.g.:\n\n",
      "\u00a0\u00a0", AMR_env$bullet_icon, " ", highlight_code(paste0(x_deparsed, " |> as.sir(", ifelse(length(warn_lacking_sir_class) == 1,
        warn_lacking_sir_class,
        paste0(warn_lacking_sir_class[1], ":", warn_lacking_sir_class[length(warn_lacking_sir_class)])
      ), ")")), "\n\n",
      "\u00a0\u00a0", AMR_env$bullet_icon, " ", highlight_code(paste0(x_deparsed, " |> mutate_if(is_sir_eligible, as.sir)")), "\n\n",
      "\u00a0\u00a0", AMR_env$bullet_icon, " ", highlight_code(paste0(x_deparsed, " |> mutate(across(where(is_sir_eligible), as.sir))"))
    )
  }

  # Return data set ---------------------------------------------------------
  if (isTRUE(verbose)) {
    as_original_data_class(verbose_info, old_attributes$class) # will remove tibble groups
  } else {
    # x was analysed with only unique rows, so join everything together again
    x <- x[, c(cols_ab, ".rowid"), drop = FALSE]
    x.bak <- x.bak[, setdiff(colnames(x.bak), cols_ab), drop = FALSE]
    x.bak <- x.bak %pm>%
      pm_left_join(x, by = ".rowid")
    x.bak <- x.bak[, old_cols, drop = FALSE]
    # reset original attributes
    attributes(x.bak) <- old_attributes
    x.bak <- as_original_data_class(x.bak, old_class = class(x.bak)) # will remove tibble groups
    x.bak
  }
}

#' @rdname interpretive_rules
#' @export
eucast_rules <- function(x,
                         col_mo = NULL,
                         info = interactive(),
                         rules = getOption("AMR_interpretive_rules", default = c("breakpoints", "expected_phenotypes")),
                         ...) {
  if (!is.null(getOption("AMR_eucastrules", default = NULL))) {
    warning_("The global option {.code AMR_eucastrules} that you have set is no longer used and was ignored, set {.code AMR_interpretive_rules} instead. See {.topic [AMR-options](AMR::AMR-options)}.")
  }
  interpretive_rules(x = x, col_mo = col_mo, info = info, rules = rules, guideline = "EUCAST", ...)
}

#' @rdname interpretive_rules
#' @export
clsi_rules <- function(x,
                       col_mo = NULL,
                       info = interactive(),
                       rules = getOption("AMR_interpretive_rules", default = c("breakpoints", "expected_phenotypes")),
                       ...) {
  interpretive_rules(x = x, col_mo = col_mo, info = info, rules = rules, guideline = "CLSI", ...)
}

# helper function for editing the table ----
edit_sir <- function(x,
                     to,
                     rule,
                     rows,
                     cols,
                     last_verbose_info,
                     original_data,
                     warned,
                     info,
                     verbose,
                     overwrite,
                     add_if_missing) {
  cols <- unique(cols[!is.na(cols) & !is.null(cols)])
  rows <- unique(rows)

  # for Verbose Mode, keep track of all changes and return them
  track_changes <- list(
    added = 0,
    changed = 0,
    output = x,
    verbose_info = last_verbose_info,
    sir_warn = character(0),
    overwrite_skipped = FALSE
  )

  txt_error <- function() {
    if (isTRUE(info)) cat("", font_red_bg(" ERROR "), "\n\n")
  }
  txt_warning <- function() {
    if (warned == FALSE) {
      if (isTRUE(info)) cat(" ", font_orange_bg(" WARNING "), sep = "")
    }
    warned <<- TRUE
  }

  if (length(rows) > 0 && length(cols) > 0) {
    new_edits <- x
    if (any(!vapply(FUN.VALUE = logical(1), x[, cols, drop = FALSE], is.sir), na.rm = TRUE)) {
      track_changes$sir_warn <- cols[!vapply(FUN.VALUE = logical(1), x[, cols, drop = FALSE], is.sir)]
    }
    # current values as a character matrix (rows x cols), to compare without the overhead of Ops.sir()
    current <- matrix(
      unlist(lapply(cols, function(col) as.character(new_edits[[col]][rows])), use.names = FALSE),
      nrow = length(rows)
    )
    isSIR <- !is.na(current) & current %in% VALID_SIR_LEVELS
    # the warning about this is thrown once by interpretive_rules(), as message_not_thrown_before() is too slow to call for every rule
    track_changes$overwrite_skipped <- isFALSE(overwrite) && any(isSIR)
    # determine which cells to modify based on overwrite and add_if_missing
    if (isTRUE(overwrite)) {
      if (isTRUE(add_if_missing)) {
        apply_mask <- rep(TRUE, length(isSIR))
      } else {
        apply_mask <- isSIR
      }
    } else {
      # overwrite = FALSE, add_if_missing = TRUE: fill missing and placeholder cells only
      apply_mask <- !isSIR
    }
    # most rules leave the data as they are, so skip writing (and tracking) if no value would change
    if (is.na(to)) {
      will_change <- apply_mask & !is.na(current)
    } else {
      will_change <- apply_mask & (is.na(current) | current != to)
    }
    if (!any(will_change)) {
      return(track_changes)
    }

    do_assign <- function() {
      subset <- new_edits[rows, cols, drop = FALSE]
      mask <- matrix(apply_mask, nrow = nrow(subset), ncol = ncol(subset))
      subset[mask] <- to
      new_edits[rows, cols] <<- subset
    }

    tryCatch(
      do_assign(),
      warning = function(w) {
        if (w$message %like% "invalid factor level") {
          vapply(FUN.VALUE = logical(1), cols, function(col) {
            new_edits[, col] <<- factor(
              x = as.character(pm_pull(new_edits, col)),
              levels = unique(c(to, levels(pm_pull(new_edits, col))))
            )
            TRUE
          })
          suppressWarnings(do_assign())
          warning_(
            "in {.help [{.fun interpretive_rules}](AMR::interpretive_rules)}: value \"", to, "\" added to the factor levels of column",
            ifelse(length(cols) == 1, "", "s"),
            " ", vector_and(cols, quotes = "`", sort = FALSE),
            " because this value was not an existing factor level."
          )
          txt_warning()
          warned <<- FALSE
        } else {
          warning_("in {.help [{.fun interpretive_rules}](AMR::interpretive_rules)}: ", w$message)
          txt_warning()
        }
      },
      error = function(e) {
        txt_error()
        stop(
          paste0(
            "In row(s) ", paste(rows[seq_len(min(length(rows), 10))], collapse = ","),
            ifelse(length(rows) > 10, "...", ""),
            " while writing value '", to,
            "' to column(s) `", paste(cols, collapse = "`, `"),
            "`:\n", conditionMessage(e)
          ),
          call. = FALSE
        )
      }
    )

    track_changes$output <- new_edits
    if (isTRUE(info) || isTRUE(verbose)) {
      # the number of rows in the original data that the (unique) row IDs stand for
      n_original_rows <- function(rowids) {
        sum(original_data$.rowid %in% rowids)
      }
      for (i in seq_len(length(cols))) {
        verbose_new <- data.frame(
          rowid = new_edits[rows, ".rowid", drop = TRUE],
          col = cols[i],
          mo_fullname = new_edits[rows, "fullname", drop = TRUE],
          old = x[rows, cols[i], drop = TRUE],
          # the value after the edit, so that cells left as they were (due to `overwrite` or `add_if_missing`) are not counted
          new = as.character(new_edits[rows, cols[i], drop = TRUE]),
          rule = font_stripstyle(rule[1]),
          rule_group = font_stripstyle(rule[2]),
          rule_name = font_stripstyle(rule[3]),
          rule_source = font_stripstyle(rule[4]),
          stringsAsFactors = FALSE
        )
        colnames(verbose_new) <- c(
          "rowid", "col", "mo_fullname", "old", "new",
          "rule", "rule_group", "rule_name", "rule_source"
        )
        is_edited <- ifelse(
          is.na(verbose_new$old) | is.na(verbose_new$new),
          xor(is.na(verbose_new$old), is.na(verbose_new$new)),
          as.character(verbose_new$old) != verbose_new$new
        )
        verbose_new <- verbose_new[which(is_edited), , drop = FALSE]
        # save changes to data set 'verbose_info'
        track_changes$verbose_info <- rbind_AMR(
          track_changes$verbose_info,
          verbose_new
        )
        # count adds and changes
        track_changes$added <- track_changes$added + n_original_rows(verbose_new$rowid[is.na(verbose_new$old)])
        track_changes$changed <- track_changes$changed + n_original_rows(verbose_new$rowid[!is.na(verbose_new$old)])
      }
    }
  }
  return(track_changes)
}

# helper functions for printing ----
# The output of interpretive_rules() is printed to stdout with cat(). It is styled with cli where available (see has_cli_rlang()), with a plain-text fallback for old R versions or when cli is not installed.

rules_console_width <- function() {
  if (has_cli_rlang()) {
    cli::console_width()
  } else {
    getOption("width", 80)
  }
}

wrap_rules_text <- function(x, width = rules_console_width(), indent = 0, exdent = 0) {
  if (has_cli_rlang()) {
    out <- cli::ansi_strwrap(x, width = width, indent = indent, exdent = exdent)
  } else {
    # strwrap() would count ANSI styling as characters, so wrap on visible width
    words <- strsplit(x, " ", fixed = TRUE)[[1]]
    words <- words[words != ""]
    out <- character(0)
    line <- strrep(" ", indent)
    line_has_words <- FALSE
    for (word in words) {
      if (line_has_words && nchar(font_stripstyle(line)) + 1 + nchar(font_stripstyle(word)) > width) {
        out <- c(out, line)
        line <- paste0(strrep(" ", exdent), word)
      } else {
        line <- paste0(line, if (line_has_words) " ", word)
      }
      line_has_words <- TRUE
    }
    out <- c(out, line)
  }
  paste0(out, collapse = "\n")
}

cat_rules_header <- function(title, subtitle = NULL) {
  if (has_cli_rlang()) {
    header <- as.character(cli::rule(left = title))
  } else {
    header <- paste0(font_grey("--"), " ", title, " ")
    header <- paste0(header, font_grey(strrep("-", max(3, rules_console_width() - nchar(font_stripstyle(header))))))
  }
  cat("\n", header, "\n", sep = "")
  if (!is.null(subtitle)) {
    cat(wrap_rules_text(subtitle), "\n", sep = "")
  }
  cat("\n")
}

cat_rule_name <- function(x) {
  # printed without line ending, since edit_sir() and txt_ok() will print on the same line
  cat(wrap_rules_text(
    paste0(AMR_env$bullet_icon, " ", x),
    indent = 2,
    exdent = 4
  ))
}

format_field <- function(x) {
  if (has_cli_rlang()) {
    cli::format_inline("{.field {x}}")
  } else {
    x
  }
}

#' @rdname interpretive_rules
#' @export
eucast_dosage <- function(ab, administration = "iv", version_breakpoints = 15) {
  meet_criteria(ab, allow_class = c("character", "numeric", "integer", "factor"))
  meet_criteria(administration, allow_class = "character", is_in = dosage$administration[!is.na(dosage$administration)], has_length = 1)
  meet_criteria(version_breakpoints, allow_class = c("numeric", "integer"), has_length = 1, is_in = as.double(names(EUCAST_VERSION_BREAKPOINTS)))

  # show used version_breakpoints number once per session (AMR_env will reload every session)
  if (missing(version_breakpoints) && message_not_thrown_before("eucast_dosage", "v", gsub("[^0-9]", "", version_breakpoints), entire_session = TRUE)) {
    message_(
      "Dosages for antimicrobial drugs, as meant for ",
      format_eucast_version_nr(version_breakpoints, markdown = FALSE), ". ",
      font_red("This note will be shown once per session.")
    )
  }

  ab <- as.ab(ab)
  lst <- vector("list", length = length(ab))
  for (i in seq_len(length(ab))) {
    df <- AMR::dosage[which(AMR::dosage$eucast_version == version_breakpoints & AMR::dosage$ab == ab[i] & AMR::dosage$administration == administration), , drop = FALSE]
    lst[[i]] <- list(
      ab = "",
      name = "",
      standard_dosage = ifelse("standard_dosage" %in% df$type,
        trimws2(df[which(df$type == "standard_dosage"), "original_txt", drop = TRUE]),
        NA_character_
      ),
      high_dosage = ifelse("high_dosage" %in% df$type,
        trimws2(df[which(df$type == "high_dosage"), "original_txt", drop = TRUE]),
        NA_character_
      ),
      eucast_version = df$eucast_version[1]
    )
  }
  out <- do.call(rbind_AMR, lapply(lst, as.data.frame, stringsAsFactors = FALSE))
  rownames(out) <- NULL
  out$ab <- ab
  out$name <- ab_name(ab, language = NULL)
  if (pkg_is_available("tibble")) {
    import_fn("as_tibble", "tibble")(out)
  } else {
    out
  }
}
