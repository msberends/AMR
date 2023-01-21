# ==================================================================== #
# TITLE                                                                #
# AMR: An R Package for Working with Antimicrobial Resistance Data     #
#                                                                      #
# SOURCE                                                               #
# https://github.com/msberends/AMR                                     #
#                                                                      #
# CITE AS                                                              #
# Berends MS, Luz CF, Friedrich AW, Sinha BNM, Albers CJ, Glasner C    #
# (2022). AMR: An R Package for Working with Antimicrobial Resistance  #
# Data. Journal of Statistical Software, 104(3), 1-31.                 #
# doi:10.18637/jss.v104.i03                                            #
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

# ====================================================== #
# || Change the EUCAST version numbers in R/globals.R || #
# ====================================================== #

format_eucast_version_nr <- function(version, markdown = TRUE) {
  # for documentation - adds title, version number, year and url in markdown language
  lst <- c(EUCAST_VERSION_BREAKPOINTS, EUCAST_VERSION_EXPERT_RULES)
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

#' Apply EUCAST Rules
#'
#' @description
#' Apply rules for clinical breakpoints and intrinsic resistance as defined by the European Committee on Antimicrobial Susceptibility Testing (EUCAST, <https://www.eucast.org>), see *Source*. Use [eucast_dosage()] to get a [data.frame] with advised dosages of a certain bug-drug combination, which is based on the [dosage] data set.
#'
#' To improve the interpretation of the antibiogram before EUCAST rules are applied, some non-EUCAST rules can applied at default, see *Details*.
#' @param x a data set with antibiotic columns, such as `amox`, `AMX` and `AMC`
#' @param info a [logical] to indicate whether progress should be printed to the console, defaults to only print while in interactive sessions
#' @param rules a [character] vector that specifies which rules should be applied. Must be one or more of `"breakpoints"`, `"expert"`, `"other"`, `"custom"`, `"all"`, and defaults to `c("breakpoints", "expert")`. The default value can be set to another value, e.g. using `options(AMR_eucastrules = "all")`. If using `"custom"`, be sure to fill in argument `custom_rules` too. Custom rules can be created with [custom_eucast_rules()].
#' @param verbose a [logical] to turn Verbose mode on and off (default is off). In Verbose mode, the function does not apply rules to the data, but instead returns a data set in logbook form with extensive info about which rows and columns would be effected and in which way. Using Verbose mode takes a lot more time.
#' @param version_breakpoints the version number to use for the EUCAST Clinical Breakpoints guideline. Can be either `r vector_or(names(EUCAST_VERSION_BREAKPOINTS), reverse = TRUE)`.
#' @param version_expertrules the version number to use for the EUCAST Expert Rules and Intrinsic Resistance guideline. Can be either `r vector_or(names(EUCAST_VERSION_EXPERT_RULES), reverse = TRUE)`.
#' @param ampc_cephalosporin_resistance a [character] value that should be applied to cefotaxime, ceftriaxone and ceftazidime for AmpC de-repressed cephalosporin-resistant mutants, defaults to `NA`. Currently only works when `version_expertrules` is `3.2` and higher; these version of '*EUCAST Expert Rules on Enterobacterales*' state that results of cefotaxime, ceftriaxone and ceftazidime should be reported with a note, or results should be suppressed (emptied) for these three drugs. A value of `NA` (the default) for this argument will remove results for these three drugs, while e.g. a value of `"R"` will make the results for these drugs resistant. Use `NULL` or `FALSE` to not alter results for these three drugs of AmpC de-repressed cephalosporin-resistant mutants. Using `TRUE` is equal to using `"R"`. \cr For *EUCAST Expert Rules* v3.2, this rule applies to: `r vector_and(gsub("[^a-zA-Z ]+", "", unlist(strsplit(EUCAST_RULES_DF[which(EUCAST_RULES_DF$reference.version %in% c(3.2, 3.3) & EUCAST_RULES_DF$reference.rule %like% "ampc"), "this_value"][1], "|", fixed = TRUE))), quotes = "*")`.
#' @param ... column name of an antibiotic, see section *Antibiotics* below
#' @param ab any (vector of) text that can be coerced to a valid antibiotic drug code with [as.ab()]
#' @param administration route of administration, either `r vector_or(dosage$administration)`
#' @param only_sir_columns a [logical] to indicate whether only antibiotic columns must be detected that were transformed to class `sir` (see [as.sir()]) on beforehand (defaults to `FALSE`)
#' @param custom_rules custom rules to apply, created with [custom_eucast_rules()]
#' @inheritParams first_isolate
#' @details
#' **Note:** This function does not translate MIC values to SIR values. Use [as.sir()] for that. \cr
#' **Note:** When ampicillin (AMP, J01CA01) is not available but amoxicillin (AMX, J01CA04) is, the latter will be used for all rules where there is a dependency on ampicillin. These drugs are interchangeable when it comes to expression of antimicrobial resistance. \cr
#'
#' The file containing all EUCAST rules is located here: <https://github.com/msberends/AMR/blob/main/data-raw/eucast_rules.tsv>.  **Note:** Old taxonomic names are replaced with the current taxonomy where applicable. For example, *Ochrobactrum anthropi* was renamed to *Brucella anthropi* in 2020; the original EUCAST rules v3.1 and v3.2 did not yet contain this new taxonomic name. The `AMR` package contains the full microbial taxonomy updated until `r documentation_date(max(TAXONOMY_VERSION$GBIF$accessed_date, TAXONOMY_VERSION$LPSN$accessed_date))`, see [microorganisms].
#'
#' ### Custom Rules
#'
#' Custom rules can be created using [custom_eucast_rules()], e.g.:
#'
#' ```r
#' x <- custom_eucast_rules(AMC == "R" & genus == "Klebsiella" ~ aminopenicillins == "R",
#'                          AMC == "I" & genus == "Klebsiella" ~ aminopenicillins == "I")
#'
#' eucast_rules(example_isolates, rules = "custom", custom_rules = x)
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
#' Since these rules are not officially approved by EUCAST, they are not applied at default. To use these rules, include `"other"` to the `rules` argument, or use `eucast_rules(..., rules = "all")`. You can also set the option `AMR_eucastrules`, i.e. run `options(AMR_eucastrules = "all")`.
#' @section Antibiotics:
#' To define antibiotics column names, leave as it is to determine it automatically with [guess_ab_col()] or input a text (case-insensitive), or use `NULL` to skip a column (e.g. `TIC = NULL` to skip ticarcillin). Manually defined but non-existing columns will be skipped with a warning.
#'
#' The following antibiotics are eligible for the functions [eucast_rules()] and [mdro()]. These are shown below in the format 'name (`antimicrobial ID`, [ATC code](https://www.whocc.no/atc/structure_and_principles/))', sorted alphabetically:
#'
#' `r create_eucast_ab_documentation()`
#' @aliases EUCAST
#' @rdname eucast_rules
#' @export
#' @return The input of `x`, possibly with edited values of antibiotics. Or, if `verbose = TRUE`, a [data.frame] with all original and new values of the affected bug-drug combinations.
#' @source
#' - EUCAST Expert Rules. Version 2.0, 2012.\cr
#'   Leclercq et al. **EUCAST expert rules in antimicrobial susceptibility testing.** *Clin Microbiol Infect.* 2013;19(2):141-60; \doi{https://doi.org/10.1111/j.1469-0691.2011.03703.x}
#' - EUCAST Expert Rules, Intrinsic Resistance and Exceptional Phenotypes Tables. Version 3.1, 2016. [(link)](https://www.eucast.org/fileadmin/src/media/PDFs/EUCAST_files/Expert_Rules/Expert_rules_intrinsic_exceptional_V3.1.pdf)
#' - EUCAST Intrinsic Resistance and Unusual Phenotypes. Version 3.2, 2020. [(link)](https://www.eucast.org/fileadmin/src/media/PDFs/EUCAST_files/Expert_Rules/2020/Intrinsic_Resistance_and_Unusual_Phenotypes_Tables_v3.2_20200225.pdf)
#' - EUCAST Intrinsic Resistance and Unusual Phenotypes. Version 3.3, 2021. [(link)](https://www.eucast.org/fileadmin/src/media/PDFs/EUCAST_files/Expert_Rules/2021/Intrinsic_Resistance_and_Unusual_Phenotypes_Tables_v3.3_20211018.pdf)
#' - EUCAST Breakpoint tables for interpretation of MICs and zone diameters. Version 9.0, 2019. [(link)](https://www.eucast.org/fileadmin/src/media/PDFs/EUCAST_files/Breakpoint_tables/v_9.0_Breakpoint_Tables.xlsx)
#' - EUCAST Breakpoint tables for interpretation of MICs and zone diameters. Version 10.0, 2020. [(link)](https://www.eucast.org/fileadmin/src/media/PDFs/EUCAST_files/Breakpoint_tables/v_10.0_Breakpoint_Tables.xlsx)
#' - EUCAST Breakpoint tables for interpretation of MICs and zone diameters. Version 11.0, 2021. [(link)](https://www.eucast.org/fileadmin/src/media/PDFs/EUCAST_files/Breakpoint_tables/v_11.0_Breakpoint_Tables.xlsx)
#' - EUCAST Breakpoint tables for interpretation of MICs and zone diameters. Version 12.0, 2022. [(link)](https://www.eucast.org/fileadmin/src/media/PDFs/EUCAST_files/Breakpoint_tables/v_12.0_Breakpoint_Tables.xlsx)
#' @inheritSection AMR Reference Data Publicly Available
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
#' b <- eucast_rules(a)
#'
#' head(b)
#'
#'
#' # do not apply EUCAST rules, but rather get a data.frame
#' # containing all details about the transformations:
#' c <- eucast_rules(a, verbose = TRUE)
#' head(c)
#' }
#'
#' # Dosage guidelines:
#'
#' eucast_dosage(c("tobra", "genta", "cipro"), "iv")
#'
#' eucast_dosage(c("tobra", "genta", "cipro"), "iv", version_breakpoints = 10)
eucast_rules <- function(x,
                         col_mo = NULL,
                         info = interactive(),
                         rules = getOption("AMR_eucastrules", default = c("breakpoints", "expert")),
                         verbose = FALSE,
                         version_breakpoints = 12.0,
                         version_expertrules = 3.3,
                         ampc_cephalosporin_resistance = NA,
                         only_sir_columns = FALSE,
                         custom_rules = NULL,
                         ...) {
  meet_criteria(x, allow_class = "data.frame")
  meet_criteria(col_mo, allow_class = "character", has_length = 1, is_in = colnames(x), allow_NULL = TRUE)
  meet_criteria(info, allow_class = "logical", has_length = 1)
  meet_criteria(rules, allow_class = "character", has_length = c(1, 2, 3, 4, 5), is_in = c("breakpoints", "expert", "other", "all", "custom"))
  meet_criteria(verbose, allow_class = "logical", has_length = 1)
  meet_criteria(version_breakpoints, allow_class = c("numeric", "integer"), has_length = 1, is_in = as.double(names(EUCAST_VERSION_BREAKPOINTS)))
  meet_criteria(version_expertrules, allow_class = c("numeric", "integer"), has_length = 1, is_in = as.double(names(EUCAST_VERSION_EXPERT_RULES)))
  meet_criteria(ampc_cephalosporin_resistance, allow_class = c("logical", "character", "sir"), has_length = 1, allow_NA = TRUE, allow_NULL = TRUE)
  meet_criteria(only_sir_columns, allow_class = "logical", has_length = 1)
  meet_criteria(custom_rules, allow_class = "custom_eucast_rules", allow_NULL = TRUE)

  add_MO_lookup_to_AMR_env()

  if ("custom" %in% rules && is.null(custom_rules)) {
    warning_("in `eucast_rules()`: no custom rules were set with the `custom_rules` argument",
      immediate = TRUE
    )
    rules <- rules[rules != "custom"]
    if (length(rules) == 0) {
      if (isTRUE(info)) {
        message_("No other rules were set, returning original data", add_fn = font_red, as_note = FALSE)
      }
      return(x)
    }
  }

  x_deparsed <- deparse(substitute(x))
  if (length(x_deparsed) > 1 || any(x_deparsed %unlike% "[a-z]+")) {
    x_deparsed <- "your_data"
  }

  breakpoints_info <- EUCAST_VERSION_BREAKPOINTS[[which(as.double(names(EUCAST_VERSION_BREAKPOINTS)) == version_breakpoints)]]
  expertrules_info <- EUCAST_VERSION_EXPERT_RULES[[which(as.double(names(EUCAST_VERSION_EXPERT_RULES)) == version_expertrules)]]

  # support old setting (until AMR v1.3.0)
  if (missing(rules) && !is.null(getOption("AMR.eucast_rules"))) {
    rules <- getOption("AMR.eucast_rules")
  }

  if (interactive() && isTRUE(verbose) && isTRUE(info)) {
    txt <- paste0(
      "WARNING: In Verbose mode, the eucast_rules() function does not apply rules to the data, but instead returns a data set in logbook form with extensive info about which rows and columns would be effected and in which way.",
      "\n\nThis may overwrite your existing data if you use e.g.:",
      "\ndata <- eucast_rules(data, verbose = TRUE)\n\nDo you want to continue?"
    )
    showQuestion <- import_fn("showQuestion", "rstudioapi", error_on_fail = FALSE)
    if (!is.null(showQuestion)) {
      q_continue <- showQuestion("Using verbose = TRUE with eucast_rules()", txt)
    } else {
      q_continue <- utils::menu(choices = c("OK", "Cancel"), graphics = FALSE, title = txt)
    }
    if (q_continue %in% c(FALSE, 2)) {
      message_("Cancelled, returning original data", add_fn = font_red, as_note = FALSE)
      return(x)
    }
  }

  # try to find columns based on type
  # -- mo
  if (is.null(col_mo)) {
    col_mo <- search_type_in_df(x = x, type = "mo", info = info)
    stop_if(is.null(col_mo), "`col_mo` must be set")
  }

  decimal.mark <- getOption("OutDec")
  big.mark <- ifelse(decimal.mark != ",", ",", ".")
  formatnr <- function(x, big = big.mark, dec = decimal.mark) {
    trimws(format(x, big.mark = big, decimal.mark = dec))
  }

  warned <- FALSE
  warn_lacking_sir_class <- character(0)
  txt_ok <- function(n_added, n_changed, warned = FALSE) {
    if (warned == FALSE) {
      if (n_added + n_changed == 0) {
        cat(font_subtle(" (no changes)\n"))
      } else {
        # opening
        if (n_added > 0 && n_changed == 0) {
          cat(font_green(" ("))
        } else if (n_added == 0 && n_changed > 0) {
          cat(font_blue(" ("))
        } else {
          cat(font_grey(" ("))
        }
        # additions
        if (n_added > 0) {
          if (n_added == 1) {
            cat(font_green("1 value added"))
          } else {
            cat(font_green(formatnr(n_added), "values added"))
          }
        }
        # separator
        if (n_added > 0 && n_changed > 0) {
          cat(font_grey(", "))
        }
        # changes
        if (n_changed > 0) {
          if (n_changed == 1) {
            cat(font_blue("1 value changed"))
          } else {
            cat(font_blue(formatnr(n_changed), "values changed"))
          }
        }
        # closing
        if (n_added > 0 && n_changed == 0) {
          cat(font_green(")\n"))
        } else if (n_added == 0 && n_changed > 0) {
          cat(font_blue(")\n"))
        } else {
          cat(font_grey(")\n"))
        }
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
      "FOX1",
      "GEN",
      "MFX",
      "NAL",
      "NOR",
      "PEN",
      "PIP",
      "TCY",
      "TIC",
      "TOB"
    ),
    hard_dependencies = NULL,
    verbose = verbose,
    info = info,
    only_sir_columns = only_sir_columns,
    fn = "eucast_rules",
    ...
  )

  if (!"AMP" %in% names(cols_ab) && "AMX" %in% names(cols_ab)) {
    # ampicillin column is missing, but amoxicillin is available
    if (isTRUE(info)) {
      message_("Using column '", cols_ab[names(cols_ab) == "AMX"], "' as input for ampicillin since many EUCAST rules depend on it.")
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
      vapply(FUN.VALUE = character(1), function(x) if (x %in% AMR::antibiotics$ab) ab_name(x, language = NULL, tolower = TRUE, fast_mode = TRUE) else x) %pm>%
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
        # like PEN,FOX,GEN S (although dependency on > 2 ABx does not exist at the moment)
        # nolint start
        # x <- paste(paste0(ab_names, collapse = " and "), "are all")
        # nolint end
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
        # like PEN,FOX,GEN S,R,R (although dependency on > 2 ABx does not exist at the moment)
        paste0(
          ab_names[1], " is '", ab_results[1], "' and ",
          ab_names[2], " is '", ab_results[2], "' and ",
          ab_names[3], " is '", ab_results[3], "'"
        )
      }
    }
  }
  as.sir_no_warning <- function(x) {
    if (is.sir(x)) {
      return(x)
    }
    suppressWarnings(as.sir(x))
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
  x$`.rowid` <- vapply(
    FUN.VALUE = character(1),
    as.list(as.data.frame(t(x[, c(col_mo, cols_ab), drop = FALSE]),
      stringsAsFactors = FALSE
    )),
    function(x) {
      x[is.na(x)] <- "."
      paste0(x, collapse = "")
    }
  )

  # save original table, with the new .rowid column
  x.bak <- x
  # keep only unique rows for MO and ABx
  x <- x %pm>%
    pm_arrange(`.rowid`) %pm>%
    # big speed gain! only analyse unique rows:
    pm_distinct(`.rowid`, .keep_all = TRUE) %pm>%
    as.data.frame(stringsAsFactors = FALSE)
  x[, col_mo] <- as.mo(as.character(x[, col_mo, drop = TRUE]), info = info)
  # rename col_mo to prevent interference with joined columns
  colnames(x)[colnames(x) == col_mo] <- ".col_mo"
  col_mo <- ".col_mo"
  # join to microorganisms data set
  x <- left_join_microorganisms(x, by = col_mo, suffix = c("_oldcols", ""))
  x$gramstain <- mo_gramstain(x[, col_mo, drop = TRUE], language = NULL, info = FALSE)
  x$genus_species <- trimws(paste(x$genus, x$species))
  if (isTRUE(info) && NROW(x) > 10000) {
    message_(" OK.", add_fn = list(font_green, font_bold), as_note = FALSE)
  }

  if (any(x$genus == "Staphylococcus", na.rm = TRUE)) {
    all_staph <- AMR_env$MO_lookup[which(AMR_env$MO_lookup$genus == "Staphylococcus"), , drop = FALSE]
    all_staph$CNS_CPS <- suppressWarnings(mo_name(all_staph$mo, Becker = "all", language = NULL, info = FALSE))
  }
  if (any(x$genus == "Streptococcus", na.rm = TRUE)) {
    all_strep <- AMR_env$MO_lookup[which(AMR_env$MO_lookup$genus == "Streptococcus"), , drop = FALSE]
    all_strep$Lancefield <- suppressWarnings(mo_name(all_strep$mo, Lancefield = TRUE, language = NULL, info = FALSE))
  }

  n_added <- 0
  n_changed <- 0

  # Other rules: enzyme inhibitors ------------------------------------------
  if (any(c("all", "other") %in% rules)) {
    if (isTRUE(info)) {
      cat("\n")
      cat(word_wrap(
        font_bold(paste0(
          "Rules by this AMR package (",
          font_red(paste0(
            "v", utils::packageDescription("AMR")$Version, ", ",
            format(as.Date(utils::packageDescription("AMR")$Date), format = "%Y")
          )), "), see ?eucast_rules\n"
        ))
      ))
    }
    ab_enzyme <- subset(AMR::antibiotics, name %like% "/")[, c("ab", "name"), drop = FALSE]
    colnames(ab_enzyme) <- c("enzyme_ab", "enzyme_name")
    ab_enzyme$base_name <- gsub("^([a-zA-Z0-9]+).*", "\\1", ab_enzyme$enzyme_name)
    ab_enzyme$base_ab <- AMR::antibiotics[match(ab_enzyme$base_name, AMR::antibiotics$name), "ab", drop = TRUE]
    ab_enzyme <- subset(ab_enzyme, !is.na(base_ab))
    # make ampicillin and amoxicillin interchangable
    ampi <- subset(ab_enzyme, base_ab == "AMX")
    ampi$base_ab <- "AMP"
    ampi$base_name <- ab_name("AMP", language = NULL)
    amox <- subset(ab_enzyme, base_ab == "AMP")
    amox$base_ab <- "AMX"
    amox$base_name <- ab_name("AMX", language = NULL)
    # merge and sort
    ab_enzyme <- rbind(ab_enzyme, ampi, amox)
    ab_enzyme <- ab_enzyme[order(ab_enzyme$enzyme_name), , drop = FALSE]

    for (i in seq_len(nrow(ab_enzyme))) {
      # check if both base and base + enzyme inhibitor are part of the data set
      if (all(c(ab_enzyme$base_ab[i], ab_enzyme$enzyme_ab[i]) %in% names(cols_ab), na.rm = TRUE)) {
        col_base <- unname(cols_ab[ab_enzyme$base_ab[i]])
        col_enzyme <- unname(cols_ab[ab_enzyme$enzyme_ab[i]])

        # Set base to R where base + enzyme inhibitor is R ----
        rule_current <- paste0(
          ab_enzyme$base_name[i], " ('", font_bold(col_base), "') = R if ",
          tolower(ab_enzyme$enzyme_name[i]), " ('", font_bold(col_enzyme), "') = R"
        )
        if (isTRUE(info)) {
          cat(word_wrap(rule_current,
            width = getOption("width") - 30,
            extra_indent = 6
          ))
        }
        run_changes <- edit_sir(
          x = x,
          to = "R",
          rule = c(
            rule_current, "Other rules", "",
            paste0("Non-EUCAST: AMR package v", utils::packageDescription("AMR")$Version)
          ),
          rows = which(as.sir_no_warning(x[, col_enzyme, drop = TRUE]) == "R"),
          cols = col_base,
          last_verbose_info = verbose_info,
          original_data = x.bak,
          warned = warned,
          info = info,
          verbose = verbose
        )
        n_added <- n_added + run_changes$added
        n_changed <- n_changed + run_changes$changed
        verbose_info <- run_changes$verbose_info
        x <- run_changes$output
        warn_lacking_sir_class <- c(warn_lacking_sir_class, run_changes$sir_warn)
        # Print number of new changes
        if (isTRUE(info)) {
          # print only on last one of rules in this group
          txt_ok(n_added = n_added, n_changed = n_changed, warned = warned)
          # and reset counters
          n_added <- 0
          n_changed <- 0
        }

        # Set base + enzyme inhibitor to S where base is S ----
        rule_current <- paste0(
          ab_enzyme$enzyme_name[i], " ('", font_bold(col_enzyme), "') = S if ",
          tolower(ab_enzyme$base_name[i]), " ('", font_bold(col_base), "') = S"
        )

        if (isTRUE(info)) {
          cat(word_wrap(rule_current,
            width = getOption("width") - 30,
            extra_indent = 6
          ))
        }
        run_changes <- edit_sir(
          x = x,
          to = "S",
          rule = c(
            rule_current, "Other rules", "",
            paste0("Non-EUCAST: AMR package v", utils::packageDescription("AMR")$Version)
          ),
          rows = which(as.sir_no_warning(x[, col_base, drop = TRUE]) == "S"),
          cols = col_enzyme,
          last_verbose_info = verbose_info,
          original_data = x.bak,
          warned = warned,
          info = info,
          verbose = verbose
        )
        n_added <- n_added + run_changes$added
        n_changed <- n_changed + run_changes$changed
        verbose_info <- run_changes$verbose_info
        x <- run_changes$output
        warn_lacking_sir_class <- c(warn_lacking_sir_class, run_changes$sir_warn)
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
      message_("Skipping inheritance rules defined by this AMR package, such as setting trimethoprim (TMP) = R where trimethoprim/sulfamethoxazole (SXT) = R. Add \"other\" or \"all\" to the `rules` argument to apply those rules.")
    }
  }

  if (!any(c("all", "custom") %in% rules) && !is.null(custom_rules)) {
    if (isTRUE(info)) {
      message_("Skipping custom EUCAST rules, since the `rules` argument does not contain \"custom\".")
    }
    custom_rules <- NULL
  }

  # Official EUCAST rules ---------------------------------------------------
  eucast_notification_shown <- FALSE
  if (!is.null(list(...)$eucast_rules_df)) {
    # this allows: eucast_rules(x, eucast_rules_df = AMR:::EUCAST_RULES_DF %>% filter(is.na(have_these_values)))
    eucast_rules_df <- list(...)$eucast_rules_df
  } else {
    # otherwise internal data file, created in data-raw/_pre_commit_hook.R
    eucast_rules_df <- EUCAST_RULES_DF
  }

  # filter on user-set guideline versions ----
  if (any(c("all", "breakpoints") %in% rules)) {
    eucast_rules_df <- subset(
      eucast_rules_df,
      reference.rule_group %unlike% "breakpoint" |
        (reference.rule_group %like% "breakpoint" & reference.version == version_breakpoints)
    )
  }
  if (any(c("all", "expert") %in% rules)) {
    eucast_rules_df <- subset(
      eucast_rules_df,
      reference.rule_group %unlike% "expert" |
        (reference.rule_group %like% "expert" & reference.version == version_expertrules)
    )
  }
  # filter out AmpC de-repressed cephalosporin-resistant mutants ----
  # no need to filter on version number here - the rules contain these version number, so are inherently filtered
  # cefotaxime, ceftriaxone, ceftazidime
  if (is.null(ampc_cephalosporin_resistance) || isFALSE(ampc_cephalosporin_resistance)) {
    eucast_rules_df <- subset(
      eucast_rules_df,
      reference.rule %unlike% "ampc"
    )
  } else {
    if (isTRUE(ampc_cephalosporin_resistance)) {
      ampc_cephalosporin_resistance <- "R"
    }
    eucast_rules_df[which(eucast_rules_df$reference.rule %like% "ampc"), "to_value"] <- as.character(ampc_cephalosporin_resistance)
  }

  # Go over all rules and apply them ----
  for (i in seq_len(nrow(eucast_rules_df))) {
    rule_previous <- eucast_rules_df[max(1, i - 1), "reference.rule", drop = TRUE]
    rule_current <- eucast_rules_df[i, "reference.rule", drop = TRUE]
    rule_next <- eucast_rules_df[min(nrow(eucast_rules_df), i + 1), "reference.rule", drop = TRUE]
    rule_group_previous <- eucast_rules_df[max(1, i - 1), "reference.rule_group", drop = TRUE]
    rule_group_current <- eucast_rules_df[i, "reference.rule_group", drop = TRUE]
    # don't apply rules if user doesn't want to apply them
    if (rule_group_current %like% "breakpoint" && !any(c("all", "breakpoints") %in% rules)) {
      next
    }
    if (rule_group_current %like% "expert" && !any(c("all", "expert") %in% rules)) {
      next
    }

    if (isFALSE(info) || isFALSE(verbose)) {
      rule_text <- ""
    } else {
      if (is.na(eucast_rules_df[i, "and_these_antibiotics", drop = TRUE])) {
        rule_text <- paste0("always report as '", eucast_rules_df[i, "to_value", drop = TRUE], "': ", get_antibiotic_names(eucast_rules_df[i, "then_change_these_antibiotics", drop = TRUE]))
      } else {
        rule_text <- paste0(
          "report as '", eucast_rules_df[i, "to_value", drop = TRUE], "' when ",
          format_antibiotic_names(
            ab_names = get_antibiotic_names(eucast_rules_df[i, "and_these_antibiotics", drop = TRUE]),
            ab_results = eucast_rules_df[i, "have_these_values", drop = TRUE]
          ), ": ",
          get_antibiotic_names(eucast_rules_df[i, "then_change_these_antibiotics", drop = TRUE])
        )
      }
    }
    if (i == 1) {
      rule_previous <- ""
      rule_group_previous <- ""
    }
    if (i == nrow(eucast_rules_df)) {
      rule_next <- ""
    }

    if (isTRUE(info)) {
      # Print EUCAST intro ------------------------------------------------------
      if (rule_group_current %unlike% "other" && eucast_notification_shown == FALSE) {
        cat(
          paste0(
            "\n", font_grey(strrep("-", 0.95 * getOption("width", 100))), "\n",
            word_wrap("Rules by the ", font_bold("European Committee on Antimicrobial Susceptibility Testing (EUCAST)")), "\n",
            font_blue("https://eucast.org/"), "\n"
          )
        )
        eucast_notification_shown <- TRUE
      }

      # Print rule (group) ------------------------------------------------------
      if (rule_group_current != rule_group_previous) {
        # is new rule group, one of Breakpoints, Expert Rules and Other
        cat(font_bold(
          ifelse(
            rule_group_current %like% "breakpoint",
            paste0(
              "\n",
              word_wrap(
                breakpoints_info$title, " (",
                font_red(paste0(breakpoints_info$version_txt, ", ", breakpoints_info$year)), ")\n"
              )
            ),
            ifelse(
              rule_group_current %like% "expert",
              paste0(
                "\n",
                word_wrap(
                  expertrules_info$title, " (",
                  font_red(paste0(expertrules_info$version_txt, ", ", expertrules_info$year)), ")\n"
                )
              ),
              ""
            )
          )
        ), "\n")
      }
      # Print rule  -------------------------------------------------------------
      if (rule_current != rule_previous) {
        # is new rule within group, print its name
        cat(italicise_taxonomy(word_wrap(rule_current,
          width = getOption("width") - 30,
          extra_indent = 6
        ),
        type = "ansi"
        ))
        warned <- FALSE
      }
    }

    # Get rule from file ------------------------------------------------------
    if_mo_property <- trimws(eucast_rules_df[i, "if_mo_property", drop = TRUE])
    like_is_one_of <- trimws(eucast_rules_df[i, "like.is.one_of", drop = TRUE])
    mo_value <- trimws(eucast_rules_df[i, "this_value", drop = TRUE])

    # be sure to comprise all coagulase-negative/-positive staphylococci when they are mentioned
    if (mo_value %like% "coagulase" && any(x$genus == "Staphylococcus", na.rm = TRUE)) {
      if (mo_value %like% "negative") {
        eucast_rules_df[i, "this_value"] <- paste0(
          "^(", paste0(all_staph[which(all_staph$CNS_CPS %like% "negative"),
            "fullname",
            drop = TRUE
          ],
          collapse = "|"
          ),
          ")$"
        )
      } else {
        eucast_rules_df[i, "this_value"] <- paste0(
          "^(", paste0(all_staph[which(all_staph$CNS_CPS %like% "positive"),
            "fullname",
            drop = TRUE
          ],
          collapse = "|"
          ),
          ")$"
        )
      }
      like_is_one_of <- "like"
    }
    # be sure to comprise all beta-haemolytic Streptococci (Lancefield groups A, B, C and G) when they are mentioned
    if (mo_value %like% "group [ABCG]" && any(x$genus == "Streptococcus", na.rm = TRUE)) {
      eucast_rules_df[i, "this_value"] <- paste0(
        "^(", paste0(all_strep[which(all_strep$Lancefield %like% "group [ABCG]"),
          "fullname",
          drop = TRUE
        ],
        collapse = "|"
        ),
        ")$"
      )
      like_is_one_of <- "like"
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
      stop("invalid value for column 'like.is.one_of'", call. = FALSE)
    }

    source_antibiotics <- eucast_rules_df[i, "and_these_antibiotics", drop = TRUE]
    source_value <- trimws(unlist(strsplit(eucast_rules_df[i, "have_these_values", drop = TRUE], ",", fixed = TRUE)))
    target_antibiotics <- eucast_rules_df[i, "then_change_these_antibiotics", drop = TRUE]
    target_value <- eucast_rules_df[i, "to_value", drop = TRUE]

    if (is.na(source_antibiotics)) {
      rows <- tryCatch(which(x[, if_mo_property, drop = TRUE] %like% mo_value),
        error = function(e) integer(0)
      )
    } else {
      source_antibiotics <- get_ab_from_namespace(source_antibiotics, cols_ab)
      if (length(source_value) == 1 && length(source_antibiotics) > 1) {
        source_value <- rep(source_value, length(source_antibiotics))
      }
      if (length(source_antibiotics) == 0) {
        rows <- integer(0)
      } else if (length(source_antibiotics) == 1) {
        rows <- tryCatch(which(x[, if_mo_property, drop = TRUE] %like% mo_value &
          as.sir_no_warning(x[, source_antibiotics[1L]]) == source_value[1L]),
        error = function(e) integer(0)
        )
      } else if (length(source_antibiotics) == 2) {
        rows <- tryCatch(which(x[, if_mo_property, drop = TRUE] %like% mo_value &
          as.sir_no_warning(x[, source_antibiotics[1L]]) == source_value[1L] &
          as.sir_no_warning(x[, source_antibiotics[2L]]) == source_value[2L]),
        error = function(e) integer(0)
        )
        # nolint start
        # } else if (length(source_antibiotics) == 3) {
        #   rows <-  tryCatch(which(x[, if_mo_property, drop = TRUE] %like% mo_value
        #                           & as.sir_no_warning(x[, source_antibiotics[1L]]) == source_value[1L]
        #                           & as.sir_no_warning(x[, source_antibiotics[2L]]) == source_value[2L]
        #                           & as.sir_no_warning(x[, source_antibiotics[3L]]) == source_value[3L]),
        #                     error = function(e) integer(0))
        # nolint end
      } else {
        stop_("only 2 antibiotics supported for source_antibiotics")
      }
    }

    cols <- get_ab_from_namespace(target_antibiotics, cols_ab)

    # Apply rule on data ------------------------------------------------------
    # this will return the unique number of changes
    run_changes <- edit_sir(
      x = x,
      to = target_value,
      rule = c(
        rule_text, rule_group_current, rule_current,
        ifelse(rule_group_current %like% "breakpoint",
          paste0(breakpoints_info$title, " ", breakpoints_info$version_txt, ", ", breakpoints_info$year),
          paste0(expertrules_info$title, " ", expertrules_info$version_txt, ", ", expertrules_info$year)
        )
      ),
      rows = rows,
      cols = cols,
      last_verbose_info = verbose_info,
      original_data = x.bak,
      warned = warned,
      info = info,
      verbose = verbose
    )
    n_added <- n_added + run_changes$added
    n_changed <- n_changed + run_changes$changed
    verbose_info <- run_changes$verbose_info
    x <- run_changes$output
    warn_lacking_sir_class <- c(warn_lacking_sir_class, run_changes$sir_warn)
    # Print number of new changes ---------------------------------------------
    if (isTRUE(info) && rule_next != rule_current) {
      # print only on last one of rules in this group
      txt_ok(n_added = n_added, n_changed = n_changed, warned = warned)
      # and reset counters
      n_added <- 0
      n_changed <- 0
    }
  } # end of going over all rules

  # Apply custom rules ----
  if (!is.null(custom_rules)) {
    if (isTRUE(info)) {
      cat("\n")
      cat(font_bold("Custom EUCAST rules, set by user"), "\n")
    }
    for (i in seq_len(length(custom_rules))) {
      rule <- custom_rules[[i]]
      rows <- which(eval(parse(text = rule$query), envir = x))
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
        # print rule
        cat(italicise_taxonomy(word_wrap(format_custom_query_rule(rule$query, colours = FALSE),
          width = getOption("width") - 30,
          extra_indent = 6
        ),
        type = "ansi"
        ))
        warned <- FALSE
      }
      run_changes <- edit_sir(
        x = x,
        to = target_value,
        rule = c(
          rule_text,
          "Custom EUCAST rules",
          paste0("Custom EUCAST rule ", i),
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
        verbose = verbose
      )
      n_added <- n_added + run_changes$added
      n_changed <- n_changed + run_changes$changed
      verbose_info <- run_changes$verbose_info
      x <- run_changes$output
      warn_lacking_sir_class <- c(warn_lacking_sir_class, run_changes$sir_warn)
      # Print number of new changes ---------------------------------------------
      if (isTRUE(info) && rule_next != rule_current) {
        # print only on last one of rules in this group
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
    rownames(verbose_info) <- NULL
  }

  if (isTRUE(info)) {
    if (isTRUE(verbose)) {
      wouldve <- "would have "
    } else {
      wouldve <- ""
    }

    cat(paste0("\n", font_grey(strrep("-", 0.95 * getOption("width", 100))), "\n"))
    cat(word_wrap(paste0(
      "The rules ", paste0(wouldve, "affected "),
      font_bold(
        formatnr(pm_n_distinct(verbose_info$row)),
        "out of", formatnr(nrow(x.bak)),
        "rows"
      ),
      ", making a total of ",
      font_bold(formatnr(nrow(verbose_info)), "edits\n")
    )))

    total_n_added <- verbose_info %pm>% pm_filter(is.na(old)) %pm>% nrow()
    total_n_changed <- verbose_info %pm>% pm_filter(!is.na(old)) %pm>% nrow()

    # print added values
    if (total_n_added == 0) {
      colour <- cat # is function
    } else {
      colour <- font_green # is function
    }
    cat(colour(paste0(
      "=> ", wouldve, "added ",
      font_bold(formatnr(verbose_info %pm>%
        pm_filter(is.na(old)) %pm>%
        nrow()), "test results"),
      "\n"
    )))
    if (total_n_added > 0) {
      added_summary <- verbose_info %pm>%
        pm_filter(is.na(old)) %pm>%
        pm_count(new, name = "n")
      cat(paste("   -",
        paste0(
          formatnr(added_summary$n), " test result", ifelse(added_summary$n > 1, "s", ""),
          " added as ", paste0('"', added_summary$new, '"')
        ),
        collapse = "\n"
      ))
    }

    # print changed values
    if (total_n_changed == 0) {
      colour <- cat # is function
    } else {
      colour <- font_blue # is function
    }
    if (total_n_added + total_n_changed > 0) {
      cat("\n")
    }
    cat(colour(paste0(
      "=> ", wouldve, "changed ",
      font_bold(formatnr(verbose_info %pm>%
        pm_filter(!is.na(old)) %pm>%
        nrow()), "test results"),
      "\n"
    )))
    if (total_n_changed > 0) {
      changed_summary <- verbose_info %pm>%
        pm_filter(!is.na(old)) %pm>%
        pm_mutate(new = ifelse(is.na(new), "NA", new)) %pm>%
        pm_count(old, new, name = "n")
      cat(paste("   -",
        paste0(
          formatnr(changed_summary$n), " test result", ifelse(changed_summary$n > 1, "s", ""), " changed from ",
          paste0('"', changed_summary$old, '"'), " to ", paste0('"', changed_summary$new, '"')
        ),
        collapse = "\n"
      ))
      cat("\n")
    }

    cat(paste0(font_grey(strrep("-", 0.95 * getOption("width", 100))), "\n"))

    if (isFALSE(verbose) && total_n_added + total_n_changed > 0) {
      cat("\n", word_wrap("Use ", font_bold("eucast_rules(..., verbose = TRUE)"), " (on your original data) to get a data.frame with all specified edits instead."), "\n\n", sep = "")
    } else if (isTRUE(verbose)) {
      cat("\n", word_wrap("Used 'Verbose mode' (", font_bold("verbose = TRUE"), "), which returns a data.frame with all specified edits.\nUse ", font_bold("verbose = FALSE"), " to apply the rules on your data."), "\n\n", sep = "")
    }
  }

  if (length(warn_lacking_sir_class) > 0) {
    warn_lacking_sir_class <- unique(warn_lacking_sir_class)
    # take order from original data set
    warn_lacking_sir_class <- warn_lacking_sir_class[order(colnames(x.bak))]
    warn_lacking_sir_class <- warn_lacking_sir_class[!is.na(warn_lacking_sir_class)]
    warning_(
      "in `eucast_rules()`: not all columns with antimicrobial results are of class 'sir'. Transform them on beforehand, with e.g.:\n",
      "  - ", x_deparsed, " %>% as.sir(", ifelse(length(warn_lacking_sir_class) == 1,
        warn_lacking_sir_class,
        paste0(warn_lacking_sir_class[1], ":", warn_lacking_sir_class[length(warn_lacking_sir_class)])
      ), ")\n",
      "  - ", x_deparsed, " %>% mutate_if(is_sir_eligible, as.sir)\n",
      "  - ", x_deparsed, " %>% mutate(across(where(is_sir_eligible), as.sir))"
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
                     verbose) {
  cols <- unique(cols[!is.na(cols) & !is.null(cols)])

  # for Verbose Mode, keep track of all changes and return them
  track_changes <- list(
    added = 0,
    changed = 0,
    output = x,
    verbose_info = last_verbose_info,
    sir_warn = character(0)
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
    tryCatch(
      # insert into original table
      new_edits[rows, cols] <- to,
      warning = function(w) {
        if (w$message %like% "invalid factor level") {
          xyz <- vapply(FUN.VALUE = logical(1), cols, function(col) {
            new_edits[, col] <<- factor(
              x = as.character(pm_pull(new_edits, col)),
              levels = unique(c(to, levels(pm_pull(new_edits, col))))
            )
            TRUE
          })
          suppressWarnings(new_edits[rows, cols] <<- to)
          warning_(
            "in `eucast_rules()`: value \"", to, "\" added to the factor levels of column",
            ifelse(length(cols) == 1, "", "s"),
            " ", vector_and(cols, quotes = "`", sort = FALSE),
            " because this value was not an existing factor level."
          )
          txt_warning()
          warned <- FALSE
        } else {
          warning_("in `eucast_rules()`: ", w$message)
          txt_warning()
        }
      },
      error = function(e) {
        txt_error()
        stop(paste0(
          "In row(s) ", paste(rows[seq_len(min(length(rows), 10))], collapse = ","),
          ifelse(length(rows) > 10, "...", ""),
          " while writing value '", to,
          "' to column(s) `", paste(cols, collapse = "`, `"),
          "`:\n", e$message
        ),
        call. = FALSE
        )
      }
    )

    track_changes$output <- new_edits
    if ((isTRUE(info) || isTRUE(verbose)) && !isTRUE(all.equal(x, track_changes$output))) {
      get_original_rows <- function(rowids) {
        as.integer(rownames(original_data[which(original_data$.rowid %in% rowids), , drop = FALSE]))
      }
      for (i in seq_len(length(cols))) {
        verbose_new <- data.frame(
          rowid = new_edits[rows, ".rowid", drop = TRUE],
          col = cols[i],
          mo_fullname = new_edits[rows, "fullname", drop = TRUE],
          old = x[rows, cols[i], drop = TRUE],
          new = to,
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
        verbose_new <- verbose_new %pm>% pm_filter(old != new | is.na(old) | is.na(new) & !is.na(old))
        # save changes to data set 'verbose_info'
        track_changes$verbose_info <- rbind(track_changes$verbose_info,
          verbose_new,
          stringsAsFactors = FALSE
        )
        # count adds and changes
        track_changes$added <- track_changes$added + verbose_new %pm>%
          pm_filter(is.na(old)) %pm>%
          pm_pull(rowid) %pm>%
          get_original_rows() %pm>%
          length()
        track_changes$changed <- track_changes$changed + verbose_new %pm>%
          pm_filter(!is.na(old)) %pm>%
          pm_pull(rowid) %pm>%
          get_original_rows() %pm>%
          length()
      }
    }
  }
  return(track_changes)
}

#' @rdname eucast_rules
#' @export
eucast_dosage <- function(ab, administration = "iv", version_breakpoints = 12.0) {
  meet_criteria(ab, allow_class = c("character", "numeric", "integer", "factor"))
  meet_criteria(administration, allow_class = "character", is_in = dosage$administration[!is.na(dosage$administration)], has_length = 1)
  meet_criteria(version_breakpoints, allow_class = c("numeric", "integer"), has_length = 1, is_in = as.double(names(EUCAST_VERSION_BREAKPOINTS)))

  # show used version_breakpoints number once per session (AMR_env will reload every session)
  if (message_not_thrown_before("eucast_dosage", "v", gsub("[^0-9]", "", version_breakpoints), entire_session = TRUE)) {
    message_(
      "Dosages for antimicrobial drugs, as meant for ",
      format_eucast_version_nr(version_breakpoints, markdown = FALSE), ". ",
      font_red("This note will be shown once per session.")
    )
  }

  ab <- as.ab(ab)
  lst <- vector("list", length = length(ab))
  for (i in seq_len(length(ab))) {
    df <- AMR::dosage[which(AMR::dosage$ab == ab[i] & AMR::dosage$administration == administration), , drop = FALSE]
    lst[[i]] <- list(
      ab = "",
      name = "",
      standard_dosage = ifelse("standard_dosage" %in% df$type,
        df[which(df$type == "standard_dosage"), "original_txt", drop = TRUE],
        NA_character_
      ),
      high_dosage = ifelse("high_dosage" %in% df$type,
        df[which(df$type == "high_dosage"), "original_txt", drop = TRUE],
        NA_character_
      )
    )
  }
  out <- do.call("rbind", lapply(lst, as.data.frame, stringsAsFactors = FALSE))
  rownames(out) <- NULL
  out$ab <- ab
  out$name <- ab_name(ab, language = NULL)
  if (pkg_is_available("tibble", also_load = FALSE)) {
    import_fn("as_tibble", "tibble")(out)
  } else {
    out
  }
}
