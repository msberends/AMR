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

#' Determine Multidrug-Resistant Organisms (MDRO)
#'
#' Determine which isolates are multidrug-resistant organisms (MDRO) according to international, national, or custom guidelines.
#' @param x A [data.frame] with antimicrobials columns, like `AMX` or `amox`. Can be left blank for automatic determination.
#' @param guideline A specific guideline to follow, see sections *Supported international / national guidelines* and *Using Custom Guidelines* below. When left empty, the publication by Magiorakos *et al.* (see below) will be followed.
#' @param esbl [logical] values, or a column name containing logical values, indicating the presence of an ESBL gene (or production of its proteins).
#' @param carbapenemase [logical] values, or a column name containing logical values, indicating the presence of a carbapenemase gene (or production of its proteins).
#' @param mecA [logical] values, or a column name containing logical values, indicating the presence of a *mecA* gene (or production of its proteins).
#' @param mecC [logical] values, or a column name containing logical values, indicating the presence of a *mecC* gene (or production of its proteins).
#' @param vanA [logical] values, or a column name containing logical values, indicating the presence of a *vanA* gene (or production of its proteins).
#' @param vanB [logical] values, or a column name containing logical values, indicating the presence of a *vanB* gene (or production of its proteins).
#' @inheritParams eucast_rules
#' @param pct_required_classes Minimal required percentage of antimicrobial classes that must be available per isolate, rounded down. For example, with the default guideline, 17 antimicrobial classes must be available for *S. aureus*. Setting this `pct_required_classes` argument to `0.5` (default) means that for every *S. aureus* isolate at least 8 different classes must be available. Any lower number of available classes will return `NA` for that isolate.
#' @param combine_SI A [logical] to indicate whether all values of S and I must be merged into one, so resistance is only considered when isolates are R, not I. As this is the default behaviour of the [mdro()] function, it follows the redefinition by EUCAST about the interpretation of I (increased exposure) in 2019, see section 'Interpretation of S, I and R' below. When using `combine_SI = FALSE`, resistance is considered when isolates are R or I.
#' @param verbose A [logical] to turn Verbose mode on and off (default is off). In Verbose mode, the function does not return the MDRO results, but instead returns a data set in logbook form with extensive info about which isolates would be MDRO-positive, or why they are not.
#' @details
#' These functions are context-aware. This means that the `x` argument can be left blank if used inside a [data.frame] call, see *Examples*.
#'
#' For the `pct_required_classes` argument, values above 1 will be divided by 100. This is to support both fractions (`0.75` or `3/4`) and percentages (`75`).
#'
#' **Note:** Every test that involves the Enterobacteriaceae family, will internally be performed using its newly named *order* Enterobacterales, since the Enterobacteriaceae family has been taxonomically reclassified by Adeolu *et al.* in 2016. Before that, Enterobacteriaceae was the only family under the Enterobacteriales (with an i) order. All species under the old Enterobacteriaceae family are still under the new Enterobacterales (without an i) order, but divided into multiple families. The way tests are performed now by this [mdro()] function makes sure that results from before 2016 and after 2016 are identical.
#'
#' ### Supported International / National Guidelines
#'
#' Please suggest to implement guidelines by [letting us know](https://github.com/msberends/AMR/issues/new?template=2-feature-request.yml&title=Add%20new%20MDRO%20guideline).
#'
#' Currently supported guidelines are (case-insensitive):
#'
#' * `guideline = "CMI 2012"` (default)
#'
#'   Magiorakos AP, Srinivasan A *et al.* "Multidrug-resistant, extensively drug-resistant and pandrug-resistant bacteria: an international expert proposal for interim standard definitions for acquired resistance." Clinical Microbiology and Infection (2012) (\doi{10.1111/j.1469-0691.2011.03570.x})
#'
#' * `guideline = "EUCAST 3.3"` (or simply `guideline = "EUCAST"`)
#'
#'   The European international guideline - EUCAST Expert Rules Version 3.3 "Intrinsic Resistance and Unusual Phenotypes" ([link](https://www.eucast.org/fileadmin/src/media/PDFs/EUCAST_files/Expert_Rules/2021/Intrinsic_Resistance_and_Unusual_Phenotypes_Tables_v3.3_20211018.pdf))
#'
#'   Also:
#'
#'   * `guideline = "EUCAST 3.2"`
#'
#'      The former European international guideline - EUCAST Expert Rules Version 3.2 "Intrinsic Resistance and Unusual Phenotypes" ([link](https://www.eucast.org/fileadmin/src/media/PDFs/EUCAST_files/Expert_Rules/2020/Intrinsic_Resistance_and_Unusual_Phenotypes_Tables_v3.2_20200225.pdf))
#'
#'   * `guideline = "EUCAST 3.1"`
#'
#'      The former European international guideline - EUCAST Expert Rules Version 3.1 "Intrinsic Resistance and Exceptional Phenotypes Tables" ([link](https://www.eucast.org/fileadmin/src/media/PDFs/EUCAST_files/Expert_Rules/Expert_rules_intrinsic_exceptional_V3.1.pdf))
#'
#' * `guideline = "TB"`
#'
#'   The international guideline for multi-drug resistant tuberculosis - World Health Organization "Companion handbook to the WHO guidelines for the programmatic management of drug-resistant tuberculosis" ([link](https://www.who.int/publications/i/item/9789241548809))
#'
#' * `guideline = "MRGN"`
#'
#'   The German national guideline - Mueller et al. (2015) Antimicrobial Resistance and Infection Control 4:7; \doi{10.1186/s13756-015-0047-6}
#'
#' * `guideline = "BRMO 2024"` (or simply `guideline = "BRMO"`)
#'
#'   The Dutch national guideline - Samenwerkingverband Richtlijnen Infectiepreventie (SRI) (2024) "Bijzonder Resistente Micro-Organismen (BRMO)" ([link](https://www.sri-richtlijnen.nl/brmo))
#'
#'   Also:
#'
#'   * `guideline = "BRMO 2017"`
#'
#'     The former Dutch national guideline - Werkgroep Infectiepreventie (WIP), RIVM, last revision as of 2017: "Bijzonder Resistente Micro-Organismen (BRMO)"
#'
#' ### Using Custom Guidelines
#'
#' Using a custom MDRO guideline is of importance if you have custom rules to determine MDROs in your hospital, e.g., rules that are dependent on ward, state of contact isolation or other variables in your data.
#'
#' Custom guidelines can be set with the [custom_mdro_guideline()] function.
#'
#' @inheritSection as.sir Interpretation of SIR
#' @return
#' - If `verbose` is set to `TRUE`:\cr
#'   A [data.frame] containing columns `row_number`, `microorganism`, `MDRO`, `reason`, `all_nonsusceptible_columns`, `guideline`
#' - CMI 2012 paper - function [mdr_cmi2012()] or [mdro()]:\cr
#'   Ordered [factor] with levels `Negative` < `Multi-drug-resistant (MDR)` < `Extensively drug-resistant (XDR)` < `Pandrug-resistant (PDR)`
#' - TB guideline - function [mdr_tb()] or [`mdro(..., guideline = "TB")`][mdro()]:\cr
#'   Ordered [factor] with levels `Negative` < `Mono-resistant` < `Poly-resistant` < `Multi-drug-resistant` < `Extensively drug-resistant`
#' - German guideline - function [mrgn()] or [`mdro(..., guideline = "MRGN")`][mdro()]:\cr
#'   Ordered [factor] with levels `Negative` < `3MRGN` < `4MRGN`
#' - Everything else, except for custom guidelines:\cr
#'   Ordered [factor] with levels `Negative` < `Positive, unconfirmed` < `Positive`. The value `"Positive, unconfirmed"` means that, according to the guideline, it is not entirely sure if the isolate is multi-drug resistant and this should be confirmed with additional (e.g. genotypic) tests
#' @rdname mdro
#' @aliases MDR XDR PDR BRMO 3MRGN 4MRGN
#' @seealso [custom_mdro_guideline()]
#' @export
#' @examples
#' out <- mdro(example_isolates)
#' str(out)
#' table(out)
#'
#' out <- mdro(example_isolates, guideline = "EUCAST 3.3")
#' table(out)
#'
#' \donttest{
#' if (require("dplyr")) {
#'   # no need to define `x` when used inside dplyr verbs:
#'   example_isolates %>%
#'     mutate(MDRO = mdro()) %>%
#'     count(MDRO)
#' }
#' }
mdro <- function(x = NULL,
                 guideline = "CMI 2012",
                 col_mo = NULL,
                 esbl = NA,
                 carbapenemase = NA,
                 mecA = NA,
                 mecC = NA,
                 vanA = NA,
                 vanB = NA,
                 info = interactive(),
                 pct_required_classes = 0.5,
                 combine_SI = TRUE,
                 verbose = FALSE,
                 only_sir_columns = any(is.sir(x)),
                 ...) {
  if (is_null_or_grouped_tbl(x)) {
    # when `x` is left blank, auto determine it (get_current_data() searches underlying data within call)
    # is also a fix for using a grouped df as input (i.e., a dot as first argument)
    x <- tryCatch(get_current_data(arg_name = "x", call = -2), error = function(e) x)
  }

  meet_criteria(x, allow_class = "data.frame") # also checks dimensions to be >0
  meet_criteria(guideline, allow_class = c("list", "character"), allow_NULL = TRUE)
  if (!is.list(guideline)) meet_criteria(guideline, allow_class = "character", has_length = 1, allow_NULL = TRUE)
  meet_criteria(esbl, allow_class = c("logical", "character"), allow_NA = TRUE)
  meet_criteria(carbapenemase, allow_class = c("logical", "character"), allow_NA = TRUE)
  meet_criteria(mecA, allow_class = c("logical", "character"), allow_NA = TRUE)
  meet_criteria(mecC, allow_class = c("logical", "character"), allow_NA = TRUE)
  meet_criteria(vanA, allow_class = c("logical", "character"), allow_NA = TRUE)
  meet_criteria(vanB, allow_class = c("logical", "character"), allow_NA = TRUE)
  meet_criteria(col_mo, allow_class = "character", has_length = 1, is_in = colnames(x), allow_NULL = TRUE)
  meet_criteria(info, allow_class = "logical", has_length = 1)
  meet_criteria(pct_required_classes, allow_class = "numeric", has_length = 1)
  meet_criteria(combine_SI, allow_class = "logical", has_length = 1)
  meet_criteria(verbose, allow_class = "logical", has_length = 1)
  meet_criteria(only_sir_columns, allow_class = "logical", has_length = 1)


  if (isTRUE(only_sir_columns) && !any(is.sir(x))) {
    stop_("There were no SIR columns found in the data set, despite `only_sir_columns` being `TRUE`. Transform columns with `as.sir()` for valid antimicrobial interpretations.")
  } else if (!isTRUE(only_sir_columns) && !any(is.sir(x)) && !any(is_sir_eligible(x))) {
    stop_("There were no eligible SIR columns found in the data set. Transform columns with `as.sir()` for valid antimicrobial interpretations.")
  }

  # get gene values as TRUE/FALSE
  if (is.character(esbl)) {
    meet_criteria(esbl, is_in = colnames(x), allow_NA = FALSE, has_length = 1)
    esbl <- x[[esbl]]
    meet_criteria(esbl, allow_class = "logical", allow_NA = TRUE)
  } else if (length(esbl) == 1) {
    esbl <- rep(esbl, NROW(x))
  }
  if (is.character(carbapenemase)) {
    meet_criteria(carbapenemase, is_in = colnames(x), allow_NA = FALSE, has_length = 1)
    carbapenemase <- x[[carbapenemase]]
    meet_criteria(carbapenemase, allow_class = "logical", allow_NA = TRUE)
  } else if (length(carbapenemase) == 1) {
    carbapenemase <- rep(carbapenemase, NROW(x))
  }
  if (is.character(mecA)) {
    meet_criteria(mecA, is_in = colnames(x), allow_NA = FALSE, has_length = 1)
    mecA <- x[[mecA]]
    meet_criteria(mecA, allow_class = "logical", allow_NA = TRUE)
  } else if (length(mecA) == 1) {
    mecA <- rep(mecA, NROW(x))
  }
  if (is.character(mecC)) {
    meet_criteria(mecC, is_in = colnames(x), allow_NA = FALSE, has_length = 1)
    mecC <- x[[mecC]]
    meet_criteria(mecC, allow_class = "logical", allow_NA = TRUE)
  } else if (length(mecC) == 1) {
    mecC <- rep(mecC, NROW(x))
  }
  if (is.character(vanA)) {
    meet_criteria(vanA, is_in = colnames(x), allow_NA = FALSE, has_length = 1)
    vanA <- x[[VanA]]
    meet_criteria(vanA, allow_class = "logical", allow_NA = TRUE)
  } else if (length(vanA) == 1) {
    vanA <- rep(vanA, NROW(x))
  }
  if (is.character(vanB)) {
    meet_criteria(vanB, is_in = colnames(x), allow_NA = FALSE, has_length = 1)
    vanB <- x[[VanB]]
    meet_criteria(vanB, allow_class = "logical", allow_NA = TRUE)
  } else if (length(vanB) == 1) {
    vanB <- rep(vanB, NROW(x))
  }

  info.bak <- info
  # don't throw info's more than once per call
  if (isTRUE(info)) {
    info <- message_not_thrown_before("mdro")
  }

  if (interactive() && isTRUE(verbose) && isTRUE(info)) {
    txt <- paste0(
      "WARNING: In Verbose mode, the mdro() function does not return the MDRO results, but instead returns a data set in logbook form with extensive info about which isolates would be MDRO-positive, or why they are not.",
      "\n\nThis may overwrite your existing data if you use e.g.:",
      "\ndata <- mdro(data, verbose = TRUE)\n\nDo you want to continue?"
    )
    showQuestion <- import_fn("showQuestion", "rstudioapi", error_on_fail = FALSE)
    if (!is.null(showQuestion)) {
      q_continue <- showQuestion("Using verbose = TRUE with mdro()", txt)
    } else {
      q_continue <- utils::menu(choices = c("OK", "Cancel"), graphics = FALSE, title = txt)
    }
    if (q_continue %in% c(FALSE, 2)) {
      message_("Cancelled, returning original data", add_fn = font_red, as_note = FALSE)
      return(x)
    }
  }

  group_msg <- ""
  if (isTRUE(info.bak)) {
    # print group name if used in dplyr::group_by()
    cur_group <- import_fn("cur_group", "dplyr", error_on_fail = FALSE)
    if (!is.null(cur_group)) {
      group_df <- tryCatch(cur_group(), error = function(e) data.frame())
      if (NCOL(group_df) > 0) {
        # transform factors to characters
        group <- vapply(FUN.VALUE = character(1), group_df, function(x) {
          if (is.numeric(x)) {
            format(x)
          } else if (is.logical(x)) {
            as.character(x)
          } else {
            paste0('"', x, '"')
          }
        })
        group_msg <- paste0("\nGroup: ", paste0(names(group), " = ", group, collapse = ", "), "\n")
      }
    }
  }

  # force regular [data.frame], not a tibble or data.table
  x <- as.data.frame(x, stringsAsFactors = FALSE)

  if (pct_required_classes > 1) {
    # allow pct_required_classes = 75 -> pct_required_classes = 0.75
    pct_required_classes <- pct_required_classes / 100
  }

  guideline.bak <- guideline
  if (is.list(guideline)) {
    # Custom MDRO guideline ---------------------------------------------------
    stop_ifnot(inherits(guideline, "custom_mdro_guideline"), "use `custom_mdro_guideline()` to create custom guidelines")
    if (isTRUE(info)) {
      txt <- paste0(
        "Determining MDROs based on custom rules",
        ifelse(isTRUE(attributes(guideline)$as_factor),
          paste0(", resulting in factor levels: ", paste0(attributes(guideline)$values, collapse = " < ")),
          ""
        ),
        "."
      )
      txt <- word_wrap(txt)
      cat(txt, "\n", sep = "")
    }
    x <- run_custom_mdro_guideline(df = x, guideline = guideline, info = info)
    if (isTRUE(info.bak)) {
      cat(group_msg)
      if (sum(!is.na(x$MDRO)) == 0) {
        cat(word_wrap(font_bold(paste0("=> Found 0 MDROs since no isolates are covered by the custom guideline"))))
      } else {
        cat(word_wrap(font_bold(paste0(
          "=> Found ", sum(x$MDRO != "Negative", na.rm = TRUE),
          " custom defined MDROs out of ", sum(!is.na(x$MDRO)),
          " isolates (",
          trimws(percentage(sum(x$MDRO != "Negative", na.rm = TRUE) / sum(!is.na(x$MDRO)))),
          ")\n"
        ))))
      }
    }

    if (isTRUE(verbose)) {
      x$reason[is.na(x$reason)] <- "not covered by guideline"
      x$microorganism <- NA_character_
      x$guideline <- "Custom guideline"
      return(x[, c(
        "row_number",
        "microorganism",
        "MDRO",
        "reason",
        "all_nonsusceptible_columns",
        "guideline"
      ),
      drop = FALSE
      ])
    } else {
      return(x$MDRO)
    }
  } # end of custom MDRO guideline

  guideline <- tolower(gsub("[^a-zA-Z0-9.]+", "", guideline))
  if (is.null(guideline)) {
    # default to the paper by Magiorakos et al. (2012)
    guideline <- "cmi2012"
  }
  if (guideline == "eucast") {
    # turn into latest EUCAST guideline
    guideline <- "eucast3.3"
  }
  if (guideline %in% c("nl", "brmo")) {
    # turn into latest BRMO guideline
    guideline <- "brmo2024"
  }
  if (guideline == "de") {
    guideline <- "mrgn"
  }
  stop_ifnot(
    guideline %in% c("brmo2017", "brmo2024", "mrgn", "eucast3.1", "eucast3.2", "eucast3.3", "tb", "cmi2012"),
    "invalid guideline: ", guideline.bak
  )
  guideline <- list(code = guideline)

  # try to find columns based on type
  # -- mo
  if (is.null(col_mo)) {
    col_mo <- search_type_in_df(x = x, type = "mo", info = info)
  }
  if (is.null(col_mo) && guideline$code == "tb") {
    message_(
      "No column found as input for `col_mo`, ",
      font_bold(paste0("assuming all rows contain ", font_italic("Mycobacterium tuberculosis"), "."))
    )
    x$mo <- as.mo("Mycobacterium tuberculosis", keep_synonyms = TRUE)
    col_mo <- "mo"
  }
  stop_if(is.null(col_mo), "`col_mo` must be set")

  if (guideline$code == "cmi2012") {
    guideline$name <- "Multidrug-resistant, extensively drug-resistant and pandrug-resistant bacteria: an international expert proposal for interim standard definitions for acquired resistance."
    guideline$author <- "Magiorakos AP, Srinivasan A, Carey RB, ..., Vatopoulos A, Weber JT, Monnet DL"
    guideline$version <- NA_character_
    guideline$source_url <- paste0("Clinical Microbiology and Infection 18:3, 2012; ", font_url("https://doi.org/10.1111/j.1469-0691.2011.03570.x", "doi: 10.1111/j.1469-0691.2011.03570.x"))
    guideline$type <- "MDRs/XDRs/PDRs"
  } else if (guideline$code == "eucast3.1") {
    guideline$name <- "EUCAST Expert Rules, \"Intrinsic Resistance and Exceptional Phenotypes Tables\""
    guideline$author <- "EUCAST (European Committee on Antimicrobial Susceptibility Testing)"
    guideline$version <- "3.1, 2016"
    guideline$source_url <- font_url("https://www.eucast.org/fileadmin/src/media/PDFs/EUCAST_files/Expert_Rules/Expert_rules_intrinsic_exceptional_V3.1.pdf", "Direct download")
    guideline$type <- "EUCAST Exceptional Phenotypes"
  } else if (guideline$code == "eucast3.2") {
    guideline$name <- "EUCAST Expert Rules, \"Intrinsic Resistance and Unusual Phenotypes\""
    guideline$author <- "EUCAST (European Committee on Antimicrobial Susceptibility Testing)"
    guideline$version <- "3.2, February 2020"
    guideline$source_url <- font_url("https://www.eucast.org/fileadmin/src/media/PDFs/EUCAST_files/Expert_Rules/2020/Intrinsic_Resistance_and_Unusual_Phenotypes_Tables_v3.2_20200225.pdf", "Direct download")
    guideline$type <- "EUCAST Unusual Phenotypes"
  } else if (guideline$code == "eucast3.3") {
    guideline$name <- "EUCAST Expert Rules, \"Intrinsic Resistance and Unusual Phenotypes\""
    guideline$author <- "EUCAST (European Committee on Antimicrobial Susceptibility Testing)"
    guideline$version <- "3.3, October 2021"
    guideline$source_url <- font_url("https://www.eucast.org/fileadmin/src/media/PDFs/EUCAST_files/Expert_Rules/2021/Intrinsic_Resistance_and_Unusual_Phenotypes_Tables_v3.3_20211018.pdf", "Direct download")
    guideline$type <- "EUCAST Unusual Phenotypes"
  } else if (guideline$code == "tb") {
    guideline$name <- "Companion handbook to the WHO guidelines for the programmatic management of drug-resistant tuberculosis"
    guideline$author <- "WHO (World Health Organization)"
    guideline$version <- "WHO/HTM/TB/2014.11, 2014"
    guideline$source_url <- font_url("https://www.who.int/publications/i/item/9789241548809", "Direct download")
    guideline$type <- "MDR-TB's"

    # support per country:
  } else if (guideline$code == "mrgn") {
    guideline$name <- "Cross-border comparison of the Dutch and German guidelines on multidrug-resistant Gram-negative microorganisms"
    guideline$author <- "M\u00fcller J, Voss A, K\u00f6ck R, ..., Kern WV, Wendt C, Friedrich AW"
    guideline$version <- NA_character_
    guideline$source_url <- paste0("Antimicrobial Resistance and Infection Control 4:7, 2015; ", font_url("https://doi.org/10.1186/s13756-015-0047-6", "doi: 10.1186/s13756-015-0047-6"))
    guideline$type <- "MRGNs"
  } else if (guideline$code == "brmo2024") {
    combine_SI <- TRUE # I must not be considered resistant
    guideline$name <- "Bijzonder Resistente Micro-organismen (BRMO)"
    guideline$author <- "Samenwerkingsverband Richtlijnen Infectiepreventie (SRI)"
    guideline$version <- "November 2024"
    guideline$source_url <- font_url("https://www.sri-richtlijnen.nl/brmo", "Direct link")
    guideline$type <- "BRMOs"
  } else if (guideline$code == "brmo2017") {
    guideline$name <- "WIP-Richtlijn Bijzonder Resistente Micro-organismen (BRMO)"
    guideline$author <- "RIVM (Rijksinstituut voor de Volksgezondheid)"
    guideline$version <- "Last revision (December 2017) - since 2024 superseded by SRI guideline"
    guideline$source_url <- NA_character_
    guideline$type <- "BRMOs"
  } else {
    stop("This guideline is currently unsupported: ", guideline$code, call. = FALSE)
  }

  if (guideline$code == "cmi2012") {
    cols_ab <- get_column_abx(
      x = x,
      soft_dependencies = c(
        # [table] 1 (S aureus)
        "GEN", "RIF", "CPT", "OXA", "CIP", "MFX", "SXT", "FUS", "VAN", "TEC", "TLV", "TGC", "CLI", "DAP", "ERY", "LNZ", "CHL", "FOS", "QDA", "TCY", "DOX", "MNO",
        # [table] 2 (Enterococcus)
        "GEH", "STH", "IPM", "MEM", "DOR", "CIP", "LVX", "MFX", "VAN", "TEC", "TGC", "DAP", "LNZ", "AMP", "QDA", "DOX", "MNO",
        # [table] 3 (Enterobacteriaceae)
        "GEN", "TOB", "AMK", "NET", "CPT", "TCC", "TZP", "ETP", "IPM", "MEM", "DOR", "CZO", "CXM", "CTX", "CAZ", "FEP", "FOX", "CTT", "CIP", "SXT", "TGC", "ATM", "AMP", "AMC", "SAM", "CHL", "FOS", "COL", "TCY", "DOX", "MNO",
        # [table] 4 (Pseudomonas)
        "GEN", "TOB", "AMK", "NET", "IPM", "MEM", "DOR", "CAZ", "FEP", "CIP", "LVX", "TCC", "TZP", "ATM", "FOS", "COL", "PLB",
        # [table] 5 (Acinetobacter)
        "GEN", "TOB", "AMK", "NET", "IPM", "MEM", "DOR", "CIP", "LVX", "TZP", "TCC", "CTX", "CRO", "CAZ", "FEP", "SXT", "SAM", "COL", "PLB", "TCY", "DOX", "MNO"
      ),
      verbose = verbose,
      info = info,
      only_sir_columns = only_sir_columns,
      fn = "mdro",
      ...
    )
  } else if (guideline$code == "eucast3.2") {
    cols_ab <- get_column_abx(
      x = x,
      soft_dependencies = c("AMP", "AMX", "CIP", "DAL", "DAP", "ERV", "FDX", "GEN", "LNZ", "MEM", "MTR", "OMC", "ORI", "PEN", "QDA", "RIF", "TEC", "TGC", "TLV", "TOB", "TZD", "VAN"),
      verbose = verbose,
      info = info,
      only_sir_columns = only_sir_columns,
      fn = "mdro",
      ...
    )
  } else if (guideline$code == "eucast3.3") {
    cols_ab <- get_column_abx(
      x = x,
      soft_dependencies = c("AMP", "AMX", "CIP", "DAL", "DAP", "ERV", "FDX", "GEN", "LNZ", "MEM", "MTR", "OMC", "ORI", "PEN", "QDA", "RIF", "TEC", "TGC", "TLV", "TOB", "TZD", "VAN"),
      verbose = verbose,
      info = info,
      only_sir_columns = only_sir_columns,
      fn = "mdro",
      ...
    )
  } else if (guideline$code == "tb") {
    cols_ab <- get_column_abx(
      x = x,
      soft_dependencies = c("CAP", "ETH", "GAT", "INH", "PZA", "RIF", "RIB", "RFP"),
      verbose = verbose,
      info = info,
      only_sir_columns = only_sir_columns,
      fn = "mdro",
      ...
    )
  } else if (guideline$code == "brmo") {
    # Dutch 2024 guideline
    cols_ab <- get_column_abx(
      x = x,
      soft_dependencies = c("SXT", "GEN", "TOB", "AMK", "IPM", "MEM", "CIP", "LVX", "NOR", "PIP", "CAZ", "VAN", "PEN", "AMX", "AMP", "FLC", "OXA", "FOX", "FOX1"),
      verbose = verbose,
      info = info,
      only_sir_columns = only_sir_columns,
      fn = "mdro",
      ...
    )
  } else if (guideline$code == "mrgn") {
    cols_ab <- get_column_abx(
      x = x,
      soft_dependencies = c("PIP", "CTX", "CAZ", "IPM", "MEM", "CIP"),
      verbose = verbose,
      info = info,
      only_sir_columns = only_sir_columns,
      fn = "mdro",
      ...
    )
  } else {
    cols_ab <- get_column_abx(
      x = x,
      verbose = verbose,
      info = info,
      only_sir_columns = only_sir_columns,
      fn = "mdro",
      ...
    )
  }
  if (!"AMP" %in% names(cols_ab) && "AMX" %in% names(cols_ab)) {
    # ampicillin column is missing, but amoxicillin is available
    if (isTRUE(info)) {
      message_("Using column '", cols_ab[names(cols_ab) == "AMX"], "' as input for ampicillin since many MDRO rules depend on it.", add_fn = font_red)
    }
    cols_ab <- c(cols_ab, c(AMP = unname(cols_ab[names(cols_ab) == "AMX"])))
  }
  cols_ab <- cols_ab[!duplicated(cols_ab)]

  # nolint start
  AMC <- cols_ab["AMC"]
  AMK <- cols_ab["AMK"]
  AMP <- cols_ab["AMP"]
  AMX <- cols_ab["AMX"]
  ATM <- cols_ab["ATM"]
  AZL <- cols_ab["AZL"]
  AZM <- cols_ab["AZM"]
  BPR <- cols_ab["BPR"]
  CAC <- cols_ab["CAC"]
  CAT <- cols_ab["CAT"]
  CAZ <- cols_ab["CAZ"]
  CCV <- cols_ab["CCV"]
  CDR <- cols_ab["CDR"]
  CDZ <- cols_ab["CDZ"]
  CEC <- cols_ab["CEC"]
  CED <- cols_ab["CED"]
  CEI <- cols_ab["CEI"]
  CEP <- cols_ab["CEP"]
  CFM <- cols_ab["CFM"]
  CFM1 <- cols_ab["CFM1"]
  CFP <- cols_ab["CFP"]
  CFR <- cols_ab["CFR"]
  CFS <- cols_ab["CFS"]
  CHL <- cols_ab["CHL"]
  CID <- cols_ab["CID"]
  CIP <- cols_ab["CIP"]
  CLI <- cols_ab["CLI"]
  CLR <- cols_ab["CLR"]
  CMX <- cols_ab["CMX"]
  CMZ <- cols_ab["CMZ"]
  CND <- cols_ab["CND"]
  COL <- cols_ab["COL"]
  CPD <- cols_ab["CPD"]
  CPM <- cols_ab["CPM"]
  CPO <- cols_ab["CPO"]
  CPR <- cols_ab["CPR"]
  CPT <- cols_ab["CPT"]
  CRD <- cols_ab["CRD"]
  CRO <- cols_ab["CRO"]
  CSL <- cols_ab["CSL"]
  CTB <- cols_ab["CTB"]
  CTF <- cols_ab["CTF"]
  CTL <- cols_ab["CTL"]
  CTT <- cols_ab["CTT"]
  CTX <- cols_ab["CTX"]
  CTZ <- cols_ab["CTZ"]
  CXM <- cols_ab["CXM"]
  CZA <- cols_ab["CZA"]
  CZD <- cols_ab["CZD"]
  CZO <- cols_ab["CZO"]
  CZX <- cols_ab["CZX"]
  DAL <- cols_ab["DAL"]
  DAP <- cols_ab["DAP"]
  DIT <- cols_ab["DIT"]
  DIZ <- cols_ab["DIZ"]
  DOR <- cols_ab["DOR"]
  DOX <- cols_ab["DOX"]
  ENX <- cols_ab["ENX"]
  ERV <- cols_ab["ERV"]
  ERY <- cols_ab["ERY"]
  ETP <- cols_ab["ETP"]
  FDX <- cols_ab["FDX"]
  FEP <- cols_ab["FEP"]
  FLC <- cols_ab["FLC"]
  FLE <- cols_ab["FLE"]
  FOS <- cols_ab["FOS"]
  FOX <- cols_ab["FOX"]
  FOX1 <- cols_ab["FOX1"]
  FUS <- cols_ab["FUS"]
  GAT <- cols_ab["GAT"]
  GEH <- cols_ab["GEH"]
  GEM <- cols_ab["GEM"]
  GEN <- cols_ab["GEN"]
  GRX <- cols_ab["GRX"]
  HAP <- cols_ab["HAP"]
  IPM <- cols_ab["IPM"]
  KAN <- cols_ab["KAN"]
  LEX <- cols_ab["LEX"]
  LIN <- cols_ab["LIN"]
  LNZ <- cols_ab["LNZ"]
  LOM <- cols_ab["LOM"]
  LOR <- cols_ab["LOR"]
  LTM <- cols_ab["LTM"]
  LVX <- cols_ab["LVX"]
  MAN <- cols_ab["MAN"]
  MEM <- cols_ab["MEM"]
  MEV <- cols_ab["MEV"]
  MEZ <- cols_ab["MEZ"]
  MFX <- cols_ab["MFX"]
  MNO <- cols_ab["MNO"]
  MTR <- cols_ab["MTR"]
  NAL <- cols_ab["NAL"]
  NEO <- cols_ab["NEO"]
  NET <- cols_ab["NET"]
  NIT <- cols_ab["NIT"]
  NOR <- cols_ab["NOR"]
  NOV <- cols_ab["NOV"]
  OFX <- cols_ab["OFX"]
  OMC <- cols_ab["OMC"]
  ORI <- cols_ab["ORI"]
  OXA <- cols_ab["OXA"]
  PAZ <- cols_ab["PAZ"]
  PEF <- cols_ab["PEF"]
  PEN <- cols_ab["PEN"]
  PIP <- cols_ab["PIP"]
  PLB <- cols_ab["PLB"]
  PRI <- cols_ab["PRI"]
  PRU <- cols_ab["PRU"]
  QDA <- cols_ab["QDA"]
  RFL <- cols_ab["RFL"]
  RID <- cols_ab["RID"]
  RIF <- cols_ab["RIF"]
  RXT <- cols_ab["RXT"]
  SAM <- cols_ab["SAM"]
  SIS <- cols_ab["SIS"]
  SPT <- cols_ab["SPT"]
  SPX <- cols_ab["SPX"]
  STH <- cols_ab["STH"]
  SXT <- cols_ab["SXT"]
  TCC <- cols_ab["TCC"]
  TCY <- cols_ab["TCY"]
  TEC <- cols_ab["TEC"]
  TGC <- cols_ab["TGC"]
  TIC <- cols_ab["TIC"]
  TLV <- cols_ab["TLV"]
  TMP <- cols_ab["TMP"]
  TMX <- cols_ab["TMX"]
  TOB <- cols_ab["TOB"]
  TVA <- cols_ab["TVA"]
  TZD <- cols_ab["TZD"]
  TZP <- cols_ab["TZP"]
  VAN <- cols_ab["VAN"]
  # additional for TB
  CAP <- cols_ab["CAP"]
  ETH <- cols_ab["ETH"]
  GAT <- cols_ab["GAT"]
  INH <- cols_ab["INH"]
  PZA <- cols_ab["PZA"]
  RIF <- cols_ab["RIF"]
  RIB <- cols_ab["RIB"]
  RFP <- cols_ab["RFP"]
  abx_tb <- c(CAP, ETH, GAT, INH, PZA, RIF, RIB, RFP)
  abx_tb <- abx_tb[!is.na(abx_tb)]
  stop_if(guideline$code == "tb" & length(abx_tb) == 0, "no antimycobacterials found in data set")
  # nolint end
  if (isTRUE(combine_SI)) {
    search_result <- "R"
  } else {
    search_result <- c("R", "I")
  }

  if (isTRUE(info)) {
    if (isTRUE(combine_SI)) {
      cat(font_red("\nOnly results with 'R' are considered as resistance. Use `combine_SI = FALSE` to also consider 'I' as resistance.\n"))
    } else {
      cat(font_red("\nResults with 'R' or 'I' are considered as resistance. Use `combine_SI = TRUE` to only consider 'R' as resistance.\n"))
    }
    cat("\n", word_wrap("Determining multidrug-resistant organisms (MDRO), according to:"), "\n",
      word_wrap(paste0(font_bold("Guideline: "), font_italic(guideline$name)), extra_indent = 11, as_note = FALSE), "\n",
      word_wrap(paste0(font_bold("Author(s): "), guideline$author), extra_indent = 11, as_note = FALSE), "\n",
      ifelse(!is.na(guideline$version),
        paste0(word_wrap(paste0(font_bold("Version:   "), guideline$version), extra_indent = 11, as_note = FALSE), "\n"),
        ""
      ),
      paste0(font_bold("Source:    "), guideline$source_url),
      "\n\n",
      sep = ""
    )
  }

  ab_missing <- function(ab) {
    isTRUE(ab %in% c(NULL, NA)) | length(ab) == 0
  }
  ab_NA <- function(x) {
    x[!is.na(x)]
  }
  try_ab <- function(expr) {
    out <- tryCatch(expr, error = function(e) FALSE)
    out[is.na(out)] <- FALSE
    out
  }
  col_values <- function(df, col, return_if_lacking = "") {
    if (col %in% colnames(df)) {
      df[[col]]
    } else {
      rep(return_if_lacking, NROW(df))
    }
  }
  NA_as_FALSE <- function(x) {
    x[is.na(x)] <- FALSE
    x
  }

  # antimicrobial classes
  # nolint start
  aminoglycosides <- c(TOB, GEN)
  cephalosporins <- c(CDZ, CAC, CEC, CFR, RID, MAN, CTZ, CZD, CZO, CDR, DIT, FEP, CAT, CFM, CMX, CMZ, DIZ, CID, CFP, CSL, CND, CTX, CTT, CTF, FOX, CPM, CPO, CPD, CPR, CRD, CFS, CPT, CAZ, CCV, CTL, CTB, CZX, BPR, CFM1, CEI, CRO, CXM, LEX, CEP, HAP, CED, LTM, LOR)
  cephalosporins_1st <- c(CAC, CFR, RID, CTZ, CZD, CZO, CRD, CTL, LEX, CEP, HAP, CED)
  cephalosporins_2nd <- c(CEC, MAN, CMZ, CID, CND, CTT, CTF, FOX, CPR, CXM, LOR)
  cephalosporins_3rd <- c(CDZ, CDR, DIT, CAT, CFM, CMX, DIZ, CFP, CSL, CTX, CPM, CPD, CFS, CAZ, CCV, CTB, CZX, CRO, LTM)
  carbapenems <- c(DOR, ETP, IPM, MEM, MEV)
  fluoroquinolones <- c(CIP, ENX, FLE, GAT, GEM, GRX, LVX, LOM, MFX, NOR, OFX, PAZ, PEF, PRU, RFL, SPX, TMX, TVA)
  # nolint end

  # helper function for editing the table
  trans_tbl <- function(to, rows, cols, any_all, reason = NULL) {
    cols.bak <- cols
    if (identical(cols, "any")) {
      cols <- unique(cols_ab)
    }
    cols <- cols[!ab_missing(cols)]
    cols <- cols[!is.na(cols)]
    if (length(rows) > 0 && length(cols) > 0) {
      x[, cols] <- as.data.frame(
        lapply(
          x[, cols, drop = FALSE],
          function(col) as.sir(col)
        ),
        stringsAsFactors = FALSE
      )
      x[rows, "all_nonsusceptible_columns"] <<- vapply(
        FUN.VALUE = character(1),
        rows,
        function(row, group_vct = cols_ab) {
          cols_nonsus <- vapply(
            FUN.VALUE = logical(1),
            x[row, group_vct, drop = FALSE],
            function(y) y %in% search_result
          )
          paste(
            unique(sort(c(
              unlist(strsplit(x[row, "all_nonsusceptible_columns", drop = TRUE], ", ", fixed = TRUE)),
              names(cols_nonsus)[cols_nonsus]
            ))),
            collapse = ", "
          )
        }
      )

      if (any_all == "any") {
        search_function <- any
      } else if (any_all == "all") {
        search_function <- all
      }
      x_transposed <- as.list(as.data.frame(t(x[, cols, drop = FALSE]),
        stringsAsFactors = FALSE
      ))
      rows_affected <- vapply(
        FUN.VALUE = logical(1),
        x_transposed,
        function(y) search_function(y %in% search_result, na.rm = TRUE) | identical(cols.bak, "any")
      )
      rows_affected <- x[which(rows_affected), "row_number", drop = TRUE]
      rows_to_change <- rows[rows %in% rows_affected]
      rows_not_to_change <- rows[!rows %in% c(rows_affected, rows_to_change)]
      rows_not_to_change <- rows_not_to_change[is.na(x[rows_not_to_change, "reason"])]
      if (is.null(reason)) {
        reason <- paste0(
          any_all,
          " of the required antimicrobials ",
          ifelse(any_all == "any", "is", "are"),
          " R",
          ifelse(!isTRUE(combine_SI), " or I", "")
        )
      }
      x[rows_to_change, "MDRO"] <<- to
      x[rows_to_change, "reason"] <<- reason
      x[rows_not_to_change, "reason"] <<- "guideline criteria not met"
    }
  }

  trans_tbl2 <- function(txt, rows, lst) {
    if (isTRUE(info)) {
      message_(txt, "...", appendLF = FALSE, as_note = FALSE)
    }
    if (length(rows) > 0) {
      # function specific for the CMI paper of 2012 (Magiorakos et al.)
      lst_vector <- unlist(lst)[!is.na(unlist(lst))]
      # keep only unique ones:
      lst_vector <- lst_vector[!duplicated(paste(lst_vector, names(lst_vector)))]

      x[, lst_vector] <- as.data.frame(
        lapply(
          x[, lst_vector, drop = FALSE],
          function(col) as.sir(col)
        ),
        stringsAsFactors = FALSE
      )
      x[rows, "classes_in_guideline"] <<- length(lst)
      x[rows, "classes_available"] <<- vapply(
        FUN.VALUE = double(1),
        rows,
        function(row, group_tbl = lst) {
          sum(vapply(
            FUN.VALUE = logical(1),
            group_tbl,
            function(group) any(unlist(x[row, group[!is.na(group)], drop = TRUE]) %in% c("S", "SDD", "I", "R"))
          ))
        }
      )

      if (isTRUE(verbose)) {
        x[rows, "all_nonsusceptible_columns"] <<- vapply(
          FUN.VALUE = character(1),
          rows,
          function(row, group_vct = lst_vector) {
            cols_nonsus <- vapply(FUN.VALUE = logical(1), x[row, group_vct, drop = FALSE], function(y) y %in% search_result)
            paste(unique(sort(names(cols_nonsus)[cols_nonsus])), collapse = ", ")
          }
        )
      }
      x[rows, "classes_affected"] <<- vapply(
        FUN.VALUE = double(1),
        rows,
        function(row, group_tbl = lst) {
          sum(
            vapply(
              FUN.VALUE = logical(1),
              group_tbl,
              function(group) {
                any(unlist(x[row, group[!is.na(group)], drop = TRUE]) %in% search_result, na.rm = TRUE)
              }
            ),
            na.rm = TRUE
          )
        }
      )
      # for PDR; all drugs are R (or I if combine_SI = FALSE)
      x_transposed <- as.list(as.data.frame(t(x[rows, lst_vector, drop = FALSE]),
        stringsAsFactors = FALSE
      ))
      row_filter <- vapply(FUN.VALUE = logical(1), x_transposed, function(y) all(y %in% search_result, na.rm = TRUE))
      x[which(row_filter), "classes_affected"] <<- 999
    }

    if (isTRUE(info)) {
      message_(" OK.", add_fn = list(font_green, font_bold), as_note = FALSE)
    }
  }

  x[, col_mo] <- as.mo(as.character(x[, col_mo, drop = TRUE]))
  # rename col_mo to prevent interference with joined columns
  colnames(x)[colnames(x) == col_mo] <- ".col_mo"
  col_mo <- ".col_mo"
  # join to microorganisms data set
  x <- left_join_microorganisms(x, by = col_mo)
  x$MDRO <- ifelse(!is.na(x$genus), 1, NA_integer_)
  x$row_number <- seq_len(nrow(x))
  x$reason <- NA_character_
  x$all_nonsusceptible_columns <- ""

  if (guideline$code == "cmi2012") {
    # CMI, 2012 ---------------------------------------------------------------
    # Non-susceptible = R and I
    # (see header 'Approaches to Creating Definitions for MDR, XDR and PDR' in paper)

    # take amoxicillin if ampicillin is unavailable
    if (is.na(AMP) && !is.na(AMX)) {
      if (isTRUE(verbose)) {
        message_("Filling ampicillin (AMP) results with amoxicillin (AMX) results")
      }
      AMP <- AMX
    }
    # take ceftriaxone if cefotaxime is unavailable and vice versa
    if (is.na(CRO) && !is.na(CTX)) {
      if (isTRUE(verbose)) {
        message_("Filling ceftriaxone (CRO) results with cefotaxime (CTX) results")
      }
      CRO <- CTX
    }
    if (is.na(CTX) && !is.na(CRO)) {
      if (isTRUE(verbose)) {
        message_("Filling cefotaxime (CTX) results with ceftriaxone (CRO) results")
      }
      CTX <- CRO
    }

    # intrinsic resistant must not be considered for the determination of MDR,
    # so let's just remove them, meticulously following the paper
    x[which(x$genus == "Enterococcus" & x$species == "faecium"), ab_NA(IPM)] <- NA
    x[which(x$genus == "Enterococcus" & x$species == "faecalis"), ab_NA(QDA)] <- NA
    x[which((x$genus == "Providencia" & x$species == "rettgeri") |
      (x$genus == "Providencia" & x$species == "stuartii")), ab_NA(c(GEN, TOB, NET))] <- NA
    x[which(x$genus == "Escherichia" & x$species == "hermannii"), ab_NA(c(TCC, TZP))] <- NA
    x[which((x$genus == "Citrobacter" & x$species == "freundii") |
      (x$genus == "Enterobacter" & x$species == "aerogenes") |
      (x$genus == "Klebsiella" & x$species == "aerogenes") # new name (2017)
    | (x$genus == "Enterobacter" & x$species == "cloacae") |
      (x$genus == "Hafnia" & x$species == "alvei") |
      (x$genus == "Morganella" & x$species == "morganii") |
      (x$genus == "Proteus" & x$species == "penneri") |
      (x$genus == "Proteus" & x$species == "vulgaris") |
      (x$genus == "Serratia" & x$species == "marcescens")), ab_NA(CZO)] <- NA
    x[which((x$genus == "Morganella" & x$species == "morganii") |
      (x$genus == "Proteus" & x$species == "penneri") |
      (x$genus == "Proteus" & x$species == "vulgaris") |
      (x$genus == "Serratia" & x$species == "marcescens")), ab_NA(CXM)] <- NA
    x[which((x$genus == "Morganella" & x$species == "morganii") |
      (x$genus == "Proteus" & x$species == "mirabilis") |
      (x$genus == "Proteus" & x$species == "penneri") |
      (x$genus == "Proteus" & x$species == "vulgaris") |
      (x$genus == "Providencia" & x$species == "rettgeri") |
      (x$genus == "Providencia" & x$species == "stuartii")), ab_NA(TGC)] <- NA
    x[which((x$genus == "Citrobacter" & x$species == "koseri") |
      (x$genus == "Citrobacter" & x$species == "freundii") |
      (x$genus == "Enterobacter" & x$species == "aerogenes") |
      (x$genus == "Klebsiella" & x$species == "aerogenes") # new name (2017)
    | (x$genus == "Enterobacter" & x$species == "cloacae") |
      (x$genus == "Escherichia" & x$species == "hermannii") |
      (x$genus == "Hafnia" & x$species == "alvei") |
      (x$genus == "Klebsiella") |
      (x$genus == "Morganella" & x$species == "morganii") |
      (x$genus == "Proteus" & x$species == "penneri") |
      (x$genus == "Proteus" & x$species == "vulgaris") |
      (x$genus == "Providencia" & x$species == "rettgeri") |
      (x$genus == "Providencia" & x$species == "stuartii") |
      (x$genus == "Serratia" & x$species == "marcescens")), ab_NA(AMP)] <- NA
    x[which((x$genus == "Citrobacter" & x$species == "freundii") |
      (x$genus == "Enterobacter" & x$species == "aerogenes") |
      (x$genus == "Klebsiella" & x$species == "aerogenes") # new name (2017)
    | (x$genus == "Enterobacter" & x$species == "cloacae") |
      (x$genus == "Hafnia" & x$species == "alvei") |
      (x$genus == "Morganella" & x$species == "morganii") |
      (x$genus == "Providencia" & x$species == "rettgeri") |
      (x$genus == "Providencia" & x$species == "stuartii") |
      (x$genus == "Serratia" & x$species == "marcescens")), ab_NA(AMC)] <- NA
    x[which((x$genus == "Citrobacter" & x$species == "freundii") |
      (x$genus == "Citrobacter" & x$species == "koseri") |
      (x$genus == "Enterobacter" & x$species == "aerogenes") |
      (x$genus == "Klebsiella" & x$species == "aerogenes") # new name (2017)
    | (x$genus == "Enterobacter" & x$species == "cloacae") |
      (x$genus == "Hafnia" & x$species == "alvei") |
      (x$genus == "Providencia" & x$species == "rettgeri") |
      (x$genus == "Serratia" & x$species == "marcescens")), ab_NA(SAM)] <- NA
    x[which((x$genus == "Morganella" & x$species == "morganii") |
      (x$genus == "Proteus" & x$species == "mirabilis") |
      (x$genus == "Proteus" & x$species == "penneri") |
      (x$genus == "Proteus" & x$species == "vulgaris") |
      (x$genus == "Providencia" & x$species == "rettgeri") |
      (x$genus == "Providencia" & x$species == "stuartii") |
      (x$genus == "Serratia" & x$species == "marcescens")), ab_NA(COL)] <- NA
    x[which((x$genus == "Morganella" & x$species == "morganii") |
      (x$genus == "Proteus" & x$species == "mirabilis") |
      (x$genus == "Proteus" & x$species == "penneri") |
      (x$genus == "Proteus" & x$species == "vulgaris") |
      (x$genus == "Providencia" & x$species == "rettgeri") |
      (x$genus == "Providencia" & x$species == "stuartii")), ab_NA(TCY)] <- NA
    x[which((x$genus == "Morganella" & x$species == "morganii") |
      (x$genus == "Proteus" & x$species == "penneri") |
      (x$genus == "Proteus" & x$species == "vulgaris") |
      (x$genus == "Providencia" & x$species == "rettgeri") |
      (x$genus == "Providencia" & x$species == "stuartii")), ab_NA(c(DOX, MNO))] <- NA

    x$classes_in_guideline <- NA_integer_
    x$classes_available <- NA_integer_
    x$classes_affected <- NA_integer_

    # now add the MDR levels to the data
    trans_tbl(
      2,
      which(x$genus == "Staphylococcus" & x$species == "aureus"),
      c(OXA, FOX),
      "any"
    )
    trans_tbl2(
      paste("Table 1 -", font_italic("Staphylococcus aureus")),
      which(x$genus == "Staphylococcus" & x$species == "aureus"),
      list(
        GEN,
        RIF,
        CPT,
        c(OXA, FOX),
        c(CIP, MFX),
        SXT,
        FUS,
        c(VAN, TEC, TLV),
        TGC,
        CLI,
        DAP,
        ERY,
        LNZ,
        CHL,
        FOS,
        QDA,
        c(TCY, DOX, MNO)
      )
    )
    trans_tbl2(
      paste("Table 2 -", font_italic("Enterococcus"), "spp."),
      which(x$genus == "Enterococcus"),
      list(
        GEH,
        STH,
        c(IPM, MEM, DOR),
        c(CIP, LVX, MFX),
        c(VAN, TEC),
        TGC,
        DAP,
        LNZ,
        AMP,
        QDA,
        c(DOX, MNO)
      )
    )
    trans_tbl2(
      paste0("Table 3 - ", font_italic("Enterobacteriaceae")),
      # this new order was previously 'Enterobacteriales' and contained only the family 'Enterobacteriaceae':
      which(x$order == "Enterobacterales"),
      list(
        c(GEN, TOB, AMK, NET),
        CPT,
        c(TCC, TZP),
        c(ETP, IPM, MEM, DOR),
        CZO,
        CXM,
        c(CTX, CAZ, FEP),
        c(FOX, CTT),
        CIP,
        SXT,
        TGC,
        ATM,
        AMP,
        c(AMC, SAM),
        CHL,
        FOS,
        COL,
        c(TCY, DOX, MNO)
      )
    )
    trans_tbl2(
      paste("Table 4 -", font_italic("Pseudomonas aeruginosa")),
      which(x$genus == "Pseudomonas" & x$species == "aeruginosa"),
      list(
        c(GEN, TOB, AMK, NET),
        c(IPM, MEM, DOR),
        c(CAZ, FEP),
        c(CIP, LVX),
        c(TCC, TZP),
        ATM,
        FOS,
        c(COL, PLB)
      )
    )
    trans_tbl2(
      paste("Table 5 -", font_italic("Acinetobacter"), "spp."),
      which(x$genus == "Acinetobacter"),
      list(
        c(GEN, TOB, AMK, NET),
        c(IPM, MEM, DOR),
        c(CIP, LVX),
        c(TZP, TCC),
        c(CTX, CRO, CAZ, FEP),
        SXT,
        SAM,
        c(COL, PLB),
        c(TCY, DOX, MNO)
      )
    )

    # now set MDROs:
    # MDR (=2): >=3 classes affected
    x[which(x$classes_affected >= 3), "MDRO"] <- 2
    if (isTRUE(verbose)) {
      x[which(x$classes_affected >= 3), "reason"] <- paste0(
        "at least 3 classes contain R",
        ifelse(!isTRUE(combine_SI), " or I", ""), ": ",
        x$classes_affected[which(x$classes_affected >= 3)],
        " out of ", x$classes_available[which(x$classes_affected >= 3)], " available classes"
      )
    }

    # XDR (=3): all but <=2 classes affected
    x[which((x$classes_in_guideline - x$classes_affected) <= 2), "MDRO"] <- 3
    if (isTRUE(verbose)) {
      x[which(x$MDRO == 3), "reason"] <- paste0(
        "less than 3 classes remain susceptible (", x$classes_in_guideline[which((x$classes_in_guideline - x$classes_affected) <= 2)] - x$classes_affected[which(x$MDRO == 3)],
        " out of ", x$classes_in_guideline[which(x$MDRO == 3)], " classes)"
      )
    }

    # PDR (=4): all drugs are R
    x[which(x$classes_affected == 999 & x$classes_in_guideline == x$classes_available), "MDRO"] <- 4
    if (isTRUE(verbose)) {
      x[which(x$MDRO == 4), "reason"] <- paste(
        "all antimicrobials in all",
        x$classes_in_guideline[which(x$MDRO == 4)],
        "classes were tested R",
        ifelse(!isTRUE(combine_SI), " or I", "")
      )
    }

    # not enough classes available
    x[which(x$MDRO %in% c(1, 3) & x$classes_available < floor(x$classes_in_guideline * pct_required_classes)), "MDRO"] <- -1
    if (isTRUE(verbose)) {
      x[which(x$MDRO == -1), "reason"] <- paste0(
        "not enough classes available: ", x$classes_available[which(x$MDRO == -1)],
        " of required ", (floor(x$classes_in_guideline * pct_required_classes))[which(x$MDRO == -1)],
        " (~", percentage(pct_required_classes), " of ", x$classes_in_guideline[which(x$MDRO == -1)], ")"
      )
    }

    # add antimicrobial names of resistant ones to verbose output
  }

  if (guideline$code == "eucast3.1") {
    # EUCAST 3.1 --------------------------------------------------------------
    # Table 5
    trans_tbl(
      3,
      which(x$order == "Enterobacterales" |
        (x$genus == "Pseudomonas" & x$species == "aeruginosa") |
        x$genus == "Acinetobacter"),
      COL,
      "all"
    )
    trans_tbl(
      3,
      which(x$genus == "Salmonella" & x$species == "Typhi"),
      c(carbapenems, fluoroquinolones),
      "any"
    )
    trans_tbl(
      3,
      which(x$genus == "Haemophilus" & x$species == "influenzae"),
      c(cephalosporins_3rd, carbapenems, fluoroquinolones),
      "any"
    )
    trans_tbl(
      3,
      which(x$genus == "Moraxella" & x$species == "catarrhalis"),
      c(cephalosporins_3rd, fluoroquinolones),
      "any"
    )
    trans_tbl(
      3,
      which(x$genus == "Neisseria" & x$species == "meningitidis"),
      c(cephalosporins_3rd, fluoroquinolones),
      "any"
    )
    trans_tbl(
      3,
      which(x$genus == "Neisseria" & x$species == "gonorrhoeae"),
      AZM,
      "any"
    )
    # Table 6
    trans_tbl(
      3,
      which(x$fullname %like% "^(Coagulase-negative|Staphylococcus (aureus|epidermidis|hominis|haemolyticus|intermedius|pseudointermedius))"),
      c(VAN, TEC, DAP, LNZ, QDA, TGC),
      "any"
    )
    trans_tbl(
      3,
      which(x$genus == "Corynebacterium"),
      c(VAN, TEC, DAP, LNZ, QDA, TGC),
      "any"
    )
    trans_tbl(
      3,
      which(x$genus == "Streptococcus" & x$species == "pneumoniae"),
      c(carbapenems, VAN, TEC, DAP, LNZ, QDA, TGC, RIF),
      "any"
    )
    trans_tbl(
      3, # Sr. groups A/B/C/G
      which(x$fullname %like% "^Streptococcus (group (A|B|C|G)|pyogenes|agalactiae|equisimilis|equi|zooepidemicus|dysgalactiae|anginosus)"),
      c(PEN, cephalosporins, VAN, TEC, DAP, LNZ, QDA, TGC),
      "any"
    )
    trans_tbl(
      3,
      which(x$genus == "Enterococcus"),
      c(DAP, LNZ, TGC, TEC),
      "any"
    )
    trans_tbl(
      3,
      which(x$genus == "Enterococcus" & x$species == "faecalis"),
      c(AMP, AMX),
      "any"
    )
    # Table 7
    trans_tbl(
      3,
      which(x$genus == "Bacteroides"),
      MTR,
      "any"
    )
    trans_tbl(
      3,
      which(x$genus %in% c("Clostridium", "Clostridioides") & x$species == "difficile"),
      c(MTR, VAN),
      "any"
    )
  }

  if (guideline$code == "eucast3.2") {
    # EUCAST 3.2 --------------------------------------------------------------
    # Table 6
    trans_tbl(
      3,
      which((x$order == "Enterobacterales" &
        !x$family == "Morganellaceae" &
        !(x$genus == "Serratia" & x$species == "marcescens")) |
        (x$genus == "Pseudomonas" & x$species == "aeruginosa") |
        x$genus == "Acinetobacter"),
      COL,
      "all"
    )
    trans_tbl(
      3,
      which(x$genus == "Salmonella" & x$species == "Typhi"),
      c(carbapenems),
      "any"
    )
    trans_tbl(
      3,
      which(x$genus == "Haemophilus" & x$species == "influenzae"),
      c(cephalosporins_3rd, carbapenems, fluoroquinolones),
      "any"
    )
    trans_tbl(
      3,
      which(x$genus == "Moraxella" & x$species == "catarrhalis"),
      c(cephalosporins_3rd, fluoroquinolones),
      "any"
    )
    trans_tbl(
      3,
      which(x$genus == "Neisseria" & x$species == "meningitidis"),
      c(cephalosporins_3rd, fluoroquinolones),
      "any"
    )
    trans_tbl(
      3,
      which(x$genus == "Neisseria" & x$species == "gonorrhoeae"),
      SPT,
      "any"
    )
    # Table 7
    trans_tbl(
      3,
      which(x$genus == "Staphylococcus" & x$species == "aureus"),
      c(VAN, TEC, TLV, DAL, ORI, DAP, LNZ, TZD, QDA, TGC, ERV, OMC),
      "any"
    )
    trans_tbl(
      3,
      which(x$mo %in% MO_CONS), # coagulase-negative Staphylococcus
      c(VAN, TLV, DAL, ORI, DAP, LNZ, TZD, QDA, TGC, ERV, OMC),
      "any"
    )
    trans_tbl(
      3,
      which(x$genus == "Corynebacterium"),
      c(VAN, TEC, TLV, DAL, ORI, DAP, LNZ, TZD, QDA, TGC),
      "any"
    )
    trans_tbl(
      3,
      which(x$genus == "Streptococcus" & x$species == "pneumoniae"),
      c(carbapenems, VAN, TEC, TLV, DAL, ORI, DAP, LNZ, TZD, QDA, TGC, ERV, OMC, RIF),
      "any"
    )
    trans_tbl(
      3, # Sr. groups A/B/C/G
      which(x$mo %in% MO_STREP_ABCG),
      c(PEN, cephalosporins, VAN, TEC, TLV, DAL, ORI, DAP, LNZ, TZD, QDA, TGC, ERV, OMC),
      "any"
    )
    trans_tbl(
      3,
      which(x$genus == "Enterococcus"),
      c(DAP, LNZ, TGC, ERV, OMC, TEC),
      "any"
    )
    trans_tbl(
      3,
      which(x$genus == "Enterococcus" & x$species == "faecalis"),
      c(AMP, AMX),
      "any"
    )
    # Table 8
    trans_tbl(
      3,
      which(x$genus == "Bacteroides"),
      MTR,
      "any"
    )
    trans_tbl(
      3,
      which(x$genus %in% c("Clostridium", "Clostridioides") & x$species == "difficile"),
      c(MTR, VAN, FDX),
      "any"
    )
  }

  if (guideline$code == "eucast3.3") {
    # EUCAST 3.3 --------------------------------------------------------------
    # note: this guideline is equal to EUCAST 3.2 - no MDRO insights changed
    # Table 6
    trans_tbl(
      3,
      which((x$order == "Enterobacterales" &
        !x$family == "Morganellaceae" &
        !(x$genus == "Serratia" & x$species == "marcescens")) |
        (x$genus == "Pseudomonas" & x$species == "aeruginosa") |
        x$genus == "Acinetobacter"),
      COL,
      "all"
    )
    trans_tbl(
      3,
      which(x$genus == "Salmonella" & x$species == "Typhi"),
      c(carbapenems),
      "any"
    )
    trans_tbl(
      3,
      which(x$genus == "Haemophilus" & x$species == "influenzae"),
      c(cephalosporins_3rd, carbapenems, fluoroquinolones),
      "any"
    )
    trans_tbl(
      3,
      which(x$genus == "Moraxella" & x$species == "catarrhalis"),
      c(cephalosporins_3rd, fluoroquinolones),
      "any"
    )
    trans_tbl(
      3,
      which(x$genus == "Neisseria" & x$species == "meningitidis"),
      c(cephalosporins_3rd, fluoroquinolones),
      "any"
    )
    trans_tbl(
      3,
      which(x$genus == "Neisseria" & x$species == "gonorrhoeae"),
      SPT,
      "any"
    )
    # Table 7
    trans_tbl(
      3,
      which(x$genus == "Staphylococcus" & x$species == "aureus"),
      c(VAN, TEC, TLV, DAL, ORI, DAP, LNZ, TZD, QDA, TGC, ERV, OMC),
      "any"
    )
    trans_tbl(
      3,
      which(x$mo %in% MO_CONS), # coagulase-negative Staphylococcus
      c(VAN, TLV, DAL, ORI, DAP, LNZ, TZD, QDA, TGC, ERV, OMC),
      "any"
    )
    trans_tbl(
      3,
      which(x$genus == "Corynebacterium"),
      c(VAN, TEC, TLV, DAL, ORI, DAP, LNZ, TZD, QDA, TGC),
      "any"
    )
    trans_tbl(
      3,
      which(x$genus == "Streptococcus" & x$species == "pneumoniae"),
      c(carbapenems, VAN, TEC, TLV, DAL, ORI, DAP, LNZ, TZD, QDA, TGC, ERV, OMC, RIF),
      "any"
    )
    trans_tbl(
      3, # Sr. groups A/B/C/G
      which(x$mo %in% MO_STREP_ABCG),
      c(PEN, cephalosporins, VAN, TEC, TLV, DAL, ORI, DAP, LNZ, TZD, QDA, TGC, ERV, OMC),
      "any"
    )
    trans_tbl(
      3,
      which(x$genus == "Enterococcus"),
      c(DAP, LNZ, TGC, ERV, OMC, TEC),
      "any"
    )
    trans_tbl(
      3,
      which(x$genus == "Enterococcus" & x$species == "faecalis"),
      c(AMP, AMX),
      "any"
    )
    # Table 8
    trans_tbl(
      3,
      which(x$genus == "Bacteroides"),
      MTR,
      "any"
    )
    trans_tbl(
      3,
      which(x$genus %in% c("Clostridium", "Clostridioides") & x$species == "difficile"),
      c(MTR, VAN, FDX),
      "any"
    )
  }

  if (guideline$code == "mrgn") {
    # Germany -----------------------------------------------------------------

    # Table 1
    trans_tbl(
      2, # 3MRGN
      which((x$order == "Enterobacterales" | # following in fact the old Enterobacteriaceae classification
        (x$genus == "Acinetobacter" & x$species == "baumannii")) &
        try_ab(x[, PIP, drop = TRUE] == "R") &
        (try_ab(x[, CTX, drop = TRUE] == "R") | try_ab(x[, CAZ, drop = TRUE] == "R")) &
        (try_ab(x[, IPM, drop = TRUE] != "R") | try_ab(x[, MEM, drop = TRUE] != "R")) &
        try_ab(x[, CIP, drop = TRUE] == "R")),
      c(PIP, CTX, CAZ, IPM, MEM, CIP),
      "any"
    )

    trans_tbl(
      3, # 4MRGN, overwrites 3MRGN if applicable
      which((x$order == "Enterobacterales" | # following in fact the old Enterobacteriaceae classification
        (x$genus == "Acinetobacter" & x$species == "baumannii")) &
        try_ab(x[, PIP, drop = TRUE] == "R") &
        (try_ab(x[, CTX, drop = TRUE] == "R") | try_ab(x[, CAZ, drop = TRUE] == "R")) &
        (try_ab(x[, IPM, drop = TRUE] == "R") | try_ab(x[, MEM, drop = TRUE] == "R")) &
        try_ab(x[, CIP, drop = TRUE] == "R")),
      c(PIP, CTX, CAZ, IPM, MEM, CIP),
      "any"
    )

    trans_tbl(
      3, # 4MRGN, overwrites 3MRGN if applicable
      which((x$order == "Enterobacterales" | # following in fact the old Enterobacteriaceae classification
        (x$genus == "Acinetobacter" & x$species == "baumannii")) &
        (try_ab(x[, IPM, drop = TRUE] == "R") | try_ab(x[, MEM, drop = TRUE] == "R"))),
      c(IPM, MEM),
      "any"
    )

    trans_tbl(
      2, # 3MRGN, if only 1 group is S
      which(x$genus == "Pseudomonas" & x$species == "aeruginosa" &
        try_ab(x[, PIP, drop = TRUE] == "S") +
          try_ab(x[, CTX, drop = TRUE] == "S") +
          try_ab(x[, CAZ, drop = TRUE] == "S") +
          try_ab(x[, IPM, drop = TRUE] == "S") +
          try_ab(x[, MEM, drop = TRUE] == "S") +
          try_ab(x[, CIP, drop = TRUE] == "S") == 1),
      c(PIP, CTX, CAZ, IPM, MEM, CIP),
      "any"
    )

    trans_tbl(
      3, # 4MRGN otherwise
      which((x$genus == "Pseudomonas" & x$species == "aeruginosa") &
        try_ab(x[, PIP, drop = TRUE] == "R") &
        (try_ab(x[, CTX, drop = TRUE] == "R") | try_ab(x[, CAZ, drop = TRUE] == "R")) &
        (try_ab(x[, IPM, drop = TRUE] == "R") | try_ab(x[, MEM, drop = TRUE] == "R")) &
        try_ab(x[, CIP, drop = TRUE] == "R")),
      c(PIP, CTX, CAZ, IPM, MEM, CIP),
      "any"
    )

    x[which(x$MDRO == 2), "reason"] <- "3MRGN"
    x[which(x$MDRO == 3), "reason"] <- "4MRGN"
  }

  if (guideline$code == "brmo2024") {
    # Netherlands 2024 --------------------------------------------------------
    aminoglycosides <- c(GEN, TOB, AMK) # note 4: gentamicin or tobramycin or amikacin
    aminoglycosides_serratia_marcescens <- GEN # note 4: TOB and AMK do not count towards S. marcescens
    fluoroquinolones <- c(CIP, NOR, LVX) # note 5:  ciprofloxacin or norfloxacin or levofloxacin
    carbapenems <- carbapenems[!is.na(carbapenems)]
    carbapenems_without_imipenem <- carbapenems[carbapenems != IPM]
    amino <- AMX %or% AMP
    third <- CAZ %or% CTX
    ESBLs <- c(amino, third)
    ESBLs <- ESBLs[!is.na(ESBLs)]
    if (length(ESBLs) != 2) {
      ESBLs <- character(0)
    }

    # Enterobacterales
    if (length(ESBLs) > 0) {
      trans_tbl(
        2, # positive, unconfirmed
        rows = which(x$order == "Enterobacterales" & col_values(x, ESBLs[1]) == "R" & col_values(x, ESBLs[2]) == "R" & is.na(esbl)),
        cols = c(AMX %or% AMP, cephalosporins_3rd),
        any_all = "all",
        reason = "Enterobacterales: potential ESBL"
      )
    }
    trans_tbl(
      3, # positive
      rows = which(x$order == "Enterobacterales" & esbl == TRUE),
      cols = "any",
      any_all = "any",
      reason = "Enterobacterales: ESBL"
    )
    trans_tbl(
      3,
      rows = which(x$order == "Enterobacterales" & (x$genus %in% c("Proteus", "Providencia") | paste(x$genus, x$species) %in% c("Serratia marcescens", "Morganella morganii"))),
      cols = carbapenems_without_imipenem,
      any_all = "any",
      reason = "Enterobacterales: carbapenem resistance"
    )
    trans_tbl(
      3,
      rows = which(x$order == "Enterobacterales" & !(x$genus %in% c("Proteus", "Providencia") | paste(x$genus, x$species) %in% c("Serratia marcescens", "Morganella morganii"))),
      cols = carbapenems,
      any_all = "any",
      reason = "Enterobacterales: carbapenem resistance"
    )
    trans_tbl(
      3,
      rows = which(x$order == "Enterobacterales" & carbapenemase == TRUE),
      cols = "any",
      any_all = "any",
      reason = "Enterobacterales: carbapenemase"
    )
    trans_tbl(
      3,
      rows = which(col_values(x, SXT) == "R" &
        (col_values(x, GEN) == "R" | col_values(x, TOB) == "R" | col_values(x, AMK) == "R") &
        (col_values(x, CIP) == "R" | col_values(x, NOR) == "R" | col_values(x, LVX) == "R") &
        (x$genus %in% c("Enterobacter", "Providencia") | paste(x$genus, x$species) %in% c("Citrobacter freundii", "Klebsiella aerogenes", "Hafnia alvei", "Morganella morganii"))),
      cols = c(SXT, aminoglycosides, fluoroquinolones),
      any_all = "any",
      reason = "Enterobacterales group II: aminoglycoside + fluoroquinolone + cotrimoxazol"
    )
    trans_tbl(
      3,
      rows = which(col_values(x, SXT) == "R" &
        col_values(x, GEN) == "R" &
        (col_values(x, CIP) == "R" | col_values(x, NOR) == "R" | col_values(x, LVX) == "R") &
        paste(x$genus, x$species) == "Serratia marcescens"),
      cols = c(SXT, aminoglycosides_serratia_marcescens, fluoroquinolones),
      any_all = "any",
      reason = "Enterobacterales group II: aminoglycoside + fluoroquinolone + cotrimoxazol"
    )

    # Acinetobacter baumannii-calcoaceticus complex
    trans_tbl(
      3,
      rows = which((col_values(x, GEN) == "R" | col_values(x, TOB) == "R" | col_values(x, AMK) == "R") &
        (col_values(x, CIP) == "R" | col_values(x, LVX) == "R") &
        x[[col_mo]] %in% AMR::microorganisms.groups$mo[AMR::microorganisms.groups$mo_group_name == "Acinetobacter baumannii complex"]),
      cols = c(aminoglycosides, CIP, LVX),
      any_all = "any",
      reason = "A. baumannii-calcoaceticus complex: aminoglycoside + ciprofloxacin or levofloxacin"
    )
    trans_tbl(
      2, # unconfirmed
      rows = which(x[[col_mo]] %in% AMR::microorganisms.groups$mo[AMR::microorganisms.groups$mo_group_name == "Acinetobacter baumannii complex"] & is.na(carbapenemase)),
      cols = carbapenems,
      any_all = "any",
      reason = "A. baumannii-calcoaceticus complex: potential carbapenemase"
    )
    trans_tbl(
      3,
      rows = which(x[[col_mo]] %in% AMR::microorganisms.groups$mo[AMR::microorganisms.groups$mo_group_name == "Acinetobacter baumannii complex"] & carbapenemase == TRUE),
      cols = carbapenems,
      any_all = "any",
      reason = "A. baumannii-calcoaceticus complex: carbapenemase"
    )

    # Pseudomonas aeruginosa
    x$psae <- 0
    x$psae <- x$psae + ifelse(NA_as_FALSE(col_values(x, TOB) == "R") | NA_as_FALSE(col_values(x, AMK) == "R"), 1, 0)
    x$psae <- x$psae + ifelse(NA_as_FALSE(col_values(x, IPM) == "R") | NA_as_FALSE(col_values(x, MEM) == "R"), 1, 0)
    x$psae <- x$psae + ifelse(NA_as_FALSE(col_values(x, PIP) == "R") | NA_as_FALSE(col_values(x, TZP) == "R"), 1, 0)
    x$psae <- x$psae + ifelse(NA_as_FALSE(col_values(x, CAZ) == "R") | NA_as_FALSE(col_values(x, CZA) == "R"), 1, 0)
    x$psae <- x$psae + ifelse(NA_as_FALSE(col_values(x, CIP) == "R") | NA_as_FALSE(col_values(x, NOR) == "R") | NA_as_FALSE(col_values(x, LVX) == "R"), 1, 0)
    trans_tbl(
      1,
      rows = which(x$genus == "Pseudomonas" & x$species == "aeruginosa"),
      cols = "any",
      any_all = "all", # this will set all negatives to "guideline criteria not met" instead of "not covered by guideline"
      reason = "guideline criteria not met"
    )
    trans_tbl(
      3,
      rows = which(x$genus == "Pseudomonas" & x$species == "aeruginosa" & x$psae >= 3),
      cols = "any",
      any_all = "any", # this is the actual one, overwriting the ones with x$psae >= 3
      reason = "P. aeruginosa: at least 3 classes contain R"
    )

    # Enterococcus faecium
    trans_tbl(
      3,
      rows = which(x$genus == "Enterococcus" & x$species == "faecium"),
      cols = c(PEN %or% AMX %or% AMP, VAN),
      any_all = "all",
      reason = "E. faecium: vancomycin + penicillin group"
    )
    trans_tbl(
      3,
      rows = which(x$genus == "Enterococcus" & x$species == "faecium" & (vanA == TRUE | vanB == TRUE)),
      cols = c(PEN, AMX, AMP, VAN),
      any_all = "any",
      reason = "E. faecium: vanA/vanB gene + penicillin group"
    )

    # Staphylococcus aureus complex (= aureus, argenteus or schweitzeri)
    trans_tbl(
      2,
      rows = which(x$genus == "Staphylococcus" & x$species %in% c("aureus", "argenteus", "schweitzeri") & (is.na(mecA) | is.na(mecC))),
      cols = c(AMC, TZP, FLC, OXA, FOX, FOX1),
      any_all = "any",
      reason = "S. aureus complex: potential MRSA"
    )
    trans_tbl(
      3,
      rows = which(x$genus == "Staphylococcus" & x$species %in% c("aureus", "argenteus", "schweitzeri") & (mecA == TRUE | mecC == TRUE)),
      cols = "any",
      any_all = "any",
      reason = "S. aureus complex: mecA/mecC gene"
    )

    # Candida auris
    trans_tbl(
      3,
      rows = which(x$genus == "Candida" & x$species == "auris"),
      cols = "any",
      any_all = "any",
      reason = "C. auris: regardless of resistance"
    )
  }

  if (guideline$code == "brmo2017") {
    # Netherlands 2017 --------------------------------------------------------
    aminoglycosides <- aminoglycosides[!is.na(aminoglycosides)]
    fluoroquinolones <- fluoroquinolones[!is.na(fluoroquinolones)]
    carbapenems <- carbapenems[!is.na(carbapenems)]
    amino <- AMX %or% AMP
    third <- CAZ %or% CTX
    ESBLs <- c(amino, third)
    ESBLs <- ESBLs[!is.na(ESBLs)]
    if (length(ESBLs) != 2) {
      ESBLs <- character(0)
    }

    # Table 1
    trans_tbl(
      3,
      which(x$order == "Enterobacterales"), # following in fact the old Enterobacteriaceae classification
      c(aminoglycosides, fluoroquinolones),
      "all"
    )

    trans_tbl(
      2,
      which(x$order == "Enterobacterales"), # following in fact the old Enterobacteriaceae classification
      carbapenems,
      "any"
    )

    trans_tbl(
      2,
      which(x$order == "Enterobacterales"), # following in fact the old Enterobacteriaceae classification
      ESBLs,
      "all"
    )

    # Table 2
    trans_tbl(
      2,
      which(x$genus == "Acinetobacter"),
      c(carbapenems),
      "any"
    )
    trans_tbl(
      3,
      which(x$genus == "Acinetobacter"),
      c(aminoglycosides, fluoroquinolones),
      "all"
    )

    trans_tbl(
      3,
      which(x$genus == "Stenotrophomonas" & x$species == "maltophilia"),
      SXT,
      "all"
    )

    if (!ab_missing(MEM) && !ab_missing(IPM) &&
      !ab_missing(GEN) && !ab_missing(TOB) &&
      !ab_missing(CIP) &&
      !ab_missing(CAZ) &&
      !ab_missing(TZP)) {
      x$psae <- 0
      x[which(x[, MEM, drop = TRUE] == "R" | x[, IPM, drop = TRUE] == "R"), "psae"] <- 1 + x[which(x[, MEM, drop = TRUE] == "R" | x[, IPM, drop = TRUE] == "R"), "psae"]
      x[which(x[, GEN, drop = TRUE] == "R" & x[, TOB, drop = TRUE] == "R"), "psae"] <- 1 + x[which(x[, GEN, drop = TRUE] == "R" & x[, TOB, drop = TRUE] == "R"), "psae"]
      x[which(x[, CIP, drop = TRUE] == "R"), "psae"] <- 1 + x[which(x[, CIP, drop = TRUE] == "R"), "psae"]
      x[which(x[, CAZ, drop = TRUE] == "R"), "psae"] <- 1 + x[which(x[, CAZ, drop = TRUE] == "R"), "psae"]
      x[which(x[, TZP, drop = TRUE] == "R"), "psae"] <- 1 + x[which(x[, TZP, drop = TRUE] == "R"), "psae"]
    } else {
      x$psae <- 0
    }
    trans_tbl(
      3,
      which(x$genus == "Pseudomonas" & x$species == "aeruginosa" & x$psae >= 3),
      c(CAZ, CIP, GEN, IPM, MEM, TOB, TZP),
      "any"
    )
    x[which(x$genus == "Pseudomonas" & x$species == "aeruginosa" & x$psae >= 3), "reason"] <- paste0("at least 3 classes contain R", ifelse(!isTRUE(combine_SI), " or I", ""))

    # Table 3
    trans_tbl(
      3,
      which(x$genus == "Streptococcus" & x$species == "pneumoniae"),
      PEN,
      "all"
    )
    trans_tbl(
      3,
      which(x$genus == "Streptococcus" & x$species == "pneumoniae"),
      VAN,
      "all"
    )
    trans_tbl(
      3,
      which(x$genus == "Enterococcus" & x$species == "faecium"),
      c(PEN, VAN),
      "all"
    )
  }

  if (guideline$code == "tb") {
    # Tuberculosis ------------------------------------------------------------
    prepare_drug <- function(ab) {
      # returns vector values of drug
      # if `ab` is a column name, looks up the values in `x`
      if (length(ab) == 1 && is.character(ab)) {
        if (ab %in% colnames(x)) {
          ab <- x[, ab, drop = TRUE]
        }
      }
      ab <- as.character(as.sir(ab))
      ab[is.na(ab)] <- ""
      ab
    }
    drug_is_R <- function(ab) {
      # returns [logical] vector
      ab <- prepare_drug(ab)
      if (length(ab) == 0) {
        rep(FALSE, NROW(x))
      } else if (length(ab) == 1) {
        rep(ab, NROW(x)) == "R"
      } else {
        ab == "R"
      }
    }
    drug_is_not_R <- function(ab) {
      # returns [logical] vector
      ab <- prepare_drug(ab)
      if (length(ab) == 0) {
        rep(TRUE, NROW(x))
      } else if (length(ab) == 1) {
        rep(ab, NROW(x)) != "R"
      } else {
        ab != "R"
      }
    }

    x$mono_count <- 0
    x[drug_is_R(INH), "mono_count"] <- x[drug_is_R(INH), "mono_count", drop = TRUE] + 1
    x[drug_is_R(RIF), "mono_count"] <- x[drug_is_R(RIF), "mono_count", drop = TRUE] + 1
    x[drug_is_R(ETH), "mono_count"] <- x[drug_is_R(ETH), "mono_count", drop = TRUE] + 1
    x[drug_is_R(PZA), "mono_count"] <- x[drug_is_R(PZA), "mono_count", drop = TRUE] + 1
    x[drug_is_R(RIB), "mono_count"] <- x[drug_is_R(RIB), "mono_count", drop = TRUE] + 1
    x[drug_is_R(RFP), "mono_count"] <- x[drug_is_R(RFP), "mono_count", drop = TRUE] + 1

    x$mono <- x$mono_count > 0
    x$poly <- x$mono_count > 1 & drug_is_not_R(RIF) & drug_is_not_R(INH)
    x$mdr <- drug_is_R(RIF) & drug_is_R(INH)
    x$xdr <- drug_is_R(LVX) | drug_is_R(MFX) | drug_is_R(GAT)
    x$second <- drug_is_R(CAP) | drug_is_R(KAN) | drug_is_R(AMK)
    x$xdr <- x$mdr & x$xdr & x$second
    x$MDRO <- ifelse(x$xdr, 5,
      ifelse(x$mdr, 4,
        ifelse(x$poly, 3,
          ifelse(x$mono, 2,
            1
          )
        )
      )
    )
    # keep all real TB, make other species NA
    x$MDRO <- ifelse(x$fullname == "Mycobacterium tuberculosis", x$MDRO, NA_real_)
    x$reason <- "PDR/MDR/XDR criteria were met"
  }

  # some more info on negative results
  if (isTRUE(verbose)) {
    if (guideline$code == "cmi2012") {
      x[which(x$MDRO == 1 & !is.na(x$classes_affected)), "reason"] <- paste0(
        x$classes_affected[which(x$MDRO == 1 & !is.na(x$classes_affected))],
        " of ",
        x$classes_available[which(x$MDRO == 1 & !is.na(x$classes_affected))],
        " available classes contain R",
        ifelse(!isTRUE(combine_SI), " or I", ""),
        " (3 required for MDR)"
      )
    } else {
      # x[which(x$MDRO == 1), "reason"] <- "too few antimicrobials are R"
    }
  }

  if (isTRUE(info.bak)) {
    cat(group_msg)
    if (sum(!is.na(x$MDRO)) == 0) {
      cat(font_bold(paste0("=> Found 0 MDROs since no isolates are covered by the guideline")))
    } else {
      cat(font_bold(paste0(
        "=> Found ", sum(x$MDRO %in% 2:5, na.rm = TRUE), " ", guideline$type, " out of ", sum(!is.na(x$MDRO)),
        " isolates (", trimws(percentage(sum(x$MDRO %in% 2:5, na.rm = TRUE) / sum(!is.na(x$MDRO)))), ")"
      )))
    }
  }

  # Fill in blanks ----
  # for rows that have no results
  x_transposed <- as.list(as.data.frame(t(x[, cols_ab, drop = FALSE]),
    stringsAsFactors = FALSE
  ))
  rows_empty <- which(vapply(
    FUN.VALUE = logical(1),
    x_transposed,
    function(y) all(is.na(y))
  ))
  if (length(rows_empty) > 0) {
    if (isTRUE(info.bak)) {
      cat(font_italic(paste0(" (", length(rows_empty), " isolates had no test results)\n")))
    }
  } else if (isTRUE(info.bak)) {
    cat("\n")
  }

  if (isTRUE(info.bak) && !isTRUE(verbose)) {
    cat("\nRerun with 'verbose = TRUE' to retrieve detailed info and reasons for every MDRO classification.\n")
  }

  # Results ----
  if (guideline$code == "cmi2012") {
    if (any(x$MDRO == -1, na.rm = TRUE)) {
      if (message_not_thrown_before("mdro", "availability")) {
        warning_(
          "in `mdro()`: NA introduced for isolates where the available percentage of antimicrobial classes was below ",
          percentage(pct_required_classes), " (set with `pct_required_classes`)"
        )
      }
      # set these -1s to NA
      x[which(x$MDRO == -1), "MDRO"] <- NA_integer_
    }
    x$MDRO <- factor(
      x = x$MDRO,
      levels = 1:4,
      labels = c(
        "Negative", "Multi-drug-resistant (MDR)",
        "Extensively drug-resistant (XDR)", "Pandrug-resistant (PDR)"
      ),
      ordered = TRUE
    )
  } else if (guideline$code == "tb") {
    x$MDRO <- factor(
      x = x$MDRO,
      levels = 1:5,
      labels = c(
        "Negative", "Mono-resistant", "Poly-resistant",
        "Multi-drug-resistant", "Extensively drug-resistant"
      ),
      ordered = TRUE
    )
  } else if (guideline$code == "mrgn") {
    x$MDRO <- factor(
      x = x$MDRO,
      levels = 1:3,
      labels = c("Negative", "3MRGN", "4MRGN"),
      ordered = TRUE
    )
  } else {
    x$MDRO <- factor(
      x = x$MDRO,
      levels = 1:3,
      labels = c("Negative", "Positive, unconfirmed", "Positive"),
      ordered = TRUE
    )
  }

  if (isTRUE(verbose)) {
    # fill in empty reasons
    x$reason[is.na(x$reason)] <- "not covered by guideline"
    x[rows_empty, "reason"] <- paste(x[rows_empty, "reason"], "(note: no available test results)")
    # format data set
    colnames(x)[colnames(x) == col_mo] <- "microorganism"
    x$microorganism <- mo_name(x$microorganism, language = NULL)
    x$guideline <- paste0(guideline$author, " - ", guideline$name, ", ", guideline$version, ")")
    x[, c(
      "row_number",
      "microorganism",
      "MDRO",
      "reason",
      "all_nonsusceptible_columns",
      "guideline"
    ),
    drop = FALSE
    ]
  } else {
    x$MDRO
  }
}

#' @rdname mdro
#' @export
brmo <- function(x = NULL, only_sir_columns = any(is.sir(x)), ...) {
  meet_criteria(x, allow_class = "data.frame", allow_NULL = TRUE)
  meet_criteria(only_sir_columns, allow_class = "logical", has_length = 1)
  stop_if(
    "guideline" %in% names(list(...)),
    "argument `guideline` must not be set since this is a guideline-specific function"
  )
  mdro(x = x, only_sir_columns = only_sir_columns, guideline = "BRMO", ...)
}


#' @rdname mdro
#' @export
mrgn <- function(x = NULL, only_sir_columns = any(is.sir(x)), verbose = FALSE, ...) {
  meet_criteria(x, allow_class = "data.frame", allow_NULL = TRUE)
  meet_criteria(only_sir_columns, allow_class = "logical", has_length = 1)
  stop_if(
    "guideline" %in% names(list(...)),
    "argument `guideline` must not be set since this is a guideline-specific function"
  )
  mdro(x = x, only_sir_columns = only_sir_columns, verbose = verbose, guideline = "MRGN", ...)
}

#' @rdname mdro
#' @export
mdr_tb <- function(x = NULL, only_sir_columns = any(is.sir(x)), verbose = FALSE, ...) {
  meet_criteria(x, allow_class = "data.frame", allow_NULL = TRUE)
  meet_criteria(only_sir_columns, allow_class = "logical", has_length = 1)
  stop_if(
    "guideline" %in% names(list(...)),
    "argument `guideline` must not be set since this is a guideline-specific function"
  )
  mdro(x = x, only_sir_columns = only_sir_columns, verbose = verbose, guideline = "TB", ...)
}

#' @rdname mdro
#' @export
mdr_cmi2012 <- function(x = NULL, only_sir_columns = any(is.sir(x)), verbose = FALSE, ...) {
  meet_criteria(x, allow_class = "data.frame", allow_NULL = TRUE)
  meet_criteria(only_sir_columns, allow_class = "logical", has_length = 1)
  stop_if(
    "guideline" %in% names(list(...)),
    "argument `guideline` must not be set since this is a guideline-specific function"
  )
  mdro(x = x, only_sir_columns = only_sir_columns, verbose = verbose, guideline = "CMI 2012", ...)
}

#' @rdname mdro
#' @export
eucast_exceptional_phenotypes <- function(x = NULL, only_sir_columns = any(is.sir(x)), verbose = FALSE, ...) {
  meet_criteria(x, allow_class = "data.frame", allow_NULL = TRUE)
  meet_criteria(only_sir_columns, allow_class = "logical", has_length = 1)
  stop_if(
    "guideline" %in% names(list(...)),
    "argument `guideline` must not be set since this is a guideline-specific function"
  )
  mdro(x = x, only_sir_columns = only_sir_columns, verbose = verbose, guideline = "EUCAST", ...)
}
