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

#' Interpret MIC and Disk Values, or Clean Raw SIR Data
#'
#' Interpret minimum inhibitory concentration (MIC) values and disk diffusion diameters according to EUCAST or CLSI, or clean up existing SIR values. This transforms the input to a new class [`sir`], which is an ordered [factor] with levels `S < I < R`.
#' @rdname as.sir
#' @param x vector of values (for class [`mic`]: MIC values in mg/L, for class [`disk`]: a disk diffusion radius in millimetres)
#' @param mo any (vector of) text that can be coerced to valid microorganism codes with [as.mo()], can be left empty to determine it automatically
#' @param ab any (vector of) text that can be coerced to a valid antimicrobial drug code with [as.ab()]
#' @param uti (Urinary Tract Infection) A vector with [logical]s (`TRUE` or `FALSE`) to specify whether a UTI specific interpretation from the guideline should be chosen. For using [as.sir()] on a [data.frame], this can also be a column containing [logical]s or when left blank, the data set will be searched for a column 'specimen', and rows within this column containing 'urin' (such as 'urine', 'urina') will be regarded isolates from a UTI. See *Examples*.
#' @inheritParams first_isolate
#' @param guideline defaults to EUCAST `r  max(as.integer(gsub("[^0-9]", "", subset(clinical_breakpoints, guideline %like% "EUCAST")$guideline)))` (the latest implemented EUCAST guideline in the [clinical_breakpoints] data set), but can be set with the [option][options()] `AMR_guideline`. Supports EUCAST (`r min(as.integer(gsub("[^0-9]", "", subset(clinical_breakpoints, guideline %like% "EUCAST")$guideline)))`-`r max(as.integer(gsub("[^0-9]", "", subset(clinical_breakpoints, guideline %like% "EUCAST")$guideline)))`) and CLSI (`r min(as.integer(gsub("[^0-9]", "", subset(clinical_breakpoints, guideline %like% "CLSI")$guideline)))`-`r max(as.integer(gsub("[^0-9]", "", subset(clinical_breakpoints, guideline %like% "CLSI")$guideline)))`), see *Details*.
#' @param conserve_capped_values a [logical] to indicate that MIC values starting with `">"` (but not `">="`) must always return "R" , and that MIC values starting with `"<"` (but not `"<="`) must always return "S"
#' @param add_intrinsic_resistance *(only useful when using a EUCAST guideline)* a [logical] to indicate whether intrinsic antibiotic resistance must also be considered for applicable bug-drug combinations, meaning that e.g. ampicillin will always return "R" in *Klebsiella* species. Determination is based on the [intrinsic_resistant] data set, that itself is based on `r format_eucast_version_nr(3.3)`.
#' @param reference_data a [data.frame] to be used for interpretation, which defaults to the [clinical_breakpoints] data set. Changing this argument allows for using own interpretation guidelines. This argument must contain a data set that is equal in structure to the [clinical_breakpoints] data set (same column names and column types). Please note that the `guideline` argument will be ignored when `reference_data` is manually set.
#' @param threshold maximum fraction of invalid antimicrobial interpretations of `x`, see *Examples*
#' @param ... for using on a [data.frame]: names of columns to apply [as.sir()] on (supports tidy selection such as `column1:column4`). Otherwise: arguments passed on to methods.
#' @details
#' ### How it Works
#'
#' The [as.sir()] function works in four ways:
#'
#' 1. For **cleaning raw / untransformed data**. The data will be cleaned to only contain values S, I and R and will try its best to determine this with some intelligence. For example, mixed values with SIR interpretations and MIC values such as `"<0.25; S"` will be coerced to `"S"`. Combined interpretations for multiple test methods (as seen in laboratory records) such as `"S; S"` will be coerced to `"S"`, but a value like `"S; I"` will return `NA` with a warning that the input is unclear.
#'
#' 2. For **interpreting minimum inhibitory concentration (MIC) values** according to EUCAST or CLSI. You must clean your MIC values first using [as.mic()], that also gives your columns the new data class [`mic`]. Also, be sure to have a column with microorganism names or codes. It will be found automatically, but can be set manually using the `mo` argument.
#'    * Using `dplyr`, SIR interpretation can be done very easily with either:
#'      ```
#'      your_data %>% mutate_if(is.mic, as.sir)
#'      your_data %>% mutate(across(where(is.mic), as.sir))
#'      ```
#'    * Operators like "<=" will be stripped before interpretation. When using `conserve_capped_values = TRUE`, an MIC value of e.g. ">2" will always return "R", even if the breakpoint according to the chosen guideline is ">=4". This is to prevent that capped values from raw laboratory data would not be treated conservatively. The default behaviour (`conserve_capped_values = FALSE`) considers ">2" to be lower than ">=4" and might in this case return "S" or "I".
#' 3. For **interpreting disk diffusion diameters** according to EUCAST or CLSI. You must clean your disk zones first using [as.disk()], that also gives your columns the new data class [`disk`]. Also, be sure to have a column with microorganism names or codes. It will be found automatically, but can be set manually using the `mo` argument.
#'    * Using `dplyr`, SIR interpretation can be done very easily with either:
#'      ```
#'      your_data %>% mutate_if(is.disk, as.sir)
#'      your_data %>% mutate(across(where(is.disk), as.sir))
#'      ```
#' 4. For **interpreting a complete data set**, with automatic determination of MIC values, disk diffusion diameters, microorganism names or codes, and antimicrobial test results. This is done very simply by running `as.sir(your_data)`.
#'
#' For points 2, 3 and 4: Use [sir_interpretation_history()] to retrieve a [data.frame] (or [tibble][tibble::tibble()] if the `tibble` package is installed) with all results of the last [as.sir()] call.
#'
#' ### Supported Guidelines
#'
#' For interpreting MIC values as well as disk diffusion diameters, currently implemented guidelines are EUCAST (`r min(as.integer(gsub("[^0-9]", "", subset(clinical_breakpoints, guideline %like% "EUCAST")$guideline)))`-`r max(as.integer(gsub("[^0-9]", "", subset(clinical_breakpoints, guideline %like% "EUCAST")$guideline)))`) and CLSI (`r min(as.integer(gsub("[^0-9]", "", subset(clinical_breakpoints, guideline %like% "CLSI")$guideline)))`-`r max(as.integer(gsub("[^0-9]", "", subset(clinical_breakpoints, guideline %like% "CLSI")$guideline)))`).
#'
#' Thus, the `guideline` argument must be set to e.g., ``r paste0('"', subset(clinical_breakpoints, guideline %like% "EUCAST")$guideline[1], '"')`` or ``r paste0('"', subset(clinical_breakpoints, guideline %like% "CLSI")$guideline[1], '"')``. By simply using `"EUCAST"` (the default) or `"CLSI"` as input, the latest included version of that guideline will automatically be selected. You can set your own data set using the `reference_data` argument. The `guideline` argument will then be ignored.
#' 
#' You can set the default guideline with the `AMR_guideline` [option][options()] (e.g. in your `.Rprofile` file), such as:
#' 
#' ```
#'   options(AMR_guideline = "CLSI")
#'   options(AMR_guideline = "CLSI 2018")
#'   options(AMR_guideline = "EUCAST 2020")
#'   # or to reset:
#'   options(AMR_guideline = NULL)
#' ```
#'
#' ### After Interpretation
#'
#' After using [as.sir()], you can use the [eucast_rules()] defined by EUCAST to (1) apply inferred susceptibility and resistance based on results of other antimicrobials and (2) apply intrinsic resistance based on taxonomic properties of a microorganism.
#'
#' ### Machine-Readable Interpretation Guidelines
#'
#' The repository of this package [contains a machine-readable version](https://github.com/msberends/AMR/blob/main/data-raw/clinical_breakpoints.txt) of all guidelines. This is a CSV file consisting of `r format(nrow(AMR::clinical_breakpoints), big.mark = ",")` rows and `r ncol(AMR::clinical_breakpoints)` columns. This file is machine-readable, since it contains one row for every unique combination of the test method (MIC or disk diffusion), the antimicrobial drug and the microorganism. **This allows for easy implementation of these rules in laboratory information systems (LIS)**. Note that it only contains interpretation guidelines for humans - interpretation guidelines from CLSI for animals were removed.
#'
#' ### Other
#'
#' The function [is.sir()] detects if the input contains class `sir`. If the input is a [data.frame], it iterates over all columns and returns a [logical] vector.
#'
#' The function [is_sir_eligible()] returns `TRUE` when a columns contains at most 5% invalid antimicrobial interpretations (not S and/or I and/or R), and `FALSE` otherwise. The threshold of 5% can be set with the `threshold` argument. If the input is a [data.frame], it iterates over all columns and returns a [logical] vector.
#' @section Interpretation of R and S/I:
#' In 2019, the European Committee on Antimicrobial Susceptibility Testing (EUCAST) has decided to change the definitions of susceptibility testing categories R and S/I as shown below (<https://www.eucast.org/newsiandr/>).
#'
#' - **R = Resistant**\cr
#'   A microorganism is categorised as *Resistant* when there is a high likelihood of therapeutic failure even when there is increased exposure. Exposure is a function of how the mode of administration, dose, dosing interval, infusion time, as well as distribution and excretion of the antimicrobial agent will influence the infecting organism at the site of infection.
#' - **S = Susceptible**\cr
#'   A microorganism is categorised as *Susceptible, standard dosing regimen*, when there is a high likelihood of therapeutic success using a standard dosing regimen of the agent.
#' - **I = Susceptible, Increased exposure**\cr
#'   A microorganism is categorised as *Susceptible, Increased exposure* when there is a high likelihood of therapeutic success because exposure to the agent is increased by adjusting the dosing regimen or by its concentration at the site of infection.
#'
#' This AMR package honours this insight. Use [susceptibility()] (equal to [proportion_SI()]) to determine antimicrobial susceptibility and [count_susceptible()] (equal to [count_SI()]) to count susceptible isolates.
#' @return Ordered [factor] with new class `sir`
#' @aliases sir
#' @export
#' @seealso [as.mic()], [as.disk()], [as.mo()]
#' @source
#' For interpretations of minimum inhibitory concentration (MIC) values and disk diffusion diameters:
#'
#' - **M39 Analysis and Presentation of Cumulative Antimicrobial Susceptibility Test Data**, `r min(as.integer(gsub("[^0-9]", "", subset(clinical_breakpoints, guideline %like% "CLSI")$guideline)))`-`r max(as.integer(gsub("[^0-9]", "", subset(clinical_breakpoints, guideline %like% "CLSI")$guideline)))`, *Clinical and Laboratory Standards Institute* (CLSI). <https://clsi.org/standards/products/microbiology/documents/m39/>.
#' - **M100 Performance Standard for Antimicrobial Susceptibility Testing**, `r min(as.integer(gsub("[^0-9]", "", subset(clinical_breakpoints, guideline %like% "CLSI")$guideline)))`-`r max(as.integer(gsub("[^0-9]", "", subset(clinical_breakpoints, guideline %like% "CLSI")$guideline)))`, *Clinical and Laboratory Standards Institute* (CLSI). <https://clsi.org/standards/products/microbiology/documents/m100/>.
#' - **Breakpoint tables for interpretation of MICs and zone diameters**, `r min(as.integer(gsub("[^0-9]", "", subset(clinical_breakpoints, guideline %like% "EUCAST")$guideline)))`-`r max(as.integer(gsub("[^0-9]", "", subset(clinical_breakpoints, guideline %like% "EUCAST")$guideline)))`, *European Committee on Antimicrobial Susceptibility Testing* (EUCAST). <https://www.eucast.org/clinical_breakpoints>.
#' @inheritSection AMR Reference Data Publicly Available
#' @examples
#' example_isolates
#' summary(example_isolates) # see all SIR results at a glance
#'
#' # For INTERPRETING disk diffusion and MIC values -----------------------
#'
#' # a whole data set, even with combined MIC values and disk zones
#' df <- data.frame(
#'   microorganism = "Escherichia coli",
#'   AMP = as.mic(8),
#'   CIP = as.mic(0.256),
#'   GEN = as.disk(18),
#'   TOB = as.disk(16),
#'   ERY = "R"
#' )
#' as.sir(df)
#'
#' # return a 'logbook' about the results:
#' sir_interpretation_history()
#'
#' # for single values
#' as.sir(
#'   x = as.mic(2),
#'   mo = as.mo("S. pneumoniae"),
#'   ab = "AMP",
#'   guideline = "EUCAST"
#' )
#'
#' as.sir(
#'   x = as.disk(18),
#'   mo = "Strep pneu", # `mo` will be coerced with as.mo()
#'   ab = "ampicillin", # and `ab` with as.ab()
#'   guideline = "EUCAST"
#' )
#'
#' \donttest{
#' # the dplyr way
#' if (require("dplyr")) {
#'   df %>% mutate_if(is.mic, as.sir)
#'   df %>% mutate_if(function(x) is.mic(x) | is.disk(x), as.sir)
#'   df %>% mutate(across(where(is.mic), as.sir))
#'   df %>% mutate_at(vars(AMP:TOB), as.sir)
#'   df %>% mutate(across(AMP:TOB, as.sir))
#'
#'   df %>%
#'     mutate_at(vars(AMP:TOB), as.sir, mo = .$microorganism)
#'
#'   # to include information about urinary tract infections (UTI)
#'   data.frame(
#'     mo = "E. coli",
#'     NIT = c("<= 2", 32),
#'     from_the_bladder = c(TRUE, FALSE)
#'   ) %>%
#'     as.sir(uti = "from_the_bladder")
#'
#'   data.frame(
#'     mo = "E. coli",
#'     NIT = c("<= 2", 32),
#'     specimen = c("urine", "blood")
#'   ) %>%
#'     as.sir() # automatically determines urine isolates
#'
#'   df %>%
#'     mutate_at(vars(AMP:TOB), as.sir, mo = "E. coli", uti = TRUE)
#' }
#'
#' # For CLEANING existing SIR values ------------------------------------
#'
#' as.sir(c("S", "I", "R", "A", "B", "C"))
#' as.sir("<= 0.002; S") # will return "S"
#' sir_data <- as.sir(c(rep("S", 474), rep("I", 36), rep("R", 370)))
#' is.sir(sir_data)
#' plot(sir_data) # for percentages
#' barplot(sir_data) # for frequencies
#'
#' # the dplyr way
#' if (require("dplyr")) {
#'   example_isolates %>%
#'     mutate_at(vars(PEN:RIF), as.sir)
#'   # same:
#'   example_isolates %>%
#'     as.sir(PEN:RIF)
#'
#'   # fastest way to transform all columns with already valid AMR results to class `sir`:
#'   example_isolates %>%
#'     mutate_if(is_sir_eligible, as.sir)
#'
#'   # since dplyr 1.0.0, this can also be:
#'   # example_isolates %>%
#'   #   mutate(across(where(is_sir_eligible), as.sir))
#' }
#' }
as.sir <- function(x, ...) {
  UseMethod("as.sir")
}

#' @rdname as.sir
#' @details `NA_sir_` is a missing value of the new `sir` class, analogous to e.g. base \R's [`NA_character_`][base::NA].
#' @export
NA_sir_ <- set_clean_class(factor(NA, levels = c("S", "I", "R"), ordered = TRUE),
  new_class = c("sir", "ordered", "factor")
)

#' @rdname as.sir
#' @export
is.sir <- function(x) {
  if (inherits(x, "data.frame")) {
    unname(vapply(FUN.VALUE = logical(1), x, is.sir))
  } else {
    inherits(x, "sir")
  }
}

#' @rdname as.sir
#' @export
is_sir_eligible <- function(x, threshold = 0.05) {
  meet_criteria(threshold, allow_class = "numeric", has_length = 1)

  if (inherits(x, "data.frame")) {
    # iterate this function over all columns
    return(unname(vapply(FUN.VALUE = logical(1), x, is_sir_eligible)))
  }

  stop_if(NCOL(x) > 1, "`x` must be a one-dimensional vector.")
  if (any(c(
    "numeric",
    "integer",
    "mo",
    "ab",
    "Date",
    "POSIXt",
    "raw",
    "hms",
    "mic",
    "disk"
  )
  %in% class(x))) {
    # no transformation needed
    return(FALSE)
  } else if (all(x %in% c("R", "S", "I", NA)) & !all(is.na(x))) {
    return(TRUE)
  } else if (!any(c("R", "S", "I") %in% x, na.rm = TRUE) & !all(is.na(x))) {
    return(FALSE)
  } else {
    x <- x[!is.na(x) & !is.null(x) & !x %in% c("", "-", "NULL")]
    if (length(x) == 0) {
      # no other values than empty
      cur_col <- get_current_column()
      if (!is.null(cur_col)) {
        ab <- suppressWarnings(as.ab(cur_col, fast_mode = TRUE, info = FALSE))
        if (!is.na(ab)) {
          # this is a valid antibiotic drug code
          message_(
            "Column '", font_bold(cur_col), "' is SIR eligible (despite only having empty values), since it seems to be ",
            ab_name(ab, language = NULL, tolower = TRUE), " (", ab, ")"
          )
          return(TRUE)
        }
      }
      # all values empty and no antibiotic col name - return FALSE
      return(FALSE)
    }
    # transform all values and see if it meets the set threshold
    checked <- suppressWarnings(as.sir(x))
    outcome <- sum(is.na(checked)) / length(x)
    outcome <= threshold
  }
}

#' @export
# extra param: warn (logical, to never throw a warning)
as.sir.default <- function(x, ...) {
  if (is.sir(x)) {
    return(x)
  }

  x.bak <- x
  x <- as.character(x) # this is needed to prevent the vctrs pkg from throwing an error
  
  if (inherits(x.bak, c("integer", "numeric", "double")) && all(x %in% c(1:3, NA))) {
    # support haven package for importing e.g., from SPSS - it adds the 'labels' attribute
    lbls <- attributes(x.bak)$labels
    if (!is.null(lbls) && all(c("R", "S", "I") %in% names(lbls)) && all(c(1:3) %in% lbls)) {
      x[x.bak == 1] <- names(lbls[lbls == 1])
      x[x.bak == 2] <- names(lbls[lbls == 2])
      x[x.bak == 3] <- names(lbls[lbls == 3])
    } else {
      x[x.bak == 1] <- "S"
      x[x.bak == 2] <- "I"
      x[x.bak == 3] <- "R"
    }
  } else if (inherits(x.bak, "character") && all(x %in% c("1", "2", "3", "S", "I", "R", NA_character_))) {
    x[x.bak == "1"] <- "S"
    x[x.bak == "2"] <- "I"
    x[x.bak == "3"] <- "R"
  } else if (!all(is.na(x)) && !identical(levels(x), c("R", "S", "I")) && !all(x %in% c("R", "S", "I", NA))) {
    if (all(x %unlike% "(R|S|I)", na.rm = TRUE)) {
      # check if they are actually MICs or disks
      if (all_valid_mics(x)) {
        warning_("in `as.sir()`: the input seems to contain MIC values. You can transform them with `as.mic()` before running `as.sir()` to interpret them.")
      } else if (all_valid_disks(x)) {
        warning_("in `as.sir()`: the input seems to contain disk diffusion values. You can transform them with `as.disk()` before running `as.sir()` to interpret them.")
      }
    }

    # trim leading and trailing spaces, new lines, etc.
    x <- trimws2(as.character(unlist(x)))
    x[x %in% c(NA, "", "-", "NULL")] <- NA_character_
    x.bak <- x
    
    na_before <- length(x[is.na(x)])

    # correct for translations
    trans_R <- unlist(TRANSLATIONS[
      which(TRANSLATIONS$pattern == "Resistant"),
      LANGUAGES_SUPPORTED[LANGUAGES_SUPPORTED %in% colnames(TRANSLATIONS)]
    ])
    trans_S <- unlist(TRANSLATIONS[
      which(TRANSLATIONS$pattern == "Susceptible"),
      LANGUAGES_SUPPORTED[LANGUAGES_SUPPORTED %in% colnames(TRANSLATIONS)]
    ])
    trans_I <- unlist(TRANSLATIONS[
      which(TRANSLATIONS$pattern %in% c("Incr. exposure", "Susceptible, incr. exp.", "Intermediate")),
      LANGUAGES_SUPPORTED[LANGUAGES_SUPPORTED %in% colnames(TRANSLATIONS)]
    ])
    x <- gsub(paste0(unique(trans_R[!is.na(trans_R)]), collapse = "|"), "R", x, ignore.case = TRUE)
    x <- gsub(paste0(unique(trans_S[!is.na(trans_S)]), collapse = "|"), "S", x, ignore.case = TRUE)
    x <- gsub(paste0(unique(trans_I[!is.na(trans_I)]), collapse = "|"), "I", x, ignore.case = TRUE)
    # replace all English textual input
    x[x %like% "([^a-z]|^)res(is(tant)?)?"] <- "R"
    x[x %like% "([^a-z]|^)sus(cep(tible)?)?"] <- "S"
    x[x %like% "([^a-z]|^)int(er(mediate)?)?|incr.*exp"] <- "I"
    # remove other invalid characters
    # set to capitals
    x <- toupper(x)
    x <- gsub("[^A-Z]+", "", x, perl = TRUE)
    # CLSI uses SDD for "susceptible dose-dependent"
    x <- gsub("SDD", "I", x, fixed = TRUE)
    # some labs now report "H" instead of "I" to not interfere with EUCAST prior to 2019
    x <- gsub("H", "I", x, fixed = TRUE)
    # MIPS uses D for Dose-dependent (which is I, but it will throw a note)
    x <- gsub("D", "I", x, fixed = TRUE)
    # MIPS uses U for "susceptible urine"
    x <- gsub("U", "S", x, fixed = TRUE)
    # in cases of "S;S" keep S, but in case of "S;I" make it NA
    x <- gsub("^S+$", "S", x)
    x <- gsub("^I+$", "I", x)
    x <- gsub("^R+$", "R", x)
    x[!x %in% c("S", "I", "R")] <- NA_character_
    na_after <- length(x[is.na(x) | x == ""])

    if (!isFALSE(list(...)$warn)) { # so as.sir(..., warn = FALSE) will never throw a warning
      if (na_before != na_after) {
        list_missing <- x.bak[is.na(x) & !is.na(x.bak) & x.bak != ""] %pm>%
          unique() %pm>%
          sort() %pm>%
          vector_and(quotes = TRUE)
        cur_col <- get_current_column()
        warning_("in `as.sir()`: ", na_after - na_before, " result",
          ifelse(na_after - na_before > 1, "s", ""),
          ifelse(is.null(cur_col), "", paste0(" in column '", cur_col, "'")),
          " truncated (",
          round(((na_after - na_before) / length(x)) * 100),
          "%) that were invalid antimicrobial interpretations: ",
          list_missing,
          call = FALSE
        )
      }
      if (any(toupper(x.bak[!is.na(x.bak)]) == "U") && message_not_thrown_before("as.sir", "U")) {
        warning_("in `as.sir()`: 'U' was interpreted as 'S', following some laboratory systems")
      }
      if (any(toupper(x.bak[!is.na(x.bak)]) == "D") && message_not_thrown_before("as.sir", "D")) {
        warning_("in `as.sir()`: 'D' (dose-dependent) was interpreted as 'I', following some laboratory systems")
      }
      if (any(toupper(x.bak[!is.na(x.bak)]) == "SDD") && message_not_thrown_before("as.sir", "SDD")) {
        warning_("in `as.sir()`: 'SDD' (susceptible dose-dependent, coined by CLSI) was interpreted as 'I' to comply with EUCAST's 'I'")
      }
      if (any(toupper(x.bak[!is.na(x.bak)]) == "H") && message_not_thrown_before("as.sir", "H")) {
        warning_("in `as.sir()`: 'H' was interpreted as 'I', following some laboratory systems")
      }
    }
  }

  set_clean_class(factor(x, levels = c("S", "I", "R"), ordered = TRUE),
    new_class = c("sir", "ordered", "factor")
  )
}

#' @rdname as.sir
#' @export
as.sir.mic <- function(x,
                       mo = NULL,
                       ab = deparse(substitute(x)),
                       guideline = getOption("AMR_guideline", "EUCAST"),
                       uti = NULL,
                       conserve_capped_values = FALSE,
                       add_intrinsic_resistance = FALSE,
                       reference_data = AMR::clinical_breakpoints,
                       ...) {
  as_sir_method(
    method_short = "mic",
    method_long = "MIC values",
    x = x,
    mo = mo,
    ab = ab,
    guideline = guideline,
    uti = uti,
    conserve_capped_values = conserve_capped_values,
    add_intrinsic_resistance = add_intrinsic_resistance,
    reference_data = reference_data,
    ...
  )
}

#' @rdname as.sir
#' @export
as.sir.disk <- function(x,
                        mo = NULL,
                        ab = deparse(substitute(x)),
                        guideline = getOption("AMR_guideline", "EUCAST"),
                        uti = NULL,
                        add_intrinsic_resistance = FALSE,
                        reference_data = AMR::clinical_breakpoints,
                        ...) {
  as_sir_method(
    method_short = "disk",
    method_long = "disk diffusion zones",
    x = x,
    mo = mo,
    ab = ab,
    guideline = guideline,
    uti = uti,
    conserve_capped_values = FALSE,
    add_intrinsic_resistance = add_intrinsic_resistance,
    reference_data = reference_data,
    ...
  )
}

#' @rdname as.sir
#' @export
as.sir.data.frame <- function(x,
                              ...,
                              col_mo = NULL,
                              guideline = getOption("AMR_guideline", "EUCAST"),
                              uti = NULL,
                              conserve_capped_values = FALSE,
                              add_intrinsic_resistance = FALSE,
                              reference_data = AMR::clinical_breakpoints) {
  meet_criteria(x, allow_class = "data.frame") # will also check for dimensions > 0
  meet_criteria(col_mo, allow_class = "character", is_in = colnames(x), allow_NULL = TRUE)
  meet_criteria(guideline, allow_class = "character", has_length = 1)
  meet_criteria(uti, allow_class = c("logical", "character"), allow_NULL = TRUE, allow_NA = TRUE)
  meet_criteria(conserve_capped_values, allow_class = "logical", has_length = 1)
  meet_criteria(add_intrinsic_resistance, allow_class = "logical", has_length = 1)
  meet_criteria(reference_data, allow_class = "data.frame")

  x.bak <- x
  for (i in seq_len(ncol(x))) {
    # don't keep factors, overwriting them is hard
    if (is.factor(x[, i, drop = TRUE])) {
      x[, i] <- as.character(x[, i, drop = TRUE])
    }
  }

  # -- MO
  col_mo.bak <- col_mo
  if (is.null(col_mo)) {
    col_mo <- search_type_in_df(x = x, type = "mo", info = FALSE)
  }

  # -- UTIs
  col_uti <- uti
  if (is.null(col_uti)) {
    col_uti <- search_type_in_df(x = x, type = "uti")
  }
  if (!is.null(col_uti)) {
    if (is.logical(col_uti)) {
      # already a [logical] vector as input
      if (length(col_uti) == 1) {
        uti <- rep(col_uti, NROW(x))
      } else {
        uti <- col_uti
      }
    } else {
      # column found, transform to logical
      stop_if(
        length(col_uti) != 1 | !col_uti %in% colnames(x),
        "argument `uti` must be a [logical] vector, of must be a single column name of `x`"
      )
      uti <- as.logical(x[, col_uti, drop = TRUE])
    }
  } else {
    # col_uti is still NULL - look for specimen column and make logicals of the urines
    col_specimen <- suppressMessages(search_type_in_df(x = x, type = "specimen"))
    if (!is.null(col_specimen)) {
      uti <- x[, col_specimen, drop = TRUE] %like% "urin"
      values <- sort(unique(x[uti, col_specimen, drop = TRUE]))
      if (length(values) > 1) {
        plural <- c("s", "", "")
      } else {
        plural <- c("", "s", "a ")
      }
      message_(
        "Assuming value", plural[1], " ",
        vector_and(values, quotes = TRUE),
        " in column '", font_bold(col_specimen),
        "' reflect", plural[2], " ", plural[3], "urinary tract infection", plural[1],
        ".\n  Use `as.sir(uti = FALSE)` to prevent this."
      )
    } else {
      # no data about UTI's found
      uti <- NULL
    }
  }

  i <- 0
  if (tryCatch(length(list(...)) > 0, error = function(e) TRUE)) {
    sel <- colnames(pm_select(x, ...))
  } else {
    sel <- colnames(x)
  }
  if (!is.null(col_mo)) {
    sel <- sel[sel != col_mo]
  }

  ab_cols <- colnames(x)[vapply(FUN.VALUE = logical(1), x, function(y) {
    i <<- i + 1
    check <- is.mic(y) | is.disk(y)
    ab <- colnames(x)[i]
    if (!is.null(col_mo) && ab == col_mo) {
      return(FALSE)
    }
    if (!is.null(col_uti) && ab == col_uti) {
      return(FALSE)
    }
    if (length(sel) == 0 || (length(sel) > 0 && ab %in% sel)) {
      ab_coerced <- suppressWarnings(as.ab(ab))
      if (is.na(ab_coerced) || (length(sel) > 0 & !ab %in% sel)) {
        # not even a valid AB code
        return(FALSE)
      } else {
        return(TRUE)
      }
    } else {
      return(FALSE)
    }
  })]

  stop_if(
    length(ab_cols) == 0,
    "no columns with MIC values, disk zones or antibiotic column names found in this data set. Use as.mic() or as.disk() to transform antimicrobial columns."
  )
  # set type per column
  types <- character(length(ab_cols))
  types[vapply(FUN.VALUE = logical(1), x.bak[, ab_cols, drop = FALSE], is.disk)] <- "disk"
  types[vapply(FUN.VALUE = logical(1), x.bak[, ab_cols, drop = FALSE], is.mic)] <- "mic"
  types[types == "" & vapply(FUN.VALUE = logical(1), x[, ab_cols, drop = FALSE], all_valid_disks)] <- "disk"
  types[types == "" & vapply(FUN.VALUE = logical(1), x[, ab_cols, drop = FALSE], all_valid_mics)] <- "mic"
  types[types == "" & !vapply(FUN.VALUE = logical(1), x.bak[, ab_cols, drop = FALSE], is.sir)] <- "sir"
  if (any(types %in% c("mic", "disk"), na.rm = TRUE)) {
    # now we need an mo column
    stop_if(is.null(col_mo), "`col_mo` must be set")
    # if not null, we already found it, now find again so a message will show
    if (is.null(col_mo.bak)) {
      col_mo <- search_type_in_df(x = x, type = "mo")
    }
    x_mo <- as.mo(x[, col_mo, drop = TRUE])
  }

  for (i in seq_len(length(ab_cols))) {
    if (types[i] == "mic") {
      x[, ab_cols[i]] <- x %pm>%
        pm_pull(ab_cols[i]) %pm>%
        as.character() %pm>%
        as.mic() %pm>%
        as.sir(
          mo = x_mo,
          mo.bak = x[, col_mo, drop = TRUE],
          ab = ab_cols[i],
          guideline = guideline,
          uti = uti,
          conserve_capped_values = conserve_capped_values,
          add_intrinsic_resistance = add_intrinsic_resistance,
          reference_data = reference_data,
          is_data.frame = TRUE
        )
    } else if (types[i] == "disk") {
      x[, ab_cols[i]] <- x %pm>%
        pm_pull(ab_cols[i]) %pm>%
        as.character() %pm>%
        as.disk() %pm>%
        as.sir(
          mo = x_mo,
          mo.bak = x[, col_mo, drop = TRUE],
          ab = ab_cols[i],
          guideline = guideline,
          uti = uti,
          add_intrinsic_resistance = add_intrinsic_resistance,
          reference_data = reference_data,
          is_data.frame = TRUE
        )
    } else if (types[i] == "sir") {
      show_message <- FALSE
      ab <- ab_cols[i]
      ab_coerced <- suppressWarnings(as.ab(ab))
      if (!all(x[, ab_cols[i], drop = TRUE] %in% c("R", "S", "I", NA), na.rm = TRUE)) {
        show_message <- TRUE
        # only print message if values are not already clean
        message_("=> Cleaning values in column '", font_bold(ab), "' (",
          ifelse(ab_coerced != toupper(ab), paste0(ab_coerced, ", "), ""),
          ab_name(ab_coerced, tolower = TRUE), ")... ",
          appendLF = FALSE,
          as_note = FALSE
        )
      } else if (!is.sir(x.bak[, ab_cols[i], drop = TRUE])) {
        show_message <- TRUE
        # only print message if class not already set
        message_("=> Assigning class 'sir' to already clean column '", font_bold(ab), "' (",
          ifelse(ab_coerced != toupper(ab), paste0(ab_coerced, ", "), ""),
          ab_name(ab_coerced, tolower = TRUE, language = NULL), ")... ",
          appendLF = FALSE,
          as_note = FALSE
        )
      }
      x[, ab_cols[i]] <- as.sir.default(x = as.character(x[, ab_cols[i], drop = TRUE]))
      if (show_message == TRUE) {
        message_(" OK.", add_fn = list(font_green), as_note = FALSE)
      }
    }
  }

  x
}

get_guideline <- function(guideline, reference_data) {
  if (!identical(reference_data, AMR::clinical_breakpoints)) {
    return(guideline)
  }
  guideline_param <- toupper(guideline)
  if (guideline_param %in% c("CLSI", "EUCAST")) {
    guideline_param <- rev(sort(subset(reference_data, guideline %like% guideline_param)$guideline))[1L]
  }
  if (guideline_param %unlike% " ") {
    # like 'EUCAST2020', should be 'EUCAST 2020'
    guideline_param <- gsub("([a-z]+)([0-9]+)", "\\1 \\2", guideline_param, ignore.case = TRUE)
  }

  stop_ifnot(guideline_param %in% reference_data$guideline,
    "invalid guideline: '", guideline,
    "'.\nValid guidelines are: ", vector_and(reference_data$guideline, quotes = TRUE, reverse = TRUE),
    call = FALSE
  )

  guideline_param
}

as_sir_method <- function(method_short,
                          method_long,
                          x,
                          mo,
                          ab,
                          guideline,
                          uti,
                          conserve_capped_values,
                          add_intrinsic_resistance,
                          reference_data,
                          ...) {
  meet_criteria(x, allow_NA = TRUE, .call_depth = -2)
  meet_criteria(mo, allow_class = c("mo", "character"), allow_NULL = TRUE, .call_depth = -2)
  meet_criteria(ab, allow_class = c("ab", "character"), has_length = 1, .call_depth = -2)
  meet_criteria(guideline, allow_class = "character", has_length = 1, .call_depth = -2)
  meet_criteria(uti, allow_class = "logical", has_length = c(1, length(x)), allow_NULL = TRUE, allow_NA = TRUE, .call_depth = -2)
  meet_criteria(conserve_capped_values, allow_class = "logical", has_length = 1, .call_depth = -2)
  meet_criteria(add_intrinsic_resistance, allow_class = "logical", has_length = 1, .call_depth = -2)
  meet_criteria(reference_data, allow_class = "data.frame", .call_depth = -2)
  check_reference_data(reference_data)

  # for dplyr's across()
  cur_column_dplyr <- import_fn("cur_column", "dplyr", error_on_fail = FALSE)
  if (!is.null(cur_column_dplyr) && tryCatch(is.data.frame(get_current_data("ab", call = 0)), error = function(e) FALSE)) {
    # try to get current column, which will only be available when in across()
    ab <- tryCatch(cur_column_dplyr(),
      error = function(e) ab
    )
  }

  # for auto-determining mo
  mo_var_found <- ""
  if (is.null(mo)) {
    tryCatch(
      {
        df <- get_current_data(arg_name = "mo", call = -3) # will return an error if not found
        mo <- NULL
        try(
          {
            mo <- suppressMessages(search_type_in_df(df, "mo"))
          },
          silent = TRUE
        )
        if (!is.null(df) && !is.null(mo) && is.data.frame(df)) {
          mo_var_found <- paste0(" based on column '", font_bold(mo), "'")
          mo <- df[, mo, drop = TRUE]
        }
      },
      error = function(e) {
        mo <- NULL
      }
    )
  }
  if (is.null(mo)) {
    stop_("No information was supplied about the microorganisms (missing argument `mo` and no column of class 'mo' found). See ?as.sir.\n\n",
      "To transform certain columns with e.g. mutate(), use `data %>% mutate(across(..., as.sir, mo = x))`, where x is your column with microorganisms.\n",
      "To tranform all ", method_long, " in a data set, use `data %>% as.sir()` or `data %>% mutate_if(is.", method_short, ", as.sir)`.",
      call = FALSE
    )
  }

  if (length(ab) == 1 && ab %like% paste0("as.", method_short)) {
    stop_("No unambiguous name was supplied about the antibiotic (argument `ab`). See ?as.sir.", call = FALSE)
  }

  ab.bak <- ab
  ab <- suppressWarnings(as.ab(ab))
  if (!is.null(list(...)$mo.bak)) {
    mo.bak <- list(...)$mo.bak
  } else {
    mo.bak <- mo
  }
  # be sure to take current taxonomy, as the clinical_breakpoints data set only contains current taxonomy
  mo <- suppressWarnings(suppressMessages(as.mo(mo, keep_synonyms = FALSE, inf0 = FALSE)))
  guideline_coerced <- get_guideline(guideline, reference_data)
  if (is.na(ab)) {
    message_("Returning NAs for unknown drug: '", font_bold(ab.bak),
      "'. Rename this column to a drug name or code, and check the output with `as.ab()`.",
      add_fn = font_red,
      as_note = FALSE
    )
    return(as.sir(rep(NA, length(x))))
  }
  if (length(mo) == 1) {
    mo <- rep(mo, length(x))
  }
  if (is.null(uti)) {
    uti <- NA
  }
  if (length(uti) == 1) {
    uti <- rep(uti, length(x))
  }
  
  if (isTRUE(add_intrinsic_resistance) && guideline_coerced %unlike% "EUCAST") {
    if (message_not_thrown_before("as.sir", "intrinsic")) {
      warning_("in `as.sir()`: using 'add_intrinsic_resistance' is only useful when using EUCAST guidelines, since the rules for intrinsic resistance are based on EUCAST.")
    }
  }
  
  agent_formatted <- paste0("'", font_bold(ab.bak), "'")
  agent_name <- ab_name(ab, tolower = TRUE, language = NULL)
  if (generalise_antibiotic_name(ab.bak) == generalise_antibiotic_name(agent_name)) {
    agent_formatted <- paste0(
      agent_formatted,
      " (", ab, ")"
    )
  } else if (generalise_antibiotic_name(ab) != generalise_antibiotic_name(agent_name)) {
    agent_formatted <- paste0(
      agent_formatted,
      " (", ifelse(ab.bak == ab, "",
        paste0(ab, ", ")
      ), agent_name, ")"
    )
  }
  message_("=> Interpreting ", method_long, " of ", ifelse(isTRUE(list(...)$is_data.frame), "column ", ""),
    agent_formatted,
    mo_var_found,
    " according to ", ifelse(identical(reference_data, AMR::clinical_breakpoints),
      font_bold(guideline_coerced),
      "manually defined 'reference_data'"
    ),
    "... ",
    appendLF = FALSE,
    as_note = FALSE
  )
  
  msg_note <- function(messages) {
    for (i in seq_len(length(messages))) {
      messages[i] <- word_wrap(extra_indent = 5, messages[i])
    }
    message(font_green(font_bold(" Note:\n")),
            paste0("   ", font_black(AMR_env$bullet_icon)," ", font_black(messages, collapse = NULL) , collapse = "\n"))
  }

  method <- method_short

  metadata_mo <- get_mo_uncertainties()

  df <- data.frame(values = x,
                   mo = mo,
                   result = NA_sir_,
                   uti = uti,
                   stringsAsFactors = FALSE)
  if (method == "mic") {
    # when as.sir.mic is called directly
    df$values <- as.mic(df$values) 
  } else if (method == "disk") {
    # when as.sir.disk is called directly
    df$values <- as.disk(df$values)
  }

  rise_warning <- FALSE
  rise_note <- FALSE
  method_coerced <- toupper(method)
  ab_coerced <- ab
  mo_coerced <- mo
  
  if (identical(reference_data, AMR::clinical_breakpoints)) {
    breakpoints <- reference_data %pm>%
      subset(guideline == guideline_coerced & method == method_coerced & ab == ab_coerced)
    if (ab_coerced == "AMX" && nrow(breakpoints) == 0) {
      ab_coerced <- "AMP"
      breakpoints <- reference_data %pm>%
        subset(guideline == guideline_coerced & method == method_coerced & ab == ab_coerced)
    }
  } else {
    breakpoints <- reference_data %pm>%
      subset(method == method_coerced & ab == ab_coerced)
  }
  
  msgs <- character(0)
  if (nrow(breakpoints) == 0) {
    # apparently no breakpoints found
    msg_note(paste0("No ", method_coerced, " breakpoints available for ",
                    suppressMessages(suppressWarnings(ab_name(ab_coerced, language = NULL, tolower = TRUE))),
                    " (", ab_coerced, ")"))
    load_mo_uncertainties(metadata_mo)
    return(rep(NA_sir_, nrow(df)))
  }
  
  if (guideline_coerced %like% "EUCAST") {
    any_is_intrinsic_resistant <- FALSE
    add_intrinsic_resistance_to_AMR_env()
  }
  
  # run the rules
  for (mo_unique in unique(df$mo)) {
    
    rows <- which(df$mo == mo_unique)
    values <- df[rows, "values", drop = TRUE]
    uti <- df[rows, "uti", drop = TRUE]
    new_sir <- rep(NA_sir_, length(rows))
    
    # find different mo properties
    mo_current_genus <- as.mo(mo_genus(mo_unique, language = NULL))
    mo_current_family <- as.mo(mo_family(mo_unique, language = NULL))
    mo_current_order <- as.mo(mo_order(mo_unique, language = NULL))
    mo_current_class <- as.mo(mo_class(mo_unique, language = NULL))
    if (mo_genus(mo_unique, language = NULL) == "Staphylococcus") {
      mo_current_becker <- as.mo(mo_unique, Becker = TRUE)
    } else {
      mo_current_becker <- mo_unique
    }
    if (mo_genus(mo_unique, language = NULL) == "Streptococcus") {
      mo_current_lancefield <- as.mo(mo_unique, Lancefield = TRUE)
    } else {
      mo_current_lancefield <- mo_unique
    }
    mo_current_other <- as.mo("UNKNOWN")
    # formatted for notes
    mo_formatted <- suppressMessages(suppressWarnings(mo_fullname(mo_unique, language = NULL, keep_synonyms = FALSE)))
    if (!mo_rank(mo_unique) %in% c("kingdom", "phylum", "class", "order")) {
      mo_formatted <- font_italic(mo_formatted)
    }
    ab_formatted <- paste0(suppressMessages(suppressWarnings(ab_name(ab_coerced, language = NULL, tolower = TRUE))), 
                           " (", ab_coerced, ")")
    
    # gather all available breakpoints for current MO and sort on taxonomic rank 
    # (this will prefer species breakpoints over order breakpoints)
    breakpoints_current <- breakpoints %pm>%
      subset(mo %in% c(mo_current_genus, mo_current_family,
                       mo_current_order, mo_current_class,
                       mo_current_becker, mo_current_lancefield,
                       mo_current_other))
    
    if (any(df[rows, "uti", drop = TRUE], na.rm = TRUE)) {
      breakpoints_current <- breakpoints_current %pm>%
        # be as specific as possible (i.e. prefer species over genus):
        # the below `pm_desc(uti)` will put `TRUE` on top and FALSE on bottom
        pm_arrange(rank_index, pm_desc(uti)) # 'uti' is a column in data set 'clinical_breakpoints'
    } else {
      breakpoints_current <- breakpoints_current %pm>%
        # sort UTI = FALSE first, then UTI = TRUE
        pm_arrange(rank_index, uti)
    }
    
    # throw notes for different body sites
    if (nrow(breakpoints_current) == 1 && all(breakpoints_current$uti == TRUE) && any(uti %in% c(FALSE, NA)) && message_not_thrown_before("as.sir", "uti", ab_coerced)) {
      # only UTI breakpoints available
      warning_("in `as.sir()`: interpretation of ", font_bold(ab_formatted), " is only available for (uncomplicated) urinary tract infections (UTI) for some microorganisms, thus assuming `uti = TRUE`. See `?as.sir`.")
      rise_warning <- TRUE
    } else if (nrow(breakpoints_current) > 1 && length(unique(breakpoints_current$site)) > 1 && any(is.na(uti)) && all(c(TRUE, FALSE) %in% breakpoints_current$uti, na.rm = TRUE) && message_not_thrown_before("as.sir", "siteUTI", mo_unique, ab_coerced)) {
      # both UTI and Non-UTI breakpoints available
      msgs <- c(msgs, paste0("Breakpoints for UTI ", font_underline("and"), " non-UTI available for ", ab_formatted, " in ", mo_formatted, " - assuming non-UTI. Use argument `uti` to set which isolates are from urine. See `?as.sir`."))
      breakpoints_current <- breakpoints_current %pm>%
        pm_filter(uti == FALSE)
    } else if (nrow(breakpoints_current) > 1 && length(unique(breakpoints_current$site)) > 1 && all(breakpoints_current$uti == FALSE, na.rm = TRUE) && message_not_thrown_before("as.sir", "siteOther", mo_unique, ab_coerced)) {
      # breakpoints for multiple body sites available
      site <- breakpoints_current[1L, "site", drop = FALSE] # this is the one we'll take
      if (is.na(site)) {
        site <- paste0("an unspecified body site")
      } else {
        site <- paste0("body site '", site, "'")
      }
      msgs <- c(msgs, paste0("Multiple breakpoints available for ", ab_formatted, " in ", mo_formatted, " - assuming ", site, "."))
    }
    
    # first check if mo is intrinsic resistant
    if (isTRUE(add_intrinsic_resistance) && guideline_coerced %like% "EUCAST" && paste(mo_unique, ab_coerced) %in% AMR_env$intrinsic_resistant) {
      msgs <- c(msgs, paste0("Intrinsic resistance applied for ", ab_formatted, " in ", mo_formatted, ""))
      new_sir <- rep(as.sir("R"), length(rows))
      
    } else {
      # then run the rules
      breakpoints_current <- breakpoints_current[1L, , drop = FALSE]
      
      if (method == "mic") {
        new_sir <- quick_case_when(
          is.na(values) ~ NA_sir_,
          values <= breakpoints_current$breakpoint_S ~ as.sir("S"),
          guideline_coerced %like% "EUCAST" & values > breakpoints_current$breakpoint_R ~ as.sir("R"),
          guideline_coerced %like% "CLSI" & values >= breakpoints_current$breakpoint_R ~ as.sir("R"),
          # return "I" when breakpoints are in the middle
          !is.na(breakpoints_current$breakpoint_S) & !is.na(breakpoints_current$breakpoint_R) ~ as.sir("I"),
          # and NA otherwise
          TRUE ~ NA_sir_
        )
        
      } else if (method == "disk") {
        new_sir <- quick_case_when(
          is.na(values) ~ NA_sir_,
          as.double(values) >= as.double(breakpoints_current$breakpoint_S) ~ as.sir("S"),
          guideline_coerced %like% "EUCAST" & as.double(values) < as.double(breakpoints_current$breakpoint_R) ~ as.sir("R"),
          guideline_coerced %like% "CLSI" & as.double(values) <= as.double(breakpoints_current$breakpoint_R) ~ as.sir("R"),
          # return "I" when breakpoints are in the middle
          !is.na(breakpoints_current$breakpoint_S) & !is.na(breakpoints_current$breakpoint_R) ~ as.sir("I"),
          # and NA otherwise
          TRUE ~ NA_sir_
        )
      }

      # write to verbose output
      AMR_env$sir_interpretation_history <- rbind(
        AMR_env$sir_interpretation_history,
        # recycling 1 to 2 rows does not seem to work, which is why rep() was added
        data.frame(
          datetime = rep(Sys.time(), length(rows)),
          index = rows,
          ab_input = rep(ab.bak, length(rows)),
          ab_guideline = rep(ab_coerced, length(rows)),
          mo_input = rep(mo.bak[match(mo_unique, df$mo)][1], length(rows)),
          mo_guideline = rep(breakpoints_current[, "mo", drop = TRUE], length(rows)),
          guideline = rep(guideline_coerced, length(rows)),
          ref_table = rep(breakpoints_current[, "ref_tbl", drop = TRUE], length(rows)),
          method = rep(method_coerced, length(rows)),
          input = as.double(values),
          outcome = as.sir(new_sir),
          breakpoint_S_R = rep(paste0(breakpoints_current[, "breakpoint_S", drop = TRUE], "-", breakpoints_current[, "breakpoint_R", drop = TRUE]), length(rows)),
          stringsAsFactors = FALSE
        )
      )
    }
    
    df[rows, "result"] <- new_sir
  }
  
  if (isTRUE(rise_warning)) {
    message(font_yellow(font_bold(" * WARNING *")))
  } else if (length(msgs) == 0) {
    message(font_green(" OK."))
  } else {
    msg_note(sort(msgs))
  }
  
  load_mo_uncertainties(metadata_mo)
  
  df$result
}

#' @rdname as.sir
#' @param clean a [logical] to indicate whether previously stored results should be forgotten after returning the 'logbook' with results
#' @export
sir_interpretation_history <- function(clean = FALSE) {
  meet_criteria(clean, allow_class = "logical", has_length = 1)

  out.bak <- AMR_env$sir_interpretation_history
  out <- out.bak
  if (NROW(out) == 0) {
    message_("No results to return. Run `as.sir()` on MIC values or disk diffusion zones first to see a 'logbook' data set here.")
    return(invisible(NULL))
  }
  out$ab_guideline <- as.ab(out$ab_guideline)
  out$mo_guideline <- as.mo(out$mo_guideline)
  out$outcome <- as.sir(out$outcome)
  # keep stored for next use
  if (isTRUE(clean)) {
    AMR_env$sir_interpretation_history <- AMR_env$sir_interpretation_history[0, , drop = FALSE]
  } else {
    AMR_env$sir_interpretation_history <- out.bak
  }

  if (pkg_is_available("tibble", also_load = FALSE)) {
    import_fn("as_tibble", "tibble")(out)
  } else {
    out
  }
}

# will be exported using s3_register() in R/zzz.R
pillar_shaft.sir <- function(x, ...) {
  out <- trimws(format(x))
  if (has_colour()) {
    # colours will anyway not work when has_colour() == FALSE,
    # but then the indentation should also not be applied
    out[is.na(x)] <- font_grey("  NA")
    out[x == "S"] <- font_green_bg("  S  ")
    out[x == "I"] <- font_orange_bg("  I  ")
    out[x == "R"] <- font_red_bg("  R  ")
  }
  create_pillar_column(out, align = "left", width = 5)
}

# will be exported using s3_register() in R/zzz.R
type_sum.sir <- function(x, ...) {
  "sir"
}

# will be exported using s3_register() in R/zzz.R
freq.sir <- function(x, ...) {
  x_name <- deparse(substitute(x))
  x_name <- gsub(".*[$]", "", x_name)
  if (x_name %in% c("x", ".")) {
    # try again going through system calls
    x_name <- stats::na.omit(vapply(
      FUN.VALUE = character(1),
      sys.calls(),
      function(call) {
        call_txt <- as.character(call)
        ifelse(call_txt[1] %like% "freq$", call_txt[length(call_txt)], character(0))
      }
    ))[1L]
  }
  ab <- suppressMessages(suppressWarnings(as.ab(x_name)))
  digits <- list(...)$digits
  if (is.null(digits)) {
    digits <- 2
  }
  if (!is.na(ab)) {
    cleaner::freq.default(
      x = x, ...,
      .add_header = list(
        Drug = paste0(ab_name(ab, language = NULL), " (", ab, ", ", paste(ab_atc(ab), collapse = "/"), ")"),
        `Drug group` = ab_group(ab, language = NULL),
        `%SI` = trimws(percentage(susceptibility(x, minimum = 0, as_percent = FALSE),
          digits = digits
        ))
      )
    )
  } else {
    cleaner::freq.default(
      x = x, ...,
      .add_header = list(
        `%SI` = trimws(percentage(susceptibility(x, minimum = 0, as_percent = FALSE),
          digits = digits
        ))
      )
    )
  }
}


# will be exported using s3_register() in R/zzz.R
get_skimmers.sir <- function(column) {
  # get the variable name 'skim_variable'
  name_call <- function(.data) {
    calls <- sys.calls()
    frms <- sys.frames()
    calls_txt <- vapply(calls, function(x) paste(deparse(x), collapse = ""), FUN.VALUE = character(1))
    if (any(calls_txt %like% "skim_variable", na.rm = TRUE)) {
      ind <- which(calls_txt %like% "skim_variable")[1L]
      vars <- tryCatch(eval(parse(text = ".data$skim_variable$sir"), envir = frms[[ind]]),
        error = function(e) NULL
      )
      tryCatch(ab_name(as.character(calls[[length(calls)]][[2]]), language = NULL),
        error = function(e) NA_character_
      )
    } else {
      NA_character_
    }
  }

  skimr::sfl(
    skim_type = "sir",
    ab_name = name_call,
    count_R = count_R,
    count_S = count_susceptible,
    count_I = count_I,
    prop_R = ~ proportion_R(., minimum = 0),
    prop_S = ~ susceptibility(., minimum = 0),
    prop_I = ~ proportion_I(., minimum = 0)
  )
}

#' @method print sir
#' @export
#' @noRd
print.sir <- function(x, ...) {
  cat("Class 'sir'\n")
  print(as.character(x), quote = FALSE)
}

#' @method droplevels sir
#' @export
#' @noRd
droplevels.sir <- function(x, exclude = if (any(is.na(levels(x)))) NULL else NA, ...) {
  x <- droplevels.factor(x, exclude = exclude, ...)
  class(x) <- c("sir", "ordered", "factor")
  x
}

#' @method summary sir
#' @export
#' @noRd
summary.sir <- function(object, ...) {
  x <- object
  n <- sum(!is.na(x))
  S <- sum(x == "S", na.rm = TRUE)
  I <- sum(x == "I", na.rm = TRUE)
  R <- sum(x == "R", na.rm = TRUE)
  pad <- function(x) {
    if (x == "0%") {
      x <- " 0.0%"
    }
    if (nchar(x) < 5) {
      x <- paste0(rep(" ", 5 - nchar(x)), x)
    }
    x
  }
  value <- c(
    "Class" = "sir",
    "%R" = paste0(pad(percentage(R / n, digits = 1)), " (n=", R, ")"),
    "%SI" = paste0(pad(percentage((S + I) / n, digits = 1)), " (n=", S + I, ")"),
    "- %S" = paste0(pad(percentage(S / n, digits = 1)), " (n=", S, ")"),
    "- %I" = paste0(pad(percentage(I / n, digits = 1)), " (n=", I, ")")
  )
  class(value) <- c("summaryDefault", "table")
  value
}

#' @method [<- sir
#' @export
#' @noRd
"[<-.sir" <- function(i, j, ..., value) {
  value <- as.sir(value)
  y <- NextMethod()
  attributes(y) <- attributes(i)
  y
}
#' @method [[<- sir
#' @export
#' @noRd
"[[<-.sir" <- function(i, j, ..., value) {
  value <- as.sir(value)
  y <- NextMethod()
  attributes(y) <- attributes(i)
  y
}
#' @method c sir
#' @export
#' @noRd
c.sir <- function(...) {
  as.sir(unlist(lapply(list(...), as.character)))
}

#' @method unique sir
#' @export
#' @noRd
unique.sir <- function(x, incomparables = FALSE, ...) {
  y <- NextMethod()
  attributes(y) <- attributes(x)
  y
}

#' @method rep sir
#' @export
#' @noRd
rep.sir <- function(x, ...) {
  y <- NextMethod()
  attributes(y) <- attributes(x)
  y
}

check_reference_data <- function(reference_data) {
  if (!identical(reference_data, AMR::clinical_breakpoints)) {
    class_sir <- vapply(FUN.VALUE = character(1), clinical_breakpoints, function(x) paste0("<", class(x), ">", collapse = " and "))
    class_ref <- vapply(FUN.VALUE = character(1), reference_data, function(x) paste0("<", class(x), ">", collapse = " and "))
    if (!all(names(class_sir) == names(class_ref))) {
      stop_("`reference_data` must have the same column names as the 'clinical_breakpoints' data set.", call = -2)
    }
    if (!all(class_sir == class_ref)) {
      class_sir[class_sir != class_ref][1]
      stop_("`reference_data` must be the same structure as the 'clinical_breakpoints' data set. Column '", names(class_ref[class_sir != class_ref][1]), "' is of class ", class_ref[class_sir != class_ref][1], ", but should be of class ", class_sir[class_sir != class_ref][1], ".", call = -2)
    }
  }
}
