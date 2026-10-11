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

#' Combine SIR Results of Multiple Antimicrobials
#'
#' @description Combine the SIR results of two or more antimicrobials into one SIR result per isolate, e.g. to create a variable for the combination therapy of piperacillin/tazobactam with gentamicin. The result is the best result of the antimicrobials, as a combination therapy covers the isolate if any of its antimicrobials is active. See *Combination Therapy* below for the exact rules.
#'
#' The same rules are used by [antibiogram()] and [wisca()] for combination antibiograms, and by the [`count_*()`][count()] and [`proportion_*()`][proportion()] functions when more than one variable is used, so the result of [sir_combine()] can be used in any further analysis.
#' @param ... Two or more vectors (or columns) with antimicrobial interpretations, or a data frame with such columns. They will be transformed internally with [as.sir()] if needed. All vectors must be of the same length.
#' @param only_all_tested A [logical] to indicate that isolates must be tested for all antimicrobials. If `TRUE`, the result is `NA` for isolates with a missing result for any of the antimicrobials. If `FALSE` (default), a missing result only leads to `NA` if a result for the missing antimicrobial could have changed the result of the combination, see *Combination Therapy* below.
#' @param guideline Either `"EUCAST"` (default) or `"CLSI"`, to determine whether the 'I' category is considered susceptible. This only affects isolates with a missing result for any of the antimicrobials when `only_all_tested = FALSE`, see *Combination Therapy* below. The default can be changed with the `AMR_guideline` option, see [AMR-options].
#' @section Combination Therapy:
#' A combination of antimicrobials is considered as good as its most active antimicrobial: `"S"` is preferred over `"WT"`, and `"R"` over `"NWT"` and `"NS"`. The value `"NI"` (non-interpretable) is considered a missing result.
#'
#' A missing result can only be ignored if the result of the combination is susceptible regardless of the missing result, which depends on the `only_all_tested` argument and on the guideline: under EUCAST, the 'I' category is considered susceptible (*susceptible, increased exposure*), and under CLSI, it is not.
#'
#' This is an example for two antimicrobials, Drug A and Drug B (the order of the antimicrobials does not matter):
#'
#' ```
#' -------------------------------------------------------------------------
#'                             Result of the combination
#'                     -----------------------------------------------------
#'  Drug A    Drug B   only_all_tested = FALSE       only_all_tested = TRUE
#' --------  --------  ----------------------------  ----------------------
#'     S     S, I, R   S                             S
#'     I       I, R    I                             I
#'     R         R     R                             R
#'     S       <NA>    S                             <NA>
#'     I       <NA>    I (EUCAST), <NA> (CLSI)       <NA>
#'     R       <NA>    <NA>                          <NA>
#'   <NA>      <NA>    <NA>                          <NA>
#' -------------------------------------------------------------------------
#' ```
#'
#' The [`count_*()`][count()] and [`proportion_*()`][proportion()] functions, such as [susceptibility()] and [resistance()], apply exactly these rules when more than one variable is used: an isolate is considered tested if the result of the combination is not `NA`, and is then counted by the result of the combination. Consequently, [count_susceptible()] `+` [count_resistant()] `=` [count_all()], and [susceptibility()] `+` [resistance()] `= 1`, for both values of `only_all_tested`. The guideline is taken from their `guideline` argument if they have one, and otherwise from the `AMR_guideline` option.
#'
#' In [antibiogram()] and [wisca()], the 'I' and 'SDD' categories are considered susceptible if `combine_SI = TRUE` (default), and the result of the combination is then reduced to S or R.
#' @return A [`sir`] vector with the same length as the input.
#' @seealso [count()] and [proportion()] for counting and calculating susceptibility of (combinations of) antimicrobials, and [antibiogram()] for combination antibiograms.
#' @export
#' @examples
#' drug_a <- as.sir(c("S", "R", "R", "I", "I", NA, "R", "SDD"))
#' drug_b <- as.sir(c("R", "S", "R", "R", NA, NA, NA, "I"))
#'
#' sir_combine(drug_a, drug_b)
#' sir_combine(drug_a, drug_b, guideline = "CLSI")
#' sir_combine(drug_a, drug_b, only_all_tested = TRUE)
#'
#' # the count and proportion functions use the same rules
#' susceptibility(drug_a, drug_b, minimum = 0)
#' susceptibility(sir_combine(drug_a, drug_b), minimum = 0)
#'
#' \donttest{
#' if (require("dplyr")) {
#'   # create a new variable for piperacillin/tazobactam + gentamicin
#'   example_isolates %>%
#'     filter(!is.na(TZP)) %>%
#'     mutate(TZP_GEN = sir_combine(TZP, GEN)) %>%
#'     select(mo, TZP, GEN, TZP_GEN)
#' }
#' if (require("dplyr")) {
#'   # combine any number of antimicrobials
#'   example_isolates %>%
#'     mutate(any_aminoglycoside = sir_combine(GEN, TOB, AMK, KAN)) %>%
#'     count(any_aminoglycoside)
#' }
#' if (require("dplyr")) {
#'   # use the new variable in further analyses
#'   # THOUGH CONSIDER wisca() FOR RELIABLE ESTIMATES
#'   example_isolates %>%
#'     mutate(TZP_GEN = sir_combine(TZP, GEN)) %>%
#'     group_by(ward) %>%
#'     summarise(
#'       TZP_mono = susceptibility(TZP),
#'       GEN_mono = susceptibility(GEN),
#'       TZP_GEN = susceptibility(TZP_GEN),
#'       # or directly, without combining first:
#'       TZP_GEN2 = susceptibility(TZP, GEN)
#'     )
#' }
#' }
sir_combine <- function(..., only_all_tested = FALSE, guideline = getOption("AMR_guideline", "EUCAST")) {
  meet_criteria(only_all_tested, allow_class = "logical", has_length = 1)
  meet_criteria(guideline, allow_class = "character", has_length = 1)
  if (guideline %like% "EUCAST") {
    guideline <- "EUCAST"
  } else if (guideline %like% "CLSI") {
    guideline <- "CLSI"
  }
  meet_criteria(guideline, allow_class = "character", is_in = c("EUCAST", "CLSI"), has_length = 1)

  # unpack data frames, so that e.g. sir_combine(df[, c("GEN", "TOB")], CIP) works
  cols <- list()
  for (dot in list(...)) {
    if (is.data.frame(dot)) {
      cols <- c(cols, unname(as.list(dot)))
    } else {
      cols <- c(cols, list(dot))
    }
  }
  stop_if(length(cols) < 2, "at least two vectors or columns must be provided to combine, but ", length(cols), " ", ifelse(length(cols) == 1, "was", "were"), " provided.")
  col_lengths <- vapply(FUN.VALUE = integer(1), cols, length)
  stop_if(length(unique(col_lengths)) > 1, "all vectors or columns must be of the same length, but they have lengths ", vector_and(sort(unique(col_lengths)), quotes = FALSE), ".")
  cols <- lapply(cols, function(x) if (is.sir(x)) x else as.sir(x))

  out <- sir_combine_engine(cols, S_values = get_susceptible_values(guideline), only_all_tested = only_all_tested)

  if (isFALSE(only_all_tested) && is.null(getOption("AMR_guideline")) && missing(guideline)) {
    # only notify about the assumed guideline if it affected the result
    any_missing <- Reduce(`|`, lapply(cols, function(x) is.na(x) | as.character(x) == "NI"))
    if (any(any_missing & as.character(out) %in% "I") && message_not_thrown_before("sir_combine", "eucast_default", entire_session = TRUE)) {
      message_("{.help [{.fun sir_combine}](AMR::sir_combine)} assumes the EUCAST guideline and thus considers the 'I' category susceptible, so a missing result does not change a combination with an 'I' result. Set the {.arg guideline} argument or the {.code AMR_guideline} option to either \"CLSI\" or \"EUCAST\", see {.topic [AMR-options](AMR::AMR-options)}.")
      message_("This message will be shown once per session.")
    }
  }
  out
}

get_susceptible_values <- function(guideline) {
  # the SIR values considered susceptible, as in count_susceptible() and susceptibility()
  if (!is.null(guideline) && guideline %like% "CLSI") {
    c("S", "SDD", "WT")
  } else {
    c("S", "SDD", "I", "WT")
  }
}

sir_combine_engine <- function(cols, S_values, only_all_tested = FALSE) {
  # this is the engine for combination therapy, used by sir_combine(), sir_calc() and antibiogram()
  # - cols: a list of vectors of equal length, with class 'sir' or with SIR values as character
  # - S_values: the values considered susceptible, so that a missing result cannot change the result of the combination
  # the order of preference: the first value is the best result, values with equal rank are resolved by this order
  preference <- c("S", "WT", "SDD", "I", "R", "NWT", "NS")
  ranks <- c(1L, 1L, 2L, 3L, 4L, 4L, 4L)
  cols <- lapply(cols, as.character)
  if (length(cols) == 0) {
    return(as_sir_structure(character(0)))
  }
  # 'NI' and NA get no rank, so they are missing results
  col_ranks <- lapply(cols, function(x) ranks[match(x, preference)])
  best <- do.call(pmin, c(col_ranks, list(na.rm = TRUE)))
  out <- rep(NA_character_, length(best))
  for (i in seq_along(preference)) {
    is_present <- Reduce(`|`, lapply(cols, function(x) x %in% preference[i]))
    out[is.na(out) & best %in% ranks[i] & is_present] <- preference[i]
  }
  any_missing <- Reduce(`|`, lapply(col_ranks, is.na))
  if (isTRUE(only_all_tested)) {
    out[any_missing] <- NA_character_
  } else {
    # a missing result could have been better, unless the combination is already susceptible
    out[any_missing & !out %in% S_values] <- NA_character_
  }
  as_sir_structure(out)
}
