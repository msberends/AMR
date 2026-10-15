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
#' @description Combine the SIR results of two or more antimicrobials into one SIR result per isolate. This can answer two different questions, set with the `method` argument:
#'
#' * `"most_susceptible"` (default): *Is at least one of these antimicrobials expected to work?* This is the result of a combination therapy, e.g. piperacillin/tazobactam with gentamicin, which covers the isolate if any of its antimicrobials is active.
#' * `"most_resistant"`: *Is resistance demonstrated to any of these antimicrobials?* This is the result for a group of antimicrobials, e.g. resistance to third-generation cephalosporins if an isolate is resistant to ceftazidime or ceftriaxone, as used in surveillance reports.
#'
#' See *Combination Therapy* below for the exact rules. The same rules are used by [antibiogram()] and [wisca()] for combination antibiograms, and by the [`count_*()`][count()] and [`proportion_*()`][proportion()] functions when more than one variable is used.
#' @param ... Two or more vectors (or columns) with antimicrobial interpretations, or a data frame with such columns. They will be transformed internally with [as.sir()] if needed. All vectors must be of the same length.
#' @param method Either `"most_susceptible"` (default) to return the most susceptible result, as for combination therapy, or `"most_resistant"` to return the most resistant result, as for resistance within a group of antimicrobials. See *Combination Therapy* below.
#' @param only_all_tested A [logical] to indicate that isolates must be tested for all antimicrobials. If `TRUE`, the result is `NA` for isolates with a missing result for any of the antimicrobials. If `FALSE` (default), a missing result only leads to `NA` if a result for the missing antimicrobial could have changed the result, see *Combination Therapy* below.
#' @param guideline A guideline name, `"EUCAST"` or `"CLSI"` (or any version thereof, e.g. `"CLSI 2025"`), to determine whether the 'I' category is considered susceptible (EUCAST) or resistant (CLSI). This only affects isolates with a missing result for any of the antimicrobials when `only_all_tested = FALSE`, see *Combination Therapy* below. The default is the package option [`AMR_guideline`][AMR-options], or `"EUCAST"` if it is not set. See section *Interpretation of SIR* below.
#' @inheritSection as.sir Interpretation of SIR
#' @section Combination Therapy:
#' Results are ordered from susceptible to resistant: S (and WT), SDD, I, R (and NWT, NS). With `method = "most_susceptible"` (in [sir_combine()]), the most susceptible result of the antimicrobials is returned, and with `method = "most_resistant"`, the most resistant result. The value `"NI"` (non-interpretable) is considered a missing result.
#'
#' A missing result can only be ignored if it cannot change the conclusion: a susceptible result for `"most_susceptible"`, and a non-susceptible result for `"most_resistant"`. Under EUCAST, the 'I' category is considered susceptible (*susceptible, increased exposure*), and under CLSI, it is not. With `only_all_tested = TRUE`, any missing result leads to `NA`.
#'
#' This is an example for two antimicrobials, Drug A and Drug B, under EUCAST (the order of the antimicrobials does not matter):
#'
#' ```
#' -------------------------------------------------------------------------
#'                       Result of method =            Result of method =
#'                      "most_susceptible"              "most_resistant"
#'                    -----------------------  -----------------------------
#'                     only_all_tested =        only_all_tested =
#'  Drug A    Drug B     FALSE      TRUE          FALSE         TRUE
#' --------  --------  ---------  ---------    -----------  -------------
#'     S        S          S          S             S             S
#'     S        I          S          S             I             I
#'     S        R          S          S             R             R
#'     I        R          I          I             R             R
#'     R        R          R          R             R             R
#'     S      <NA>         S        <NA>          <NA>          <NA>
#'     I      <NA>         I        <NA>          <NA>          <NA>
#'     R      <NA>       <NA>       <NA>            R           <NA>
#'   <NA>     <NA>       <NA>       <NA>          <NA>          <NA>
#' -------------------------------------------------------------------------
#' ```
#'
#' Under CLSI, the combination of I with a missing result is `<NA>` for `"most_susceptible"` and I for `"most_resistant"`.
#'
#' The [`count_*()`][count()] and [`proportion_*()`][proportion()] functions apply these rules when more than one variable is used. With `only_all_tested = FALSE` (default), the functions for susceptibility, such as [susceptibility()] and [count_susceptible()], count an isolate as susceptible if any of the antimicrobials is susceptible (`"most_susceptible"`), and the functions for resistance, such as [resistance()] and [count_resistant()], count an isolate as resistant if any of the antimicrobials is resistant (`"most_resistant"`). An isolate is considered tested if it has such a result, or if it was tested for all antimicrobials. Therefore, for combinations with `only_all_tested = FALSE`:
#'
#' ```
#'     count_S()    +   count_I()    +   count_R()    >= count_all()
#'   proportion_S() + proportion_I() + proportion_R() >= 1
#' ```
#'
#' With `only_all_tested = TRUE`, all these functions count the result of the combination therapy (`"most_susceptible"`) of isolates tested for all antimicrobials, so that:
#'
#' ```
#'     count_S()    +   count_I()    +   count_R()    = count_all()
#'   proportion_S() + proportion_I() + proportion_R() = 1
#' ```
#'
#' In [antibiogram()] and [wisca()], the result of the combination therapy (`"most_susceptible"`) is used, the 'I' category is considered susceptible or not following the `guideline` argument ('SDD' is always considered susceptible), and the result is then reduced to S or R.
#' @return A [`sir`] vector with the same length as the input.
#' @seealso [count()] and [proportion()] for counting and calculating susceptibility of (combinations of) antimicrobials, and [antibiogram()] for combination antibiograms.
#' @export
#' @examples
#' drug_a <- as.sir(c("S", "R", "R", "I", "I", NA, "R", "SDD"))
#' drug_b <- as.sir(c("R", "S", "R", "R", NA, NA, NA, "I"))
#'
#' # combination therapy: is at least one of them expected to work?
#' sir_combine(drug_a, drug_b)
#' sir_combine(drug_a, drug_b, guideline = "CLSI")
#' sir_combine(drug_a, drug_b, only_all_tested = TRUE)
#'
#' # resistance within a group: is resistance demonstrated to any of them?
#' sir_combine(drug_a, drug_b, method = "most_resistant")
#'
#' # the count and proportion functions use the same rules
#' susceptibility(drug_a, drug_b, minimum = 0)
#' susceptibility(sir_combine(drug_a, drug_b), minimum = 0)
#' resistance(drug_a, drug_b, minimum = 0)
#' resistance(sir_combine(drug_a, drug_b, method = "most_resistant"), minimum = 0)
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
#'   # resistance to third-generation cephalosporins
#'   example_isolates %>%
#'     mutate(cephalosporins_3rd = sir_combine(CAZ, CRO, method = "most_resistant")) %>%
#'     count(cephalosporins_3rd)
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
sir_combine <- function(..., method = c("most_susceptible", "most_resistant"), only_all_tested = FALSE, guideline = getOption("AMR_guideline", "EUCAST")) {
  method <- match.arg(method)
  meet_criteria(only_all_tested, allow_class = "logical", has_length = 1)
  guideline_missing <- missing(guideline)
  guideline <- validate_sir_guideline(guideline, fn = "sir_combine")

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

  out <- sir_combine_engine(cols, S_values = get_susceptible_values(guideline), only_all_tested = only_all_tested, method = method)

  if (isFALSE(only_all_tested) && is.null(getOption("AMR_guideline")) && guideline_missing) {
    # only notify about the assumed guideline if it affected the result
    out_clsi <- sir_combine_engine(cols, S_values = get_susceptible_values("CLSI"), only_all_tested = only_all_tested, method = method)
    if (!identical(out, out_clsi)) {
      message_guideline_assumed("sir_combine", txt = ", which affected the result for 'I' results combined with a missing result. Use {.code guideline = \"CLSI\"} to consider it resistant.")
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

sir_combine_engine <- function(cols, S_values, only_all_tested = FALSE, method = "most_susceptible") {
  # this is the engine for combining SIR results, used by sir_combine(), sir_calc() and antibiogram()
  # - cols: a list of vectors of equal length, with class 'sir' or with SIR values as character
  # - S_values: the values considered susceptible; all other values are considered non-susceptible
  # - method: "most_susceptible" (combination therapy) or "most_resistant" (resistance within a group)
  # the order of preference: values with equal rank are resolved by this order
  preference <- c("S", "WT", "SDD", "I", "R", "NWT", "NS")
  ranks <- c(1L, 1L, 2L, 3L, 4L, 4L, 4L)
  cols <- lapply(cols, as.character)
  if (length(cols) == 0) {
    return(as_sir_structure(character(0)))
  }
  # 'NI' and NA get no rank, so they are missing results
  col_ranks <- lapply(cols, function(x) ranks[match(x, preference)])
  if (method == "most_susceptible") {
    selected <- do.call(pmin, c(col_ranks, list(na.rm = TRUE)))
    # a missing result cannot change a susceptible result
    decisive_values <- S_values
  } else {
    selected <- do.call(pmax, c(col_ranks, list(na.rm = TRUE)))
    # a missing result cannot change a non-susceptible result
    decisive_values <- setdiff(preference, S_values)
  }
  out <- rep(NA_character_, length(selected))
  for (i in seq_along(preference)) {
    is_present <- Reduce(`|`, lapply(cols, function(x) x %in% preference[i]))
    out[is.na(out) & selected %in% ranks[i] & is_present] <- preference[i]
  }
  any_missing <- Reduce(`|`, lapply(col_ranks, is.na))
  if (isTRUE(only_all_tested)) {
    out[any_missing] <- NA_character_
  } else {
    out[any_missing & !out %in% decisive_values] <- NA_character_
  }
  as_sir_structure(out)
}
