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
# how to conduct AMR data analysis: https://msberends.github.io/AMR/   #
# ==================================================================== #

#' Calculate Antimicrobial Resistance
#'
#' @description These functions can be used to calculate the (co-)resistance or susceptibility of microbial isolates (i.e. percentage of S, SI, I, IR or R). All functions support quasiquotation with pipes, can be used in `summarise()` from the `dplyr` package and also support grouped variables, see *Examples*.
#'
#' [resistance()] should be used to calculate resistance, [susceptibility()] should be used to calculate susceptibility.\cr
#' @param ... One or more vectors (or columns) with antibiotic interpretations. They will be transformed internally with [as.sir()] if needed. Use multiple columns to calculate (the lack of) co-resistance: the probability where one of two drugs have a resistant or susceptible result. See *Examples*.
#' @param minimum The minimum allowed number of available (tested) isolates. Any isolate count lower than `minimum` will return `NA` with a warning. The default number of `30` isolates is advised by the Clinical and Laboratory Standards Institute (CLSI) as best practice, see *Source*.
#' @param as_percent A [logical] to indicate whether the output must be returned as a hundred fold with % sign (a character). A value of `0.123456` will then be returned as `"12.3%"`.
#' @param only_all_tested (for combination therapies, i.e. using more than one variable for `...`): a [logical] to indicate that isolates must be tested for all antimicrobials, see section *Combination Therapy* below
#' @param data A [data.frame] containing columns with class [`sir`] (see [as.sir()])
#' @param translate_ab A column name of the [antimicrobials] data set to translate the antibiotic abbreviations to, using [ab_property()]
#' @inheritParams ab_property
#' @param combine_SI A [logical] to indicate whether all values of S, SDD, and I must be merged into one, so the output only consists of S+SDD+I vs. R (susceptible vs. resistant) - the default is `TRUE`
#' @param ab_result Antibiotic results to test against, must be one or more values of "S", "SDD", "I", or "R"
#' @param confidence_level The confidence level for the returned confidence interval. For the calculation, the number of S or SI isolates, and R isolates are compared with the total number of available isolates with R, S, or I by using [binom.test()], i.e., the Clopper-Pearson method.
#' @param side The side of the confidence interval to return. The default is `"both"` for a length 2 vector, but can also be (abbreviated as) `"min"`/`"left"`/`"lower"`/`"less"` or `"max"`/`"right"`/`"higher"`/`"greater"`.
#' @param collapse A [logical] to indicate whether the output values should be 'collapsed', i.e. be merged together into one value, or a character value to use for collapsing
#' @inheritSection as.sir Interpretation of SIR
#' @details
#' For a more automated and comprehensive analysis, consider using [antibiogram()] or [wisca()], which streamline many aspects of susceptibility reporting and, importantly, also support WISCA. The functions described here offer a more hands-on, manual approach for greater customisation.
#'
#' **Remember that you should filter your data to let it contain only first isolates!** This is needed to exclude duplicates and to reduce selection bias. Use [first_isolate()] to determine them in your data set with one of the four available algorithms.
#'
#' The function [resistance()] is equal to the function [proportion_R()]. The function [susceptibility()] is equal to the function [proportion_SI()]. Since AMR v3.0, [proportion_SI()] and [proportion_I()] include dose-dependent susceptibility ('SDD').
#'
#' Use [sir_confidence_interval()] to calculate the confidence interval, which relies on [binom.test()], i.e., the Clopper-Pearson method. This function returns a vector of length 2 at default for antimicrobial *resistance*. Change the `side` argument to "left"/"min" or "right"/"max" to return a single value, and change the `ab_result` argument to e.g. `c("S", "I")` to test for antimicrobial *susceptibility*, see Examples.
#'
#' These functions are not meant to count isolates, but to calculate the proportion of resistance/susceptibility. Use the [`count_*()`][AMR::count()] functions to count isolates. The function [susceptibility()] is essentially equal to [count_susceptible()]` / `[count_all()]. *Low counts can influence the outcome - the `proportion_*()` functions may camouflage this, since they only return the proportion (albeit dependent on the `minimum` argument).*
#'
#' The function [proportion_df()] takes any variable from `data` that has an [`sir`] class (created with [as.sir()]) and calculates the proportions S, I, and R. It also supports grouped variables. The function [sir_df()] works exactly like [proportion_df()], but adds the number of isolates.
#' @section Combination Therapy:
#' When using more than one variable for `...` (= combination therapy), use `only_all_tested` to only count isolates that are tested for all antimicrobials/variables that you test them for. See this example for two antimicrobials, Drug A and Drug B, about how [susceptibility()] works to calculate the %SI:
#'
#'
#' ```
#' --------------------------------------------------------------------
#'                     only_all_tested = FALSE  only_all_tested = TRUE
#'                     -----------------------  -----------------------
#'  Drug A    Drug B   considered   considered  considered   considered
#'                     susceptible    tested    susceptible    tested
#' --------  --------  -----------  ----------  -----------  ----------
#'  S or I    S or I        X            X           X            X
#'    R       S or I        X            X           X            X
#'   <NA>     S or I        X            X           -            -
#'  S or I      R           X            X           X            X
#'    R         R           -            X           -            X
#'   <NA>       R           -            -           -            -
#'  S or I     <NA>         X            X           -            -
#'    R        <NA>         -            -           -            -
#'   <NA>      <NA>         -            -           -            -
#' --------------------------------------------------------------------
#' ```
#'
#' Please note that, in combination therapies, for `only_all_tested = TRUE` applies that:
#'
#' ```
#'     count_S()    +   count_I()    +   count_R()    = count_all()
#'   proportion_S() + proportion_I() + proportion_R() = 1
#' ```
#'
#' and that, in combination therapies, for `only_all_tested = FALSE` applies that:
#'
#' ```
#'     count_S()    +   count_I()    +   count_R()    >= count_all()
#'   proportion_S() + proportion_I() + proportion_R() >= 1
#' ```
#'
#' Using `only_all_tested` has no impact when only using one antibiotic as input.
#' @source **M39 Analysis and Presentation of Cumulative Antimicrobial Susceptibility Test Data, 5th Edition**, 2022, *Clinical and Laboratory Standards Institute (CLSI)*. <https://clsi.org/standards/products/microbiology/documents/m39/>.
#' @seealso [AMR::count()] to count resistant and susceptible isolates.
#' @return A [double] or, when `as_percent = TRUE`, a [character].
#' @rdname proportion
#' @aliases portion
#' @name proportion
#' @export
#' @examples
#' # example_isolates is a data set available in the AMR package.
#' # run ?example_isolates for more info.
#' example_isolates
#'
#'
#' # base R ------------------------------------------------------------
#' # determines %R
#' resistance(example_isolates$AMX)
#' sir_confidence_interval(example_isolates$AMX)
#' sir_confidence_interval(example_isolates$AMX,
#'   confidence_level = 0.975
#' )
#' sir_confidence_interval(example_isolates$AMX,
#'   confidence_level = 0.975,
#'   collapse = ", "
#' )
#'
#' # determines %S+I:
#' susceptibility(example_isolates$AMX)
#' sir_confidence_interval(example_isolates$AMX,
#'   ab_result = c("S", "I")
#' )
#'
#' # be more specific
#' proportion_S(example_isolates$AMX)
#' proportion_SI(example_isolates$AMX)
#' proportion_I(example_isolates$AMX)
#' proportion_IR(example_isolates$AMX)
#' proportion_R(example_isolates$AMX)
#'
#' # dplyr -------------------------------------------------------------
#' \donttest{
#' if (require("dplyr")) {
#'   example_isolates %>%
#'     group_by(ward) %>%
#'     summarise(
#'       r = resistance(CIP),
#'       n = n_sir(CIP)
#'     ) # n_sir works like n_distinct in dplyr, see ?n_sir
#' }
#' if (require("dplyr")) {
#'   example_isolates %>%
#'     group_by(ward) %>%
#'     summarise(
#'       cipro_R = resistance(CIP),
#'       ci_min = sir_confidence_interval(CIP, side = "min"),
#'       ci_max = sir_confidence_interval(CIP, side = "max"),
#'     )
#' }
#' if (require("dplyr")) {
#'   # scoped dplyr verbs with antimicrobial selectors
#'   # (you could also use across() of course)
#'   example_isolates %>%
#'     group_by(ward) %>%
#'     summarise_at(
#'       c(aminoglycosides(), carbapenems()),
#'       resistance
#'     )
#' }
#' if (require("dplyr")) {
#'   example_isolates %>%
#'     group_by(ward) %>%
#'     summarise(
#'       R = resistance(CIP, as_percent = TRUE),
#'       SI = susceptibility(CIP, as_percent = TRUE),
#'       n1 = count_all(CIP), # the actual total; sum of all three
#'       n2 = n_sir(CIP), # same - analogous to n_distinct
#'       total = n()
#'     ) # NOT the number of tested isolates!
#'
#'   # Calculate co-resistance between amoxicillin/clav acid and gentamicin,
#'   # so we can see that combination therapy does a lot more than mono therapy:
#'   example_isolates %>% susceptibility(AMC) # %SI = 76.3%
#'   example_isolates %>% count_all(AMC) #   n = 1879
#'
#'   example_isolates %>% susceptibility(GEN) # %SI = 75.4%
#'   example_isolates %>% count_all(GEN) #   n = 1855
#'
#'   example_isolates %>% susceptibility(AMC, GEN) # %SI = 94.1%
#'   example_isolates %>% count_all(AMC, GEN) #   n = 1939
#'
#'
#'   # See Details on how `only_all_tested` works. Example:
#'   example_isolates %>%
#'     summarise(
#'       numerator = count_susceptible(AMC, GEN),
#'       denominator = count_all(AMC, GEN),
#'       proportion = susceptibility(AMC, GEN)
#'     )
#'
#'   example_isolates %>%
#'     summarise(
#'       numerator = count_susceptible(AMC, GEN, only_all_tested = TRUE),
#'       denominator = count_all(AMC, GEN, only_all_tested = TRUE),
#'       proportion = susceptibility(AMC, GEN, only_all_tested = TRUE)
#'     )
#'
#'
#'   example_isolates %>%
#'     group_by(ward) %>%
#'     summarise(
#'       cipro_p = susceptibility(CIP, as_percent = TRUE),
#'       cipro_n = count_all(CIP),
#'       genta_p = susceptibility(GEN, as_percent = TRUE),
#'       genta_n = count_all(GEN),
#'       combination_p = susceptibility(CIP, GEN, as_percent = TRUE),
#'       combination_n = count_all(CIP, GEN)
#'     )
#'
#'   # Get proportions S/I/R immediately of all sir columns
#'   example_isolates %>%
#'     select(AMX, CIP) %>%
#'     proportion_df(translate = FALSE)
#'
#'   # It also supports grouping variables
#'   # (use sir_df to also include the count)
#'   example_isolates %>%
#'     select(ward, AMX, CIP) %>%
#'     group_by(ward) %>%
#'     sir_df(translate = FALSE)
#' }
#' }
resistance <- function(...,
                       minimum = 30,
                       as_percent = FALSE,
                       only_all_tested = FALSE) {
  tryCatch(
    sir_calc(...,
      ab_result = "R",
      minimum = minimum,
      as_percent = as_percent,
      only_all_tested = only_all_tested,
      only_count = FALSE
    ),
    error = function(e) stop_(gsub("in sir_calc(): ", "", e$message, fixed = TRUE), call = -5)
  )
}

#' @rdname proportion
#' @export
susceptibility <- function(...,
                           minimum = 30,
                           as_percent = FALSE,
                           only_all_tested = FALSE) {
  tryCatch(
    sir_calc(...,
      ab_result = c("S", "SDD", "I"),
      minimum = minimum,
      as_percent = as_percent,
      only_all_tested = only_all_tested,
      only_count = FALSE
    ),
    error = function(e) stop_(gsub("in sir_calc(): ", "", e$message, fixed = TRUE), call = -5)
  )
}

#' @rdname proportion
#' @export
sir_confidence_interval <- function(...,
                                    ab_result = "R",
                                    minimum = 30,
                                    as_percent = FALSE,
                                    only_all_tested = FALSE,
                                    confidence_level = 0.95,
                                    side = "both",
                                    collapse = FALSE) {
  meet_criteria(ab_result, allow_class = c("character", "sir"), has_length = c(1:5), is_in = c("S", "SDD", "I", "R", "NI"))
  meet_criteria(minimum, allow_class = c("numeric", "integer"), has_length = 1, is_positive_or_zero = TRUE, is_finite = TRUE)
  meet_criteria(as_percent, allow_class = "logical", has_length = 1)
  meet_criteria(only_all_tested, allow_class = "logical", has_length = 1)
  meet_criteria(confidence_level, allow_class = "numeric", is_positive = TRUE, has_length = 1)
  meet_criteria(side, allow_class = "character", has_length = 1, is_in = c("both", "b", "left", "l", "lower", "lowest", "less", "min", "right", "r", "higher", "highest", "greater", "g", "max"))
  meet_criteria(collapse, allow_class = c("logical", "character"), has_length = 1)

  x <- tryCatch(
    sir_calc(...,
      ab_result = ab_result,
      only_all_tested = only_all_tested,
      only_count = TRUE
    ),
    error = function(e) stop_(gsub("in sir_calc(): ", "", e$message, fixed = TRUE), call = -5)
  )
  n <- tryCatch(
    sir_calc(...,
      ab_result = c("S", "SDD", "I", "R", "NI"),
      only_all_tested = only_all_tested,
      only_count = TRUE
    ),
    error = function(e) stop_(gsub("in sir_calc(): ", "", e$message, fixed = TRUE), call = -5)
  )

  if (x == 0) {
    out <- c(0, 0)
  } else {
    # this applies the Clopper-Pearson method
    out <- stats::binom.test(x = x, n = n, conf.level = confidence_level)$conf.int
  }
  out <- set_clean_class(out, "numeric")

  if (side %in% c("left", "l", "lower", "lowest", "less", "min")) {
    out <- out[1]
  } else if (side %in% c("right", "r", "higher", "highest", "greater", "g", "max")) {
    out <- out[2]
  }
  if (isTRUE(as_percent)) {
    out <- trimws(percentage(out, digits = 1))
  }
  if (!isFALSE(collapse) && length(out) > 1) {
    if (is.numeric(out)) {
      out <- round(out, digits = 3)
    }
    # out[is.na(out)] <- 0
    out <- paste(out, collapse = ifelse(isTRUE(collapse), "-", collapse))
  }

  if (n < minimum) {
    warning_("Introducing NA: ",
      ifelse(n == 0, "no", paste("only", n)),
      " results available for `sir_confidence_interval()` (`minimum` = ", minimum, ").",
      call = FALSE
    )
    if (is.character(out)) {
      return(NA_character_)
    } else {
      return(NA_real_)
    }
  }
  out
}

#' @rdname proportion
#' @export
proportion_R <- function(...,
                         minimum = 30,
                         as_percent = FALSE,
                         only_all_tested = FALSE) {
  tryCatch(
    sir_calc(...,
      ab_result = "R",
      minimum = minimum,
      as_percent = as_percent,
      only_all_tested = only_all_tested,
      only_count = FALSE
    ),
    error = function(e) stop_(gsub("in sir_calc(): ", "", e$message, fixed = TRUE), call = -5)
  )
}

#' @rdname proportion
#' @export
proportion_IR <- function(...,
                          minimum = 30,
                          as_percent = FALSE,
                          only_all_tested = FALSE) {
  if (message_not_thrown_before("proportion_IR", entire_session = TRUE)) {
    message_("Note that `proportion_IR()` will also include dose-dependent susceptibility, 'SDD'. This note will be shown once for this session.", as_note = FALSE)
  }
  tryCatch(
    sir_calc(...,
      ab_result = c("I", "SDD", "R"),
      minimum = minimum,
      as_percent = as_percent,
      only_all_tested = only_all_tested,
      only_count = FALSE
    ),
    error = function(e) stop_(gsub("in sir_calc(): ", "", e$message, fixed = TRUE), call = -5)
  )
}

#' @rdname proportion
#' @export
proportion_I <- function(...,
                         minimum = 30,
                         as_percent = FALSE,
                         only_all_tested = FALSE) {
  if (message_not_thrown_before("proportion_I", entire_session = TRUE)) {
    message_("Note that `proportion_I()` will also include dose-dependent susceptibility, 'SDD'. This note will be shown once for this session.", as_note = FALSE)
  }
  tryCatch(
    sir_calc(...,
      ab_result = c("I", "SDD"),
      minimum = minimum,
      as_percent = as_percent,
      only_all_tested = only_all_tested,
      only_count = FALSE
    ),
    error = function(e) stop_(gsub("in sir_calc(): ", "", e$message, fixed = TRUE), call = -5)
  )
}

#' @rdname proportion
#' @export
proportion_SI <- function(...,
                          minimum = 30,
                          as_percent = FALSE,
                          only_all_tested = FALSE) {
  if (message_not_thrown_before("proportion_SI", entire_session = TRUE)) {
    message_("Note that `proportion_SI()` will also include dose-dependent susceptibility, 'SDD'. This note will be shown once for this session.", as_note = FALSE)
  }
  tryCatch(
    sir_calc(...,
      ab_result = c("S", "I", "SDD"),
      minimum = minimum,
      as_percent = as_percent,
      only_all_tested = only_all_tested,
      only_count = FALSE
    ),
    error = function(e) stop_(gsub("in sir_calc(): ", "", e$message, fixed = TRUE), call = -5)
  )
}

#' @rdname proportion
#' @export
proportion_S <- function(...,
                         minimum = 30,
                         as_percent = FALSE,
                         only_all_tested = FALSE) {
  tryCatch(
    sir_calc(...,
      ab_result = "S",
      minimum = minimum,
      as_percent = as_percent,
      only_all_tested = only_all_tested,
      only_count = FALSE
    ),
    error = function(e) stop_(gsub("in sir_calc(): ", "", e$message, fixed = TRUE), call = -5)
  )
}

#' @rdname proportion
#' @export
proportion_df <- function(data,
                          translate_ab = "name",
                          language = get_AMR_locale(),
                          minimum = 30,
                          as_percent = FALSE,
                          combine_SI = TRUE,
                          confidence_level = 0.95) {
  tryCatch(
    sir_calc_df(
      type = "proportion",
      data = data,
      translate_ab = translate_ab,
      language = language,
      minimum = minimum,
      as_percent = as_percent,
      combine_SI = combine_SI,
      confidence_level = confidence_level
    ),
    error = function(e) stop_(gsub("in sir_calc_df(): ", "", e$message, fixed = TRUE), call = -5)
  )
}
