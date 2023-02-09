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

#' Calculate Microbial Resistance
#'
#' @description These functions can be used to calculate the (co-)resistance or susceptibility of microbial isolates (i.e. percentage of S, SI, I, IR or R). All functions support quasiquotation with pipes, can be used in `summarise()` from the `dplyr` package and also support grouped variables, see *Examples*.
#'
#' [resistance()] should be used to calculate resistance, [susceptibility()] should be used to calculate susceptibility.\cr
#' @param ... one or more vectors (or columns) with antibiotic interpretations. They will be transformed internally with [as.sir()] if needed. Use multiple columns to calculate (the lack of) co-resistance: the probability where one of two drugs have a resistant or susceptible result. See *Examples*.
#' @param minimum the minimum allowed number of available (tested) isolates. Any isolate count lower than `minimum` will return `NA` with a warning. The default number of `30` isolates is advised by the Clinical and Laboratory Standards Institute (CLSI) as best practice, see *Source*.
#' @param as_percent a [logical] to indicate whether the output must be returned as a hundred fold with % sign (a character). A value of `0.123456` will then be returned as `"12.3%"`.
#' @param only_all_tested (for combination therapies, i.e. using more than one variable for `...`): a [logical] to indicate that isolates must be tested for all antibiotics, see section *Combination Therapy* below
#' @param data a [data.frame] containing columns with class [`sir`] (see [as.sir()])
#' @param translate_ab a column name of the [antibiotics] data set to translate the antibiotic abbreviations to, using [ab_property()]
#' @inheritParams ab_property
#' @param combine_SI a [logical] to indicate whether all values of S and I must be merged into one, so the output only consists of S+I vs. R (susceptible vs. resistant), defaults to `TRUE`
#' @param ab_result antibiotic results to test against, must be one or more values of "S", "I", or "R"
#' @param confidence_level the confidence level for the returned confidence interval. For the calculation, the number of S or SI isolates, and R isolates are compared with the total number of available isolates with R, S, or I by using [binom.test()], i.e., the Clopper-Pearson method.
#' @param side the side of the confidence interval to return. Defaults to `"both"` for a length 2 vector, but can also be (abbreviated as) `"min"`/`"left"`/`"lower"`/`"less"` or `"max"`/`"right"`/`"higher"`/`"greater"`.
#' @inheritSection as.sir Interpretation of SIR
#' @details
#' The function [resistance()] is equal to the function [proportion_R()]. The function [susceptibility()] is equal to the function [proportion_SI()].
#'
#' Use [sir_confidence_interval()] to calculate the confidence interval, which relies on [binom.test()], i.e., the Clopper-Pearson method. This function returns a vector of length 2 at default for antimicrobial *resistance*. Change the `side` argument to "left"/"min" or "right"/"max" to return a single value, and change the `ab_result` argument to e.g. `c("S", "I")` to test for antimicrobial *susceptibility*, see Examples.
#'
#' **Remember that you should filter your data to let it contain only first isolates!** This is needed to exclude duplicates and to reduce selection bias. Use [first_isolate()] to determine them in your data set.
#'
#' These functions are not meant to count isolates, but to calculate the proportion of resistance/susceptibility. Use the [`count()`][AMR::count()] functions to count isolates. The function [susceptibility()] is essentially equal to `count_susceptible() / count_all()`. *Low counts can influence the outcome - the `proportion` functions may camouflage this, since they only return the proportion (albeit being dependent on the `minimum` argument).*
#'
#' The function [proportion_df()] takes any variable from `data` that has an [`sir`] class (created with [as.sir()]) and calculates the proportions S, I, and R. It also supports grouped variables. The function [sir_df()] works exactly like [proportion_df()], but adds the number of isolates.
#' @section Combination Therapy:
#' When using more than one variable for `...` (= combination therapy), use `only_all_tested` to only count isolates that are tested for all antibiotics/variables that you test them for. See this example for two antibiotics, Drug A and Drug B, about how [susceptibility()] works to calculate the %SI:
#'
#' ```
#' --------------------------------------------------------------------
#'                     only_all_tested = FALSE  only_all_tested = TRUE
#'                     -----------------------  -----------------------
#'  Drug A    Drug B   include as  include as   include as  include as
#'                     numerator   denominator  numerator   denominator
#' --------  --------  ----------  -----------  ----------  -----------
#'  S or I    S or I       X            X            X            X
#'    R       S or I       X            X            X            X
#'   <NA>     S or I       X            X            -            -
#'  S or I      R          X            X            X            X
#'    R         R          -            X            -            X
#'   <NA>       R          -            -            -            -
#'  S or I     <NA>        X            X            -            -
#'    R        <NA>        -            -            -            -
#'   <NA>      <NA>        -            -            -            -
#' --------------------------------------------------------------------
#' ```
#'
#' Please note that, in combination therapies, for `only_all_tested = TRUE` applies that:
#' ```
#'     count_S()    +   count_I()    +   count_R()    = count_all()
#'   proportion_S() + proportion_I() + proportion_R() = 1
#' ```
#' and that, in combination therapies, for `only_all_tested = FALSE` applies that:
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
#'
#' # base R ------------------------------------------------------------
#' # determines %R
#' resistance(example_isolates$AMX)
#' sir_confidence_interval(example_isolates$AMX)
#' sir_confidence_interval(example_isolates$AMX,
#'   confidence_level = 0.975
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
#'   # scoped dplyr verbs with antibiotic selectors
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
      ab_result = c("S", "I"),
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
                                    side = "both") {
  meet_criteria(ab_result, allow_class = c("character", "sir"), has_length = c(1, 2, 3), is_in = c("S", "I", "R"))
  meet_criteria(confidence_level, allow_class = "numeric", is_positive = TRUE, has_length = 1)
  meet_criteria(side, allow_class = "character", has_length = 1, is_in = c("both", "b", "left", "l", "lower", "lowest", "less", "min", "right", "r", "higher", "highest", "greater", "g", "max"))
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
      ab_result = c("S", "I", "R"),
      only_all_tested = only_all_tested,
      only_count = TRUE
    ),
    error = function(e) stop_(gsub("in sir_calc(): ", "", e$message, fixed = TRUE), call = -5)
  )

  if (n < minimum) {
    warning_("Introducing NA: ",
      ifelse(n == 0, "no", paste("only", n)),
      " results available for `sir_confidence_interval()` (`minimum` = ", minimum, ").",
      call = FALSE
    )
    if (as_percent == TRUE) {
      return(NA_character_)
    } else {
      return(NA_real_)
    }
  }

  out <- stats::binom.test(x = x, n = n, conf.level = confidence_level)$conf.int
  out <- set_clean_class(out, "double")

  if (side %in% c("left", "l", "lower", "lowest", "less", "min")) {
    out <- out[1]
  } else if (side %in% c("right", "r", "higher", "highest", "greater", "g", "max")) {
    out <- out[2]
  }
  if (as_percent == TRUE) {
    percentage(out, digits = 1)
  } else {
    out
  }
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
  tryCatch(
    sir_calc(...,
      ab_result = c("I", "R"),
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
  tryCatch(
    sir_calc(...,
      ab_result = "I",
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
  tryCatch(
    sir_calc(...,
      ab_result = c("S", "I"),
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
