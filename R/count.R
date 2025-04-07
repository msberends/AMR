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

#' Count Available Isolates
#'
#' @description These functions can be used to count resistant/susceptible microbial isolates. All functions support quasiquotation with pipes, can be used in `summarise()` from the `dplyr` package and also support grouped variables, see *Examples*.
#'
#' [count_resistant()] should be used to count resistant isolates, [count_susceptible()] should be used to count susceptible isolates.
#' @param ... One or more vectors (or columns) with antibiotic interpretations. They will be transformed internally with [as.sir()] if needed.
#' @inheritParams proportion
#' @inheritSection as.sir Interpretation of SIR
#' @details These functions are meant to count isolates. Use the [resistance()]/[susceptibility()] functions to calculate microbial resistance/susceptibility.
#'
#' The function [count_resistant()] is equal to the function [count_R()]. The function [count_susceptible()] is equal to the function [count_SI()].
#'
#' The function [n_sir()] is an alias of [count_all()]. They can be used to count all available isolates, i.e. where all input antimicrobials have an available result (S, I or R). Their use is equal to `n_distinct()`. Their function is equal to `count_susceptible(...) + count_resistant(...)`.
#'
#' The function [count_df()] takes any variable from `data` that has an [`sir`] class (created with [as.sir()]) and counts the number of S's, I's and R's. It also supports grouped variables. The function [sir_df()] works exactly like [count_df()], but adds the percentage of S, I and R.
#' @inheritSection proportion Combination Therapy
#' @seealso [`proportion_*`][proportion] to calculate microbial resistance and susceptibility.
#' @return An [integer]
#' @rdname count
#' @name count
#' @export
#' @examples
#' # example_isolates is a data set available in the AMR package.
#' # run ?example_isolates for more info.
#'
#' # base R ------------------------------------------------------------
#' count_resistant(example_isolates$AMX) # counts "R"
#' count_susceptible(example_isolates$AMX) # counts "S" and "I"
#' count_all(example_isolates$AMX) # counts "S", "I" and "R"
#'
#' # be more specific
#' count_S(example_isolates$AMX)
#' count_SI(example_isolates$AMX)
#' count_I(example_isolates$AMX)
#' count_IR(example_isolates$AMX)
#' count_R(example_isolates$AMX)
#'
#' # Count all available isolates
#' count_all(example_isolates$AMX)
#' n_sir(example_isolates$AMX)
#'
#' # n_sir() is an alias of count_all().
#' # Since it counts all available isolates, you can
#' # calculate back to count e.g. susceptible isolates.
#' # These results are the same:
#' count_susceptible(example_isolates$AMX)
#' susceptibility(example_isolates$AMX) * n_sir(example_isolates$AMX)
#'
#' # dplyr -------------------------------------------------------------
#' \donttest{
#' if (require("dplyr")) {
#'   example_isolates %>%
#'     group_by(ward) %>%
#'     summarise(
#'       R = count_R(CIP),
#'       I = count_I(CIP),
#'       S = count_S(CIP),
#'       n1 = count_all(CIP), # the actual total; sum of all three
#'       n2 = n_sir(CIP), # same - analogous to n_distinct
#'       total = n()
#'     ) # NOT the number of tested isolates!
#'
#'   # Number of available isolates for a whole antibiotic class
#'   # (i.e., in this data set columns GEN, TOB, AMK, KAN)
#'   example_isolates %>%
#'     group_by(ward) %>%
#'     summarise(across(aminoglycosides(), n_sir))
#'
#'   # Count co-resistance between amoxicillin/clav acid and gentamicin,
#'   # so we can see that combination therapy does a lot more than mono therapy.
#'   # Please mind that `susceptibility()` calculates percentages right away instead.
#'   example_isolates %>% count_susceptible(AMC) # 1433
#'   example_isolates %>% count_all(AMC) # 1879
#'
#'   example_isolates %>% count_susceptible(GEN) # 1399
#'   example_isolates %>% count_all(GEN) # 1855
#'
#'   example_isolates %>% count_susceptible(AMC, GEN) # 1764
#'   example_isolates %>% count_all(AMC, GEN) # 1936
#'
#'   # Get number of S+I vs. R immediately of selected columns
#'   example_isolates %>%
#'     select(AMX, CIP) %>%
#'     count_df(translate = FALSE)
#'
#'   # It also supports grouping variables
#'   example_isolates %>%
#'     select(ward, AMX, CIP) %>%
#'     group_by(ward) %>%
#'     count_df(translate = FALSE)
#' }
#' }
count_resistant <- function(..., only_all_tested = FALSE) {
  tryCatch(
    sir_calc(...,
      ab_result = "R",
      only_all_tested = only_all_tested,
      only_count = TRUE
    ),
    error = function(e) stop_(gsub("in sir_calc(): ", "", e$message, fixed = TRUE), call = -5)
  )
}

#' @rdname count
#' @export
count_susceptible <- function(..., only_all_tested = FALSE) {
  tryCatch(
    sir_calc(...,
      ab_result = c("S", "SDD", "I"),
      only_all_tested = only_all_tested,
      only_count = TRUE
    ),
    error = function(e) stop_(gsub("in sir_calc(): ", "", e$message, fixed = TRUE), call = -5)
  )
}

#' @rdname count
#' @export
count_S <- function(..., only_all_tested = FALSE) {
  tryCatch(
    sir_calc(...,
      ab_result = "S",
      only_all_tested = only_all_tested,
      only_count = TRUE
    ),
    error = function(e) stop_(gsub("in sir_calc(): ", "", e$message, fixed = TRUE), call = -5)
  )
}

#' @rdname count
#' @export
count_SI <- function(..., only_all_tested = FALSE) {
  tryCatch(
    sir_calc(...,
      ab_result = c("S", "SDD", "I"),
      only_all_tested = only_all_tested,
      only_count = TRUE
    ),
    error = function(e) stop_(gsub("in sir_calc(): ", "", e$message, fixed = TRUE), call = -5)
  )
}

#' @rdname count
#' @export
count_I <- function(..., only_all_tested = FALSE) {
  tryCatch(
    sir_calc(...,
      ab_result = c("I", "SDD"),
      only_all_tested = only_all_tested,
      only_count = TRUE
    ),
    error = function(e) stop_(gsub("in sir_calc(): ", "", e$message, fixed = TRUE), call = -5)
  )
}

#' @rdname count
#' @export
count_IR <- function(..., only_all_tested = FALSE) {
  tryCatch(
    sir_calc(...,
      ab_result = c("I", "SDD", "R"),
      only_all_tested = only_all_tested,
      only_count = TRUE
    ),
    error = function(e) stop_(gsub("in sir_calc(): ", "", e$message, fixed = TRUE), call = -5)
  )
}

#' @rdname count
#' @export
count_R <- function(..., only_all_tested = FALSE) {
  tryCatch(
    sir_calc(...,
      ab_result = "R",
      only_all_tested = only_all_tested,
      only_count = TRUE
    ),
    error = function(e) stop_(gsub("in sir_calc(): ", "", e$message, fixed = TRUE), call = -5)
  )
}

#' @rdname count
#' @export
count_all <- function(..., only_all_tested = FALSE) {
  tryCatch(
    sir_calc(...,
      ab_result = c("S", "SDD", "I", "R", "NI"),
      only_all_tested = only_all_tested,
      only_count = TRUE
    ),
    error = function(e) stop_(gsub("in sir_calc(): ", "", e$message, fixed = TRUE), call = -5)
  )
}

#' @rdname count
#' @export
n_sir <- count_all

#' @rdname count
#' @export
count_df <- function(data,
                     translate_ab = "name",
                     language = get_AMR_locale(),
                     combine_SI = TRUE) {
  tryCatch(
    sir_calc_df(
      type = "count",
      data = data,
      translate_ab = translate_ab,
      language = language,
      combine_SI = combine_SI,
      confidence_level = 0.95 # doesn't matter, will be removed
    ),
    error = function(e) stop_(gsub("in sir_calc_df(): ", "", e$message, fixed = TRUE), call = -5)
  )
}
