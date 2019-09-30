# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Analysis                              #
#                                                                      #
# SOURCE                                                               #
# https://gitlab.com/msberends/AMR                                     #
#                                                                      #
# LICENCE                                                              #
# (c) 2019 Berends MS (m.s.berends@umcg.nl), Luz CF (c.f.luz@umcg.nl)  #
#                                                                      #
# This R package is free software; you can freely use and distribute   #
# it for both personal and commercial purposes under the terms of the  #
# GNU General Public License version 2.0 (GNU GPL-2), as published by  #
# the Free Software Foundation.                                        #
#                                                                      #
# This R package was created for academic research and was publicly    #
# released in the hope that it will be useful, but it comes WITHOUT    #
# ANY WARRANTY OR LIABILITY.                                           #
# Visit our website for more info: https://msberends.gitlab.io/AMR.    #
# ==================================================================== #

#' Calculate resistance of isolates
#'
#' @description These functions can be used to calculate the (co-)resistance of microbial isolates (i.e. percentage of S, SI, I, IR or R). All functions support quasiquotation with pipes, can be used in \code{dplyr}s \code{\link[dplyr]{summarise}} and support grouped variables, see \emph{Examples}.
#'
#' \code{portion_R} and \code{portion_IR} can be used to calculate resistance, \code{portion_S} and \code{portion_SI} can be used to calculate susceptibility.\cr
#' @param ... one or more vectors (or columns) with antibiotic interpretations. They will be transformed internally with \code{\link{as.rsi}} if needed. Use multiple columns to calculate (the lack of) co-resistance: the probability where one of two drugs have a resistant or susceptible result. See Examples.
#' @param minimum the minimum allowed number of available (tested) isolates. Any isolate count lower than \code{minimum} will return \code{NA} with a warning. The default number of \code{30} isolates is advised by the Clinical and Laboratory Standards Institute (CLSI) as best practice, see Source.
#' @param as_percent a logical to indicate whether the output must be returned as a hundred fold with \% sign (a character) using\code{\link[clean]{percentage}}. A value of \code{0.123456} will then be returned as \code{"12.3\%"}.
#' @param only_all_tested (for combination therapies, i.e. using more than one variable for \code{...}) a logical to indicate that isolates must be tested for all antibiotics, see section \emph{Combination therapy} below
#' @param data a \code{data.frame} containing columns with class \code{rsi} (see \code{\link{as.rsi}})
#' @param translate_ab a column name of the \code{\link{antibiotics}} data set to translate the antibiotic abbreviations to, using \code{\link{ab_property}}
#' @inheritParams ab_property
#' @param combine_SI a logical to indicate whether all values of S and I must be merged into one, so the output only consists of S+I vs. R (susceptible vs. resistant). This used to be the parameter \code{combine_IR}, but this now follows the redefinition by EUCAST about the interpretion of I (increased exposure) in 2019, see section 'Interpretation of S, I and R' below. Default is \code{TRUE}.
#' @param combine_IR a logical to indicate whether all values of I and R must be merged into one, so the output only consists of S vs. I+R (susceptible vs. non-susceptible). This is outdated, see parameter \code{combine_SI}.
#' @inheritSection as.rsi Interpretation of S, I and R
#' @details \strong{Remember that you should filter your table to let it contain only first isolates!} This is needed to exclude duplicates and to reduce selection bias. Use \code{\link{first_isolate}} to determine them in your data set.
#'
#' These functions are not meant to count isolates, but to calculate the portion of resistance/susceptibility. Use the \code{\link[AMR]{count}} functions to count isolates. The function \code{portion_SI()} is essentially equal to \code{count_SI() / count_all()}. \emph{Low counts can infuence the outcome - the \code{portion} functions may camouflage this, since they only return the portion (albeit being dependent on the \code{minimum} parameter).}
#'
#' The function \code{portion_df} takes any variable from \code{data} that has an \code{"rsi"} class (created with \code{\link{as.rsi}}) and calculates the portions R, I and S. The resulting \emph{tidy data} (see Source) \code{data.frame} will have three rows (S/I/R) and a column for each group and each variable with class \code{"rsi"}.
#'
#' The function \code{rsi_df} works exactly like \code{portion_df}, but adds the number of isolates.
#' @section Combination therapy:
#' When using more than one variable for \code{...} (= combination therapy)), use \code{only_all_tested} to only count isolates that are tested for all antibiotics/variables that you test them for. See this example for two antibiotics, Antibiotic A and Antibiotic B, about how \code{portion_SI} works to calculate the \%SI:
#'
#' \preformatted{
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
#' }
#'
#' Please note that, in combination therapies, for \code{only_all_tested = TRUE} applies that:
#' \preformatted{
#'    count_S()  +  count_I()  +  count_R()  == count_all()
#'   portion_S() + portion_I() + portion_R() == 1
#' }
#' and that, in combination therapies, for \code{only_all_tested = FALSE} applies that:
#' \preformatted{
#'    count_S()  +  count_I()  +  count_R()  >= count_all()
#'   portion_S() + portion_I() + portion_R() >= 1
#' }
#'
#' Using \code{only_all_tested} has no impact when only using one antibiotic as input.
#' @source \strong{M39 Analysis and Presentation of Cumulative Antimicrobial Susceptibility Test Data, 4th Edition}, 2014, \emph{Clinical and Laboratory Standards Institute (CLSI)}. \url{https://clsi.org/standards/products/microbiology/documents/m39/}.
#'
#' Wickham H. \strong{Tidy Data.} The Journal of Statistical Software, vol. 59, 2014. \url{http://vita.had.co.nz/papers/tidy-data.html}
#' @seealso \code{\link[AMR]{count}_*} to count resistant and susceptible isolates.
#' @keywords resistance susceptibility rsi_df rsi antibiotics isolate isolates
#' @return Double or, when \code{as_percent = TRUE}, a character.
#' @rdname portion
#' @name portion
#' @export
#' @inheritSection AMR Read more on our website!
#' @examples
#' # example_isolates is a data set available in the AMR package.
#' ?example_isolates
#'
#' # Calculate resistance
#' portion_R(example_isolates$AMX)
#' portion_IR(example_isolates$AMX)
#'
#' # Or susceptibility
#' portion_S(example_isolates$AMX)
#' portion_SI(example_isolates$AMX)
#'

#' # Do the above with pipes:
#' library(dplyr)
#' example_isolates %>% portion_R(AMX)
#' example_isolates %>% portion_IR(AMX)
#' example_isolates %>% portion_S(AMX)
#' example_isolates %>% portion_SI(AMX)
#'
#' example_isolates %>%
#'   group_by(hospital_id) %>%
#'   summarise(p = portion_SI(CIP),
#'             n = n_rsi(CIP)) # n_rsi works like n_distinct in dplyr
#'
#' example_isolates %>%
#'   group_by(hospital_id) %>%
#'   summarise(R = portion_R(CIP, as_percent = TRUE),
#'             I = portion_I(CIP, as_percent = TRUE),
#'             S = portion_S(CIP, as_percent = TRUE),
#'             n1 = count_all(CIP),  # the actual total; sum of all three
#'             n2 = n_rsi(CIP),      # same - analogous to n_distinct
#'             total = n())          # NOT the number of tested isolates!
#'
#' # Calculate co-resistance between amoxicillin/clav acid and gentamicin,
#' # so we can see that combination therapy does a lot more than mono therapy:
#' example_isolates %>% portion_SI(AMC)      # %SI = 76.3%
#' example_isolates %>% count_all(AMC)       #   n = 1879
#'
#' example_isolates %>% portion_SI(GEN)      # %SI = 75.4%
#' example_isolates %>% count_all(GEN)       #   n = 1855
#'
#' example_isolates %>% portion_SI(AMC, GEN) # %SI = 94.1%
#' example_isolates %>% count_all(AMC, GEN)  #   n = 1939
#'
#'
#' # See Details on how `only_all_tested` works. Example:
#' example_isolates %>%
#'   summarise(numerator = count_SI(AMC, GEN),
#'             denominator = count_all(AMC, GEN),
#'             portion = portion_SI(AMC, GEN))
#' #   numerator denominator portion
#' #        1764        1936  0.9408
#' example_isolates %>%
#'   summarise(numerator = count_SI(AMC, GEN, only_all_tested = TRUE),
#'             denominator = count_all(AMC, GEN, only_all_tested = TRUE),
#'             portion = portion_SI(AMC, GEN, only_all_tested = TRUE))
#' #   numerator denominator portion
#' #       1687        1798   0.9383
#'
#'
#' example_isolates %>%
#'   group_by(hospital_id) %>%
#'   summarise(cipro_p = portion_SI(CIP, as_percent = TRUE),
#'             cipro_n = count_all(CIP),
#'             genta_p = portion_SI(GEN, as_percent = TRUE),
#'             genta_n = count_all(GEN),
#'             combination_p = portion_SI(CIP, GEN, as_percent = TRUE),
#'             combination_n = count_all(CIP, GEN))
#'
#' # Get portions S/I/R immediately of all rsi columns
#' example_isolates %>%
#'   select(AMX, CIP) %>%
#'   portion_df(translate = FALSE)
#'
#' # It also supports grouping variables
#' example_isolates %>%
#'   select(hospital_id, AMX, CIP) %>%
#'   group_by(hospital_id) %>%
#'   portion_df(translate = FALSE)
#'
#'
#' \dontrun{
#'
#' # calculate current empiric combination therapy of Helicobacter gastritis:
#' my_table %>%
#'   filter(first_isolate == TRUE,
#'          genus == "Helicobacter") %>%
#'   summarise(p = portion_S(AMX, MTR),  # amoxicillin with metronidazole
#'             n = count_all(AMX, MTR))
#' }
portion_R <- function(...,
                      minimum = 30,
                      as_percent = FALSE,
                      only_all_tested = FALSE) {
  rsi_calc(...,
           ab_result = "R",
           minimum = minimum,
           as_percent = as_percent,
           only_all_tested = only_all_tested,
           only_count = FALSE)
}

#' @rdname portion
#' @export
portion_IR <- function(...,
                       minimum = 30,
                       as_percent = FALSE,
                       only_all_tested = FALSE) {
  rsi_calc(...,
           ab_result = c("I", "R"),
           minimum = minimum,
           as_percent = as_percent,
           only_all_tested = only_all_tested,
           only_count = FALSE)
}

#' @rdname portion
#' @export
portion_I <- function(...,
                      minimum = 30,
                      as_percent = FALSE,
                      only_all_tested = FALSE) {
  rsi_calc(...,
           ab_result = "I",
           minimum = minimum,
           as_percent = as_percent,
           only_all_tested = only_all_tested,
           only_count = FALSE)
}

#' @rdname portion
#' @export
portion_SI <- function(...,
                       minimum = 30,
                       as_percent = FALSE,
                       only_all_tested = FALSE) {
  rsi_calc(...,
           ab_result = c("S", "I"),
           minimum = minimum,
           as_percent = as_percent,
           only_all_tested = only_all_tested,
           only_count = FALSE)
}

#' @rdname portion
#' @export
portion_S <- function(...,
                      minimum = 30,
                      as_percent = FALSE,
                      only_all_tested = FALSE) {
  rsi_calc(...,
           ab_result = "S",
           minimum = minimum,
           as_percent = as_percent,
           only_all_tested = only_all_tested,
           only_count = FALSE)
}

#' @rdname portion
#' @importFrom dplyr %>% select_if bind_rows summarise_if mutate group_vars select everything
#' @export
portion_df <- function(data,
                       translate_ab = "name",
                       language = get_locale(),
                       minimum = 30,
                       as_percent = FALSE,
                       combine_SI = TRUE,
                       combine_IR = FALSE) {

  rsi_calc_df(type = "portion",
              data = data,
              translate_ab = translate_ab,
              language = language,
              minimum = minimum,
              as_percent = as_percent,
              combine_SI = combine_SI,
              combine_IR = combine_IR,
              combine_SI_missing = missing(combine_SI))
}
