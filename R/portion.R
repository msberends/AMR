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
#' @param as_percent a logical to indicate whether the output must be returned as a hundred fold with \% sign (a character). A value of \code{0.123456} will then be returned as \code{"12.3\%"}.
#' @param also_single_tested a logical to indicate whether (in combination therapies) also observations should be included where not all antibiotics were tested, but at least one of the tested antibiotics contains a target interpretation (e.g. S in case of \code{portion_S} and R in case of \code{portion_R}). \strong{This would lead to selection bias in almost all cases.}
#' @param data a \code{data.frame} containing columns with class \code{rsi} (see \code{\link{as.rsi}})
#' @param translate_ab a column name of the \code{\link{antibiotics}} data set to translate the antibiotic abbreviations to, using \code{\link{ab_property}}
#' @inheritParams ab_property
#' @param combine_SI a logical to indicate whether all values of S and I must be merged into one, so the output only consists of S+I vs. R (susceptible vs. resistant). This used to be the parameter \code{combine_IR}, but this now follows the redefinition by EUCAST about the interpretion of I (increased exposure) in 2019, see section 'Interpretation of S, I and R' below. Default is \code{TRUE}.
#' @param combine_IR a logical to indicate whether all values of I and R must be merged into one, so the output only consists of S vs. I+R (susceptible vs. non-susceptible). This is outdated, see parameter \code{combine_SI}.
#' @inheritSection as.rsi Interpretation of S, I and R
#' @details \strong{Remember that you should filter your table to let it contain only first isolates!} Use \code{\link{first_isolate}} to determine them in your data set.
#'
#' These functions are not meant to count isolates, but to calculate the portion of resistance/susceptibility. Use the \code{\link[AMR]{count}} functions to count isolates. \emph{Low counts can infuence the outcome - these \code{portion} functions may camouflage this, since they only return the portion albeit being dependent on the \code{minimum} parameter.}
#'
#' \code{portion_df} takes any variable from \code{data} that has an \code{"rsi"} class (created with \code{\link{as.rsi}}) and calculates the portions R, I and S. The resulting \emph{tidy data} (see Source) \code{data.frame} will have three rows (S/I/R) and a column for each variable with class \code{"rsi"}.
#' \if{html}{
#    (created with https://www.latex4technics.com/)
#'   \cr\cr
#'   To calculate the probability (\emph{p}) of susceptibility of one antibiotic, we use this formula:
#'   \out{<div style="text-align: center;">}\figure{combi_therapy_2.png}\out{</div>}
#'   To calculate the probability (\emph{p}) of susceptibility of more antibiotics (i.e. combination therapy), we need to check whether one of them has a susceptible result (as numerator) and count all cases where all antibiotics were tested (as denominator). \cr
#'   \cr
#'   For two antibiotics:
#'   \out{<div style="text-align: center;">}\figure{combi_therapy_2.png}\out{</div>}
#'   \cr
#'   For three antibiotics:
#'   \out{<div style="text-align: center;">}\figure{combi_therapy_2.png}\out{</div>}
#'   \cr
#'   And so on.
#' }
#'
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
#' # septic_patients is a data set available in the AMR package. It is true, genuine data.
#' ?septic_patients
#'
#' # Calculate resistance
#' portion_R(septic_patients$AMX)
#' portion_IR(septic_patients$AMX)
#'
#' # Or susceptibility
#' portion_S(septic_patients$AMX)
#' portion_SI(septic_patients$AMX)
#'

#' # Do the above with pipes:
#' library(dplyr)
#' septic_patients %>% portion_R(AMX)
#' septic_patients %>% portion_IR(AMX)
#' septic_patients %>% portion_S(AMX)
#' septic_patients %>% portion_SI(AMX)
#'
#' septic_patients %>%
#'   group_by(hospital_id) %>%
#'   summarise(p = portion_S(CIP),
#'             n = n_rsi(CIP)) # n_rsi works like n_distinct in dplyr
#'
#' septic_patients %>%
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
#' septic_patients %>% portion_S(AMC)       # S = 71.4%
#' septic_patients %>% count_all(AMC)       # n = 1879
#'
#' septic_patients %>% portion_S(GEN)       # S = 74.0%
#' septic_patients %>% count_all(GEN)       # n = 1855
#'
#' septic_patients %>% portion_S(AMC, GEN)  # S = 92.3%
#' septic_patients %>% count_all(AMC, GEN)  # n = 1798
#'
#'
#' septic_patients %>%
#'   group_by(hospital_id) %>%
#'   summarise(cipro_p = portion_S(CIP, as_percent = TRUE),
#'             cipro_n = count_all(CIP),
#'             genta_p = portion_S(GEN, as_percent = TRUE),
#'             genta_n = count_all(GEN),
#'             combination_p = portion_S(CIP, GEN, as_percent = TRUE),
#'             combination_n = count_all(CIP, GEN))
#'
#' # Get portions S/I/R immediately of all rsi columns
#' septic_patients %>%
#'   select(AMX, CIP) %>%
#'   portion_df(translate = FALSE)
#'
#' # It also supports grouping variables
#' septic_patients %>%
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
                      also_single_tested = FALSE) {
  rsi_calc(...,
           type = "R",
           include_I = FALSE,
           minimum = minimum,
           as_percent = as_percent,
           also_single_tested = also_single_tested,
           only_count = FALSE)
}

#' @rdname portion
#' @export
portion_IR <- function(...,
                       minimum = 30,
                       as_percent = FALSE,
                       also_single_tested = FALSE) {
  rsi_calc(...,
           type = "R",
           include_I = TRUE,
           minimum = minimum,
           as_percent = as_percent,
           also_single_tested = also_single_tested,
           only_count = FALSE)
}

#' @rdname portion
#' @export
portion_I <- function(...,
                      minimum = 30,
                      as_percent = FALSE,
                      also_single_tested = FALSE) {
  rsi_calc(...,
           type = "I",
           include_I = FALSE,
           minimum = minimum,
           as_percent = as_percent,
           also_single_tested = also_single_tested,
           only_count = FALSE)
}

#' @rdname portion
#' @export
portion_SI <- function(...,
                       minimum = 30,
                       as_percent = FALSE,
                       also_single_tested = FALSE) {
  rsi_calc(...,
           type = "S",
           include_I = TRUE,
           minimum = minimum,
           as_percent = as_percent,
           also_single_tested = also_single_tested,
           only_count = FALSE)
}

#' @rdname portion
#' @export
portion_S <- function(...,
                      minimum = 30,
                      as_percent = FALSE,
                      also_single_tested = FALSE) {
  rsi_calc(...,
           type = "S",
           include_I = FALSE,
           minimum = minimum,
           as_percent = as_percent,
           also_single_tested = also_single_tested,
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
