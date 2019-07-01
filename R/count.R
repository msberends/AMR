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

#' Count isolates
#'
#' @description These functions can be used to count resistant/susceptible microbial isolates. All functions support quasiquotation with pipes, can be used in \code{dplyr}s \code{\link[dplyr]{summarise}} and support grouped variables, see \emph{Examples}.
#'
#' \code{count_R} and \code{count_IR} can be used to count resistant isolates, \code{count_S} and \code{count_SI} can be used to count susceptible isolates.\cr
#' @param ... one or more vectors (or columns) with antibiotic interpretations. They will be transformed internally with \code{\link{as.rsi}} if needed.
#' @inheritParams portion
#' @inheritSection as.rsi Interpretation of S, I and R
#' @details These functions are meant to count isolates. Use the \code{\link{portion}_*} functions to calculate microbial resistance.
#'
#' The function \code{n_rsi} is an alias of \code{count_all}. They can be used to count all available isolates, i.e. where all input antibiotics have an available result (S, I or R). Their use is equal to \code{\link{n_distinct}}. Their function is equal to \code{count_S(...) + count_IR(...)}.
#'
#' The function \code{count_df} takes any variable from \code{data} that has an \code{"rsi"} class (created with \code{\link{as.rsi}}) and counts the amounts of S, I and R. The resulting \emph{tidy data} (see Source) \code{data.frame} will have three rows (S/I/R) and a column for each variable with class \code{"rsi"}.
#'
#' The function \code{rsi_df} works exactly like \code{count_df}, but adds the percentage of S, I and R.
#' @inheritSection portion Combination therapy
#' @source Wickham H. \strong{Tidy Data.} The Journal of Statistical Software, vol. 59, 2014. \url{http://vita.had.co.nz/papers/tidy-data.html}
#' @seealso \code{\link{portion}_*} to calculate microbial resistance and susceptibility.
#' @keywords resistance susceptibility rsi antibiotics isolate isolates
#' @return Integer
#' @rdname count
#' @name count
#' @export
#' @inheritSection AMR Read more on our website!
#' @examples
#' # septic_patients is a data set available in the AMR package. It is true, genuine data.
#' ?septic_patients
#'
#' # Count resistant isolates
#' count_R(septic_patients$AMX)
#' count_IR(septic_patients$AMX)
#'
#' # Or susceptible isolates
#' count_S(septic_patients$AMX)
#' count_SI(septic_patients$AMX)
#'
#' # Count all available isolates
#' count_all(septic_patients$AMX)
#' n_rsi(septic_patients$AMX)
#'
#' # Since n_rsi counts available isolates, you can
#' # calculate back to count e.g. non-susceptible isolates.
#' # This results in the same:
#' count_SI(septic_patients$AMX)
#' portion_SI(septic_patients$AMX) * n_rsi(septic_patients$AMX)
#'
#' library(dplyr)
#' septic_patients %>%
#'   group_by(hospital_id) %>%
#'   summarise(R  = count_R(CIP),
#'             I  = count_I(CIP),
#'             S  = count_S(CIP),
#'             n1 = count_all(CIP),  # the actual total; sum of all three
#'             n2 = n_rsi(CIP),      # same - analogous to n_distinct
#'             total = n())          # NOT the number of tested isolates!
#'
#' # Count co-resistance between amoxicillin/clav acid and gentamicin,
#' # so we can see that combination therapy does a lot more than mono therapy.
#' # Please mind that `portion_SI` calculates percentages right away instead.
#' count_SI(septic_patients$AMC)  # 1433
#' count_all(septic_patients$AMC) # 1879
#'
#' count_SI(septic_patients$GEN)  # 1399
#' count_all(septic_patients$GEN) # 1855
#'
#' with(septic_patients,
#'      count_SI(AMC, GEN))        # 1764
#' with(septic_patients,
#'      n_rsi(AMC, GEN))           # 1936
#'
#' # Get portions S/I/R immediately of all rsi columns
#' septic_patients %>%
#'   select(AMX, CIP) %>%
#'   count_df(translate = FALSE)
#'
#' # It also supports grouping variables
#' septic_patients %>%
#'   select(hospital_id, AMX, CIP) %>%
#'   group_by(hospital_id) %>%
#'   count_df(translate = FALSE)
#'
count_R <- function(..., only_all_tested = FALSE) {
  rsi_calc(...,
           ab_result = "R",
           only_all_tested = only_all_tested,
           only_count = TRUE)
}

#' @rdname count
#' @export
count_IR <- function(..., only_all_tested = FALSE) {
  rsi_calc(...,
           ab_result = c("I", "R"),
           only_all_tested = only_all_tested,
           only_count = TRUE)
}

#' @rdname count
#' @export
count_I <- function(..., only_all_tested = FALSE) {
  rsi_calc(...,
           ab_result = "I",
           only_all_tested = only_all_tested,
           only_count = TRUE)
}

#' @rdname count
#' @export
count_SI <- function(..., only_all_tested = FALSE) {
  rsi_calc(...,
           ab_result = c("S", "I"),
           only_all_tested = only_all_tested,
           only_count = TRUE)
}

#' @rdname count
#' @export
count_S <- function(..., only_all_tested = FALSE) {
  rsi_calc(...,
           ab_result = "S",
           only_all_tested = only_all_tested,
           only_count = TRUE)
}

#' @rdname count
#' @export
count_all <- function(..., only_all_tested = FALSE) {
  rsi_calc(...,
           ab_result = c("S", "I", "R"),
           only_all_tested = only_all_tested,
           only_count = TRUE)
}

#' @rdname count
#' @export
n_rsi<- count_all

#' @rdname count
#' @export
count_df <- function(data,
                     translate_ab = "name",
                     language = get_locale(),
                     combine_SI = TRUE,
                     combine_IR = FALSE) {

  rsi_calc_df(type = "count",
              data = data,
              translate_ab = translate_ab,
              language = language,
              combine_SI = combine_SI,
              combine_IR = combine_IR,
              combine_SI_missing = missing(combine_SI))
}
