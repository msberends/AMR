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
#' @details These functions are meant to count isolates. Use the \code{\link{portion}_*} functions to calculate microbial resistance.
#'
#' \code{n_rsi} is an alias of \code{count_all}. They can be used to count all available isolates, i.e. where all input antibiotics have an available result (S, I or R). Their use is equal to \code{\link{n_distinct}}. Their function is equal to \code{count_S(...) + count_IR(...)}.
#'
#' \code{count_df} takes any variable from \code{data} that has an \code{"rsi"} class (created with \code{\link{as.rsi}}) and counts the amounts of R, I and S. The resulting \emph{tidy data} (see Source) \code{data.frame} will have three rows (S/I/R) and a column for each variable with class \code{"rsi"}.
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
#' count_IR(septic_patients$AMX)
#' portion_IR(septic_patients$AMX) * n_rsi(septic_patients$AMX)
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
#' # Please mind that `portion_S` calculates percentages right away instead.
#' count_S(septic_patients$AMC)   # S = 1342 (71.4%)
#' count_all(septic_patients$AMC) # n = 1879
#'
#' count_S(septic_patients$GEN)   # S = 1372 (74.0%)
#' count_all(septic_patients$GEN) # n = 1855
#'
#' with(septic_patients,
#'      count_S(AMC, GEN))         # S = 1660 (92.3%)
#' with(septic_patients,           # n = 1798
#'      n_rsi(AMC, GEN))
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
count_R <- function(..., also_single_tested = FALSE) {
  rsi_calc(...,
           type = "R",
           include_I = FALSE,
           minimum = 0,
           as_percent = FALSE,
           also_single_tested = FALSE,
           only_count = TRUE)
}

#' @rdname count
#' @export
count_IR <- function(..., also_single_tested = FALSE) {
  rsi_calc(...,
           type = "R",
           include_I = TRUE,
           minimum = 0,
           as_percent = FALSE,
           also_single_tested = FALSE,
           only_count = TRUE)
}

#' @rdname count
#' @export
count_I <- function(..., also_single_tested = FALSE) {
  rsi_calc(...,
           type = "I",
           include_I = FALSE,
           minimum = 0,
           as_percent = FALSE,
           also_single_tested = FALSE,
           only_count = TRUE)
}

#' @rdname count
#' @export
count_SI <- function(..., also_single_tested = FALSE) {
  rsi_calc(...,
           type = "S",
           include_I = TRUE,
           minimum = 0,
           as_percent = FALSE,
           also_single_tested = FALSE,
           only_count = TRUE)
}

#' @rdname count
#' @export
count_S <- function(..., also_single_tested = FALSE) {
  rsi_calc(...,
           type = "S",
           include_I = FALSE,
           minimum = 0,
           as_percent = FALSE,
           also_single_tested = FALSE,
           only_count = TRUE)
}

#' @rdname count
#' @export
count_all <- function(...) {
  # only print warnings once, if needed
  count_S(...) + suppressWarnings(count_IR(...))
}

#' @rdname count
#' @export
n_rsi <- function(...) {
  # only print warnings once, if needed
  count_S(...) + suppressWarnings(count_IR(...))
}

#' @rdname count
#' @importFrom dplyr %>% select_if bind_rows summarise_if mutate group_vars select everything
#' @export
count_df <- function(data,
                     translate_ab = "name",
                     language = get_locale(),
                     combine_IR = FALSE) {

  if (!"data.frame" %in% class(data)) {
    stop("`count_df` must be called on a data.frame")
  }

  if (data %>% select_if(is.rsi) %>% ncol() == 0) {
    stop("No columns with class 'rsi' found. See ?as.rsi.")
  }

  if (as.character(translate_ab) %in% c("TRUE", "official")) {
    translate_ab <- "name"
  }

  resS <- summarise_if(.tbl = data,
                       .predicate = is.rsi,
                       .funs = count_S) %>%
    mutate(Interpretation = "S") %>%
    select(Interpretation, everything())

  if (combine_IR == FALSE) {
    resI <- summarise_if(.tbl = data,
                         .predicate = is.rsi,
                         .funs = count_I) %>%
      mutate(Interpretation = "I") %>%
      select(Interpretation, everything())

    resR <- summarise_if(.tbl = data,
                         .predicate = is.rsi,
                         .funs = count_R) %>%
      mutate(Interpretation = "R") %>%
      select(Interpretation, everything())

    data.groups <- group_vars(data)

    res <- bind_rows(resS, resI, resR) %>%
      mutate(Interpretation = factor(Interpretation, levels = c("R", "I", "S"), ordered = TRUE)) %>%
      tidyr::gather(Antibiotic, Value, -Interpretation, -data.groups)
  } else {
    resIR <- summarise_if(.tbl = data,
                          .predicate = is.rsi,
                          .funs = count_IR) %>%
      mutate(Interpretation = "IR") %>%
      select(Interpretation, everything())

    data.groups <- group_vars(data)

    res <- bind_rows(resS, resIR) %>%
      mutate(Interpretation = factor(Interpretation, levels = c("IR", "S"), ordered = TRUE)) %>%
      tidyr::gather(Antibiotic, Value, -Interpretation, -data.groups)
  }

  if (!translate_ab == FALSE) {
    res <- res %>% mutate(Antibiotic = ab_property(Antibiotic, property = translate_ab, language = language))
  }

  res
}
