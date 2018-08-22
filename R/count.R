# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Analysis                              #
#                                                                      #
# AUTHORS                                                              #
# Berends MS (m.s.berends@umcg.nl), Luz CF (c.f.luz@umcg.nl)           #
#                                                                      #
# LICENCE                                                              #
# This program is free software; you can redistribute it and/or modify #
# it under the terms of the GNU General Public License version 2.0,    #
# as published by the Free Software Foundation.                        #
#                                                                      #
# This program is distributed in the hope that it will be useful,      #
# but WITHOUT ANY WARRANTY; without even the implied warranty of       #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        #
# GNU General Public License for more details.                         #
# ==================================================================== #

#' Count isolates
#'
#' @description These functions can be used to count resistant/susceptible microbial isolates. All functions support quasiquotation with pipes, can be used in \code{dplyr}s \code{\link[dplyr]{summarise}} and support grouped variables, see \emph{Examples}.
#'
#' \code{count_R} and \code{count_IR} can be used to count resistant isolates, \code{count_S} and \code{count_SI} can be used to count susceptible isolates.\cr
#' @inheritParams portion
#' @details \strong{Remember that you should filter your table to let it contain only first isolates!} Use \code{\link{first_isolate}} to determine them in your data set.
#'
#' These functions are meant to count isolates. Use the \code{\link{portion}_*} functions to calculate microbial resistance.
#'
#' \code{count_df} takes any variable from \code{data} that has an \code{"rsi"} class (created with \code{\link{as.rsi}}) and counts the amounts of R, I and S. The resulting \emph{tidy data} (see Source) \code{data.frame} will have three rows (S/I/R) and a column for each variable with class \code{"rsi"}.
#' @source Wickham H. \strong{Tidy Data.} The Journal of Statistical Software, vol. 59, 2014. \url{http://vita.had.co.nz/papers/tidy-data.html}
#' @seealso \code{\link{portion}_*} to calculate microbial resistance and susceptibility.\cr
#' \code{\link{n_rsi}} to count all cases where antimicrobial results are available.
#' @keywords resistance susceptibility rsi antibiotics isolate isolates
#' @return Integer
#' @rdname count
#' @name count
#' @export
#' @examples
#' # septic_patients is a data set available in the AMR package. It is true, genuine data.
#' ?septic_patients
#'
#' # Count resistant isolates
#' count_R(septic_patients$amox)
#' count_IR(septic_patients$amox)
#'
#' # Or susceptibile isolates
#' count_S(septic_patients$amox)
#' count_SI(septic_patients$amox)
#'
#' # Since n_rsi counts available isolates, you can
#' # calculate back to count e.g. non-susceptible isolates.
#' # This results in the same:
#' count_IR(septic_patients$amox)
#' portion_IR(septic_patients$amox) * n_rsi(septic_patients$amox)
#'
#' library(dplyr)
#' septic_patients %>%
#'   group_by(hospital_id) %>%
#'   summarise(R = count_R(cipr),
#'             I = count_I(cipr),
#'             S = count_S(cipr),
#'             n = n_rsi(cipr), # the actual total; sum of all three
#'             total = n())     # NOT the amount of tested isolates!
#'
#' # Count co-resistance between amoxicillin/clav acid and gentamicin,
#' # so we can see that combination therapy does a lot more than mono therapy.
#' # Please mind that `portion_S` calculates percentages right away instead.
#' count_S(septic_patients$amcl)   # S = 1056 (67.3%)
#' n_rsi(septic_patients$amcl)     # n = 1570
#'
#' count_S(septic_patients$gent)   # S = 1363 (74.0%)
#' n_rsi(septic_patients$gent)     # n = 1842
#'
#' with(septic_patients,
#'      count_S(amcl, gent))       # S = 1385 (92.1%)
#' with(septic_patients,           # n = 1504
#'      n_rsi(amcl, gent))
#'
#' # Get portions S/I/R immediately of all rsi columns
#' septic_patients %>%
#'   select(amox, cipr) %>%
#'   count_df(translate = FALSE)
#'
#' # It also supports grouping variables
#' septic_patients %>%
#'   select(hospital_id, amox, cipr) %>%
#'   group_by(hospital_id) %>%
#'   count_df(translate = FALSE)
#'
count_R <- function(...) {
  rsi_calc(...,
           type = "R",
           include_I = FALSE,
           minimum = 0,
           as_percent = FALSE,
           only_count = TRUE)
}

#' @rdname count
#' @export
count_IR <- function(...) {
  rsi_calc(...,
           type = "R",
           include_I = TRUE,
           minimum = 0,
           as_percent = FALSE,
           only_count = TRUE)
}

#' @rdname count
#' @export
count_I <- function(...) {
  rsi_calc(...,
           type = "I",
           include_I = FALSE,
           minimum = 0,
           as_percent = FALSE,
           only_count = TRUE)
}

#' @rdname count
#' @export
count_SI <- function(...) {
  rsi_calc(...,
           type = "S",
           include_I = TRUE,
           minimum = 0,
           as_percent = FALSE,
           only_count = TRUE)
}

#' @rdname count
#' @export
count_S <- function(...) {
  rsi_calc(...,
           type = "S",
           include_I = FALSE,
           minimum = 0,
           as_percent = FALSE,
           only_count = TRUE)
}

#' @rdname count
#' @importFrom dplyr %>% select_if bind_rows summarise_if mutate group_vars select everything
#' @export
count_df <- function(data,
                     translate_ab = getOption("get_antibiotic_names", "official")) {

  if (data %>% select_if(is.rsi) %>% ncol() == 0) {
    stop("No columns with class 'rsi' found. See ?as.rsi.")
  }

  if (as.character(translate_ab) == "TRUE") {
    translate_ab <- "official"
  }
  options(get_antibiotic_names = translate_ab)

  resS <- summarise_if(.tbl = data,
                       .predicate = is.rsi,
                       .funs = count_S) %>%
    mutate(Interpretation = "S") %>%
    select(Interpretation, everything())

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
    tidyr::gather(Antibiotic, Count, -Interpretation, -data.groups)

  if (!translate_ab == FALSE) {
    if (!tolower(translate_ab) %in% tolower(colnames(AMR::antibiotics))) {
      stop("Parameter `translate_ab` does not occur in the `antibiotics` data set.", call. = FALSE)
    }
    res <- res %>% mutate(Antibiotic = abname(Antibiotic, from = "guess", to = translate_ab))
  }

  res
}
