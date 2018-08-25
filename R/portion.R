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

#' Calculate resistance of isolates
#'
#' @description These functions can be used to calculate the (co-)resistance of microbial isolates (i.e. percentage S, SI, I, IR or R). All functions support quasiquotation with pipes, can be used in \code{dplyr}s \code{\link[dplyr]{summarise}} and support grouped variables, see \emph{Examples}.
#'
#' \code{portion_R} and \code{portion_IR} can be used to calculate resistance, \code{portion_S} and \code{portion_SI} can be used to calculate susceptibility.\cr
#' @param ... one or more vectors (or columns) with antibiotic interpretations. They will be transformed internally with \code{\link{as.rsi}} if needed. Use multiple columns to calculate (the lack of) co-resistance: the probability where one of two drugs have a resistant or susceptible result. See Examples.
#' @param minimum minimal amount of available isolates. Any number lower than \code{minimum} will return \code{NA}. The default number of \code{30} isolates is advised by the CLSI as best practice, see Source.
#' @param as_percent logical to indicate whether the output must be returned as a hundred fold with \% sign (a character). A value of \code{0.123456} will then be returned as \code{"12.3\%"}.
#' @param data a \code{data.frame} containing columns with class \code{rsi} (see \code{\link{as.rsi}})
#' @param translate_ab a column name of the \code{\link{antibiotics}} data set to translate the antibiotic abbreviations to, using \code{\link{abname}}. This can be set with \code{\link{getOption}("get_antibiotic_names")}.
#' @details \strong{Remember that you should filter your table to let it contain only first isolates!} Use \code{\link{first_isolate}} to determine them in your data set.
#'
#' These functions are not meant to count isolates, but to calculate the portion of resistance/susceptibility. If a column has been transformed with \code{\link{as.rsi}}, just use e.g. \code{isolates[isolates == "R"]} to get the resistant ones. You could then calculate the \code{\link{length}} of it.
#'
#' \code{portion_df} takes any variable from \code{data} that has an \code{"rsi"} class (created with \code{\link{as.rsi}}) and calculates the portions R, I and S. The resulting \emph{tidy data} (see Source) \code{data.frame} will have three rows (S/I/R) and a column for each variable with class \code{"rsi"}.
#'
#' The old \code{\link{rsi}} function is still available for backwards compatibility but is deprecated.
#' \if{html}{
#'   \cr\cr
#'   To calculate the probability (\emph{p}) of susceptibility of one antibiotic, we use this formula:
#'   \out{<div style="text-align: center">}\figure{mono_therapy.png}\out{</div>}
#'   To calculate the probability (\emph{p}) of susceptibility of more antibiotics (i.e. combination therapy), we need to check whether one of them has a susceptible result (as numerator) and count all cases where all antibiotics were tested (as denominator). \cr
#'   \cr
#'   For two antibiotics:
#'   \out{<div style="text-align: center">}\figure{combi_therapy_2.png}\out{</div>}
#'   \cr
#'   For three antibiotics:
#'   \out{<div style="text-align: center">}\figure{combi_therapy_3.png}\out{</div>}
#'   \cr
#'   And so on.
#' }
#' @source \strong{M39 Analysis and Presentation of Cumulative Antimicrobial Susceptibility Test Data, 4th Edition}, 2014, \emph{Clinical and Laboratory Standards Institute (CLSI)}. \url{https://clsi.org/standards/products/microbiology/documents/m39/}.
#'
#' Wickham H. \strong{Tidy Data.} The Journal of Statistical Software, vol. 59, 2014. \url{http://vita.had.co.nz/papers/tidy-data.html}
#' @seealso \code{\link[AMR]{count}_*} to count resistant and susceptibile isolates.\cr
#' \code{\link{n_rsi}} to count all cases where antimicrobial results are available.
#' @keywords resistance susceptibility rsi_df rsi antibiotics isolate isolates
#' @return Double or, when \code{as_percent = TRUE}, a character.
#' @rdname portion
#' @name portion
#' @export
#' @examples
#' # septic_patients is a data set available in the AMR package. It is true, genuine data.
#' ?septic_patients
#'
#' # Calculate resistance
#' portion_R(septic_patients$amox)
#' portion_IR(septic_patients$amox)
#'
#' # Or susceptibility
#' portion_S(septic_patients$amox)
#' portion_SI(septic_patients$amox)
#'

#' # Do the above with pipes:
#' library(dplyr)
#' septic_patients %>% portion_R(amox)
#' septic_patients %>% portion_IR(amox)
#' septic_patients %>% portion_S(amox)
#' septic_patients %>% portion_SI(amox)
#'
#' septic_patients %>%
#'   group_by(hospital_id) %>%
#'   summarise(p = portion_S(cipr),
#'             n = n_rsi(cipr)) # n_rsi works like n_distinct in dplyr
#'
#' septic_patients %>%
#'   group_by(hospital_id) %>%
#'   summarise(R = portion_R(cipr, as_percent = TRUE),
#'             I = portion_I(cipr, as_percent = TRUE),
#'             S = portion_S(cipr, as_percent = TRUE),
#'             n = n_rsi(cipr), # works like n_distinct in dplyr
#'             total = n())     # NOT the amount of tested isolates!
#'
#' # Calculate co-resistance between amoxicillin/clav acid and gentamicin,
#' # so we can see that combination therapy does a lot more than mono therapy:
#' septic_patients %>% portion_S(amcl)       # S = 67.3%
#' septic_patients %>% n_rsi(amcl)           # n = 1570
#'
#' septic_patients %>% portion_S(gent)       # S = 74.0%
#' septic_patients %>% n_rsi(gent)           # n = 1842
#'
#' septic_patients %>% portion_S(amcl, gent) # S = 92.1%
#' septic_patients %>% n_rsi(amcl, gent)     # n = 1504
#'
#'
#' septic_patients %>%
#'   group_by(hospital_id) %>%
#'   summarise(cipro_p = portion_S(cipr, as_percent = TRUE),
#'             cipro_n = n_rsi(cipr),
#'             genta_p = portion_S(gent, as_percent = TRUE),
#'             genta_n = n_rsi(gent),
#'             combination_p = portion_S(cipr, gent, as_percent = TRUE),
#'             combination_n = n_rsi(cipr, gent))
#'
#' # Get portions S/I/R immediately of all rsi columns
#' septic_patients %>%
#'   select(amox, cipr) %>%
#'   portion_df(translate = FALSE)
#'
#' # It also supports grouping variables
#' septic_patients %>%
#'   select(hospital_id, amox, cipr) %>%
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
#'   summarise(p = portion_S(amox, metr),  # amoxicillin with metronidazole
#'             n = n_rsi(amox, metr))
#' }
portion_R <- function(...,
                      minimum = 30,
                      as_percent = FALSE) {
  rsi_calc(...,
           type = "R",
           include_I = FALSE,
           minimum = minimum,
           as_percent = as_percent,
           only_count = FALSE)
}

#' @rdname portion
#' @export
portion_IR <- function(...,
                       minimum = 30,
                       as_percent = FALSE) {
  rsi_calc(...,
           type = "R",
           include_I = TRUE,
           minimum = minimum,
           as_percent = as_percent,
           only_count = FALSE)
}

#' @rdname portion
#' @export
portion_I <- function(...,
                      minimum = 30,
                      as_percent = FALSE) {
  rsi_calc(...,
           type = "I",
           include_I = FALSE,
           minimum = minimum,
           as_percent = as_percent,
           only_count = FALSE)
}

#' @rdname portion
#' @export
portion_SI <- function(...,
                       minimum = 30,
                       as_percent = FALSE) {
  rsi_calc(...,
           type = "S",
           include_I = TRUE,
           minimum = minimum,
           as_percent = as_percent,
           only_count = FALSE)
}

#' @rdname portion
#' @export
portion_S <- function(...,
                      minimum = 30,
                      as_percent = FALSE) {
  rsi_calc(...,
           type = "S",
           include_I = FALSE,
           minimum = minimum,
           as_percent = as_percent,
           only_count = FALSE)
}

#' @rdname portion
#' @importFrom dplyr %>% select_if bind_rows summarise_if mutate group_vars select everything
#' @export
portion_df <- function(data,
                       translate_ab = getOption("get_antibiotic_names", "official"),
                       minimum = 30,
                       as_percent = FALSE) {

  if (data %>% select_if(is.rsi) %>% ncol() == 0) {
    stop("No columns with class 'rsi' found. See ?as.rsi.")
  }

  if (as.character(translate_ab) == "TRUE") {
    translate_ab <- "official"
  }
  options(get_antibiotic_names = translate_ab)

  resS <- summarise_if(.tbl = data,
                       .predicate = is.rsi,
                       .funs = portion_S,
                       minimum = minimum,
                       as_percent = as_percent) %>%
    mutate(Interpretation = "S") %>%
    select(Interpretation, everything())

  resI <- summarise_if(.tbl = data,
                       .predicate = is.rsi,
                       .funs = portion_I,
                       minimum = minimum,
                       as_percent = as_percent) %>%
    mutate(Interpretation = "I") %>%
    select(Interpretation, everything())

  resR <- summarise_if(.tbl = data,
                       .predicate = is.rsi,
                       .funs = portion_R,
                       minimum = minimum,
                       as_percent = as_percent) %>%
    mutate(Interpretation = "R") %>%
    select(Interpretation, everything())

  data.groups <- group_vars(data)

  res <- bind_rows(resS, resI, resR) %>%
    mutate(Interpretation = factor(Interpretation, levels = c("R", "I", "S"), ordered = TRUE)) %>%
    tidyr::gather(Antibiotic, Percentage, -Interpretation, -data.groups)

  if (!translate_ab == FALSE) {
    if (!tolower(translate_ab) %in% tolower(colnames(AMR::antibiotics))) {
      stop("Parameter `translate_ab` does not occur in the `antibiotics` data set.", call. = FALSE)
    }
    res <- res %>% mutate(Antibiotic = abname(Antibiotic, from = "guess", to = translate_ab))
  }

  res
}


#' Calculate resistance of isolates
#'
#' This function is deprecated. Use the \code{\link{portion}} functions instead.
#' @inheritParams portion
#' @param ab1,ab2  vector (or column) with antibiotic interpretations. It will be transformed internally with \code{\link{as.rsi}} if needed.
#' @param interpretation antimicrobial interpretation to check for
#' @param ... deprecated parameters to support usage on older versions
#' @importFrom dplyr tibble case_when
#' @export
rsi <- function(ab1,
                ab2 = NULL,
                interpretation = "IR",
                minimum = 30,
                as_percent = FALSE,
                ...) {

  if (all(is.null(ab2))) {
    df <- tibble(ab1 = ab1)
  } else {
    df <- tibble(ab1 = ab1,
                 ab2 = ab2)
  }

  result <- case_when(
    interpretation == "S"             ~ portion_S(df, minimum = minimum, as_percent = FALSE),
    interpretation %in% c("SI", "IS") ~ portion_SI(df, minimum = minimum, as_percent = FALSE),
    interpretation == "I"             ~ portion_I(df, minimum = minimum, as_percent = FALSE),
    interpretation %in% c("RI", "IR") ~ portion_IR(df, minimum = minimum, as_percent = FALSE),
    interpretation == "R"             ~ portion_R(df, minimum = minimum, as_percent = FALSE),
    TRUE ~ -1
  )
  if (result == -1) {
    stop("invalid interpretation")
  }

  .Deprecated(new = paste0("portion_", interpretation))

  if (as_percent == TRUE) {
    percent(result, force_zero = TRUE)
  } else {
    result
  }
}
