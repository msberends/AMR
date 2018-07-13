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
#' These functions can be used to calculate the (co-)resistance of microbial isolates (i.e. percentage S, SI, I, IR or R). All functions can be used in \code{dplyr}s \code{\link[dplyr]{summarise}} and support grouped variables, see \emph{Examples}.
#' @param ab,ab1,ab2 vector of antibiotic interpretations, they will be transformed internally with \code{\link{as.rsi}}
#' @param include_I logical to indicate whether antimicrobial interpretations of "I" should be included
#' @param minimum minimal amount of available isolates. Any number lower than \code{minimum} will return \code{NA}.
#' @param as_percent logical to indicate whether the output must be returned as percent (text), will else be a double
#' @param interpretation antimicrobial interpretation
#' @details \strong{Remember that you should filter your table to let it contain only first isolates!} Use \code{\link{first_isolate}} to determine them in your data set.
#'
#' All return values are calculated using hybrid evaluation (i.e. using C++), which makes these functions 60-65 times faster than in \code{AMR} v0.2.0 and below. The \code{rsi} function is available for backwards compatibility and deprecated. It now uses the \code{resistance} and \code{susceptibility} functions internally, based on the \code{interpretation} parameter.
#' \if{html}{
#'   \cr
#'   To calculate the probability (\emph{p}) of susceptibility of one antibiotic, we use this formula:
#'   \out{<div style="text-align: center">}\figure{mono_therapy.png}\out{</div>}
#'   To calculate the probability (\emph{p}) of susceptibility of more antibiotics (i.e. combination therapy), we need to check whether one of them has a susceptible result (as numerator) and count all cases where all antibiotics were tested (as denominator). \cr
#'   \cr
#'   For two antibiotics:
#'   \out{<div style="text-align: center">}\figure{combi_therapy_2.png}\out{</div>}
#'   \cr
#'   Theoretically for three antibiotics:
#'   \out{<div style="text-align: center">}\figure{combi_therapy_3.png}\out{</div>}
#' }
#' @keywords resistance susceptibility rsi_df antibiotics isolate isolates
#' @return Double or, when \code{as_percent = TRUE}, a character.
#' @rdname resistance
#' @export
#' @examples
#' library(dplyr)
#'
#' septic_patients %>%
#'   group_by(hospital_id) %>%
#'   summarise(p = susceptibility(cipr),
#'             n = n_rsi(cipr)) # n_rsi works like n_distinct in dplyr
#'
#' septic_patients %>%
#'   group_by(hospital_id) %>%
#'   summarise(cipro_p = susceptibility(cipr, as_percent = TRUE),
#'             cipro_n = n_rsi(cipr),
#'             genta_p = susceptibility(gent, as_percent = TRUE),
#'             genta_n = n_rsi(gent),
#'             combination_p = susceptibility(cipr, gent, as_percent = TRUE),
#'             combination_n = n_rsi(cipr, gent))
#'
#'
#' # Calculate resistance
#' resistance(septic_patients$amox)
#' rsi(septic_patients$amox, interpretation = "IR") # deprecated
#'
#' # Or susceptibility
#' susceptibility(septic_patients$amox)
#' rsi(septic_patients$amox, interpretation = "S") # deprecated
#'
#'
#' # Calculate co-resistance between amoxicillin/clav acid and gentamicin,
#' # so we can see that combination therapy does a lot more than mono therapy:
#' susceptibility(septic_patients$amcl) # p = 67.8%
#' n_rsi(septic_patients$amcl)          # n = 1641
#'
#' susceptibility(septic_patients$gent) # p = 69.1%
#' n_rsi(septic_patients$gent)          # n = 1863
#'
#' with(septic_patients,
#'      susceptibility(amcl, gent))     # p = 90.6%
#' with(septic_patients,
#'      n_rsi(amcl, gent))              # n = 1580
#'
#' \dontrun{
#' # calculate current empiric combination therapy of Helicobacter gastritis:
#' my_table %>%
#'   filter(first_isolate == TRUE,
#'          genus == "Helicobacter") %>%
#'   summarise(p = susceptibility(amox, metr), # amoxicillin with metronidazole
#'             n = n_rsi(amox, metr))
#' }
resistance <- function(ab,
                       include_I = TRUE,
                       minimum = 30,
                       as_percent = FALSE) {

  if (NCOL(ab) > 1) {
    stop('`ab` must be a vector of antimicrobial interpretations', call. = FALSE)
  }
  if (!is.logical(include_I)) {
    stop('`include_I` must be logical', call. = FALSE)
  }
  if (!is.numeric(minimum)) {
    stop('`minimum` must be numeric', call. = FALSE)
  }
  if (!is.logical(as_percent)) {
    stop('`as_percent` must be logical', call. = FALSE)
  }

  x <- as.integer(as.rsi(ab))
  total <-  .Call(`_AMR_rsi_calc_total`, x)
  if (total < minimum) {
    return(NA)
  }
  found <- .Call(`_AMR_rsi_calc_R`, x, include_I)

  if (as_percent == TRUE) {
    percent(found / total, force_zero = TRUE)
  } else {
    found / total
  }
}

#' @rdname resistance
#' @export
susceptibility <- function(ab1,
                           ab2 = NULL,
                           include_I = FALSE,
                           minimum = 30,
                           as_percent = FALSE) {

  if (NCOL(ab1) > 1) {
    stop('`ab1` must be a vector of antimicrobial interpretations', call. = FALSE)
  }
  if (!is.logical(include_I)) {
    stop('`include_I` must be logical', call. = FALSE)
  }
  if (!is.numeric(minimum)) {
    stop('`minimum` must be numeric', call. = FALSE)
  }
  if (!is.logical(as_percent)) {
    stop('`as_percent` must be logical', call. = FALSE)
  }

  if (!is.null(ab2)) {
    if (NCOL(ab2) > 1) {
      stop('`ab2` must be a vector of antimicrobial interpretations', call. = FALSE)
    }
    x <- apply(X = data.frame(ab1 = as.integer(as.rsi(ab1)),
                              ab2 = as.integer(as.rsi(ab2))),
               MARGIN = 1,
               FUN = min)
  } else {
    x <- as.integer(as.rsi(ab1))
  }
  total <-  .Call(`_AMR_rsi_calc_total`, x)
  if (total < minimum) {
    return(NA)
  }
  found <- .Call(`_AMR_rsi_calc_S`, x, include_I)

  if (as_percent == TRUE) {
    percent(found / total, force_zero = TRUE)
  } else {
    found / total
  }
}

#' @rdname resistance
#' @export
n_rsi <- function(ab1, ab2 = NULL) {
  if (NCOL(ab1) > 1) {
    stop('`ab1` must be a vector of antimicrobial interpretations', call. = FALSE)
  }
  if (!is.null(ab2)) {
    if (NCOL(ab2) > 1) {
      stop('`ab2` must be a vector of antimicrobial interpretations', call. = FALSE)
    }
    x <- apply(X = data.frame(ab1 = as.integer(as.rsi(ab1)),
                              ab2 = as.integer(as.rsi(ab2))),
               MARGIN = 1,
               FUN = min)
  } else {
    x <- as.integer(as.rsi(ab1))
  }
  .Call(`_AMR_rsi_calc_total`, x)
}


#' @rdname resistance
#' @export
rsi <- function(ab1,
                ab2 = NULL,
                interpretation = "IR",
                minimum = 30,
                as_percent = FALSE) {
  warning("'rsi' is deprecated. Use 'resistance' or 'susceptibility' instead.", call. = FALSE)
  if (interpretation %in% c('IR', 'RI')) {
    resistance(ab = ab1, include_I = TRUE, minimum = minimum, as_percent = as_percent)
  } else if (interpretation == 'R') {
    resistance(ab = ab1, include_I = FALSE, minimum = minimum, as_percent = as_percent)
  } else  if (interpretation %in% c('IS', 'SI')) {
    susceptibility(ab1 = ab1, ab2 = ab2, include_I = TRUE, minimum = minimum, as_percent = as_percent)
  } else if (interpretation == 'S') {
    susceptibility(ab1 = ab1, ab2 = ab2, include_I = FALSE, minimum = minimum, as_percent = as_percent)
  } else {
    stop('invalid `interpretation`')
  }
}

#' Predict antimicrobial resistance
#'
#' Create a prediction model to predict antimicrobial resistance for the next years on statistical solid ground. Standard errors (SE) will be returned as columns \code{se_min} and \code{se_max}. See Examples for a real live example.
#' @param tbl table that contains columns \code{col_ab} and \code{col_date}
#' @param col_ab column name of \code{tbl} with antimicrobial interpretations (\code{R}, \code{I} and \code{S}), supports tidyverse-like quotation
#' @param col_date column name of the date, will be used to calculate years if this column doesn't consist of years already, supports tidyverse-like quotation
#' @param year_max highest year to use in the prediction model, deafults to 15 years after today
#' @param year_every unit of sequence between lowest year found in the data and \code{year_max}
#' @param model the statistical model of choice. Valid values are \code{"binomial"} (or \code{"binom"} or \code{"logit"}) or \code{"loglin"} or \code{"linear"} (or \code{"lin"}).
#' @param I_as_R treat \code{I} as \code{R}
#' @param preserve_measurements overwrite predictions of years that are actually available in the data, with the original data. The standard errors of those years will be \code{NA}.
#' @param info print textual analysis with the name and \code{\link{summary}} of the model.
#' @return \code{data.frame} with columns \code{year}, \code{probR}, \code{se_min} and \code{se_max}.
#' @seealso \code{\link{resistance}} \cr \code{\link{lm}} \code{\link{glm}}
#' @rdname resistance_predict
#' @export
#' @importFrom dplyr %>% pull mutate group_by_at summarise filter
#' @importFrom reshape2 dcast
#' @examples
#' \dontrun{
#' # use it directly:
#' rsi_predict(tbl = tbl[which(first_isolate == TRUE & genus == "Haemophilus"),],
#'             col_ab = "amcl", col_date = "date")
#'
#' # or with dplyr so you can actually read it:
#' library(dplyr)
#' tbl %>%
#'   filter(first_isolate == TRUE,
#'          genus == "Haemophilus") %>%
#'   rsi_predict(amcl, date)
#' }
#'
#'
#' # real live example:
#' library(dplyr)
#' septic_patients %>%
#'   # get bacteria properties like genus and species
#'   left_join_microorganisms("bactid") %>%
#'   # calculate first isolates
#'   mutate(first_isolate =
#'            first_isolate(.,
#'                          "date",
#'                          "patient_id",
#'                          "bactid",
#'                          col_specimen = NA,
#'                          col_icu = NA)) %>%
#'   # filter on first E. coli isolates
#'   filter(genus == "Escherichia",
#'          species == "coli",
#'          first_isolate == TRUE) %>%
#'   # predict resistance of cefotaxime for next years
#'   rsi_predict(col_ab = "cfot",
#'               col_date = "date",
#'               year_max = 2025,
#'               preserve_measurements = FALSE)
#'
resistance_predict <- function(tbl,
                        col_ab,
                        col_date,
                        year_max = as.integer(format(as.Date(Sys.Date()), '%Y')) + 15,
                        year_every = 1,
                        model = 'binomial',
                        I_as_R = TRUE,
                        preserve_measurements = TRUE,
                        info = TRUE) {

  if (nrow(tbl) == 0) {
    stop('This table does not contain any observations.')
  }

  if (!col_ab %in% colnames(tbl)) {
    stop('Column ', col_ab, ' not found.')
  }

  if (!col_date %in% colnames(tbl)) {
    stop('Column ', col_date, ' not found.')
  }
  if ('grouped_df' %in% class(tbl)) {
    # no grouped tibbles please, mutate will throw errors
    tbl <- base::as.data.frame(tbl, stringsAsFactors = FALSE)
  }

  if (I_as_R == TRUE) {
    tbl[, col_ab] <- gsub('I', 'R', tbl %>% pull(col_ab))
  }

  if (!all(tbl %>% pull(col_ab) %>% as.rsi() %in% c(NA, 'S', 'I', 'R'))) {
    stop('Column ', col_ab, ' must contain antimicrobial interpretations (S, I, R).')
  }

  year <- function(x) {
    if (all(grepl('^[0-9]{4}$', x))) {
      x
    } else {
      as.integer(format(as.Date(x), '%Y'))
    }
  }

  years_predict <- seq(from = min(year(tbl %>% pull(col_date))), to = year_max, by = year_every)

  df <- tbl %>%
    mutate(year = year(tbl %>% pull(col_date))) %>%
    group_by_at(c('year', col_ab)) %>%
    summarise(n())
  colnames(df) <- c('year', 'antibiotic', 'count')
  df <- df %>%
    reshape2::dcast(year ~ antibiotic, value.var = 'count')

  if (model %in% c('binomial', 'binom', 'logit')) {
    logitmodel <- with(df, glm(cbind(R, S) ~ year, family = binomial))
    if (info == TRUE) {
      cat('\nLogistic regression model (logit) with binomial distribution')
      cat('\n------------------------------------------------------------\n')
      print(summary(logitmodel))
    }

    predictmodel <- stats::predict(logitmodel, newdata = with(df, list(year = years_predict)), type = "response", se.fit = TRUE)
    prediction <- predictmodel$fit
    se <- predictmodel$se.fit

  } else if (model == 'loglin') {
    loglinmodel <- with(df, glm(R ~ year, family = poisson))
    if (info == TRUE) {
      cat('\nLog-linear regression model (loglin) with poisson distribution')
      cat('\n--------------------------------------------------------------\n')
      print(summary(loglinmodel))
    }

    predictmodel <- stats::predict(loglinmodel, newdata = with(df, list(year = years_predict)), type = "response", se.fit = TRUE)
    prediction <- predictmodel$fit
    se <- predictmodel$se.fit

  } else if (model %in% c('lin', 'linear')) {
    linmodel <- with(df, lm((R / (R + S)) ~ year))
    if (info == TRUE) {
      cat('\nLinear regression model')
      cat('\n-----------------------\n')
      print(summary(linmodel))
    }

    predictmodel <- stats::predict(linmodel, newdata = with(df, list(year = years_predict)), se.fit = TRUE)
    prediction <- predictmodel$fit
    se <- predictmodel$se.fit

  } else {
    stop('No valid model selected.')
  }

  # prepare the output dataframe
  prediction <- data.frame(year = years_predict, probR = prediction, stringsAsFactors = FALSE)

  prediction$se_min <- prediction$probR - se
  prediction$se_max <- prediction$probR + se

  if (model == 'loglin') {
    prediction$probR <- prediction$probR %>%
      format(scientific = FALSE) %>%
      as.integer()
    prediction$se_min <- prediction$se_min %>% as.integer()
    prediction$se_max <- prediction$se_max %>% as.integer()

    colnames(prediction) <- c('year', 'amountR', 'se_max', 'se_min')
  } else {
    prediction$se_max[which(prediction$se_max > 1)] <- 1
  }
  prediction$se_min[which(prediction$se_min < 0)] <- 0

  total <- prediction

  if (preserve_measurements == TRUE) {
    # geschatte data vervangen door gemeten data
    if (I_as_R == TRUE) {
      if (!'I' %in% colnames(df)) {
        df$I <- 0
      }
      df$probR <- df$R / rowSums(df[, c('R', 'S', 'I')])
    } else {
      df$probR <- df$R / rowSums(df[, c('R', 'S')])
    }
    measurements <- data.frame(year = df$year,
                           probR = df$probR,
                           se_min = NA,
                           se_max = NA,
                           stringsAsFactors = FALSE)
    colnames(measurements) <- colnames(prediction)
    prediction <- prediction %>% filter(!year %in% df$year)

    total <- rbind(measurements, prediction)
  }

  total

}

#' @rdname resistance_predict
#' @export
rsi_predict <- resistance_predict
