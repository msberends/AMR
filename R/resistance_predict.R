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

#' Predict antimicrobial resistance
#'
#' Create a prediction model to predict antimicrobial resistance for the next years on statistical solid ground. Standard errors (SE) will be returned as columns \code{se_min} and \code{se_max}. See Examples for a real live example.
#' @inheritParams first_isolate
#' @param col_ab column name of \code{tbl} with antimicrobial interpretations (\code{R}, \code{I} and \code{S})
#' @param col_date column name of the date, will be used to calculate years if this column doesn't consist of years already
#' @param year_min lowest year to use in the prediction model, dafaults the lowest year in \code{col_date}
#' @param year_max highest year to use in the prediction model, defaults to 15 years after today
#' @param year_every unit of sequence between lowest year found in the data and \code{year_max}
#' @param minimum minimal amount of available isolates per year to include. Years containing less observations will be estimated by the model.
#' @param model the statistical model of choice. Valid values are \code{"binomial"} (or \code{"binom"} or \code{"logit"}) or \code{"loglin"} or \code{"linear"} (or \code{"lin"}).
#' @param I_as_R treat \code{I} as \code{R}
#' @param preserve_measurements logical to indicate whether predictions of years that are actually available in the data should be overwritten with the original data. The standard errors of those years will be \code{NA}.
#' @param info print textual analysis with the name and \code{\link{summary}} of the model.
#' @return \code{data.frame} with columns:
#' \itemize{
#'   \item{\code{year}}
#'   \item{\code{value}, the same as \code{estimated} when \code{preserve_measurements = FALSE}, and a combination of \code{observed} and \code{estimated} otherwise}
#'   \item{\code{se_min}, the lower bound of the standard error with a minimum of \code{0}}
#'   \item{\code{se_max} the upper bound of the standard error with a maximum of \code{1}}
#'   \item{\code{observations}, the total number of observations, i.e. S + I + R}
#'   \item{\code{observed}, the original observed values}
#'   \item{\code{estimated}, the estimated values, calculated by the model}
#' }
#' @seealso The \code{\link{portion}} function to calculate resistance, \cr \code{\link{lm}} \code{\link{glm}}
#' @rdname resistance_predict
#' @export
#' @importFrom stats predict glm lm
#' @importFrom dplyr %>% pull mutate group_by_at summarise filter n_distinct arrange case_when
# @importFrom tidyr spread
#' @examples
#' \dontrun{
#' # use it with base R:
#' resistance_predict(tbl = tbl[which(first_isolate == TRUE & genus == "Haemophilus"),],
#'                    col_ab = "amcl", col_date = "date")
#'
#' # or use dplyr so you can actually read it:
#' library(dplyr)
#' tbl %>%
#'   filter(first_isolate == TRUE,
#'          genus == "Haemophilus") %>%
#'   resistance_predict(amcl, date)
#' }
#'
#'
#' # real live example:
#' library(dplyr)
#' septic_patients %>%
#'   # get bacteria properties like genus and species
#'   left_join_microorganisms("mo") %>%
#'   # calculate first isolates
#'   mutate(first_isolate =
#'            first_isolate(.,
#'                          "date",
#'                          "patient_id",
#'                          "mo",
#'                          col_specimen = NA,
#'                          col_icu = NA)) %>%
#'   # filter on first E. coli isolates
#'   filter(genus == "Escherichia",
#'          species == "coli",
#'          first_isolate == TRUE) %>%
#'   # predict resistance of cefotaxime for next years
#'   resistance_predict(col_ab = "cfot",
#'                      col_date = "date",
#'                      year_max = 2025,
#'                      preserve_measurements = TRUE,
#'                      minimum = 0)
#'
#' # create nice plots with ggplot
#' if (!require(ggplot2)) {
#'
#'   data <- septic_patients %>%
#'     filter(mo == "ESCCOL") %>%
#'     resistance_predict(col_ab = "amox",
#'                       col_date = "date",
#'                       info = FALSE,
#'                       minimum = 15)
#'
#'   ggplot(data,
#'          aes(x = year)) +
#'     geom_col(aes(y = value),
#'              fill = "grey75") +
#'     geom_errorbar(aes(ymin = se_min,
#'                       ymax = se_max),
#'                   colour = "grey50") +
#'     scale_y_continuous(limits = c(0, 1),
#'                        breaks = seq(0, 1, 0.1),
#'                        labels = paste0(seq(0, 100, 10), "%")) +
#'     labs(title = expression(paste("Forecast of amoxicillin resistance in ",
#'                                   italic("E. coli"))),
#'          y = "%IR",
#'          x = "Year") +
#'     theme_minimal(base_size = 13)
#' }
resistance_predict <- function(tbl,
                               col_ab,
                               col_date,
                               year_min = NULL,
                               year_max = NULL,
                               year_every = 1,
                               minimum = 30,
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

  if (!tbl %>% pull(col_ab) %>% is.rsi()) {
    tbl[, col_ab] <- tbl %>% pull(col_ab) %>% as.rsi()
  }

  year <- function(x) {
    if (all(grepl('^[0-9]{4}$', x))) {
      x
    } else {
      as.integer(format(as.Date(x), '%Y'))
    }
  }

  df <- tbl %>%
    mutate(year = tbl %>% pull(col_date) %>% year()) %>%
    group_by_at(c('year', col_ab)) %>%
    summarise(n())

  if (df %>% pull(col_ab) %>% n_distinct(na.rm = TRUE) < 2) {
    stop("No variety in antimicrobial interpretations - all isolates are '",
         df %>% pull(col_ab) %>% unique() %>% .[!is.na(.)], "'.",
         call. = FALSE)
  }

  colnames(df) <- c('year', 'antibiotic', 'observations')
  df <- df %>%
    filter(!is.na(antibiotic)) %>%
    tidyr::spread(antibiotic, observations, fill = 0) %>%
    mutate(total = R + S) %>%
    filter(total >= minimum)

  if (NROW(df) == 0) {
    stop('There are no observations.')
  }

  year_lowest <- min(df$year)
  if (is.null(year_min)) {
    year_min <- year_lowest
  } else {
    year_min <- max(year_min, year_lowest, na.rm = TRUE)
  }
  if (is.null(year_max)) {
    year_max <- year(Sys.Date()) + 15
  }

  years_predict <- seq(from = year_min, to = year_max, by = year_every)

  if (model %in% c('binomial', 'binom', 'logit')) {
    logitmodel <- with(df, glm(cbind(R, S) ~ year, family = binomial))
    if (info == TRUE) {
      cat('\nLogistic regression model (logit) with binomial distribution')
      cat('\n------------------------------------------------------------\n')
      print(summary(logitmodel))
    }

    predictmodel <- predict(logitmodel, newdata = with(df, list(year = years_predict)), type = "response", se.fit = TRUE)
    prediction <- predictmodel$fit
    se <- predictmodel$se.fit

  } else if (model == 'loglin') {
    loglinmodel <- with(df, glm(R ~ year, family = poisson))
    if (info == TRUE) {
      cat('\nLog-linear regression model (loglin) with poisson distribution')
      cat('\n--------------------------------------------------------------\n')
      print(summary(loglinmodel))
    }

    predictmodel <- predict(loglinmodel, newdata = with(df, list(year = years_predict)), type = "response", se.fit = TRUE)
    prediction <- predictmodel$fit
    se <- predictmodel$se.fit

  } else if (model %in% c('lin', 'linear')) {
    linmodel <- with(df, lm((R / (R + S)) ~ year))
    if (info == TRUE) {
      cat('\nLinear regression model')
      cat('\n-----------------------\n')
      print(summary(linmodel))
    }

    predictmodel <- predict(linmodel, newdata = with(df, list(year = years_predict)), se.fit = TRUE)
    prediction <- predictmodel$fit
    se <- predictmodel$se.fit

  } else {
    stop('No valid model selected.')
  }

  # prepare the output dataframe
  prediction <- data.frame(year = years_predict, value = prediction, stringsAsFactors = FALSE)

  prediction$se_min <- prediction$value - se
  prediction$se_max <- prediction$value + se

  if (model == 'loglin') {
    prediction$value <- prediction$value %>%
      format(scientific = FALSE) %>%
      as.integer()
    prediction$se_min <- prediction$se_min %>% as.integer()
    prediction$se_max <- prediction$se_max %>% as.integer()

    colnames(prediction) <- c('year', 'amountR', 'se_max', 'se_min')
  } else {
    prediction$se_max[which(prediction$se_max > 1)] <- 1
  }
  prediction$se_min[which(prediction$se_min < 0)] <- 0
  prediction$observations = NA

  total <- prediction

  if (preserve_measurements == TRUE) {
    # replace estimated data by observed data
    if (I_as_R == TRUE) {
      if (!'I' %in% colnames(df)) {
        df$I <- 0
      }
      df$value <- df$R / rowSums(df[, c('R', 'S', 'I')])
    } else {
      df$value <- df$R / rowSums(df[, c('R', 'S')])
    }
    measurements <- data.frame(year = df$year,
                               value = df$value,
                               se_min = NA,
                               se_max = NA,
                               observations = df$total,
                               stringsAsFactors = FALSE)
    colnames(measurements) <- colnames(prediction)

    total <- rbind(measurements,
                   prediction %>% filter(!year %in% df$year))
    if (model %in% c('binomial', 'binom', 'logit')) {
      total <- total %>% mutate(observed = ifelse(is.na(observations), NA, value),
                                estimated = prediction$value)
    }
  }

  if ("value" %in% colnames(total)) {
    total <- total %>%
      mutate(value = case_when(value > 1 ~ 1,
                               value < 0 ~ 0,
                               TRUE ~ value))
  }

  total %>% arrange(year)

}

#' @rdname resistance_predict
#' @export
rsi_predict <- resistance_predict
