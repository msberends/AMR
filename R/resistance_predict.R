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
# Visit our website for more info: https://msberends.gitab.io/AMR.     #
# ==================================================================== #

#' Predict antimicrobial resistance
#'
#' Create a prediction model to predict antimicrobial resistance for the next years on statistical solid ground. Standard errors (SE) will be returned as columns \code{se_min} and \code{se_max}. See Examples for a real live example.
#' @inheritParams first_isolate
#' @inheritParams graphics::plot
#' @param col_ab column name of \code{tbl} with antimicrobial interpretations (\code{R}, \code{I} and \code{S})
#' @param col_date column name of the date, will be used to calculate years if this column doesn't consist of years already, defaults to the first column of with a date class
#' @param year_min lowest year to use in the prediction model, dafaults to the lowest year in \code{col_date}
#' @param year_max highest year to use in the prediction model, defaults to 10 years after today
#' @param year_every unit of sequence between lowest year found in the data and \code{year_max}
#' @param minimum minimal amount of available isolates per year to include. Years containing less observations will be estimated by the model.
#' @param model the statistical model of choice. Valid values are \code{"binomial"} (or \code{"binom"} or \code{"logit"}) or \code{"loglin"} (or \code{"poisson"}) or \code{"linear"} (or \code{"lin"}).
#' @param I_as_R a logical to indicate whether values \code{I} should be treated as \code{R}
#' @param preserve_measurements a logical to indicate whether predictions of years that are actually available in the data should be overwritten by the original data. The standard errors of those years will be \code{NA}.
#' @param info a logical to indicate whether textual analysis should be printed with the name and \code{\link{summary}} of the statistical model.
#' @return \code{data.frame} with extra class \code{"resistance_predict"} with columns:
#' \itemize{
#'   \item{\code{year}}
#'   \item{\code{value}, the same as \code{estimated} when \code{preserve_measurements = FALSE}, and a combination of \code{observed} and \code{estimated} otherwise}
#'   \item{\code{se_min}, the lower bound of the standard error with a minimum of \code{0} (so the standard error will never go below 0\%)}
#'   \item{\code{se_max} the upper bound of the standard error with a maximum of \code{1} (so the standard error will never go above 100\%)}
#'   \item{\code{observations}, the total number of available observations in that year, i.e. S + I + R}
#'   \item{\code{observed}, the original observed resistant percentages}
#'   \item{\code{estimated}, the estimated resistant percentages, calculated by the model}
#' }
#' @seealso The \code{\link{portion}} function to calculate resistance, \cr \code{\link{lm}} \code{\link{glm}}
#' @rdname resistance_predict
#' @export
#' @importFrom stats predict glm lm
#' @importFrom dplyr %>% pull mutate mutate_at n group_by_at summarise filter filter_at all_vars n_distinct arrange case_when n_groups
#' @inheritSection AMR Read more on our website!
#' @examples
#' x <- resistance_predict(septic_patients, col_ab = "amox", year_min = 2010)
#' plot(x)
#' ggplot_rsi_predict(x)
#'
#' # use dplyr so you can actually read it:
#' library(dplyr)
#' x <- septic_patients %>%
#'   filter_first_isolate() %>%
#'   filter(mo_genus(mo) == "Staphylococcus") %>%
#'   resistance_predict("peni")
#' plot(x)
#'
#'
#' # create nice plots with ggplot yourself
#' if (!require(ggplot2)) {
#'
#'   data <- septic_patients %>%
#'     filter(mo == as.mo("E. coli")) %>%
#'     resistance_predict(col_ab = "amox",
#'                        col_date = "date",
#'                        info = FALSE,
#'                        minimum = 15)
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
                               col_date = NULL,
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

  # -- date
  if (is.null(col_date)) {
    col_date <- search_type_in_df(tbl = tbl, type = "date")
  }
  if (is.null(col_date)) {
    stop("`col_date` must be set.", call. = FALSE)
  }

  if (!col_date %in% colnames(tbl)) {
    stop('Column ', col_date, ' not found.')
  }

  if (n_groups(tbl) > 1) {
    # no grouped tibbles please, mutate will throw errors
    tbl <- base::as.data.frame(tbl, stringsAsFactors = FALSE)
  }

  year <- function(x) {
    if (all(grepl('^[0-9]{4}$', x))) {
      x
    } else {
      as.integer(format(as.Date(x), '%Y'))
    }
  }

  df <- tbl %>%
    mutate_at(col_ab, as.rsi) %>%
    mutate_at(col_ab, droplevels) %>%
    mutate_at(col_ab, funs(
      if (I_as_R == TRUE) {
        gsub("I", "R", .)
      } else {
        gsub("I", "S", .)
      }
      )) %>%
    filter_at(col_ab, all_vars(!is.na(.))) %>%
    mutate(year = pull(., col_date) %>% year()) %>%
    group_by_at(c('year', col_ab)) %>%
    summarise(n())

  if (df %>% pull(col_ab) %>% n_distinct(na.rm = TRUE) < 2) {
    stop("No variety in antimicrobial interpretations - all isolates are '",
         df %>% pull(col_ab) %>% unique(), "'.",
         call. = FALSE)
  }

  colnames(df) <- c('year', 'antibiotic', 'observations')
  df <- df %>%
    filter(!is.na(antibiotic)) %>%
    tidyr::spread(antibiotic, observations, fill = 0) %>%
    filter((R + S) >= minimum)
  df_matrix <- df %>%
    ungroup() %>%
    select(R, S) %>%
    as.matrix()

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
    year_max <- year(Sys.Date()) + 10
  }

  years <- list(year = seq(from = year_min, to = year_max, by = year_every))

  if (model %in% c('binomial', 'binom', 'logit')) {
    model <- "binomial"
    model_lm <- with(df, glm(df_matrix ~ year, family = binomial))
    if (info == TRUE) {
      cat('\nLogistic regression model (logit) with binomial distribution')
      cat('\n------------------------------------------------------------\n')
      print(summary(model_lm))
    }

    predictmodel <- predict(model_lm, newdata = years, type = "response", se.fit = TRUE)
    prediction <- predictmodel$fit
    se <- predictmodel$se.fit

  } else if (model %in% c('loglin', 'poisson')) {
    model <- "poisson"
    model_lm <- with(df, glm(R ~ year, family = poisson))
    if (info == TRUE) {
      cat('\nLog-linear regression model (loglin) with poisson distribution')
      cat('\n--------------------------------------------------------------\n')
      print(summary(model_lm))
    }

    predictmodel <- predict(model_lm, newdata = years, type = "response", se.fit = TRUE)
    prediction <- predictmodel$fit
    se <- predictmodel$se.fit

  } else if (model %in% c('lin', 'linear')) {
    model <- "linear"
    model_lm <- with(df, lm((R / (R + S)) ~ year))
    if (info == TRUE) {
      cat('\nLinear regression model')
      cat('\n-----------------------\n')
      print(summary(model_lm))
    }

    predictmodel <- predict(model_lm, newdata = years, se.fit = TRUE)
    prediction <- predictmodel$fit
    se <- predictmodel$se.fit

  } else {
    stop('No valid model selected.')
  }

  # prepare the output dataframe
  df_prediction <- data.frame(year = unlist(years),
                              value = prediction,
                              stringsAsFactors = FALSE) %>%

    mutate(se_min = value - se,
           se_max = value + se)

  if (model == 'poisson') {
    df_prediction <- df_prediction %>%
      mutate(value = value %>%
               format(scientific = FALSE) %>%
               as.integer(),
             se_min = as.integer(se_min),
             se_max = as.integer(se_max))
  } else {
    df_prediction <- df_prediction %>%
      # se_max not above 1
      mutate(se_max = ifelse(se_max > 1, 1, se_max))
  }
  df_prediction <- df_prediction %>%
    # se_min not below 0
    mutate(se_min = ifelse(se_min < 0, 0, se_min))

  df_observations <- df %>%
    ungroup() %>%
    transmute(year,
              observations = R + S,
              observed = R / (R + S))
  df_prediction <- df_prediction %>%
    left_join(df_observations, by = "year") %>%
    mutate(estimated = value)

  if (preserve_measurements == TRUE) {
    # replace estimated data by observed data
    df_prediction <- df_prediction %>%
      mutate(value = ifelse(!is.na(observed), observed, value),
             se_min = ifelse(!is.na(observed), NA, se_min),
             se_max = ifelse(!is.na(observed), NA, se_max))
  }

  df_prediction <- df_prediction %>%
    mutate(value = case_when(value > 1 ~ 1,
                             value < 0 ~ 0,
                             TRUE ~ value)) %>%
    arrange(year)

  structure(
    .Data = df_prediction,
    class = c("resistance_predict", "data.frame"),
    I_as_R = I_as_R,
    model_title = model,
    model = model_lm,
    ab = col_ab
  )
}

#' @rdname resistance_predict
#' @export
rsi_predict <- resistance_predict

#' @exportMethod plot.mic
#' @export
#' @importFrom dplyr %>% group_by summarise
#' @importFrom graphics plot axis arrows
#' @rdname resistance_predict
plot.resistance_predict <- function(x, main = paste("Resistance prediction of", attributes(x)$ab), ...) {
  if (attributes(x)$I_as_R == TRUE) {
    ylab <- "%IR"
  } else {
    ylab <- "%R"
  }
  plot(x = x$year,
       y = x$value,
       ylim = c(0, 1),
       yaxt = "n", # no y labels
       pch = 19, # closed dots
       ylab = paste0("Percentage (", ylab, ")"),
       xlab = "Year",
       main = main,
       sub = paste0("(model: ", attributes(x)$model_title, ")"))

  axis(side = 2, at = seq(0, 1, 0.1), labels = paste0(0:10 * 10, "%"))

  # arrows hack: https://stackoverflow.com/a/22037078/4575331
  arrows(x0 = x$year,
         y0 = x$se_min,
         x1 = x$year,
         y1 = x$se_max, length = 0.05, angle = 90, code = 3)
}

#' @rdname resistance_predict
#' @export
ggplot_rsi_predict <- function(x, main = paste("Resistance prediction of", attributes(x)$ab), ...) {

  if (!"resistance_predict" %in% class(x)) {
    stop("`x` must be a resistance prediction model created with resistance_predict().")
  }

  if (attributes(x)$I_as_R == TRUE) {
    ylab <- "%IR"
  } else {
    ylab <- "%R"
  }
  suppressWarnings(
    ggplot2::ggplot(x, ggplot2::aes(x = year, y = value)) +
      ggplot2::geom_col() +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = se_min, ymax = se_max)) +
      scale_y_percent() +
      labs(title = main,
           y = paste0("Percentage (", ylab, ")"),
           x = "Year",
           caption = paste0("(model: ", attributes(x)$model_title, ")"))
  )
}
