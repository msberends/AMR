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

#' Resistance of isolates in data.frame
#'
#' \strong{NOTE: use \code{\link{rsi}} in dplyr functions like \code{\link[dplyr]{summarise}}.} \cr Calculate the percentage of S, SI, I, IR or R of a \code{data.frame} containing isolates.
#' @param tbl \code{data.frame} containing columns with antibiotic interpretations.
#' @param ab character vector with 1, 2 or 3 antibiotics that occur as column names in \code{tbl}, like \code{ab = c("amox", "amcl")}
#' @param interpretation antimicrobial interpretation of which the portion must be calculated. Valid values are \code{"S"}, \code{"SI"}, \code{"I"}, \code{"IR"} or \code{"R"}.
#' @param minimum minimal amount of available isolates. Any number lower than \code{minimum} will return \code{NA} with a warning (when \code{warning = TRUE}).
#' @param percent return output as percent (text), will else (at default) be a double
#' @param info calculate the amount of available isolates and print it, like \code{n = 423}
#' @param warning show a warning when the available amount of isolates is below \code{minimum}
#' @details Remember that you should filter your table to let it contain \strong{only first isolates}!
#' @keywords rsi antibiotics isolate isolates
#' @return Double or, when \code{percent = TRUE}, a character.
#' @export
#' @importFrom dplyr %>% n_distinct filter filter_at pull vars all_vars any_vars
#' @seealso \code{\link{rsi}} for the function that can be used with \code{\link[dplyr]{summarise}} directly.
#' @examples
#' \dontrun{
#' rsi_df(tbl_with_bloodcultures, 'amcl')
#'
#' rsi_df(tbl_with_bloodcultures, c('amcl', 'gent'), interpretation = 'IR')
#'
#' library(dplyr)
#' # calculate current empiric therapy of Helicobacter gastritis:
#' my_table %>%
#'   filter(first_isolate == TRUE,
#'          genus == "Helicobacter") %>%
#'   rsi_df(ab = c("amox", "metr"))
#' }
rsi_df <- function(tbl,
                   ab,
                   interpretation = 'IR',
                   minimum = 30,
                   percent = FALSE,
                   info = TRUE,
                   warning = TRUE) {

  # in case tbl$interpretation already exists:
  interpretations_to_check <- paste(interpretation, collapse = "")

  # validate:
  if (min(grepl('^[a-z]{3,4}$', ab)) == 0 &
      min(grepl('^rsi[1-2]$', ab)) == 0) {
    for (i in 1:length(ab)) {
      ab[i] <- paste0('rsi', i)
    }
  }
  if (!grepl('^(S|SI|IS|I|IR|RI|R){1}$', interpretations_to_check)) {
    stop('Invalid `interpretation`; must be "S", "SI", "I", "IR", or "R".')
  }
  if ('is_ic' %in% colnames(tbl)) {
    if (n_distinct(tbl$is_ic) > 1 & warning == TRUE) {
      warning('Dataset contains isolates from the Intensive Care. Exclude them from proper epidemiological analysis.')
    }
  }

  # transform when checking for different results
  if (interpretations_to_check %in% c('SI', 'IS')) {
    for (i in 1:length(ab)) {
      lijst <- tbl[, ab[i]]
      if ('I' %in% lijst) {
        tbl[which(tbl[ab[i]] == 'I'), ][ab[i]] <- 'S'
      }
    }
    interpretations_to_check <- 'S'
  }
  if (interpretations_to_check %in% c('RI', 'IR')) {
    for (i in 1:length(ab)) {
      lijst <- tbl[, ab[i]]
      if ('I' %in% lijst) {
        tbl[which(tbl[ab[i]] == 'I'), ][ab[i]] <- 'R'
      }
    }
    interpretations_to_check <- 'R'
  }

  # get fraction
  if (length(ab) == 1) {
    numerator <- tbl %>%
      filter(pull(., ab[1]) == interpretations_to_check) %>%
      nrow()

    denominator <- tbl %>%
      filter(pull(., ab[1]) %in% c("S", "I", "R")) %>%
      nrow()

  } else if (length(ab) == 2) {
    numerator <- tbl %>%
      filter_at(vars(ab[1], ab[2]),
                any_vars(. == interpretations_to_check)) %>%
      filter_at(vars(ab[1], ab[2]),
                all_vars(. %in% c("S", "R", "I"))) %>%
      nrow()

    denominator <- tbl %>%
      filter_at(vars(ab[1], ab[2]),
                all_vars(. %in% c("S", "R", "I"))) %>%
      nrow()

  } else if (length(ab) == 3) {
    numerator <- tbl %>%
      filter_at(vars(ab[1], ab[2], ab[3]),
                any_vars(. == interpretations_to_check)) %>%
      filter_at(vars(ab[1], ab[2], ab[3]),
                all_vars(. %in% c("S", "R", "I"))) %>%
      nrow()

    denominator <- tbl %>%
      filter_at(vars(ab[1], ab[2], ab[3]),
                all_vars(. %in% c("S", "R", "I"))) %>%
      nrow()

  } else {
    stop('Maximum of 3 drugs allowed.')
  }

  # build text part
  if (info == TRUE) {
    cat('n =', denominator)
    info.txt1 <- percent(denominator / nrow(tbl))
    if (denominator == 0) {
      info.txt1 <- 'none'
    }
    info.txt2 <- gsub(',', ' and',
                      ab %>%
                        abname(tolower = TRUE) %>%
                        toString(), fixed = TRUE)
    info.txt2 <- gsub('rsi1 and rsi2', 'these two drugs', info.txt2, fixed = TRUE)
    info.txt2 <- gsub('rsi1', 'this drug', info.txt2, fixed = TRUE)
    cat(paste0(' (of ', nrow(tbl), ' in total; ', info.txt1, ' tested on ', info.txt2, ')\n'))
  }

  # calculate and format
  y <- numerator / denominator
  if (percent == TRUE) {
    y <- percent(y)
  }
  if (denominator < minimum) {
    if (warning == TRUE) {
      warning(paste0('TOO FEW ISOLATES OF ', toString(ab), ' (n = ', denominator, ', n < ', minimum, '); NO RESULT.'))
    }
    y <- NA
  }

  # output
  y
}

#' Resistance of isolates
#'
#' This function can be used in \code{dplyr}s \code{\link[dplyr]{summarise}}, see \emph{Examples}. Calculate the percentage S, SI, I, IR or R of a vector of isolates.
#' @param ab1,ab2 list with interpretations of an antibiotic
#' @inheritParams rsi_df
#' @details This function uses the \code{\link{rsi_df}} function internally.
#' @keywords rsi antibiotics isolate isolates
#' @return Double or, when \code{percent = TRUE}, a character.
#' @export
#' @examples
#' \dontrun{
#' tbl %>%
#'   group_by(hospital) %>%
#'   summarise(cipr = rsi(cipr))
#'
#' tbl %>%
#'   group_by(year, hospital) %>%
#'   summarise(
#'     isolates = n(),
#'     cipro = rsi(cipr %>% as.rsi(), percent = TRUE),
#'     amoxi = rsi(amox %>% as.rsi(), percent = TRUE))
#'
#' rsi(as.rsi(isolates$amox))
#'
#' rsi(as.rsi(isolates$amcl), interpretation = "S")
#' }
rsi <- function(ab1, ab2 = NA, interpretation = 'IR', minimum = 30, percent = FALSE, info = FALSE, warning = FALSE) {
  ab1.name <- deparse(substitute(ab1))
  if (ab1.name %like% '.[$].') {
    ab1.name <- unlist(strsplit(ab1.name, "$", fixed = TRUE))
    ab1.name <- ab1.name[length(ab1.name)]
  }
  if (!ab1.name %like% '^[a-z]{3,4}$') {
    ab1.name <- 'rsi1'
  }
  ab2.name <- deparse(substitute(ab2))
  if (ab2.name %like% '.[$].') {
    ab2.name <- unlist(strsplit(ab2.name, "$", fixed = TRUE))
    ab2.name <- ab2.name[length(ab2.name)]
  }
  if (!ab2.name %like% '^[a-z]{3,4}$') {
    ab2.name <- 'rsi2'
  }

  interpretation <- paste(interpretation, collapse = "")

  tbl <- tibble(rsi1 = ab1, rsi2 = ab2)
  colnames(tbl) <- c(ab1.name, ab2.name)

  if (length(ab2) == 1) {
    return(rsi_df(tbl = tbl,
                  ab = ab1.name,
                  interpretation = interpretation,
                  minimum = minimum,
                  percent = percent,
                  info = info,
                  warning = warning))
  } else {
    if (length(ab1) != length(ab2)) {
      stop('`ab1` (n = ', length(ab1), ') and `ab2` (n = ', length(ab2), ') must be of same length.', call. = FALSE)
    }
    if (interpretation != 'S') {
      warning('`interpretation` is not set to S, albeit analysing a combination therapy.')
    }
    return(rsi_df(tbl = tbl,
                  ab = c(ab1.name, ab2.name),
                  interpretation = interpretation,
                  minimum = minimum,
                  percent = percent,
                  info = info,
                  warning = warning))
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
#' @seealso \code{\link{lm}} \cr \code{\link{glm}}
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
rsi_predict <- function(tbl,
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
