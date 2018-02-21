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
#' @param antibiotics character vector with 1, 2 or 3 antibiotics that occur as column names in \code{tbl}, like \code{antibiotics = c("amox", "amcl")}
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
#'   rsi_df(antibiotics = c("amox", "metr"))
#' }
rsi_df <- function(tbl,
                   antibiotics,
                   interpretation = 'IR',
                   minimum = 30,
                   percent = FALSE,
                   info = TRUE,
                   warning = TRUE) {
  
  # we willen niet dat tbl$interpretation toevallig ook bestaat, dus:
  te_testen_uitslag_ab <- interpretation
  
  # validatie:
  if (min(grepl('^[a-z]{3,4}$', antibiotics)) == 0 &
      min(grepl('^rsi[1-2]$', antibiotics)) == 0) {
    for (i in 1:length(antibiotics)) {
      antibiotics[i] <- paste0('rsi', i)
    }
  }
  if (!grepl('^(S|SI|IS|I|IR|RI|R){1}$', te_testen_uitslag_ab)) {
    stop('Invalid `interpretation`; must be "S", "SI", "I", "IR", or "R".')
  }
  if ('is_ic' %in% colnames(tbl)) {
    if (n_distinct(tbl$is_ic) > 1) {
      warning('Dataset contains isolates from the Intensive Care. Exclude them from proper epidemiological analysis.')
    }
  }
  
  # transformeren wanneer gezocht wordt op verschillende uitslagen
  if (te_testen_uitslag_ab %in% c('SI', 'IS')) {
    for (i in 1:length(antibiotics)) {
      lijst <- tbl[, antibiotics[i]]
      if ('I' %in% lijst) {
        tbl[which(tbl[antibiotics[i]] == 'I'), ][antibiotics[i]] <- 'S'
      }
    }
    te_testen_uitslag_ab <- 'S'
  }
  if (te_testen_uitslag_ab %in% c('RI', 'IR')) {
    for (i in 1:length(antibiotics)) {
      lijst <- tbl[, antibiotics[i]]
      if ('I' %in% lijst) {
        tbl[which(tbl[antibiotics[i]] == 'I'), ][antibiotics[i]] <- 'R'
      }
    }
    te_testen_uitslag_ab <- 'R'
  }
  
  # breuk samenstellen
  if (length(antibiotics) == 1) {
    numerator <- tbl %>%
      filter(pull(., antibiotics[1]) == te_testen_uitslag_ab) %>%
      nrow()
    
    denominator <- tbl %>%
      filter(pull(., antibiotics[1]) %in% c("S", "I", "R")) %>%
      nrow()
    
  } else if (length(antibiotics) == 2) {
    numerator <- tbl %>%
      filter_at(vars(antibiotics[1], antibiotics[2]),
                any_vars(. == te_testen_uitslag_ab)) %>%
      filter_at(vars(antibiotics[1], antibiotics[2]),
                all_vars(. %in% c("S", "R", "I"))) %>%
      nrow()
    
    denominator <- tbl %>%
      filter_at(vars(antibiotics[1], antibiotics[2]),
                all_vars(. %in% c("S", "R", "I"))) %>%
      nrow()
    
  } else if (length(antibiotics) == 3) {
    numerator <- tbl %>%
      filter_at(vars(antibiotics[1], antibiotics[2], antibiotics[3]),
                any_vars(. == te_testen_uitslag_ab)) %>%
      filter_at(vars(antibiotics[1], antibiotics[2], antibiotics[3]),
                all_vars(. %in% c("S", "R", "I"))) %>%
      nrow()
    
    denominator <- tbl %>%
      filter_at(vars(antibiotics[1], antibiotics[2], antibiotics[3]),
                all_vars(. %in% c("S", "R", "I"))) %>%
      nrow()
    
  } else {
    stop('Maximum of 3 drugs allowed.')
  }
  
  # tekstdeel opbouwen
  if (info == TRUE) {
    cat('n =', denominator)
    info.txt1 <- percent(denominator / nrow(tbl))
    if (denominator == 0) {
      info.txt1 <- 'none'
    }
    info.txt2 <- gsub(',', ' and',
                      antibiotics %>%
                        abname(to = 'trivial',
                               tolower = TRUE) %>%
                        toString(), fixed = TRUE)
    info.txt2 <- gsub('rsi1 and rsi2', 'these two drugs', info.txt2, fixed = TRUE)
    info.txt2 <- gsub('rsi1', 'this drug', info.txt2, fixed = TRUE)
    cat(paste0(' (of ', nrow(tbl), ' in total; ', info.txt1, ' tested on ', info.txt2, ')\n'))
  }
  
  # rekenen en opmaken
  y <- numerator / denominator
  if (percent == TRUE) {
    y <- percent(y)
  }
  if (denominator < minimum) {
    if (warning == TRUE) {
      warning(paste0('TOO FEW ISOLATES OF ', toString(antibiotics), ' (n = ', denominator, ', n < ', minimum, '); NO RESULT.'))
    }
    y <- NA
  }
  
  # output
  y
}

#' Resistance of isolates
#'
#' This function can be used in \code{\link[dplyr]{summarise}}, see \emph{Examples}. CaBerekent het percentage S, SI, I, IR of R van een lijst isolaten. 
#' @param ab1,ab2 list with interpretations of an antibiotic
#' @inheritParams rsi_df
#' @details This function uses the \code{\link{rsi_df}} function internally.
#' @keywords rsi antibiotics isolate isolates
#' @return Double or, when \code{percent = TRUE}, a character.
#' @export
#' @examples
#' \dontrun{
#' tbl %>%
#'   group_by(year, hospital) %>%
#'   summarise(
#'     isolates = n(),
#'     cipro = rsi(cipr, percent = TRUE),
#'     amoxi = rsi(amox, percent = TRUE)
#'   )
#'
#' tbl %>%
#'   group_by(hospital) %>%
#'   summarise(cipr = rsi(cipr))
#'
#' rsi(isolates$amox)
#'
#' rsi(isolates$amcl, interpretation = "S")
#' }
rsi <- function(ab1, ab2 = NA, interpretation = 'IR', minimum = 30, percent = FALSE, info = FALSE, warning = FALSE) {
  functietekst <- as.character(match.call())
  # param 1 = functienaam
  # param 2 = ab1
  # param 3 = ab2
  ab1.naam <- functietekst[2]
  if (!grepl('^[a-z]{3,4}$', ab1.naam)) {
    ab1.naam <- 'rsi1'
  }
  ab2.naam <- functietekst[3]
  if (!grepl('^[a-z]{3,4}$', ab2.naam)) {
    ab2.naam <- 'rsi2'
  }
  
  tbl <- tibble(rsi1 = ab1, rsi2 = ab2)
  
  colnames(tbl) <- c(ab1.naam, ab2.naam)
  
  if (length(ab2) == 1) {
    return(rsi_df(tbl = tbl,
                  antibiotics = ab1.naam,
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
                  antibiotics = c(ab1.naam, ab2.naam),
                  interpretation = interpretation,
                  minimum = minimum,
                  percent = percent,
                  info = info,
                  warning = warning))
  }
}

#' Predict antimicrobial resistance
#'
#' Create a prediction model to predict antimicrobial resistance for the next years on statistical solid ground. Standard errors (SE) will be returned as columns \code{se_min} and \code{se_max}.
#' @param tbl table that contains columns \code{col_ab} and \code{col_date}
#' @param col_ab column name of \code{tbl} with antimicrobial interpretations (\code{R}, \code{I} and \code{S})
#' @param col_date column name of the date, will be used to calculate years
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
#' rsi_predict(tbl[which(first_isolate == TRUE & genus == "Haemophilus"),], "amcl")
#'   
#' # or with dplyr so you can actually read it:
#' tbl %>%
#'   filter(first_isolate == TRUE,
#'          genus == "Haemophilus") %>%
#'   rsi_predict("amcl")
#'
#' tbl %>%
#'   filter(first_isolate_weighted == TRUE,
#'          genus == "Haemophilus") %>%
#'   rsi_predict(col_ab = "amcl",
#'               year_max = 2050,
#'               year_every = 5)
#'
#' }
rsi_predict <- function(tbl,
                        col_ab,
                        col_date = 'ontvangstdatum',
                        year_max = as.integer(format(as.Date(Sys.Date()), '%Y')) + 15,
                        year_every = 1,
                        model = 'binomial',
                        I_as_R = TRUE,
                        preserve_measurements = TRUE,
                        info = TRUE) {
  
  if (I_as_R == TRUE) {
    tbl[, col_ab] <- gsub('I', 'R', tbl %>% pull(col_ab))
  }
  
  year <- function(x) {
    as.integer(format(as.Date(x), '%Y'))
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
