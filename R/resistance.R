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
#' @param info \emph{DEPRECATED} calculate the amount of available isolates and print it, like \code{n = 423}
#' @param warning \emph{DEPRECATED} show a warning when the available amount of isolates is below \code{minimum}
#' @details \strong{Remember that you should filter your table to let it contain only first isolates!} Use \code{\link{first_isolate}} to determine them in your data set.
#'
#' The functions \code{resistance}, \code{susceptibility} and \code{n_rsi} calculate using hybrid evaluation (i.e. using C++), which makes these functions 25-30 times faster than the old \code{rsi} function. This function is still available for backwards compatibility but is deprecated.
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
#'
#'
#' # How fast is this hybrid evaluation in C++ compared to R?
#' # In other words: how is the speed improvement of the new `resistance` compared to old `rsi`?
#'
#' library(microbenchmark)
#' df <- septic_patients %>% group_by(hospital_id, bactid) # 317 groups with sizes 1 to 167
#'
#' microbenchmark(old_IR = df %>% summarise(p = rsi(amox, minimum = 0, interpretation = "IR")),
#'                new_IR = df %>% summarise(p = resistance(amox, minimum = 0)),
#'                old_S = df %>% summarise(p = rsi(amox, minimum = 0, interpretation = "S")),
#'                new_S = df %>% summarise(p = susceptibility(amox, minimum = 0)),
#'                times = 5,
#'                unit = "s")
#'
#' # Unit: seconds
#' #    expr          min         lq       mean     median         uq        max neval
#' #  old_IR   1.95600230 1.96096857 1.97981537 1.96823318 2.00645711 2.00741568     5
#' #  new_IR   0.06872808 0.06984932 0.07162866 0.06987306 0.07050094 0.07919192     5
#' #   old_S   1.68893579 1.69024888 1.72461867 1.69785934 1.70428796 1.84176137     5
#' #   new_S   0.06737037 0.06838167 0.07431906 0.07745364 0.07827224 0.08011738     5
#'
#' # The old function took roughly 2 seconds, the new ones take 0.07 seconds.
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

  if (!is.rsi(ab)) {
    x <- as.rsi(ab)
    warning("Increase speed by transforming to class `rsi` on beforehand: df %>% mutate_at(vars(col10:col20), as.rsi)")
  } else {
    x <- ab
  }
  total <- length(x) - sum(is.na(x)) # faster than C++
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

  print_warning <- FALSE
  if (!is.rsi(ab1)) {
    ab1 <- as.rsi(ab1)
    print_warning <- TRUE
  }
  if (!is.null(ab2)) {
    if (NCOL(ab2) > 1) {
      stop('`ab2` must be a vector of antimicrobial interpretations', call. = FALSE)
    }
    if (!is.rsi(ab2)) {
      ab2 <- as.rsi(ab2)
      print_warning <- TRUE
    }
    x <- apply(X = data.frame(ab1 = as.integer(ab1),
                              ab2 = as.integer(ab2)),
               MARGIN = 1,
               FUN = min)
  } else {
    x <- ab1
  }
  total <- length(x) - sum(is.na(x))
  if (total < minimum) {
    return(NA)
  }
  found <- .Call(`_AMR_rsi_calc_S`, x, include_I)

  if (print_warning == TRUE) {
    warning("Increase speed by transforming to class `rsi` on beforehand: df %>% mutate_at(vars(col10:col20), as.rsi)")
  }

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
  if (!is.rsi(ab1)) {
    ab1 <- as.rsi(ab1)
  }
  if (!is.null(ab2)) {
    if (NCOL(ab2) > 1) {
      stop('`ab2` must be a vector of antimicrobial interpretations', call. = FALSE)
    }
    if (!is.rsi(ab2)) {
      ab2 <- as.rsi(ab2)
    }
    sum(!is.na(ab1) & !is.na(ab2))
  } else {
    sum(!is.na(ab1))
  }
}

#' @rdname resistance
#' @export
rsi <- function(ab1,
                ab2 = NA,
                interpretation = 'IR',
                minimum = 30,
                as_percent = FALSE,
                info = FALSE,
                warning = TRUE) {
  ab1.name <- deparse(substitute(ab1))
  if (ab1.name %like% '.[$].') {
    ab1.name <- unlist(strsplit(ab1.name, "$", fixed = TRUE))
    ab1.name <- ab1.name[length(ab1.name)]
  }
  if (!ab1.name %like% '^[a-z]{3,4}$') {
    ab1.name <- 'rsi1'
  }
  if (length(ab1) == 1 & is.character(ab1)) {
    stop('`ab1` must be a vector of antibiotic interpretations.',
         '\n  Try rsi(', ab1, ', ...) instead of rsi("', ab1, '", ...)', call. = FALSE)
  }
  ab2.name <- deparse(substitute(ab2))
  if (ab2.name %like% '.[$].') {
    ab2.name <- unlist(strsplit(ab2.name, "$", fixed = TRUE))
    ab2.name <- ab2.name[length(ab2.name)]
  }
  if (!ab2.name %like% '^[a-z]{3,4}$') {
    ab2.name <- 'rsi2'
  }
  if (length(ab2) == 1 & is.character(ab2)) {
    stop('`ab2` must be a vector of antibiotic interpretations.',
         '\n  Try rsi(', ab2, ', ...) instead of rsi("', ab2, '", ...)', call. = FALSE)
  }

  interpretation <- paste(interpretation, collapse = "")

  ab1 <- as.rsi(ab1)
  ab2 <- as.rsi(ab2)

  tbl <- tibble(rsi1 = ab1, rsi2 = ab2)
  colnames(tbl) <- c(ab1.name, ab2.name)

  if (length(ab2) == 1) {
    r <- rsi_df(tbl = tbl,
                ab = ab1.name,
                interpretation = interpretation,
                minimum = minimum,
                as_percent = FALSE,
                info = info,
                warning = warning)
  } else {
    if (length(ab1) != length(ab2)) {
      stop('`ab1` (n = ', length(ab1), ') and `ab2` (n = ', length(ab2), ') must be of same length.', call. = FALSE)
    }
    if (!interpretation %in% c('S', 'IS', 'SI')) {
      warning('`interpretation` not set to S or I/S, albeit analysing a combination therapy.', call. = FALSE)
    }
    r <- rsi_df(tbl = tbl,
                ab = c(ab1.name, ab2.name),
                interpretation = interpretation,
                minimum = minimum,
                as_percent = FALSE,
                info = info,
                warning = warning)
  }
  if (as_percent == TRUE) {
    percent(r, force_zero = TRUE)
  } else {
    r
  }
}

#' @importFrom dplyr %>% filter_at vars any_vars all_vars
#' @noRd
rsi_df <- function(tbl,
                   ab,
                   interpretation = 'IR',
                   minimum = 30,
                   as_percent = FALSE,
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
      tbl[which(tbl[, ab[i]] == 'I'), ab[i]] <- 'S'
    }
    interpretations_to_check <- 'S'
  }
  if (interpretations_to_check %in% c('RI', 'IR')) {
    for (i in 1:length(ab)) {
      tbl[which(tbl[, ab[i]] == 'I'), ab[i]] <- 'R'
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
    if (interpretations_to_check != 'S') {
      warning('`interpretation` not set to S or I/S, albeit analysing a combination therapy.', call. = FALSE)
    }
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

  } else {
    stop('Maximum of 2 drugs allowed.')
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
  if (as_percent == TRUE) {
    y <- percent(y, force_zero = TRUE)
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
#'   \item{\code{resistance}, the same as \code{estimated} when \code{preserve_measurements = FALSE}, and a combination of \code{observed} and \code{estimated} otherwise}
#'   \item{\code{se_min}, the lower bound of the standard error with a minimum of \code{0}}
#'   \item{\code{se_max} the upper bound of the standard error with a maximum of \code{1}}
#'   \item{\code{observations}, the total number of observations, i.e. S + I + R}
#'   \item{\code{observed}, the original observed values}
#'   \item{\code{estimated}, the estimated values, calculated by the model}
#' }
#' @seealso \code{\link{resistance}} \cr \code{\link{lm}} \code{\link{glm}}
#' @rdname resistance_predict
#' @export
#' @importFrom stats predict glm lm
#' @importFrom dplyr %>% pull mutate group_by_at summarise filter n_distinct arrange
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
#'     filter(bactid == "ESCCOL") %>%
#'     resistance_predict(col_ab = "amox",
#'                       col_date = "date",
#'                       info = FALSE,
#'                       minimum = 15)
#'
#'   ggplot(data,
#'          aes(x = year)) +
#'     geom_col(aes(y = resistance),
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
  prediction <- data.frame(year = years_predict, resistance = prediction, stringsAsFactors = FALSE)

  prediction$se_min <- prediction$resistance - se
  prediction$se_max <- prediction$resistance + se

  if (model == 'loglin') {
    prediction$resistance <- prediction$resistance %>%
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
      df$resistance <- df$R / rowSums(df[, c('R', 'S', 'I')])
    } else {
      df$resistance <- df$R / rowSums(df[, c('R', 'S')])
    }
    measurements <- data.frame(year = df$year,
                               resistance = df$resistance,
                               se_min = NA,
                               se_max = NA,
                               observations = df$total,
                               stringsAsFactors = FALSE)
    colnames(measurements) <- colnames(prediction)

    total <- rbind(measurements,
                   prediction %>% filter(!year %in% df$year))
    if (model %in% c('binomial', 'binom', 'logit')) {
      total <- total %>% mutate(observed = ifelse(is.na(observations), NA, resistance),
                                estimated = prediction$resistance)
    }
  }

  try(
    total$resistance[which(total$resistance > 1)] <- 1,
    total$resistance[which(total$resistance < 0)] <- 0,
    silent = TRUE
  )
  total %>% arrange(year)

}

#' @rdname resistance_predict
#' @export
rsi_predict <- resistance_predict
