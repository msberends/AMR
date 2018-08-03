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
#' These functions can be used to calculate the (co-)resistance of microbial isolates (i.e. percentage S, SI, I, IR or R). All functions can be used in \code{dplyr}s \code{\link[dplyr]{summarise}} and support grouped variables, see \emph{Examples}. \cr\cr
#' \code{rsi_R} and \code{rsi_IR} can be used to calculate resistance, \code{rsi_S} and \code{rsi_SI} can be used to calculate susceptibility.\cr
#' \code{rsi_n} counts all cases where antimicrobial interpretations are available.
#' @param ab1 vector of antibiotic interpretations, they will be transformed internally with \code{\link{as.rsi}}
#' @param ab2 like \code{ab}, a vector of antibiotic interpretations. Use this to calculate (the lack of) co-resistance: the probability where one of two drugs have a susceptible result. See Examples.
#' @param include_I logical to indicate whether antimicrobial interpretations of "I" should be included
#' @param minimum minimal amount of available isolates. Any number lower than \code{minimum} will return \code{NA}. The default number of \code{30} isolates is advised by the CLSI as best practice, see Source.
#' @param as_percent logical to indicate whether the output must be returned as percent (text), will else be a double
#' @details \strong{Remember that you should filter your table to let it contain only first isolates!} Use \code{\link{first_isolate}} to determine them in your data set.
#'
#' The functions \code{resistance} and \code{susceptibility} are wrappers around \code{rsi_IR} and \code{rsi_S}, respectively. All functions use hybrid evaluation (i.e. using C++), which makes these functions 20-30 times faster than the old \code{\link{rsi}} function. This latter function is still available for backwards compatibility but is deprecated.
#' \if{html}{
#'   \cr\cr
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
#' @source \strong{M39 Analysis and Presentation of Cumulative Antimicrobial Susceptibility Test Data, 4th Edition}, 2014, \emph{Clinical and Laboratory Standards Institute (CLSI)}. \url{https://clsi.org/standards/products/microbiology/documents/m39/}.
#' @keywords resistance susceptibility rsi_df rsi antibiotics isolate isolates
#' @return Double or, when \code{as_percent = TRUE}, a character.
#' @rdname rsi_IR
#' @name rsi_IR
#' @export
#' @examples
#' # Calculate resistance
#' rsi_R(septic_patients$amox)
#' rsi_IR(septic_patients$amox)
#'
#' # Or susceptibility
#' rsi_S(septic_patients$amox)
#' rsi_SI(septic_patients$amox)
#'
#' # Since n_rsi counts available isolates (and is used as denominator),
#' # you can calculate back to e.g. count resistant isolates:
#' rsi_IR(septic_patients$amox) * n_rsi(septic_patients$amox)
#'
#' library(dplyr)
#' septic_patients %>%
#'   group_by(hospital_id) %>%
#'   summarise(p = rsi_S(cipr),
#'             n = rsi_n(cipr)) # n_rsi works like n_distinct in dplyr
#'
#' septic_patients %>%
#'   group_by(hospital_id) %>%
#'   summarise(R = rsi_R(cipr, as_percent = TRUE),
#'             I = rsi_I(cipr, as_percent = TRUE),
#'             S = rsi_S(cipr, as_percent = TRUE),
#'             n = rsi_n(cipr), # also: n_rsi, works like n_distinct in dplyr
#'             total = n()) # this is the length, NOT the amount of tested isolates
#'
#' # Calculate co-resistance between amoxicillin/clav acid and gentamicin,
#' # so we can see that combination therapy does a lot more than mono therapy:
#' rsi_S(septic_patients$amcl) # S = 67.3%
#' rsi_n(septic_patients$amcl) # n = 1570
#'
#' rsi_S(septic_patients$gent) # S = 74.0%
#' rsi_n(septic_patients$gent) # n = 1842
#'
#' with(septic_patients,
#'      rsi_S(amcl, gent))     # S = 92.1%
#' with(septic_patients,       # n = 1504
#'      rsi_n(amcl, gent))
#'
#' septic_patients %>%
#'   group_by(hospital_id) %>%
#'   summarise(cipro_p = rsi_S(cipr, as_percent = TRUE),
#'             cipro_n = rsi_n(cipr),
#'             genta_p = rsi_S(gent, as_percent = TRUE),
#'             genta_n = rsi_n(gent),
#'             combination_p = rsi_S(cipr, gent, as_percent = TRUE),
#'             combination_n = rsi_n(cipr, gent))
#'
#' \dontrun{
#'
#' # calculate current empiric combination therapy of Helicobacter gastritis:
#' my_table %>%
#'   filter(first_isolate == TRUE,
#'          genus == "Helicobacter") %>%
#'   summarise(p = rsi_S(amox, metr),  # amoxicillin with metronidazole
#'             n = rsi_n(amox, metr))
#' }
rsi_R <- function(ab1,
                  minimum = 30,
                  as_percent = FALSE) {
  resistance(ab1 = ab1,
             include_I = FALSE,
             minimum = minimum,
             as_percent = as_percent)
}

#' @rdname rsi_IR
#' @export
rsi_IR <- function(ab1,
                   minimum = 30,
                   as_percent = FALSE) {
  resistance(ab1 = ab1,
             include_I = TRUE,
             minimum = minimum,
             as_percent = as_percent)
}

#' @rdname rsi_IR
#' @export
rsi_I <- function(ab1,
                  minimum = 30,
                  as_percent = FALSE) {
  intermediate(ab1 = ab1,
               minimum = minimum,
               as_percent = as_percent)
}

#' @rdname rsi_IR
#' @export
rsi_SI <- function(ab1,
                   ab2 = NULL,
                   minimum = 30,
                   as_percent = FALSE) {
  susceptibility(ab1 = ab1,
                 ab2 = ab2,
                 include_I = TRUE,
                 minimum = minimum,
                 as_percent = as_percent)
}

#' @rdname rsi_IR
#' @export
rsi_S <- function(ab1,
                  ab2 = NULL,
                  minimum = 30,
                  as_percent = FALSE) {
  susceptibility(ab1 = ab1,
                 ab2 = ab2,
                 include_I = FALSE,
                 minimum = minimum,
                 as_percent = as_percent)
}

#' @rdname rsi_IR
#' @export
resistance <- function(ab1,
                       include_I = TRUE,
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

  # ab_name <- deparse(substitute(ab))

  if (!is.rsi(ab1)) {
    x <- as.rsi(ab1)
    warning("Increase speed by transforming to class `rsi` on beforehand: df %>% mutate_at(vars(col10:col20), as.rsi)")
  } else {
    x <- ab1
  }
  total <- length(x) - sum(is.na(x)) # faster than C++
  if (total < minimum) {
    # warning("Too few isolates available for ", ab_name, ": ", total, " < ", minimum, "; returning NA.", call. = FALSE)
    return(NA)
  }
  found <- .Call(`_AMR_rsi_calc_R`, x, include_I)

  if (as_percent == TRUE) {
    percent(found / total, force_zero = TRUE)
  } else {
    found / total
  }
}

#' @rdname rsi_IR
#' @export
intermediate <- function(ab1,
                         minimum = 30,
                         as_percent = FALSE) {

  if (NCOL(ab1) > 1) {
    stop('`ab1` must be a vector of antimicrobial interpretations', call. = FALSE)
  }
  if (!is.numeric(minimum)) {
    stop('`minimum` must be numeric', call. = FALSE)
  }
  if (!is.logical(as_percent)) {
    stop('`as_percent` must be logical', call. = FALSE)
  }

  # ab_name <- deparse(substitute(ab))

  if (!is.rsi(ab1)) {
    x <- as.rsi(ab1)
    warning("Increase speed by transforming to class `rsi` on beforehand: df %>% mutate_at(vars(col10:col20), as.rsi)")
  } else {
    x <- ab1
  }
  total <- length(x) - sum(is.na(x)) # faster than C++
  if (total < minimum) {
    # warning("Too few isolates available for ", ab_name, ": ", total, " < ", minimum, "; returning NA.", call. = FALSE)
    return(NA)
  }
  found <- .Call(`_AMR_rsi_calc_I`, x)

  if (as_percent == TRUE) {
    percent(found / total, force_zero = TRUE)
  } else {
    found / total
  }
}

#' @rdname rsi_IR
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
    # ab_name <- paste(deparse(substitute(ab1)), "and", deparse(substitute(ab2)))
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
    # ab_name <- deparse(substitute(ab1))
  }
  total <- length(x) - sum(is.na(x))
  if (total < minimum) {
    # warning("Too few isolates available for ", ab_name, ": ", total, " < ", minimum, "; returning NA.", call. = FALSE)
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

#' @rdname rsi_IR
#' @export
rsi_n <- function(ab1, ab2 = NULL) {
  if (NCOL(ab1) > 1) {
    stop('`ab` must be a vector of antimicrobial interpretations', call. = FALSE)
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

#' @rdname rsi_IR
#' @export
n_rsi <- rsi_n

#' @inherit resistance
#' @description This function is deprecated. Use \code{\link{rsi_IR}} instead.
#' @param info calculate the amount of available isolates and print it, like \code{n = 423}
#' @param warning show a warning when the available amount of isolates is below \code{minimum}
#' @param interpretation antimicrobial interpretation
#' @export
rsi <- function(ab1,
                ab2 = NA,
                interpretation = 'IR',
                minimum = 30,
                as_percent = FALSE,
                info = FALSE,
                warning = TRUE) {

  .Deprecated()

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
#'   \item{\code{value}, the same as \code{estimated} when \code{preserve_measurements = FALSE}, and a combination of \code{observed} and \code{estimated} otherwise}
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
