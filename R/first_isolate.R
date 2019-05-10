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

#' Determine first (weighted) isolates
#'
#' Determine first (weighted) isolates of all microorganisms of every patient per episode and (if needed) per specimen type.
#' @param tbl a \code{data.frame} containing isolates.
#' @param col_date column name of the result date (or date that is was received on the lab), defaults to the first column of with a date class
#' @param col_patient_id column name of the unique IDs of the patients, defaults to the first column that starts with 'patient' or 'patid' (case insensitive)
#' @param col_mo column name of the unique IDs of the microorganisms (see \code{\link{mo}}), defaults to the first column of class \code{mo}. Values will be coerced using \code{\link{as.mo}}.
#' @param col_testcode column name of the test codes. Use \code{col_testcode = NULL} to \strong{not} exclude certain test codes (like test codes for screening). In that case \code{testcodes_exclude} will be ignored.
#' @param col_specimen column name of the specimen type or group
#' @param col_icu column name of the logicals (\code{TRUE}/\code{FALSE}) whether a ward or department is an Intensive Care Unit (ICU)
#' @param col_keyantibiotics column name of the key antibiotics to determine first \emph{weighted} isolates, see \code{\link{key_antibiotics}}. Defaults to the first column that starts with 'key' followed by 'ab' or 'antibiotics' (case insensitive). Use \code{col_keyantibiotics = FALSE} to prevent this.
#' @param episode_days episode in days after which a genus/species combination will be determined as 'first isolate' again
#' @param testcodes_exclude character vector with test codes that should be excluded (case-insensitive)
#' @param icu_exclude logical whether ICU isolates should be excluded (rows with value \code{TRUE} in column \code{col_icu})
#' @param specimen_group value in column \code{col_specimen} to filter on
#' @param type type to determine weighed isolates; can be \code{"keyantibiotics"} or \code{"points"}, see Details
#' @param ignore_I logical to determine whether antibiotic interpretations with \code{"I"} will be ignored when \code{type = "keyantibiotics"}, see Details
#' @param points_threshold points until the comparison of key antibiotics will lead to inclusion of an isolate when \code{type = "points"}, see Details
#' @param info print progress
#' @param ... parameters passed on to the \code{first_isolate} function
#' @details \strong{WHY THIS IS SO IMPORTANT} \cr
#' To conduct an analysis of antimicrobial resistance, you should only include the first isolate of every patient per episode \href{https://www.ncbi.nlm.nih.gov/pubmed/17304462}{[1]}. If you would not do this, you could easily get an overestimate or underestimate of the resistance of an antibiotic. Imagine that a patient was admitted with an MRSA and that it was found in 5 different blood cultures the following week. The resistance percentage of oxacillin of all \emph{S. aureus} isolates would be overestimated, because you included this MRSA more than once. It would be \href{https://en.wikipedia.org/wiki/Selection_bias}{selection bias}.
#'
#' The functions \code{filter_first_isolate} and \code{filter_first_weighted_isolate} are helper functions to quickly filter on first isolates. The function \code{filter_first_isolate} is essentially equal to:
#' \preformatted{
#'  tbl \%>\%
#'    mutate(only_firsts = first_isolate(tbl, ...)) \%>\%
#'    filter(only_firsts == TRUE) \%>\%
#'    select(-only_firsts)
#' }
#' The function \code{filter_first_weighted_isolate} is essentially equal to:
#' \preformatted{
#'  tbl \%>\%
#'    mutate(keyab = key_antibiotics(.)) \%>\%
#'    mutate(only_weighted_firsts = first_isolate(tbl,
#'                                                col_keyantibiotics = "keyab", ...)) \%>\%
#'    filter(only_weighted_firsts == TRUE) \%>\%
#'    select(-only_weighted_firsts)
#' }
#' @section Key antibiotics:
#'     There are two ways to determine whether isolates can be included as first \emph{weighted} isolates which will give generally the same results: \cr
#'
#'     \strong{1. Using} \code{type = "keyantibiotics"} \strong{and parameter} \code{ignore_I} \cr
#'     Any difference from S to R (or vice versa) will (re)select an isolate as a first weighted isolate. With \code{ignore_I = FALSE}, also differences from I to S|R (or vice versa) will lead to this. This is a reliable method and 30-35 times faster than method 2. Read more about this in the \code{\link{key_antibiotics}} function. \cr
#'
#'     \strong{2. Using} \code{type = "points"} \strong{and parameter} \code{points_threshold} \cr
#'     A difference from I to S|R (or vice versa) means 0.5 points, a difference from S to R (or vice versa) means 1 point. When the sum of points exceeds \code{points_threshold}, which default to \code{2}, an isolate will be (re)selected as a first weighted isolate.
#' @rdname first_isolate
#' @keywords isolate isolates first
#' @seealso \code{\link{key_antibiotics}}
#' @export
#' @importFrom dplyr arrange_at lag between row_number filter mutate arrange pull
#' @importFrom crayon blue bold silver
#' @return Logical vector
#' @source Methodology of this function is based on: \strong{M39 Analysis and Presentation of Cumulative Antimicrobial Susceptibility Test Data, 4th Edition}, 2014, \emph{Clinical and Laboratory Standards Institute (CLSI)}. \url{https://clsi.org/standards/products/microbiology/documents/m39/}.
#' @inheritSection AMR Read more on our website!
#' @examples
#' # septic_patients is a dataset available in the AMR package. It is true, genuine data.
#' ?septic_patients
#'
#' library(dplyr)
#' # Filter on first isolates:
#' septic_patients %>%
#'   mutate(first_isolate = first_isolate(.,
#'                                        col_date = "date",
#'                                        col_patient_id = "patient_id",
#'                                        col_mo = "mo")) %>%
#'   filter(first_isolate == TRUE)
#'
#' # Which can be shortened to:
#' septic_patients %>%
#'   filter_first_isolate()
#' # or for first weighted isolates:
#' septic_patients %>%
#'   filter_first_weighted_isolate()
#'
#' # Now let's see if first isolates matter:
#' A <- septic_patients %>%
#'   group_by(hospital_id) %>%
#'   summarise(count = n_rsi(GEN),            # gentamicin availability
#'             resistance = portion_IR(GEN))  # gentamicin resistance
#'
#' B <- septic_patients %>%
#'   filter_first_weighted_isolate() %>%      # the 1st isolate filter
#'   group_by(hospital_id) %>%
#'   summarise(count = n_rsi(GEN),            # gentamicin availability
#'             resistance = portion_IR(GEN))  # gentamicin resistance
#'
#' # Have a look at A and B.
#' # B is more reliable because every isolate is only counted once.
#' # Gentamicin resitance in hospital D appears to be 3.1% higher than
#' # when you (erroneously) would have used all isolates for analysis.
#'
#'
#' ## OTHER EXAMPLES:
#'
#' \dontrun{
#'
#' # set key antibiotics to a new variable
#' tbl$keyab <- key_antibiotics(tbl)
#'
#' tbl$first_isolate <-
#'   first_isolate(tbl)
#'
#' tbl$first_isolate_weighed <-
#'   first_isolate(tbl,
#'                 col_keyantibiotics = 'keyab')
#'
#' tbl$first_blood_isolate <-
#'   first_isolate(tbl,
#'                 specimen_group = 'Blood')
#'
#' tbl$first_blood_isolate_weighed <-
#'   first_isolate(tbl,
#'                 specimen_group = 'Blood',
#'                 col_keyantibiotics = 'keyab')
#'
#' tbl$first_urine_isolate <-
#'   first_isolate(tbl,
#'                 specimen_group = 'Urine')
#'
#' tbl$first_urine_isolate_weighed <-
#'   first_isolate(tbl,
#'                 specimen_group = 'Urine',
#'                 col_keyantibiotics = 'keyab')
#'
#' tbl$first_resp_isolate <-
#'   first_isolate(tbl,
#'                 specimen_group = 'Respiratory')
#'
#' tbl$first_resp_isolate_weighed <-
#'   first_isolate(tbl,
#'                 specimen_group = 'Respiratory',
#'                 col_keyantibiotics = 'keyab')
#' }
first_isolate <- function(tbl,
                          col_date = NULL,
                          col_patient_id = NULL,
                          col_mo = NULL,
                          col_testcode = NULL,
                          col_specimen = NULL,
                          col_icu = NULL,
                          col_keyantibiotics = NULL,
                          episode_days = 365,
                          testcodes_exclude = NULL,
                          icu_exclude = FALSE,
                          specimen_group = NULL,
                          type = "keyantibiotics",
                          ignore_I = TRUE,
                          points_threshold = 2,
                          info = TRUE,
                          ...) {

  if (!is.data.frame(tbl)) {
    stop("`tbl` must be a data.frame.", call. = FALSE)
  }

  dots <- unlist(list(...))
  if (length(dots) != 0) {
    # backwards compatibility with old parameters
    dots.names <- dots %>% names()
    if ('filter_specimen' %in% dots.names) {
      specimen_group <- dots[which(dots.names == 'filter_specimen')]
    }
  }

  # try to find columns based on type
  # -- mo
  if (is.null(col_mo)) {
    col_mo <- search_type_in_df(tbl = tbl, type = "mo")
  }
  if (is.null(col_mo)) {
    stop("`col_mo` must be set.", call. = FALSE)
  }

  # -- date
  if (is.null(col_date)) {
    col_date <- search_type_in_df(tbl = tbl, type = "date")
  }
  if (is.null(col_date)) {
    stop("`col_date` must be set.", call. = FALSE)
  }
  # convert to Date (pipes/pull for supporting tibbles too)
  tbl[, col_date] <- tbl %>% pull(col_date) %>% as.Date()

  # -- patient id
  if (is.null(col_patient_id)) {
    if (all(c("First name", "Last name", "Sex", "Identification number") %in% colnames(tbl))) {
      # WHONET support
      tbl <- tbl %>% mutate(patient_id = paste(`First name`, `Last name`, Sex))
      col_patient_id <- "patient_id"
      message(blue(paste0("NOTE: Using combined columns ", bold("`First name`, `Last name` and `Sex`"), " as input for `col_patient_id`.")))
    } else {
    col_patient_id <- search_type_in_df(tbl = tbl, type = "patient_id")
    }
  }
  if (is.null(col_patient_id)) {
    stop("`col_patient_id` must be set.", call. = FALSE)
  }

  # -- key antibiotics
  if (is.null(col_keyantibiotics)) {
    col_keyantibiotics <- search_type_in_df(tbl = tbl, type = "keyantibiotics")
  }
  if (isFALSE(col_keyantibiotics)) {
    col_keyantibiotics <- NULL
  }

  # -- specimen
  if (is.null(col_specimen)) {
    col_specimen <- search_type_in_df(tbl = tbl, type = "specimen")
  }
  if (isFALSE(col_specimen)) {
    col_specimen <- NULL
  }

  # check if columns exist
  check_columns_existance <- function(column, tblname = tbl) {
    if (NROW(tblname) <= 1 | NCOL(tblname) <= 1) {
      stop('Please check tbl for existance.')
    }

    if (!is.null(column)) {
      if (!(column %in% colnames(tblname))) {
        stop('Column `', column, '` not found.')
      }
    }
  }

  check_columns_existance(col_date)
  check_columns_existance(col_patient_id)
  check_columns_existance(col_mo)
  check_columns_existance(col_testcode)
  check_columns_existance(col_icu)
  check_columns_existance(col_keyantibiotics)

  # join to microorganisms data set
  tbl <- tbl %>%
    mutate_at(vars(col_mo), as.mo) %>%
    left_join_microorganisms(by = col_mo)
  col_genus <- "genus"
  col_species <- "species"

  if (is.null(col_testcode)) {
    testcodes_exclude <- NULL
  }
  # remove testcodes
  if (!is.null(testcodes_exclude) & info == TRUE) {
    cat('[Criterion] Excluded test codes:\n', toString(testcodes_exclude), '\n')
  }

  if (is.null(col_icu)) {
    icu_exclude <- FALSE
  } else {
    tbl <- tbl %>%
      mutate(col_icu = tbl %>% pull(col_icu) %>% as.logical())
  }

  if (is.null(col_specimen)) {
    specimen_group <- NULL
  }

  # filter on specimen group and keyantibiotics when they are filled in
  if (!is.null(specimen_group)) {
    check_columns_existance(col_specimen, tbl)
    if (info == TRUE) {
      cat('[Criterion] Excluded other than specimen group \'', specimen_group, '\'\n', sep = '')
    }
  }
  if (!is.null(col_keyantibiotics)) {
    tbl <- tbl %>% mutate(key_ab = tbl %>% pull(col_keyantibiotics))
  }

  if (is.null(testcodes_exclude)) {
    testcodes_exclude <- ''
  }

  # create new dataframe with original row index and right sorting
  tbl <- tbl %>%
    mutate(first_isolate_row_index = 1:nrow(tbl),
           date_lab = tbl %>% pull(col_date),
           patient_id = tbl %>% pull(col_patient_id),
           species = tbl %>% pull(col_species),
           genus = tbl %>% pull(col_genus)) %>%
    mutate(species = if_else(is.na(species) | species == "(no MO)", "", species),
           genus = if_else(is.na(genus) | genus == "(no MO)", "", genus))

  if (is.null(specimen_group)) {
    # not filtering on specimen
    if (icu_exclude == FALSE) {
      if (info == TRUE & !is.null(col_icu)) {
        cat('[Criterion] Included isolates from ICU.\n')
      }
      tbl <- tbl %>%
        arrange_at(c(col_patient_id,
                     col_genus,
                     col_species,
                     col_date))
      row.start <- 1
      row.end <- nrow(tbl)
    } else {
      if (info == TRUE) {
        cat('[Criterion] Excluded isolates from ICU.\n')
      }
      tbl <- tbl %>%
        arrange_at(c(col_icu,
                     col_patient_id,
                     col_genus,
                     col_species,
                     col_date))

      suppressWarnings(
        row.start <- which(tbl %>% pull(col_icu) == FALSE) %>% min(na.rm = TRUE)
      )
      suppressWarnings(
        row.end <- which(tbl %>% pull(col_icu) == FALSE) %>% max(na.rm = TRUE)
      )
    }

  } else {
    # filtering on specimen and only analyse these row to save time
    if (icu_exclude == FALSE) {
      if (info == TRUE & !is.null(col_icu)) {
        cat('[Criterion] Included isolates from ICU.\n')
      }
      tbl <- tbl %>%
        arrange_at(c(col_specimen,
                     col_patient_id,
                     col_genus,
                     col_species,
                     col_date))
      suppressWarnings(
        row.start <- which(tbl %>% pull(col_specimen) == specimen_group) %>% min(na.rm = TRUE)
      )
      suppressWarnings(
        row.end <- which(tbl %>% pull(col_specimen) == specimen_group) %>% max(na.rm = TRUE)
      )
    } else {
      if (info == TRUE) {
        cat('[Criterion] Excluded isolates from ICU.\n')
      }
      tbl <- tbl %>%
        arrange_at(c(col_icu,
                     col_specimen,
                     col_patient_id,
                     col_genus,
                     col_species,
                     col_date))
      suppressWarnings(
        row.start <- which(tbl %>% pull(col_specimen) == specimen_group
                           & tbl %>% pull(col_icu) == FALSE) %>% min(na.rm = TRUE)
      )
      suppressWarnings(
        row.end <- which(tbl %>% pull(col_specimen) == specimen_group
                         & tbl %>% pull(col_icu) == FALSE) %>% max(na.rm = TRUE)
      )
    }

  }

  if (abs(row.start) == Inf | abs(row.end) == Inf) {
    if (info == TRUE) {
      message(paste("=> Found", bold("no isolates")))
    }
    # NAs where genus is unavailable
    return(tbl %>%
             mutate(real_first_isolate = if_else(genus == '', NA, FALSE)) %>%
             pull(real_first_isolate)
    )
  }

  # suppress warnings because dplyr wants us to use library(dplyr) when using filter(row_number())
  suppressWarnings(
    scope.size <- tbl %>%
      filter(
        row_number() %>% between(row.start,
                                 row.end),
        genus != "",
        species != "") %>%
      nrow()
  )

  identify_new_year = function(x, episode_days) {
    # I asked on StackOverflow:
    # https://stackoverflow.com/questions/42122245/filter-one-row-every-year
    if (length(x) == 1) {
      return(TRUE)
    }
    indices = integer(0)
    start = x[1]
    ind = 1
    indices[ind] = ind
    for (i in 2:length(x)) {
      if (as.numeric(x[i] - start >= episode_days)) {
        ind = ind + 1
        indices[ind] = i
        start = x[i]
      }
    }
    result <- rep(FALSE, length(x))
    result[indices] <- TRUE
    return(result)
  }

  # Analysis of first isolate ----
  all_first <- tbl %>%
    mutate(other_pat_or_mo = if_else(patient_id == lag(patient_id)
                                     & genus == lag(genus)
                                     & species == lag(species),
                                     FALSE,
                                     TRUE)) %>%
    group_by_at(vars(patient_id,
                     genus,
                     species)) %>%
    mutate(more_than_episode_ago = identify_new_year(x = date_lab,
                                                     episode_days = episode_days)) %>%
    ungroup()

  weighted.notice <- ''
  if (!is.null(col_keyantibiotics)) {
    weighted.notice <- 'weighted '
    if (info == TRUE) {
      if (type == 'keyantibiotics') {
        cat('[Criterion] Inclusion based on key antibiotics, ')
        if (ignore_I == FALSE) {
          cat('not ')
        }
        cat('ignoring I.\n')
      }
      if (type == 'points') {
        cat(paste0('[Criterion] Inclusion based on key antibiotics, using points threshold of '
                   , points_threshold, '.\n'))
      }
    }
    type_param <- type
    # suppress warnings because dplyr want us to use library(dplyr) when using filter(row_number())
    suppressWarnings(
      all_first <- all_first %>%
        mutate(key_ab_lag = lag(key_ab)) %>%
        mutate(key_ab_other = !key_antibiotics_equal(x = key_ab,
                                                     y = key_ab_lag,
                                                     type = type_param,
                                                     ignore_I = ignore_I,
                                                     points_threshold = points_threshold,
                                                     info = info)) %>%
        mutate(
          real_first_isolate =
            if_else(
              between(row_number(), row.start, row.end)
              & genus != ""
              & species != ""
              & (other_pat_or_mo | more_than_episode_ago | key_ab_other),
              TRUE,
              FALSE))
    )
  } else {
    # no key antibiotics
    # suppress warnings because dplyr want us to use library(dplyr) when using filter(row_number())
    suppressWarnings(
      all_first <- all_first %>%
        mutate(
          real_first_isolate =
            if_else(
              between(row_number(), row.start, row.end)
              & genus != ""
              & species != ""
              & (other_pat_or_mo | more_than_episode_ago),
              TRUE,
              FALSE))
    )
  }

  # first one as TRUE
  all_first[row.start, 'real_first_isolate'] <- TRUE
  # no tests that should be included, or ICU
  if (!is.null(col_testcode)) {
    all_first[which(all_first[, col_testcode] %in% tolower(testcodes_exclude)), 'real_first_isolate'] <- FALSE
  }
  if (icu_exclude == TRUE) {
    all_first[which(all_first[, col_icu] == TRUE), 'real_first_isolate'] <- FALSE
  }

  # NAs where genus is unavailable
  all_first <- all_first %>%
    mutate(real_first_isolate = if_else(genus %in% c('', '(no MO)', NA), NA, real_first_isolate))

  all_first <- all_first %>%
    arrange(first_isolate_row_index) %>%
    pull(real_first_isolate)

  if (info == TRUE) {
    decimal.mark <- getOption("OutDec")
    big.mark <- ifelse(decimal.mark != ",", ",", ".")
    n_found <- base::sum(all_first, na.rm = TRUE)
    p_found_total <- percent(n_found / nrow(tbl), force_zero = TRUE)
    p_found_scope <- percent(n_found / scope.size, force_zero = TRUE)
    # mark up number of found
    n_found <- base::format(n_found, big.mark = big.mark, decimal.mark = decimal.mark)
    if (p_found_total != p_found_scope) {
      msg_txt <- paste0("=> Found ",
                        bold(paste0(n_found, " first ", weighted.notice, "isolates")),
                        " (", p_found_scope, " within scope and ", p_found_total, " of total)")
    } else {
      msg_txt <- paste0("=> Found ",
                        bold(paste0(n_found, " first ", weighted.notice, "isolates")),
                        " (", p_found_total, " of total)")
    }
    base::message(msg_txt)
  }

  all_first

}

#' @rdname first_isolate
#' @importFrom dplyr filter
#' @export
filter_first_isolate <- function(tbl,
                                 col_date = NULL,
                                 col_patient_id = NULL,
                                 col_mo = NULL,
                                 ...) {
  filter(tbl, first_isolate(tbl = tbl,
                            col_date = col_date,
                            col_patient_id = col_patient_id,
                            col_mo = col_mo,
                            ...))
}

#' @rdname first_isolate
#' @importFrom dplyr %>% mutate filter
#' @export
filter_first_weighted_isolate <- function(tbl,
                                          col_date = NULL,
                                          col_patient_id = NULL,
                                          col_mo = NULL,
                                          col_keyantibiotics = NULL,
                                          ...) {
  tbl_keyab <- tbl %>%
    mutate(keyab = suppressMessages(key_antibiotics(.,
                                                    col_mo = col_mo,
                                                    ...))) %>%
    mutate(firsts = first_isolate(.,
                                  col_date = col_date,
                                  col_patient_id = col_patient_id,
                                  col_mo = col_mo,
                                  col_keyantibiotics = "keyab",
                                  ...))
  tbl[which(tbl_keyab$firsts == TRUE),]
}
