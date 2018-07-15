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

#' Determine first (weighted) isolates
#'
#' Determine first (weighted) isolates of all microorganisms of every patient per episode and (if needed) per specimen type.
#' @param tbl a \code{data.frame} containing isolates.
#' @param col_date column name of the result date (or date that is was received on the lab)
#' @param col_patient_id column name of the unique IDs of the patients
#' @param col_bactid column name of the unique IDs of the microorganisms (should occur in the \code{\link{microorganisms}} dataset). Get your bactid's with the function \code{\link{guess_bactid}}, that takes microorganism names as input.
#' @param col_testcode column name of the test codes. Use \code{col_testcode = NA} to \strong{not} exclude certain test codes (like test codes for screening). In that case \code{testcodes_exclude} will be ignored. Supports tidyverse-like quotation.
#' @param col_specimen column name of the specimen type or group
#' @param col_icu column name of the logicals (\code{TRUE}/\code{FALSE}) whether a ward or department is an Intensive Care Unit (ICU)
#' @param col_keyantibiotics column name of the key antibiotics to determine first \emph{weighted} isolates, see \code{\link{key_antibiotics}}. Supports tidyverse-like quotation.
#' @param episode_days episode in days after which a genus/species combination will be determined as 'first isolate' again
#' @param testcodes_exclude character vector with test codes that should be excluded (case-insensitive)
#' @param icu_exclude logical whether ICU isolates should be excluded
#' @param filter_specimen specimen group or type that should be excluded
#' @param output_logical return output as \code{logical} (will else be the values \code{0} or \code{1})
#' @param type type to determine weighed isolates; can be \code{"keyantibiotics"} or \code{"points"}, see Details
#' @param ignore_I logical to determine whether antibiotic interpretations with \code{"I"} will be ignored when \code{type = "keyantibiotics"}, see Details
#' @param points_threshold points until the comparison of key antibiotics will lead to inclusion of an isolate when \code{type = "points"}, see Details
#' @param info print progress
#' @param col_genus (deprecated, use \code{col_bactid} instead) column name of the genus of the microorganisms
#' @param col_species (deprecated, use \code{col_bactid} instead) column name of the species of the microorganisms
#' @details \strong{WHY THIS IS SO IMPORTANT} \cr
#'     To conduct an analysis of antimicrobial resistance, you should only include the first isolate of every patient per episode \href{https://www.ncbi.nlm.nih.gov/pubmed/17304462}{[1]}. If you would not do this, you could easily get an overestimate or underestimate of the resistance of an antibiotic. Imagine that a patient was admitted with an MRSA and that it was found in 5 different blood cultures the following week. The resistance percentage of oxacillin of all \emph{S. aureus} isolates would be overestimated, because you included this MRSA more than once. It would be \href{https://en.wikipedia.org/wiki/Selection_bias}{selection bias}.
#'
#'     \strong{DETERMINING WEIGHTED ISOLATES} \cr
#'     \strong{1. Using} \code{type = "keyantibiotics"} \strong{and parameter} \code{ignore_I} \cr
#'     To determine weighted isolates, the difference between key antibiotics will be checked. Any difference from S to R (or vice versa) will (re)select an isolate as a first weighted isolate. With \code{ignore_I = FALSE}, also differences from I to S|R (or vice versa) will lead to this. This is a reliable method and 30-35 times faster than method 2. \cr
#'     \strong{2. Using} \code{type = "points"} \strong{and parameter} \code{points_threshold} \cr
#'     To determine weighted isolates, difference between antimicrobial interpretations will be measured with points. A difference from I to S|R (or vice versa) means 0.5 points, a difference from S to R (or vice versa) means 1 point. When the sum of points exceeds \code{points_threshold}, an isolate will be (re)selected as a first weighted isolate. This method is being used by the Infection Prevention department (Dr M. Lokate) of the University Medical Center Groningen (UMCG).
#' @keywords isolate isolates first
#' @export
#' @importFrom dplyr arrange_at lag between row_number filter mutate arrange
#' @return A vector to add to table, see Examples.
#' @source Methodology of this function is based on: "M39 Analysis and Presentation of Cumulative Antimicrobial Susceptibility Test Data, 4th Edition", 2014, Clinical and Laboratory Standards Institute. \url{https://clsi.org/standards/products/microbiology/documents/m39/}.
#' @examples
#' # septic_patients is a dataset available in the AMR package
#' ?septic_patients
#' my_patients <- septic_patients
#'
#' library(dplyr)
#' my_patients$first_isolate <- my_patients %>%
#'   first_isolate(col_date = "date",
#'                 col_patient_id = "patient_id",
#'                 col_bactid = "bactid")
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
#'                 filter_specimen = 'Blood')
#'
#' tbl$first_blood_isolate_weighed <-
#'   first_isolate(tbl,
#'                 filter_specimen = 'Blood',
#'                 col_keyantibiotics = 'keyab')
#'
#' tbl$first_urine_isolate <-
#'   first_isolate(tbl,
#'                 filter_specimen = 'Urine')
#'
#' tbl$first_urine_isolate_weighed <-
#'   first_isolate(tbl,
#'                 filter_specimen = 'Urine',
#'                 col_keyantibiotics = 'keyab')
#'
#' tbl$first_resp_isolate <-
#'   first_isolate(tbl,
#'                 filter_specimen = 'Respiratory')
#'
#' tbl$first_resp_isolate_weighed <-
#'   first_isolate(tbl,
#'                 filter_specimen = 'Respiratory',
#'                 col_keyantibiotics = 'keyab')
#' }
first_isolate <- function(tbl,
                          col_date,
                          col_patient_id,
                          col_bactid = NA,
                          col_testcode = NA,
                          col_specimen = NA,
                          col_icu = NA,
                          col_keyantibiotics = NA,
                          episode_days = 365,
                          testcodes_exclude = '',
                          icu_exclude = FALSE,
                          filter_specimen = NA,
                          output_logical = TRUE,
                          type = "keyantibiotics",
                          ignore_I = TRUE,
                          points_threshold = 2,
                          info = TRUE,
                          col_genus = NA,
                          col_species = NA) {

  # bactid OR genus+species must be available
  if (is.na(col_bactid) & (is.na(col_genus) | is.na(col_species))) {
    stop('`col_bactid or both `col_genus` and `col_species` must be available.')
  }

  # check if columns exist
  check_columns_existance <- function(column, tblname = tbl) {
    if (NROW(tblname) <= 1 | NCOL(tblname) <= 1) {
      stop('Please check tbl for existance.')
    }

    if (!is.na(column)) {
      if (!(column %in% colnames(tblname))) {
        stop('Column `', column, '` not found.')
      }
    }
  }

  check_columns_existance(col_date)
  check_columns_existance(col_patient_id)
  check_columns_existance(col_bactid)
  check_columns_existance(col_genus)
  check_columns_existance(col_species)
  check_columns_existance(col_testcode)
  check_columns_existance(col_icu)
  check_columns_existance(col_keyantibiotics)

  if (!is.na(col_bactid)) {
    tbl <- tbl %>% left_join_microorganisms(by = col_bactid)
    col_genus <- "genus"
    col_species <- "species"
  }

  if (is.na(col_testcode)) {
    testcodes_exclude <- NA
  }
  # remove testcodes
  if (!is.na(testcodes_exclude[1]) & testcodes_exclude[1] != '' & info == TRUE) {
    cat('Isolates from these test codes will be ignored:\n', toString(testcodes_exclude), '\n')
  }

  if (is.na(col_icu)) {
    icu_exclude <- FALSE
  } else {
    tbl <- tbl %>%
      mutate(col_icu = tbl %>% pull(col_icu) %>% as.logical())
  }

  if (is.na(col_specimen)) {
    filter_specimen <- ''
  }

  weighted.notice <- ''
  # filter on specimen group and keyantibiotics when they are filled in
  if (!is.na(filter_specimen) & filter_specimen != '') {
    check_columns_existance(col_specimen, tbl)
    if (info == TRUE) {
      cat('Isolates other than of specimen group \'', filter_specimen, '\' will be ignored. ', sep = '')
    }
  } else {
    filter_specimen <- ''
  }
  if (col_keyantibiotics %in% c(NA, '')) {
    col_keyantibiotics <- ''
  } else {
    tbl <- tbl %>% mutate(key_ab = tbl %>% pull(col_keyantibiotics))
  }

  if (is.na(testcodes_exclude[1])) {
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

  if (filter_specimen == '') {

    if (icu_exclude == FALSE) {
      if (info == TRUE) {
        cat('Isolates from ICU will *NOT* be ignored.\n')
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
        cat('Isolates from ICU will be ignored.\n')
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
    # sort on specimen and only analyse these row to save time
    if (icu_exclude == FALSE) {
      if (info == TRUE) {
        cat('Isolates from ICU will *NOT* be ignored.\n')
      }
      tbl <- tbl %>%
        arrange_at(c(col_specimen,
                     col_patient_id,
                     col_genus,
                     col_species,
                     col_date))
      suppressWarnings(
        row.start <- which(tbl %>% pull(col_specimen) == filter_specimen) %>% min(na.rm = TRUE)
      )
      suppressWarnings(
        row.end <- which(tbl %>% pull(col_specimen) == filter_specimen) %>% max(na.rm = TRUE)
      )
    } else {
      if (info == TRUE) {
        cat('Isolates from ICU will be ignored.\n')
      }
      tbl <- tbl %>%
        arrange_at(c(col_icu,
                     col_specimen,
                     col_patient_id,
                     col_genus,
                     col_species,
                     col_date))
      suppressWarnings(
        row.start <- which(tbl %>% pull(col_specimen) == filter_specimen
                           & tbl %>% pull(col_icu) == FALSE) %>% min(na.rm = TRUE)
      )
      suppressWarnings(
        row.end <- which(tbl %>% pull(col_specimen) == filter_specimen
                         & tbl %>% pull(col_icu) == FALSE) %>% max(na.rm = TRUE)
      )
    }

  }

  if (abs(row.start) == Inf | abs(row.end) == Inf) {
    if (info == TRUE) {
      cat('No isolates found.\n')
    }
    # NA's where genus is unavailable
    tbl <- tbl %>%
      mutate(real_first_isolate = if_else(genus == '', NA, FALSE))
    if (output_logical == FALSE) {
      tbl$real_first_isolate <- tbl %>% pull(real_first_isolate) %>% as.integer()
    }
    return(tbl %>% pull(real_first_isolate))
  }

  # suppress warnings because dplyr want us to use library(dplyr) when using filter(row_number())
  suppressWarnings(
    scope.size <- tbl %>%
      filter(
        row_number() %>% between(row.start,
                                 row.end),
        genus != '') %>%
      nrow()
  )

  # Analysis of first isolate ----
  all_first <- tbl %>%
    mutate(other_pat_or_mo = if_else(patient_id == lag(patient_id)
                                     & genus == lag(genus)
                                     & species == lag(species),
                                     FALSE,
                                     TRUE),
           days_diff = 0) %>%
    mutate(days_diff = if_else(other_pat_or_mo == FALSE,
                               (date_lab - lag(date_lab)) + lag(days_diff),
                               0))

  if (col_keyantibiotics != '') {
    if (info == TRUE) {
      if (type == 'keyantibiotics') {
        cat('Key antibiotics for first weighted isolates will be compared (')
        if (ignore_I == FALSE) {
          cat('NOT ')
        }
        cat('ignoring I).')
      }
      if (type == 'points') {
        cat(paste0('Comparing antibiotics for first weighted isolates (using points threshold of '
                   , points_threshold, ')...\n'))
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
              & genus != ''
              & (other_pat_or_mo
                 | days_diff >= episode_days
                 | key_ab_other),
              TRUE,
              FALSE))
    )
    if (info == TRUE) {
      cat('\n')
    }
  } else {
    # suppress warnings because dplyr want us to use library(dplyr) when using filter(row_number())
    suppressWarnings(
      all_first <- all_first %>%
        mutate(
          real_first_isolate =
            if_else(
              between(row_number(), row.start, row.end)
              & genus != ''
              & (other_pat_or_mo
                 | days_diff >= episode_days),
              TRUE,
              FALSE))
    )
  }

  # first one as TRUE
  all_first[row.start, 'real_first_isolate'] <- TRUE
  # no tests that should be included, or ICU
  if (!is.na(col_testcode)) {
    all_first[which(all_first[, col_testcode] %in% tolower(testcodes_exclude)), 'real_first_isolate'] <- FALSE
  }
  if (icu_exclude == TRUE) {
    all_first[which(all_first[, col_icu] == TRUE), 'real_first_isolate'] <- FALSE
  }

  # NA's where genus is unavailable
  all_first <- all_first %>%
    mutate(real_first_isolate = if_else(genus %in% c('', '(no MO)', NA), NA, real_first_isolate))

  all_first <- all_first %>%
    arrange(first_isolate_row_index) %>%
    pull(real_first_isolate)

  if (info == TRUE) {
    cat(paste0('\nFound ',
               all_first %>% sum(na.rm = TRUE),
               ' first ', weighted.notice, 'isolates (',
               (all_first %>% sum(na.rm = TRUE) / scope.size) %>% percent(),
               ' of isolates in scope [where genus was not empty] and ',
               (all_first %>% sum(na.rm = TRUE) / tbl %>% nrow()) %>% percent(),
               ' of total)\n'))
  }

  if (output_logical == FALSE) {
    all_first <- all_first %>% as.integer()
  }

  all_first

}

#' Key antibiotics based on bacteria ID
#'
#' @param tbl table with antibiotics coloms, like \code{amox} and \code{amcl}.
#' @inheritParams first_isolate
#' @param amcl,amox,cfot,cfta,cftr,cfur,cipr,clar,clin,clox,doxy,gent,line,mero,peni,pita,rifa,teic,trsu,vanc column names of antibiotics, case-insensitive
#' @export
#' @importFrom dplyr %>% mutate if_else
#' @return Character of length 1.
#' @seealso \code{\link{mo_property}} \code{\link{antibiotics}}
#' @examples
#' \donttest{
#' #' # set key antibiotics to a new variable
#' tbl$keyab <- key_antibiotics(tbl)
#' }
key_antibiotics <- function(tbl,
                            col_bactid = 'bactid',
                            info = TRUE,
                            amcl = 'amcl',
                            amox = 'amox',
                            cfot = 'cfot',
                            cfta = 'cfta',
                            cftr = 'cftr',
                            cfur = 'cfur',
                            cipr = 'cipr',
                            clar = 'clar',
                            clin = 'clin',
                            clox = 'clox',
                            doxy = 'doxy',
                            gent = 'gent',
                            line = 'line',
                            mero = 'mero',
                            peni = 'peni',
                            pita = 'pita',
                            rifa = 'rifa',
                            teic = 'teic',
                            trsu = 'trsu',
                            vanc = 'vanc') {

  keylist <- character(length = nrow(tbl))

  if (!col_bactid %in% colnames(tbl)) {
    stop('Column ', col_bactid, ' not found.', call. = FALSE)
  }

  # check columns
  col.list <- c(amox, cfot, cfta, cftr, cfur, cipr, clar,
                clin, clox, doxy, gent, line, mero, peni,
                pita, rifa, teic, trsu, vanc)
  col.list <- check_available_columns(tbl = tbl, col.list = col.list, info = info)
  amox <- col.list[amox]
  cfot <- col.list[cfot]
  cfta <- col.list[cfta]
  cftr <- col.list[cftr]
  cfur <- col.list[cfur]
  cipr <- col.list[cipr]
  clar <- col.list[clar]
  clin <- col.list[clin]
  clox <- col.list[clox]
  doxy <- col.list[doxy]
  gent <- col.list[gent]
  line <- col.list[line]
  mero <- col.list[mero]
  peni <- col.list[peni]
  pita <- col.list[pita]
  rifa <- col.list[rifa]
  teic <- col.list[teic]
  trsu <- col.list[trsu]
  vanc <- col.list[vanc]

  # join microorganisms
  tbl <- tbl %>% left_join_microorganisms(col_bactid)

  tbl$key_ab <- NA_character_

  # Staphylococcus
  list_ab <- c(clox, trsu, teic, vanc, doxy, line, clar, rifa)
  list_ab <- list_ab[list_ab %in% colnames(tbl)]
  tbl <- tbl %>% mutate(key_ab =
                          if_else(genus == 'Staphylococcus',
                                  apply(X = tbl[, list_ab],
                                        MARGIN = 1,
                                        FUN = function(x) paste(x, collapse = "")),
                                  key_ab))

  # Rest of Gram +
  list_ab <- c(peni, amox, teic, vanc, clin, line, clar, trsu)
  list_ab <- list_ab[list_ab %in% colnames(tbl)]
  tbl <- tbl %>% mutate(key_ab =
                          if_else(gramstain %like% '^Positive ',
                                  apply(X = tbl[, list_ab],
                                        MARGIN = 1,
                                        FUN = function(x) paste(x, collapse = "")),
                                  key_ab))

  # Gram -
  list_ab <- c(amox, amcl, pita, cfur, cfot, cfta, cftr, mero, cipr, trsu, gent)
  list_ab <- list_ab[list_ab %in% colnames(tbl)]
  tbl <- tbl %>% mutate(key_ab =
                          if_else(gramstain %like% '^Negative ',
                                  apply(X = tbl[, list_ab],
                                        MARGIN = 1,
                                        FUN = function(x) paste(x, collapse = "")),
                                  key_ab))

  # format
  tbl <- tbl %>%
    mutate(key_ab = gsub('(NA|NULL)', '-', key_ab) %>% toupper())

  tbl$key_ab

}

#' @importFrom dplyr progress_estimated %>%
#' @noRd
key_antibiotics_equal <- function(x,
                                  y,
                                  type = c("keyantibiotics", "points"),
                                  ignore_I = TRUE,
                                  points_threshold = 2,
                                  info = FALSE) {
  # x is active row, y is lag
  type <- type[1]

  if (length(x) != length(y)) {
    stop('Length of `x` and `y` must be equal.')
  }

  result <- logical(length(x))

  if (type == "keyantibiotics") {
    if (ignore_I == TRUE) {
      # evaluation using regular expression will treat '?' as any character
      # so I is actually ignored then
      x <- gsub('I', '?', x, ignore.case = TRUE)
      y <- gsub('I', '?', y, ignore.case = TRUE)
    }
    for (i in 1:length(x)) {
      result[i] <- grepl(x = x[i],
                         pattern = y[i],
                         ignore.case = TRUE) |
        grepl(x = y[i],
              pattern = x[i],
              ignore.case = TRUE)
    }
    return(result)
  } else {

    if (info == TRUE) {
      p <- dplyr::progress_estimated(length(x))
    }

    for (i in 1:length(x)) {

      if (info == TRUE) {
        p$tick()$print()
      }

      if (is.na(x[i])) {
        x[i] <- ''
      }
      if (is.na(y[i])) {
        y[i] <- ''
      }

      if (nchar(x[i]) != nchar(y[i])) {

        result[i] <- FALSE

      } else if (x[i] == '' & y[i] == '') {

        result[i] <- TRUE

      } else {

        x2 <- strsplit(x[i], "")[[1]]
        y2 <- strsplit(y[i], "")[[1]]

        if (type == 'points') {
          # count points for every single character:
          # - no change is 0 points
          # - I <-> S|R is 0.5 point
          # - S|R <-> R|S is 1 point
          # use the levels of as.rsi (S = 1, I = 2, R = 3)

          suppressWarnings(x2 <- x2 %>% as.rsi() %>% as.double())
          suppressWarnings(y2 <- y2 %>% as.rsi() %>% as.double())

          points <- (x2 - y2) %>% abs() %>% sum(na.rm = TRUE)
          result[i] <- ((points / 2) >= points_threshold)

        } else {
          stop('`', type, '` is not a valid value for type, must be "points" or "keyantibiotics". See ?first_isolate.')
        }
      }
    }
    if (info == TRUE) {
      cat('\n')
    }
    result
  }
}
