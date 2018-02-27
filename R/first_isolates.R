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
#' @param col_genus column name of the genus of the microorganisms
#' @param col_species column name of the species of the microorganisms
#' @param col_testcode column name of the test codes. Use \code{col_testcode = NA} to \strong{not} exclude certain test codes (like test codes for screening). In that case \code{testcodes_exclude} will be ignored.
#' @param col_specimen column name of the specimen type or group
#' @param col_icu column name of the logicals (\code{TRUE}/\code{FALSE}) whether a ward or department is an Intensive Care Unit (ICU)
#' @param col_keyantibiotics column name of the key antibiotics to determine first \emph{weighted} isolates, see \code{\link{key_antibiotics}}.
#' @param episode_days episode in days after which a genus/species combination will be determined as 'first isolate' again
#' @param testcodes_exclude character vector with test codes that should be excluded (caseINsensitive)
#' @param icu_exclude logical whether ICU isolates should be excluded
#' @param filter_specimen specimen group or type that should be excluded
#' @param output_logical return output as \code{logical} (will else the values \code{0} or \code{1})
#' @param points_threshold points until the comparison of key antibiotics will lead to inclusion of an isolate, see Details
#' @param info print progress
#' @details \strong{Why this is so important} \cr
#'     To conduct an analysis of antimicrobial resistance, you should only include the first isolate of every patient per episode \href{https://www.ncbi.nlm.nih.gov/pubmed/17304462}{[1]}. If you would not do this, you could easily get an overestimate or underestimate of the resistance of an antibiotic. Imagine that a patient was admitted with an MRSA and that it was found in 5 different blood cultures the following week. The resistance percentage of oxacillin of all \emph{S. aureus} isolates would be overestimated, because you included this MRSA more than once. It would be \href{https://en.wikipedia.org/wiki/Selection_bias}{selection bias}.

#'     \strong{\code{points_threshold}} \cr
#'     To compare key antibiotics, the difference between antimicrobial interpretations will be measured. A difference from I to S|R (or vice versa) means 0.5 points. A difference from S to R (or vice versa) means 1 point. When the sum of points exceeds \code{points_threshold}, an isolate will be (re)selected as a first weighted isolate.
#' @keywords isolate isolates first
#' @export
#' @importFrom dplyr arrange_at lag between row_number filter mutate arrange
#' @return A vector to add to table, see Examples.
#' @examples
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
                          col_genus,
                          col_species,
                          col_testcode = NA,
                          col_specimen,
                          col_icu,
                          col_keyantibiotics = NA,
                          episode_days = 365,
                          testcodes_exclude = '',
                          icu_exclude = FALSE,
                          filter_specimen = NA,
                          output_logical = TRUE,
                          points_threshold = 2,
                          info = TRUE) {
  
  # controleren of kolommen wel bestaan
  check_columns_existance <- function(column, tblname = tbl) {
    if (NROW(tblname) <= 1 | NCOL(tblname) <= 1) {
      stop('Please check tbl for existance.')
    }
    
    if (!is.na(column)) {
      if (!(column %in% colnames(tblname))) {
        stop('Column ', column, ' not found.')
      }
    }
  }
  
  check_columns_existance(col_date)
  check_columns_existance(col_patient_id)
  check_columns_existance(col_genus)
  check_columns_existance(col_species)
  check_columns_existance(col_testcode)
  check_columns_existance(col_icu)
  check_columns_existance(col_keyantibiotics)
  
  if (is.na(col_testcode)) {
    testcodes_exclude <- NA
  }
  # testcodes verwijderen die ingevuld zijn
  if (!is.na(testcodes_exclude[1]) & testcodes_exclude[1] != '' & info == TRUE) {
    cat('Isolates from these test codes will be ignored:\n', toString(testcodes_exclude), '\n')
  }
  
  if (is.na(col_icu)) {
    icu_exclude <- FALSE
  } else {
    tbl <- tbl %>%
      mutate(col_icu = tbl %>% pull(col_icu) %>% as.logical())
  }
  
  specgroup.notice <- ''
  weighted.notice <- ''
  # filteren op materiaalgroep en sleutelantibiotica gebruiken wanneer deze ingevuld zijn
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
  
  # nieuwe dataframe maken met de oorspronkelijke rij-index, 0-bepaling en juiste sortering
  #cat('Sorting table...')
  tbl <- tbl %>%
    mutate(first_isolate_row_index = 1:nrow(tbl),
           eersteisolaatbepaling = 0,
           date_lab = tbl %>% pull(col_date),
           patient_id = tbl %>% pull(col_patient_id),
           species = tbl %>% pull(col_species),
           genus = tbl %>% pull(col_genus)) %>%
    mutate(species = if_else(is.na(species), '', species),
           genus = if_else(is.na(genus), '', genus))
  
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
    # sorteren op materiaal en alleen die rijen analyseren om tijd te besparen
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
    # NA's maken waar genus niet beschikbaar is
    tbl <- tbl %>%
      mutate(real_first_isolate = if_else(genus == '', NA, FALSE))
    if (output_logical == FALSE) {
      tbl$real_first_isolate <- tbl %>% pull(real_first_isolate) %>% as.integer()
    }
    return(tbl %>% pull(real_first_isolate))
  }
  
  scope.size <- tbl %>%
    filter(row_number() %>%
             between(row.start,
                     row.end),
           genus != '') %>%
    nrow()
  
  # Analyse van eerste isolaat ----
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
      cat(paste0('Comparing key antibiotics for first weighted isolates (using points threshold of '
                 , points_threshold, ')...\n'))
    }
    all_first <- all_first %>%
      mutate(key_ab_lag = lag(key_ab)) %>%
      mutate(key_ab_other = !key_antibiotics_equal(x = key_ab,
                                                   y = key_ab_lag,
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
    if (info == TRUE) {
      cat('\n')
    }
  } else {
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
  }
  
  # allereerst isolaat als TRUE
  all_first[row.start, 'real_first_isolate'] <- TRUE
  # geen testen die uitgesloten moeten worden, of ICU
  if (!is.na(col_testcode)) {
    all_first[which(all_first[, col_testcode] %in% tolower(testcodes_exclude)), 'real_first_isolate'] <- FALSE
  }
  if (icu_exclude == TRUE) {
    all_first[which(all_first[, col_icu] == TRUE), 'real_first_isolate'] <- FALSE
  }
  
  # NA's maken waar genus niet beschikbaar is
  all_first <- all_first %>%
    mutate(real_first_isolate = if_else(genus == '', NA, real_first_isolate))
  
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
#' @param col_bactcode column of bacteria IDs in \code{tbl}; these should occur in \code{bactlist$bactid}, see \code{\link{bactlist}}
#' @param info print warnings
#' @param amcl,amox,cfot,cfta,cftr,cfur,cipr,clar,clin,clox,doxy,gent,line,mero,peni,pita,rifa,teic,trsu,vanc column names of antibiotics.
#' @export
#' @importFrom dplyr %>% mutate if_else 
#' @return Character of length 1.
#' @seealso \code{\link{mo_property}} \code{\link{ablist}}
#' @examples 
#' \donttest{
#' #' # set key antibiotics to a new variable
#' tbl$keyab <- key_antibiotics(tbl)
#' }
key_antibiotics <- function(tbl,
                            col_bactcode = 'bactid',
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
  
  # check columns
  col.list <- c(amox, cfot, cfta, cftr, cfur, cipr, clar,
                clin, clox, doxy, gent, line, mero, peni,
                pita, rifa, teic, trsu, vanc)
  col.list <- col.list[!is.na(col.list)]
  if (!all(col.list %in% colnames(tbl))) {
    if (info == TRUE) {
      warning('These columns do not exist and will be ignored:\n',
              col.list[!(col.list %in% colnames(tbl))] %>% toString(),
              immediate. = TRUE,
              call. = FALSE)
    }
  }
  
  # bactlist aan vastknopen
  tbl <- tbl %>% left_join_bactlist(col_bactcode)
  
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
                          if_else(gramstain %like% '^Positi[e]?ve',
                                  apply(X = tbl[, list_ab],
                                        MARGIN = 1,
                                        FUN = function(x) paste(x, collapse = "")),
                                  key_ab))
  
  # Gram -
  list_ab <- c(amox, amcl, pita, cfur, cfot, cfta, cftr, mero, cipr, trsu, gent)
  list_ab <- list_ab[list_ab %in% colnames(tbl)]
  tbl <- tbl %>% mutate(key_ab =
                          if_else(gramstain %like% '^Negati[e]?ve',
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
key_antibiotics_equal <- function(x, y, points_threshold = 2, info = FALSE) {
  # x is active row, y is lag

  if (length(x) != length(y)) {
    stop('Length of `x` and `y` must be equal.')
  }
  
  result <- logical(length(x))
  
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
      
      # count points for every single character:
      # - no change is 0 points
      # - I <-> S|R is 0.5 point
      # - S|R <-> R|S is 1 point
      # use the levels of as.rsi (S = 1, I = 2, R = 3)

      x2 <- strsplit(x[i], "")[[1]] %>% as.rsi() %>% as.double()
      y2 <- strsplit(y[i], "")[[1]] %>% as.rsi() %>% as.double()
      
      points <- (x2 - y2) %>% abs() %>% sum(na.rm = TRUE)
      result[i] <- ((points / 2) >= points_threshold)
    }
  }
  if (info == TRUE) {
    cat('\n')
  }
  result
}
