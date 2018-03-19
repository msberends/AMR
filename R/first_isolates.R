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
#' @param col_date column name of the result date (or date that is was received on the lab), supports tidyverse-like quotation
#' @param col_patient_id column name of the unique IDs of the patients, supports tidyverse-like quotation
#' @param col_genus column name of the genus of the microorganisms, supports tidyverse-like quotation
#' @param col_species column name of the species of the microorganisms, supports tidyverse-like quotation
#' @param col_testcode column name of the test codes. Use \code{col_testcode = NA} to \strong{not} exclude certain test codes (like test codes for screening). In that case \code{testcodes_exclude} will be ignored. Supports tidyverse-like quotation.
#' @param col_specimen column name of the specimen type or group, supports tidyverse-like quotation
#' @param col_icu column name of the logicals (\code{TRUE}/\code{FALSE}) whether a ward or department is an Intensive Care Unit (ICU), supports tidyverse-like quotation
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
#' @details \strong{WHY THIS IS SO IMPORTANT} \cr
#'     To conduct an analysis of antimicrobial resistance, you should only include the first isolate of every patient per episode \href{https://www.ncbi.nlm.nih.gov/pubmed/17304462}{[1]}. If you would not do this, you could easily get an overestimate or underestimate of the resistance of an antibiotic. Imagine that a patient was admitted with an MRSA and that it was found in 5 different blood cultures the following week. The resistance percentage of oxacillin of all \emph{S. aureus} isolates would be overestimated, because you included this MRSA more than once. It would be \href{https://en.wikipedia.org/wiki/Selection_bias}{selection bias}.
#'
#'     \strong{DETERMINING WEIGHTED ISOLATES} \cr
#'     \strong{1. Using \code{type = "keyantibiotics"} and parameter \code{ignore_I}} \cr
#'     To determine weighted isolates, the difference between key antibiotics will be checked. Any difference from S to R (or vice versa) will (re)select an isolate as a first weighted isolate. With \code{ignore_I == FALSE}, also differences from I to S|R (or vice versa) will lead to this. This is a reliable and fast method. \cr
#'     \strong{2. Using \code{type = "points"} and parameter \code{points_threshold}} \cr
#'     To determine weighted isolates, difference between antimicrobial interpretations will be measured with points. A difference from I to S|R (or vice versa) means 0.5 points. A difference from S to R (or vice versa) means 1 point. When the sum of points exceeds \code{points_threshold}, an isolate will be (re)selected as a first weighted isolate. This method is being used by the Infection Prevention department (Dr M. Lokate) of the University Medical Center Groningen (UMCG).
#' @keywords isolate isolates first
#' @export
#' @importFrom dplyr arrange_at lag between row_number filter mutate arrange
#' @return A vector to add to table, see Examples.
#' @examples
#' # septic_patients is a dataset available in the AMR package
#' ?septic_patients
#' my_patients <- septic_patients
#' 
#' library(dplyr)
#' my_patients$first_isolate <- my_patients %>%
#'   left_join_bactlist() %>%
#'   first_isolate(col_date = date,
#'                 col_patient_id = patient_id,
#'                 col_genus = genus,
#'                 col_species = species)
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
                          col_genus,
                          col_species,
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
                          info = TRUE) {
  
  # support tidyverse-like quotation
  col_date <- quasiquotate(deparse(substitute(col_date)), col_date)
  col_patient_id <- quasiquotate(deparse(substitute(col_patient_id)), col_patient_id)
  col_genus <- quasiquotate(deparse(substitute(col_genus)), col_genus)
  col_species <- quasiquotate(deparse(substitute(col_species)), col_species)
  col_testcode <- quasiquotate(deparse(substitute(col_testcode)), col_testcode)
  col_specimen <- quasiquotate(deparse(substitute(col_specimen)), col_specimen)
  col_icu <- quasiquotate(deparse(substitute(col_icu)), col_icu)
  col_keyantibiotics <- quasiquotate(deparse(substitute(col_keyantibiotics)), col_keyantibiotics)
  
  # check if columns exist
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
  
  specgroup.notice <- ''
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
  
  scope.size <- tbl %>%
    filter(row_number() %>%
             between(row.start,
                     row.end),
           genus != '') %>%
    nrow()
  
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
        cat('Comparing key antibiotics for first weighted isolates (')
        if (ignore_I == FALSE) {
          cat('NOT ')
        }
        cat('ignoring I)...\n')
      }
      if (type == 'points') {
        cat(paste0('Comparing antibiotics for first weighted isolates (using points threshold of '
                   , points_threshold, ')...\n'))
      }
    }
    type_param <- type
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
  for (i in 1:length(col.list)) {
    if (toupper(col.list[i]) %in% colnames(tbl)) {
      col.list[i] <- toupper(col.list[i])
    } else if (tolower(col.list[i]) %in% colnames(tbl)) {
      col.list[i] <- tolower(col.list[i])
    }
  }
  if (!all(col.list %in% colnames(tbl))) {
    if (info == TRUE) {
      warning('These columns do not exist and will be ignored:\n',
              col.list[!(col.list %in% colnames(tbl))] %>% toString(),
              immediate. = TRUE,
              call. = FALSE)
    }
  }
  amox <- col.list[1]
  cfot <- col.list[2]
  cfta <- col.list[3]
  cftr <- col.list[4]
  cfur <- col.list[5]
  cipr <- col.list[6]
  clar <- col.list[7]
  clin <- col.list[8]
  clox <- col.list[9]
  doxy <- col.list[10]
  gent <- col.list[11]
  line <- col.list[12]
  mero <- col.list[13]
  peni <- col.list[14]
  pita <- col.list[15]
  rifa <- col.list[16]
  teic <- col.list[17]
  trsu <- col.list[18]
  vanc <- col.list[19]
  
  # join bactlist
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
        
      } else if (type == 'keyantibiotics') {
        # check if key antibiotics are exactly the same
        # also possible to ignore I, so only S <-> R and S <-> R are counted
        if (ignore_I == TRUE) {
          valid_chars <- c('S', 's', 'R', 'r')
        } else {
          valid_chars <- c('S', 's', 'I', 'i', 'R', 'r')
        }
        
        # remove invalid values (like "-", NA) on both locations
        x2[which(!x2 %in% valid_chars)] <- '?'
        x2[which(!y2 %in% valid_chars)] <- '?'
        y2[which(!x2 %in% valid_chars)] <- '?'
        y2[which(!y2 %in% valid_chars)] <- '?'
        
        result[i] <- all(x2 == y2)
        
      } else {
        stop('`', type, '` is not a valid value for type, must be `points` or `keyantibiotics`. See ?first_isolate.')
      }
    }
  }
  if (info == TRUE) {
    cat('\n')
  }
  result
}

#' Find bacteria ID based on genus/species
#'
#' Use this function to determine a valid ID based on a genus (and species). This input could be a full name (like \code{"Staphylococcus aureus"}), an abbreviated name (like \code{"S. aureus"}), or just a genus. You could also use a \code{\link{paste}} of a genus and species column to use the full name as input: \code{x = paste(df$genus, df$species)}, where \code{df} is your dataframe.
#' @param x character vector to determine \code{bactid}
#' @export
#' @importFrom dplyr %>% filter slice pull
#' @return Character (vector).
#' @seealso \code{\link{bactlist}} for the dataframe that is being used to determine ID's.
#' @examples 
#' # These examples all return "STAAUR", the ID of S. aureus:
#' guess_bactid("stau")
#' guess_bactid("STAU")
#' guess_bactid("staaur")
#' guess_bactid("S. aureus")
#' guess_bactid("S aureus")
#' guess_bactid("Staphylococcus aureus")
#' guess_bactid("MRSA") # Methicillin-resistant S. aureus
#' guess_bactid("VISA") # Vancomycin Intermediate S. aureus
guess_bactid <- function(x) {
  # remove dots and other non-text in case of "E. coli" except spaces
  x <- gsub("[^a-zA-Z ]+", "", x)
  x.bak <- x
  # replace space by regex sign
  x <- gsub(" ", ".*", x, fixed = TRUE)
  # add start and stop
  x_species <- paste(x, 'species')
  x <- paste0('^', x, '$')
  
  for (i in 1:length(x)) {
    if (tolower(x[i]) == '^e.*coli$') {
      # avoid detection of Entamoeba coli in case of Escherichia coli
      x[i] <- 'Escherichia coli'
    }
    if (tolower(x[i]) == '^st.*au$'
        | tolower(x[i]) == '^stau$'
        | tolower(x[i]) == '^staaur$') {
      # avoid detection of Staphylococcus auricularis in case of S. aureus
      x[i] <- 'Staphylococcus aureus'
    }
    if (tolower(x[i]) == '^p.*aer$') {
      # avoid detection of Pasteurella aerogenes in case of Pseudomonas aeruginosa
      x[i] <- 'Pseudomonas aeruginosa'
    }
    # translate known trivial names to genus+species
    if (toupper(x.bak[i]) == 'MRSA'
        | toupper(x.bak[i]) == 'VISA'
        | toupper(x.bak[i]) == 'VRSA') {
      x[i] <- 'Staphylococcus aureus'
    }
    if (toupper(x.bak[i]) == 'MRSE') {
      x[i] <- 'Staphylococcus epidermidis'
    }
    if (toupper(x.bak[i]) == 'VRE') {
      x[i] <- 'Enterococcus'
    }

    # let's try the ID's first
    found <- AMR::bactlist %>% filter(bactid == x.bak[i])
    
    if (nrow(found) == 0) {
      # now try exact match
      found <- AMR::bactlist %>% filter(fullname == x[i])
    }
    if (nrow(found) == 0) {
      # try any match
      found <- AMR::bactlist %>% filter(fullname %like% x[i])
    }
    if (nrow(found) == 0) {
      # try only genus, with 'species' attached
      found <- AMR::bactlist %>% filter(fullname %like% x_species[i])
    }
    if (nrow(found) == 0) {
      # try splitting of characters and then find ID
      # like esco = E. coli, klpn = K. pneumoniae, stau = S. aureus
      x_length <- nchar(x.bak[i])
      x[i] <- paste0(x.bak[i] %>% substr(1, x_length / 2) %>% trimws(),
                     '.* ',
                     x.bak[i] %>% substr((x_length / 2) + 1, x_length) %>% trimws())
      found <- AMR::bactlist %>% filter(fullname %like% paste0('^', x[i]))
    }
    
    if (nrow(found) != 0) {
      x[i] <- found %>% 
        slice(1) %>% 
        pull(bactid)
    } else {
      x[i] <- ""
    }
  }
  x
}
