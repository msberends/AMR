# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Analysis                              #
#                                                                      #
# SOURCE                                                               #
# https://gitlab.com/msberends/AMR                                     #
#                                                                      #
# LICENCE                                                              #
# (c) 2018-2020 Berends MS, Luz CF et al.                              #
#                                                                      #
# This R package is free software; you can freely use and distribute   #
# it for both personal and commercial purposes under the terms of the  #
# GNU General Public License version 2.0 (GNU GPL-2), as published by  #
# the Free Software Foundation.                                        #
#                                                                      #
# We created this package for both routine data analysis and academic  #
# research and it was publicly released in the hope that it will be    #
# useful, but it comes WITHOUT ANY WARRANTY OR LIABILITY.              #
# Visit our website for more info: https://msberends.gitlab.io/AMR.    #
# ==================================================================== #

#' Determine first (weighted) isolates
#'
#' Determine first (weighted) isolates of all microorganisms of every patient per episode and (if needed) per specimen type.
#' @inheritSection lifecycle Stable lifecycle
#' @param x a [`data.frame`] containing isolates.
#' @param col_date column name of the result date (or date that is was received on the lab), defaults to the first column of with a date class
#' @param col_patient_id column name of the unique IDs of the patients, defaults to the first column that starts with 'patient' or 'patid' (case insensitive)
#' @param col_mo column name of the IDs of the microorganisms (see [as.mo()]), defaults to the first column of class [`mo`]. Values will be coerced using [as.mo()].
#' @param col_testcode column name of the test codes. Use `col_testcode = NULL` to **not** exclude certain test codes (like test codes for screening). In that case `testcodes_exclude` will be ignored.
#' @param col_specimen column name of the specimen type or group
#' @param col_icu column name of the logicals (`TRUE`/`FALSE`) whether a ward or department is an Intensive Care Unit (ICU)
#' @param col_keyantibiotics column name of the key antibiotics to determine first *weighted* isolates, see [key_antibiotics()]. Defaults to the first column that starts with 'key' followed by 'ab' or 'antibiotics' (case insensitive). Use `col_keyantibiotics = FALSE` to prevent this.
#' @param episode_days episode in days after which a genus/species combination will be determined as 'first isolate' again. The default of 365 days is based on the guideline by CLSI, see Source. 
#' @param testcodes_exclude character vector with test codes that should be excluded (case-insensitive)
#' @param icu_exclude logical whether ICU isolates should be excluded (rows with value `TRUE` in column `col_icu`)
#' @param specimen_group value in column `col_specimen` to filter on
#' @param type type to determine weighed isolates; can be `"keyantibiotics"` or `"points"`, see Details
#' @param ignore_I logical to determine whether antibiotic interpretations with `"I"` will be ignored when `type = "keyantibiotics"`, see Details
#' @param points_threshold points until the comparison of key antibiotics will lead to inclusion of an isolate when `type = "points"`, see Details
#' @param info print progress
#' @param include_unknown logical to determine whether 'unknown' microorganisms should be included too, i.e. microbial code `"UNKNOWN"`, which defaults to `FALSE`. For WHONET users, this means that all records with organism code `"con"` (*contamination*) will be excluded at default. Isolates with a microbial ID of `NA` will always be excluded as first isolate.
#' @param ... parameters passed on to the [first_isolate()] function
#' @details **WHY THIS IS SO IMPORTANT** \cr
#' To conduct an analysis of antimicrobial resistance, you should only include the first isolate of every patient per episode [(ref)](https://www.ncbi.nlm.nih.gov/pubmed/17304462). If you would not do this, you could easily get an overestimate or underestimate of the resistance of an antibiotic. Imagine that a patient was admitted with an MRSA and that it was found in 5 different blood cultures the following week. The resistance percentage of oxacillin of all *S. aureus* isolates would be overestimated, because you included this MRSA more than once. It would be [selection bias](https://en.wikipedia.org/wiki/Selection_bias).
#'
#' All isolates with a microbial ID of `NA` will be excluded as first isolate.
#'
#' The functions [filter_first_isolate()] and [filter_first_weighted_isolate()] are helper functions to quickly filter on first isolates. The function [filter_first_isolate()] is essentially equal to:
#' ```
#'  x %>%
#'    mutate(only_firsts = first_isolate(x, ...)) %>%
#'    filter(only_firsts == TRUE) %>%
#'    select(-only_firsts)
#' ```
#' The function [filter_first_weighted_isolate()] is essentially equal to:
#' ```
#'  x %>%
#'    mutate(keyab = key_antibiotics(.)) %>%
#'    mutate(only_weighted_firsts = first_isolate(x,
#'                                                col_keyantibiotics = "keyab", ...)) %>%
#'    filter(only_weighted_firsts == TRUE) %>%
#'    select(-only_weighted_firsts)
#' ```
#' @section Key antibiotics:
#' There are two ways to determine whether isolates can be included as first *weighted* isolates which will give generally the same results:
#'
#' 1. Using `type = "keyantibiotics"` and parameter `ignore_I`
#' 
#'    Any difference from S to R (or vice versa) will (re)select an isolate as a first weighted isolate. With `ignore_I = FALSE`, also differences from I to S|R (or vice versa) will lead to this. This is a reliable method and 30-35 times faster than method 2. Read more about this in the [key_antibiotics()] function.
#'    
#' 2. Using `type = "points"` and parameter `points_threshold`
#' 
#'    A difference from I to S|R (or vice versa) means 0.5 points, a difference from S to R (or vice versa) means 1 point. When the sum of points exceeds `points_threshold`, which default to `2`, an isolate will be (re)selected as a first weighted isolate.
#' @rdname first_isolate
#' @seealso [key_antibiotics()]
#' @export
#' @importFrom dplyr arrange_at lag between row_number filter mutate arrange pull ungroup
#' @importFrom crayon blue bold silver
# @importFrom clean percentage
#' @return A [`logical`] vector
#' @source Methodology of this function is based on:
#' 
#' **M39 Analysis and Presentation of Cumulative Antimicrobial Susceptibility Test Data, 4th Edition**, 2014, *Clinical and Laboratory Standards Institute (CLSI)*. <https://clsi.org/standards/products/microbiology/documents/m39/>.
#' @inheritSection AMR Read more on our website!
#' @examples
#' # `example_isolates` is a dataset available in the AMR package.
#' # See ?example_isolates.
#'
#' library(dplyr)
#' # Filter on first isolates:
#' example_isolates %>%
#'   mutate(first_isolate = first_isolate(.)) %>%
#'   filter(first_isolate == TRUE)
#' 
#' # Now let's see if first isolates matter:
#' A <- example_isolates %>%
#'   group_by(hospital_id) %>%
#'   summarise(count = n_rsi(GEN),            # gentamicin availability
#'             resistance = resistance(GEN))  # gentamicin resistance
#'
#' B <- example_isolates %>%
#'   filter_first_weighted_isolate() %>%      # the 1st isolate filter
#'   group_by(hospital_id) %>%
#'   summarise(count = n_rsi(GEN),            # gentamicin availability
#'             resistance = resistance(GEN))  # gentamicin resistance
#'
#' # Have a look at A and B.
#' # B is more reliable because every isolate is counted only once.
#' # Gentamicin resitance in hospital D appears to be 3.7% higher than
#' # when you (erroneously) would have used all isolates for analysis.
#'
#'
#' ## OTHER EXAMPLES:
#'
#' \dontrun{
#' 
#' # Short-hand versions:
#' example_isolates %>%
#'   filter_first_isolate()
#'   
#' example_isolates %>%
#'   filter_first_weighted_isolate()
#'
#'
#' # set key antibiotics to a new variable
#' x$keyab <- key_antibiotics(x)
#'
#' x$first_isolate <- first_isolate(x)
#'
#' x$first_isolate_weighed <- first_isolate(x, col_keyantibiotics = 'keyab')
#'
#' x$first_blood_isolate <- first_isolate(x, specimen_group = "Blood")
#' }
first_isolate <- function(x,
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
                          info = interactive(),
                          include_unknown = FALSE,
                          ...) {
  
  if (!is.data.frame(x)) {
    stop("`x` must be a data.frame.", call. = FALSE)
  }
  
  dots <- unlist(list(...))
  if (length(dots) != 0) {
    # backwards compatibility with old parameters
    dots.names <- dots %>% names()
    if ("filter_specimen" %in% dots.names) {
      specimen_group <- dots[which(dots.names == "filter_specimen")]
    }
    if ("tbl" %in% dots.names) {
      x <- dots[which(dots.names == "tbl")]
    }
  }
  
  # try to find columns based on type
  # -- mo
  if (is.null(col_mo)) {
    col_mo <- search_type_in_df(x = x, type = "mo")
  }
  if (is.null(col_mo)) {
    stop("`col_mo` must be set.", call. = FALSE)
  }
  
  # -- date
  if (is.null(col_date)) {
    col_date <- search_type_in_df(x = x, type = "date")
  }
  if (is.null(col_date)) {
    stop("`col_date` must be set.", call. = FALSE)
  }
  # convert to Date (pipes/pull for supporting tibbles too)
  dates <- x %>% pull(col_date) %>% as.Date()
  dates[is.na(dates)] <- as.Date("1970-01-01")
  x[, col_date] <- dates
  
  # -- patient id
  if (is.null(col_patient_id)) {
    if (all(c("First name", "Last name", "Sex") %in% colnames(x))) {
      # WHONET support
      x <- x %>% mutate(patient_id = paste(`First name`, `Last name`, Sex))
      col_patient_id <- "patient_id"
      message(blue(paste0("NOTE: Using combined columns `", bold("First name"), "`, `", bold("Last name"), "` and `", bold("Sex"), "` as input for `col_patient_id`")))
    } else {
      col_patient_id <- search_type_in_df(x = x, type = "patient_id")
    }
  }
  if (is.null(col_patient_id)) {
    stop("`col_patient_id` must be set.", call. = FALSE)
  }
  
  # -- key antibiotics
  if (is.null(col_keyantibiotics)) {
    col_keyantibiotics <- search_type_in_df(x = x, type = "keyantibiotics")
  }
  if (isFALSE(col_keyantibiotics)) {
    col_keyantibiotics <- NULL
  }
  
  # -- specimen
  if (is.null(col_specimen) & !is.null(specimen_group)) {
    col_specimen <- search_type_in_df(x = x, type = "specimen")
  }
  if (isFALSE(col_specimen)) {
    col_specimen <- NULL
  }
  
  # check if columns exist
  check_columns_existance <- function(column, tblname = x) {
    if (NROW(tblname) <= 1 | NCOL(tblname) <= 1) {
      stop("Please check tbl for existance.")
    }
    
    if (!is.null(column)) {
      if (!(column %in% colnames(tblname))) {
        stop("Column `", column, "` not found.")
      }
    }
  }
  
  check_columns_existance(col_date)
  check_columns_existance(col_patient_id)
  check_columns_existance(col_mo)
  check_columns_existance(col_testcode)
  check_columns_existance(col_icu)
  check_columns_existance(col_keyantibiotics)
  
  # create new dataframe with original row index
  x <- x %>%
    mutate(newvar_row_index = seq_len(nrow(x)),
           newvar_mo = x %>% pull(col_mo) %>% as.mo(),
           newvar_genus_species = paste(mo_genus(newvar_mo), mo_species(newvar_mo)),
           newvar_date = x %>% pull(col_date),
           newvar_patient_id = x %>% pull(col_patient_id))
  
  if (is.null(col_testcode)) {
    testcodes_exclude <- NULL
  }
  # remove testcodes
  if (!is.null(testcodes_exclude) & info == TRUE) {
    message(blue(paste0("[Criterion] Excluded test codes: ", toString(testcodes_exclude))))
  }
  
  if (is.null(col_icu)) {
    icu_exclude <- FALSE
  } else {
    x <- x %>%
      mutate(col_icu = x %>% pull(col_icu) %>% as.logical())
  }
  
  if (is.null(col_specimen)) {
    specimen_group <- NULL
  }
  
  # filter on specimen group and keyantibiotics when they are filled in
  if (!is.null(specimen_group)) {
    check_columns_existance(col_specimen, x)
    if (info == TRUE) {
      message(blue(paste0("[Criterion] Excluded other than specimen group '", specimen_group, "'")))
    }
  }
  if (!is.null(col_keyantibiotics)) {
    x <- x %>% mutate(key_ab = x %>% pull(col_keyantibiotics))
  }
  
  if (is.null(testcodes_exclude)) {
    testcodes_exclude <- ""
  }
  
  # arrange data to the right sorting
  if (is.null(specimen_group)) {
    # not filtering on specimen
    if (icu_exclude == FALSE) {
      if (info == TRUE & !is.null(col_icu)) {
        message(blue("[Criterion] Included isolates from ICU"))
      }
      x <- x %>%
        arrange(newvar_patient_id,
                newvar_genus_species,
                newvar_date)
      row.start <- 1
      row.end <- nrow(x)
    } else {
      if (info == TRUE) {
        message(blue("[Criterion] Excluded isolates from ICU"))
      }
      x <- x %>%
        arrange_at(c(col_icu,
                     "newvar_patient_id",
                     "newvar_genus_species",
                     "newvar_date"))
      
      suppressWarnings(
        row.start <- which(x %>% pull(col_icu) == FALSE) %>% min(na.rm = TRUE)
      )
      suppressWarnings(
        row.end <- which(x %>% pull(col_icu) == FALSE) %>% max(na.rm = TRUE)
      )
    }
    
  } else {
    # filtering on specimen and only analyse these row to save time
    if (icu_exclude == FALSE) {
      if (info == TRUE & !is.null(col_icu)) {
        message(blue("[Criterion] Included isolates from ICU.\n"))
      }
      x <- x %>%
        arrange_at(c(col_specimen,
                     "newvar_patient_id",
                     "newvar_genus_species",
                     "newvar_date"))
      suppressWarnings(
        row.start <- which(x %>% pull(col_specimen) == specimen_group) %>% min(na.rm = TRUE)
      )
      suppressWarnings(
        row.end <- which(x %>% pull(col_specimen) == specimen_group) %>% max(na.rm = TRUE)
      )
    } else {
      if (info == TRUE) {
        message(blue("[Criterion] Excluded isolates from ICU"))
      }
      x <- x %>%
        arrange_at(c(col_icu,
                     col_specimen,
                     "newvar_patient_id",
                     "newvar_genus_species",
                     "newvar_date"))
      suppressWarnings(
        row.start <- min(which(x %>% pull(col_specimen) == specimen_group
                               & x %>% pull(col_icu) == FALSE), 
                         na.rm = TRUE)
      )
      suppressWarnings(
        row.end <- max(which(x %>% pull(col_specimen) == specimen_group & 
                               x %>% pull(col_icu) == FALSE),
                       na.rm = TRUE)
      )
    }
    
  }
  
  # no isolates found
  if (abs(row.start) == Inf | abs(row.end) == Inf) {
    if (info == TRUE) {
      message(paste("=> Found", bold("no isolates")))
    }
    return(rep(FALSE, nrow(x)))
  }
  
  # did find some isolates - add new index numbers of rows
  x <- x %>% mutate(newvar_row_index_sorted = seq_len(nrow(.)))
  
  scope.size <- row.end - row.start + 1
  
  identify_new_year <- function(x, episode_days) {
    # I asked on StackOverflow:
    # https://stackoverflow.com/questions/42122245/filter-one-row-every-year
    if (length(x) == 1) {
      return(TRUE)
    }
    indices <- integer(0)
    start <- x[1]
    ind <- 1
    indices[ind] <- ind
    for (i in 2:length(x)) {
      if (isTRUE(as.numeric(x[i] - start) >= episode_days)) {
        ind <- ind + 1
        indices[ind] <- i
        start <- x[i]
      }
    }
    result <- rep(FALSE, length(x))
    result[indices] <- TRUE
    return(result)
  }
  
  # Analysis of first isolate ----
  all_first <- x %>%
    mutate(other_pat_or_mo = if_else(newvar_patient_id == lag(newvar_patient_id)
                                     & newvar_genus_species == lag(newvar_genus_species),
                                     FALSE,
                                     TRUE)) %>%
    group_by(newvar_patient_id,
             newvar_genus_species) %>%
    mutate(more_than_episode_ago = identify_new_year(x = newvar_date,
                                                     episode_days = episode_days)) %>%
    ungroup()
  
  weighted.notice <- ""
  if (!is.null(col_keyantibiotics)) {
    weighted.notice <- "weighted "
    if (info == TRUE) {
      if (type == "keyantibiotics") {
        message(blue(paste0("[Criterion] Inclusion based on key antibiotics, ",
                            ifelse(ignore_I == FALSE, "not ", ""),
                            "ignoring I")))
      }
      if (type == "points") {
        message(blue(paste0("[Criterion] Inclusion based on key antibiotics, using points threshold of "
                            , points_threshold)))
      }
    }
    type_param <- type
    
    all_first <- all_first %>%
      mutate(key_ab_lag = lag(key_ab)) %>%
      mutate(key_ab_other = !key_antibiotics_equal(y = key_ab,
                                                   z = key_ab_lag,
                                                   type = type_param,
                                                   ignore_I = ignore_I,
                                                   points_threshold = points_threshold,
                                                   info = info)) %>%
      mutate(
        real_first_isolate =
          if_else(
            newvar_row_index_sorted %>% between(row.start, row.end)
            & newvar_genus_species != ""
            & (other_pat_or_mo | more_than_episode_ago | key_ab_other),
            TRUE,
            FALSE))
    
  } else {
    # no key antibiotics
    all_first <- all_first %>%
      mutate(
        real_first_isolate =
          if_else(
            newvar_row_index_sorted %>% between(row.start, row.end)
            & newvar_genus_species != ""
            & (other_pat_or_mo | more_than_episode_ago),
            TRUE,
            FALSE))
    
  }
  
  # first one as TRUE
  all_first[row.start, "real_first_isolate"] <- TRUE
  # no tests that should be included, or ICU
  if (!is.null(col_testcode)) {
    all_first[which(all_first[, col_testcode] %in% tolower(testcodes_exclude)), "real_first_isolate"] <- FALSE
  }
  if (icu_exclude == TRUE) {
    all_first[which(all_first[, col_icu] == TRUE), "real_first_isolate"] <- FALSE
  }
  
  decimal.mark <- getOption("OutDec")
  big.mark <- ifelse(decimal.mark != ",", ",", ".")
  
  # handle empty microorganisms
  if (any(all_first$newvar_mo == "UNKNOWN", na.rm = TRUE) & info == TRUE) {
    message(blue(paste0("NOTE: ", ifelse(include_unknown == TRUE, "Included ", "Excluded "), 
                        format(sum(all_first$newvar_mo == "UNKNOWN"),
                               decimal.mark = decimal.mark, big.mark = big.mark), 
                        " isolates with a microbial ID 'UNKNOWN' (column `", bold(col_mo), "`)")))
  }
  all_first[which(all_first$newvar_mo == "UNKNOWN"), "real_first_isolate"] <- include_unknown
  
  # exclude all NAs
  if (any(is.na(all_first$newvar_mo)) & info == TRUE) {
    message(blue(paste0("NOTE: Excluded ", format(sum(is.na(all_first$newvar_mo)),
                                                  decimal.mark = decimal.mark, big.mark = big.mark), 
                        " isolates with a microbial ID 'NA' (column `", bold(col_mo), "`)")))
  }
  all_first[which(is.na(all_first$newvar_mo)), "real_first_isolate"] <- FALSE
  
  # arrange back according to original sorting again
  all_first <- all_first %>%
    arrange(newvar_row_index) %>%
    pull(real_first_isolate)
  
  if (info == TRUE) {
    n_found <- base::sum(all_first, na.rm = TRUE)
    p_found_total <- percentage(n_found / nrow(x))
    p_found_scope <- percentage(n_found / scope.size)
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
filter_first_isolate <- function(x,
                                 col_date = NULL,
                                 col_patient_id = NULL,
                                 col_mo = NULL,
                                 ...) {
  filter(x, first_isolate(x = x,
                          col_date = col_date,
                          col_patient_id = col_patient_id,
                          col_mo = col_mo,
                          ...))
}

#' @rdname first_isolate
#' @importFrom dplyr %>% mutate filter
#' @export
filter_first_weighted_isolate <- function(x,
                                          col_date = NULL,
                                          col_patient_id = NULL,
                                          col_mo = NULL,
                                          col_keyantibiotics = NULL,
                                          ...) {
  tbl_keyab <- x %>%
    mutate(keyab = suppressMessages(key_antibiotics(.,
                                                    col_mo = col_mo,
                                                    ...))) %>%
    mutate(firsts = first_isolate(.,
                                  col_date = col_date,
                                  col_patient_id = col_patient_id,
                                  col_mo = col_mo,
                                  col_keyantibiotics = "keyab",
                                  ...))
  x[which(tbl_keyab$firsts == TRUE), ]
}
