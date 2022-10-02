# ==================================================================== #
# TITLE                                                                #
# AMR: An R Package for Working with Antimicrobial Resistance Data     #
#                                                                      #
# SOURCE                                                               #
# https://github.com/msberends/AMR                                     #
#                                                                      #
# CITE AS                                                              #
# Berends MS, Luz CF, Friedrich AW, Sinha BNM, Albers CJ, Glasner C    #
# (2022). AMR: An R Package for Working with Antimicrobial Resistance  #
# Data. Journal of Statistical Software, 104(3), 1-31.                 #
# doi:10.18637/jss.v104.i03                                            #
#                                                                      #
# Developed at the University of Groningen, the Netherlands, in        #
# collaboration with non-profit organisations Certe Medical            #
# Diagnostics & Advice, and University Medical Center Groningen.       #
#                                                                      #
# This R package is free software; you can freely use and distribute   #
# it for both personal and commercial purposes under the terms of the  #
# GNU General Public License version 2.0 (GNU GPL-2), as published by  #
# the Free Software Foundation.                                        #
# We created this package for both routine data analysis and academic  #
# research and it was publicly released in the hope that it will be    #
# useful, but it comes WITHOUT ANY WARRANTY OR LIABILITY.              #
#                                                                      #
# Visit our website for the full manual and a complete tutorial about  #
# how to conduct AMR data analysis: https://msberends.github.io/AMR/   #
# ==================================================================== #

#' Determine First Isolates
#'
#' Determine first isolates of all microorganisms of every patient per episode and (if needed) per specimen type. These functions support all four methods as summarised by Hindler *et al.* in 2007 (\doi{10.1086/511864}). To determine patient episodes not necessarily based on microorganisms, use [is_new_episode()] that also supports grouping with the `dplyr` package.
#' @param x a [data.frame] containing isolates. Can be left blank for automatic determination, see *Examples*.
#' @param col_date column name of the result date (or date that is was received on the lab), defaults to the first column with a date class
#' @param col_patient_id column name of the unique IDs of the patients, defaults to the first column that starts with 'patient' or 'patid' (case insensitive)
#' @param col_mo column name of the IDs of the microorganisms (see [as.mo()]), defaults to the first column of class [`mo`]. Values will be coerced using [as.mo()].
#' @param col_testcode column name of the test codes. Use `col_testcode = NULL` to **not** exclude certain test codes (such as test codes for screening). In that case `testcodes_exclude` will be ignored.
#' @param col_specimen column name of the specimen type or group
#' @param col_icu column name of the logicals (`TRUE`/`FALSE`) whether a ward or department is an Intensive Care Unit (ICU). This can also be a [logical] vector with the same length as rows in `x`.
#' @param col_keyantimicrobials (only useful when `method = "phenotype-based"`) column name of the key antimicrobials to determine first isolates, see [key_antimicrobials()]. Defaults to the first column that starts with 'key' followed by 'ab' or 'antibiotics' or 'antimicrobials' (case insensitive). Use `col_keyantimicrobials = FALSE` to prevent this. Can also be the output of [key_antimicrobials()].
#' @param episode_days episode in days after which a genus/species combination will be determined as 'first isolate' again. The default of 365 days is based on the guideline by CLSI, see *Source*.
#' @param testcodes_exclude a [character] vector with test codes that should be excluded (case-insensitive)
#' @param icu_exclude a [logical] to indicate whether ICU isolates should be excluded (rows with value `TRUE` in the column set with `col_icu`)
#' @param specimen_group value in the column set with `col_specimen` to filter on
#' @param type type to determine weighed isolates; can be `"keyantimicrobials"` or `"points"`, see *Details*
#' @param method the method to apply, either `"phenotype-based"`, `"episode-based"`, `"patient-based"` or `"isolate-based"` (can be abbreviated), see *Details*. The default is `"phenotype-based"` if antimicrobial test results are present in the data, and `"episode-based"` otherwise.
#' @param ignore_I [logical] to indicate whether antibiotic interpretations with `"I"` will be ignored when `type = "keyantimicrobials"`, see *Details*
#' @param points_threshold minimum number of points to require before differences in the antibiogram will lead to inclusion of an isolate when `type = "points"`, see *Details*
#' @param info a [logical] to indicate info should be printed, defaults to `TRUE` only in interactive mode
#' @param include_unknown a [logical] to indicate whether 'unknown' microorganisms should be included too, i.e. microbial code `"UNKNOWN"`, which defaults to `FALSE`. For WHONET users, this means that all records with organism code `"con"` (*contamination*) will be excluded at default. Isolates with a microbial ID of `NA` will always be excluded as first isolate.
#' @param include_untested_rsi a [logical] to indicate whether also rows without antibiotic results are still eligible for becoming a first isolate. Use `include_untested_rsi = FALSE` to always return `FALSE` for such rows. This checks the data set for columns of class `<rsi>` and consequently requires transforming columns with antibiotic results using [as.rsi()] first.
#' @param ... arguments passed on to [first_isolate()] when using [filter_first_isolate()], otherwise arguments passed on to [key_antimicrobials()] (such as `universal`, `gram_negative`, `gram_positive`)
#' @details
#' To conduct epidemiological analyses on antimicrobial resistance data, only so-called first isolates should be included to prevent overestimation and underestimation of antimicrobial resistance. Different methods can be used to do so, see below.
#'
#' These functions are context-aware. This means that the `x` argument can be left blank if used inside a [data.frame] call, see *Examples*.
#'
#' The [first_isolate()] function is a wrapper around the [is_new_episode()] function, but more efficient for data sets containing microorganism codes or names.
#'
#' All isolates with a microbial ID of `NA` will be excluded as first isolate.
#'
#' ## Different methods
#'
#' According to Hindler *et al.* (2007, \doi{10.1086/511864}), there are different methods (algorithms) to select first isolates with increasing reliability: isolate-based, patient-based, episode-based and phenotype-based. All methods select on a combination of the taxonomic genus and species (not subspecies).
#'
#' All mentioned methods are covered in the [first_isolate()] function:
#'
#'
#' | **Method**                                       | **Function to apply**                                 |
#' |--------------------------------------------------|-------------------------------------------------------|
#' | **Isolate-based**                                | `first_isolate(x, method = "isolate-based")`          |
#' | *(= all isolates)*                               |                                                       |
#' |                                                  |                                                       |
#' |                                                  |                                                       |
#' | **Patient-based**                                | `first_isolate(x, method = "patient-based")`          |
#' | *(= first isolate per patient)*                  |                                                       |
#' |                                                  |                                                       |
#' |                                                  |                                                       |
#' | **Episode-based**                                | `first_isolate(x, method = "episode-based")`, or:     |
#' | *(= first isolate per episode)*                  |                                                       |
#' | - 7-Day interval from initial isolate            | - `first_isolate(x, method = "e", episode_days = 7)`  |
#' | - 30-Day interval from initial isolate           | - `first_isolate(x, method = "e", episode_days = 30)` |
#' |                                                  |                                                       |
#' |                                                  |                                                       |
#' | **Phenotype-based**                              | `first_isolate(x, method = "phenotype-based")`, or:   |
#' | *(= first isolate per phenotype)*                |                                                       |
#' | - Major difference in any antimicrobial result   | - `first_isolate(x, type = "points")`                 |
#' | - Any difference in key antimicrobial results    | - `first_isolate(x, type = "keyantimicrobials")`      |
#'
#' ### Isolate-based
#'
#' This method does not require any selection, as all isolates should be included. It does, however, respect all arguments set in the [first_isolate()] function. For example, the default setting for `include_unknown` (`FALSE`) will omit selection of rows without a microbial ID.
#'
#' ### Patient-based
#'
#' To include every genus-species combination per patient once, set the `episode_days` to `Inf`. Although often inappropriate, this method makes sure that no duplicate isolates are selected from the same patient. In a large longitudinal data set, this could mean that isolates are *excluded* that were found years after the initial isolate.
#'
#' ### Episode-based
#'
#' To include every genus-species combination per patient episode once, set the `episode_days` to a sensible number of days. Depending on the type of analysis, this could be 14, 30, 60 or 365. Short episodes are common for analysing specific hospital or ward data, long episodes are common for analysing regional and national data.
#'
#' This is the most common method to correct for duplicate isolates. Patients are categorised into episodes based on their ID and dates (e.g., the date of specimen receipt or laboratory result). While this is a common method, it does not take into account antimicrobial test results. This means that e.g. a methicillin-resistant *Staphylococcus aureus* (MRSA) isolate cannot be differentiated from a wildtype *Staphylococcus aureus* isolate.
#'
#' ### Phenotype-based
#'
#' This is a more reliable method, since it also *weighs* the antibiogram (antimicrobial test results) yielding so-called 'first weighted isolates'. There are two different methods to weigh the antibiogram:
#'
#' 1. Using `type = "points"` and argument `points_threshold` (default)
#'
#'    This method weighs *all* antimicrobial agents available in the data set. Any difference from I to S or R (or vice versa) counts as `0.5` points, a difference from S to R (or vice versa) counts as `1` point. When the sum of points exceeds `points_threshold`, which defaults to `2`, an isolate will be selected as a first weighted isolate.
#'
#'    All antimicrobials are internally selected using the [all_antimicrobials()] function. The output of this function does not need to be passed to the [first_isolate()] function.
#'
#'
#' 2. Using `type = "keyantimicrobials"` and argument `ignore_I`
#'
#'    This method only weighs specific antimicrobial agents, called *key antimicrobials*. Any difference from S to R (or vice versa) in these key antimicrobials will select an isolate as a first weighted isolate. With `ignore_I = FALSE`, also differences from I to S or R (or vice versa) will lead to this.
#'
#'    Key antimicrobials are internally selected using the [key_antimicrobials()] function, but can also be added manually as a variable to the data and set in the `col_keyantimicrobials` argument. Another option is to pass the output of the [key_antimicrobials()] function directly to the `col_keyantimicrobials` argument.
#'
#'
#' The default method is phenotype-based (using `type = "points"`) and episode-based (using `episode_days = 365`). This makes sure that every genus-species combination is selected per patient once per year, while taking into account all antimicrobial test results. If no antimicrobial test results are available in the data set, only the episode-based method is applied at default.
#' @rdname first_isolate
#' @seealso [key_antimicrobials()]
#' @export
#' @return A [logical] vector
#' @source Methodology of this function is strictly based on:
#'
#' - **M39 Analysis and Presentation of Cumulative Antimicrobial Susceptibility Test Data, 4th Edition**, 2014, *Clinical and Laboratory Standards Institute (CLSI)*. <https://clsi.org/standards/products/microbiology/documents/m39/>.
#'
#' - Hindler JF and Stelling J (2007). **Analysis and Presentation of Cumulative Antibiograms: A New Consensus Guideline from the Clinical and Laboratory Standards Institute.** Clinical Infectious Diseases, 44(6), 867-873. \doi{10.1086/511864}
#' @examples
#' # `example_isolates` is a data set available in the AMR package.
#' # See ?example_isolates.
#'
#' example_isolates[first_isolate(), ]
#' \donttest{
#' # get all first Gram-negatives
#' example_isolates[which(first_isolate(info = FALSE) & mo_is_gram_negative()), ]
#'
#' if (require("dplyr")) {
#'   # filter on first isolates using dplyr:
#'   example_isolates %>%
#'     filter(first_isolate())
#' }
#' if (require("dplyr")) {
#'
#'   # short-hand version:
#'   example_isolates %>%
#'     filter_first_isolate(info = FALSE)
#' }
#' if (require("dplyr")) {
#'
#'   # flag the first isolates per group:
#'   example_isolates %>%
#'     group_by(ward) %>%
#'     mutate(first = first_isolate()) %>%
#'     select(ward, date, patient, mo, first)
#' }
#' }
first_isolate <- function(x = NULL,
                          col_date = NULL,
                          col_patient_id = NULL,
                          col_mo = NULL,
                          col_testcode = NULL,
                          col_specimen = NULL,
                          col_icu = NULL,
                          col_keyantimicrobials = NULL,
                          episode_days = 365,
                          testcodes_exclude = NULL,
                          icu_exclude = FALSE,
                          specimen_group = NULL,
                          type = "points",
                          method = c("phenotype-based", "episode-based", "patient-based", "isolate-based"),
                          ignore_I = TRUE,
                          points_threshold = 2,
                          info = interactive(),
                          include_unknown = FALSE,
                          include_untested_rsi = TRUE,
                          ...) {
  if (is_null_or_grouped_tbl(x)) {
    # when `x` is left blank, auto determine it (get_current_data() also contains dplyr::cur_data_all())
    # is also fix for using a grouped df as input (a dot as first argument)
    x <- tryCatch(get_current_data(arg_name = "x", call = -2), error = function(e) x)
  }
  meet_criteria(x, allow_class = "data.frame") # also checks dimensions to be >0
  meet_criteria(col_date, allow_class = "character", has_length = 1, allow_NULL = TRUE, is_in = colnames(x))
  meet_criteria(col_patient_id, allow_class = "character", has_length = 1, allow_NULL = TRUE, is_in = colnames(x))
  meet_criteria(col_mo, allow_class = "character", has_length = 1, allow_NULL = TRUE, is_in = colnames(x))
  meet_criteria(col_testcode, allow_class = "character", has_length = 1, allow_NULL = TRUE, is_in = colnames(x))
  if (isFALSE(col_specimen)) {
    col_specimen <- NULL
  }
  meet_criteria(col_specimen, allow_class = "character", has_length = 1, allow_NULL = TRUE, is_in = colnames(x))
  if (is.logical(col_icu)) {
    meet_criteria(col_icu, allow_class = "logical", has_length = c(1, nrow(x)), allow_NULL = TRUE)
    if (length(col_icu) == 1) {
      col_icu <- rep(col_icu, nrow(x))
    }
  } else {
    meet_criteria(col_icu, allow_class = c("character", "logical"), has_length = 1, allow_NULL = TRUE, is_in = colnames(x))
    col_icu <- x[, col_icu, drop = TRUE]
  }
  # method
  method <- coerce_method(method)
  meet_criteria(method, allow_class = "character", has_length = 1, is_in = c("phenotype-based", "episode-based", "patient-based", "isolate-based"))
  # key antimicrobials
  if (length(col_keyantimicrobials) > 1) {
    meet_criteria(col_keyantimicrobials, allow_class = "character", has_length = nrow(x))
    x$keyabcol <- col_keyantimicrobials
    col_keyantimicrobials <- "keyabcol"
  } else {
    if (isFALSE(col_keyantimicrobials)) {
      col_keyantimicrobials <- NULL
      # method cannot be phenotype-based anymore
      if (method == "phenotype-based") {
        method <- "episode-based"
      }
    }
    meet_criteria(col_keyantimicrobials, allow_class = "character", has_length = 1, allow_NULL = TRUE, is_in = colnames(x))
  }
  meet_criteria(episode_days, allow_class = c("numeric", "integer"), has_length = 1, is_positive = TRUE, is_finite = FALSE)
  meet_criteria(testcodes_exclude, allow_class = "character", allow_NULL = TRUE)
  meet_criteria(icu_exclude, allow_class = "logical", has_length = 1)
  meet_criteria(specimen_group, allow_class = "character", has_length = 1, allow_NULL = TRUE)
  meet_criteria(type, allow_class = "character", has_length = 1, is_in = c("points", "keyantimicrobials"))
  meet_criteria(ignore_I, allow_class = "logical", has_length = 1)
  meet_criteria(points_threshold, allow_class = c("numeric", "integer"), has_length = 1, is_positive = TRUE, is_finite = TRUE)
  meet_criteria(info, allow_class = "logical", has_length = 1)
  meet_criteria(include_unknown, allow_class = "logical", has_length = 1)
  meet_criteria(include_untested_rsi, allow_class = "logical", has_length = 1)

  # remove data.table, grouping from tibbles, etc.
  x <- as.data.frame(x, stringsAsFactors = FALSE)

  any_col_contains_rsi <- any(vapply(
    FUN.VALUE = logical(1),
    X = x,
    # check only first 10,000 rows
    FUN = function(x) any(as.character(x[1:10000]) %in% c("R", "S", "I"), na.rm = TRUE),
    USE.NAMES = FALSE
  ))
  if (method == "phenotype-based" && !any_col_contains_rsi) {
    method <- "episode-based"
  }
  if (info == TRUE && message_not_thrown_before("first_isolate", "method")) {
    message_(paste0(
      "Determining first isolates ",
      ifelse(method %in% c("episode-based", "phenotype-based"),
        ifelse(is.infinite(episode_days),
          "without a specified episode length",
          paste("using an episode length of", episode_days, "days")
        ),
        ""
      )
    ),
    as_note = FALSE,
    add_fn = font_black
    )
  }

  # try to find columns based on type
  # -- mo
  if (is.null(col_mo)) {
    col_mo <- search_type_in_df(x = x, type = "mo", info = info)
    stop_if(is.null(col_mo), "`col_mo` must be set")
  }

  # methods ----
  if (method == "isolate-based") {
    episode_days <- Inf
    col_keyantimicrobials <- NULL
    x$dummy_dates <- Sys.Date()
    col_date <- "dummy_dates"
    x$dummy_patients <- paste("dummy", seq_len(nrow(x))) # all 'patients' must be unique
    col_patient_id <- "dummy_patients"
  } else if (method == "patient-based") {
    episode_days <- Inf
    col_keyantimicrobials <- NULL
  } else if (method == "episode-based") {
    col_keyantimicrobials <- NULL
  } else if (method == "phenotype-based") {
    if (missing(type) && !is.null(col_keyantimicrobials)) {
      # type = "points" is default, but not set explicitly, while col_keyantimicrobials is
      type <- "keyantimicrobials"
    }
    if (type == "points") {
      x$keyantimicrobials <- all_antimicrobials(x, only_rsi_columns = FALSE)
      col_keyantimicrobials <- "keyantimicrobials"
    } else if (type == "keyantimicrobials" && is.null(col_keyantimicrobials)) {
      col_keyantimicrobials <- search_type_in_df(x = x, type = "keyantimicrobials", info = info)
      if (is.null(col_keyantimicrobials)) {
        # still not found as a column, create it ourselves
        x$keyantimicrobials <- key_antimicrobials(x, only_rsi_columns = FALSE, col_mo = col_mo, ...)
        col_keyantimicrobials <- "keyantimicrobials"
      }
    }
  }

  # -- date
  if (is.null(col_date)) {
    col_date <- search_type_in_df(x = x, type = "date", info = info)
    stop_if(is.null(col_date), "`col_date` must be set")
  }

  # -- patient id
  if (is.null(col_patient_id)) {
    if (all(c("First name", "Last name", "Sex") %in% colnames(x))) {
      # WHONET support
      x$patient_id <- paste(x$`First name`, x$`Last name`, x$Sex)
      col_patient_id <- "patient_id"
      message_("Using combined columns '", font_bold("First name"), "', '", font_bold("Last name"), "' and '", font_bold("Sex"), "' as input for `col_patient_id`")
    } else {
      col_patient_id <- search_type_in_df(x = x, type = "patient_id", info = info)
    }
    stop_if(is.null(col_patient_id), "`col_patient_id` must be set")
  }

  # -- specimen
  if (is.null(col_specimen) && !is.null(specimen_group)) {
    col_specimen <- search_type_in_df(x = x, type = "specimen", info = info)
  }

  # check if columns exist
  check_columns_existance <- function(column, tblname = x) {
    if (!is.null(column)) {
      stop_ifnot(column %in% colnames(tblname),
        "Column '", column, "' not found.",
        call = FALSE
      )
    }
  }

  check_columns_existance(col_date)
  check_columns_existance(col_patient_id)
  check_columns_existance(col_mo)
  check_columns_existance(col_testcode)
  check_columns_existance(col_keyantimicrobials)

  # convert dates to Date
  dates <- as.Date(x[, col_date, drop = TRUE])
  dates[is.na(dates)] <- as.Date("1970-01-01")
  x[, col_date] <- dates

  # create original row index
  x$newvar_row_index <- seq_len(nrow(x))
  x$newvar_mo <- as.mo(x[, col_mo, drop = TRUE])
  x$newvar_genus_species <- paste(mo_genus(x$newvar_mo), mo_species(x$newvar_mo))
  x$newvar_date <- x[, col_date, drop = TRUE]
  x$newvar_patient_id <- x[, col_patient_id, drop = TRUE]

  if (is.null(col_testcode)) {
    testcodes_exclude <- NULL
  }
  # remove testcodes
  if (!is.null(testcodes_exclude) && info == TRUE && message_not_thrown_before("first_isolate", "excludingtestcodes")) {
    message_("Excluding test codes: ", vector_and(testcodes_exclude, quotes = TRUE),
      add_fn = font_black,
      as_note = FALSE
    )
  }

  if (is.null(col_specimen)) {
    specimen_group <- NULL
  }

  # filter on specimen group and keyantibiotics when they are filled in
  if (!is.null(specimen_group)) {
    check_columns_existance(col_specimen, x)
    if (info == TRUE && message_not_thrown_before("first_isolate", "excludingspecimen")) {
      message_("Excluding other than specimen group '", specimen_group, "'",
        add_fn = font_black,
        as_note = FALSE
      )
    }
  }
  if (!is.null(col_keyantimicrobials)) {
    x$newvar_key_ab <- x[, col_keyantimicrobials, drop = TRUE]
  }

  if (is.null(testcodes_exclude)) {
    testcodes_exclude <- ""
  }

  # arrange data to the right sorting
  if (is.null(specimen_group)) {
    x <- x[order(
      x$newvar_patient_id,
      x$newvar_genus_species,
      x$newvar_date
    ), ]
    rownames(x) <- NULL
    row.start <- 1
    row.end <- nrow(x)
  } else {
    # filtering on specimen and only analyse these rows to save time
    x <- x[order(
      pm_pull(x, col_specimen),
      x$newvar_patient_id,
      x$newvar_genus_species,
      x$newvar_date
    ), ]
    rownames(x) <- NULL
    suppressWarnings(
      row.start <- which(x %pm>% pm_pull(col_specimen) == specimen_group) %pm>% min(na.rm = TRUE)
    )
    suppressWarnings(
      row.end <- which(x %pm>% pm_pull(col_specimen) == specimen_group) %pm>% max(na.rm = TRUE)
    )
  }

  # speed up - return immediately if obvious
  if (abs(row.start) == Inf || abs(row.end) == Inf) {
    if (info == TRUE) {
      message_("=> Found ", font_bold("no isolates"),
        add_fn = font_black,
        as_note = FALSE
      )
    }
    return(rep(FALSE, nrow(x)))
  }
  if (row.start == row.end) {
    if (info == TRUE) {
      message_("=> Found ", font_bold("1 first isolate"), ", as the data only contained 1 row",
        add_fn = font_black,
        as_note = FALSE
      )
    }
    return(TRUE)
  }
  if (length(c(row.start:row.end)) == pm_n_distinct(x[c(row.start:row.end), col_mo, drop = TRUE])) {
    if (info == TRUE) {
      message_("=> Found ", font_bold(paste(length(c(row.start:row.end)), "first isolates")),
        ", as all isolates were different microbial species",
        add_fn = font_black,
        as_note = FALSE
      )
    }
    return(rep(TRUE, length(c(row.start:row.end))))
  }

  # did find some isolates - add new index numbers of rows
  x$newvar_row_index_sorted <- seq_len(nrow(x))

  scope.size <- nrow(x[which(x$newvar_row_index_sorted %in% c(row.start + 1:row.end) &
    !is.na(x$newvar_mo)), , drop = FALSE])

  # Analysis of first isolate ----
  if (!is.null(col_keyantimicrobials)) {
    if (info == TRUE && message_not_thrown_before("first_isolate", "type")) {
      if (type == "keyantimicrobials") {
        message_("Basing inclusion on key antimicrobials, ",
          ifelse(ignore_I == FALSE, "not ", ""),
          "ignoring I",
          add_fn = font_black,
          as_note = FALSE
        )
      }
      if (type == "points") {
        message_("Basing inclusion on all antimicrobial results, using a points threshold of ",
          points_threshold,
          add_fn = font_black,
          as_note = FALSE
        )
      }
    }
  }

  x$other_pat_or_mo <- !(x$newvar_patient_id == pm_lag(x$newvar_patient_id) & x$newvar_genus_species == pm_lag(x$newvar_genus_species))

  x$episode_group <- paste(x$newvar_patient_id, x$newvar_genus_species)
  x$more_than_episode_ago <- unlist(lapply(split(
    x$newvar_date,
    x$episode_group
  ),
  exec_episode, # this will skip meet_criteria() in is_new_episode(), saving time
  type = "logical",
  episode_days = episode_days
  ),
  use.names = FALSE
  )

  if (!is.null(col_keyantimicrobials)) {
    # with key antibiotics
    x$other_key_ab <- !antimicrobials_equal(
      y = x$newvar_key_ab,
      z = pm_lag(x$newvar_key_ab),
      type = type,
      ignore_I = ignore_I,
      points_threshold = points_threshold
    )
    x$newvar_first_isolate <- pm_if_else(
      x$newvar_row_index_sorted >= row.start &
        x$newvar_row_index_sorted <= row.end &
        x$newvar_genus_species != "" &
        (x$other_pat_or_mo | x$more_than_episode_ago | x$other_key_ab),
      TRUE,
      FALSE
    )
  } else {
    # no key antibiotics
    x$newvar_first_isolate <- pm_if_else(
      x$newvar_row_index_sorted >= row.start &
        x$newvar_row_index_sorted <= row.end &
        x$newvar_genus_species != "" &
        (x$other_pat_or_mo | x$more_than_episode_ago),
      TRUE,
      FALSE
    )
  }

  # first one as TRUE
  x[row.start, "newvar_first_isolate"] <- TRUE
  # no tests that should be included, or ICU
  if (!is.null(col_testcode)) {
    x[which(x[, col_testcode] %in% tolower(testcodes_exclude)), "newvar_first_isolate"] <- FALSE
  }
  if (!is.null(col_icu)) {
    if (icu_exclude == TRUE) {
      message_("Excluding ", format(sum(col_icu, na.rm = TRUE), big.mark = ","), " isolates from ICU.",
        add_fn = font_black,
        as_note = FALSE
      )
      x[which(col_icu), "newvar_first_isolate"] <- FALSE
    } else {
      message_("Including isolates from ICU.",
        add_fn = font_black,
        as_note = FALSE
      )
    }
  }

  decimal.mark <- getOption("OutDec")
  big.mark <- ifelse(decimal.mark != ",", ",", ".")

  if (info == TRUE) {
    # print group name if used in dplyr::group_by()
    cur_group <- import_fn("cur_group", "dplyr", error_on_fail = FALSE)
    if (!is.null(cur_group)) {
      group_df <- tryCatch(cur_group(), error = function(e) data.frame())
      if (NCOL(group_df) > 0) {
        # transform factors to characters
        group <- vapply(FUN.VALUE = character(1), group_df, function(x) {
          if (is.numeric(x)) {
            format(x)
          } else if (is.logical(x)) {
            as.character(x)
          } else {
            paste0('"', x, '"')
          }
        })
        message_("\nGroup: ", paste0(names(group), " = ", group, collapse = ", "), "\n",
          as_note = FALSE,
          add_fn = font_red
        )
      }
    }
  }

  # handle empty microorganisms
  if (any(x$newvar_mo == "UNKNOWN", na.rm = TRUE) && info == TRUE) {
    message_(
      ifelse(include_unknown == TRUE, "Included ", "Excluded "),
      format(sum(x$newvar_mo == "UNKNOWN", na.rm = TRUE),
        decimal.mark = decimal.mark, big.mark = big.mark
      ),
      " isolates with a microbial ID 'UNKNOWN' (in column '", font_bold(col_mo), "')"
    )
  }
  x[which(x$newvar_mo == "UNKNOWN"), "newvar_first_isolate"] <- include_unknown

  # exclude all NAs
  if (anyNA(x$newvar_mo) && info == TRUE) {
    message_(
      "Excluded ", format(sum(is.na(x$newvar_mo), na.rm = TRUE),
        decimal.mark = decimal.mark, big.mark = big.mark
      ),
      " isolates with a microbial ID 'NA' (in column '", font_bold(col_mo), "')"
    )
  }
  x[which(is.na(x$newvar_mo)), "newvar_first_isolate"] <- FALSE

  # handle isolates without antibiogram
  if (include_untested_rsi == FALSE && any(is.rsi(x))) {
    rsi_all_NA <- which(unname(vapply(
      FUN.VALUE = logical(1),
      as.data.frame(t(x[, is.rsi(x), drop = FALSE])),
      function(rsi_values) all(is.na(rsi_values))
    )))
    x[rsi_all_NA, "newvar_first_isolate"] <- FALSE
  }

  # arrange back according to original sorting again
  x <- x[order(x$newvar_row_index), , drop = FALSE]
  rownames(x) <- NULL

  if (info == TRUE) {
    n_found <- sum(x$newvar_first_isolate, na.rm = TRUE)
    p_found_total <- percentage(n_found / nrow(x[which(!is.na(x$newvar_mo)), , drop = FALSE]), digits = 1)
    p_found_scope <- percentage(n_found / scope.size, digits = 1)
    if (p_found_total %unlike% "[.]") {
      p_found_total <- gsub("%", ".0%", p_found_total, fixed = TRUE)
    }
    if (p_found_scope %unlike% "[.]") {
      p_found_scope <- gsub("%", ".0%", p_found_scope, fixed = TRUE)
    }
    # mark up number of found
    n_found <- format(n_found, big.mark = big.mark, decimal.mark = decimal.mark)
    message_(paste0(
      "=> Found ",
      font_bold(paste0(
        n_found,
        ifelse(method == "isolate-based", "", paste0(" '", method, "'")),
        " first isolates"
      )),
      " (",
      ifelse(p_found_total != p_found_scope,
        paste0(p_found_scope, " within scope and "),
        ""
      ),
      p_found_total, " of total where a microbial ID was available)"
    ),
    add_fn = font_black, as_note = FALSE
    )
  }

  x$newvar_first_isolate
}

#' @rdname first_isolate
#' @export
filter_first_isolate <- function(x = NULL,
                                 col_date = NULL,
                                 col_patient_id = NULL,
                                 col_mo = NULL,
                                 episode_days = 365,
                                 method = c("phenotype-based", "episode-based", "patient-based", "isolate-based"),
                                 ...) {
  if (is_null_or_grouped_tbl(x)) {
    # when `x` is left blank, auto determine it (get_current_data() also contains dplyr::cur_data_all())
    # is also fix for using a grouped df as input (a dot as first argument)
    x <- tryCatch(get_current_data(arg_name = "x", call = -2), error = function(e) x)
  }
  meet_criteria(x, allow_class = "data.frame") # also checks dimensions to be >0
  meet_criteria(col_date, allow_class = "character", has_length = 1, allow_NULL = TRUE, is_in = colnames(x))
  meet_criteria(col_patient_id, allow_class = "character", has_length = 1, allow_NULL = TRUE, is_in = colnames(x))
  meet_criteria(col_mo, allow_class = "character", has_length = 1, allow_NULL = TRUE, is_in = colnames(x))
  meet_criteria(episode_days, allow_class = c("numeric", "integer"), has_length = 1, is_positive = TRUE, is_finite = FALSE)
  method <- coerce_method(method)
  meet_criteria(method, allow_class = "character", has_length = 1, is_in = c("phenotype-based", "episode-based", "patient-based", "isolate-based"))

  subset(x, first_isolate(
    x = x,
    col_date = col_date,
    col_patient_id = col_patient_id,
    col_mo = col_mo,
    episode_days = episode_days,
    method = method,
    ...
  ))
}

coerce_method <- function(method) {
  if (is.null(method)) {
    return(method)
  }
  method <- tolower(as.character(method[1L]))
  method[method %like% "^(p$|pheno)"] <- "phenotype-based"
  method[method %like% "^(e$|episode)"] <- "episode-based"
  method[method %like% "^pat"] <- "patient-based"
  method[method %like% "^(i$|iso)"] <- "isolate-based"
  method
}
