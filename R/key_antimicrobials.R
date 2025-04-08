# ==================================================================== #
# TITLE:                                                               #
# AMR: An R Package for Working with Antimicrobial Resistance Data     #
#                                                                      #
# SOURCE CODE:                                                         #
# https://github.com/msberends/AMR                                     #
#                                                                      #
# PLEASE CITE THIS SOFTWARE AS:                                        #
# Berends MS, Luz CF, Friedrich AW, et al. (2022).                     #
# AMR: An R Package for Working with Antimicrobial Resistance Data.    #
# Journal of Statistical Software, 104(3), 1-31.                       #
# https://doi.org/10.18637/jss.v104.i03                                #
#                                                                      #
# Developed at the University of Groningen and the University Medical  #
# Center Groningen in The Netherlands, in collaboration with many      #
# colleagues from around the world, see our website.                   #
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
# how to conduct AMR data analysis: https://amr-for-r.org/             #
# ==================================================================== #

#' (Key) Antimicrobials for First Weighted Isolates
#'
#' These functions can be used to determine first weighted isolates by considering the phenotype for isolate selection (see [first_isolate()]). Using a phenotype-based method to determine first isolates is more reliable than methods that disregard phenotypes.
#' @param x A [data.frame] with antimicrobials columns, like `AMX` or `amox`. Can be left blank to determine automatically.
#' @param y,z [character] vectors to compare.
#' @inheritParams first_isolate
#' @param universal Names of **broad-spectrum** antimicrobial drugs, case-insensitive. Set to `NULL` to ignore. See *Details* for the default antimicrobial drugs.
#' @param gram_negative Names of antibiotic drugs for **Gram-positives**, case-insensitive. Set to `NULL` to ignore. See *Details* for the default antibiotic drugs.
#' @param gram_positive Names of antibiotic drugs for **Gram-negatives**, case-insensitive. Set to `NULL` to ignore. See *Details* for the default antibiotic drugs.
#' @param antifungal Names of antifungal drugs for **fungi**, case-insensitive. Set to `NULL` to ignore. See *Details* for the default antifungal drugs.
#' @param only_sir_columns A [logical] to indicate whether only columns must be included that were transformed to class `sir` (see [as.sir()]) on beforehand (default is `FALSE`).
#' @param ... Ignored, only in place to allow future extensions.
#' @details
#' The [key_antimicrobials()] and [all_antimicrobials()] functions are context-aware. This means that the `x` argument can be left blank if used inside a [data.frame] call, see *Examples*.
#'
#' The function [key_antimicrobials()] returns a [character] vector with 12 antimicrobial results for every isolate. The function [all_antimicrobials()] returns a [character] vector with all antimicrobial drug results for every isolate. These vectors can then be compared using [antimicrobials_equal()], to check if two isolates have generally the same antibiogram. Missing and invalid values are replaced with a dot (`"."`) by [key_antimicrobials()] and ignored by [antimicrobials_equal()].
#'
#' Please see the [first_isolate()] function how these important functions enable the 'phenotype-based' method for determination of first isolates.
#'
#' The default antimicrobial drugs used for **all rows** (set in `universal`) are:
#'
#' - Ampicillin
#' - Amoxicillin/clavulanic acid
#' - Cefuroxime
#' - Ciprofloxacin
#' - Piperacillin/tazobactam
#' - Trimethoprim/sulfamethoxazole
#'
#' The default antimicrobial drugs used for **Gram-negative bacteria** (set in `gram_negative`) are:
#'
#' - Cefotaxime
#' - Ceftazidime
#' - Colistin
#' - Gentamicin
#' - Meropenem
#' - Tobramycin
#'
#' The default antimicrobial drugs used for **Gram-positive bacteria** (set in `gram_positive`) are:
#'
#' - Erythromycin
#' - Oxacillin
#' - Rifampin
#' - Teicoplanin
#' - Tetracycline
#' - Vancomycin
#'
#' The default antimicrobial drugs used for **fungi** (set in `antifungal`) are:
#'
#' - Anidulafungin
#' - Caspofungin
#' - Fluconazole
#' - Miconazole
#' - Nystatin
#' - Voriconazole
#' @rdname key_antimicrobials
#' @export
#' @seealso [first_isolate()]
#' @examples
#' # `example_isolates` is a data set available in the AMR package.
#' # See ?example_isolates.
#'
#' # output of the `key_antimicrobials()` function could be like this:
#' strainA <- "SSSRR.S.R..S"
#' strainB <- "SSSIRSSSRSSS"
#'
#' # those strings can be compared with:
#' antimicrobials_equal(strainA, strainB, type = "keyantimicrobials")
#' # TRUE, because I is ignored (as well as missing values)
#'
#' antimicrobials_equal(strainA, strainB, type = "keyantimicrobials", ignore_I = FALSE)
#' # FALSE, because I is not ignored and so the 4th [character] differs
#'
#' \donttest{
#' if (require("dplyr")) {
#'   # set key antimicrobials to a new variable
#'   my_patients <- example_isolates %>%
#'     mutate(keyab = key_antimicrobials(antifungal = NULL)) %>% # no need to define `x`
#'     mutate(
#'       # now calculate first isolates
#'       first_regular = first_isolate(col_keyantimicrobials = FALSE),
#'       # and first WEIGHTED isolates
#'       first_weighted = first_isolate(col_keyantimicrobials = "keyab")
#'     )
#'
#'   # Check the difference in this data set, 'weighted' results in more isolates:
#'   sum(my_patients$first_regular, na.rm = TRUE)
#'   sum(my_patients$first_weighted, na.rm = TRUE)
#' }
#' }
key_antimicrobials <- function(x = NULL,
                               col_mo = NULL,
                               universal = c(
                                 "ampicillin", "amoxicillin/clavulanic acid", "cefuroxime",
                                 "piperacillin/tazobactam", "ciprofloxacin", "trimethoprim/sulfamethoxazole"
                               ),
                               gram_negative = c(
                                 "gentamicin", "tobramycin", "colistin",
                                 "cefotaxime", "ceftazidime", "meropenem"
                               ),
                               gram_positive = c(
                                 "vancomycin", "teicoplanin", "tetracycline",
                                 "erythromycin", "oxacillin", "rifampin"
                               ),
                               antifungal = c(
                                 "anidulafungin", "caspofungin", "fluconazole",
                                 "miconazole", "nystatin", "voriconazole"
                               ),
                               only_sir_columns = FALSE,
                               ...) {
  if (is_null_or_grouped_tbl(x)) {
    # when `x` is left blank, auto determine it (get_current_data() searches underlying data within call)
    # is also fix for using a grouped df as input (a dot as first argument)
    x <- tryCatch(get_current_data(arg_name = "x", call = -2), error = function(e) x)
  }
  meet_criteria(x, allow_class = "data.frame") # also checks dimensions to be >0
  meet_criteria(col_mo, allow_class = "character", has_length = 1, allow_NULL = TRUE, allow_NA = TRUE, is_in = colnames(x))
  meet_criteria(universal, allow_class = "character", allow_NULL = TRUE)
  meet_criteria(gram_negative, allow_class = "character", allow_NULL = TRUE)
  meet_criteria(gram_positive, allow_class = "character", allow_NULL = TRUE)
  meet_criteria(antifungal, allow_class = "character", allow_NULL = TRUE)
  meet_criteria(only_sir_columns, allow_class = "logical", has_length = 1)

  # force regular data.frame, not a tibble or data.table
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  cols <- get_column_abx(x, info = FALSE, only_sir_columns = only_sir_columns, fn = "key_antimicrobials")

  # try to find columns based on type
  # -- mo
  if (is.null(col_mo)) {
    col_mo <- search_type_in_df(x = x, type = "mo", info = FALSE)
  }
  if (is.null(col_mo)) {
    warning_("in `key_antimicrobials()`: no column found for `col_mo`, ignoring antibiotics set in `gram_negative` and `gram_positive`, and antimycotics set in `antifungal`")
    gramstain <- NA_character_
    kingdom <- NA_character_
  } else {
    x.mo <- as.mo(x[, col_mo, drop = TRUE])
    gramstain <- mo_gramstain(x.mo, language = NULL)
    kingdom <- mo_kingdom(x.mo, language = NULL)
  }

  AMR_string <- function(x, values, name, filter, cols = cols) {
    if (is.null(values)) {
      return(rep(NA_character_, length(which(filter))))
    }

    values_old_length <- length(values)
    values <- as.ab(values, flag_multiple_results = FALSE, info = FALSE)
    values <- cols[names(cols) %in% values]
    values_new_length <- length(values)

    if (values_new_length < values_old_length &&
      any(filter, na.rm = TRUE) &&
      message_not_thrown_before("key_antimicrobials", name)) {
      warning_(
        "in `key_antimicrobials()`: ",
        ifelse(values_new_length == 0,
          "No columns available ",
          paste0("Only using ", values_new_length, " out of ", values_old_length, " defined columns ")
        ),
        "as key antimicrobials for ", name, "s. See `?key_antimicrobials`."
      )
    }

    generate_antimicrobials_string(x[which(filter), c(universal, values), drop = FALSE])
  }

  if (is.null(universal)) {
    universal <- character(0)
  } else {
    universal <- as.ab(universal, flag_multiple_results = FALSE, info = FALSE)
    universal <- cols[names(cols) %in% universal]
  }

  key_ab <- rep(NA_character_, nrow(x))

  key_ab[which(gramstain == "Gram-negative")] <- AMR_string(
    x = x,
    values = gram_negative,
    name = "Gram-negative",
    filter = gramstain == "Gram-negative",
    cols = cols
  )

  key_ab[which(gramstain == "Gram-positive")] <- AMR_string(
    x = x,
    values = gram_positive,
    name = "Gram-positive",
    filter = gramstain == "Gram-positive",
    cols = cols
  )

  key_ab[which(kingdom == "Fungi")] <- AMR_string(
    x = x,
    values = antifungal,
    name = "antifungal",
    filter = kingdom == "Fungi",
    cols = cols
  )

  # back-up - only use `universal`
  key_ab[which(is.na(key_ab))] <- AMR_string(
    x = x,
    values = character(0),
    name = "",
    filter = is.na(key_ab),
    cols = cols
  )

  if (length(unique(key_ab)) == 1) {
    warning_("in `key_antimicrobials()`: no distinct key antibiotics determined.")
  }

  key_ab
}

#' @rdname key_antimicrobials
#' @export
all_antimicrobials <- function(x = NULL,
                               only_sir_columns = FALSE,
                               ...) {
  if (is_null_or_grouped_tbl(x)) {
    # when `x` is left blank, auto determine it (get_current_data() searches underlying data within call)
    # is also fix for using a grouped df as input (a dot as first argument)
    x <- tryCatch(get_current_data(arg_name = "x", call = -2), error = function(e) x)
  }
  meet_criteria(x, allow_class = "data.frame") # also checks dimensions to be >0
  meet_criteria(only_sir_columns, allow_class = "logical", has_length = 1)

  # force regular data.frame, not a tibble or data.table
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  cols <- get_column_abx(x,
    only_sir_columns = only_sir_columns, info = FALSE,
    sort = FALSE, fn = "all_antimicrobials"
  )

  generate_antimicrobials_string(x[, cols, drop = FALSE])
}

generate_antimicrobials_string <- function(df) {
  if (NCOL(df) == 0) {
    return(rep("", NROW(df)))
  }
  if (NROW(df) == 0) {
    return(character(0))
  }
  tryCatch(
    {
      do.call(
        paste0,
        lapply(
          as.list(df),
          function(x) {
            x <- toupper(as.character(x))
            x[x == "SDD"] <- "I"
            # ignore "NI" here, no use for determining first isolates
            x[!x %in% c("S", "I", "R")] <- "."
            paste(x)
          }
        )
      )
    },
    error = function(e) rep(strrep(".", NCOL(df)), NROW(df))
  )
}

#' @rdname key_antimicrobials
#' @export
antimicrobials_equal <- function(y,
                                 z,
                                 type = c("points", "keyantimicrobials"),
                                 ignore_I = TRUE,
                                 points_threshold = 2,
                                 ...) {
  meet_criteria(y, allow_class = "character")
  meet_criteria(z, allow_class = "character")
  stop_if(missing(type), "argument \"type\" is missing, with no default")
  meet_criteria(type, allow_class = "character", has_length = 1, is_in = c("points", "keyantimicrobials"))
  meet_criteria(ignore_I, allow_class = "logical", has_length = 1)
  meet_criteria(points_threshold, allow_class = c("numeric", "integer"), has_length = 1, is_positive = TRUE, is_finite = TRUE)
  stop_ifnot(length(y) == length(z), "length of `y` and `z` must be equal")

  key2sir <- function(val) {
    val <- strsplit(val, "", fixed = TRUE)[[1L]]
    val.int <- rep(NA_real_, length(val))
    val.int[val == "S"] <- 1
    val.int[val %in% c("I", "SDD")] <- 2
    val.int[val == "R"] <- 3
    val.int
  }
  # only run on uniques
  uniq <- unique(c(y, z))
  uniq_list <- lapply(uniq, key2sir)
  names(uniq_list) <- uniq

  y <- uniq_list[match(y, names(uniq_list))]
  z <- uniq_list[match(z, names(uniq_list))]

  determine_equality <- function(a, b, type, points_threshold, ignore_I) {
    if (length(a) != length(b)) {
      # incomparable, so not equal
      return(FALSE)
    }
    # ignore NAs on both sides
    NA_ind <- which(is.na(a) | is.na(b))
    a[NA_ind] <- NA_real_
    b[NA_ind] <- NA_real_

    if (type == "points") {
      # count points for every single character:
      # - no change is 0 points
      # - I <-> S|R is 0.5 point
      # - S|R <-> R|S is 1 point
      # use the levels of as.sir (S = 1, I = 2, R = 3)
      # and divide by 2 (S = 0.5, I = 1, R = 1.5)
      (sum(abs(a - b), na.rm = TRUE) / 2) < points_threshold
    } else {
      if (ignore_I == TRUE) {
        ind <- which(a == 2 | b == 2) # since as.double(as.sir("I")) == 2
        a[ind] <- NA_real_
        b[ind] <- NA_real_
      }
      all(a == b, na.rm = TRUE)
    }
  }
  out <- unlist(Map(
    f = determine_equality,
    y,
    z,
    MoreArgs = list(
      type = type,
      points_threshold = points_threshold,
      ignore_I = ignore_I
    ),
    USE.NAMES = FALSE
  ))
  out[is.na(y) | is.na(z)] <- NA
  out
}
