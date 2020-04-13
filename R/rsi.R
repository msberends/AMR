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

#' Class 'rsi'
#'
#' Interpret MIC values and disk diffusion diameters according to EUCAST or CLSI, or clean up existing R/SI values. This transforms the input to a new class [`rsi`], which is an ordered factor with levels `S < I < R`. Invalid antimicrobial interpretations will be translated as `NA` with a warning.
#' @inheritSection lifecycle Stable lifecycle
#' @rdname as.rsi
#' @param x vector of values (for class [`mic`]: an MIC value in mg/L, for class [`disk`]: a disk diffusion radius in millimetres)
#' @param mo any (vector of) text that can be coerced to a valid microorganism code with [as.mo()]
#' @param ab any (vector of) text that can be coerced to a valid antimicrobial code with [as.ab()]
#' @param uti (Urinary Tract Infection) A vector with [logical]s (`TRUE` or `FALSE`) to specify whether a UTI specific interpretation from the guideline should be chosen. For using [as.rsi()] on a [data.frame], this can also be a column containing [logical]s or when left blank, the data set will be search for a 'specimen' and rows containing 'urin' in that column will be regarded isolates from a UTI. See *Examples*.
#' @inheritParams first_isolate
#' @param guideline defaults to the latest included EUCAST guideline, run `unique(rsi_translation$guideline)` for all options
#' @param threshold maximum fraction of invalid antimicrobial interpretations of `x`, please see *Examples*
#' @param ... parameters passed on to methods
#' @details Run `unique(rsi_translation$guideline)` for a list of all supported guidelines. The repository of this package contains [this machine readable version](https://gitlab.com/msberends/AMR/blob/master/data-raw/rsi_translation.txt) of these guidelines.
#' 
#' These guidelines are machine readable, since [](https://gitlab.com/msberends/AMR/blob/master/data-raw/rsi_translation.txt).
#'
#' After using [as.rsi()], you can use [eucast_rules()] to (1) apply inferred susceptibility and resistance based on results of other antimicrobials and (2) apply intrinsic resistance based on taxonomic properties of a microorganism.
#'
#' The function [is.rsi.eligible()] returns `TRUE` when a columns contains at most 5% invalid antimicrobial interpretations (not S and/or I and/or R), and `FALSE` otherwise. The threshold of 5% can be set with the `threshold` parameter.
#' @section Interpretation of R and S/I:
#' In 2019, the European Committee on Antimicrobial Susceptibility Testing (EUCAST) has decided to change the definitions of susceptibility testing categories R and S/I as shown below (<http://www.eucast.org/newsiandr/>).
#'
#' - **R = Resistant**\cr
#'   A microorganism is categorised as *Resistant* when there is a high likelihood of therapeutic failure even when there is increased exposure. Exposure is a function of how the mode of administration, dose, dosing interval, infusion time, as well as distribution and excretion of the antimicrobial agent will influence the infecting organism at the site of infection.
#' - **S = Susceptible**\cr
#'   A microorganism is categorised as *Susceptible, standard dosing regimen*, when there is a high likelihood of therapeutic success using a standard dosing regimen of the agent.
#' - **I = Increased exposure, but still susceptible**\cr
#'   A microorganism is categorised as *Susceptible, Increased exposure* when there is a high likelihood of therapeutic success because exposure to the agent is increased by adjusting the dosing regimen or by its concentration at the site of infection.
#'
#' This AMR package honours this new insight. Use [susceptibility()] (equal to [proportion_SI()]) to determine antimicrobial susceptibility and [count_susceptible()] (equal to [count_SI()]) to count susceptible isolates.
#' @return Ordered factor with new class [`rsi`]
#' @aliases rsi
#' @export
#' @seealso [as.mic()]
#' @inheritSection AMR Read more on our website!
#' @examples
#' # For INTERPRETING disk diffusion and MIC values -----------------------
#'        
#' # a whole data set, even with combined MIC values and disk zones
#' df <- data.frame(microorganism = "E. coli",
#'                  AMP = as.mic(8),
#'                  CIP = as.mic(0.256),
#'                  GEN = as.disk(18),
#'                  TOB = as.disk(16),
#'                  NIT = as.mic(32))
#' as.rsi(df)
#' 
#' \donttest{
#' 
#' # the dplyr way
#' library(dplyr)
#' df %>%
#'   mutate_at(vars(AMP:TOB), as.rsi, mo = "E. coli")
#'   
#' df %>%
#'   mutate_at(vars(AMP:TOB), as.rsi, mo = .$microorganism)
#'   
#' # to include information about urinary tract infections (UTI)
#' data.frame(mo = "E. coli",
#'            NIT = c("<= 2", 32),
#'            from_the_bladder = c(TRUE, FALSE)) %>%
#'   as.rsi(uti = "from_the_bladder")
#'   
#' data.frame(mo = "E. coli",
#'            NIT = c("<= 2", 32),
#'            specimen = c("urine", "blood")) %>%
#'   as.rsi() # automatically determines urine isolates
#' 
#' df %>%
#'   mutate_at(vars(AMP:NIT), as.rsi, mo = "E. coli", uti = TRUE)  
#' }
#' 
#' # for single values
#' as.rsi(x = as.mic(2),
#'        mo = as.mo("S. pneumoniae"),
#'        ab = "AMP",
#'        guideline = "EUCAST")
#'
#' as.rsi(x = as.disk(18),
#'        mo = "Strep pneu",  # `mo` will be coerced with as.mo()
#'        ab = "ampicillin",  # and `ab` with as.ab()
#'        guideline = "EUCAST")
#'
#'
#' # For CLEANING existing R/SI values ------------------------------------
#' 
#' as.rsi(c("S", "I", "R", "A", "B", "C"))
#' as.rsi("<= 0.002; S") # will return "S"
#' 
#' rsi_data <- as.rsi(c(rep("S", 474), rep("I", 36), rep("R", 370)))
#' is.rsi(rsi_data)
#' plot(rsi_data)    # for percentages
#' barplot(rsi_data) # for frequencies
#' freq(rsi_data)    # frequency table with informative header
#'
#' library(dplyr)
#' example_isolates %>%
#'   mutate_at(vars(PEN:RIF), as.rsi)
#'
#' # fastest way to transform all columns with already valid AMR results to class `rsi`:
#' example_isolates %>%
#'   mutate_if(is.rsi.eligible, as.rsi)
#'   
#' # note: from dplyr 1.0.0 on, this will be: 
#' # example_isolates %>%
#' #   mutate(across(is.rsi.eligible, as.rsi))
#'
#' # default threshold of `is.rsi.eligible` is 5%.
#' is.rsi.eligible(WHONET$`First name`) # fails, >80% is invalid
#' is.rsi.eligible(WHONET$`First name`, threshold = 0.99) # succeeds
as.rsi <- function(x, ...) {
  UseMethod("as.rsi")
}

#' @export
as.rsi.default <- function(x, ...) {
  if (is.rsi(x)) {
    x
  } else if (identical(levels(x), c("S", "I", "R"))) {
    structure(x, class = c("rsi", "ordered", "factor"))
  } else if (inherits(x, "integer") & all(x %in% c(1:3, NA))) {
    x[x == 1] <- "S"
    x[x == 2] <- "I"
    x[x == 3] <- "R"
    structure(.Data = factor(x, levels = c("S", "I", "R"), ordered = TRUE),
              class =  c("rsi", "ordered", "factor"))
  } else {

    ab <- deparse(substitute(x))
    if (!any(x %like% "(R|S|I)", na.rm = TRUE)) {
      if (!is.na(suppressWarnings(as.ab(ab)))) {
        # check if they are actually MICs or disks now that the antibiotic name is valid
        if (all_valid_mics(x)) {
          as.rsi(as.mic(x), ab = ab, ...)
        } else if (all_valid_disks(x)) {
          as.rsi(as.disk(x), ab = ab, ...)
        }
      }
    }
    
    x <- x %>% unlist()
    x.bak <- x
    
    na_before <- x[is.na(x) | x == ""] %>% length()
    # remove all spaces
    x <- gsub(" +", "", x)
    # remove all MIC-like values: numbers, operators and periods
    x <- gsub("[0-9.,;:<=>]+", "", x)
    # remove everything between brackets, and 'high' and 'low'
    x <- gsub("([(].*[)])", "", x)
    x <- gsub("(high|low)", "", x, ignore.case = TRUE)
    # disallow more than 3 characters
    x[nchar(x) > 3] <- NA
    # set to capitals
    x <- toupper(x)
    # remove all invalid characters
    x <- gsub("[^RSI]+", "", x)
    # in cases of "S;S" keep S, but in case of "S;I" make it NA
    x <- gsub("^S+$", "S", x)
    x <- gsub("^I+$", "I", x)
    x <- gsub("^R+$", "R", x)
    x[!x %in% c("S", "I", "R")] <- NA
    na_after <- x[is.na(x) | x == ""] %>% length()
    
    if (!isFALSE(list(...)$warn)) { # so as.rsi(..., warn = FALSE) will never throw a warning
      if (na_before != na_after) {
        list_missing <- x.bak[is.na(x) & !is.na(x.bak) & x.bak != ""] %>%
          unique() %>%
          sort()
        list_missing <- paste0('"', list_missing, '"', collapse = ", ")
        warning(na_after - na_before, " results truncated (",
                round(((na_after - na_before) / length(x)) * 100),
                "%) that were invalid antimicrobial interpretations: ",
                list_missing, call. = FALSE)
      }
    }
    
    structure(.Data = factor(x, levels = c("S", "I", "R"), ordered = TRUE),
              class =  c("rsi", "ordered", "factor"))
  }
}

#' @rdname as.rsi
#' @export
as.rsi.mic <- function(x, mo, ab = deparse(substitute(x)), guideline = "EUCAST", uti = FALSE, ...) {
  if (missing(mo)) {
    stop('No information was supplied about the microorganisms (missing parameter "mo"). See ?as.rsi.\n\n',
         "To transform certain columns with e.g. mutate_at(), use\n",
         "`data %>% mutate_at(vars(...), as.rsi, mo = .$x)`, where x is your column with microorganisms.\n\n",
         "To tranform all MIC variables in a data set, use `as.rsi(data)` or `data %>% as.rsi()`.", call. = FALSE)
  }
  
  ab_coerced <- suppressWarnings(as.ab(ab))
  mo_coerced <- suppressWarnings(as.mo(mo))
  guideline_coerced <- get_guideline(guideline)
  if (is.na(ab_coerced)) {
    message(red(paste0("Unknown drug: `", bold(ab), "`. Rename this column to a drug name or code, and check the output with as.ab().")))
    return(as.rsi(rep(NA, length(x))))
  }
  if (length(mo_coerced) == 1) {
    mo_coerced <- rep(mo_coerced, length(x))
  }
  if (length(uti) == 1) {
    uti <- rep(uti, length(x))
  }
  
  message(blue(paste0("=> Interpreting MIC values of `", bold(ab), "` (",
                      ifelse(ab_coerced != ab, paste0(ab_coerced, ", "), ""),
                      ab_name(ab_coerced, tolower = TRUE), ") using guideline ", bold(guideline_coerced), " ... ")),
          appendLF = FALSE)
  result <- exec_as.rsi(method = "mic",
                        x = x,
                        mo = mo_coerced,
                        ab = ab_coerced,
                        guideline = guideline_coerced,
                        uti = uti) # exec_as.rsi will return message(blue(" OK."))
  result
}

#' @rdname as.rsi
#' @export
as.rsi.disk <- function(x, mo, ab = deparse(substitute(x)), guideline = "EUCAST", uti = FALSE, ...) {
  if (missing(mo)) {
    stop('No information was supplied about the microorganisms (missing parameter "mo"). See ?as.rsi.\n\n',
         "To transform certain columns with e.g. mutate_at(), use\n",
         "`data %>% mutate_at(vars(...), as.rsi, mo = .$x)`, where x is your column with microorganisms.\n\n",
         "To tranform all disk diffusion zones in a data set, use `as.rsi(data)` or `data %>% as.rsi()`.", call. = FALSE)
  }
  
  ab_coerced <- suppressWarnings(as.ab(ab))
  mo_coerced <- suppressWarnings(as.mo(mo))
  guideline_coerced <- get_guideline(guideline)
  if (is.na(ab_coerced)) {
    message(red(paste0("Unknown drug: `", bold(ab), "`. Rename this column to a drug name or code, and check the output with as.ab().")))
    return(as.rsi(rep(NA, length(x))))
  }
  if (length(mo_coerced) == 1) {
    mo_coerced <- rep(mo_coerced, length(x))
  }
  if (length(uti) == 1) {
    uti <- rep(uti, length(x))
  }
  
  message(blue(paste0("=> Interpreting disk zones of `", bold(ab), "` (",
                      ifelse(ab_coerced != ab, paste0(ab_coerced, ", "), ""),
                      ab_name(ab_coerced, tolower = TRUE), ") using guideline ", bold(guideline_coerced), " ... ")),
          appendLF = FALSE)
  result <- exec_as.rsi(method = "disk",
                        x = x,
                        mo = mo_coerced,
                        ab = ab_coerced,
                        guideline = guideline_coerced,
                        uti = uti) # exec_as.rsi will return message(blue(" OK."))
  result
}

#' @rdname as.rsi
#' @importFrom crayon red blue bold
#' @export
as.rsi.data.frame <- function(x, col_mo = NULL, guideline = "EUCAST", uti = NULL, ...) {
  # try to find columns based on type
  # -- mo
  if (is.null(col_mo)) {
    col_mo <- search_type_in_df(x = x, type = "mo")
  }
  if (is.null(col_mo)) {
    stop("`col_mo` must be set.", call. = FALSE)
  }
  # -- UTIs
  col_uti <- uti
  if (is.null(col_uti)) {
    col_uti <- search_type_in_df(x = x, type = "uti")
  }
  if (!is.null(col_uti)) {
    if (is.logical(col_uti)) {
      # already a logical vector as input
      if (length(col_uti) == 1) {
        uti <- rep(col_uti, NROW(x))
      } else {
        uti <- col_uti
      }
    } else {
      # column found, transform to logical
      uti <- as.logical(x[, col_uti, drop = TRUE])
    }
  } else {
    # look for specimen column and make logicals of the urines
    col_specimen <- suppressMessages(search_type_in_df(x = x, type = "specimen"))
    if (!is.null(col_specimen)) {
      uti <- x[, col_specimen, drop = TRUE] %like% "urin"
      values <- sort(unique(x[uti, col_specimen, drop = TRUE]))
      if (length(values) > 1) {
        plural <- c("s", "", "")
      } else {
        plural <- c("", "s", "a ")
      }
      message(blue(paste0("NOTE: Assuming value", plural[1], " ", 
              paste(paste0('"', values, '"'), collapse = ", "),
              " in column `", bold(col_specimen),
              "` reflect", plural[2], " ", plural[3], "urinary tract infection", plural[1], ".\n  Use `as.rsi(uti = FALSE)` to prevent this.")))
    } else {
      # no data about UTI's found
      uti <- FALSE
    }
  }

  i <- 0
  ab_cols <- colnames(x)[sapply(x, function(y) {
    i <<- i + 1
    check <- is.mic(y) | is.disk(y)
    ab <- colnames(x)[i]
    ab_coerced <- suppressWarnings(as.ab(ab))
    if (is.na(ab_coerced)) {
      # not even a valid AB code
      return(FALSE)
    } else if (!check & all_valid_mics(y)) {
      message(blue(paste0("NOTE: Assuming column `", ab, "` (",
                          ifelse(ab_coerced != ab, paste0(ab_coerced, ", "), ""),
                          ab_name(ab_coerced, tolower = TRUE), ") contains MIC values.")))
      return(TRUE)
    } else if (!check & all_valid_disks(y)) {
      message(blue(paste0("NOTE: Assuming column `", ab, "` (",
                          ifelse(ab_coerced != ab, paste0(ab_coerced, ", "), ""),
                          ab_name(ab_coerced, tolower = TRUE), ") contains disk zones.")))
      return(TRUE)
    } else {
      return(check)
    }
  })]
  
  if (length(ab_cols) == 0) {
    stop("No columns with MIC values or disk zones found in this data set. Use as.mic() or as.disk() to transform antimicrobial columns.", call. = FALSE)
  }
  
  # set type per column
  types <- character(length(ab_cols))
  types[sapply(x[, ab_cols], is.mic)] <- "mic"
  types[types == "" & sapply(x[, ab_cols], all_valid_mics)] <- "mic"
  types[sapply(x[, ab_cols], is.disk)] <- "disk"
  types[types == "" & sapply(x[, ab_cols], all_valid_disks)] <- "disk"
  
  for (i in seq_len(length(ab_cols))) {
    if (types[i] == "mic") {
      x[, ab_cols[i]] <- as.rsi.mic(x = x %>% pull(ab_cols[i]),
                                    mo = x %>% pull(col_mo),
                                    ab = ab_cols[i],
                                    guideline = guideline,
                                    uti = uti)
    } else if (types[i] == "disk") {
      x[, ab_cols[i]] <- as.rsi.disk(x = x %>% pull(ab_cols[i]),
                                     mo = x %>% pull(col_mo),
                                     ab = ab_cols[i],
                                     guideline = guideline,
                                     uti = uti)
    }
  }
  
  x
}

#' @importFrom dplyr %>% filter pull
get_guideline <- function(guideline) {
  guideline_param <- toupper(guideline)
  if (guideline_param %in% c("CLSI", "EUCAST")) {
    guideline_param <- rsi_translation %>%
      filter(guideline %like% guideline_param) %>%
      pull(guideline) %>%
      sort() %>%
      rev() %>%
      .[1]
  }
  if (!guideline_param %like% " ") {
    # like 'EUCAST2020', should be 'EUCAST 2020'
    guideline_param <- gsub("([a-z]+)([0-9]+)", "\\1 \\2", guideline_param, ignore.case = TRUE)
  }
  
  if (!guideline_param %in% rsi_translation$guideline) {
    stop(paste0("invalid guideline: '", guideline,
                "'.\nValid guidelines are: ", paste0("'", unique(rsi_translation$guideline), "'", collapse = ", "), "."),
         call. = FALSE)
  }
  
  guideline_param
  
}

#' @importFrom dplyr %>% case_when desc arrange filter n_distinct
#' @importFrom crayon green red bold
exec_as.rsi <- function(method, x, mo, ab, guideline, uti) {
  if (method == "mic") {
    x <- as.mic(x) # when as.rsi.mic is called directly
  } else if (method == "disk") {
    x <- as.disk(x) # when as.rsi.disk is called directly
  }
  
  warned <- FALSE
  method_param <- toupper(method)
  
  mo_genus <- as.mo(mo_genus(mo))
  mo_family <- as.mo(mo_family(mo))
  mo_order <- as.mo(mo_order(mo))
  mo_becker <- as.mo(mo, Becker = TRUE)
  mo_lancefield <- as.mo(mo, Lancefield = TRUE)
  mo_other <- as.mo("other")
  
  guideline_coerced <- get_guideline(guideline)
  if (guideline_coerced != guideline) {
    message(blue(paste0("Note: Using guideline ", bold(guideline_coerced), " as input for `guideline`.")))
  }
  
  new_rsi <- rep(NA_character_, length(x))
  ab_param <- ab
  trans <- rsi_translation %>%
    filter(guideline == guideline_coerced & method == method_param & ab == ab_param) %>%
    mutate(lookup = paste(mo, ab))
  
  lookup_mo <- paste(mo, ab)
  lookup_genus <- paste(mo_genus, ab)
  lookup_family <- paste(mo_family, ab)
  lookup_order <- paste(mo_order, ab)
  lookup_becker <- paste(mo_becker, ab)
  lookup_lancefield <- paste(mo_lancefield, ab)
  lookup_other <- paste(mo_other, ab)
  
  if (all(trans$uti == TRUE, na.rm = TRUE) & all(uti == FALSE)) {
    message(red("WARNING."))
    warning("Interpretation of ", bold(ab_name(ab, tolower = TRUE)), " for some microorganisms is only available for (uncomplicated) urinary tract infections (UTI).\n  Use parameter 'uti' to set which isolates are from urine. See ?as.rsi.", call. = FALSE)
    warned <- TRUE
  }
  
  for (i in seq_len(length(x))) {
    get_record <- trans %>%
      # no UTI for now
      filter(lookup %in% c(lookup_mo[i],
                           lookup_genus[i],
                           lookup_family[i],
                           lookup_order[i],
                           lookup_becker[i],
                           lookup_lancefield[i],
                           lookup_other[i]))
    
    if (isTRUE(uti[i])) {
      get_record <- get_record %>% 
        # be as specific as possible (i.e. prefer species over genus):
        # desc(uti) = TRUE on top and FALSE on bottom
        arrange(desc(uti), desc(nchar(mo))) %>% # 'uti' is a column in rsi_translation
        .[1L, ]
    } else {
      get_record <- get_record %>% 
        filter(uti == FALSE) %>% # 'uti' is a column in rsi_translation
        arrange(desc(nchar(mo))) %>%
        .[1L, ]
    }
    
    if (NROW(get_record) > 0) {
      if (is.na(x[i])) {
        new_rsi[i] <- NA_character_
      } else if (method == "mic") {
        mic_input <- x[i]
        mic_S <- as.mic(get_record$breakpoint_S)
        mic_R <- as.mic(get_record$breakpoint_R)
        new_rsi[i] <- case_when(isTRUE(which(levels(mic_input) == mic_input) <= which(levels(mic_S) == mic_S)) ~ "S",
                                isTRUE(which(levels(mic_input) == mic_input) >= which(levels(mic_R) == mic_R)) ~ "R",
                                !is.na(get_record$breakpoint_S) & !is.na(get_record$breakpoint_R) ~ "I",
                                TRUE ~ NA_character_)
      } else if (method == "disk") {
        new_rsi[i] <- case_when(isTRUE(as.double(x[i]) >= as.double(get_record$breakpoint_S)) ~ "S",
                                isTRUE(as.double(x[i]) <= as.double(get_record$breakpoint_R)) ~ "R",
                                !is.na(get_record$breakpoint_S) & !is.na(get_record$breakpoint_R) ~ "I",
                                TRUE ~ NA_character_)
      }
    }
  }
  if (warned == FALSE) {
    message(green("OK."))
  }
  structure(.Data = factor(new_rsi, levels = c("S", "I", "R"), ordered = TRUE),
            class =  c("rsi", "ordered", "factor"))
}

#' @rdname as.rsi
#' @export
is.rsi <- function(x) {
  inherits(x, "rsi")
}

#' @rdname as.rsi
#' @export
is.rsi.eligible <- function(x, threshold = 0.05) {
  if (NCOL(x) > 1) {
    stop("`x` must be a one-dimensional vector.")
  }
  if (any(c("logical",
            "numeric",
            "integer",
            "mo",
            "Date",
            "POSIXct",
            "rsi",
            "raw",
            "hms")
          %in% class(x))) {
    # no transformation needed
    FALSE
  } else {
    x <- x[!is.na(x) & !is.null(x) & !identical(x, "")]
    if (length(x) == 0) {
      return(FALSE)
    }
    checked <- suppressWarnings(as.rsi(x))
    outcome <- sum(is.na(checked)) / length(x)
    outcome <= threshold
  }
}

#' @exportMethod print.rsi
#' @export
#' @importFrom dplyr %>%
#' @noRd
print.rsi <- function(x, ...) {
  cat("Class 'rsi'\n")
  print(as.character(x), quote = FALSE)
}

#' @exportMethod droplevels.rsi
#' @export
#' @noRd
droplevels.rsi <- function(x, exclude = if (anyNA(levels(x))) NULL else NA, ...) {
  x <- droplevels.factor(x, exclude = exclude, ...)
  class(x) <- c("rsi", "ordered", "factor")
  x
}

#' @exportMethod summary.rsi
#' @export
#' @noRd
summary.rsi <- function(object, ...) {
  x <- object
  c(
    "Class" = "rsi",
    "<NA>" = sum(is.na(x)),
    "Sum S" = sum(x == "S", na.rm = TRUE),
    "Sum IR" = sum(x %in% c("I", "R"), na.rm = TRUE),
    "-Sum R" = sum(x == "R", na.rm = TRUE),
    "-Sum I" = sum(x == "I", na.rm = TRUE)
  )
}

#' @exportMethod plot.rsi
#' @export
#' @importFrom dplyr %>% group_by summarise filter mutate if_else n_distinct
#' @importFrom graphics plot text
#' @noRd
plot.rsi <- function(x,
                     lwd = 2,
                     ylim = NULL,
                     ylab = "Percentage",
                     xlab = "Antimicrobial Interpretation",
                     main = paste("Susceptibility Analysis of", deparse(substitute(x))),
                     axes = FALSE,
                     ...) {
  suppressWarnings(
    data <- data.frame(x = x,
                       y = 1,
                       stringsAsFactors = TRUE) %>%
      group_by(x) %>%
      summarise(n = sum(y)) %>%
      filter(!is.na(x)) %>%
      mutate(s = round((n / sum(n)) * 100, 1))
  )
  if (!"S" %in% data$x) {
    data <- rbind(data, data.frame(x = "S", n = 0, s = 0))
  }
  if (!"I" %in% data$x) {
    data <- rbind(data, data.frame(x = "I", n = 0, s = 0))
  }
  if (!"R" %in% data$x) {
    data <- rbind(data, data.frame(x = "R", n = 0, s = 0))
  }
  
  data$x <- factor(data$x, levels = c("S", "I", "R"), ordered = TRUE)
  
  ymax <- if_else(max(data$s) > 95, 105, 100)
  
  plot(x = data$x,
       y = data$s,
       lwd = lwd,
       ylim = c(0, ymax),
       ylab = ylab,
       xlab = xlab,
       main = main,
       axes = axes,
       ...)
  # x axis
  axis(side = 1, at = 1:n_distinct(data$x), labels = levels(data$x), lwd = 0)
  # y axis, 0-100%
  axis(side = 2, at = seq(0, 100, 5))
  
  text(x = data$x,
       y = data$s + 4,
       labels = paste0(data$s, "% (n = ", data$n, ")"))
}


#' @exportMethod barplot.rsi
#' @export
#' @importFrom dplyr %>% group_by summarise
#' @importFrom graphics barplot axis par
#' @noRd
barplot.rsi <- function(height,
                        col = c("chartreuse4", "chartreuse3", "brown3"),
                        xlab = ifelse(beside, "Antimicrobial Interpretation", ""),
                        main = paste("Antimicrobial resistance of", deparse(substitute(height))),
                        ylab = "Frequency",
                        beside = TRUE,
                        axes = beside,
                        ...) {
  
  if (axes == TRUE) {
    par(mar =  c(5, 4, 4, 2) + 0.1)
  } else {
    par(mar =  c(2, 4, 4, 2) + 0.1)
  }
  
  barplot(as.matrix(table(height)),
          col = col,
          xlab = xlab,
          main = main,
          ylab = ylab,
          beside = beside,
          axes = FALSE,
          ...)
  # y axis, 0-100%
  axis(side = 2, at = seq(0, max(table(height)) + max(table(height)) * 1.1, by = 25))
  if (axes == TRUE && beside == TRUE) {
    axis(side = 1, labels = levels(height), at = c(1, 2, 3) + 0.5, lwd = 0)
  }
}

#' @importFrom vctrs vec_ptype_abbr
#' @export
vec_ptype_abbr.rsi <- function(x, ...) {
  "rsi"
}

#' @importFrom vctrs vec_ptype_full
#' @export
vec_ptype_full.rsi <- function(x, ...) {
  "rsi"
}

#' @importFrom pillar pillar_shaft
#' @importFrom crayon bgGreen bgYellow bgRed black white
#' @export 
pillar_shaft.rsi <- function(x, ...) {
  out <- trimws(format(x))
  out[is.na(x)] <- pillar::style_subtle(" NA")
  out[x == "S"] <- bgGreen(white(" S "))
  out[x == "I"] <- bgYellow(black(" I "))
  out[x == "R"] <- bgRed(white(" R "))
  pillar::new_pillar_shaft_simple(out, align = "left", width = 3)
}

#' @exportMethod [<-.rsi
#' @export
#' @noRd
"[<-.rsi" <- function(i, j, ..., value) {
  value <- as.rsi(value)
  y <- NextMethod()
  attributes(y) <- attributes(i)
  y
}
#' @exportMethod [[<-.rsi
#' @export
#' @noRd
"[[<-.rsi" <- function(i, j, ..., value) {
  value <- as.rsi(value)
  y <- NextMethod()
  attributes(y) <- attributes(i)
  y
}
#' @exportMethod c.rsi
#' @export
#' @noRd
c.rsi <- function(x, ...) {
  y <- unlist(lapply(list(...), as.character))
  x <- as.character(x)
  as.rsi(c(x, y))
}
