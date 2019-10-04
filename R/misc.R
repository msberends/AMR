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

# No export, no Rd
addin_insert_in <- function() {
  rstudioapi::insertText(" %in% ")
}

# No export, no Rd
addin_insert_like <- function() {
  rstudioapi::insertText(" %like% ")
}

#' @importFrom crayon blue bold red
#' @importFrom dplyr %>% pull
search_type_in_df <- function(x, type) {
  # try to find columns based on type
  found <- NULL

  colnames(x) <- trimws(colnames(x))

  # -- mo
  if (type == "mo") {
    if ("mo" %in% lapply(x, class)) {
      found <- colnames(x)[lapply(x, class) == "mo"][1]
    } else if (any(colnames(x) %like% "^(mo|microorganism|organism|bacteria|bacterie)s?$")) {
      found <- colnames(x)[colnames(x) %like% "^(mo|microorganism|organism|bacteria|bacterie)s?$"][1]
    } else if (any(colnames(x) %like% "^(microorganism|organism|bacteria|bacterie)")) {
      found <- colnames(x)[colnames(x) %like% "^(microorganism|organism|bacteria|bacterie)"][1]
    } else if (any(colnames(x) %like% "species")) {
      found <- colnames(x)[colnames(x) %like% "species"][1]
    }

  }
  # -- key antibiotics
  if (type == "keyantibiotics") {
    if (any(colnames(x) %like% "^key.*(ab|antibiotics)")) {
      found <- colnames(x)[colnames(x) %like% "^key.*(ab|antibiotics)"][1]
    }
  }
  # -- date
  if (type == "date") {
    if (any(colnames(x) %like% "^(specimen date|specimen_date|spec_date)")) {
      # WHONET support
      found <- colnames(x)[colnames(x) %like% "^(specimen date|specimen_date|spec_date)"][1]
      if (!any(class(x %>% pull(found)) %in% c("Date", "POSIXct"))) {
        stop(red(paste0("ERROR: Found column `", bold(found), "` to be used as input for `col_", type,
                        "`, but this column contains no valid dates. Transform its values to valid dates first.")),
             call. = FALSE)
      }
    } else {
      for (i in 1:ncol(x)) {
        if (any(class(x %>% pull(i)) %in% c("Date", "POSIXct"))) {
          found <- colnames(x)[i]
          break
        }
      }
    }
  }
  # -- patient id
  if (type == "patient_id") {
    if (any(colnames(x) %like% "^(identification |patient|patid)")) {
      found <- colnames(x)[colnames(x) %like% "^(identification |patient|patid)"][1]
    }
  }
  # -- specimen
  if (type == "specimen") {
    if (any(colnames(x) %like% "(specimen type|spec_type)")) {
      found <- colnames(x)[colnames(x) %like% "(specimen type|spec_type)"][1]
    } else if (any(colnames(x) %like% "^(specimen)")) {
      found <- colnames(x)[colnames(x) %like% "^(specimen)"][1]
    }
  }

  if (!is.null(found)) {
    msg <- paste0("NOTE: Using column `", bold(found), "` as input for `col_", type, "`.")
    if (type %in% c("keyantibiotics", "specimen")) {
      msg <- paste(msg, "Use", bold(paste0("col_", type), "= FALSE"), "to prevent this.")
    }
    message(blue(msg))
  }
  found
}

stopifnot_installed_package <- function(package) {
  # no "utils::installed.packages()" since it requires non-staged install since R 3.6.0
  # https://developer.r-project.org/Blog/public/2019/02/14/staged-install/index.html
  get(".packageName", envir = asNamespace(package))
  return(invisible())
}

"%or%" <- function(x, y) {
  if (is.null(x) | is.null(y)) {
    if (is.null(x)) {
      return(y)
    } else {
      return(x)
    }
  }
  ifelse(!is.na(x),
         x,
         ifelse(!is.na(y), y, NA))
}

class_integrity_check <- function(value, type, check_vector) {
  if (!all(value[!is.na(value)] %in% check_vector)) {
    warning(paste0("invalid ", type, ", NA generated"), call. = FALSE)
    value[!value %in% check_vector] <- NA
  }
  value
}




# Percentages -------------------------------------------------------------
# Can all be removed when clean 1.2.0 is on CRAN

getdecimalplaces <- function(x, minimum = 0, maximum = 3) {
  if (maximum < minimum) {
    maximum <- minimum
  }
  if (minimum > maximum) {
    minimum <- maximum
  }
  max_places <- max(unlist(lapply(strsplit(sub('0+$', '', 
                                               as.character(x * 100)), ".", fixed = TRUE),
                                  function(y) ifelse(length(y) == 2, nchar(y[2]), 0))), na.rm = TRUE)
  max(min(max_places,
          maximum, na.rm = TRUE),
      minimum, na.rm = TRUE)
}

round2 <- function(x, digits = 0, force_zero = TRUE) {
  # https://stackoverflow.com/a/12688836/4575331
  val <- (trunc((abs(x) * 10 ^ digits) + 0.5) / 10 ^ digits) * sign(x)
  if (digits > 0 & force_zero == TRUE) {
    val[val != as.integer(val)] <- paste0(val[val != as.integer(val)],
                                          strrep("0", max(0, digits - nchar(gsub(".*[.](.*)$", "\\1", val[val != as.integer(val)])))))
  }
  val
}

percentage <- function(x, digits = NULL, ...) {
  if (is.null(digits)) {
    digits <- getdecimalplaces(x, minimum = 0, maximum = 1)
  }
  # round right: percentage(0.4455) should return "44.6%", not "44.5%"
  x <- as.numeric(round2(x, digits = digits + 2))
  x_formatted <- format(as.double(x) * 100, scientific = FALSE, digits = digits, nsmall = digits, ...)
  x_formatted[!is.na(x)] <- paste0(x_formatted[!is.na(x)], "%")
  x_formatted
}
