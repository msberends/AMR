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


#' Transform to ATC code
#'
#' Use this function to determine the ATC code of one or more antibiotics. The data set \code{\link{antibiotics}} will be searched for abbreviations, official names and trade names.
#' @param x character vector to determine \code{ATC} code
#' @rdname as.atc
#' @aliases atc
#' @keywords atc
#' @export
#' @importFrom dplyr %>% filter slice pull
#' @details Use the \code{\link{ab_property}} functions to get properties based on the returned ATC code, see Examples.
#'
#' In the ATC classification system, the active substances are classified in a hierarchy with five different levels.  The system has fourteen main anatomical/pharmacological groups or 1st levels. Each ATC main group is divided into 2nd levels which could be either pharmacological or therapeutic groups.  The 3rd and 4th levels are chemical, pharmacological or therapeutic subgroups and the 5th level is the chemical substance.  The 2nd, 3rd and 4th levels are often used to identify pharmacological subgroups when that is considered more appropriate than therapeutic or chemical subgroups.
#'   Source: \url{https://www.whocc.no/atc/structure_and_principles/}
#' @return Character (vector) with class \code{"act"}. Unknown values will return \code{NA}.
#' @seealso \code{\link{antibiotics}} for the dataframe that is being used to determine ATCs.
#' @examples
#' # These examples all return "J01FA01", the ATC code of Erythromycin:
#' as.atc("J01FA01")
#' as.atc("Erythromycin")
#' as.atc("eryt")
#' as.atc("ERYT")
#' as.atc("ERY")
#' as.atc("Erythrocin") # Trade name
#' as.atc("Eryzole")    # Trade name
#' as.atc("Pediamycin") # Trade name
#'
#' # Use ab_* functions to get a specific property based on an ATC code
#' Cipro <- as.atc("cipro") # returns `J01MA02`
#' ab_official(Cipro)       # returns "Ciprofloxacin"
#' ab_umcg(Cipro)           # returns "CIPR", the code used in the UMCG
as.atc <- function(x) {

  x.new <- rep(NA_character_, length(x))
  x.bak <- x
  x <- unique(x[!is.na(x)])
  failures <- character(0)

  for (i in 1:length(x)) {
    fail <- TRUE

    # first try atc
    found <- AMR::antibiotics[which(AMR::antibiotics$atc == x[i]),]$atc
    if (length(found) > 0) {
      fail <- FALSE
      x.new[is.na(x.new) & x.bak == x[i]] <- found[1L]
    }

    # try ATC in code form, even if it does not exist in the antibiotics data set YET
    if (length(found) == 0 & x[i] %like% '[A-Z][0-9][0-9][A-Z][A-Z][0-9][0-9]') {
      warning("ATC code ", x[i], " is not yet in the `antibiotics` data set.")
      fail <- FALSE
      x.new[is.na(x.new) & x.bak == x[i]] <- x[i]
    }

    # try abbreviation of certe and glims
    found <- AMR::antibiotics[which(tolower(AMR::antibiotics$certe) == tolower(x[i])),]$atc
    if (length(found) > 0) {
      fail <- FALSE
      x.new[is.na(x.new) & x.bak == x[i]] <- found[1L]
    }
    found <- AMR::antibiotics[which(tolower(AMR::antibiotics$umcg) == tolower(x[i])),]$atc
    if (length(found) > 0) {
      fail <- FALSE
      x.new[is.na(x.new) & x.bak == x[i]] <- found[1L]
    }

    # try exact official name
    found <- AMR::antibiotics[which(tolower(AMR::antibiotics$official) == tolower(x[i])),]$atc
    if (length(found) > 0) {
      fail <- FALSE
      x.new[is.na(x.new) & x.bak == x[i]] <- found[1L]
    }

    # try exact official Dutch
    found <- AMR::antibiotics[which(tolower(AMR::antibiotics$official_nl) == tolower(x[i])),]$atc
    if (length(found) > 0) {
      fail <- FALSE
      x.new[is.na(x.new) & x.bak == x[i]] <- found[1L]
    }

    # try trade name
    found <- AMR::antibiotics[which(paste0("(", AMR::antibiotics$trade_name, ")") %like% x[i]),]$atc
    if (length(found) > 0) {
      fail <- FALSE
      x.new[is.na(x.new) & x.bak == x[i]] <- found[1L]
    }

    # try abbreviation
    found <- AMR::antibiotics[which(paste0("(", AMR::antibiotics$abbr, ")") %like% x[i]),]$atc
    if (length(found) > 0) {
      fail <- FALSE
      x.new[is.na(x.new) & x.bak == x[i]] <- found[1L]
    }

    # nothing helped, try first chars of official name, but only if nchar > 4 (cipro, nitro, fosfo)
    if (nchar(x[i]) > 4) {
      found <- AMR::antibiotics[which(AMR::antibiotics$official %like% paste0("^", substr(x[i], 1, 5))),]$atc
      if (length(found) > 0) {
        fail <- FALSE
        x.new[is.na(x.new) & x.bak == x[i]] <- found[1L]
      }
    }

    # not found
    if (fail == TRUE) {
      failures <- c(failures, x[i])
    }
  }

  failures <- failures[!failures %in% c(NA, NULL, NaN)]
  if (length(failures) > 0) {
    warning("These values could not be coerced to a valid atc: ",
            paste('"', unique(failures), '"', sep = "", collapse = ', '),
            ".",
            call. = FALSE)
  }
  class(x.new) <- "atc"
  attr(x.new, 'package') <- 'AMR'
  x.new
}

#' @rdname as.atc
#' @export
guess_atc <- as.atc

#' @rdname as.atc
#' @export
is.atc <- function(x) {
  identical(class(x), "atc")
}


#' @exportMethod print.atc
#' @export
#' @noRd
print.atc <- function(x, ...) {
  cat("Class 'atc'\n")
  print.default(as.character(x), quote = FALSE)
}

#' @exportMethod as.data.frame.atc
#' @export
#' @noRd
as.data.frame.atc <- function (x, ...) {
  # same as as.data.frame.character but with removed stringsAsFactors
  nm <- paste(deparse(substitute(x), width.cutoff = 500L),
              collapse = " ")
  if (!"nm" %in% names(list(...))) {
    as.data.frame.vector(x, ..., nm = nm)
  } else {
    as.data.frame.vector(x, ...)
  }
}

#' @exportMethod pull.atc
#' @export
#' @importFrom dplyr pull
#' @noRd
pull.atc <- function(.data, ...) {
  pull(as.data.frame(.data), ...)
}

#' Properties of an ATC code
#'
#' Gets data from the WHO to determine properties of an ATC (e.g. an antibiotic) like name, defined daily dose (DDD) or standard unit. \cr \strong{This function requires an internet connection.}
#' @param atc_code a character or character vector with ATC code(s) of antibiotic(s)
#' @param property property of an ATC code. Valid values are \code{"ATC"}, \code{"Name"}, \code{"DDD"}, \code{"U"} (\code{"unit"}), \code{"Adm.R"}, \code{"Note"} and \code{groups}. For this last option, all hierarchical groups of an ATC code will be returned, see Examples.
#' @param administration type of administration when using \code{property = "Adm.R"}, see Details
#' @param url url of website of the WHO. The sign \code{\%s} can be used as a placeholder for ATC codes.
#' @param ... parameters to pass on to \code{atc_property}
#' @details
#' Options for parameter \code{administration}:
#' \itemize{
#'   \item{\code{"Implant"}}{ = Implant}
#'   \item{\code{"Inhal"}}{ = Inhalation}
#'   \item{\code{"Instill"}}{ = Instillation}
#'   \item{\code{"N"}}{ = nasal}
#'   \item{\code{"O"}}{ = oral}
#'   \item{\code{"P"}}{ = parenteral}
#'   \item{\code{"R"}}{ = rectal}
#'   \item{\code{"SL"}}{ = sublingual/buccal}
#'   \item{\code{"TD"}}{ = transdermal}
#'   \item{\code{"V"}}{ = vaginal}
#' }
#'
#' Abbreviations of return values when using \code{property = "U"} (unit):
#' \itemize{
#'   \item{\code{"g"}}{ = gram}
#'   \item{\code{"mg"}}{ = milligram}
#'   \item{\code{"mcg"}}{ = microgram}
#'   \item{\code{"U"}}{ = unit}
#'   \item{\code{"TU"}}{ = thousand units}
#'   \item{\code{"MU"}}{ = million units}
#'   \item{\code{"mmol"}}{ = millimole}
#'   \item{\code{"ml"}}{ = milliliter (e.g. eyedrops)}
#' }
#' @export
#' @rdname atc_property
#' @importFrom dplyr %>% progress_estimated
#' @importFrom xml2 read_html
#' @importFrom rvest html_children html_node html_nodes html_table
#' @importFrom curl nslookup
#' @source \url{https://www.whocc.no/atc_ddd_alterations__cumulative/ddd_alterations/abbrevations/}
#' @examples
#' \donttest{
#' # What's the ATC of amoxicillin?
#' guess_atc("Amoxicillin")
#' # [1] "J01CA04"
#'
#' # oral DDD (Defined Daily Dose) of amoxicillin
#' atc_property("J01CA04", "DDD", "O")
#' # parenteral DDD (Defined Daily Dose) of amoxicillin
#' atc_property("J01CA04", "DDD", "P")
#'
#' atc_property("J01CA04", property = "groups") # search hierarchical groups of amoxicillin
#' # [1] "ANTIINFECTIVES FOR SYSTEMIC USE"
#' # [2] "ANTIBACTERIALS FOR SYSTEMIC USE"
#' # [3] "BETA-LACTAM ANTIBACTERIALS, PENICILLINS"
#' # [4] "Penicillins with extended spectrum"
#' }
atc_property <- function(atc_code,
                         property,
                         administration = 'O',
                         url = 'https://www.whocc.no/atc_ddd_index/?code=%s&showdescription=no') {

  # check active network interface, from https://stackoverflow.com/a/5078002/4575331
  has_internet <- function(url) {
    # extract host from given url
    # https://www.whocc.no/atc_ddd_index/ -> www.whocc.no
    url <- url %>%
      gsub("^(http://|https://)", "", .) %>%
      strsplit('/', fixed = TRUE) %>%
      unlist() %>%
      .[1]
    !is.null(curl::nslookup(url, error = FALSE))
  }
  # check for connection using the ATC of amoxicillin
  if (!has_internet(url = url)) {
    message("The URL could not be reached.")
    return(rep(NA, length(atc_code)))
  }

  if (length(property) != 1L) {
    stop('`property` must be of length 1', call. = FALSE)
  }
  if (length(administration) != 1L) {
    stop('`administration` must be of length 1', call. = FALSE)
  }

  # also allow unit as property
  if (property %like% 'unit') {
    property <- 'U'
  }

  # validation of properties
  valid_properties <- c("ATC", "Name", "DDD", "U", "Adm.R", "Note", "groups")
  valid_properties.bak <- valid_properties

  property <- tolower(property)
  valid_properties <- tolower(valid_properties)

  if (!property %in% valid_properties) {
    stop('Invalid `property`, use one of ', paste(valid_properties.bak, collapse = ", "), '.')
  }

  if (property == 'ddd') {
    returnvalue <- rep(NA_real_, length(atc_code))
  } else if (property == 'groups') {
    returnvalue <- list()
  } else {
    returnvalue <- rep(NA_character_, length(atc_code))
  }

  progress <- progress_estimated(n = length(atc_code))

  for (i in 1:length(atc_code)) {

    progress$tick()$print()

    atc_url <- sub('%s', atc_code[i], url, fixed = TRUE)

    if (property == "groups") {
      tbl <- xml2::read_html(atc_url) %>%
        rvest::html_node("#content") %>%
        rvest::html_children() %>%
        rvest::html_node("a")

      # get URLS of items
      hrefs <- tbl %>% rvest::html_attr("href")
      # get text of items
      texts <- tbl %>% rvest::html_text()
      # select only text items where URL like "code="
      texts <- texts[grepl("?code=", tolower(hrefs), fixed = TRUE)]
      # last one is antibiotics, skip it
      texts <- texts[1:length(texts) - 1]
      returnvalue <- c(list(texts), returnvalue)

    } else {
      tbl <- xml2::read_html(atc_url) %>%
        rvest::html_nodes('table') %>%
        rvest::html_table(header = TRUE) %>%
        as.data.frame(stringsAsFactors = FALSE)

      # case insensitive column names
      colnames(tbl) <- tolower(colnames(tbl)) %>% gsub('^atc.*', 'atc', .)

      if (length(tbl) == 0) {
        warning('ATC not found: ', atc_code[i], '. Please check ', atc_url, '.', call. = FALSE)
        returnvalue[i] <- NA
        next
      }

      if (property %in% c('atc', 'name')) {
        # ATC and name are only in first row
        returnvalue[i] <- tbl[1, property]
      } else {
        if (!'adm.r' %in% colnames(tbl) | is.na(tbl[1, 'adm.r'])) {
          returnvalue[i] <- NA
          next
        } else {
          for (j in 1:nrow(tbl)) {
            if (tbl[j, 'adm.r'] == administration) {
              returnvalue[i] <- tbl[j, property]
            }
          }
        }
      }
    }
  }

  if (property == "groups" & length(returnvalue) == 1) {
    returnvalue <- returnvalue[[1]]
  }

  returnvalue
}

#' @rdname atc_property
#' @export
atc_groups <- function(atc_code, ...) {
  atc_property(atc_code = atc_code, property = "groups", ...)
}

#' @rdname atc_property
#' @export
atc_ddd <- function(atc_code, ...) {
  atc_property(atc_code = atc_code, property = "ddd", ...)
}

