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

#' Get ATC properties from WHOCC website
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
#' @rdname atc_online
#' @importFrom dplyr %>% progress_estimated
#' @inheritSection AMR Read more on our website!
#' @source \url{https://www.whocc.no/atc_ddd_alterations__cumulative/ddd_alterations/abbrevations/}
#' @examples
#' \donttest{
#' # oral DDD (Defined Daily Dose) of amoxicillin
#' atc_online_property("J01CA04", "DDD", "O")
#' # parenteral DDD (Defined Daily Dose) of amoxicillin
#' atc_online_property("J01CA04", "DDD", "P")
#'
#' atc_online_property("J01CA04", property = "groups") # search hierarchical groups of amoxicillin
#' # [1] "ANTIINFECTIVES FOR SYSTEMIC USE"
#' # [2] "ANTIBACTERIALS FOR SYSTEMIC USE"
#' # [3] "BETA-LACTAM ANTIBACTERIALS, PENICILLINS"
#' # [4] "Penicillins with extended spectrum"
#' }
atc_online_property <- function(atc_code,
                                property,
                                administration = "O",
                                url = "https://www.whocc.no/atc_ddd_index/?code=%s&showdescription=no") {

  if (!all(c("curl", "rvest", "xml2") %in% rownames(utils::installed.packages()))) {
    stop("Packages 'xml2', 'rvest' and 'curl' are required for this function")
  }

  if (!all(atc_code %in% AMR::antibiotics)) {
    atc_code <- as.character(ab_atc(atc_code))
  }

  if (!curl::has_internet()) {
    message("There appears to be no internet connection.")
    return(rep(NA, length(atc_code)))
  }

  if (length(property) != 1L) {
    stop("`property` must be of length 1", call. = FALSE)
  }
  if (length(administration) != 1L) {
    stop("`administration` must be of length 1", call. = FALSE)
  }

  # also allow unit as property
  if (property %like% "unit") {
    property <- "U"
  }

  # validation of properties
  valid_properties <- c("ATC", "Name", "DDD", "U", "Adm.R", "Note", "groups")
  valid_properties.bak <- valid_properties

  property <- tolower(property)
  valid_properties <- tolower(valid_properties)

  if (!property %in% valid_properties) {
    stop("Invalid `property`, use one of ", paste(valid_properties.bak, collapse = ", "), ".")
  }

  if (property == "ddd") {
    returnvalue <- rep(NA_real_, length(atc_code))
  } else if (property == "groups") {
    returnvalue <- list()
  } else {
    returnvalue <- rep(NA_character_, length(atc_code))
  }

  progress <- progress_estimated(n = length(atc_code))

  for (i in seq_len(length(atc_code))) {

    progress$tick()$print()

    atc_url <- sub("%s", atc_code[i], url, fixed = TRUE)

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
      texts <- texts[seq_len(length(texts)) - 1]
      returnvalue <- c(list(texts), returnvalue)

    } else {
      tbl <- xml2::read_html(atc_url) %>%
        rvest::html_nodes("table") %>%
        rvest::html_table(header = TRUE) %>%
        as.data.frame(stringsAsFactors = FALSE)

      # case insensitive column names
      colnames(tbl) <- tolower(colnames(tbl)) %>% gsub("^atc.*", "atc", .)

      if (length(tbl) == 0) {
        warning("ATC not found: ", atc_code[i], ". Please check ", atc_url, ".", call. = FALSE)
        returnvalue[i] <- NA
        next
      }

      if (property %in% c("atc", "name")) {
        # ATC and name are only in first row
        returnvalue[i] <- tbl[1, property]
      } else {
        if (!"adm.r" %in% colnames(tbl) | is.na(tbl[1, "adm.r"])) {
          returnvalue[i] <- NA
          next
        } else {
          for (j in seq_len(nrow(tbl))) {
            if (tbl[j, "adm.r"] == administration) {
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

#' @rdname atc_online
#' @export
atc_online_groups <- function(atc_code, ...) {
  atc_online_property(atc_code = atc_code, property = "groups", ...)
}

#' @rdname atc_online
#' @export
atc_online_ddd <- function(atc_code, ...) {
  atc_online_property(atc_code = atc_code, property = "ddd", ...)
}
