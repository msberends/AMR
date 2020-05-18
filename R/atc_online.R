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

#' Get ATC properties from WHOCC website
#'
#' @inheritSection lifecycle Questioning lifecycle
#' @description Gets data from the WHO to determine properties of an ATC (e.g. an antibiotic) like name, defined daily dose (DDD) or standard unit.
#' 
#' **This function requires an internet connection.**
#' @param atc_code a character or character vector with ATC code(s) of antibiotic(s)
#' @param property property of an ATC code. Valid values are `"ATC"`, `"Name"`, `"DDD"`, `"U"` (`"unit"`), `"Adm.R"`, `"Note"` and `groups`. For this last option, all hierarchical groups of an ATC code will be returned, see Examples.
#' @param administration type of administration when using `property = "Adm.R"`, see Details
#' @param url url of website of the WHO. The sign `%s` can be used as a placeholder for ATC codes.
#' @param ... parameters to pass on to `atc_property`
#' @details
#' Options for parameter `administration`:
#' 
#' - `"Implant"` = Implant
#' - `"Inhal"` = Inhalation
#' - `"Instill"` = Instillation
#' - `"N"` = nasal
#' - `"O"` = oral
#' - `"P"` = parenteral
#' - `"R"` = rectal
#' - `"SL"` = sublingual/buccal
#' - `"TD"` = transdermal
#' - `"V"` = vaginal
#'
#' Abbreviations of return values when using `property = "U"` (unit):
#' 
#' - `"g"` = gram
#' - `"mg"` = milligram
#' - `"mcg"`` = microgram
#' - `"U"` = unit
#' - `"TU"` = thousand units
#' - `"MU"` = million units
#' - `"mmol"` = millimole
#' - `"ml"` = milliliter (e.g. eyedrops)
#' @export
#' @rdname atc_online
#' @inheritSection AMR Read more on our website!
#' @source <https://www.whocc.no/atc_ddd_alterations__cumulative/ddd_alterations/abbrevations/>
#' @examples
#' \dontrun{
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
  
  stopifnot_installed_package(c("curl", "rvest", "xml2"))
  has_internet <- get("has_internet", envir = asNamespace("curl"))
  html_attr <- get("html_attr", envir = asNamespace("rvest"))
  html_children <- get("html_children", envir = asNamespace("rvest"))
  html_node <- get("html_node", envir = asNamespace("rvest"))
  html_nodes <- get("html_nodes", envir = asNamespace("rvest"))
  html_table <- get("html_table", envir = asNamespace("rvest"))
  html_text <- get("html_text", envir = asNamespace("rvest"))
  read_html <- get("read_html", envir = asNamespace("xml2"))

  check_dataset_integrity()
  
  if (!all(atc_code %in% antibiotics)) {
    atc_code <- as.character(ab_atc(atc_code))
  }
  
  if (!has_internet()) {
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

    progress$tick()

    atc_url <- sub("%s", atc_code[i], url, fixed = TRUE)

    if (property == "groups") {
      tbl <- read_html(atc_url) %>%
        html_node("#content") %>%
        html_children() %>%
        html_node("a")

      # get URLS of items
      hrefs <- tbl %>% html_attr("href")
      # get text of items
      texts <- tbl %>% html_text()
      # select only text items where URL like "code="
      texts <- texts[grepl("?code=", tolower(hrefs), fixed = TRUE)]
      # last one is antibiotics, skip it
      texts <- texts[seq_len(length(texts)) - 1]
      returnvalue <- c(list(texts), returnvalue)

    } else {
      tbl <- read_html(atc_url) %>%
        html_nodes("table") %>%
        html_table(header = TRUE) %>%
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
