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
# how to conduct AMR data analysis: https://amr-for-r.org              #
# ==================================================================== #

#' Get ATC Properties from WHOCC Website
#'
#' Gets data from the WHOCC website to determine properties of an Anatomical Therapeutic Chemical (ATC) (e.g. an antimicrobial), such as the name, defined daily dose (DDD) or standard unit.
#' @param atc_code A [character] (vector) with ATC code(s) of antimicrobials, will be coerced with [as.ab()] and [ab_atc()] internally if not a valid ATC code.
#' @param property Property of an ATC code. Valid values are `"ATC"`, `"Name"`, `"DDD"`, `"U"` (`"unit"`), `"Adm.R"`, `"Note"` and `groups`. For this last option, all hierarchical groups of an ATC code will be returned, see *Examples*.
#' @param administration Type of administration when using `property = "Adm.R"`, see *Details*.
#' @param url URL of website of the WHOCC. The sign `%s` can be used as a placeholder for ATC codes.
#' @param url_vet URL of website of the WHOCC for veterinary medicine. The sign `%s` can be used as a placeholder for ATC_vet codes (that all start with "Q").
#' @param ... Arguments to pass on to `atc_property`.
#' @details
#' Options for argument `administration`:
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
#' - `"mcg"` = microgram
#' - `"U"` = unit
#' - `"TU"` = thousand units
#' - `"MU"` = million units
#' - `"mmol"` = millimole
#' - `"ml"` = millilitre (e.g. eyedrops)
#'
#' **N.B. This function requires an internet connection and only works if the following packages are installed: `curl`, `rvest`, `xml2`.**
#' @export
#' @rdname atc_online
#' @source <https://atcddd.fhi.no/atc_ddd_alterations__cumulative/ddd_alterations/abbrevations/>
#' @examples
#' \donttest{
#' if (requireNamespace("curl") && requireNamespace("rvest") && requireNamespace("xml2")) {
#'   # oral DDD (Defined Daily Dose) of amoxicillin
#'   atc_online_property("J01CA04", "DDD", "O")
#'   atc_online_ddd(ab_atc("amox"))
#'
#'   # parenteral DDD (Defined Daily Dose) of amoxicillin
#'   atc_online_property("J01CA04", "DDD", "P")
#'
#'   atc_online_property("J01CA04", property = "groups") # search hierarchical groups of amoxicillin
#' }
#' }
atc_online_property <- function(atc_code,
                                property,
                                administration = "O",
                                url = "https://atcddd.fhi.no/atc_ddd_index/?code=%s&showdescription=no",
                                url_vet = "https://atcddd.fhi.no/atcvet/atcvet_index/?code=%s&showdescription=no") {
  meet_criteria(atc_code, allow_class = "character", allow_NA = TRUE)
  meet_criteria(property, allow_class = "character", has_length = 1, is_in = c("ATC", "Name", "DDD", "U", "unit", "Adm.R", "Note", "groups"), ignore.case = TRUE)
  meet_criteria(administration, allow_class = "character", has_length = 1)
  meet_criteria(url, allow_class = "character", has_length = 1, looks_like = "https?://")
  meet_criteria(url_vet, allow_class = "character", has_length = 1, looks_like = "https?://")

  has_internet <- import_fn("has_internet", "curl")
  html_attr <- import_fn("html_attr", "rvest")
  html_children <- import_fn("html_children", "rvest")
  html_node <- import_fn("html_node", "rvest")
  html_nodes <- import_fn("html_nodes", "rvest")
  html_table <- import_fn("html_table", "rvest")
  html_text <- import_fn("html_text", "rvest")
  read_html <- import_fn("read_html", "xml2")

  if (!all(atc_code %in% unlist(AMR::antimicrobials$atc))) {
    missing <- atc_code %unlike% "[A-Z][0-9][0-9][A-Z][A-Z][0-9][0-9]"
    atc_code[missing] <- as.character(ab_atc(atc_code[missing], only_first = TRUE))
  }

  if (!has_internet()) {
    message_("There appears to be no internet connection, returning NA.",
      add_fn = font_red,
      as_note = FALSE
    )
    return(rep(NA, length(atc_code)))
  }

  property <- tolower(property)
  # also allow unit as property
  if (property == "unit") {
    property <- "u"
  }
  if (property == "ddd") {
    returnvalue <- rep(NA_real_, length(atc_code))
  } else if (property == "groups") {
    returnvalue <- list()
  } else {
    returnvalue <- rep(NA_character_, length(atc_code))
  }

  progress <- progress_ticker(n = length(atc_code), 3)
  on.exit(close(progress))

  for (i in seq_len(length(atc_code))) {
    progress$tick()

    if (is.na(atc_code[i])) {
      next
    }

    if (atc_code[i] %like% "^Q") {
      # veterinary drugs, ATC_vet codes start with a "Q"
      atc_url <- url_vet
    } else {
      atc_url <- url
    }
    atc_url <- sub("%s", atc_code[i], atc_url, fixed = TRUE)

    if (property == "groups") {
      out <- tryCatch(
        read_html(atc_url) %pm>%
          html_node("#content") %pm>%
          html_children() %pm>%
          html_node("a"),
        error = function(e) NULL
      )
      if (is.null(out)) {
        message_("Connection to ", atc_url, " failed.")
        return(rep(NA, length(atc_code)))
      }

      # get URLS of items
      hrefs <- out %pm>% html_attr("href")
      # get text of items
      texts <- out %pm>% html_text()
      # select only text items where URL like "code="
      texts <- texts[grepl("?code=", tolower(hrefs), fixed = TRUE)]
      # last one is antibiotics, skip it
      texts <- texts[seq_len(length(texts)) - 1]
      returnvalue <- c(list(texts), returnvalue)
    } else {
      out <- tryCatch(
        read_html(atc_url) %pm>%
          html_nodes("table") %pm>%
          html_table(header = TRUE) %pm>%
          as.data.frame(stringsAsFactors = FALSE),
        error = function(e) NULL
      )
      if (is.null(out)) {
        message_("Connection to ", atc_url, " failed.")
        return(rep(NA, length(atc_code)))
      }

      # case insensitive column names
      colnames(out) <- gsub("^atc.*", "atc", tolower(colnames(out)))

      if (length(out) == 0) {
        message_("in `atc_online_property()`: no properties found for ATC ", atc_code[i], ". Please check ", font_url(atc_url, "this WHOCC webpage"), ".")
        returnvalue[i] <- NA
        next
      }

      if (property %in% c("atc", "name")) {
        # ATC and name are only in first row
        returnvalue[i] <- out[1, property, drop = TRUE]
      } else {
        if (!"adm.r" %in% colnames(out) || is.na(out[1, "adm.r", drop = TRUE])) {
          returnvalue[i] <- NA
          next
        } else {
          for (j in seq_len(nrow(out))) {
            if (out[j, "adm.r"] == administration) {
              returnvalue[i] <- out[j, property, drop = TRUE]
            }
          }
        }
      }
    }
  }

  if (property == "groups" && length(returnvalue) == 1) {
    returnvalue <- returnvalue[[1]]
  }

  returnvalue
}

#' @rdname atc_online
#' @export
atc_online_groups <- function(atc_code, ...) {
  meet_criteria(atc_code, allow_class = "character", allow_NA = TRUE)
  atc_online_property(atc_code = atc_code, property = "groups", ...)
}

#' @rdname atc_online
#' @export
atc_online_ddd <- function(atc_code, ...) {
  meet_criteria(atc_code, allow_class = "character", allow_NA = TRUE)
  atc_online_property(atc_code = atc_code, property = "ddd", ...)
}

#' @rdname atc_online
#' @export
atc_online_ddd_units <- function(atc_code, ...) {
  meet_criteria(atc_code, allow_class = "character", allow_NA = TRUE)
  atc_online_property(atc_code = atc_code, property = "unit", ...)
}
