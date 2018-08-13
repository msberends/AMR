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


#' Find ATC code based on antibiotic property
#'
#' Use this function to determine the ATC code of one or more antibiotics. The dataset \code{\link{antibiotics}} will be searched for abbreviations, official names and trade names.
#' @param x character vector to determine \code{ATC} code
#' @export
#' @importFrom dplyr %>% filter slice pull
#' @details In the ATC classification system, the active substances are classified in a hierarchy with five different levels.  The system has fourteen main anatomical/pharmacological groups or 1st levels.  Each ATC main group is divided into 2nd levels which could be either pharmacological or therapeutic groups.  The 3rd and 4th levels are chemical, pharmacological or therapeutic subgroups and the 5th level is the chemical substance.  The 2nd, 3rd and 4th levels are often used to identify pharmacological subgroups when that is considered more appropriate than therapeutic or chemical subgroups.
#'   Source: \url{https://www.whocc.no/atc/structure_and_principles/}
#' @return Character (vector).
#' @seealso \code{\link{antibiotics}} for the dataframe that is being used to determine ATC's.
#' @examples
#' # These examples all return "J01FA01", the ATC code of Erythromycin:
#' guess_atc("J01FA01")
#' guess_atc("Erythromycin")
#' guess_atc("eryt")
#' guess_atc("ERYT")
#' guess_atc("ERY")
#' guess_atc("Erythrocin") # Trade name
#' guess_atc("Eryzole")    # Trade name
#' guess_atc("Pediamycin") # Trade name
guess_atc <- function(x) {

  # use this later to further fill AMR::antibiotics

  # drug <- "Ciprofloxacin"
  # url <- xml2::read_html(paste0("https://www.ncbi.nlm.nih.gov/pccompound?term=", drug)) %>%
  #   html_nodes(".rslt") %>%
  #   .[[1]] %>%
  #   html_nodes(".title a") %>%
  #   html_attr("href") %>%
  #   gsub("/compound/", "/rest/pug_view/data/compound/", ., fixed = TRUE) %>%
  #   paste0("/XML/?response_type=display")
  # synonyms <- url %>%
  #   read_xml() %>%
  #   xml_contents() %>% .[[6]] %>%
  #   xml_contents() %>% .[[8]] %>%
  #   xml_contents() %>% .[[3]] %>%
  #   xml_contents() %>% .[[3]] %>%
  #   xml_contents() %>%
  #   paste() %>%
  #   .[. %like% "StringValueList"] %>%
  #   gsub("[</]+StringValueList[>]", "", .)


  for (i in 1:length(x)) {

    # first try atc
    found <- AMR::antibiotics %>% filter(atc == x[i])

    if (nrow(found) == 0) {
      # try abbreviation of molis and glims
      found <- AMR::antibiotics %>% filter(tolower(molis) == tolower(x[i]) | tolower(umcg) == tolower(x[i]))
    }

    if (nrow(found) == 0) {
      # try exact official name
      found <- AMR::antibiotics[which(tolower(AMR::antibiotics$official) == tolower(x[i])),]
    }

    if (nrow(found) == 0) {
      # try trade name
      found <- AMR::antibiotics[which(paste0("(", AMR::antibiotics$trade_name, ")") %like% x[i]),]
    }

    if (nrow(found) == 0) {
      # try abbreviation
      found <- AMR::antibiotics[which(paste0("(", AMR::antibiotics$abbr, ")") %like% x[i]),]
    }
    # if (nrow(found) == 0) {
    #   # loosely try official name
    #   found <- AMR::antibiotics[which(AMR::antibiotics$official %like% x[i]),]
    # }

    if (nrow(found) != 0) {
      x[i] <- found %>%
        slice(1) %>%
        pull(atc)
    } else {
      x[i] <- NA
    }
  }
  x
}
