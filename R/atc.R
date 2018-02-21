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
#' Gets data from the WHO to determine properties of an ATC of e.g. an antibiotic.
#' @param atc_code a character or character vector with ATC code(s) of antibiotic(s)
#' @param property property of an ATC code. Valid values are \code{"ATC code"}, \code{"Name"}, \code{"DDD"}, \code{"U"} (\code{"unit"}), \code{"Adm.R"} en \code{"Note"}.
#' @param administration type of administration, see \emph{Details}
#' @param url url of website of the WHO. The sign \code{\%s} can be used as a placeholder for ATC codes.
#' @details
#' Abbreviations for the property \code{"Adm.R"} (parameter \code{administration}):
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
#' Abbreviations for the property \code{"U"} (unit):
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
#' @importFrom dplyr %>% progress_estimated
#' @importFrom xml2 read_html
#' @importFrom rvest html_nodes html_table
#' @source \url{https://www.whocc.no/atc_ddd_alterations__cumulative/ddd_alterations/abbrevations/}
atc_property <- function(atc_code,
                         property,
                         administration = 'O',
                         url = 'https://www.whocc.no/atc_ddd_index/?code=%s&showdescription=no') {
  
  # property <- property %>% tolower()
  #
  if (property %like% 'unit') {
    property <- 'U'
  }
  
  # validation of properties
  valid_properties.bak <- c("ATC code", "Name", "DDD", "U", "Adm.R", "Note")
  valid_properties <- valid_properties.bak #%>% tolower()
  if (!property %in% valid_properties) {
    stop('Invalid `property`, use one of ', paste(valid_properties, collapse = ", "), '.')
  }
  
  returnvalue <- rep(NA_character_, length(atc_code))
  if (property == 'DDD') {
    returnvalue <- rep(NA_real_, length(atc_code))
  }
  
  progress <- progress_estimated(n = length(atc_code))
  
  for (i in 1:length(atc_code)) {
    
    progress$tick()$print()
    
    atc_url <- sub('%s', atc_code[i], url, fixed = TRUE)
    tbl <- xml2::read_html(atc_url) %>%
      rvest::html_nodes('table') %>%
      rvest::html_table(header = TRUE)
    
    if (length(tbl) == 0) {
      warning('ATC not found: ', atc_code[i], '. Please check ', atc_url, '.', call. = FALSE)
      returnvalue[i] <- NA
      next
    }
    
    tbl <- tbl[[1]]
    
    if (property == 'Name') {
      returnvalue[i] <- tbl[1, 2]
    } else {
      
      names(returnvalue)[i] <- tbl[1, 2] %>% as.character()
      
      if (!'Adm.R' %in% colnames(tbl) | is.na(tbl[1, 'Adm.R'])) {
        returnvalue[i] <- NA
        next
      } else {
        for (j in 1:nrow(tbl)) {
          if (tbl[j, 'Adm.R'] == administration) {
            returnvalue[i] <- tbl[j, property]
          }
        }
      }
    }
  }
  
  cat('\n')
  returnvalue
  
}

#' Name of an antibiotic
#'
#' Convert antibiotic codes (from a laboratory information system like MOLIS or GLIMS) to a (trivial) antibiotic name or ATC code, or vice versa. This uses the data from \code{\link{ablist}}.
#' @param abcode a code or name, like \code{"amox"}, \code{"cftr"} or \code{"J01CA04"}
#' @param from,to type to transform from and to. See \code{\link{ablist}} for its column names.
#' @param textbetween text to put between multiple returned texts
#' @param tolower return output as lower case with function \code{\link{tolower}}.
#' @keywords ab antibiotics
#' @export
#' @importFrom dplyr %>% filter select slice 
#' @examples
#' abname("AMCL")
#' # "amoxicillin and enzyme inhibitor"
#'
#' abname("AMCL+GENT")
#' # "amoxicillin and enzyme inhibitor + gentamicin"
#'
#' abname(c("AMCL", "GENT"))
#' # "amoxicillin and enzyme inhibitor" "gentamicin" 
#'
#' abname("AMCL", to = "trivial")
#' # "Amoxicilline/clavulaanzuur"
#'
#' abname("AMCL", to = "atc")
#' # "J01CR02"
#'
#' abname("J01CR02", from = "atc", to = "umcg")
#' # "AMCL"
#'
#' @source \code{\link{ablist}}
abname <- function(abcode, from = 'umcg', to = 'official', textbetween = ' + ', tolower = FALSE) {
  
  ablist <- AMR::ablist
  colnames(ablist) <- colnames(ablist) %>% tolower()
  from <- from %>% tolower()
  to <- to %>% tolower()
  
  if (!from %in% colnames(ablist) |
      !to %in% colnames(ablist)) {
    stop(paste0('Invalid `from` or `to`. Choose one of ', 
                colnames(ablist) %>% paste(collapse = ","), '.'), call. = FALSE)
  }
  
  abcode <- as.character(abcode)
  
  for (i in 1:length(abcode)) {
    drug <- abcode[i]
    if (!grepl('+', drug, fixed = TRUE) & !grepl(' en ', drug, fixed = TRUE)) {
      # bestaat maar uit 1 middel
      if (any(ablist[, from] == drug)) {
        abcode[i] <-
          ablist %>%
          filter(.[, from] == drug) %>%
          select(to) %>%
          slice(1) %>%
          as.character()
      } else {
        # niet gevonden
        warning('Code "', drug, '" not found in antibiotics list.', call. = FALSE)
        abcode[i] <- NA
      }
    } else {
      # meerdere middelen
      if (grepl('+', drug, fixed = TRUE)) {
        drug.group <-
          strsplit(drug, '+', fixed = TRUE) %>%
          unlist() %>%
          trimws('both')
      } else if (grepl(' en ', drug, fixed = TRUE)) {
        drug.group <-
          strsplit(drug, ' en ', fixed = TRUE) %>%
          unlist() %>%
          trimws('both')
      } else {
        warning('Invalid concat.')
        abcode[i] <- NA
        next
      }
      
      for (j in 1:length(drug.group)) {
        drug.group[j] <-
          ablist %>%
          filter(.[, from] == drug.group[j]) %>%
          select(to) %>%
          slice(1) %>%
          as.character()
        if (j > 1 & to %in% c('official', 'trivial')) {
          drug.group[j] <- drug.group[j] %>% tolower()
        }
      }
      abcode[i] <- paste(drug.group, collapse = textbetween)
    }
  }
  
  if (tolower == TRUE) {
    abcode <- abcode %>% tolower()
  }
  
  abcode
}
