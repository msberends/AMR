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
#' Gets data from the WHO to determine properties of an ATC (e.g. an antibiotic) like name, defined daily dose (DDD) or standard unit. \strong{This function requires an internet connection.}
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
#' @examples
#' \donttest{
#' atc_property("J01CA04", "DDD", "O") # oral DDD (Defined Daily Dose) of amoxicillin
#' atc_property("J01CA04", "DDD", "P") # parenteral DDD (Defined Daily Dose) of amoxicillin
#' }
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
#' Convert antibiotic codes (from a laboratory information system like MOLIS or GLIMS) to a (trivial) antibiotic name or ATC code, or vice versa. This uses the data from \code{\link{antibiotics}}.
#' @param abcode a code or name, like \code{"AMOX"}, \code{"AMCL"} or \code{"J01CA04"}
#' @param from,to type to transform from and to. See \code{\link{antibiotics}} for its column names. WIth \code{from = "guess"} the from will be guessed from \code{"atc"}, \code{"molis"} and \code{"umcg"}.
#' @param textbetween text to put between multiple returned texts
#' @param tolower return output as lower case with function \code{\link{tolower}}.
#' @keywords ab antibiotics
#' @source \code{\link{antibiotics}}
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
#' abname("AMCL", to = "trivial_nl")
#' # "Amoxicilline/clavulaanzuur"
#'
#' abname("AMCL", to = "atc")
#' # "J01CR02"
#'
#' abname("J01CR02", from = "atc", to = "umcg")
#' # "AMCL"
abname <- function(abcode, from = c("guess", "atc", "molis", "umcg"), to = 'official', textbetween = ' + ', tolower = FALSE) {

  antibiotics <- AMR::antibiotics

  from <- from[1]
  if (from == "guess") {
    for (i in 1:3) {
      if (abcode[1] %in% (antibiotics %>% pull(i))) {
        from <- colnames(antibiotics)[i]
      }
    }
    if (from == "guess") {
      from <- "umcg"
    }
  }

  colnames(antibiotics) <- colnames(antibiotics) %>% tolower()
  from <- from %>% tolower()
  to <- to %>% tolower()

  if (!from %in% colnames(antibiotics) |
      !to %in% colnames(antibiotics)) {
    stop(paste0('Invalid `from` or `to`. Choose one of ',
                colnames(antibiotics) %>% paste(collapse = ","), '.'), call. = FALSE)
  }

  abcode <- as.character(abcode)

  for (i in 1:length(abcode)) {
    drug <- abcode[i]
    if (!grepl('+', drug, fixed = TRUE) & !grepl(' en ', drug, fixed = TRUE)) {
      # only 1 drug
      if (drug %in% (antibiotics %>% pull(from))) {
        abcode[i] <-
          antibiotics %>%
          filter(.[, from] == drug) %>%
          select(to) %>%
          slice(1) %>%
          as.character()
      } else {
        # not found
        warning('Code "', drug, '" not found in antibiotics list.', call. = FALSE)
        abcode[i] <- NA
      }
    } else {
      # more than 1 drug
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
          antibiotics %>%
          filter(.[, from] == drug.group[j]) %>%
          select(to) %>%
          slice(1) %>%
          as.character()
        if (j > 1 & to %in% c('official', 'trivial_nl')) {
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

#' Find bacteria ID based on genus/species
#'
#' Use this function to determine a valid ID based on a genus (and species). This input could be a full name (like \code{"Staphylococcus aureus"}), an abbreviated name (like \code{"S. aureus"}), or just a genus. You could also use a \code{\link{paste}} of a genus and species column to use the full name as input: \code{x = paste(df$genus, df$species)}, where \code{df} is your dataframe.
#' @param x character vector to determine \code{bactid}
#' @export
#' @importFrom dplyr %>% filter slice pull
#' @return Character (vector).
#' @seealso \code{\link{microorganisms}} for the dataframe that is being used to determine ID's.
#' @examples
#' # These examples all return "STAAUR", the ID of S. aureus:
#' guess_bactid("stau")
#' guess_bactid("STAU")
#' guess_bactid("staaur")
#' guess_bactid("S. aureus")
#' guess_bactid("S aureus")
#' guess_bactid("Staphylococcus aureus")
#' guess_bactid("MRSA") # Methicillin-resistant S. aureus
#' guess_bactid("VISA") # Vancomycin Intermediate S. aureus
guess_bactid <- function(x) {
  # remove dots and other non-text in case of "E. coli" except spaces
  x <- gsub("[^a-zA-Z ]+", "", x)
  # but spaces before and after should be omitted
  x <- trimws(x, which = "both")
  x.bak <- x
  # replace space by regex sign
  x <- gsub(" ", ".*", x, fixed = TRUE)
  # add start and stop
  x_species <- paste(x, 'species')
  x <- paste0('^', x, '$')

  for (i in 1:length(x)) {
    if (tolower(x[i]) == '^e.*coli$') {
      # avoid detection of Entamoeba coli in case of E. coli
      x[i] <- 'Escherichia coli'
    }
    if (tolower(x[i]) == '^h.*influenzae$') {
      # avoid detection of Haematobacter influenzae in case of H. influenzae
      x[i] <- 'Haemophilus influenzae'
    }
    if (tolower(x[i]) == '^st.*au$'
        | tolower(x[i]) == '^stau$'
        | tolower(x[i]) == '^staaur$') {
      # avoid detection of Staphylococcus auricularis in case of S. aureus
      x[i] <- 'Staphylococcus aureus'
    }
    if (tolower(x[i]) == '^p.*aer$') {
      # avoid detection of Pasteurella aerogenes in case of Pseudomonas aeruginosa
      x[i] <- 'Pseudomonas aeruginosa'
    }

    # translate known trivial names to genus+species
    if (toupper(x.bak[i]) == 'MRSA'
        | toupper(x.bak[i]) == 'VISA'
        | toupper(x.bak[i]) == 'VRSA') {
      x[i] <- 'Staphylococcus aureus'
    }
    if (toupper(x.bak[i]) == 'MRSE') {
      x[i] <- 'Staphylococcus epidermidis'
    }
    if (toupper(x.bak[i]) == 'VRE') {
      x[i] <- 'Enterococcus'
    }
    if (toupper(x.bak[i]) == 'MRPA') {
      # multi resistant P. aeruginosa
      x[i] <- 'Pseudomonas aeruginosa'
    }
    if (toupper(x.bak[i]) == 'PISP'
        | toupper(x.bak[i]) == 'PRSP') {
      # peni resistant S. pneumoniae
      x[i] <- 'Streptococcus pneumoniae'
    }
    if (toupper(x.bak[i]) == 'VISP'
        | toupper(x.bak[i]) == 'VRSP') {
      # vanco resistant S. pneumoniae
      x[i] <- 'Streptococcus pneumoniae'
    }

    # let's try the ID's first
    found <- AMR::microorganisms %>% filter(bactid == x.bak[i])

    if (nrow(found) == 0) {
      # now try exact match
      found <- AMR::microorganisms %>% filter(fullname == x[i])
    }
    if (nrow(found) == 0) {
      # try any match
      found <- AMR::microorganisms %>% filter(fullname %like% x[i])
    }
    if (nrow(found) == 0) {
      # try only genus, with 'species' attached
      found <- AMR::microorganisms %>% filter(fullname %like% x_species[i])
    }
    if (nrow(found) == 0) {
      # search for GLIMS code
      if (toupper(x.bak[i]) %in% toupper(AMR::microorganisms.umcg$mocode)) {
        found <- AMR::microorganisms.umcg %>% filter(toupper(mocode) == toupper(x.bak[i]))
      }
    }
    if (nrow(found) == 0) {
      # try splitting of characters and then find ID
      # like esco = E. coli, klpn = K. pneumoniae, stau = S. aureus
      x_split <- x
      x_length <- nchar(x.bak[i])
      x_split[i] <- paste0(x.bak[i] %>% substr(1, x_length / 2) %>% trimws(),
                     '.* ',
                     x.bak[i] %>% substr((x_length / 2) + 1, x_length) %>% trimws())
      found <- AMR::microorganisms %>% filter(fullname %like% paste0('^', x_split[i]))
    }
    if (nrow(found) == 0) {
      # try any match with text before and after original search string
      # so "negative rods" will be "GNR"
      if (x.bak[i] %like% "^Gram") {
        x.bak[i] <- gsub("^Gram", "", x.bak[i], ignore.case = TRUE)
        # remove leading and trailing spaces again
        x.bak[i] <- trimws(x.bak[i], which = "both")
      }
      found <- AMR::microorganisms %>% filter(fullname %like% x.bak[i])
    }

    if (nrow(found) != 0) {
      x[i] <- found %>%
        slice(1) %>%
        pull(bactid)
    } else {
      x[i] <- ""
    }
  }
  x
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
      found <- AMR::antibiotics %>% filter(molis == x[i] | umcg == x[i])
    }

    if (nrow(found) == 0) {
      # try exact official name
      found <- AMR::antibiotics[which(AMR::antibiotics$official == x[i]),]
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
