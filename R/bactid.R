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

#' Transform to bacteria ID
#'
#' Use this function to determine a valid ID based on a genus (and species). This input can be a full name (like \code{"Staphylococcus aureus"}), an abbreviated name (like \code{"S. aureus"}), or just a genus. You could also \code{\link{select}} a genus and species column, zie Examples.
#' @param x a character vector or a dataframe with one or two columns
#' @rdname as.bactid
#' @details Some exceptions have been built in to get more logical results, based on prevalence of human pathogens. For example:
#' \itemize{
#'   \item{\code{"E. coli"} will return the ID of \emph{Escherichia coli} and not \emph{Entamoeba coli}, although the latter would alphabetically come first}
#'   \item{\code{"H. influenzae"} will return the ID of \emph{Haemophilus influenzae} and not \emph{Haematobacter influenzae}}
#'   \item{Something like \code{"p aer"} will return the ID of \emph{Pseudomonas aeruginosa} and not \emph{Pasteurella aerogenes}}
#'   \item{Something like \code{"stau"} or \code{"staaur"} will return the ID of \emph{Staphylococcus aureus} and not \emph{Staphylococcus auricularis}}
#' }
#' Moreover, this function also supports ID's based on only Gram stain, when the species is not known. \cr
#' For example, \code{"Gram negative rods"} and \code{"GNR"} will both return the ID of a Gram negative rod: \code{GNR}.
#' @export
#' @importFrom dplyr %>% filter pull
#' @return Character (vector) with class \code{"bactid"}. Unknown values will return \code{NA}.
#' @seealso \code{\link{microorganisms}} for the dataframe that is being used to determine ID's.
#' @examples
#' # These examples all return "STAAUR", the ID of S. aureus:
#' as.bactid("stau")
#' as.bactid("STAU")
#' as.bactid("staaur")
#' as.bactid("S. aureus")
#' as.bactid("S aureus")
#' as.bactid("Staphylococcus aureus")
#' as.bactid("MRSA") # Methicillin Resistant S. aureus
#' as.bactid("VISA") # Vancomycin Intermediate S. aureus
#' as.bactid("VRSA") # Vancomycin Resistant S. aureus
#'
#' \dontrun{
#' df$bactid <- as.bactid(df$microorganism_name)
#'
#' # the select function of tidyverse is also supported:
#' library(dplyr)
#' df$bactid <- df %>%
#'   select(microorganism_name) %>%
#'   as.bactid()
#'
#' # and can even contain 2 columns, which is convenient for genus/species combinations:
#' df$bactid <- df %>%
#'   select(genus, species) %>%
#'   as.bactid()
#'
#' # same result:
#' df <- df %>%
#'   mutate(bactid = paste(genus, species) %>%
#'                     as.bactid())
#' }
as.bactid <- function(x) {

  failures <- character(0)

  if (NCOL(x) == 2) {
    # support tidyverse selection like: df %>% select(colA, colB)
    # paste these columns together
    x_vector <- vector("character", NROW(x))
    for (i in 1:NROW(x)) {
      x_vector[i] <- paste(pull(x[i,], 1), pull(x[i,], 2), sep = " ")
    }
    x <- x_vector
  } else {
    if (NCOL(x) > 2) {
      stop('`x` can be 2 columns at most', call. = FALSE)
    }

    # support tidyverse selection like: df %>% select(colA)
    if (!is.vector(x)) {
      x <- pull(x, 1)
    }
  }

  x.fullbackup <- x
  # remove dots and other non-text in case of "E. coli" except spaces
  x <- gsub("[^a-zA-Z0-9 ]+", "", x)
  # but spaces before and after should be omitted
  x <- trimws(x, which = "both")
  x.backup <- x
  # replace space by regex sign
  x <- gsub(" ", ".*", x, fixed = TRUE)
  # add start and stop
  x_species <- paste(x, 'species')
  x <- paste0('^', x, '$')

  for (i in 1:length(x)) {

    if (x.fullbackup[i] %in% AMR::microorganisms$bactid) {
      # is already a valid bactid
      x[i] <- x.fullbackup[i]
      next
    }
    if (x.backup[i] %in% AMR::microorganisms$bactid) {
      # is already a valid bactid
      x[i] <- x.backup[i]
      next
    }

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
    if (tolower(x[i]) %like% 'coagulase'
        | tolower(x[i]) %like% 'cns'
        | tolower(x[i]) %like% 'cons') {
      # coerce S. coagulase negative, also as CNS and CoNS
      x[i] <- 'Coagulase Negative Staphylococcus (CNS)'
    }

    # translate known trivial names to genus+species
    if (!is.na(x.backup[i])) {
      if (toupper(x.backup[i]) == 'MRSA'
          | toupper(x.backup[i]) == 'VISA'
          | toupper(x.backup[i]) == 'VRSA') {
        x[i] <- 'Staphylococcus aureus'
      }
      if (toupper(x.backup[i]) == 'MRSE') {
        x[i] <- 'Staphylococcus epidermidis'
      }
      if (toupper(x.backup[i]) == 'VRE') {
        x[i] <- 'Enterococcus'
      }
      if (toupper(x.backup[i]) == 'MRPA') {
        # multi resistant P. aeruginosa
        x[i] <- 'Pseudomonas aeruginosa'
      }
      if (toupper(x.backup[i]) == 'PISP'
          | toupper(x.backup[i]) == 'PRSP') {
        # peni resistant S. pneumoniae
        x[i] <- 'Streptococcus pneumoniae'
      }
      if (toupper(x.backup[i]) == 'VISP'
          | toupper(x.backup[i]) == 'VRSP') {
        # vanco resistant S. pneumoniae
        x[i] <- 'Streptococcus pneumoniae'
      }
    }

    # let's try the ID's first
    found <- AMR::microorganisms %>% filter(bactid == x.backup[i])

    if (nrow(found) == 0) {
      # now try exact match
      found <- AMR::microorganisms %>% filter(fullname == x[i])
    }
    if (nrow(found) == 0) {
      # try any match
      found <- AMR::microorganisms %>% filter(fullname %like% x[i])
    }
    if (nrow(found) == 0) {
      # try exact match of only genus, with 'species' attached
      # (e.g. this prevents Streptococcus for becoming Peptostreptococcus, since "p" < "s")
      found <- AMR::microorganisms %>% filter(fullname == x_species[i])
    }
    if (nrow(found) == 0) {
      # try any match of only genus, with 'species' attached
      found <- AMR::microorganisms %>% filter(fullname %like% x_species[i])
    }
    if (nrow(found) == 0) {
      # search for GLIMS code
      if (toupper(x.backup[i]) %in% toupper(AMR::microorganisms.umcg$mocode)) {
        found <- AMR::microorganisms.umcg %>% filter(toupper(mocode) == toupper(x.backup[i]))
      }
    }
    if (nrow(found) == 0) {
      # try splitting of characters and then find ID
      # like esco = E. coli, klpn = K. pneumoniae, stau = S. aureus
      x_split <- x
      x_length <- nchar(x.backup[i])
      x_split[i] <- paste0(x.backup[i] %>% substr(1, x_length / 2) %>% trimws(),
                           '.* ',
                           x.backup[i] %>% substr((x_length / 2) + 1, x_length) %>% trimws())
      found <- AMR::microorganisms %>% filter(fullname %like% paste0('^', x_split[i]))
    }
    if (nrow(found) == 0) {
      # try any match with text before and after original search string
      # so "negative rods" will be "GNR"
      if (x.backup[i] %like% "^Gram") {
        x.backup[i] <- gsub("^Gram", "", x.backup[i], ignore.case = TRUE)
        # remove leading and trailing spaces again
        x.backup[i] <- trimws(x.backup[i], which = "both")
      }
      if (!is.na(x.backup[i])) {
        found <- AMR::microorganisms %>% filter(fullname %like% x.backup[i])
      }
    }

    if (nrow(found) != 0 & x.backup[i] != "") {
      x[i] <- as.character(found[1, 'bactid'])
    } else {
      x[i] <- NA_character_
      failures <- c(failures, x.fullbackup[i])
    }
  }

  failures <- failures[!failures %in% c(NA, NULL, NaN)]
  if (length(failures) > 0) {
    warning("These values could not be coerced to a valid bactid: ",
            paste('"', unique(failures), '"', sep = "", collapse = ', '),
            ".",
            call. = FALSE)
  }
  class(x) <- "bactid"
  attr(x, 'package') <- 'AMR'
  x
}

#' @rdname as.bactid
#' @export
guess_bactid <- as.bactid

#' @rdname as.bactid
#' @export
is.bactid <- function(x) {
  identical(class(x), "bactid")
}

#' @exportMethod print.bactid
#' @export
#' @noRd
print.bactid <- function(x, ...) {
  cat("Class 'bactid'\n")
  print.default(as.character(x), quote = FALSE)
}

#' @exportMethod as.data.frame.bactid
#' @export
#' @noRd
as.data.frame.bactid <- function (x, ...) {
  # same as as.data.frame.character but with removed stringsAsFactors
  nm <- paste(deparse(substitute(x), width.cutoff = 500L),
              collapse = " ")
  if (!"nm" %in% names(list(...))) {
    as.data.frame.vector(x, ..., nm = nm)
  } else {
    as.data.frame.vector(x, ...)
  }
}

#' @exportMethod pull.bactid
#' @export
#' @importFrom dplyr pull
#' @noRd
pull.bactid <- function(.data, ...) {
  pull(as.data.frame(.data), ...)
}
