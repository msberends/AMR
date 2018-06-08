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

#' Find bacteria ID based on genus/species
#'
#' Use this function to determine a valid ID based on a genus (and species). This input could be a full name (like \code{"Staphylococcus aureus"}), an abbreviated name (like \code{"S. aureus"}), or just a genus. You could also \code{\link{select}} a genus and species column, zie Examples.
#' @param x character vector or a dataframe with one or two columns
#' @export
#' @importFrom dplyr %>% filter pull
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
#'
#' \dontrun{
#' df$bactid <- guess_bactid(df$microorganism_name)
#'
#' # the select function of tidyverse is also supported:
#' df$bactid <- df %>% select(microorganism_name) %>% guess_bactid()
#'
#' # and can even contain 2 columns, which is convenient for genus/species combinations:
#' df$bactid <- df %>% select(genus, species) %>% guess_bactid()
#' # same result:
#' df <- df %>% mutate(bactid = paste(genus, species) %>% guess_bactid())
#' }
guess_bactid <- function(x) {

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
    if (!is.na(x.bak[i])) {
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
      if (!is.na(x.bak[i])) {
        found <- AMR::microorganisms %>% filter(fullname %like% x.bak[i])
      }
    }

    if (nrow(found) != 0) {
      x[i] <- as.character(found[1, 'bactid'])
    } else {
      x[i] <- ""
    }
  }
  x
}
