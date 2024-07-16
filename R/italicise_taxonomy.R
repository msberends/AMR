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
# how to conduct AMR data analysis: https://msberends.github.io/AMR/   #
# ==================================================================== #

#' Italicise Taxonomic Families, Genera, Species, Subspecies
#'
#' According to the binomial nomenclature, the lowest four taxonomic levels (family, genus, species, subspecies) should be printed in italics. This function finds taxonomic names within strings and makes them italic.
#' @param string a [character] (vector)
#' @param type type of conversion of the taxonomic names, either "markdown", "html" or "ansi", see *Details*
#' @details
#' This function finds the taxonomic names and makes them italic based on the [microorganisms] data set.
#'
#' The taxonomic names can be italicised using markdown (the default) by adding `*` before and after the taxonomic names, or `<i>` and `</i>` when using html. When using 'ansi', ANSI colours will be added using `\033[3m` before and `\033[23m` after the taxonomic names. If multiple ANSI colours are not available, no conversion will occur.
#'
#' This function also supports abbreviation of the genus if it is followed by a species, such as "E. coli" and "K. pneumoniae ozaenae".
#' @export
#' @examples
#' italicise_taxonomy("An overview of Staphylococcus aureus isolates")
#' italicise_taxonomy("An overview of S. aureus isolates")
#'
#' cat(italicise_taxonomy("An overview of S. aureus isolates", type = "ansi"))
italicise_taxonomy <- function(string, type = c("markdown", "ansi", "html")) {
  if (missing(type)) {
    type <- "markdown"
  }
  meet_criteria(string, allow_class = "character")
  meet_criteria(type, allow_class = "character", has_length = 1, is_in = c("markdown", "ansi", "html"))

  add_MO_lookup_to_AMR_env()

  if (type == "markdown") {
    before <- "*"
    after <- "*"
  } else if (type == "html") {
    before <- "<i>"
    after <- "</i>"
  } else if (type == "ansi") {
    if (!has_colour() && !identical(Sys.getenv("IN_PKGDOWN"), "true")) {
      return(string)
    }
    before <- "\033[3m"
    after <- "\033[23m"
  }

  vapply(
    FUN.VALUE = character(1),
    string,
    function(s) {
      s_split <- unlist(strsplit(s, " ", fixed = TRUE))

      search_strings <- gsub("[^a-zA-Z-]", "", s_split)

      ind_species <- search_strings != "" &
        search_strings %in% AMR_env$MO_lookup[
          which(AMR_env$MO_lookup$rank %in% c(
            "family",
            "genus",
            "species",
            "subspecies",
            "infraspecies",
            "subsp."
          )),
          "species",
          drop = TRUE
        ]

      ind_fullname <- search_strings != "" &
        search_strings %in% c(
          AMR_env$MO_lookup[
            which(AMR_env$MO_lookup$rank %in% c(
              "family",
              "genus",
              "species",
              "subspecies",
              "infraspecies",
              "subsp."
            )),
            "fullname",
            drop = TRUE
          ],
          AMR_env$MO_lookup[
            which(AMR_env$MO_lookup$rank %in% c(
              "family",
              "genus",
              "species",
              "subspecies",
              "infraspecies",
              "subsp."
            )),
            "subspecies",
            drop = TRUE
          ]
        )

      # also support E. coli, add "E." to indices
      has_previous_genera_abbr <- s_split[which(ind_species) - 1] %like_case% "^[A-Z][.]?$"
      ind_species <- c(which(ind_species), which(ind_species)[has_previous_genera_abbr] - 1)

      ind <- c(ind_species, which(ind_fullname))

      s_split[ind] <- paste0(before, s_split[ind], after)
      s_paste <- paste(s_split, collapse = " ")

      # clean up a bit
      s_paste <- gsub(paste0(after, " ", before), " ", s_paste, fixed = TRUE)

      s_paste
    },
    USE.NAMES = FALSE
  )
}

#' @rdname italicise_taxonomy
#' @export
italicize_taxonomy <- function(string, type = c("markdown", "ansi", "html")) {
  if (missing(type)) {
    type <- "markdown"
  }
  italicise_taxonomy(string = string, type = type)
}
