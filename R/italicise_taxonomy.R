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
# Center Groningen in the Netherlands, in collaboration with many      #
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

#' Italicise Taxonomic Families, Genera, Species, and Subspecies
#'
#' According to the binomial nomenclature, the lowest four taxonomic levels (family, genus, species, subspecies) should be printed in italics. This function finds taxonomic names within strings and makes them italic.
#' @param string A [character] (vector).
#' @param type Type of conversion of the taxonomic names, either "markdown", "html", "ansi" or "plotmath", see *Details*.
#' @details
#' This function finds the taxonomic names and makes them italic based on the [microorganisms] data set. Family and genus names are italicised anywhere in the text. Species and subspecies names are only italicised when they follow a genus, an abbreviated genus (such as "E." in "E. coli"), or another species or subspecies name, so that ordinary words that are also species names (such as "coli" or "major") stay upright. Subspecies may also follow "subsp.", "ssp." or "var.", which stay upright themselves, e.g., "*K. pneumoniae* subsp. *ozaenae*". *Salmonella* serovars (such as Typhimurium) and groups are not italicised, as is the convention. Punctuation around the names, such as brackets and commas, stays outside the italics.
#'
#' The taxonomic names can be italicised using markdown (the default) by adding `*` before and after the taxonomic names, or `<i>` and `</i>` when using html. When using 'ansi', ANSI colours will be added using `\033[3m` before and `\033[23m` after the taxonomic names. If multiple ANSI colours are not available, no conversion will occur.
#'
#' Use `type = "plotmath"` to get an [expression] vector for plots (see [plotmath][grDevices::plotmath]), which works for titles, axes and legends in both base \R and `ggplot2`. For axis and legend labels, the function [labels_taxonomy()] does exactly this and can directly be used as labeller, e.g., `ggplot2::scale_x_discrete(labels = labels_taxonomy)`.
#' @return For "markdown", "html" and "ansi", a [character] vector. For "plotmath", an [expression] vector. Both have the same length as `string` and keep `NA` values.
#' @seealso [labels_taxonomy()] to italicise taxonomic names in plot labels.
#' @export
#' @examples
#' italicise_taxonomy("An overview of Staphylococcus aureus isolates")
#' italicise_taxonomy("An overview of S. aureus isolates")
#' italicise_taxonomy("Infections with E. coli (and coli bacteria in general)")
#'
#' cat(italicise_taxonomy("An overview of S. aureus isolates", type = "ansi"))
#'
#' # for plots:
#' italicise_taxonomy("An overview of S. aureus isolates", type = "plotmath")
#' plot(1:10, main = italicise_taxonomy("Growth of S. aureus", type = "plotmath"))
italicise_taxonomy <- function(string, type = c("markdown", "ansi", "html", "plotmath")) {
  if (missing(type)) {
    type <- "markdown"
  }
  meet_criteria(string, allow_class = "character")
  meet_criteria(type, allow_class = "character", has_length = 1, is_in = c("markdown", "ansi", "html", "plotmath"))

  pieces <- taxonomy_pieces(string)

  if (type == "plotmath") {
    return(as.expression(lapply(pieces, function(p) {
      if (is.null(p)) {
        return(NA_character_)
      }
      if (!any(p$italic)) {
        return(paste0(p$text, collapse = ""))
      }
      parts <- lapply(seq_along(p$text), function(i) {
        if (p$italic[i]) call("italic", p$text[i]) else p$text[i]
      })
      if (length(parts) == 1) {
        parts[[1]]
      } else {
        as.call(c(as.name("paste"), parts))
      }
    })))
  }

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
    pieces,
    function(p) {
      if (is.null(p)) {
        return(NA_character_)
      }
      p$text[p$italic] <- paste0(before, p$text[p$italic], after)
      paste0(p$text, collapse = "")
    },
    USE.NAMES = FALSE
  )
}

contains_taxonomy <- function(x) {
  # whether any value holds a taxonomic name that would be italicised
  x <- unique(as.character(x))
  any(vapply(taxonomy_pieces(x[!is.na(x)]), function(p) any(p$italic), logical(1)))
}

get_taxonomy_names <- function() {
  # cached for the session, reset by reset_mo_cache() when the lookup table changes
  if (is.null(AMR_env$taxonomy_names)) {
    add_MO_lookup_to_AMR_env()
    lookup <- AMR_env$MO_lookup
    families <- unique(lookup$family)
    genera <- unique(lookup$genus)
    epithets <- unique(c(
      lookup$species[lookup$rank %in% c("species", "subspecies", "species group")],
      lookup$subspecies[lookup$rank == "subspecies"]
    ))
    AMR_env$taxonomy_names <- list(
      # only real names, so no "(unknown genus)" or similar
      families = families[families %like_case% "^[A-Z][a-z-]+$"],
      genera = genera[genera %like_case% "^[A-Z][a-z-]+$"],
      # only lowercase epithets, so no Salmonella serovars ("Typhimurium") or groups ("Group A")
      epithets = epithets[epithets %like_case% "^[a-z][a-z-]*$"]
    )
  }
  AMR_env$taxonomy_names
}

taxonomy_pieces <- function(string) {
  # returns a list with for every string NULL (if NA) or a list with `text` and `italic`
  nms <- get_taxonomy_names()
  out <- vector("list", length(string))
  if (length(string) == 0) {
    return(out)
  }
  # split all strings at once, so the lookups below run only once
  words <- strsplit(string, " ", fixed = TRUE)
  # strsplit() drops a trailing empty string
  trailing <- which(!is.na(string) & string %like_case% " $")
  words[trailing] <- lapply(words[trailing], c, "")
  n_words <- lengths(words)
  all_words <- unlist(words, use.names = FALSE)
  all_words[is.na(all_words)] <- ""

  # leading punctuation, the word itself, trailing punctuation
  prefix <- gsub("^([^A-Za-z]*).*$", "\\1", all_words, perl = TRUE)
  rest <- substr(all_words, nchar(prefix) + 1, nchar(all_words))
  core <- gsub("^([A-Za-z]([A-Za-z-]*[A-Za-z])?).*$", "\\1", rest, perl = TRUE)
  core[rest %unlike_case% "^[A-Za-z]"] <- ""
  suffix <- substr(rest, nchar(core) + 1, nchar(rest))
  # abbreviated genus, e.g. "E." in "E. coli"
  is_abbr <- core %like_case% "^[A-Z]$" & suffix %like_case% "^[.]"
  core[is_abbr] <- paste0(core[is_abbr], ".")
  suffix[is_abbr] <- substr(suffix[is_abbr], 2, nchar(suffix[is_abbr]))

  is_genus <- core != "" & core %in% nms$genera
  is_higher <- is_genus | (core != "" & core %in% nms$families)
  is_epithet <- core != "" & core %in% nms$epithets
  is_connector <- all_words %in% c("subsp.", "ssp.", "var.")
  ends <- cumsum(n_words)
  starts <- ends - n_words + 1

  for (s in which(!is.na(string))) {
    ind <- seq.int(starts[s], length.out = n_words[s])
    pre <- prefix[ind]
    cr <- core[ind]
    suf <- suffix[ind]
    higher <- is_higher[ind]
    genus <- is_genus[ind]
    abbr <- is_abbr[ind]
    epithet <- is_epithet[ind]
    connector <- is_connector[ind]
    n <- length(ind)
    italic <- higher
    # whether the next word can be a species or subspecies name
    can_follow <- logical(n)
    for (i in seq_len(n)) {
      if (i > 1 && epithet[i] && pre[i] == "" && can_follow[i - 1]) {
        italic[i] <- TRUE
        if (abbr[i - 1]) {
          italic[i - 1] <- TRUE
        }
      }
      if (suf[i] == "" || connector[i]) {
        can_follow[i] <- genus[i] ||
          abbr[i] ||
          (italic[i] && epithet[i]) ||
          (connector[i] && i > 1 && italic[i - 1] && epithet[i - 1])
      }
    }
    # connectors such as "subsp." are not part of the core, so keep them upright
    cr[connector] <- all_words[ind][connector]
    suf[connector] <- ""
    italic[connector] <- FALSE

    text <- character(0)
    ital <- logical(0)
    for (i in seq_len(n)) {
      text <- c(text, pre[i], cr[i], suf[i])
      ital <- c(ital, FALSE, italic[i], FALSE)
      if (i < n) {
        # the space between two italic words belongs to the italics
        text <- c(text, " ")
        ital <- c(ital, italic[i] && italic[i + 1] && suf[i] == "" && pre[i + 1] == "")
      }
    }
    keep <- text != ""
    text <- text[keep]
    ital <- ital[keep]
    if (length(text) == 0) {
      out[[s]] <- list(text = "", italic = FALSE)
      next
    }
    # merge adjacent pieces of the same style
    grp <- cumsum(c(TRUE, ital[-1] != ital[-length(ital)]))
    out[s] <- list(list(
      text = vapply(split(text, grp), paste0, character(1), collapse = "", USE.NAMES = FALSE),
      italic = vapply(split(ital, grp), function(x) x[1], logical(1), USE.NAMES = FALSE)
    ))
  }
  out
}

#' @rdname italicise_taxonomy
#' @export
italicize_taxonomy <- function(string, type = c("markdown", "ansi", "html", "plotmath")) {
  if (missing(type)) {
    type <- "markdown"
  }
  italicise_taxonomy(string = string, type = type)
}
