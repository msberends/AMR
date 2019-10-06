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

#' Transform to antibiotic ID
#'
#' Use this function to determine the antibiotic code of one or more antibiotics. The data set \code{\link{antibiotics}} will be searched for abbreviations, official names and synonyms (brand names).
#' @param x character vector to determine to antibiotic ID
#' @param ... arguments passed on to internal functions
#' @rdname as.ab
#' @keywords atc
#' @inheritSection WHOCC WHOCC
#' @export
#' @importFrom dplyr %>% filter slice pull
#' @details All entries in the \code{\link{antibiotics}} data set have three different identifiers: a human readable EARS-Net code (column \code{ab}, used by ECDC and WHONET), an ATC code (column \code{atc}, used by WHO), and a CID code (column \code{cid}, Compound ID, used by PubChem). The data set contains more than 5,000 official brand names from many different countries, as found in PubChem.
#'
#' Use the \code{\link{ab_property}} functions to get properties based on the returned antibiotic ID, see Examples.
#' @section Source:
#' World Health Organization (WHO) Collaborating Centre for Drug Statistics Methodology: \url{https://www.whocc.no/atc_ddd_index/}
#'
#' WHONET 2019 software: \url{http://www.whonet.org/software.html}
#'
#' European Commission Public Health PHARMACEUTICALS - COMMUNITY REGISTER: \url{http://ec.europa.eu/health/documents/community-register/html/atc.htm}
#' @return Character (vector) with class \code{"ab"}. Unknown values will return \code{NA}.
#' @seealso \code{\link{antibiotics}} for the dataframe that is being used to determine ATCs.
#' @inheritSection AMR Read more on our website!
#' @examples
#' # These examples all return "ERY", the ID of Erythromycin:
#' as.ab("J01FA01")
#' as.ab("J 01 FA 01")
#' as.ab("Erythromycin")
#' as.ab("eryt")
#' as.ab("   eryt 123")
#' as.ab("ERYT")
#' as.ab("ERY")
#' as.ab("eritromicine") # spelled wrong, yet works
#' as.ab("Erythrocin")   # trade name
#' as.ab("Romycin")      # trade name
#'
#' # Use ab_* functions to get a specific properties (see ?ab_property);
#' # they use as.ab() internally:
#' ab_name("J01FA01")    # "Erythromycin"
#' ab_name("eryt")       # "Erythromycin"
as.ab <- function(x, ...) {
  if (is.ab(x)) {
    return(x)
  }

  if (all(toupper(x) %in% AMR::antibiotics$ab)) {
    # valid AB code, but not yet right class
    return(structure(.Data = toupper(x),
                     class = "ab"))
  }

  x_bak <- x
  # remove diacritics
  x <- iconv(x, from = "UTF-8", to = "ASCII//TRANSLIT")
  x <- gsub('"', "", x, fixed = TRUE)
  # remove suffices
  x_bak_clean <- gsub("_(mic|rsi|dis[ck])$", "", x, ignore.case = TRUE)
  # remove disk concentrations, like LVX_NM -> LVX
  x_bak_clean <- gsub("_[A-Z]{2}[0-9_]{0,3}$", "", x_bak_clean, ignore.case = TRUE)
  # remove part between brackets if that's followed by another string
  x_bak_clean <- gsub("(.*)+ [(].*[)]", "\\1", x_bak_clean)
  # keep only a-Z, 0-9, space, slash and dash
  # x_bak_clean <- gsub("[^A-Z0-9 /-]", "", x_bak_clean, ignore.case = TRUE)
  # keep only max 1 space
  x_bak_clean <- trimws(gsub(" +", " ", x_bak_clean, ignore.case = TRUE))
  # non-character, space or number should be a slash
  x_bak_clean <- gsub("[^A-Za-z0-9 -]", "/", x_bak_clean)
  # spaces around non-characters must be removed: amox + clav -> amox/clav
  x_bak_clean <- gsub("(.*[a-zA-Z0-9]) ([^a-zA-Z0-9].*)", "\\1\\2", x_bak_clean)
  x_bak_clean <- gsub("(.*[^a-zA-Z0-9]) ([a-zA-Z0-9].*)", "\\1\\2", x_bak_clean)

  x <- unique(x_bak_clean)
  x_new <- rep(NA_character_, length(x))
  x_unknown <- character(0)

  for (i in 1:length(x)) {
    if (is.na(x[i]) | is.null(x[i])) {
      next
    }
    if (identical(x[i], "")) {
      x_unknown <- c(x_unknown, x_bak[x[i] == x_bak_clean][1])
      next
    }
    # prevent "bacteria" from coercing to TMP, since Bacterial is a brand name of it
    if (identical(tolower(x[i]), "bacteria")) {
      x_unknown <- c(x_unknown, x_bak[x[i] == x_bak_clean][1])
      next
    }

    # exact AB code
    found <- AMR::antibiotics[which(AMR::antibiotics$ab == toupper(x[i])),]$ab
    if (length(found) > 0) {
      x_new[i] <- found[1L]
      next
    }

    # exact ATC code
    found <- AMR::antibiotics[which(AMR::antibiotics$atc == toupper(x[i])),]$ab
    if (length(found) > 0) {
      x_new[i] <- found[1L]
      next
    }

    # exact CID code
    found <- AMR::antibiotics[which(AMR::antibiotics$cid == x[i]),]$ab
    if (length(found) > 0) {
      x_new[i] <- found[1L]
      next
    }

    # exact name
    found <- AMR::antibiotics[which(toupper(AMR::antibiotics$name) == toupper(x[i])),]$ab
    if (length(found) > 0) {
      x_new[i] <- found[1L]
      next
    }

    # exact synonym
    synonym_found <- unlist(lapply(AMR::antibiotics$synonyms,
                                   function(s) if (toupper(x[i]) %in% toupper(s)) {
                                     TRUE
                                   } else {
                                     FALSE
                                   }))
    found <- AMR::antibiotics$ab[synonym_found == TRUE]
    if (length(found) > 0) {
      x_new[i] <- found[1L]
      next
    }

    # exact abbreviation
    abbr_found <- unlist(lapply(AMR::antibiotics$abbreviations,
                                function(a) if (toupper(x[i]) %in% toupper(a)) {
                                  TRUE
                                } else {
                                  FALSE
                                }))
    found <- AMR::antibiotics$ab[abbr_found == TRUE]
    if (length(found) > 0) {
      x_new[i] <- found[1L]
      next
    }

    # first >=4 characters of name
    if (nchar(x[i]) >= 4) {
      found <- AMR::antibiotics[which(toupper(AMR::antibiotics$name) %like% paste0("^", x[i])),]$ab
      if (length(found) > 0) {
        x_new[i] <- found[1L]
        next
      }
    }

    # allow characters that resemble others, but only continue when having more than 3 characters
    if (nchar(x[i]) <= 3) {
      x_unknown <- c(x_unknown, x_bak[x[i] == x_bak_clean][1])
      next
    }
    x_spelling <- tolower(x[i])
    x_spelling <- gsub("[iy]+", "[iy]+", x_spelling)
    x_spelling <- gsub("(c|k|q|qu|s|z|x|ks)+", "(c|k|q|qu|s|z|x|ks)+", x_spelling)
    x_spelling <- gsub("(ph|f|v)+", "(ph|f|v)+", x_spelling)
    x_spelling <- gsub("(th|t)+", "(th|t)+", x_spelling)
    x_spelling <- gsub("a+", "a+", x_spelling)
    x_spelling <- gsub("e+", "e+", x_spelling)
    x_spelling <- gsub("o+", "o+", x_spelling)
    # allow any ending of -in/-ine and -im/-ime
    x_spelling <- gsub("(\\[iy\\]\\+(n|m)|\\[iy\\]\\+(n|m)e\\+)$", "[iy]+(n|m)e*", x_spelling)
    # allow any ending of -ol/-ole
    x_spelling <- gsub("(o\\+l|o\\+le\\+)$", "o+le*", x_spelling)
    # allow any ending of -on/-one
    x_spelling <- gsub("(o\\+n|o\\+ne\\+)$", "o+ne*", x_spelling)
    # replace multiple same characters to single one with '+', like "ll" -> "l+"
    x_spelling <- gsub("(.)\\1+", "\\1+", x_spelling)
  
    # try if name starts with it
    found <- AMR::antibiotics[which(AMR::antibiotics$name %like% paste0("^", x_spelling)),]$ab
    if (length(found) > 0) {
      x_new[i] <- found[1L]
      next
    }
    # and try if any synonym starts with it
    synonym_found <- unlist(lapply(AMR::antibiotics$synonyms,
                                   function(s) if (any(s %like% paste0("^", x_spelling))) {
                                     TRUE
                                   } else {
                                     FALSE
                                   }))
    found <- AMR::antibiotics$ab[synonym_found == TRUE]
    if (length(found) > 0) {
      x_new[i] <- found[1L]
      next
    }

    # try by removing all spaces
    if (x[i] %like% " ") {
      found <- suppressWarnings(as.ab(gsub(" +", "", x[i])))
      if (length(found) > 0 & !is.na(found)) {
        x_new[i] <- found[1L]
        next
      }
    }

    # try by removing all spaces and numbers
    if (x[i] %like% " " | x[i] %like% "[0-9]") {
      found <- suppressWarnings(as.ab(gsub("[ 0-9]", "", x[i])))
      if (length(found) > 0 & !is.na(found)) {
        x_new[i] <- found[1L]
        next
      }
    }
    
    if (!isFALSE(list(...)$initial_search)) {
      # transform back from other languages and try again
      x_translated <- paste(lapply(strsplit(x[i], "[^a-zA-Z0-9 ]"),
                                   function(y) {
                                     for (i in 1:length(y)) {
                                       y[i] <- ifelse(tolower(y[i]) %in% tolower(translations_file$replacement),
                                                      translations_file[which(tolower(translations_file$replacement) == tolower(y[i]) &
                                                                                !isFALSE(translations_file$fixed)), "pattern"],
                                                      y[i])
                                     }
                                     y
                                   })[[1]],
                            collapse = "/")
      x_translated_guess <- suppressWarnings(as.ab(x_translated, initial_search = FALSE))
      if (!is.na(x_translated_guess)) {
        x_new[i] <- x_translated_guess
        next
      }
      
      if (!isFALSE(list(...)$initial_search2)) {
        # now also try to coerce brandname combinations like "Amoxy/clavulanic acid"
        x_translated <- paste(lapply(strsplit(x_translated, "[^a-zA-Z0-9 ]"),
                                     function(y) {
                                       for (i in 1:length(y)) {
                                         y_name <- suppressWarnings(ab_name(y[i], language = NULL, initial_search = FALSE, initial_search2 = FALSE))
                                         y[i] <- ifelse(!is.na(y_name),
                                                        y_name,
                                                        y[i])
                                       }
                                       y
                                     })[[1]],
                              collapse = "/")
        x_translated_guess <- suppressWarnings(as.ab(x_translated, initial_search = FALSE))
        if (!is.na(x_translated_guess)) {
          x_new[i] <- x_translated_guess
          next
        }
      }
    }
    
    # not found
    x_unknown <- c(x_unknown, x_bak[x[i] == x_bak_clean][1])
  }

  # take failed ATC codes apart from rest
  x_unknown_ATCs <- x_unknown[x_unknown %like% "[A-Z][0-9][0-9][A-Z][A-Z][0-9][0-9]"]
  x_unknown <- x_unknown[!x_unknown %in% x_unknown_ATCs]
  if (length(x_unknown_ATCs) > 0) {
    warning("These ATC codes are not (yet) in the antibiotics data set: ",
            paste('"', sort(unique(x_unknown_ATCs)), '"', sep = "", collapse = ', '),
            ".",
            call. = FALSE)
  }

  if (length(x_unknown) > 0) {
    warning("These values could not be coerced to a valid antimicrobial ID: ",
            paste('"', sort(unique(x_unknown)), '"', sep = "", collapse = ', '),
            ".",
            call. = FALSE)
  }

  x_result <- data.frame(x = x_bak_clean, stringsAsFactors = FALSE) %>%
    left_join(data.frame(x = x, x_new = x_new, stringsAsFactors = FALSE), by = "x") %>%
    pull(x_new)

  if (length(x_result) == 0) {
    x_result <- NA_character_
  }

  structure(.Data = x_result,
            class = "ab")
}

#' @rdname as.ab
#' @export
is.ab <- function(x) {
  identical(class(x), "ab")
}

#' @exportMethod print.ab
#' @export
#' @noRd
print.ab <- function(x, ...) {
  cat("Class 'ab'\n")
  print(as.character(x), quote = FALSE)
}

#' @exportMethod as.data.frame.ab
#' @export
#' @noRd
as.data.frame.ab <- function (x, ...) {
  # same as as.data.frame.character but with removed stringsAsFactors
  nm <- paste(deparse(substitute(x), width.cutoff = 500L),
              collapse = " ")
  if (!"nm" %in% names(list(...))) {
    as.data.frame.vector(x, ..., nm = nm)
  } else {
    as.data.frame.vector(x, ...)
  }
}

#' @exportMethod [.ab
#' @export
#' @noRd
"[.ab" <- function(x, ...) {
  y <- NextMethod()
  attributes(y) <- attributes(x)
  y
}
#' @exportMethod [[.ab
#' @export
#' @noRd
"[[.ab" <- function(x, ...) {
  y <- NextMethod()
  attributes(y) <- attributes(x)
  y
}
#' @exportMethod [<-.ab
#' @export
#' @noRd
"[<-.ab" <- function(i, j, ..., value) {
  y <- NextMethod()
  attributes(y) <- attributes(i)
  class_integrity_check(y, "antimicrobial code", AMR::antibiotics$ab)
}
#' @exportMethod [[<-.ab
#' @export
#' @noRd
"[[<-.ab" <- function(i, j, ..., value) {
  y <- NextMethod()
  attributes(y) <- attributes(i)
  class_integrity_check(y, "antimicrobial code", AMR::antibiotics$ab)
}
#' @exportMethod c.ab
#' @export
#' @noRd
c.ab <- function(x, ...) {
  y <- NextMethod()
  attributes(y) <- attributes(x)
  class_integrity_check(y, "antimicrobial code", AMR::antibiotics$ab)
}

#' @importFrom pillar type_sum
#' @export
type_sum.ab <- function(x) {
  "ab"
}

#' @importFrom pillar pillar_shaft
#' @export
pillar_shaft.ab <- function(x, ...) {
  out <- format(x)
  out[is.na(x)] <- pillar::style_na("NA")
  pillar::new_pillar_shaft_simple(out, align = "left", min_width = 4)
}
