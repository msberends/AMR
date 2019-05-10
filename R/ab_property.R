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

#' Property of an antibiotic
#'
#' Use these functions to return a specific property of an antibiotic from the \code{\link{antibiotics}} data set. All input values will be evaluated internally with \code{\link{as.ab}}.
#' @param x any (vector of) text that can be coerced to a valid microorganism code with \code{\link{as.ab}}
#' @param tolower logical to indicate whether the first character of every output should be transformed to a lower case character. This will lead to e.g. "polymyxin B" and not "polymyxin b".
#' @param property one of the column names of one of the \code{\link{antibiotics}} data set
#' @param language language of the returned text, defaults to system language (see \code{\link{get_locale}}) and can also be set with \code{\link{getOption}("AMR_locale")}. Use \code{language = NULL} or \code{language = ""} to prevent translation.
#' @param administration way of administration, either \code{"oral"} or \code{"iv"}
#' @param units a logical to indicate whether the units instead of the DDDs itself must be returned, see Examples
#' @param ... other parameters passed on to \code{\link{as.ab}}
#' @details All output will be \link{translate}d where possible.
#' @inheritSection as.ab Source
#' @rdname ab_property
#' @name ab_property
#' @return \itemize{
#'   \item{An \code{integer} in case of \code{ab_cid}}
#'   \item{A named \code{list} in case of multiple \code{ab_synonyms}}
#'   \item{A \code{double} in case of \code{ab_ddd}}
#'   \item{A \code{character} in all other cases}
#' }
#' @export
#' @seealso \code{\link{antibiotics}}
#' @inheritSection AMR Read more on our website!
#' @examples
#' # all properties:
#' ab_name("AMX")       # "Amoxicillin"
#' ab_atc("AMX")        # J01CA04 (ATC code from the WHO)
#' ab_cid("AMX")        # 33613 (Compound ID from PubChem)
#'
#' ab_synonyms("AMX")   # a list with brand names of amoxicillin
#' ab_tradenames("AMX") # same
#'
#' ab_group("AMX")      # "Beta-lactams/penicillins"
#' ab_atc_group1("AMX") # "Beta-lactam antibacterials, penicillins"
#' ab_atc_group2("AMX") # "Penicillins with extended spectrum"
#'
#' ab_name(x = c("AMC", "PLB"))  # "Amoxicillin/clavulanic acid" "Polymyxin B"
#' ab_name(x = c("AMC", "PLB"),
#'         tolower = TRUE)       # "amoxicillin/clavulanic acid" "polymyxin B"
#'
#' ab_ddd("AMX", "oral")               #  1
#' ab_ddd("AMX", "oral", units = TRUE) # "g"
#' ab_ddd("AMX", "iv")                 #  1
#' ab_ddd("AMX", "iv", units = TRUE)   # "g"
#'
#' # all ab_* functions use as.ab() internally:
#' ab_name("Fluclox")   # "Flucloxacillin"
#' ab_name("fluklox")   # "Flucloxacillin"
#' ab_name("floxapen")  # "Flucloxacillin"
#' ab_name(21319)       # "Flucloxacillin" (using CID)
#' ab_name("J01CF05")   # "Flucloxacillin" (using ATC)
ab_name <- function(x, language = get_locale(), tolower = FALSE, ...) {
  x <- ab_validate(x = x, property = "name", ...)
  res <- t(x, language = language)
  if (tolower == TRUE) {
    # use perl to only transform the first character
    # as we want "polymyxin B", not "polymyxin b"
    res <- gsub("^([A-Z])", "\\L\\1", res, perl = TRUE)
  }
  res
}

#' @rdname ab_property
#' @export
ab_atc <- function(x, ...) {
  ab_validate(x = x, property = "atc", ...)
}

#' @rdname ab_property
#' @export
ab_cid <- function(x, ...) {
  ab_validate(x = x, property = "cid", ...)
}

#' @rdname ab_property
#' @export
ab_synonyms <- function(x, ...) {
  syns <- ab_validate(x = x, property = "synonyms", ...)
  names(syns) <- x
  if (length(syns) == 1) {
    unname(unlist(syns))
  } else {
    syns
  }
}

#' @rdname ab_property
#' @export
ab_tradenames <- function(x, ...) {
  ab_synonyms(x, ...)
}

#' @rdname ab_property
#' @export
ab_group <- function(x, ...) {
  ab_validate(x = x, property = "group", ...)
}

#' @rdname ab_property
#' @export
ab_atc_group1 <- function(x, ...) {
  ab_validate(x = x, property = "atc_group1", ...)
}

#' @rdname ab_property
#' @export
ab_atc_group2 <- function(x, ...) {
  ab_validate(x = x, property = "atc_group2", ...)
}

#' @rdname ab_property
#' @export
ab_ddd <- function(x, administration = "oral", units = FALSE, ...) {
  if (!administration %in% c("oral", "iv")) {
    stop("`administration` must be 'oral' or 'iv'", call. = FALSE)
  }
  ddd_prop <- administration
  if (units == TRUE) {
    ddd_prop <- paste0(ddd_prop, "_units")
  } else {
    ddd_prop <- paste0(ddd_prop, "_ddd")
  }
  ab_validate(x = x, property = ddd_prop, ...)
}

#' @rdname ab_property
#' @export
ab_property <- function(x, property = 'name', language = get_locale(), ...) {
  if (length(property) != 1L) {
    stop("'property' must be of length 1.")
  }
  if (!property %in% colnames(AMR::antibiotics)) {
    stop("invalid property: '", property, "' - use a column name of the `antibiotics` data set")
  }

  t(ab_validate(x = x, property = property, ...), language = language)
}

ab_validate <- function(x, property, ...) {
  if (!"AMR" %in% base::.packages()) {
    library("AMR")
    # check onLoad() in R/zzz.R: data tables are created there.
  }

  # try to catch an error when inputting an invalid parameter
  # so the 'call.' can be set to FALSE
  tryCatch(x[1L] %in% AMR::antibiotics[1, property],
           error = function(e) stop(e$message, call. = FALSE))

  if (!all(x %in% AMR::antibiotics[, property])) {
    x <- data.frame(ab = as.ab(x), stringsAsFactors = FALSE) %>%
      left_join(antibiotics %>% select(c("ab", property)), by = "ab") %>%
      pull(property)
  }
  if (property %in% c("ab", "atc")) {
    return(structure(x, class = property))
  } else if (property == "cid") {
    return(as.integer(x))
  } else if (property %like% "ddd") {
    return(as.double(x))
  } else {
    return(x)
  }
}
