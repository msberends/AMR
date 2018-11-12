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

#' Get language for AMR
#'
#' Determines the system language to be used for language-dependent output of AMR functions, like \code{\link{mo_gramstain}} and \code{\link{mo_type}}.
#' @details The system language can be overwritten with \code{\link{getOption}("AMR_locale")}.
#' @section Supported languages:
#' Supported languages are \code{"en"} (English), \code{"de"} (German), \code{"nl"} (Dutch), \code{"es"} (Spanish), \code{"it"} (Italian), \code{"fr"} (French), and \code{"pt"} (Portuguese).
#' @export
get_locale <- function() {
  if (!is.null(getOption("AMR_locale"))) {
    if (getOption("AMR_locale") %in% c("en", "de", "nl", "es", "it", "fr", "pt")) {
      return(getOption("AMR_locale"))
    }
  }
  lang <- Sys.getlocale("LC_COLLATE")
  # grepl with case = FALSE is faster than like
  if (grepl("^(English|en_|EN_)", lang, ignore.case = FALSE)) {
    # as first option to optimise speed
    "en"
  } else if (grepl("^(German|Deutsch|de_|DE_)", lang, ignore.case = FALSE)) {
    "de"
  } else if (grepl("^(Dutch|Nederlands|nl_|NL_)", lang, ignore.case = FALSE)) {
    "nl"
  } else if (grepl("^(Spanish|Espa.ol|es_|ES_)", lang, ignore.case = FALSE)) {
    "es"
  } else if (grepl("^(Italian|Italiano|it_|IT_)", lang, ignore.case = FALSE)) {
    "it"
  } else if (grepl("^(French|Fran.ais|fr_|FR_)", lang, ignore.case = FALSE)) {
    "fr"
  } else if (grepl("^(Portuguese|Portugu.s|pt_|PT_)", lang, ignore.case = FALSE)) {
    "pt"
  } else {
    # other language, set to English
    "en"
  }
}
