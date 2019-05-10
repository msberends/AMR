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

#' Translate strings from AMR package
#'
#' For language-dependent output of AMR functions, like \code{\link{mo_fullname}} and \code{\link{mo_type}}.
#' @details Strings will be translated to foreign languages if they are defined in a local translation file. This file comes with this package and can be found when running:
#'
#' \code{system.file("translations.tsv", package = "AMR")}
#'
#' This file will be read by all functions where a translated output can be desired, like all \code{\link{mo_property}} functions (\code{\link{mo_fullname}}, \code{\link{mo_type}}, etc.). Please suggest your own translations \href{https://gitlab.com/msberends/AMR/issues/new?issue[title]=Translation suggestion}{by creating a new issue on our repository}.
#'
#' The system language will be used at default, if supported, using \code{\link{get_locale}}. The system language can be overwritten with \code{\link{getOption}("AMR_locale")}.
#' @inheritSection AMR Read more on our website!
#' @rdname translate
#' @name translate
#' @export
#' @examples
#' # The 'language' parameter of below functions
#' # will be set automatically to your system language
#' # with get_locale()
#'
#' # English
#' mo_fullname("CoNS", language = "en")
#' #> "Coagulase-negative Staphylococcus (CoNS)"
#'
#' # German
#' mo_fullname("CoNS", language = "de")
#' #> "Koagulase-negative Staphylococcus (KNS)"
#'
#' # Dutch
#' mo_fullname("CoNS", language = "nl")
#' #> "Coagulase-negatieve Staphylococcus (CNS)"
#'
#' # Spanish
#' mo_fullname("CoNS", language = "es")
#' #> "Staphylococcus coagulasa negativo (SCN)"
#'
#' # Italian
#' mo_fullname("CoNS", language = "it")
#' #> "Staphylococcus negativo coagulasi (CoNS)"
#'
#' # Portuguese
#' mo_fullname("CoNS", language = "pt")
#' #> "Staphylococcus coagulase negativo (CoNS)"
get_locale <- function() {
  if (getOption("AMR_locale", "en") != "en") {
    return(getOption("AMR_locale"))
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
    # other language -> set to English
    "en"
  }
}
