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

#' WHOCC: WHO Collaborating Centre for Drug Statistics Methodology
#'
#' All antimicrobial drugs and their official names, ATC codes, ATC groups and defined daily dose (DDD) are included in this package, using the WHO Collaborating Centre for Drug Statistics Methodology.
#' @section WHOCC:
#' \if{html}{\figure{logo_who.png}{options: height=60px style=margin-bottom:5px} \cr}
#' This package contains \strong{all ~450 antimicrobial drugs} and their Anatomical Therapeutic Chemical (ATC) codes, ATC groups and Defined Daily Dose (DDD) from the World Health Organization Collaborating Centre for Drug Statistics Methodology (WHOCC, \url{https://www.whocc.no}) and the Pharmaceuticals Community Register of the European Commission (\url{http://ec.europa.eu/health/documents/community-register/html/atc.htm}). 
#'
#' These have become the gold standard for international drug utilisation monitoring and research.
#'
#' The WHOCC is located in Oslo at the Norwegian Institute of Public Health and funded by the Norwegian government. The European Commission is the executive of the European Union and promotes its general interest.
#' 
#' \strong{NOTE: The WHOCC copyright does not allow use for commercial purposes, unlike any other info from this package. See \url{https://www.whocc.no/copyright_disclaimer/}.}
#' @inheritSection AMR Read more on our website!
#' @name WHOCC
#' @rdname WHOCC
#' @examples
#' as.ab("meropenem")
#' ab_name("J01DH02")
#'
#' ab_tradenames("flucloxacillin")
NULL
