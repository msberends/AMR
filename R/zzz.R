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

#' The \code{AMR} Package
#'
#' Welcome to the \code{AMR} package. This page gives some additional contact information about the authors.
#' @details
#' This package was intended to simplify the analysis and prediction of Antimicrobial Resistance (AMR) and work with antibiotic properties by using evidence-based methods.
#'
#' This package was created for academic research by PhD students of the Faculty of Medical Sciences of the University of Groningen and the Medical Microbiology & Infection Prevention (MMBI) department of the University Medical Center Groningen (UMCG).
#' @section Authors:
#' Matthijs S. Berends[1,2] Christian F. Luz[1], Erwin E.A. Hassing[2],  Corinna Glasner[1],  Alex W. Friedrich[1],  Bhanu Sinha[1] \cr
#'
#' [1] Department of Medical Microbiology, University of Groningen, University Medical Center Groningen, Groningen, the Netherlands - \url{rug.nl} \url{umcg.nl} \cr
#' [2] Certe Medical Diagnostics & Advice, Groningen, the Netherlands - \url{certe.nl}
#' @section Contact us:
#' For suggestions, comments or questions, please contact us at:
#'
#' Matthijs S. Berends \cr
#' m.s.berends [at] umcg [dot] nl \cr
#' Department of Medical Microbiology, University of Groningen \cr
#' University Medical Center Groningen \cr
#' Post Office Box 30001 \cr
#' 9700 RB Groningen
#'
#' If you have found a bug, please file a new issue at: \cr
#' \url{https://gitlab.com/msberends/AMR/issues}
#' @name AMR
#' @rdname AMR
NULL

.onLoad <- function(libname, pkgname) {
  backports::import(pkgname)
}
