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
# Visit our website for more info: https://msberends.gitab.io/AMR.     #
# ==================================================================== #

#' The \code{AMR} Package
#'
#' Welcome to the \code{AMR} package. This page gives some additional contact information about the authors.
#' @details
#' This package was intended to simplify the analysis and prediction of Antimicrobial Resistance (AMR) and to work with microbial and antimicrobial properties by using evidence-based methods.
#'
#' This package was created for both academic research and routine analysis by PhD students of the Faculty of Medical Sciences of the University of Groningen and the Medical Microbiology & Infection Prevention (MMBI) department of the University Medical Center Groningen (UMCG).
#' @section Read more on our website!:
#' \if{html}{\figure{logo.png}{options: height=40px style=margin-bottom:5px} \cr}
#' On our website \url{https://msberends.gitlab.io/AMR} you can find \href{https://msberends.gitlab.io/AMR/articles/AMR.html}{a omprehensive tutorial} about how to conduct AMR analysis and find \href{https://msberends.gitlab.io/AMR/reference}{the complete documentation of all functions}, which reads a lot easier than in R.
#' @section Authors:
#' Matthijs S. Berends[1,2] Christian F. Luz[1], Erwin E.A. Hassing[2],  Corinna Glasner[1],  Alex W. Friedrich[1],  Bhanu N.M. Sinha[1] \cr
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
