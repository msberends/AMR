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

#' The \code{AMR} Package
#'
#' Welcome to the \code{AMR} package.
#' @details
#' \code{AMR} is a free and open-source R package to simplify the analysis and prediction of Antimicrobial Resistance (AMR) and to work with microbial and antimicrobial properties by using evidence-based methods. It supports any table format, including WHONET/EARS-Net data.
#'
#' We created this package for both academic research and routine analysis at the Faculty of Medical Sciences of the University of Groningen and the Medical Microbiology & Infection Prevention (MMBI) department of the University Medical Center Groningen (UMCG). This R package is actively maintained and free software; you can freely use and distribute it for both personal and commercial (but not patent) purposes under the terms of the GNU General Public License version 2.0 (GPL-2), as published by the Free Software Foundation.
#'
#' This package can be used for:
#' \itemize{
#'   \item{Reference for microorganisms, since it contains all microbial (sub)species from the Catalogue of Life}
#'   \item{Interpreting raw MIC and disk diffusion values, based on the latest CLSI or EUCAST guidelines}
#'   \item{Calculating antimicrobial resistance}
#'   \item{Determining multi-drug resistance (MDR) / multi-drug resistant organisms (MDRO)}
#'   \item{Calculating empirical susceptibility of both mono therapy and combination therapy}
#'   \item{Predicting future antimicrobial resistance using regression models}
#'   \item{Getting properties for any microorganism (like Gram stain, species, genus or family)}
#'   \item{Getting properties for any antibiotic (like name, ATC code, defined daily dose or trade name)}
#'   \item{Plotting antimicrobial resistance}
#'   \item{Determining first isolates to be used for AMR analysis}
#'   \item{Applying EUCAST expert rules}
#'   \item{Descriptive statistics: frequency tables, kurtosis and skewness}
#' }
#' @section Authors:
#' Matthijs S. Berends[1,2] Christian F. Luz[1], Erwin E.A. Hassing[2],  Corinna Glasner[1],  Alex W. Friedrich[1],  Bhanu N.M. Sinha[1] \cr
#'
#' [1] Department of Medical Microbiology, University of Groningen, University Medical Center Groningen, Groningen, the Netherlands - \url{https://www.rug.nl} \url{https://www.umcg.nl} \cr
#' [2] Certe Medical Diagnostics & Advice, Groningen, the Netherlands - \url{https://www.certe.nl}

#' @section Read more on our website!:
#' On our website \url{https://msberends.gitlab.io/AMR} you can find \href{https://msberends.gitlab.io/AMR/articles/AMR.html}{a tutorial} about how to conduct AMR analysis, the \href{https://msberends.gitlab.io/AMR/reference}{complete documentation of all functions} (which reads a lot easier than here in R) and \href{https://msberends.gitlab.io/AMR/articles/WHONET.html}{an example analysis using WHONET data}.

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
#  # prevent NOTE on R >= 3.6
#' @importFrom microbenchmark microbenchmark
NULL
