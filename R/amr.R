# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Analysis                              #
#                                                                      #
# SOURCE                                                               #
# https://gitlab.com/msberends/AMR                                     #
#                                                                      #
# LICENCE                                                              #
# (c) 2018-2020 Berends MS, Luz CF et al.                              #
#                                                                      #
# This R package is free software; you can freely use and distribute   #
# it for both personal and commercial purposes under the terms of the  #
# GNU General Public License version 2.0 (GNU GPL-2), as published by  #
# the Free Software Foundation.                                        #
#                                                                      #
# We created this package for both routine data analysis and academic  #
# research and it was publicly released in the hope that it will be    #
# useful, but it comes WITHOUT ANY WARRANTY OR LIABILITY.              #
# Visit our website for more info: https://msberends.gitlab.io/AMR.    #
# ==================================================================== #

#' The `AMR` Package
#'
#' Welcome to the `AMR` package.
#' @details
#' `AMR` is a free and open-source R package to simplify the analysis and prediction of Antimicrobial Resistance (AMR) and to work with microbial and antimicrobial properties by using evidence-based methods. It supports any table format, including WHONET/EARS-Net data.
#'
#' We created this package for both academic research and routine analysis at the Faculty of Medical Sciences of the University of Groningen and the Medical Microbiology & Infection Prevention (MMBI) department of the University Medical Center Groningen (UMCG). This R package is actively maintained and free software; you can freely use and distribute it for both personal and commercial (but not patent) purposes under the terms of the GNU General Public License version 2.0 (GPL-2), as published by the Free Software Foundation.
#'
#' This package can be used for:
#' - Reference for the taxonomy of microorganisms, since the package contains all microbial (sub)species from the [Catalogue of Life](http://www.catalogueoflife.org)
#' - Interpreting raw MIC and disk diffusion values, based on the latest CLSI or EUCAST guidelines
#' - Determining first isolates to be used for AMR analysis
#' - Calculating antimicrobial resistance
#' - Determining multi-drug resistance (MDR) / multi-drug resistant organisms (MDRO)
#' - Calculating (empirical) susceptibility of both mono therapy and combination therapies
#' - Predicting future antimicrobial resistance using regression models
#' - Getting properties for any microorganism (like Gram stain, species, genus or family)
#' - Getting properties for any antibiotic (like name, EARS-Net code, ATC code, PubChem code, defined daily dose or trade name)
#' - Plotting antimicrobial resistance
#' - Applying EUCAST expert rules

#' @section Read more on our website!:
#' On our website <https://msberends.gitlab.io/AMR> you can find [a tutorial](https://msberends.gitlab.io/AMR/articles/AMR.html) about how to conduct AMR analysis, the [complete documentation of all functions](https://msberends.gitlab.io/AMR/reference) (which reads a lot easier than here in R) and [an example analysis using WHONET data](https://msberends.gitlab.io/AMR/articles/WHONET.html).
#' @section Contact us:
#' For suggestions, comments or questions, please contact us at:
#'
#' Matthijs S. Berends \cr
#' m.s.berends at umcg dot nl \cr
#' Department of Medical Microbiology, University of Groningen \cr
#' University Medical Center Groningen \cr
#' Post Office Box 30001 \cr
#' 9700 RB Groningen \cr
#' The Netherlands
#'
#' If you have found a bug, please file a new issue at: \cr
#' <https://gitlab.com/msberends/AMR/issues>
#' @name AMR
#' @rdname AMR
#' @importFrom microbenchmark microbenchmark
#' @importFrom knitr kable
NULL
