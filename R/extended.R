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

#' Extended functions
#'
#' These functions are extensions of functions in other packages.
#' @inheritSection lifecycle Stable lifecycle
#' @inheritSection AMR Read more on our website!
#' @export
#' @keywords internal
#' @name extended-functions
#' @rdname extended-functions
#' @exportMethod scale_type.mo
#' @export
scale_type.mo <- function(x) {
  # fix for:
  # "Don't know how to automatically pick scale for object of type mo. Defaulting to continuous."
  # "Error: Discrete value supplied to continuous scale"
  "discrete"
}

#' @rdname extended-functions
#' @exportMethod scale_type.ab
#' @export
scale_type.ab <- function(x) {
  # fix for:
  # "Don't know how to automatically pick scale for object of type mo. Defaulting to continuous."
  # "Error: Discrete value supplied to continuous scale"
  "discrete"
}
