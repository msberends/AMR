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

#' Extended functions
#'
#' These functions are extensions of functions in other packages.
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
