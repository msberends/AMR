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

#' `vctrs` methods
#' 
#' These methods are needed to support methods used by the tidyverse, like joining and transforming data, with new classes that come with this package.
#' @inheritSection lifecycle Stable lifecycle
#' @inheritSection AMR Read more on our website!
#' @keywords internal
#' @name AMR-vctrs
NULL


# Class mo ----------------------------------------------------------------


#' @exportMethod vec_ptype_abbr.mo
#' @importFrom vctrs vec_ptype_abbr
#' @export
vec_ptype_abbr.mo <- function(x, ...) {
  "mo"
}

#' @exportMethod vec_ptype_full.mo
#' @importFrom vctrs vec_ptype_full
#' @export
vec_ptype_full.mo <- function(x, ...) {
  "mo"
}

#' @rdname AMR-vctrs
#' @export
vec_ptype2.mo <- function(x, y, ...) {
  UseMethod("vec_ptype2.mo", y)
}

#' @method vec_ptype2.mo default
#' @export
vec_ptype2.mo.default <- function(x, y, ..., x_arg = "x", y_arg = "y") {
  vctrs::vec_default_ptype2(x, y, x_arg = x_arg, y_arg = y_arg)
}

#' @method vec_ptype2.mo character
#' @export
vec_ptype2.mo.character <- function(x, y, ...) {
  x
}

#' @method vec_ptype2.character mo
#' @importFrom vctrs vec_ptype2.character
#' @export
vec_ptype2.character.mo <- function(x, y, ...) {
  y
}

#' @rdname AMR-vctrs
#' @export
vec_cast.mo <- function(x, to, ...) {
  UseMethod("vec_cast.mo")
}

#' @method vec_cast.mo mo
#' @export
vec_cast.mo.mo <- function(x, to, ...) {
  as.mo(x)
}

#' @method vec_cast.mo character
#' @export
vec_cast.mo.character <- function(x, to, ...) {
  as.mo(x)
}

#' @method vec_cast.mo default
#' @export
vec_cast.mo.default <- function(x, to, ...) {
  vec_default_cast(x, to)
}

# @method vec_cast.character mo
#' @exportMethod vec_cast.character.mo
#' @importFrom vctrs vec_cast
#' @export
vec_cast.character.mo <- function(x, to, ...) {
  # purrr::map_chr(x, stringr::str_c, collapse = " ")
  unclass(x)
}


# Class ab ----------------------------------------------------------------


#' @exportMethod vec_ptype_abbr.ab
#' @importFrom vctrs vec_ptype_abbr
#' @export
vec_ptype_abbr.ab <- function(x, ...) {
  "ab"
}

#' @exportMethod vec_ptype_full.ab
#' @importFrom vctrs vec_ptype_full
#' @export
vec_ptype_full.ab <- function(x, ...) {
  "ab"
}

#' @rdname AMR-vctrs
#' @export
vec_ptype2.ab <- function(x, y, ...) {
  UseMethod("vec_ptype2.ab", y)
}

#' @method vec_ptype2.ab default
#' @export
vec_ptype2.ab.default <- function(x, y, ..., x_arg = "x", y_arg = "y") {
  vctrs::vec_default_ptype2(x, y, x_arg = x_arg, y_arg = y_arg)
}

#' @method vec_ptype2.ab character
#' @export
vec_ptype2.ab.character <- function(x, y, ...) {
  x
}

#' @method vec_ptype2.character ab
#' @importFrom vctrs vec_ptype2.character
#' @export
vec_ptype2.character.ab <- function(x, y, ...) {
  y
}

#' @rdname AMR-vctrs
#' @export
vec_cast.ab <- function(x, to, ...) {
  UseMethod("vec_cast.ab")
}

#' @method vec_cast.ab ab
#' @export
vec_cast.ab.ab <- function(x, to, ...) {
  as.ab(x)
}

#' @method vec_cast.ab character
#' @export
vec_cast.ab.character <- function(x, to, ...) {
  as.ab(x)
}

#' @method vec_cast.ab default
#' @export
vec_cast.ab.default <- function(x, to, ...) {
  vec_default_cast(x, to)
}

# @method vec_cast.character ab
#' @exportMethod vec_cast.character.ab
#' @importFrom vctrs vec_cast
#' @export
vec_cast.character.ab <- function(x, to, ...) {
  # purrr::map_chr(x, stringr::str_c, collapse = " ")
  unclass(x)
}


# Class disk --------------------------------------------------------------


#' @exportMethod vec_ptype_abbr.disk
#' @importFrom vctrs vec_ptype_abbr
#' @export
vec_ptype_abbr.disk <- function(x, ...) {
  "disk"
}

#' @exportMethod vec_ptype_full.disk
#' @importFrom vctrs vec_ptype_full
#' @export
vec_ptype_full.disk <- function(x, ...) {
  "disk"
}


# Class rsi --------------------------------------------------------------


#' @exportMethod vec_ptype_abbr.rsi
#' @importFrom vctrs vec_ptype_abbr
#' @export
vec_ptype_abbr.rsi <- function(x, ...) {
  "rsi"
}

#' @exportMethod vec_ptype_full.rsi
#' @importFrom vctrs vec_ptype_full
#' @export
vec_ptype_full.rsi <- function(x, ...) {
  "rsi"
}


# Class mic --------------------------------------------------------------


#' @exportMethod vec_ptype_abbr.mic
#' @importFrom vctrs vec_ptype_abbr
#' @export
vec_ptype_abbr.mic <- function(x, ...) {
  "mic"
}

#' @exportMethod vec_ptype_full.mic
#' @importFrom vctrs vec_ptype_full
#' @export
vec_ptype_full.mic <- function(x, ...) {
  "mic"
}
