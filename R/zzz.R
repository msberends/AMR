.onLoad <- function(libname, pkgname) {
  backports::import(pkgname)
}

#' @importFrom Rcpp evalCpp
#' @useDynLib AMR, .registration = TRUE
NULL
