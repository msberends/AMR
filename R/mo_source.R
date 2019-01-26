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

#' Use predefined reference data set
#'
#' These functions can be used to predefine your own reference to be used in \code{\link{as.mo}} and consequently all \code{mo_*} functions like \code{\link{mo_genus}} and \code{\link{mo_gramstain}}.
#' @param path location of your reference file, see Details
#' @rdname mo_source
#' @name mo_source
#' @aliases set_mo_source get_mo_source
#' @details The reference file can be a text file seperated with commas (CSV) or pipes, an Excel file (old 'xls' format or new 'xlsx' format) or an R object file (extension '.rds'). To use an Excel file, you need to have the \code{readxl} package installed.
#'
#' \code{set_mo_source} will check the file for validity: it must be a \code{data.frame}, must have a column named \code{"mo"} which contains values from \code{microorganisms$mo} and must have a reference column with your own defined values. If all tests pass, \code{set_mo_source} will read the file into R and export it to \code{"~/.mo_source.rds"}. This compressed data file will then be used at default for MO determination (function \code{\link{as.mo}} and consequently all \code{mo_*} functions like \code{\link{mo_genus}} and \code{\link{mo_gramstain}}). The location of the original file will be saved as option with \code{\link{options}(mo_source = path)}. Its timestamp will be saved with \code{\link{options}(mo_source_datetime = ...)}.
#'
#' \code{get_mo_source} will return the data set by reading \code{"~/.mo_source.rds"} with \code{\link{readRDS}}. If the original file has changed (the file defined with \code{path}), it will call \code{set_mo_source} to update the data file automatically.
#'
#' Reading an Excel file (\code{.xlsx}) with only one row has a size of 8-9 kB. The compressed file will have a size of 0.1 kB and can be read by \code{get_mo_source} in only a couple of microseconds (a millionth of a second).
#' @importFrom dplyr select everything
#' @export
#' @inheritSection AMR Read more on our website!
#' @examples
#' \dontrun{
#'
#' # imagine this Excel file (mo codes looked up in `microorganisms` data set):
#' #            A                B
#' # 1  our code            mo
#' # 2  lab_mo_ecoli        B_ESCHR_COL
#' # 3  lab_mo_kpneumoniae  B_KLBSL_PNE
#'
#' # 1. We save it as 'home/me/ourcodes.xlsx'
#'
#' # 2. We use it for input:
#' set_mo_source("C:\path\ourcodes.xlsx")
#' #> Created mo_source file '~/.mo_source.rds' from 'home/me/ourcodes.xlsx'.
#'
#' # 3. And use it in our functions:
#' as.mo("lab_mo_ecoli")
#' #> B_ESCHR_COL
#'
#' mo_genus("lab_mo_kpneumoniae")
#' #> "Klebsiella"
#'
#' # 4. It will look for changes itself:
#' # (add new row to the Excel file and save it)
#'
#' mo_genus("lab_mo_kpneumoniae")
#' #> Updated mo_source file '~/.mo_source.rds' from 'home/me/ourcodes.xlsx'.
#' #> "Klebsiella"
#' }
set_mo_source <- function(path) {

  if (!is.character(path) | length(path) > 1) {
    stop("`path` must be a character of length 1.")
  }

  if (path == "") {
    options(mo_source = NULL)
    options(mo_source_timestamp = NULL)
    if (file.exists("~/.mo_source.rds")) {
      unlink("~/.mo_source.rds")
      message("Removed mo_source file '~/.mo_source.rds'.")
    }
    return(invisible())
  }

  if (!file.exists(path)) {
    stop("File not found: ", path)
  }

  is_valid <- function(df) {
    valid <- TRUE
    if (!is.data.frame(df)) {
      valid <- FALSE
    } else if (!"mo" %in% colnames(df)) {
      valid <- FALSE
    } else if (!all(df$mo %in% AMR::microorganisms$mo)) {
      valid <- FALSE
    } else if (NCOL(df) < 2) {
      valid <- FALSE
    }
    valid
  }

  if (path %like% '[.]rds$') {
    df <- readRDS(path)

  } else if (path %like% '[.]xlsx?$') {
    # is Excel file (old or new)
    if (!"readxl" %in% utils::installed.packages()) {
      stop("Install the 'readxl' package first.")
    }
    if (path %like% '[.]xlsx$') {
      df <- readxl::read_xlsx(path)
    } else {
      df <- readxl::read_xls(path)
    }

  } else {
    # try comma first
    try(
      df <- utils::read.table(header = TRUE, sep = ",", stringsAsFactors = FALSE),
      silent = TRUE)
    if (!is_valid(df)) {
      # try pipe
      try(
        df <- utils::read.table(header = TRUE, sep = "|", stringsAsFactors = FALSE),
        silent = TRUE)
    }
  }

  if (!is_valid(df)) {
    stop("File must contain a column with self-defined values and a reference column `mo` with valid values from the `microorganisms` data set.")
  }

  if (colnames(df)[1] == "mo") {
    # put mo to the end
    df <- df %>% select(-"mo", everything(), "mo")
  }

  df <- as.data.frame(df, stringAsFactors = FALSE)

  # success
  if (file.exists("~/.mo_source.rds")) {
    action <- "Updated"
  } else {
    action <- "Created"
  }
  saveRDS(df, "~/.mo_source.rds")
  options(mo_source = path)
  options(mo_source_timestamp = as.character(file.info(path)$mtime))
  message(action, " mo_source file '~/.mo_source.rds' from '", path, "'.")
}

#' @rdname mo_source
#' @export
get_mo_source <- function() {

  if (is.null(getOption("mo_source", NULL))) {
    return(NULL)
  } else {
    old_time <- as.POSIXct(getOption("mo_source_timestamp"))
    new_time <- as.POSIXct(as.character(file.info(getOption("mo_source", ""))$mtime))

    if (is.na(new_time)) {
      # source file was deleted, remove reference too
      set_mo_source("")
      return(NULL)
    }
    if (new_time != old_time) {
      # set updated source
      set_mo_source(getOption("mo_source"))
    }
  }

  readRDS("~/.mo_source.rds")
}