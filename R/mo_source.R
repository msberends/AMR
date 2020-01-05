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

#' Use predefined reference data set
#'
#' @description These functions can be used to predefine your own reference to be used in [as.mo()] and consequently all `mo_*` functions like [mo_genus()] and [mo_gramstain()].
#'
#' This is **the fastest way** to have your organisation (or analysis) specific codes picked up and translated by this package.
#' @inheritSection lifecycle Stable lifecycle
#' @param path location of your reference file, see Details
#' @rdname mo_source
#' @name mo_source
#' @aliases set_mo_source get_mo_source
#' @details The reference file can be a text file seperated with commas (CSV) or tabs or pipes, an Excel file (either 'xls' or 'xlsx' format) or an R object file (extension '.rds'). To use an Excel file, you need to have the `readxl` package installed.
#'
#' [set_mo_source()] will check the file for validity: it must be a [`data.frame`], must have a column named `"mo"` which contains values from [`microorganisms$mo`][microorganisms] and must have a reference column with your own defined values. If all tests pass, [set_mo_source()] will read the file into R and export it to `"~/.mo_source.rds"`. This compressed data file will then be used at default for MO determination (function [as.mo()] and consequently all `mo_*` functions like [mo_genus()] and [mo_gramstain()]). The location of the original file will be saved as option with `options(mo_source = path)`. Its timestamp will be saved with `options(mo_source_datetime = ...)`.
#'
#' [get_mo_source()] will return the data set by reading `"~/.mo_source.rds"` with [readRDS()]. If the original file has changed (the file defined with `path`), it will call [set_mo_source()] to update the data file automatically.
#'
#' Reading an Excel file (`.xlsx`) with only one row has a size of 8-9 kB. The compressed file used by this package will have a size of 0.1 kB and can be read by [get_mo_source()] in only a couple of microseconds (a millionth of a second).
#' 
#' ## How it works
#' 
#' Imagine this data on a sheet of an Excel file (mo codes were looked up in the `microorganisms` data set). The first column contains the organisation specific codes, the second column contains an MO code from this package:
#' ```
#'   |         A          |       B      |
#' --|--------------------|--------------|
#' 1 | Organisation XYZ   | mo           |
#' 2 | lab_mo_ecoli       | B_ESCHR_COLI |
#' 3 | lab_mo_kpneumoniae | B_KLBSL_PNMN |
#' 4 |                    |              |
#' ```
#'
#' We save it as `"home/me/ourcodes.xlsx"`. Now we have to set it as a source:
#' ```
#' set_mo_source("home/me/ourcodes.xlsx")
#' # Created mo_source file '~/.mo_source.rds' from 'home/me/ourcodes.xlsx'.
#' ```
#'
#' It has now created a file `"~/.mo_source.rds"` with the contents of our Excel file, but only the first column with foreign values and the 'mo' column will be kept.
#'
#' And now we can use it in our functions:
#' ```
#' as.mo("lab_mo_ecoli")
#' [1] B_ESCHR_COLI
#'
#' mo_genus("lab_mo_kpneumoniae")
#' [1] "Klebsiella"
#'
#' # other input values still work too
#' as.mo(c("Escherichia coli", "E. coli", "lab_mo_ecoli"))
#' [1] B_ESCHR_COLI B_ESCHR_COLI B_ESCHR_COLI
#' ```
#'
#' If we edit the Excel file to, let's say, by adding row 4 like this:
#' ```
#'   |         A          |       B      |
#' --|--------------------|--------------|
#' 1 | Organisation XYZ   | mo           |
#' 2 | lab_mo_ecoli       | B_ESCHR_COLI |
#' 3 | lab_mo_kpneumoniae | B_KLBSL_PNMN |
#' 4 | lab_Staph_aureus   | B_STPHY_AURS |
#' 5 |                    |              |
#' ```
#'
#' ...any new usage of an MO function in this package will update your data file:
#' ```
#' as.mo("lab_mo_ecoli")
#' # Updated mo_source file '~/.mo_source.rds' from 'home/me/ourcodes.xlsx'.
#' [1] B_ESCHR_COLI
#'
#' mo_genus("lab_Staph_aureus")
#' [1] "Staphylococcus"
#' ```
#'
#' To remove the reference data file completely, just use `""` or `NULL` as input for `[set_mo_source()]`:
#' ```
#' set_mo_source(NULL)
#' # Removed mo_source file '~/.mo_source.rds'.
#' ```
#' @importFrom dplyr select everything
#' @export
#' @inheritSection AMR Read more on our website!
set_mo_source <- function(path) {

  file_location <- path.expand("~/mo_source.rds")

  if (!is.character(path) | length(path) > 1) {
    stop("`path` must be a character of length 1.")
  }

  if (path %in% c(NULL, "")) {
    options(mo_source = NULL)
    options(mo_source_timestamp = NULL)
    if (file.exists(file_location)) {
      unlink(file_location)
      message("Removed mo_source file '", file_location, "'.")
    }
    return(invisible())
  }

  if (!file.exists(path)) {
    stop("File not found: ", path)
  }

  if (path %like% "[.]rds$") {
    df <- readRDS(path)

  } else if (path %like% "[.]xlsx?$") {
    # is Excel file (old or new)
    stopifnot_installed_package("readxl")
    df <- readxl::read_excel(path)

  } else if (path %like% "[.]tsv$") {
    df <- utils::read.table(header = TRUE, sep = "\t", stringsAsFactors = FALSE)

  } else {
    # try comma first
    try(
      df <- utils::read.table(header = TRUE, sep = ",", stringsAsFactors = FALSE),
      silent = TRUE)
    if (!mo_source_isvalid(df)) {
      # try tab
      try(
        df <- utils::read.table(header = TRUE, sep = "\t", stringsAsFactors = FALSE),
        silent = TRUE)
    }
    if (!mo_source_isvalid(df)) {
      # try pipe
      try(
        df <- utils::read.table(header = TRUE, sep = "|", stringsAsFactors = FALSE),
        silent = TRUE)
    }
  }

  if (!mo_source_isvalid(df)) {
    stop("File must contain a column with self-defined values and a reference column `mo` with valid values from the `microorganisms` data set.")
  }

  df <- df %>% filter(!is.na(mo))

  # keep only first two columns, second must be mo
  if (colnames(df)[1] == "mo") {
    df <- df[, c(2, 1)]
  } else {
    df <- df[, c(1, 2)]
  }

  df <- as.data.frame(df, stringAsFactors = FALSE)

  # success
  if (file.exists(file_location)) {
    action <- "Updated"
  } else {
    action <- "Created"
  }
  saveRDS(df, file_location)
  options(mo_source = path)
  options(mo_source_timestamp = as.character(file.info(path)$mtime))
  message(action, " mo_source file '", file_location, "' from '", path, "'.")
}

#' @rdname mo_source
#' @export
get_mo_source <- function() {
  if (is.null(getOption("mo_source", NULL))) {
    NULL
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
    file_location <- path.expand("~/mo_source.rds")
    readRDS(file_location)
  }
}

mo_source_isvalid <- function(x) {
  if (deparse(substitute(x)) == "get_mo_source()") {
    return(TRUE)
  }
  if (identical(x, get_mo_source())) {
    return(TRUE)
  }
  if (is.null(x)) {
    return(TRUE)
  }
  if (!is.data.frame(x)) {
    return(FALSE)
  }
  if (!"mo" %in% colnames(x)) {
    return(FALSE)
  }
  all(x$mo %in% c("", AMR::microorganisms$mo, microorganisms.translation$mo_old), na.rm = TRUE)
}
