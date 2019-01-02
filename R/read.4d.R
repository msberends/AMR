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

#' Read data from 4D database
#'
#' This function is only useful for the MMB department of the UMCG. Use this function to \strong{import data by just defining the \code{file} parameter}. It will automatically transform birth dates and calculate patients age, translate the column names to English, transform the MO codes with \code{\link{as.mo}} and transform all antimicrobial columns with \code{\link{as.rsi}}.
#' @inheritParams utils::read.table
#' @param info a logical to indicate whether info about the import should be printed, defaults to \code{TRUE} in interactive sessions
#' @details Column names will be transformed, but the original column names are set as a "label" attribute and can be seen in e.g. RStudio Viewer.
#' @inheritSection AMR Read more on our website!
#' @export
read.4D <- function(file,
                    info = interactive(),
                    header = TRUE,
                    row.names = NULL,
                    sep = "\t",
                    quote = "\"'",
                    dec = ",",
                    na.strings = c("NA", "", "."),
                    skip = 2,
                    check.names = TRUE,
                    strip.white = TRUE,
                    fill = TRUE,
                    blank.lines.skip = TRUE,
                    stringsAsFactors = FALSE,
                    fileEncoding = "UTF-8",
                    encoding = "UTF-8") {

  if (info == TRUE) {
    message("Importing ", file, "... ", appendLF = FALSE)
  }
  data_4D <- utils::read.table(file = file,
                               row.names = row.names,
                               header = header,
                               sep = sep,
                               quote = quote,
                               dec = dec,
                               na.strings = na.strings,
                               skip = skip,
                               check.names = check.names,
                               strip.white = strip.white,
                               fill = fill,
                               blank.lines.skip = blank.lines.skip,
                               stringsAsFactors = stringsAsFactors,
                               fileEncoding = fileEncoding,
                               encoding = encoding)

  # helper function for dates
  to_date_4D <- function(x) {
    date_regular <- as.Date(x, format = "%d-%m-%y")
    posixlt <- as.POSIXlt(date_regular)
    # born after today will be born 100 years ago
    # based on https://stackoverflow.com/a/3312971/4575331
    posixlt[date_regular > Sys.Date() & !is.na(posixlt)]$year <- posixlt[date_regular > Sys.Date() & !is.na(posixlt)]$year - 100
    as.Date(posixlt)
  }

  if (info == TRUE) {
    message("OK\nTransforming column names... ", appendLF = FALSE)
  }
  if ("row.names" %in% colnames(data_4D) & all(is.na(data_4D[, ncol(data_4D)]))) {
    # remove first column name "row.names" and remove last empty column
    colnames(data_4D) <- c(colnames(data_4D)[2:ncol(data_4D)], "_skip_last")
    data_4D <- data_4D[, -ncol(data_4D)]
  }

  colnames(data_4D) <- tolower(colnames(data_4D))
  if (all(c("afnamedat", "gebdatum") %in% colnames(data_4D))) {
    # add age
    data_4D$age <- NA_integer_
  }
  cols_wanted <- c("patientnr", "gebdatum", "age", "mv", "monsternr", "afnamedat", "bepaling",
                   "afd.", "spec", "mat", "matbijz.",  "mocode",
                   "amfo", "amox", "anid", "azit", "casp", "cecl", "cefe", "cfcl",
                   "cfot", "cfox", "cfta", "cftr", "cfur", "chlo", "cipr", "clin",
                   "cocl", "ctta", "dapt", "doxy", "eryt", "fluo", "fluz", "fosf",
                   "fusi", "gehi", "gent", "imip", "kana", "levo", "line", "mero",
                   "metr", "mico", "mino", "moxi", "mupi", "nali", "nitr", "norf",
                   "oxac", "peni", "pipe", "pita", "poly", "posa", "quda", "rifa",
                   "spat", "teic", "tige", "tobr", "trim", "trsu", "vana", "vanb",
                   "vanc", "vori")
  # this ones actually exist
  cols_wanted <- cols_wanted[cols_wanted %in% colnames(data_4D)]
  # order of columns
  data_4D <- data_4D[, cols_wanted]

  # backup original column names
  colnames.bak <- toupper(colnames(data_4D))
  colnames.bak[colnames.bak == "AGE"] <- NA_character_

  # rename of columns
  colnames(data_4D) <- gsub("patientnr", "patient_id", colnames(data_4D), fixed = TRUE)
  colnames(data_4D) <- gsub("gebdatum", "date_birth", colnames(data_4D), fixed = TRUE)
  colnames(data_4D) <- gsub("mv", "gender", colnames(data_4D), fixed = TRUE)
  colnames(data_4D) <- gsub("monsternr", "sample_id", colnames(data_4D), fixed = TRUE)
  colnames(data_4D) <- gsub("afnamedat", "date_received", colnames(data_4D), fixed = TRUE)
  colnames(data_4D) <- gsub("bepaling", "sample_test", colnames(data_4D), fixed = TRUE)
  colnames(data_4D) <- gsub("afd.", "department", colnames(data_4D), fixed = TRUE)
  colnames(data_4D) <- gsub("spec", "specialty", colnames(data_4D), fixed = TRUE)
  colnames(data_4D) <- gsub("matbijz.", "specimen_type", colnames(data_4D), fixed = TRUE)
  colnames(data_4D) <- gsub("mat", "specimen_group", colnames(data_4D), fixed = TRUE)
  colnames(data_4D) <- gsub("mocode", "mo", colnames(data_4D), fixed = TRUE)

  if (info == TRUE) {
    message("OK\nTransforming dates and age... ", appendLF = FALSE)
  }
  if ("date_birth" %in% colnames(data_4D)) {
    data_4D$date_birth <- to_date_4D(data_4D$date_birth)

  }
  if ("date_received" %in% colnames(data_4D)) {
    data_4D$date_received <- to_date_4D(data_4D$date_received)
  }
  if ("age" %in% colnames(data_4D)) {
    data_4D$age <- age(data_4D$date_birth, data_4D$date_received)
  }
  if ("gender" %in% colnames(data_4D)) {
    data_4D$gender[data_4D$gender == "V"] <- "F"
  }

  if (info == TRUE) {
    message("OK\nTransforming MO codes... ", appendLF = FALSE)
  }
  if ("mo" %in% colnames(data_4D)) {
    data_4D$mo <- as.mo(data_4D$mo)
    # column right of mo is:
    drug1 <- colnames(data_4D)[grep("^mo$", colnames(data_4D)) + 1]
    if (!is.na(drug1)) {
      # and last is:
      drug_last <- colnames(data_4D)[length(data_4D)]
      # transform those to rsi:
      data_4D <- suppressWarnings(mutate_at(data_4D, vars(drug1:drug_last), as.rsi))
    }
  }

  # set original column names as label (can be seen in RStudio Viewer)
  if (info == TRUE) {
    message("OK\nSetting original column names as label... ", appendLF = FALSE)
  }
  for (i in 1:ncol(data_4D)) {
    if (!is.na(colnames.bak[i])) {
      attr(data_4D[, i], "label") <- colnames.bak[i]
    }
  }

  if (info == TRUE) {
    message("OK\nSetting query as label to data.frame... ", appendLF = FALSE)
  }
  qry <- readLines(con <- file(file, open="r"))[1]
  close(con)
  attr(data_4D, "label") <- qry

  if (info == TRUE) {
    message("OK")
  }

  data_4D
}

