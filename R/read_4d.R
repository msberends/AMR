# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Analysis                              #
#                                                                      #
# AUTHORS                                                              #
# Berends MS (m.s.berends@umcg.nl), Luz CF (c.f.luz@umcg.nl)           #
#                                                                      #
# LICENCE                                                              #
# This program is free software; you can redistribute it and/or modify #
# it under the terms of the GNU General Public License version 2.0,    #
# as published by the Free Software Foundation.                        #
#                                                                      #
# This program is distributed in the hope that it will be useful,      #
# but WITHOUT ANY WARRANTY; without even the implied warranty of       #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        #
# GNU General Public License for more details.                         #
# ==================================================================== #

#' Read data from 4D database
#'
#' This function is only useful for the MMB department of the UMCG. Use this function to \strong{import data by just defining the \code{file} parameter}. It will automatically transform birth dates and calculate patients age, translate the data set to English, transform the \code{mo} with \code{\link{as.mo}} and transform all antimicrobial columns with \code{\link{as.rsi}}.
#' @inheritParams utils::read.table
#' @export
read_4D <- function(file,
                    header = TRUE,
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

  data_4D <- utils::read.table(file = file,
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

  # helper functions
  to_date_4D <- function(x) {
    date_regular <- as.Date(x, format = "%d-%m-%y")
    posixlt <- as.POSIXlt(date_regular)
    # born after today will be born 100 years ago
    # based on https://stackoverflow.com/a/3312971/4575331
    posixlt[date_regular > Sys.Date()]$year <- posixlt[date_regular > Sys.Date()]$year - 100
    as.Date(posixlt)
  }
  to_age_4D <- function(from, to) {
    from_lt = as.POSIXlt(from)
    to_lt = as.POSIXlt(to)

    age = to_lt$year - from_lt$year

    ifelse(to_lt$mon < from_lt$mon |
             (to_lt$mon == from_lt$mon & to_lt$mday < from_lt$mday),
           age - 1, age)
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

  if ("date_birth" %in% colnames(data_4D)) {
    data_4D$date_birth <- to_date_4D(data_4D$date_birth)

  }
  if ("date_received" %in% colnames(data_4D)) {
    data_4D$date_received <- to_date_4D(data_4D$date_received)
  }
  if ("age" %in% colnames(data_4D)) {
    data_4D$age <- to_age_4D(data_4D$date_birth, data_4D$date_received)
  }
  if ("gender" %in% colnames(data_4D)) {
    data_4D$gender[data_4D$gender == "V"] <- "F"
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

  data_4D
}

