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

#' Key antibiotics for first \emph{weighted} isolates
#'
#' These function can be used to determine first isolates (see \code{\link{first_isolate}}). Using key antibiotics to determine first isolates is more reliable than without key antibiotics. These selected isolates will then be called first \emph{weighted} isolates.
#' @param tbl table with antibiotics coloms, like \code{amox} and \code{amcl}.
#' @param x,y characters to compare
#' @inheritParams first_isolate
#' @param universal_1,universal_2,universal_3,universal_4,universal_5,universal_6 column names of \strong{broad-spectrum} antibiotics, case-insensitive
#' @param GramPos_1,GramPos_2,GramPos_3,GramPos_4,GramPos_5,GramPos_6 column names of antibiotics for \strong{Gram positives}, case-insensitive
#' @param GramNeg_1,GramNeg_2,GramNeg_3,GramNeg_4,GramNeg_5,GramNeg_6 column names of antibiotics for \strong{Gram negatives}, case-insensitive
#' @param warnings give warning about missing antibiotic columns, they will anyway be ignored
#' @details The function \code{key_antibiotics} returns a character vector with 12 antibiotic results for every isolate. These isolates can then be compared using \code{key_antibiotics_equal}, to check if two isolates have generally the same antibiogram. Missing and invalid values are replaced with a dot (\code{"."}). The \code{\link{first_isolate}} function only uses this function on the same microbial species from the same patient. Using this, an MRSA will be included after a susceptible \emph{S. aureus} (MSSA) found within the same episode (see \code{episode} parameter of \code{\link{first_isolate}}). Without key antibiotic comparison it wouldn't.
#'
#'   At default, the antibiotics that are used for \strong{Gram positive bacteria} are (colum names): \cr
#'   \code{"amox"}, \code{"amcl"}, \code{"cfur"}, \code{"pita"}, \code{"cipr"}, \code{"trsu"} (until here is universal), \code{"vanc"}, \code{"teic"}, \code{"tetr"}, \code{"eryt"}, \code{"oxac"}, \code{"rifa"}.
#'
#'   At default, the antibiotics that are used for \strong{Gram negative bacteria} are (colum names): \cr
#'   \code{"amox"}, \code{"amcl"}, \code{"cfur"}, \code{"pita"}, \code{"cipr"}, \code{"trsu"} (until here is universal), \code{"gent"}, \code{"tobr"}, \code{"coli"}, \code{"cfot"}, \code{"cfta"}, \code{"mero"}.
#'
#'
#'   The function \code{key_antibiotics_equal} checks the characters returned by \code{key_antibiotics} for equality, and returns a logical vector.
#' @inheritSection first_isolate Key antibiotics
#' @rdname key_antibiotics
#' @export
#' @importFrom dplyr %>% mutate if_else
#' @seealso \code{\link{first_isolate}}
#' @examples
#' \dontrun{
#' # set key antibiotics to a new variable
#' tbl$keyab <- key_antibiotics(tbl)
#'
#' # add regular first isolates
#' tbl$first_isolate <-
#'   first_isolate(tbl)
#'
#' # add first WEIGHTED isolates using key antibiotics
#' tbl$first_isolate_weighed <-
#'   first_isolate(tbl,
#'                 col_keyantibiotics = 'keyab')
#' }
#'
#' # output of the `key_antibiotics` function could be like this:
#' strainA <- "SSSRR.S.R..S"
#' strainB <- "SSSIRSSSRSSS"
#'
#' key_antibiotics_equal(strainA, strainB)
#' # TRUE, because I is ignored (as are missing values)
#'
#' key_antibiotics_equal(strainA, strainB, ignore_I = FALSE)
#' # FALSE, because I is not ignored and so the 4th value differs
key_antibiotics <- function(tbl,
                            col_bactid = "bactid",
                            universal_1 = "amox",
                            universal_2 = "amcl",
                            universal_3 = "cfur",
                            universal_4 = "pita",
                            universal_5 = "cipr",
                            universal_6 = "trsu",
                            GramPos_1 = "vanc",
                            GramPos_2 = "teic",
                            GramPos_3 = "tetr",
                            GramPos_4 = "eryt",
                            GramPos_5 = "oxac",
                            GramPos_6 = "rifa",
                            GramNeg_1 = "gent",
                            GramNeg_2 = "tobr",
                            GramNeg_3 = "coli",
                            GramNeg_4 = "cfot",
                            GramNeg_5 = "cfta",
                            GramNeg_6 = "mero",
                            warnings = TRUE) {

  if (!col_bactid %in% colnames(tbl)) {
    stop('Column ', col_bactid, ' not found.', call. = FALSE)
  }

  # check columns
  col.list <- c(universal_1, universal_2, universal_3, universal_4, universal_5, universal_6,
                GramPos_1, GramPos_2, GramPos_3, GramPos_4, GramPos_5, GramPos_6,
                GramNeg_1, GramNeg_2, GramNeg_3, GramNeg_4, GramNeg_5, GramNeg_6)
  col.list <- check_available_columns(tbl = tbl, col.list = col.list, info = warnings)
  universal_1 <- col.list[universal_1]
  universal_2 <- col.list[universal_2]
  universal_3 <- col.list[universal_3]
  universal_4 <- col.list[universal_4]
  universal_5 <- col.list[universal_5]
  universal_6 <- col.list[universal_6]
  GramPos_1 <- col.list[GramPos_1]
  GramPos_2 <- col.list[GramPos_2]
  GramPos_3 <- col.list[GramPos_3]
  GramPos_4 <- col.list[GramPos_4]
  GramPos_5 <- col.list[GramPos_5]
  GramPos_6 <- col.list[GramPos_6]
  GramNeg_1 <- col.list[GramNeg_1]
  GramNeg_2 <- col.list[GramNeg_2]
  GramNeg_3 <- col.list[GramNeg_3]
  GramNeg_4 <- col.list[GramNeg_4]
  GramNeg_5 <- col.list[GramNeg_5]
  GramNeg_6 <- col.list[GramNeg_6]

  universal <- c(universal_1, universal_2, universal_3,
                 universal_4, universal_5, universal_6)

  gram_positive = c(universal,
                    GramPos_1, GramPos_2, GramPos_3,
                    GramPos_4, GramPos_5, GramPos_6)
  gram_positive <- gram_positive[!is.na(gram_positive)]

  gram_negative = c(universal,
                    GramNeg_1, GramNeg_2, GramNeg_3,
                    GramNeg_4, GramNeg_5, GramNeg_6)
  gram_negative <- gram_negative[!is.na(gram_negative)]

  # join microorganisms
  tbl <- tbl %>% left_join_microorganisms(col_bactid)

  tbl$key_ab <- NA_character_

  # Gram +
  tbl <- tbl %>% mutate(key_ab =
                          if_else(gramstain %like% '^Positive ',
                                  apply(X = tbl[, gram_positive],
                                        MARGIN = 1,
                                        FUN = function(x) paste(x, collapse = "")),
                                  key_ab))

  # Gram -
  tbl <- tbl %>% mutate(key_ab =
                          if_else(gramstain %like% '^Negative ',
                                  apply(X = tbl[, gram_negative],
                                        MARGIN = 1,
                                        FUN = function(x) paste(x, collapse = "")),
                                  key_ab))

  # format
  key_abs <- tbl %>%
    pull(key_ab) %>%
    gsub('(NA|NULL)', '.', .) %>%
    gsub('[^SIR]', '.', ., ignore.case = TRUE)

  key_abs

}

#' @importFrom dplyr progress_estimated %>%
#' @rdname key_antibiotics
#' @export
key_antibiotics_equal <- function(x,
                                  y,
                                  type = c("keyantibiotics", "points"),
                                  ignore_I = TRUE,
                                  points_threshold = 2,
                                  info = FALSE) {
  # x is active row, y is lag
  type <- type[1]

  if (length(x) != length(y)) {
    stop('Length of `x` and `y` must be equal.')
  }

  result <- logical(length(x))

  if (type == "keyantibiotics") {
    if (ignore_I == TRUE) {
      # evaluation using regular expression will treat '.' as any character
      # so I is actually ignored then
      x <- gsub('I', '.', x, ignore.case = TRUE)
      y <- gsub('I', '.', y, ignore.case = TRUE)
    }

    for (i in 1:length(x)) {
      x_split <- strsplit(x[i], "")[[1]]
      y_split <- strsplit(y[i], "")[[1]]
      y_split[x_split == "."] <- "."
      x_split[y_split == "."] <- "."
      x_checkfor <- paste(x_split, collapse = "")
      y_checkfor <- paste(y_split, collapse = "")
      result[i] <- nchar(x[i]) == nchar(y[i]) &
        (x_checkfor %like% y_checkfor |
           y_checkfor %like% x_checkfor)
    }
    return(result)

  } else {

    if (type != 'points') {
      stop('`', type, '` is not a valid value for type, must be "points" or "keyantibiotics". See ?first_isolate.')
    }

    if (info == TRUE) {
      p <- dplyr::progress_estimated(length(x))
    }

    for (i in 1:length(x)) {

      if (info == TRUE) {
        p$tick()$print()
      }

      if (is.na(x[i])) {
        x[i] <- ''
      }
      if (is.na(y[i])) {
        y[i] <- ''
      }

      if (nchar(x[i]) != nchar(y[i])) {

        result[i] <- FALSE

      } else if (x[i] == '' & y[i] == '') {

        result[i] <- TRUE

      } else {

        x2 <- strsplit(x[i], "")[[1]]
        y2 <- strsplit(y[i], "")[[1]]

        # count points for every single character:
        # - no change is 0 points
        # - I <-> S|R is 0.5 point
        # - S|R <-> R|S is 1 point
        # use the levels of as.rsi (S = 1, I = 2, R = 3)

        suppressWarnings(x2 <- x2 %>% as.rsi() %>% as.double())
        suppressWarnings(y2 <- y2 %>% as.rsi() %>% as.double())

        points <- (x2 - y2) %>% abs() %>% sum(na.rm = TRUE)
        result[i] <- ((points / 2) >= points_threshold)
      }
    }
    if (info == TRUE) {
      cat('\n')
    }
    result
  }
}
