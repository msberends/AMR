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
#' @inheritParams first_isolate
#' @param amcl,amox,cfot,cfta,cfur,cipr,coli,eryt,gent,mero,oxac,pita,rifa,teic,tetr,tobr,trsu,vanc column names of antibiotics, case-insensitive
#' @details The function \code{key_antibiotics} returns a character vector with antibiotic results.
#'
#'   The antibiotics that are used for \strong{Gram positive bacteria} are (colum names): \cr
#'   amox, amcl, cfur, pita, cipr, trsu, vanc, teic, tetr, eryt, oxac, rifa.
#'
#'   The antibiotics that are used for \strong{Gram negative bacteria} are (colum names): \cr
#'   amox, amcl, cfur, pita, cipr, trsu, gent, tobr, coli, cfot, cfta, mero.
#'
#'
#'   The function \code{key_antibiotics_equal} checks the characters returned by \code{key_antibiotics} for equality, and returns a logical value.
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
key_antibiotics <- function(tbl,
                            col_bactid = "bactid",
                            amcl = "amcl",
                            amox = "amox",
                            cfot = "cfot",
                            cfta = "cfta",
                            cfur = "cfur",
                            cipr = "cipr",
                            coli = "coli",
                            eryt = "eryt",
                            gent = "gent",
                            mero = "mero",
                            oxac = "oxac",
                            pita = "pita",
                            rifa = "rifa",
                            teic = "teic",
                            tetr = "tetr",
                            tobr = "tobr",
                            trsu = "trsu",
                            vanc = "vanc",
                            info = TRUE) {

  if (!col_bactid %in% colnames(tbl)) {
    stop('Column ', col_bactid, ' not found.', call. = FALSE)
  }

  # check columns
  col.list <- c(amcl, amox, cfot, cfta, cfur, cipr,
                coli, eryt, gent, mero, oxac, pita,
                rifa, teic, tetr, tobr, trsu, vanc)
  col.list <- check_available_columns(tbl = tbl, col.list = col.list, info = info)
  amcl <- col.list[amcl]
  amox <- col.list[amox]
  cfot <- col.list[cfot]
  cfta <- col.list[cfta]
  cfur <- col.list[cfur]
  cipr <- col.list[cipr]
  coli <- col.list[coli]
  eryt <- col.list[eryt]
  gent <- col.list[gent]
  mero <- col.list[mero]
  oxac <- col.list[oxac]
  pita <- col.list[pita]
  rifa <- col.list[rifa]
  teic <- col.list[teic]
  tetr <- col.list[tetr]
  tobr <- col.list[tobr]
  trsu <- col.list[trsu]
  vanc <- col.list[vanc]

  gram_positive = c(amox, amcl, cfur, pita, cipr, trsu,
                    # specific for G+:
                    vanc, teic, tetr, eryt, oxac, rifa)
  gram_positive <- gram_positive[!is.na(gram_positive)]

  gram_negative = c(amox, amcl, cfur, pita, cipr, trsu,
                    # specific for G-:
                    gent, tobr, coli, cfot, cfta, mero)
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
    gsub('(NA|NULL)', '-', .)

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
      # evaluation using regular expression will treat '?' as any character
      # so I is actually ignored then
      x <- gsub('I', '?', x, ignore.case = TRUE)
      y <- gsub('I', '?', y, ignore.case = TRUE)
    }
    for (i in 1:length(x)) {
      result[i] <- grepl(x = x[i],
                         pattern = y[i],
                         ignore.case = TRUE) |
        grepl(x = y[i],
              pattern = x[i],
              ignore.case = TRUE)
    }
    return(result)
  } else {

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

        if (type == 'points') {
          # count points for every single character:
          # - no change is 0 points
          # - I <-> S|R is 0.5 point
          # - S|R <-> R|S is 1 point
          # use the levels of as.rsi (S = 1, I = 2, R = 3)

          suppressWarnings(x2 <- x2 %>% as.rsi() %>% as.double())
          suppressWarnings(y2 <- y2 %>% as.rsi() %>% as.double())

          points <- (x2 - y2) %>% abs() %>% sum(na.rm = TRUE)
          result[i] <- ((points / 2) >= points_threshold)

        } else {
          stop('`', type, '` is not a valid value for type, must be "points" or "keyantibiotics". See ?first_isolate.')
        }
      }
    }
    if (info == TRUE) {
      cat('\n')
    }
    result
  }
}
