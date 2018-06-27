#' Import/export from clipboard
#'
#' @description These are helper functions around \code{\link{read.table}} and \code{\link{write.table}} to import from and export to clipboard with support for Windows, Linux and macOS.
#'
#' The data will be read and written as tab-separated by default, which makes it possible to copy and paste from other software like Excel and SPSS without further transformation.
#'
#' Supports automatic column type transformation and supports new classes \code{\link{as.rsi}} and \code{\link{as.mic}}.
#' @rdname clipboard
#' @name clipboard
#' @inheritParams utils::read.table
#' @inheritParams utils::write.table
#' @inheritParams readr::locale
#' @param startrow \emph{n}th row to start importing from. When \code{header = TRUE}, the import will start on row \code{startrow} \emph{below} the header.
#' @param as_vector a logical value indicating whether data consisting of only one column should be imported as vector using \code{\link[dplyr]{pull}}. This will strip off the header.
#' @param guess_col_types a logical value indicating whether column types should be guessed with the \code{readr} package.
#' @param info print info about copying
#' @keywords clipboard clipboard_import clipboard_export import export
#' @importFrom dplyr %>% pull as_tibble
#' @importFrom clipr read_clip_tbl write_clip
#' @importFrom utils read.delim write.table object.size
#' @importFrom readr parse_guess locale
#' @details
#' When using \code{guess_col_types = TRUE}, all column types will be determined automatically with \code{\link[readr]{parse_guess}} from the \code{readr} package. Besides, the antimicrobial classes in this AMR package (\code{\link{as.rsi}} and \code{\link{as.mic}}) are also supported.
#'
#'   \if{html}{
#'     \strong{Example for copying from Excel:}
#'     \out{<div style="text-align: left">}\figure{clipboard_copy.png}\out{</div>}
#'     \cr
#'     \strong{And pasting in R:} \cr
#'     \cr
#'     \code{> data <- clipboard_import()} \cr
#'     \code{> data} \cr
#'     \out{<div style="text-align: left">}\figure{clipboard_paste.png}\out{</div>}
#'     \cr
#'     \strong{The resulting data contains the right RSI-classes:} \cr
#'     \cr
#'     \code{> data$amox} \cr
#'     \out{<div style="text-align: left">}\figure{clipboard_rsi.png}\out{</div>}
#'   }
#' @export
clipboard_import <- function(sep = '\t',
                             header = TRUE,
                             dec = ".",
                             na = c("", "NA", "NULL"),
                             startrow = 1,
                             as_vector = TRUE,
                             guess_col_types = TRUE,
                             date_names = 'en',
                             date_format = '%Y-%m-%d',
                             time_format = '%H:%M',
                             tz = Sys.timezone(),
                             encoding = "UTF-8",
                             info = FALSE) {

  if (clipr::clipr_available() & info == TRUE) {
    cat('Importing from clipboard...')
  }

  # this will fail when clipr is not available
  import_tbl <- clipr::read_clip_tbl(file = file,
                                     sep = sep,
                                     header = header,
                                     strip.white = TRUE,
                                     dec = dec,
                                     na.strings = na,
                                     encoding = 'UTF-8',
                                     stringsAsFactors = FALSE)

  if (info == TRUE) {
    cat('OK\n')
  }

  # use tibble, so column types will be translated correctly
  import_tbl <- as_tibble(import_tbl)

  if (startrow > 1) {
    # would else lose column headers
    import_tbl <- import_tbl[startrow:nrow(import_tbl),]
  }

  colnames(import_tbl) <- gsub('[.]+', '_', colnames(import_tbl))

  if (guess_col_types == TRUE) {
    if (info == TRUE) {
      cat('Transforming columns with readr::parse_guess...')
    }
    import_tbl <- clipboard_format(tbl = import_tbl,
                                   date_names = date_names,
                                   date_format = date_format,
                                   time_format = time_format,
                                   decimal_mark = dec,
                                   tz = tz,
                                   encoding = encoding,
                                   na = na)
    if (info == TRUE) {
      cat('OK\n')
    }
  }

  if (NCOL(import_tbl) == 1 & as_vector == TRUE) {
    import_tbl <-  import_tbl %>% pull(1)
  }

  if (info == TRUE) {
    cat("Successfully imported from clipboard:", NROW(import_tbl), "obs. of", NCOL(import_tbl), "variables.\n")
  }

  import_tbl

}

#' @rdname clipboard
#' @importFrom dplyr %>% pull as_tibble
#' @export
clipboard_export <- function(x,
                             sep = '\t',
                             dec = ".",
                             na = "",
                             header = TRUE,
                             info = TRUE) {

  clipr::write_clip(content = x,
                    na = na,
                    sep = sep,
                    row.names = FALSE,
                    col.names = header,
                    dec = dec,
                    quote = FALSE)

  if (info == TRUE) {
    cat("Successfully exported to clipboard:", NROW(x), "obs. of", NCOL(x), "variables.\n")
  }

}

clipboard_format <- function(tbl,
                             date_names = 'en',
                             date_format = '%Y-%m-%d',
                             time_format = '%H:%M',
                             decimal_mark = '.',
                             tz = Sys.timezone(),
                             encoding = "UTF-8",
                             na = c("", "NA", "NULL")) {

  date_format <- date_generic(date_format)
  time_format <- date_generic(time_format)
  # set col types with readr
  for (i in 1:ncol(tbl)) {
    if (!all(tbl %>% pull(i) %>% class() %in% c('list', 'matrix'))) {
      tbl[, i] <- readr::parse_guess(x = tbl %>% pull(i) %>% as.character(),
                                     na = na,
                                     locale = readr::locale(date_names = date_names,
                                                            date_format = date_format,
                                                            time_format = time_format,
                                                            decimal_mark = decimal_mark,
                                                            encoding = encoding,
                                                            tz = tz,
                                                            asciify = FALSE))
    }
    if (any(tbl %>% pull(i) %>% class() %in% c('factor', 'character'))) {
      # get values
      distinct_val <- tbl %>% pull(i) %>% unique() %>% sort()
      # remove ASCII escape character: https://en.wikipedia.org/wiki/Escape_character#ASCII_escape_character
      tbl[, i] <- tbl %>% pull(i) %>% gsub('\033', ' ', ., fixed = TRUE)
      # look for RSI, shouldn't all be "" and must be valid antibiotic interpretations
      if (!all(distinct_val[!is.na(distinct_val)] == '')
          & all(distinct_val[!is.na(distinct_val)] %in% c('', 'I', 'I;I', 'R', 'R;R', 'S', 'S;S'))) {
        tbl[, i] <- tbl %>% pull(i) %>% as.rsi()
      }
    }
    # convert to MIC class
    if (colnames(tbl)[i] %like% '_mic$') {
      tbl[, i] <- tbl %>% pull(i) %>% as.mic()
    }
  }
  tbl
}

