#' Import/export from clipboard
#'
#' @description These are helper functions around \code{\link{read.table}} and \code{\link{write.table}} to import from and export to clipboard with support for Windows, Linux and macOS.
#'
#' The data will be read and written as tab-separated by default, which makes it possible to copy and paste from other software like Excel and SPSS without further transformation.
#'
#' This also supports automatic column type transformation, with AMR classes \code{\link{as.rsi}} and \code{\link{as.mic}}.
#' @rdname clipboard
#' @name clipboard
#' @inheritParams base::data.frame
#' @inheritParams utils::read.table
#' @inheritParams utils::write.table
#' @inheritParams readr::locale
#' @param startrow \emph{n}th row to start importing from. When \code{header = TRUE}, the import will start on row \code{startrow} \emph{below} the header.
#' @param as_vector a logical value indicating whether data consisting of only one column should be imported as vector using \code{\link[dplyr]{pull}}. This will strip off the header.
#' @param guess_col_types a logical value indicating whether column types should be guessed and transformed automatically with \code{\link[readr]{parse_guess}} from the \code{readr} package. Besides, the antimicrobial classes in this AMR package (\code{\link{as.rsi}} and \code{\link{as.mic}}) are also supported.
#' @param remove_ASCII_escape_char remove ASCII escape character
#' @param info print info to console
#' @keywords clipboard clipboard_import clipboard_export import export
#' @importFrom dplyr %>% pull as_tibble
#' @importFrom clipr read_clip_tbl write_clip
#' @importFrom utils read.delim write.table object.size
#' @importFrom readr parse_guess locale
#' @details
#' The parameter \code{stringsAsFactors} defaults to \code{FALSE}, as opposed to most base \R methods.
#'
#' The parameters \code{date_format} and \code{time_format} also support generic date and time formats like \code{"dd-mm-yyyy"} like Excel.
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
#' @examples
#' \dontrun{
#'
#' df1 <- data.frame(a = letters[1:12],
#'                   b = runif(n = 12, min = 1000, max = 2000),
#'                   stringsAsFactors = FALSE)
#' clipboard_export(df1)
#' df2 <- clipboard_import()
#' identical(df1, df2)
#'
#' # send frequency table to clipboard (e.g. for pasting in Excel)
#' septic_patients %>%
#'   freq(age) %>%
#'   format() %>%       # this will format the percentages
#'   clipboard_export()
#' }
clipboard_import <- function(sep = '\t',
                             quote = "",
                             header = TRUE,
                             dec = ".",
                             na = c("", "NA", "NULL"),
                             stringsAsFactors = FALSE,
                             startrow = 1,
                             as_vector = TRUE,
                             guess_col_types = TRUE,
                             date_names = 'en',
                             date_format = '%Y-%m-%d',
                             time_format = '%H:%M',
                             remove_ASCII_escape_char = FALSE,
                             tz = Sys.timezone(),
                             encoding = "UTF-8",
                             info = TRUE) {

  if (!clipr::clipr_available() & Sys.info()['sysname'] == "Linux") {
    # try to support on X11, by setting the R variable DISPLAY
    Sys.setenv(DISPLAY = "localhost:10.0")
  }

  # this will fail when clipr is (still) not available
  import_tbl <- clipr::read_clip_tbl(file = file,
                                     sep = sep,
                                     quote = quote,
                                     header = header,
                                     strip.white = TRUE,
                                     dec = dec,
                                     na.strings = na,
                                     encoding = 'UTF-8',
                                     stringsAsFactors = stringsAsFactors)

  # use tibble, so column types will be translated correctly
  import_tbl <- as_tibble(import_tbl)

  if (startrow > 1) {
    # would else lose column headers
    import_tbl <- import_tbl[startrow:NROW(import_tbl),]
  }

  colnames(import_tbl) <- gsub('[.]+', '_', colnames(import_tbl))

  if (guess_col_types == TRUE) {
    if (info == TRUE) {
      cat('Transforming data by guessing column types...')
    }
    import_tbl <- tbl_parse_guess(tbl = import_tbl,
                                  date_names = date_names,
                                  date_format = date_format,
                                  time_format = time_format,
                                  decimal_mark = dec,
                                  tz = tz,
                                  encoding = encoding,
                                  remove_ASCII_escape_char = remove_ASCII_escape_char,
                                  na = na)
    if (info == TRUE) {
      cat('OK\n')
    }
  }

  if (NCOL(import_tbl) == 1 & as_vector == TRUE) {
    import_tbl <-  import_tbl %>% pull(1)
  }

  # and transform back to data.frame
  import_tbl <- as.data.frame(import_tbl, stringsAsFactors = stringsAsFactors)

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

  if (!clipr::clipr_available() & Sys.info()['sysname'] == "Linux") {
    # try to support on X11, by setting the R variable DISPLAY
    Sys.setenv(DISPLAY = "localhost:10.0")
  }

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

