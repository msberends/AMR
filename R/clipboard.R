#' Import/export from clipboard
#'
#' These are helper functions around \code{\link{read.table}} and \code{\link{write.table}} to import from and export to clipboard. The data will be read and written as tab-separated by default, which makes it possible to copy and paste from other software like Excel and SPSS without further transformation.
#' @rdname clipboard
#' @name clipboard
#' @inheritParams utils::read.table
#' @inheritParams utils::write.table
#' @param startrow \emph{n}th row to start importing from. For \code{clipboard_import}, when \code{header = TRUE} the import will start on row \code{startrow} \emph{below} the header.
#' @param as_vector a logical value indicating whether data consisting of only one column should be imported as vector using \code{\link[dplyr]{pull}}. This will strip off the header.
#' @keywords clipboard clipboard_import clipboard_export import export
#' @importFrom dplyr %>% pull as_tibble
#' @importFrom utils read.delim write.table object.size writeClipboard
#' @details For \code{clipboard_export}, the reserved clipboard size for exporting will be set automatically to 125\% of the object size of \code{x}. This way, it is possible to export data with thousands of rows as the only limit will be your systems RAM.
#' @export
#' @return data.frame
clipboard_import <- function(sep = '\t',
                             header = TRUE,
                             dec = ".",
                             na = c("", "NA", "NULL"),
                             startrow = 1,
                             as_vector = TRUE) {
  
  import_tbl <- read.delim(file = 'clipboard',
                           sep = sep,
                           header = header,
                           strip.white = TRUE,
                           dec = dec,
                           na.strings = na,
                           fileEncoding = 'UTF-8',
                           encoding = 'UTF-8',
                           stringsAsFactors = FALSE)
  
  # use tibble, so column types will be translated correctly
  import_tbl <- as_tibble(import_tbl)
  
  if (startrow > 1) {
    # would else lose column headers
    import_tbl <- import_tbl[startrow:nrow(import_tbl),]
  }
  
  colnames(import_tbl) <- gsub('[.]+', '_', colnames(import_tbl))
  
  if (NCOL(import_tbl) == 1 & as_vector == TRUE) {
    import_tbl %>% pull(1)
  } else {
    import_tbl
  }
}

#' @rdname clipboard
#' @importFrom dplyr %>% pull as_tibble
#' @export
clipboard_export <- function(x,
                             sep = '\t',
                             dec = ".",
                             na = "",
                             header = TRUE) {
  
  x <- deparse(substitute(x))
  size <- x %>%
    get() %>% 
    object.size() %>%
    formatC(format = 'd') %>%
    as.integer()
  
  x <- get(x)
  
  if (size > 25 * 1024 * 1024) {
    # above 25 MB use a hacker function
    writeClipboard(knitr::kable(x))
  } else {
    # set size of clipboard to 125% of the object size of x
    write.table(x = x,
                file = paste0("clipboard-", size * 1.25),
                sep = sep,
                na = na,
                row.names = FALSE,
                col.names = header,
                dec = dec,
                quote = FALSE)
  }
  cat("Successfully exported to clipboard:", NROW(x), "obs. of", NCOL(x), "variables.\n")
  
}
