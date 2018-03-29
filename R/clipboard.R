#' Import/export from clipboard
#'
#' These are helper functions around \code{\link{read.table}} and \code{\link{write.table}} to import from and export to clipboard, with support for Windows, Linux and macOS. The data will be read and written as tab-separated by default, which makes it possible to copy and paste from other software like Excel and SPSS without further transformation.
#' @rdname clipboard
#' @name clipboard
#' @inheritParams utils::read.table
#' @inheritParams utils::write.table
#' @param startrow \emph{n}th row to start importing from. For \code{clipboard_import}, when \code{header = TRUE} the import will start on row \code{startrow} \emph{below} the header.
#' @param as_vector a logical value indicating whether data consisting of only one column should be imported as vector using \code{\link[dplyr]{pull}}. This will strip off the header.
#' @param info print info about copying
#' @keywords clipboard clipboard_import clipboard_export import export
#' @importFrom dplyr %>% pull as_tibble
#' @importFrom utils read.delim write.table object.size
#' @details For \code{clipboard_export}, the reserved clipboard size for exporting will be set automatically to 125\% of the object size of \code{x}. This way, it is possible to export data with thousands of rows as the only limit will be your systems RAM.
#' @export
#' @return data.frame
clipboard_import <- function(sep = '\t',
                             header = TRUE,
                             dec = ".",
                             na = c("", "NA", "NULL"),
                             startrow = 1,
                             as_vector = TRUE) {

  if (is_Windows() == TRUE) {
    file <- 'clipboard'
  } else {
    # use xclip package
    check_xclip()
    file <- pipe("xclip -o -selection c", "r")
    on.exit(close(file))
  }
  
  import_tbl <- read.delim(file = file,
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
                             header = TRUE,
                             info = TRUE) {
  
  x <- deparse(substitute(x))
  size <- x %>%
    get() %>% 
    object.size() %>%
    formatC(format = 'd') %>%
    as.integer()
  
  x <- get(x)

  if (is_Windows() == TRUE) {
    # set size of clipboard to 125% of the object size of x
    file <- paste0("clipboard-", size * 1.25)
  } else {
    # use xclip package
    check_xclip()
    file <- pipe("xclip -i", "w")
    on.exit(close(file))
  }

  write.table(x = x,
              file = file,
              sep = sep,
              na = na,
              row.names = FALSE,
              col.names = header,
              dec = dec,
              quote = FALSE)

  if (info == TRUE) {
    cat("Successfully exported to clipboard:", NROW(x), "obs. of", NCOL(x), "variables.\n")
  }
}

is_Windows <- function() {
  Sys.info()['sysname'] %like% "Windows"
}
check_xclip <- function() {
  if (!isTRUE(file.exists(Sys.which("xclip")[1L]))) {
      stop("Please install Linux package xclip first.")
  }
}
