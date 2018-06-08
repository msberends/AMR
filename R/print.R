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

#' Printing Data Tables and Tibbles
#'
#' Print a data table or tibble. It prints: \cr- The \strong{first and last rows} like \code{data.table}s are printed by the \code{data.table} package,\cr- A \strong{header} and \strong{left aligned text} like \code{tibble}s are printed by the \code{tibble} package with info about grouped variables,\cr- \strong{Unchanged values} and \strong{support for row names} like \code{data.frame}s are printed by the \code{base} package.
#' @inheritParams base::print.data.frame
#' @param nmax amount of rows to print in total. When the total amount of rows exceeds this limit, the first and last \code{nmax / 2} rows will be printed. Use \code{nmax = NA} to print all rows.
#' @param header print header with information about data size and tibble grouping
#' @param print.keys print keys for \code{data.table}
#' @param na value to print instead of NA
#' @param width amount of white spaces to keep between columns, must be at least 1
#' @rdname print
#' @name print
#' @importFrom dplyr %>% n_groups group_vars group_size filter pull select
#' @importFrom data.table data.table
#' @importFrom utils object.size
#' @exportMethod print.tbl_df
#' @export
#' @examples
#' # more reliable data view:
#' library(dplyr)
#' starwars
#' print(starwars, width = 3)
#'
#' # This is how the tibble package prints since v1.4.0:
#' # (mind the quite unfamiliar underscores and ending dots)
#' tibble(now_what = c(1.2345, 2345.67, 321.456)) %>% tibble:::print.tbl_df()
#'
#' # This is how this AMR package prints:
#' # (every number shown as you would expect)
#' tibble(now_what = c(1.2345, 2345.67, 321.456))
#'
#' # also supports info about groups (look at header)
#' starwars %>% group_by(homeworld, gender)
print.tbl_df <- function(x,
                         nmax = 10,
                         header = TRUE,
                         row.names = TRUE,
                         right = FALSE,
                         width = 1,
                         na = "<NA>",
                         ...) {
  prettyprint_df(x = x,
                 nmax = nmax,
                 header = header,
                 row.names = row.names,
                 print.keys = FALSE,
                 right = right,
                 width = width,
                 na = na,
                 ...)
}

#' @rdname print
#' @exportMethod print.tbl
#' @export
print.tbl <- function(x, ...) {
  prettyprint_df(x, ...)
}

#' @rdname print
#' @exportMethod print.data.table
#' @export
print.data.table <- function(x,
                             print.keys = FALSE,
                             ...) {
  prettyprint_df(x = x,
                 print.keys = print.keys,
                 ...)
}

printDT <- data.table:::print.data.table
prettyprint_df <- function(x,
                           nmax = 10,
                           header = TRUE,
                           row.names = TRUE,
                           print.keys = FALSE,
                           right = FALSE,
                           width = 1,
                           na = "<NA>",
                           ...) {

  ansi_reset <- "\u001B[0m"
  ansi_black <- "\u001B[30m"
  ansi_red <- "\u001B[31m"
  ansi_green <- "\u001B[32m"
  ansi_yellow <- "\u001B[33m"
  ansi_blue <- "\u001B[34m"
  ansi_purple <- "\u001B[35m"
  ansi_cyan <- "\u001B[36m"
  ansi_white <- "\u001B[37m"
  ansi_gray <- "\u001B[38;5;246m"

  if (width < 1) {
    stop('`width` must be at least 1.', call. = FALSE)
  }

  if (is.na(nmax)) {
    nmax <- NROW(x)
  }
  n <- nmax
  if (n %% 2 == 1) {
    # odd number; add 1
    n <- n + 1
  }

  width <- width - 1

  if ('tbl_df' %in% class(x)) {
    type <- 'tibble'
  } else if ('data.table' %in% class(x)) {
    type <- 'data.table'
  } else {
    type <- 'data.frame'
  }

  if (header == TRUE) {
    if (NCOL(x) == 1) {
      vars <- 'variable'
    } else {
      vars <- 'variables'
    }

    size <- object.size(x) %>% as.double() %>% size_humanreadable()

    cat(paste0("A ", type,": ",
               format(NROW(x)),
               " obs. of ",
               format(NCOL(x)),
               " ", vars,
               ansi_gray, " (", size, ")\n", ansi_reset))
    if ('grouped_df' %in% class(x) & n_groups(x) > 0) {
      cat(paste0("Grouped by ",
                 x %>% group_vars() %>% paste0(ansi_red, ., ansi_reset) %>% paste0(collapse = " and "),
                 ansi_gray,
                 " (",
                 x %>% n_groups(),
                 " groups with sizes between ",
                 x %>% group_size() %>% min(),
                 " and ",
                 x %>% group_size() %>% max(),
                 ")\n",
                 ansi_reset))
    }
    if (!is.null(attributes(x)$qry)) {
      cat(paste0(ansi_gray, "This data contains a query. Use qry() to view it.\n", ansi_reset))
    }
    cat("\n")
  }

  # data.table where keys should be printed
  if (print.keys == TRUE) {
    printDT(x,
            class = header,
            row.names = row.names,
            print.keys = TRUE,
            right = right,
            ...
    )
    return(invisible())
  }

  # tibbles give warning when setting column names
  x <- x %>% base::as.data.frame(stringsAsFactors = FALSE)

  # extra space of 3 chars, right to row name or number
  if (NROW(x) > 0) {
    maxrowchars <- rownames(x) %>% nchar() %>% max() + 3
    rownames(x) <- paste0(rownames(x), strrep(" ", maxrowchars - nchar(rownames(x))))
  } else {
    maxrowchars <- 0
  }

  x.bak <- x

  if (n + 1 < nrow(x)) {
    # remove in between part, 1 extra for ~~~~ between first and last part
    rows_list <- c(1:(n / 2 + 1), (nrow(x) - (n / 2) + 1):nrow(x))
    x <- as.data.frame(x.bak[rows_list,], stringsAsFactors = FALSE)
    colnames(x) <- colnames(x.bak)
    rownames(x) <- rownames(x.bak)[rows_list]
    # set inbetweener between parts
    rownames(x)[n / 2 + 1] <- strrep("~", maxrowchars)
  }

  if (header == TRUE) {
    # add 1 row for classes
    # class will be marked up per column
    if (NROW(x.bak) > 0) {
      rownames.x <- rownames(x)
      # suppress warnings because dplyr want us to use library(dplyr) when using filter(row_number())
      suppressWarnings(
        x <- x %>%
          filter(row_number() == 1) %>%
          rbind(x, stringsAsFactors = FALSE)
      )
      rownames(x) <- c('*', rownames.x)
    }

    # select 1st class per column and abbreviate
    classes <- x.bak %>%
      sapply(class) %>%
      lapply(
        function(c) {
          # do print all POSIX classes like "POSct/t"
          if ('POSIXct' %in% c) {
            paste0('POS',
                   c %>%
                     gsub('POSIX', '', .) %>%
                     paste0(collapse = '/'))
          } else {
            if (NCOL(.) > 1) {
              .[1,]
            } else {
              c[[1]]
            }
          }
        }) %>%
      unlist() %>%
      gsub("character",  "chr", ., fixed = TRUE) %>%
      gsub("complex",    "cplx", ., fixed = TRUE) %>%
      gsub("Date",       "Date", ., fixed = TRUE) %>%
      gsub("double",     "dbl", ., fixed = TRUE) %>%
      gsub("expression", "expr", ., fixed = TRUE) %>%
      gsub("factor",     "fct", ., fixed = TRUE) %>%
      gsub("IDate",      "IDat", ., fixed = TRUE) %>%
      gsub("integer",    "int", ., fixed = TRUE) %>%
      gsub("integer64",  "i64", ., fixed = TRUE) %>%
      gsub("list",       "list", ., fixed = TRUE) %>%
      gsub("logical",    "lgl", ., fixed = TRUE) %>%
      gsub("numeric",    "dbl", ., fixed = TRUE) %>%
      gsub("ordered",    "ord", ., fixed = TRUE) %>%
      gsub("percent",    "pct", ., fixed = TRUE) %>%
      gsub("single",     "sgl", ., fixed = TRUE) %>%
      paste0("<", ., ">")
  }

  # markup cols

  for (i in 1:ncol(x)) {
    if (all(!class(x[, i]) %in% class(x.bak[, i]))) {
      class(x[, i]) <- class(x.bak[, i])
    }
    try(x[, i] <- format(x %>% pull(i)), silent = TRUE)
    # replace NAs
    if (nchar(na) < 2) {
      # make as long as the text "NA"
      na <- paste0(na, strrep(" ", 2 - nchar(na)))
    }
    try(x[, i] <- gsub("^NA$", na, trimws(x[, i], 'both')), silent = TRUE)
    # place class into 1st row
    if (header == TRUE) {
      x[1, i] <- classes[i]
    }
    # dashes between two parts when exceeding nmax
    maxvalchars <- max(colnames(x)[i] %>% nchar(), x[, i] %>% nchar() %>% max())
    if (n + 1 < nrow(x.bak)) {
      x[n / 2 + if_else(header == TRUE, 2, 1), i] <- strrep("~", maxvalchars)
    }

    # align according to `right` parameter, but only factors, logicals text, but not MICs
    if (any(x.bak %>% pull(i) %>% class() %in% c('factor', 'character', 'logical'))
        & !("mic" %in% (x.bak %>% pull(i) %>% class()))) {
      vals <- x %>% pull(i) %>% trimws('both')
      colname <- colnames(x)[i] %>% trimws('both')
      if (right == FALSE) {
        vals <- paste0(vals, strrep(" ", maxvalchars - nchar(vals)))
        colname <- paste0(colname, strrep(" ", maxvalchars - nchar(colname)))
      } else {
        vals <- paste0(strrep(" ", maxvalchars - nchar(vals)), vals)
        colname <- paste0(strrep(" ", maxvalchars - nchar(colname)), colname)
      }
      x[, i] <- vals
      colnames(x)[i] <- colname
    }

    # add left padding according to `width` parameter
    # but not in 1st col when row names are off
    if (row.names == TRUE | i > 1) {
      x[, i] <- paste0(strrep(" ", width), x[, i])
      colnames(x)[i] <- paste0(strrep(" ", width), colnames(x)[i])
    }

    # strip columns that do not fit (width + 2 extra chars as margin)
    width_console <- options()$width
    width_until_col <- x %>%
      select(1:i) %>%
      apply(1, paste, collapse = strrep(" ", width + 2)) %>%
      nchar() %>%
      max()
    width_until_col_before <- x %>%
      select(1:(max(i, 2) - 1)) %>%
      apply(1, paste, collapse = strrep(" ", width + 2)) %>%
      nchar() %>%
      max()
    extraspace <- maxrowchars + nchar(rownames(x)[length(rownames(x))])
    width_until_colnames <- colnames(x)[1:i] %>% paste0(collapse = strrep(" ", width + 1)) %>% nchar() + extraspace
    width_until_colnames_before <- colnames(x)[1:(max(i, 2) - 1)] %>% paste0(collapse = strrep(" ", width + 1)) %>% nchar() + extraspace

    if (i > 1 &
        (width_until_col > width_console
         | width_until_colnames > width_console)) {
      if (width_until_col_before > width_console
          | width_until_colnames_before > width_console) {
        x <- x[, 1:(i - 2)]
      } else {
        x <- x[, 1:(i - 1)]
      }
      break
    }
  }

  # empty table, row name of header should be "*"
  if (NROW(x.bak) == 0) {
    rownames(x) <- '* '
  }

  # and here it is...
  suppressWarnings(
    base::print.data.frame(x, row.names = row.names, ...)
  )

  # print rest of col names when they were stripped
  if (ncol(x) < ncol(x.bak)) {
    x.notshown <- x.bak %>% select((ncol(x) + 1):ncol(x.bak))
    if (ncol(x.notshown) == 1) {
      cat('... and 1 more column: ')
    } else {
      cat('... and', ncol(x.notshown), 'more columns: ')
    }
    cat(x.notshown %>%
          colnames() %>%
          paste0(' ', ansi_gray, classes[(ncol(x) + 1):ncol(x.bak)], ansi_reset) %>%
          paste0(collapse = ", "), '\n')
  }
}
