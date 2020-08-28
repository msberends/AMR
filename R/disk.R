# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Analysis                              #
#                                                                      #
# SOURCE                                                               #
# https://github.com/msberends/AMR                                     #
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
# Visit our website for more info: https://msberends.github.io/AMR.    #
# ==================================================================== #

#' Class 'disk'
#'
#' This transforms a vector to a new class [`disk`], which is a growth zone size (around an antibiotic disk) in millimetres between 6 and 50.
#' @inheritSection lifecycle Stable lifecycle
#' @rdname as.disk
#' @param x vector
#' @param na.rm a logical indicating whether missing values should be removed
#' @details Interpret disk values as RSI values with [as.rsi()]. It supports guidelines from EUCAST and CLSI.
#' @return An [`integer`] with additional new class [`disk`]
#' @aliases disk
#' @export
#' @seealso [as.rsi()]
#' @inheritSection AMR Read more on our website!
#' @examples
#' \dontrun{
#' # transform existing disk zones to the `disk` class
#' library(dplyr)
#' df <- data.frame(microorganism = "E. coli",
#'                  AMP = 20,
#'                  CIP = 14,
#'                  GEN = 18,
#'                  TOB = 16)
#' df <- df %>% mutate_at(vars(AMP:TOB), as.disk)
#' df
#' 
#' # interpret disk values, see ?as.rsi
#' as.rsi(x = as.disk(18),
#'        mo = "Strep pneu",  # `mo` will be coerced with as.mo()
#'        ab = "ampicillin",  # and `ab` with as.ab()
#'        guideline = "EUCAST")
#'        
#' as.rsi(df)
#' }
as.disk <- function(x, na.rm = FALSE) {
  if (!is.disk(x)) {
    x <- x %>% unlist()
    if (na.rm == TRUE) {
      x <- x[!is.na(x)]
    }
    x.bak <- x
    
    na_before <- length(x[is.na(x)])
    
    # heavily based on the function from our cleaner package:
    clean_double2 <- function(x, remove = "[^0-9.,-]", fixed = FALSE) {
      x <- gsub(",", ".", x)
      # remove ending dot/comma
      x <- gsub("[,.]$", "", x)
      # only keep last dot/comma
      reverse <- function(x) sapply(lapply(strsplit(x, NULL), rev), paste, collapse = "")
      x <- sub("{{dot}}", ".", 
               gsub(".", "",
                    reverse(sub(".", "}}tod{{",
                                reverse(x), 
                                fixed = TRUE)),
                    fixed = TRUE), 
               fixed = TRUE)
      x_clean <- gsub(remove, "", x, ignore.case = TRUE, fixed = fixed)
      # remove everything that is not a number or dot
      as.numeric(gsub("[^0-9.]+", "", x_clean))
    }
    
    # round up and make it an integer
    x <- as.integer(ceiling(clean_double2(x)))
    
    # disks can never be less than 6 mm (size of smallest disk) or more than 50 mm
    x[x < 6 | x > 50] <- NA_integer_
    na_after <- length(x[is.na(x)])
    
    if (na_before != na_after) {
      list_missing <- x.bak[is.na(x) & !is.na(x.bak)] %>%
        unique() %>%
        sort()
      list_missing <- paste0('"', list_missing, '"', collapse = ", ")
      warning(na_after - na_before, " results truncated (",
              round(((na_after - na_before) / length(x)) * 100),
              "%) that were invalid disk zones: ",
              list_missing, call. = FALSE)
    }
  }
  structure(as.integer(x),
            class = c("disk", "integer"))
}

all_valid_disks <- function(x) {
  x_disk <- suppressWarnings(as.disk(x[!is.na(x)]))
  !any(is.na(x_disk)) & !all(is.na(x))
}

#' @rdname as.disk
#' @export
is.disk <- function(x) {
  inherits(x, "disk")
}

# will be exported using s3_register() in R/zzz.R
pillar_shaft.disk <- function(x, ...) {
  out <- trimws(format(x))
  out[is.na(x)] <- font_na(NA)
  create_pillar_column(out, align = "right", width = 2)
}

# will be exported using s3_register() in R/zzz.R
type_sum.disk <- function(x, ...) {
  "disk"
}

#' @method print disk
#' @export
#' @noRd
print.disk <- function(x, ...) {
  cat("Class <disk>\n")
  print(as.integer(x), quote = FALSE)
}

#' @method [ disk
#' @export
#' @noRd
"[.disk" <- function(x, ...) {
  y <- NextMethod()
  attributes(y) <- attributes(x)
  y
}
#' @method [[ disk
#' @export
#' @noRd
"[[.disk" <- function(x, ...) {
  y <- NextMethod()
  attributes(y) <- attributes(x)
  y
}
#' @method [<- disk
#' @export
#' @noRd
"[<-.disk" <- function(i, j, ..., value) {
  value <- as.disk(value)
  y <- NextMethod()
  attributes(y) <- attributes(i)
  y
}
#' @method [[<- disk
#' @export
#' @noRd
"[[<-.disk" <- function(i, j, ..., value) {
  value <- as.disk(value)
  y <- NextMethod()
  attributes(y) <- attributes(i)
  y
}
#' @method c disk
#' @export
#' @noRd
c.disk <- function(x, ...) {
  y <- NextMethod()
  y <- as.disk(y)
  attributes(y) <- attributes(x)
  y
}
