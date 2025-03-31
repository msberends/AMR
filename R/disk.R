# ==================================================================== #
# TITLE:                                                               #
# AMR: An R Package for Working with Antimicrobial Resistance Data     #
#                                                                      #
# SOURCE CODE:                                                         #
# https://github.com/msberends/AMR                                     #
#                                                                      #
# PLEASE CITE THIS SOFTWARE AS:                                        #
# Berends MS, Luz CF, Friedrich AW, et al. (2022).                     #
# AMR: An R Package for Working with Antimicrobial Resistance Data.    #
# Journal of Statistical Software, 104(3), 1-31.                       #
# https://doi.org/10.18637/jss.v104.i03                                #
#                                                                      #
# Developed at the University of Groningen and the University Medical  #
# Center Groningen in The Netherlands, in collaboration with many      #
# colleagues from around the world, see our website.                   #
#                                                                      #
# This R package is free software; you can freely use and distribute   #
# it for both personal and commercial purposes under the terms of the  #
# GNU General Public License version 2.0 (GNU GPL-2), as published by  #
# the Free Software Foundation.                                        #
# We created this package for both routine data analysis and academic  #
# research and it was publicly released in the hope that it will be    #
# useful, but it comes WITHOUT ANY WARRANTY OR LIABILITY.              #
#                                                                      #
# Visit our website for the full manual and a complete tutorial about  #
# how to conduct AMR data analysis: https://msberends.github.io/AMR/   #
# ==================================================================== #

#' Transform Input to Disk Diffusion Diameters
#'
#' This transforms a vector to a new class [`disk`], which is a disk diffusion growth zone size (around an antibiotic disk) in millimetres between 0 and 50.
#' @rdname as.disk
#' @param x Vector
#' @param na.rm A [logical] indicating whether missing values should be removed
#' @details Interpret disk values as SIR values with [as.sir()]. It supports guidelines from EUCAST and CLSI.
#'
#' Disk diffusion growth zone sizes must be between 0 and 50 millimetres. Values higher than 50 but lower than 100 will be maximised to 50. All others input values outside the 0-50 range will return `NA`.
#' @return An [integer] with additional class [`disk`]
#' @aliases disk
#' @export
#' @seealso [as.sir()]
#' @examples
#' # transform existing disk zones to the `disk` class (using base R)
#' df <- data.frame(
#'   microorganism = "Escherichia coli",
#'   AMP = 20,
#'   CIP = 14,
#'   GEN = 18,
#'   TOB = 16
#' )
#' df[, 2:5] <- lapply(df[, 2:5], as.disk)
#' str(df)
#'
#' \donttest{
#' # transforming is easier with dplyr:
#' if (require("dplyr")) {
#'   df %>% mutate(across(AMP:TOB, as.disk))
#' }
#' }
#'
#' # interpret disk values, see ?as.sir
#' as.sir(
#'   x = as.disk(18),
#'   mo = "Strep pneu", # `mo` will be coerced with as.mo()
#'   ab = "ampicillin", # and `ab` with as.ab()
#'   guideline = "EUCAST"
#' )
#'
#' # interpret whole data set, pretend to be all from urinary tract infections:
#' as.sir(df, uti = TRUE)
as.disk <- function(x, na.rm = FALSE) {
  meet_criteria(x, allow_NA = TRUE)
  meet_criteria(na.rm, allow_class = "logical", has_length = 1)

  if (!is.disk(x)) {
    x <- unlist(x)
    if (isTRUE(na.rm)) {
      x <- x[!is.na(x)]
    }
    x[trimws2(x) == ""] <- NA
    x.bak <- x

    na_before <- length(x[is.na(x)])

    # heavily based on cleaner::clean_double():
    clean_double2 <- function(x, remove = "[^0-9.,-]", fixed = FALSE) {
      x <- gsub(",", ".", x, fixed = TRUE)
      # remove ending dot/comma
      x <- gsub("[,.]$", "", x)
      # only keep last dot/comma
      reverse <- function(x) vapply(FUN.VALUE = character(1), lapply(strsplit(x, NULL), rev), paste, collapse = "")
      x <- sub("{{dot}}", ".",
        gsub(".", "",
          reverse(sub(".", "}}tod{{",
            reverse(x),
            fixed = TRUE
          )),
          fixed = TRUE
        ),
        fixed = TRUE
      )
      x_clean <- gsub(remove, "", x, ignore.case = TRUE, fixed = fixed)
      # remove everything that is not a number or dot
      as.double(gsub("[^0-9.]+", "", x_clean))
    }

    # round up and make it an integer
    x <- as.integer(ceiling(clean_double2(x)))

    # disks can never be less than 0 mm or more than 50 mm
    x[x < 0 | x > 99] <- NA_integer_
    x[x > 50] <- 50L
    na_after <- length(x[is.na(x)])

    if (na_before != na_after) {
      list_missing <- x.bak[is.na(x) & !is.na(x.bak)] %pm>%
        unique() %pm>%
        sort() %pm>%
        vector_and(quotes = TRUE)
      cur_col <- get_current_column()
      warning_("in `as.disk()`: ", na_after - na_before, " result",
        ifelse(na_after - na_before > 1, "s", ""),
        ifelse(is.null(cur_col), "", paste0(" in index '", cur_col, "'")),
        " truncated (",
        round(((na_after - na_before) / length(x)) * 100),
        "%) that were invalid disk zones: ",
        list_missing,
        call = FALSE
      )
    }
  }
  set_clean_class(as.integer(x),
    new_class = c("disk", "integer")
  )
}

all_valid_disks <- function(x) {
  if (!inherits(x, c("disk", "character", "numeric", "integer"))) {
    return(FALSE)
  }
  x_disk <- tryCatch(suppressWarnings(as.disk(x[!is.na(x)])),
    error = function(e) NA
  )
  !anyNA(x_disk) && !all(is.na(x))
}

#' @rdname as.disk
#' @details `NA_disk_` is a missing value of the new `disk` class.
#' @export
NA_disk_ <- set_clean_class(as.integer(NA_real_),
  new_class = c("disk", "integer")
)

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

#' @method print disk
#' @export
#' @noRd
print.disk <- function(x, ...) {
  cat("Class 'disk'\n")
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
c.disk <- function(...) {
  as.disk(unlist(lapply(list(...), as.character)))
}

#' @method unique disk
#' @export
#' @noRd
unique.disk <- function(x, incomparables = FALSE, ...) {
  y <- NextMethod()
  attributes(y) <- attributes(x)
  y
}

#' @method rep disk
#' @export
#' @noRd
rep.disk <- function(x, ...) {
  y <- NextMethod()
  attributes(y) <- attributes(x)
  y
}

# will be exported using s3_register() in R/zzz.R
get_skimmers.disk <- function(column) {
  skimr::sfl(
    skim_type = "disk",
    min = ~ min(as.double(.), na.rm = TRUE),
    max = ~ max(as.double(.), na.rm = TRUE),
    median = ~ stats::median(as.double(.), na.rm = TRUE),
    n_unique = ~ length(unique(stats::na.omit(.))),
    hist = ~ skimr::inline_hist(stats::na.omit(as.double(.)))
  )
}
