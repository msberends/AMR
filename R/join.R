#' Join a table with \code{microorganisms}
#'
#' Join the dataset \code{\link{microorganisms}} easily to an existing table or character vector.
#' @rdname join
#' @name join
#' @aliases join inner_join
#' @param x existing table to join, also supports character vectors
#' @param by a variable to join by - could be a column name of \code{x} with values that exist in \code{microorganisms$bactid} (like \code{by = "bacteria_id"}), or another column in \code{\link{microorganisms}} (but then it should be named, like \code{by = c("my_genus_species" = "fullname")})
#' @param suffix if there are non-joined duplicate variables in \code{x} and \code{y}, these suffixes will be added to the output to disambiguate them. Should be a character vector of length 2.
#' @param ... other parameters to pass on to \code{dplyr::\link[dplyr]{join}}.
#' @details As opposed to the \code{\link[dplyr]{join}} functions of \code{dplyr}, characters vectors are supported and at default existing columns will get a suffix \code{"2"} and the newly joined columns will not get a suffix. See \code{\link[dplyr]{join}} for more information.
#' @export
#' @examples 
#' left_join_microorganisms("STAAUR")
#' 
#' library(dplyr)
#' septic_patients %>% left_join_microorganisms()
#' 
#' df <- data.frame(date = seq(from = as.Date("2018-01-01"),
#'                             to = as.Date("2018-01-07"),
#'                             by = 1),
#'                  bacteria_id = c("STAAUR", "STAAUR", "STAAUR", "STAAUR",
#'                                  "ESCCOL", "ESCCOL", "ESCCOL"),
#'                  stringsAsFactors = FALSE)
#' colnames(df)
#' df2 <- left_join_microorganisms(df, "bacteria_id")
#' colnames(df2)
inner_join_microorganisms <- function(x, by = 'bactid', suffix = c("2", ""), ...) {
  if (any(class(x) %in% c('character', 'factor'))) {
    x <- data.frame(bactid = x, stringsAsFactors = FALSE)
  }
  # no name set to `by` parameter
  if (is.null(names(by))) {
    joinby <- colnames(AMR::microorganisms)[1]
    names(joinby) <- by
  } else {
    joinby <- by
  }
  join <- dplyr::inner_join(x = x, y = AMR::microorganisms, by = joinby, suffix = c("2", ""), ...)
  if (nrow(join) > nrow(x)) {
    warning('the newly joined tbl contains ', nrow(join) - nrow(x), ' rows more that its original')
  }
  join
}

#' @rdname join
#' @export
left_join_microorganisms <- function(x, by = 'bactid', suffix = c("2", ""), ...) {
  if (any(class(x) %in% c('character', 'factor'))) {
    x <- data.frame(bactid = x, stringsAsFactors = FALSE)
  }
  # no name set to `by` parameter
  if (is.null(names(by))) {
    joinby <- colnames(AMR::microorganisms)[1]
    names(joinby) <- by
  } else {
    joinby <- by
  }
  join <- dplyr::left_join(x = x, y = AMR::microorganisms, by = joinby, suffix = c("2", ""), ...)
  if (nrow(join) > nrow(x)) {
    warning('the newly joined tbl contains ', nrow(join) - nrow(x), ' rows more that its original')
  }
  join
}

#' @rdname join
#' @export
right_join_microorganisms <- function(x, by = 'bactid', suffix = c("2", ""), ...) {
  if (any(class(x) %in% c('character', 'factor'))) {
    x <- data.frame(bactid = x, stringsAsFactors = FALSE)
  }
  # no name set to `by` parameter
  if (is.null(names(by))) {
    joinby <- colnames(AMR::microorganisms)[1]
    names(joinby) <- by
  } else {
    joinby <- by
  }
  join <- dplyr::right_join(x = x, y = AMR::microorganisms, by = joinby, suffix = c("2", ""), ...)
  if (nrow(join) > nrow(x)) {
    warning('the newly joined tbl contains ', nrow(join) - nrow(x), ' rows more that its original')
  }
  join
}

#' @rdname join
#' @export
full_join_microorganisms <- function(x, by = 'bactid', suffix = c("2", ""), ...) {
  if (any(class(x) %in% c('character', 'factor'))) {
    x <- data.frame(bactid = x, stringsAsFactors = FALSE)
  }
  # no name set to `by` parameter
  if (is.null(names(by))) {
    joinby <- colnames(AMR::microorganisms)[1]
    names(joinby) <- by
  } else {
    joinby <- by
  }
  join <- dplyr::full_join(x = x, y = AMR::microorganisms, by = joinby, suffix = c("2", ""), ...)
  if (nrow(join) > nrow(x)) {
    warning('the newly joined tbl contains ', nrow(join) - nrow(x), ' rows more that its original')
  }
  join
}

#' @rdname join
#' @export
semi_join_microorganisms <- function(x, by = 'bactid', ...) {
  if (any(class(x) %in% c('character', 'factor'))) {
    x <- data.frame(bactid = x, stringsAsFactors = FALSE)
  }
  # no name set to `by` parameter
  if (is.null(names(by))) {
    joinby <- colnames(AMR::microorganisms)[1]
    names(joinby) <- by
  } else {
    joinby <- by
  }
  dplyr::semi_join(x = x, y = AMR::microorganisms, by = joinby, ...)
}

#' @rdname join
#' @export
anti_join_microorganisms <- function(x, by = 'bactid', ...) {
  if (any(class(x) %in% c('character', 'factor'))) {
    x <- data.frame(bactid = x, stringsAsFactors = FALSE)
  }
  # no name set to `by` parameter
  if (is.null(names(by))) {
    joinby <- colnames(AMR::microorganisms)[1]
    names(joinby) <- by
  } else {
    joinby <- by
  }
  dplyr::anti_join(x = x, y = AMR::microorganisms, by = joinby, ...)
}
