#' Join a table with \code{bactlist}
#'
#' Join the list of microorganisms \code{\link{bactlist}} easily to an existing table.
#' @rdname join
#' @name join
#' @aliases join inner_join
#' @param x existing table to join
#' @param by a variable to join by - could be a column name of \code{x} with values that exist in \code{bactlist$bactid} (like \code{by = "bacteria_id"}), or another column in \code{\link{bactlist}} (but then it should be named, like \code{by = c("my_genus_species" = "fullname")})
#' @param ... other parameters to pass on to \code{dplyr::\link[dplyr]{join}}.
#' @details As opposed to the \code{\link[dplyr]{join}} functions of \code{dplyr}, at default existing columns will get a suffix \code{"2"} and the newly joined columns will not get a suffix. See \code{\link[dplyr]{join}} for more information.
#' @export
#' @examples 
#' df <- data.frame(date = seq(from = as.Date("2018-01-01"),
#'                             to = as.Date("2018-01-07"),
#'                             by = 1),
#'                  bacteria_id = c("STAAUR", "STAAUR", "STAAUR", "STAAUR",
#'                                  "ESCCOL", "ESCCOL", "ESCCOL"),
#'                  stringsAsFactors = FALSE)
#'                  
#' colnames(df)
#' df2 <- left_join_bactlist(df, "bacteria_id")
#' colnames(df2)
inner_join_bactlist <- function(x, by = 'bactid', ...) {
  # no name set to `by` parameter
  if (is.null(names(by))) {
    joinby <- colnames(AMR::bactlist)[1]
    names(joinby) <- by
  } else {
    joinby <- by
  }
  join <- dplyr::inner_join(x = x, y = AMR::bactlist, by = joinby, suffix = c("2", ""), ...)
  if (nrow(join) > nrow(x)) {
    warning('the newly joined tbl contains ', nrow(join) - nrow(x), ' rows more that its original')
  }
  join
}

#' @rdname join
#' @export
left_join_bactlist <- function(x, by = 'bactid', ...) {
  # no name set to `by` parameter
  if (is.null(names(by))) {
    joinby <- colnames(AMR::bactlist)[1]
    names(joinby) <- by
  } else {
    joinby <- by
  }
  join <- dplyr::left_join(x = x, y = AMR::bactlist, by = joinby, suffix = c("2", ""), ...)
  if (nrow(join) > nrow(x)) {
    warning('the newly joined tbl contains ', nrow(join) - nrow(x), ' rows more that its original')
  }
  join
}

#' @rdname join
#' @export
right_join_bactlist <- function(x, by = 'bactid', ...) {
  # no name set to `by` parameter
  if (is.null(names(by))) {
    joinby <- colnames(AMR::bactlist)[1]
    names(joinby) <- by
  } else {
    joinby <- by
  }
  join <- dplyr::right_join(x = x, y = AMR::bactlist, by = joinby, suffix = c("2", ""), ...)
  if (nrow(join) > nrow(x)) {
    warning('the newly joined tbl contains ', nrow(join) - nrow(x), ' rows more that its original')
  }
  join
}

#' @rdname join
#' @export
full_join_bactlist <- function(x, by = 'bactid', ...) {
  # no name set to `by` parameter
  if (is.null(names(by))) {
    joinby <- colnames(AMR::bactlist)[1]
    names(joinby) <- by
  } else {
    joinby <- by
  }
  dplyr::full_join(x = x, y = AMR::bactlist, by = joinby, suffix = c("2", ""), ...)
}

#' @rdname join
#' @export
semi_join_bactlist <- function(x, by = 'bactid', ...) {
  # no name set to `by` parameter
  if (is.null(names(by))) {
    joinby <- colnames(AMR::bactlist)[1]
    names(joinby) <- by
  } else {
    joinby <- by
  }
  dplyr::semi_join(x = x, y = AMR::bactlist, by = joinby, ...)
}

#' @rdname join
#' @export
anti_join_bactlist <- function(x, by = 'bactid', ...) {
  # no name set to `by` parameter
  if (is.null(names(by))) {
    joinby <- colnames(AMR::bactlist)[1]
    names(joinby) <- by
  } else {
    joinby <- by
  }
  dplyr::anti_join(x = x, y = AMR::bactlist, by = joinby, ...)
}
