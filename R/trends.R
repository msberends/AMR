#' Detect trends using Machine Learning
#'
#' Test text
#' @param data a \code{data.frame}
#' @param threshold_unique do not analyse more unique \code{threshold_unique} items per variable
#' @param na.rm a logical value indicating whether \code{NA} values should be stripped before the computation proceeds.
#' @param info print relevant combinations to console
#' @return A \code{list} with class \code{"trends"}
#' @importFrom stats na.omit
#' @importFrom broom tidy
# @export
trends <- function(data, threshold_unique = 30, na.rm = TRUE, info = TRUE) {

  cols <- colnames(data)
  relevant <- list()
  count <- 0
  for (x in 1:length(cols)) {
    for (y in 1:length(cols)) {
      if (x == y) {
        next
      }
      if (n_distinct(data[, x]) > threshold_unique | n_distinct(data[, y]) > threshold_unique) {
        next
      }
      count <- count + 1
      df <- data %>%
        group_by_at(c(cols[x], cols[y])) %>%
        summarise(n = n())
      n <- df %>% pull(n)
      # linear regression model
      lin <- stats::lm(1:length(n) ~ n, na.action = ifelse(na.rm == TRUE, na.omit, NULL))

      res <- list(
        df = df,
        x = cols[x],
        y = cols[y],
        m = base::mean(n, na.rm = na.rm),
        sd = stats::sd(n, na.rm = na.rm),
        cv = cv(n, na.rm = na.rm),
        cqv = cqv(n, na.rm = na.rm),
        kurtosis = kurtosis(n, na.rm = na.rm),
        skewness = skewness(n, na.rm = na.rm),
        lin.p = broom::tidy(lin)[2, 'p.value']
        #binom.p <- broom::tidy(binom)[2, 'p.value']
      )

      include <- TRUE
      # ML part
      if (res$cv > 0.25) {
        res$reason <- "cv > 0.25"
      } else if (res$cqv > 0.75) {
        res$reason <- "cqv > 0.75"
      } else {
        include <- FALSE
      }

      if (include == TRUE) {
        relevant <- c(relevant, list(res))
        if (info == TRUE) {
          # minus one because the whole data will be added later
          cat(paste0("[", length(relevant), "]"), "Relevant:", cols[x], "vs.", cols[y], "\n")
        }
      }

    }
  }

  cat("Total of", count, "combinations analysed;", length(relevant), "seem relevant.\n")
  class(relevant) <- 'trends'
  relevant <- c(relevant, list(data = data))
  relevant

}

# @exportMethod print.trends
# @export
#' @noRd
print.trends <- function(x, ...) {
  cat(length(x) - 1, "relevant trends, out of", length(x$data)^2, "\n")
}

# @exportMethod plot.trends
# @export
#' @noRd
# plot.trends <- function(x, n = NULL, ...) {
#   if (is.null(n)) {
#     oask <- devAskNewPage(TRUE)
#     on.exit(devAskNewPage(oask))
#     n <- c(1:(length(x) - 1))
#   } else {
#     if (n > length(x) - 1) {
#       stop('trend unavailable, max is ', length(x) - 1, call. = FALSE)
#     }
#     oask <- NULL
#   }
#   for (i in n) {
#     data <- x[[i]]$df
#     if (as.character(i) %like% '1$') {
#       suffix <- "st"
#     } else if (as.character(i) %like% '2$') {
#       suffix <- "nd"
#     } else if (as.character(i) %like% '3$') {
#       suffix <- "rd"
#     } else {
#       suffix <- "th"
#     }
#     if (!is.null(oask)) {
#       cat(paste("Coming up:", colnames(data)[1], "vs.", colnames(data)[2]), "\n")
#     }
#     print(
#       ggplot(
#         data,
#         aes_string(x = colnames(data)[1],
#                    y = colnames(data)[3],
#                    group = colnames(data)[2],
#                    fill = colnames(data)[2])) +
#         geom_col(position = "dodge") +
#         theme_minimal() +
#         labs(title = paste(colnames(data)[1], "vs.", colnames(data)[2]),
#              subtitle = paste0(i, suffix, " trend"))
#     )
#   }
# }
