# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Data Analysis for R                   #
#                                                                      #
# SOURCE                                                               #
# https://github.com/msberends/AMR                                     #
#                                                                      #
# LICENCE                                                              #
# (c) 2018-2022 Berends MS, Luz CF et al.                              #
# Developed at the University of Groningen, the Netherlands, in        #
# collaboration with non-profit organisations Certe Medical            #
# Diagnostics & Advice, and University Medical Center Groningen.       #
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

#' Plotting for Classes `rsi`, `mic` and `disk`
#'
#' Functions to plot classes `rsi`, `mic` and `disk`, with support for base \R and `ggplot2`.

#' @param x,object values created with [as.mic()], [as.disk()] or [as.rsi()] (or their `random_*` variants, such as [random_mic()])
#' @param mo any (vector of) text that can be coerced to a valid microorganism code with [as.mo()]
#' @param ab any (vector of) text that can be coerced to a valid antimicrobial code with [as.ab()]
#' @param guideline interpretation guideline to use, defaults to the latest included EUCAST guideline, see *Details*
#' @param main,title title of the plot
#' @param xlab,ylab axis title
#' @param colours_RSI colours to use for filling in the bars, must be a vector of three values (in the order R, S and I). The default colours are colour-blind friendly.
#' @param language language to be used to translate 'Susceptible', 'Increased exposure'/'Intermediate' and 'Resistant', defaults to system language (see [get_AMR_locale()]) and can be overwritten by setting the option `AMR_locale`, e.g. `options(AMR_locale = "de")`, see [translate]. Use `language = NULL` or `language = ""` to prevent translation.
#' @param expand a [logical] to indicate whether the range on the x axis should be expanded between the lowest and highest value. For MIC values, intermediate values will be factors of 2 starting from the highest MIC value. For disk diameters, the whole diameter range will be filled.
#' @details
#' The interpretation of "I" will be named "Increased exposure" for all EUCAST guidelines since 2019, and will be named "Intermediate" in all other cases.
#'
#' For interpreting MIC values as well as disk diffusion diameters, supported guidelines to be used as input for the `guideline` argument are: `r vector_and(AMR::rsi_translation$guideline, quotes = TRUE, reverse = TRUE)`.
#'
#' Simply using `"CLSI"` or `"EUCAST"` as input will automatically select the latest version of that guideline.
#' @name plot
#' @rdname plot
#' @return The `autoplot()` functions return a [`ggplot`][ggplot2::ggplot()] model that is extendible with any `ggplot2` function.
#'
#' The `fortify()` functions return a [data.frame] as an extension for usage in the [ggplot2::ggplot()] function.
#' @param ... arguments passed on to methods
#' @examples
#' some_mic_values <- random_mic(size = 100)
#' some_disk_values <- random_disk(size = 100, mo = "Escherichia coli", ab = "cipro")
#' some_rsi_values <- random_rsi(50, prob_RSI = c(0.30, 0.55, 0.05))
#'
#' plot(some_mic_values)
#' plot(some_disk_values)
#' plot(some_rsi_values)
#'
#' # when providing the microorganism and antibiotic, colours will show interpretations:
#' plot(some_mic_values, mo = "S. aureus", ab = "ampicillin")
#' plot(some_disk_values, mo = "Escherichia coli", ab = "cipro")
#' plot(some_disk_values, mo = "Escherichia coli", ab = "cipro", language = "uk")
#'
#' \donttest{
#' if (require("ggplot2")) {
#'   autoplot(some_mic_values)
#'   autoplot(some_disk_values, mo = "Escherichia coli", ab = "cipro")
#'   autoplot(some_rsi_values)
#' }
#' }
NULL

#' @method plot mic
#' @importFrom graphics barplot axis mtext legend
#' @export
#' @rdname plot
plot.mic <- function(x,
                     mo = NULL,
                     ab = NULL,
                     guideline = "EUCAST",
                     main = deparse(substitute(x)),
                     ylab = "Frequency",
                     xlab = "Minimum Inhibitory Concentration (mg/L)",
                     colours_RSI = c("#ED553B", "#3CAEA3", "#F6D55C"),
                     language = get_AMR_locale(),
                     expand = TRUE,
                     ...) {
  meet_criteria(mo, allow_class = c("mo", "character"), allow_NULL = TRUE)
  meet_criteria(ab, allow_class = c("ab", "character"), allow_NULL = TRUE)
  meet_criteria(guideline, allow_class = "character", has_length = 1)
  meet_criteria(main, allow_class = "character", has_length = 1, allow_NULL = TRUE)
  meet_criteria(ylab, allow_class = "character", has_length = 1)
  meet_criteria(xlab, allow_class = "character", has_length = 1)
  meet_criteria(colours_RSI, allow_class = "character", has_length = c(1, 3))
  meet_criteria(language, has_length = 1, is_in = c(LANGUAGES_SUPPORTED, ""), allow_NULL = TRUE, allow_NA = TRUE)
  meet_criteria(expand, allow_class = "logical", has_length = 1)

  # translate if not specifically set
  if (missing(ylab)) {
    ylab <- translate_into_language(ylab, language = language)
  }
  if (missing(xlab)) {
    xlab <- translate_into_language(xlab, language = language)
  }

  if (length(colours_RSI) == 1) {
    colours_RSI <- rep(colours_RSI, 3)
  }
  main <- gsub(" +", " ", paste0(main, collapse = " "))

  x <- plot_prepare_table(x, expand = expand)

  cols_sub <- plot_colours_subtitle_guideline(
    x = x,
    mo = mo,
    ab = ab,
    guideline = guideline,
    colours_RSI = colours_RSI,
    fn = as.mic,
    language = language,
    ...
  )
  barplot(x,
    col = cols_sub$cols,
    main = main,
    ylim = c(0, max(x) * ifelse(any(colours_RSI %in% cols_sub$cols), 1.1, 1)),
    ylab = ylab,
    xlab = xlab,
    axes = FALSE
  )
  axis(2, seq(0, max(x)))
  if (!is.null(cols_sub$sub)) {
    mtext(side = 3, line = 0.5, adj = 0.5, cex = 0.75, cols_sub$sub)
  }

  if (any(colours_RSI %in% cols_sub$cols)) {
    legend_txt <- character(0)
    legend_col <- character(0)
    if (any(cols_sub$cols == colours_RSI[2] & cols_sub$count > 0)) {
      legend_txt <- "Susceptible"
      legend_col <- colours_RSI[2]
    }
    if (any(cols_sub$cols == colours_RSI[3] & cols_sub$count > 0)) {
      legend_txt <- c(legend_txt, plot_name_of_I(cols_sub$guideline))
      legend_col <- c(legend_col, colours_RSI[3])
    }
    if (any(cols_sub$cols == colours_RSI[1] & cols_sub$count > 0)) {
      legend_txt <- c(legend_txt, "Resistant")
      legend_col <- c(legend_col, colours_RSI[1])
    }

    legend("top",
      x.intersp = 0.5,
      legend = translate_into_language(legend_txt, language = language),
      fill = legend_col,
      horiz = TRUE,
      cex = 0.75,
      box.lwd = 0,
      box.col = "#FFFFFF55",
      bg = "#FFFFFF55"
    )
  }
}

#' @method barplot mic
#' @export
#' @noRd
barplot.mic <- function(height,
                        mo = NULL,
                        ab = NULL,
                        guideline = "EUCAST",
                        main = deparse(substitute(height)),
                        ylab = "Frequency",
                        xlab = "Minimum Inhibitory Concentration (mg/L)",
                        colours_RSI = c("#ED553B", "#3CAEA3", "#F6D55C"),
                        language = get_AMR_locale(),
                        expand = TRUE,
                        ...) {
  meet_criteria(main, allow_class = "character", has_length = 1, allow_NULL = TRUE)
  meet_criteria(ylab, allow_class = "character", has_length = 1)
  meet_criteria(xlab, allow_class = "character", has_length = 1)
  meet_criteria(mo, allow_class = c("mo", "character"), allow_NULL = TRUE)
  meet_criteria(ab, allow_class = c("ab", "character"), allow_NULL = TRUE)
  meet_criteria(guideline, allow_class = "character", has_length = 1)
  meet_criteria(colours_RSI, allow_class = "character", has_length = c(1, 3))
  meet_criteria(language, has_length = 1, is_in = c(LANGUAGES_SUPPORTED, ""), allow_NULL = TRUE, allow_NA = TRUE)
  meet_criteria(expand, allow_class = "logical", has_length = 1)

  # translate if not specifically set
  if (missing(ylab)) {
    ylab <- translate_into_language(ylab, language = language)
  }
  if (missing(xlab)) {
    xlab <- translate_into_language(xlab, language = language)
  }

  main <- gsub(" +", " ", paste0(main, collapse = " "))

  plot(
    x = height,
    main = main,
    ylab = ylab,
    xlab = xlab,
    mo = mo,
    ab = ab,
    guideline = guideline,
    colours_RSI = colours_RSI,
    ...
  )
}

#' @method autoplot mic
#' @rdname plot
# will be exported using s3_register() in R/zzz.R
autoplot.mic <- function(object,
                         mo = NULL,
                         ab = NULL,
                         guideline = "EUCAST",
                         title = deparse(substitute(object)),
                         ylab = "Frequency",
                         xlab = "Minimum Inhibitory Concentration (mg/L)",
                         colours_RSI = c("#ED553B", "#3CAEA3", "#F6D55C"),
                         language = get_AMR_locale(),
                         expand = TRUE,
                         ...) {
  stop_ifnot_installed("ggplot2")
  meet_criteria(mo, allow_class = c("mo", "character"), allow_NULL = TRUE)
  meet_criteria(ab, allow_class = c("ab", "character"), allow_NULL = TRUE)
  meet_criteria(guideline, allow_class = "character", has_length = 1)
  meet_criteria(title, allow_class = "character", allow_NULL = TRUE)
  meet_criteria(ylab, allow_class = "character", has_length = 1)
  meet_criteria(xlab, allow_class = "character", has_length = 1)
  meet_criteria(colours_RSI, allow_class = "character", has_length = c(1, 3))
  meet_criteria(language, has_length = 1, is_in = c(LANGUAGES_SUPPORTED, ""), allow_NULL = TRUE, allow_NA = TRUE)
  meet_criteria(expand, allow_class = "logical", has_length = 1)

  # translate if not specifically set
  if (missing(ylab)) {
    ylab <- translate_into_language(ylab, language = language)
  }
  if (missing(xlab)) {
    xlab <- translate_into_language(xlab, language = language)
  }

  if ("main" %in% names(list(...))) {
    title <- list(...)$main
  }
  if (!is.null(title)) {
    title <- gsub(" +", " ", paste0(title, collapse = " "))
  }

  x <- plot_prepare_table(object, expand = expand)
  cols_sub <- plot_colours_subtitle_guideline(
    x = x,
    mo = mo,
    ab = ab,
    guideline = guideline,
    colours_RSI = colours_RSI,
    fn = as.mic,
    language = language,
    ...
  )
  df <- as.data.frame(x, stringsAsFactors = TRUE)
  colnames(df) <- c("mic", "count")
  df$cols <- cols_sub$cols
  df$cols[df$cols == colours_RSI[1]] <- "Resistant"
  df$cols[df$cols == colours_RSI[2]] <- "Susceptible"
  df$cols[df$cols == colours_RSI[3]] <- plot_name_of_I(cols_sub$guideline)
  df$cols <- factor(translate_into_language(df$cols, language = language),
    levels = translate_into_language(c("Susceptible", plot_name_of_I(cols_sub$guideline), "Resistant"),
      language = language
    ),
    ordered = TRUE
  )
  p <- ggplot2::ggplot(df)

  if (any(colours_RSI %in% cols_sub$cols)) {
    vals <- c(
      "Resistant" = colours_RSI[1],
      "Susceptible" = colours_RSI[2],
      "Susceptible, incr. exp." = colours_RSI[3],
      "Intermediate" = colours_RSI[3]
    )
    names(vals) <- translate_into_language(names(vals), language = language)
    p <- p +
      ggplot2::geom_col(ggplot2::aes(x = mic, y = count, fill = cols)) +
      # limits = force is needed because of a ggplot2 >= 3.3.4 bug (#4511)
      ggplot2::scale_fill_manual(
        values = vals,
        name = NULL,
        limits = force
      )
  } else {
    p <- p +
      ggplot2::geom_col(ggplot2::aes(x = mic, y = count))
  }

  p +
    ggplot2::labs(title = title, x = xlab, y = ylab, subtitle = cols_sub$sub)
}

#' @method fortify mic
#' @rdname plot
# will be exported using s3_register() in R/zzz.R
fortify.mic <- function(object, ...) {
  stats::setNames(
    as.data.frame(plot_prepare_table(object, expand = FALSE)),
    c("x", "y")
  )
}

#' @method plot disk
#' @export
#' @importFrom graphics barplot axis mtext legend
#' @rdname plot
plot.disk <- function(x,
                      main = deparse(substitute(x)),
                      ylab = "Frequency",
                      xlab = "Disk diffusion diameter (mm)",
                      mo = NULL,
                      ab = NULL,
                      guideline = "EUCAST",
                      colours_RSI = c("#ED553B", "#3CAEA3", "#F6D55C"),
                      language = get_AMR_locale(),
                      expand = TRUE,
                      ...) {
  meet_criteria(main, allow_class = "character", has_length = 1, allow_NULL = TRUE)
  meet_criteria(ylab, allow_class = "character", has_length = 1)
  meet_criteria(xlab, allow_class = "character", has_length = 1)
  meet_criteria(mo, allow_class = c("mo", "character"), allow_NULL = TRUE)
  meet_criteria(ab, allow_class = c("ab", "character"), allow_NULL = TRUE)
  meet_criteria(guideline, allow_class = "character", has_length = 1)
  meet_criteria(colours_RSI, allow_class = "character", has_length = c(1, 3))
  meet_criteria(language, has_length = 1, is_in = c(LANGUAGES_SUPPORTED, ""), allow_NULL = TRUE, allow_NA = TRUE)
  meet_criteria(expand, allow_class = "logical", has_length = 1)

  # translate if not specifically set
  if (missing(ylab)) {
    ylab <- translate_into_language(ylab, language = language)
  }
  if (missing(xlab)) {
    xlab <- translate_into_language(xlab, language = language)
  }

  if (length(colours_RSI) == 1) {
    colours_RSI <- rep(colours_RSI, 3)
  }
  main <- gsub(" +", " ", paste0(main, collapse = " "))

  x <- plot_prepare_table(x, expand = expand)

  cols_sub <- plot_colours_subtitle_guideline(
    x = x,
    mo = mo,
    ab = ab,
    guideline = guideline,
    colours_RSI = colours_RSI,
    fn = as.disk,
    language = language,
    ...
  )

  barplot(x,
    col = cols_sub$cols,
    main = main,
    ylim = c(0, max(x) * ifelse(any(colours_RSI %in% cols_sub$cols), 1.1, 1)),
    ylab = ylab,
    xlab = xlab,
    axes = FALSE
  )
  axis(2, seq(0, max(x)))
  if (!is.null(cols_sub$sub)) {
    mtext(side = 3, line = 0.5, adj = 0.5, cex = 0.75, cols_sub$sub)
  }

  if (any(colours_RSI %in% cols_sub$cols)) {
    legend_txt <- character(0)
    legend_col <- character(0)
    if (any(cols_sub$cols == colours_RSI[1] & cols_sub$count > 0)) {
      legend_txt <- "Resistant"
      legend_col <- colours_RSI[1]
    }
    if (any(cols_sub$cols == colours_RSI[3] & cols_sub$count > 0)) {
      legend_txt <- c(legend_txt, plot_name_of_I(cols_sub$guideline))
      legend_col <- c(legend_col, colours_RSI[3])
    }
    if (any(cols_sub$cols == colours_RSI[2] & cols_sub$count > 0)) {
      legend_txt <- c(legend_txt, "Susceptible")
      legend_col <- c(legend_col, colours_RSI[2])
    }
    legend("top",
      x.intersp = 0.5,
      legend = translate_into_language(legend_txt, language = language),
      fill = legend_col,
      horiz = TRUE,
      cex = 0.75,
      box.lwd = 0,
      box.col = "#FFFFFF55",
      bg = "#FFFFFF55"
    )
  }
}

#' @method barplot disk
#' @export
#' @noRd
barplot.disk <- function(height,
                         main = deparse(substitute(height)),
                         ylab = "Frequency",
                         xlab = "Disk diffusion diameter (mm)",
                         mo = NULL,
                         ab = NULL,
                         guideline = "EUCAST",
                         colours_RSI = c("#ED553B", "#3CAEA3", "#F6D55C"),
                         language = get_AMR_locale(),
                         expand = TRUE,
                         ...) {
  meet_criteria(main, allow_class = "character", has_length = 1, allow_NULL = TRUE)
  meet_criteria(ylab, allow_class = "character", has_length = 1)
  meet_criteria(xlab, allow_class = "character", has_length = 1)
  meet_criteria(mo, allow_class = c("mo", "character"), allow_NULL = TRUE)
  meet_criteria(ab, allow_class = c("ab", "character"), allow_NULL = TRUE)
  meet_criteria(guideline, allow_class = "character", has_length = 1)
  meet_criteria(colours_RSI, allow_class = "character", has_length = c(1, 3))
  meet_criteria(language, has_length = 1, is_in = c(LANGUAGES_SUPPORTED, ""), allow_NULL = TRUE, allow_NA = TRUE)
  meet_criteria(expand, allow_class = "logical", has_length = 1)

  # translate if not specifically set
  if (missing(ylab)) {
    ylab <- translate_into_language(ylab, language = language)
  }
  if (missing(xlab)) {
    xlab <- translate_into_language(xlab, language = language)
  }

  main <- gsub(" +", " ", paste0(main, collapse = " "))

  plot(
    x = height,
    main = main,
    ylab = ylab,
    xlab = xlab,
    mo = mo,
    ab = ab,
    guideline = guideline,
    colours_RSI = colours_RSI,
    ...
  )
}

#' @method autoplot disk
#' @rdname plot
# will be exported using s3_register() in R/zzz.R
autoplot.disk <- function(object,
                          mo = NULL,
                          ab = NULL,
                          title = deparse(substitute(object)),
                          ylab = "Frequency",
                          xlab = "Disk diffusion diameter (mm)",
                          guideline = "EUCAST",
                          colours_RSI = c("#ED553B", "#3CAEA3", "#F6D55C"),
                          language = get_AMR_locale(),
                          expand = TRUE,
                          ...) {
  stop_ifnot_installed("ggplot2")
  meet_criteria(title, allow_class = "character", allow_NULL = TRUE)
  meet_criteria(ylab, allow_class = "character", has_length = 1)
  meet_criteria(xlab, allow_class = "character", has_length = 1)
  meet_criteria(mo, allow_class = c("mo", "character"), allow_NULL = TRUE)
  meet_criteria(ab, allow_class = c("ab", "character"), allow_NULL = TRUE)
  meet_criteria(guideline, allow_class = "character", has_length = 1)
  meet_criteria(colours_RSI, allow_class = "character", has_length = c(1, 3))
  meet_criteria(language, has_length = 1, is_in = c(LANGUAGES_SUPPORTED, ""), allow_NULL = TRUE, allow_NA = TRUE)
  meet_criteria(expand, allow_class = "logical", has_length = 1)

  # translate if not specifically set
  if (missing(ylab)) {
    ylab <- translate_into_language(ylab, language = language)
  }
  if (missing(xlab)) {
    xlab <- translate_into_language(xlab, language = language)
  }

  if ("main" %in% names(list(...))) {
    title <- list(...)$main
  }
  if (!is.null(title)) {
    title <- gsub(" +", " ", paste0(title, collapse = " "))
  }

  x <- plot_prepare_table(object, expand = expand)
  cols_sub <- plot_colours_subtitle_guideline(
    x = x,
    mo = mo,
    ab = ab,
    guideline = guideline,
    colours_RSI = colours_RSI,
    fn = as.disk,
    language = language,
    ...
  )
  df <- as.data.frame(x, stringsAsFactors = TRUE)
  colnames(df) <- c("disk", "count")
  df$cols <- cols_sub$cols

  df$cols[df$cols == colours_RSI[1]] <- "Resistant"
  df$cols[df$cols == colours_RSI[2]] <- "Susceptible"
  df$cols[df$cols == colours_RSI[3]] <- plot_name_of_I(cols_sub$guideline)
  df$cols <- factor(translate_into_language(df$cols, language = language),
    levels = translate_into_language(c("Susceptible", plot_name_of_I(cols_sub$guideline), "Resistant"),
      language = language
    ),
    ordered = TRUE
  )
  p <- ggplot2::ggplot(df)

  if (any(colours_RSI %in% cols_sub$cols)) {
    vals <- c(
      "Resistant" = colours_RSI[1],
      "Susceptible" = colours_RSI[2],
      "Susceptible, incr. exp." = colours_RSI[3],
      "Intermediate" = colours_RSI[3]
    )
    names(vals) <- translate_into_language(names(vals), language = language)
    p <- p +
      ggplot2::geom_col(ggplot2::aes(x = disk, y = count, fill = cols)) +
      # limits = force is needed because of a ggplot2 >= 3.3.4 bug (#4511)
      ggplot2::scale_fill_manual(
        values = vals,
        name = NULL,
        limits = force
      )
  } else {
    p <- p +
      ggplot2::geom_col(ggplot2::aes(x = disk, y = count))
  }

  p +
    ggplot2::labs(title = title, x = xlab, y = ylab, subtitle = cols_sub$sub)
}

#' @method fortify disk
#' @rdname plot
# will be exported using s3_register() in R/zzz.R
fortify.disk <- function(object, ...) {
  stats::setNames(
    as.data.frame(plot_prepare_table(object, expand = FALSE)),
    c("x", "y")
  )
}

#' @method plot rsi
#' @export
#' @importFrom graphics plot text axis
#' @rdname plot
plot.rsi <- function(x,
                     ylab = "Percentage",
                     xlab = "Antimicrobial Interpretation",
                     main = deparse(substitute(x)),
                     language = get_AMR_locale(),
                     ...) {
  meet_criteria(ylab, allow_class = "character", has_length = 1)
  meet_criteria(xlab, allow_class = "character", has_length = 1)
  meet_criteria(main, allow_class = "character", has_length = 1, allow_NULL = TRUE)

  # translate if not specifically set
  if (missing(ylab)) {
    ylab <- translate_into_language(ylab, language = language)
  }
  if (missing(xlab)) {
    xlab <- translate_into_language(xlab, language = language)
  }

  data <- as.data.frame(table(x), stringsAsFactors = FALSE)
  colnames(data) <- c("x", "n")
  data$s <- round((data$n / sum(data$n)) * 100, 1)

  if (!"S" %in% data$x) {
    data <- rbind(data, data.frame(x = "S", n = 0, s = 0, stringsAsFactors = FALSE),
      stringsAsFactors = FALSE
    )
  }
  if (!"I" %in% data$x) {
    data <- rbind(data, data.frame(x = "I", n = 0, s = 0, stringsAsFactors = FALSE),
      stringsAsFactors = FALSE
    )
  }
  if (!"R" %in% data$x) {
    data <- rbind(data, data.frame(x = "R", n = 0, s = 0, stringsAsFactors = FALSE),
      stringsAsFactors = FALSE
    )
  }

  data$x <- factor(data$x, levels = c("S", "I", "R"), ordered = TRUE)

  ymax <- pm_if_else(max(data$s) > 95, 105, 100)

  plot(
    x = data$x,
    y = data$s,
    lwd = 2,
    ylim = c(0, ymax),
    ylab = ylab,
    xlab = xlab,
    main = main,
    axes = FALSE
  )
  # x axis
  axis(side = 1, at = 1:pm_n_distinct(data$x), labels = levels(data$x), lwd = 0)
  # y axis, 0-100%
  axis(side = 2, at = seq(0, 100, 5))

  text(
    x = data$x,
    y = data$s + 4,
    labels = paste0(data$s, "% (n = ", data$n, ")")
  )
}


#' @method barplot rsi
#' @importFrom graphics barplot axis
#' @export
#' @noRd
barplot.rsi <- function(height,
                        main = deparse(substitute(height)),
                        xlab = "Antimicrobial Interpretation",
                        ylab = "Frequency",
                        colours_RSI = c("#ED553B", "#3CAEA3", "#F6D55C"),
                        language = get_AMR_locale(),
                        expand = TRUE,
                        ...) {
  meet_criteria(xlab, allow_class = "character", has_length = 1)
  meet_criteria(main, allow_class = "character", has_length = 1, allow_NULL = TRUE)
  meet_criteria(ylab, allow_class = "character", has_length = 1)
  meet_criteria(colours_RSI, allow_class = "character", has_length = c(1, 3))
  meet_criteria(language, has_length = 1, is_in = c(LANGUAGES_SUPPORTED, ""), allow_NULL = TRUE, allow_NA = TRUE)
  meet_criteria(expand, allow_class = "logical", has_length = 1)

  # translate if not specifically set
  if (missing(ylab)) {
    ylab <- translate_into_language(ylab, language = language)
  }
  if (missing(xlab)) {
    xlab <- translate_into_language(xlab, language = language)
  }

  if (length(colours_RSI) == 1) {
    colours_RSI <- rep(colours_RSI, 3)
  } else {
    colours_RSI <- c(colours_RSI[2], colours_RSI[3], colours_RSI[1])
  }
  main <- gsub(" +", " ", paste0(main, collapse = " "))

  x <- table(height)
  x <- x[c(1, 2, 3)]
  barplot(x,
    col = colours_RSI,
    xlab = xlab,
    main = main,
    ylab = ylab,
    axes = FALSE
  )
  axis(2, seq(0, max(x)))
}

#' @method autoplot rsi
#' @rdname plot
# will be exported using s3_register() in R/zzz.R
autoplot.rsi <- function(object,
                         title = deparse(substitute(object)),
                         xlab = "Antimicrobial Interpretation",
                         ylab = "Frequency",
                         colours_RSI = c("#ED553B", "#3CAEA3", "#F6D55C"),
                         language = get_AMR_locale(),
                         ...) {
  stop_ifnot_installed("ggplot2")
  meet_criteria(title, allow_class = "character", allow_NULL = TRUE)
  meet_criteria(ylab, allow_class = "character", has_length = 1)
  meet_criteria(xlab, allow_class = "character", has_length = 1)
  meet_criteria(colours_RSI, allow_class = "character", has_length = c(1, 3))

  # translate if not specifically set
  if (missing(ylab)) {
    ylab <- translate_into_language(ylab, language = language)
  }
  if (missing(xlab)) {
    xlab <- translate_into_language(xlab, language = language)
  }

  if ("main" %in% names(list(...))) {
    title <- list(...)$main
  }
  if (!is.null(title)) {
    title <- gsub(" +", " ", paste0(title, collapse = " "))
  }

  if (length(colours_RSI) == 1) {
    colours_RSI <- rep(colours_RSI, 3)
  }

  df <- as.data.frame(table(object), stringsAsFactors = TRUE)
  colnames(df) <- c("rsi", "count")
  ggplot2::ggplot(df) +
    ggplot2::geom_col(ggplot2::aes(x = rsi, y = count, fill = rsi)) +
    # limits = force is needed because of a ggplot2 >= 3.3.4 bug (#4511)
    ggplot2::scale_fill_manual(
      values = c(
        "R" = colours_RSI[1],
        "S" = colours_RSI[2],
        "I" = colours_RSI[3]
      ),
      limits = force
    ) +
    ggplot2::labs(title = title, x = xlab, y = ylab) +
    ggplot2::theme(legend.position = "none")
}

#' @method fortify rsi
#' @rdname plot
# will be exported using s3_register() in R/zzz.R
fortify.rsi <- function(object, ...) {
  stats::setNames(
    as.data.frame(table(object)),
    c("x", "y")
  )
}

plot_prepare_table <- function(x, expand) {
  x <- x[!is.na(x)]
  stop_if(length(x) == 0, "no observations to plot", call = FALSE)
  if (is.mic(x)) {
    if (expand == TRUE) {
      # expand range for MIC by adding factors of 2 from lowest to highest so all MICs in between also print
      valid_lvls <- levels(x)
      extra_range <- max(x) / 2
      while (min(extra_range) / 2 > min(x)) {
        extra_range <- c(min(extra_range) / 2, extra_range)
      }
      nms <- extra_range
      extra_range <- rep(0, length(extra_range))
      names(extra_range) <- nms
      x <- table(droplevels(x, as.mic = FALSE))
      extra_range <- extra_range[!names(extra_range) %in% names(x) & names(extra_range) %in% valid_lvls]
      x <- as.table(c(x, extra_range))
    } else {
      x <- table(droplevels(x, as.mic = FALSE))
    }
    x <- x[order(as.double(as.mic(names(x))))]
  } else if (is.disk(x)) {
    if (expand == TRUE) {
      # expand range for disks from lowest to highest so all mm's in between also print
      extra_range <- rep(0, max(x) - min(x) - 1)
      names(extra_range) <- seq(min(x) + 1, max(x) - 1)
      x <- table(x)
      extra_range <- extra_range[!names(extra_range) %in% names(x)]
      x <- as.table(c(x, extra_range))
    } else {
      x <- table(x)
    }
    x <- x[order(as.double(names(x)))]
  }
  as.table(x)
}

plot_name_of_I <- function(guideline) {
  if (guideline %unlike% "CLSI" && as.double(gsub("[^0-9]+", "", guideline)) >= 2019) {
    # interpretation since 2019
    "Susceptible, incr. exp."
  } else {
    # interpretation until 2019
    "Intermediate"
  }
}

plot_colours_subtitle_guideline <- function(x, mo, ab, guideline, colours_RSI, fn, language, ...) {
  guideline <- get_guideline(guideline, AMR::rsi_translation)
  if (!is.null(mo) && !is.null(ab)) {
    # interpret and give colour based on MIC values
    mo <- as.mo(mo)
    ab <- as.ab(ab)
    rsi <- suppressWarnings(suppressMessages(as.rsi(fn(names(x)), mo = mo, ab = ab, guideline = guideline, ...)))
    cols <- character(length = length(rsi))
    cols[is.na(rsi)] <- "#BEBEBE"
    cols[rsi == "R"] <- colours_RSI[1]
    cols[rsi == "S"] <- colours_RSI[2]
    cols[rsi == "I"] <- colours_RSI[3]
    moname <- mo_name(mo, language = language)
    abname <- ab_name(ab, language = language)
    if (all(cols == "#BEBEBE")) {
      message_(
        "No ", guideline, " interpretations found for ",
        ab_name(ab, language = NULL, tolower = TRUE), " in ", moname
      )
      guideline_txt <- ""
    } else {
      guideline_txt <- guideline
      if (isTRUE(list(...)$uti)) {
        guideline_txt <- paste("UTIs,", guideline_txt)
      }
      guideline_txt <- paste0("(", guideline_txt, ")")
    }
    sub <- bquote(.(abname) ~ "-" ~ italic(.(moname)) ~ .(guideline_txt))
  } else {
    cols <- "#BEBEBE"
    sub <- NULL
  }
  list(cols = cols, count = as.double(x), sub = sub, guideline = guideline)
}
