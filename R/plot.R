# ==================================================================== #
# TITLE:                                                               #
# AMR: An R Package for Working with Antimicrobial Resistance Data     #
#                                                                      #
# SOURCE CODE:                                                         #
# https://github.com/msberends/AMR                                     #
#                                                                      #
# PLEASE CITE THIS SOFTWARE AS:                                        #
# Berends MS, Luz CF, Friedrich AW, Sinha BNM, Albers CJ, Glasner C    #
# (2022). AMR: An R Package for Working with Antimicrobial Resistance  #
# Data. Journal of Statistical Software, 104(3), 1-31.                 #
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

#' Plotting for Classes `sir`, `mic` and `disk`
#'
#' @description
#' Functions to plot classes `sir`, `mic` and `disk`, with support for base \R and `ggplot2`.
#' 
#' Especially the `scale_*_mic()` functions are relevant wrappers to plot MIC values for `ggplot2`. They allows custom MIC ranges and to plot intermediate log2 levels for missing MIC values.
#' @param x,object values created with [as.mic()], [as.disk()] or [as.sir()] (or their `random_*` variants, such as [random_mic()])
#' @param mo any (vector of) text that can be coerced to a valid microorganism code with [as.mo()]
#' @param ab any (vector of) text that can be coerced to a valid antimicrobial drug code with [as.ab()]
#' @param guideline interpretation guideline to use - the default is the latest included EUCAST guideline, see *Details*
#' @param main,title title of the plot
#' @param xlab,ylab axis title
#' @param colours_SIR colours to use for filling in the bars, must be a vector of three values (in the order S, I and R). The default colours are colour-blind friendly.
#' @param language language to be used to translate 'Susceptible', 'Increased exposure'/'Intermediate' and 'Resistant' - the default is system language (see [get_AMR_locale()]) and can be overwritten by setting the [package option][AMR-options] [`AMR_locale`][AMR-options], e.g. `options(AMR_locale = "de")`, see [translate]. Use `language = NULL` or `language = ""` to prevent translation.
#' @param expand a [logical] to indicate whether the range on the x axis should be expanded between the lowest and highest value. For MIC values, intermediate values will be factors of 2 starting from the highest MIC value. For disk diameters, the whole diameter range will be filled.
#' @inheritParams as.sir
#' @details
#' The interpretation of "I" will be named "Increased exposure" for all EUCAST guidelines since 2019, and will be named "Intermediate" in all other cases.
#'
#' For interpreting MIC values as well as disk diffusion diameters, supported guidelines to be used as input for the `guideline` argument are: `r vector_and(AMR::clinical_breakpoints$guideline, quotes = TRUE, reverse = TRUE)`.
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
#' some_sir_values <- random_sir(50, prob_SIR = c(0.55, 0.05, 0.30))
#'
#' plot(some_mic_values)
#' plot(some_disk_values)
#' plot(some_sir_values)
#'
#' # when providing the microorganism and antibiotic, colours will show interpretations:
#' plot(some_mic_values, mo = "S. aureus", ab = "ampicillin")
#' plot(some_disk_values, mo = "Escherichia coli", ab = "cipro")
#' plot(some_disk_values, mo = "Escherichia coli", ab = "cipro", language = "nl")
#' 
#' 
#' # Plotting using scale_x_mic()
#' \donttest{
#' if (require("ggplot2")) {
#'   mic_plot <- ggplot(data.frame(mics = as.mic(c(0.125, "<=4", 4, 8, 32, ">=32")),
#'                                 counts = c(1, 1, 2, 2, 3, 3)),
#'                      aes(mics, counts)) +
#'     geom_col()
#'   mic_plot +
#'     labs(title = "without scale_x_mic()")
#' }
#' if (require("ggplot2")) {
#'   mic_plot +
#'     scale_x_mic() +
#'     labs(title = "with scale_x_mic()")
#' }
#' if (require("ggplot2")) {
#'   mic_plot +
#'     scale_x_mic(keep_operators = "all") +
#'     labs(title = "with scale_x_mic() keeping all operators")
#' }
#' if (require("ggplot2")) {
#'   mic_plot +
#'     scale_x_mic(mic_range = c(1, 128)) +
#'     labs(title = "with scale_x_mic() using a manual range")
#' }
#' 
#' if (require("ggplot2")) {
#'   autoplot(some_mic_values)
#' }
#' if (require("ggplot2")) {
#'   autoplot(some_disk_values, mo = "Escherichia coli", ab = "cipro")
#' }
#' if (require("ggplot2")) {
#'   autoplot(some_sir_values)
#' }
#' }
NULL

#' @export
#' @inheritParams as.mic
#' @param drop a [logical] to remove intermediate MIC values, defaults to `FALSE`
#' @rdname plot
scale_x_mic <- function(keep_operators = "edges", mic_range = NULL, drop = FALSE, ...) {
  stop_ifnot_installed("ggplot2")
  meet_criteria(drop, allow_class = "logical", has_length = 1)
  scale <- ggplot2::scale_x_discrete(drop = drop, ...)
  scale$transform <- function(x, keep_ops = keep_operators, mic_rng = mic_range) {
    rescale_mic(x = x, keep_operators = keep_ops, mic_range = mic_rng, as.mic = FALSE)
  }
  scale
}

#' @export
#' @inheritParams as.mic
#' @rdname plot
scale_y_mic <- function(keep_operators = "edges", mic_range = NULL, drop = FALSE, ...) {
  stop_ifnot_installed("ggplot2")
  meet_criteria(drop, allow_class = "logical", has_length = 1)
  scale <- ggplot2::scale_y_discrete(drop = drop, ...)
  scale$transform <- function(x, keep_ops = keep_operators, mic_rng = mic_range) {
    rescale_mic(x = x, keep_operators = keep_ops, mic_range = mic_rng, as.mic = FALSE)
  }
  scale
}

#' @export
#' @inheritParams as.mic
#' @rdname plot
scale_colour_mic <- function(keep_operators = "edges", mic_range = NULL, drop = FALSE, ...) {
  stop_ifnot_installed("ggplot2")
  meet_criteria(drop, allow_class = "logical", has_length = 1)
  scale <- ggplot2::scale_colour_discrete(drop = drop, ...)
  scale$transform <- function(x, keep_ops = keep_operators, mic_rng = mic_range) {
    rescale_mic(x = x, keep_operators = keep_ops, mic_range = mic_rng, as.mic = FALSE)
  }
  scale
}

#' @export
#' @inheritParams as.mic
#' @rdname plot
scale_fill_mic <- function(keep_operators = "edges", mic_range = NULL, drop = FALSE, ...) {
  stop_ifnot_installed("ggplot2")
  meet_criteria(drop, allow_class = "logical", has_length = 1)
  scale <- ggplot2::scale_fill_discrete(drop = drop, ...)
  scale$transform <- function(x, keep_ops = keep_operators, mic_rng = mic_range) {
    rescale_mic(x = x, keep_operators = keep_ops, mic_range = mic_rng, as.mic = FALSE)
  }
  scale
}

#' @method plot mic
#' @importFrom graphics barplot axis mtext legend
#' @export
#' @rdname plot
plot.mic <- function(x,
                     mo = NULL,
                     ab = NULL,
                     guideline = "EUCAST",
                     main = deparse(substitute(x)),
                     ylab = translate_AMR("Frequency", language = language),
                     xlab = translate_AMR("Minimum Inhibitory Concentration (mg/L)", language = language),
                     colours_SIR = c("#3CAEA3", "#F6D55C", "#ED553B"),
                     language = get_AMR_locale(),
                     expand = TRUE,
                     include_PKPD = getOption("AMR_include_PKPD", TRUE),
                     breakpoint_type = getOption("AMR_breakpoint_type", "human"),
                     ...) {
  meet_criteria(mo, allow_class = c("mo", "character"), allow_NULL = TRUE)
  meet_criteria(ab, allow_class = c("ab", "character"), allow_NULL = TRUE)
  meet_criteria(guideline, allow_class = "character", has_length = 1)
  meet_criteria(main, allow_class = "character", has_length = 1, allow_NULL = TRUE)
  meet_criteria(ylab, allow_class = "character", has_length = 1)
  meet_criteria(xlab, allow_class = "character", has_length = 1)
  meet_criteria(colours_SIR, allow_class = "character", has_length = c(1, 3))
  language <- validate_language(language)
  meet_criteria(expand, allow_class = "logical", has_length = 1)

  x <- as.mic(x) # make sure that currently implemented MIC levels are used
  
  if (length(colours_SIR) == 1) {
    colours_SIR <- rep(colours_SIR, 3)
  }
  main <- gsub(" +", " ", paste0(main, collapse = " "))

  x <- plotrange_as_table(x, expand = expand)
  cols_sub <- plot_colours_subtitle_guideline(
    x = x,
    mo = mo,
    ab = ab,
    guideline = guideline,
    colours_SIR = colours_SIR,
    fn = as.mic,
    language = language,
    method = "MIC",
    include_PKPD = include_PKPD,
    breakpoint_type = breakpoint_type,
    ...
  )
  barplot(x,
    col = cols_sub$cols,
    main = main,
    ylim = c(0, max(x) * ifelse(any(colours_SIR %in% cols_sub$cols), 1.1, 1)),
    ylab = ylab,
    xlab = xlab,
    axes = FALSE
  )
  axis(2, seq(0, max(x)))
  if (!is.null(cols_sub$sub)) {
    mtext(side = 3, line = 0.5, adj = 0.5, cex = 0.75, cols_sub$sub)
  }

  if (any(colours_SIR %in% cols_sub$cols)) {
    legend_txt <- character(0)
    legend_col <- character(0)
    if (any(cols_sub$cols == colours_SIR[1] & cols_sub$count > 0)) {
      legend_txt <- c(legend_txt, "(S) Susceptible")
      legend_col <- colours_SIR[1]
    }
    if (any(cols_sub$cols == colours_SIR[2] & cols_sub$count > 0)) {
      legend_txt <- c(legend_txt, paste("(I)", plot_name_of_I(cols_sub$guideline)))
      legend_col <- c(legend_col, colours_SIR[2])
    }
    if (any(cols_sub$cols == colours_SIR[3] & cols_sub$count > 0)) {
      legend_txt <- c(legend_txt, "(R) Resistant")
      legend_col <- c(legend_col, colours_SIR[3])
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
                        ylab = translate_AMR("Frequency", language = language),
                        xlab = translate_AMR("Minimum Inhibitory Concentration (mg/L)", language = language),
                        colours_SIR = c("#3CAEA3", "#F6D55C", "#ED553B"),
                        language = get_AMR_locale(),
                        expand = TRUE,
                        ...) {
  meet_criteria(main, allow_class = "character", has_length = 1, allow_NULL = TRUE)
  meet_criteria(ylab, allow_class = "character", has_length = 1)
  meet_criteria(xlab, allow_class = "character", has_length = 1)
  meet_criteria(mo, allow_class = c("mo", "character"), allow_NULL = TRUE)
  meet_criteria(ab, allow_class = c("ab", "character"), allow_NULL = TRUE)
  meet_criteria(guideline, allow_class = "character", has_length = 1)
  meet_criteria(colours_SIR, allow_class = "character", has_length = c(1, 3))
  language <- validate_language(language)
  meet_criteria(expand, allow_class = "logical", has_length = 1)

  main <- gsub(" +", " ", paste0(main, collapse = " "))
  
  height <- as.mic(height) # make sure that currently implemented MIC levels are used

  plot(
    x = height,
    main = main,
    ylab = ylab,
    xlab = xlab,
    mo = mo,
    ab = ab,
    guideline = guideline,
    colours_SIR = colours_SIR,
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
                         ylab = translate_AMR("Frequency", language = language),
                         xlab = translate_AMR("Minimum Inhibitory Concentration (mg/L)", language = language),
                         colours_SIR = c("#3CAEA3", "#F6D55C", "#ED553B"),
                         language = get_AMR_locale(),
                         expand = TRUE,
                         include_PKPD = getOption("AMR_include_PKPD", TRUE),
                         breakpoint_type = getOption("AMR_breakpoint_type", "human"),
                         ...) {
  stop_ifnot_installed("ggplot2")
  meet_criteria(mo, allow_class = c("mo", "character"), allow_NULL = TRUE)
  meet_criteria(ab, allow_class = c("ab", "character"), allow_NULL = TRUE)
  meet_criteria(guideline, allow_class = "character", has_length = 1)
  meet_criteria(title, allow_class = "character", allow_NULL = TRUE)
  meet_criteria(ylab, allow_class = "character", has_length = 1)
  meet_criteria(xlab, allow_class = "character", has_length = 1)
  meet_criteria(colours_SIR, allow_class = "character", has_length = c(1, 3))
  language <- validate_language(language)
  meet_criteria(expand, allow_class = "logical", has_length = 1)

  if ("main" %in% names(list(...))) {
    title <- list(...)$main
  }
  if (!is.null(title)) {
    title <- gsub(" +", " ", paste0(title, collapse = " "))
  }

  object <- as.mic(object) # make sure that currently implemented MIC levels are used
  x <- plotrange_as_table(object, expand = expand)
  cols_sub <- plot_colours_subtitle_guideline(
    x = x,
    mo = mo,
    ab = ab,
    guideline = guideline,
    colours_SIR = colours_SIR,
    fn = as.mic,
    language = language,
    method = "MIC",
    include_PKPD = include_PKPD,
    breakpoint_type = breakpoint_type,
    ...
  )
  df <- as.data.frame(x, stringsAsFactors = TRUE)
  colnames(df) <- c("mic", "count")
  df$cols <- cols_sub$cols
  df$cols[df$cols == colours_SIR[1]] <- "(S) Susceptible"
  df$cols[df$cols == colours_SIR[2]] <- paste("(I)", plot_name_of_I(cols_sub$guideline))
  df$cols[df$cols == colours_SIR[3]] <- "(R) Resistant"
  df$cols <- factor(translate_into_language(df$cols, language = language),
    levels = translate_into_language(
      c(
        "(S) Susceptible",
        paste("(I)", plot_name_of_I(cols_sub$guideline)),
        "(R) Resistant"
      ),
      language = language
    ),
    ordered = TRUE
  )
  p <- ggplot2::ggplot(df)

  if (any(colours_SIR %in% cols_sub$cols)) {
    vals <- c(
      "(S) Susceptible" = colours_SIR[1],
      "(SDD) Susceptible dose-dependent" = colours_SIR[2],
      "(I) Susceptible, incr. exp." = colours_SIR[2],
      "(I) Intermediate" = colours_SIR[2],
      "(R) Resistant" = colours_SIR[3]
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
  object <- as.mic(object) # make sure that currently implemented MIC levels are used
  stats::setNames(
    as.data.frame(plotrange_as_table(object, expand = FALSE)),
    c("x", "y")
  )
}


#' @method plot disk
#' @export
#' @importFrom graphics barplot axis mtext legend
#' @rdname plot
plot.disk <- function(x,
                      main = deparse(substitute(x)),
                      ylab = translate_AMR("Frequency", language = language),
                      xlab = translate_AMR("Disk diffusion diameter (mm)", language = language),
                      mo = NULL,
                      ab = NULL,
                      guideline = "EUCAST",
                      colours_SIR = c("#3CAEA3", "#F6D55C", "#ED553B"),
                      language = get_AMR_locale(),
                      expand = TRUE,
                      include_PKPD = getOption("AMR_include_PKPD", TRUE),
                      breakpoint_type = getOption("AMR_breakpoint_type", "human"),
                      ...) {
  meet_criteria(main, allow_class = "character", has_length = 1, allow_NULL = TRUE)
  meet_criteria(ylab, allow_class = "character", has_length = 1)
  meet_criteria(xlab, allow_class = "character", has_length = 1)
  meet_criteria(mo, allow_class = c("mo", "character"), allow_NULL = TRUE)
  meet_criteria(ab, allow_class = c("ab", "character"), allow_NULL = TRUE)
  meet_criteria(guideline, allow_class = "character", has_length = 1)
  meet_criteria(colours_SIR, allow_class = "character", has_length = c(1, 3))
  language <- validate_language(language)
  meet_criteria(expand, allow_class = "logical", has_length = 1)

  if (length(colours_SIR) == 1) {
    colours_SIR <- rep(colours_SIR, 3)
  }
  main <- gsub(" +", " ", paste0(main, collapse = " "))

  x <- plotrange_as_table(x, expand = expand)
  cols_sub <- plot_colours_subtitle_guideline(
    x = x,
    mo = mo,
    ab = ab,
    guideline = guideline,
    colours_SIR = colours_SIR,
    fn = as.disk,
    language = language,
    method = "disk",
    include_PKPD = include_PKPD,
    breakpoint_type = breakpoint_type,
    ...
  )

  barplot(x,
    col = cols_sub$cols,
    main = main,
    ylim = c(0, max(x) * ifelse(any(colours_SIR %in% cols_sub$cols), 1.1, 1)),
    ylab = ylab,
    xlab = xlab,
    axes = FALSE
  )
  axis(2, seq(0, max(x)))
  if (!is.null(cols_sub$sub)) {
    mtext(side = 3, line = 0.5, adj = 0.5, cex = 0.75, cols_sub$sub)
  }

  if (any(colours_SIR %in% cols_sub$cols)) {
    legend_txt <- character(0)
    legend_col <- character(0)
    if (any(cols_sub$cols == colours_SIR[3] & cols_sub$count > 0)) {
      legend_txt <- "(R) Resistant"
      legend_col <- colours_SIR[3]
    }
    if (any(cols_sub$cols == colours_SIR[2] & cols_sub$count > 0)) {
      legend_txt <- c(legend_txt, paste("(I)", plot_name_of_I(cols_sub$guideline)))
      legend_col <- c(legend_col, colours_SIR[2])
    }
    if (any(cols_sub$cols == colours_SIR[1] & cols_sub$count > 0)) {
      legend_txt <- c(legend_txt, "(S) Susceptible")
      legend_col <- c(legend_col, colours_SIR[1])
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
                         ylab = translate_AMR("Frequency", language = language),
                         xlab = translate_AMR("Disk diffusion diameter (mm)", language = language),
                         mo = NULL,
                         ab = NULL,
                         guideline = "EUCAST",
                         colours_SIR = c("#3CAEA3", "#F6D55C", "#ED553B"),
                         language = get_AMR_locale(),
                         expand = TRUE,
                         ...) {
  meet_criteria(main, allow_class = "character", has_length = 1, allow_NULL = TRUE)
  meet_criteria(ylab, allow_class = "character", has_length = 1)
  meet_criteria(xlab, allow_class = "character", has_length = 1)
  meet_criteria(mo, allow_class = c("mo", "character"), allow_NULL = TRUE)
  meet_criteria(ab, allow_class = c("ab", "character"), allow_NULL = TRUE)
  meet_criteria(guideline, allow_class = "character", has_length = 1)
  meet_criteria(colours_SIR, allow_class = "character", has_length = c(1, 3))
  language <- validate_language(language)
  meet_criteria(expand, allow_class = "logical", has_length = 1)

  main <- gsub(" +", " ", paste0(main, collapse = " "))

  plot(
    x = height,
    main = main,
    ylab = ylab,
    xlab = xlab,
    mo = mo,
    ab = ab,
    guideline = guideline,
    colours_SIR = colours_SIR,
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
                          ylab = translate_AMR("Frequency", language = language),
                          xlab = translate_AMR("Disk diffusion diameter (mm)", language = language),
                          guideline = "EUCAST",
                          colours_SIR = c("#3CAEA3", "#F6D55C", "#ED553B"),
                          language = get_AMR_locale(),
                          expand = TRUE,
                          include_PKPD = getOption("AMR_include_PKPD", TRUE),
                          breakpoint_type = getOption("AMR_breakpoint_type", "human"),
                          ...) {
  stop_ifnot_installed("ggplot2")
  meet_criteria(title, allow_class = "character", allow_NULL = TRUE)
  meet_criteria(ylab, allow_class = "character", has_length = 1)
  meet_criteria(xlab, allow_class = "character", has_length = 1)
  meet_criteria(mo, allow_class = c("mo", "character"), allow_NULL = TRUE)
  meet_criteria(ab, allow_class = c("ab", "character"), allow_NULL = TRUE)
  meet_criteria(guideline, allow_class = "character", has_length = 1)
  meet_criteria(colours_SIR, allow_class = "character", has_length = c(1, 3))
  language <- validate_language(language)
  meet_criteria(expand, allow_class = "logical", has_length = 1)

  if ("main" %in% names(list(...))) {
    title <- list(...)$main
  }
  if (!is.null(title)) {
    title <- gsub(" +", " ", paste0(title, collapse = " "))
  }

  x <- plotrange_as_table(object, expand = expand)
  cols_sub <- plot_colours_subtitle_guideline(
    x = x,
    mo = mo,
    ab = ab,
    guideline = guideline,
    colours_SIR = colours_SIR,
    fn = as.disk,
    language = language,
    method = "disk",
    include_PKPD = include_PKPD,
    breakpoint_type = breakpoint_type,
    ...
  )
  df <- as.data.frame(x, stringsAsFactors = TRUE)
  colnames(df) <- c("disk", "count")
  df$cols <- cols_sub$cols

  df$cols[df$cols == colours_SIR[1]] <- "(S) Susceptible"
  df$cols[df$cols == colours_SIR[2]] <- paste("(I)", plot_name_of_I(cols_sub$guideline))
  df$cols[df$cols == colours_SIR[3]] <- "(R) Resistant"
  df$cols <- factor(translate_into_language(df$cols, language = language),
    levels = translate_into_language(
      c(
        "(S) Susceptible",
        paste("(I)", plot_name_of_I(cols_sub$guideline)),
        "(R) Resistant"
      ),
      language = language
    ),
    ordered = TRUE
  )
  p <- ggplot2::ggplot(df)

  if (any(colours_SIR %in% cols_sub$cols)) {
    vals <- c(
      "(S) Susceptible" = colours_SIR[1],
      "(SDD) Susceptible dose-dependent" = colours_SIR[2],
      "(I) Susceptible, incr. exp." = colours_SIR[2],
      "(I) Intermediate" = colours_SIR[2],
      "(R) Resistant" = colours_SIR[3]
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
    as.data.frame(plotrange_as_table(object, expand = FALSE)),
    c("x", "y")
  )
}

#' @method plot sir
#' @export
#' @importFrom graphics plot text axis
#' @rdname plot
plot.sir <- function(x,
                     ylab = translate_AMR("Percentage", language = language),
                     xlab = translate_AMR("Antimicrobial Interpretation", language = language),
                     main = deparse(substitute(x)),
                     language = get_AMR_locale(),
                     ...) {
  meet_criteria(ylab, allow_class = "character", has_length = 1)
  meet_criteria(xlab, allow_class = "character", has_length = 1)
  meet_criteria(main, allow_class = "character", has_length = 1, allow_NULL = TRUE)

  data <- as.data.frame(table(x), stringsAsFactors = FALSE)
  colnames(data) <- c("x", "n")
  data$s <- round((data$n / sum(data$n)) * 100, 1)

  if (!"S" %in% data$x) {
    data <- rbind_AMR(data, data.frame(x = "S", n = 0, s = 0, stringsAsFactors = FALSE))
  }
  if (!"SDD" %in% data$x) {
    data <- rbind_AMR(data, data.frame(x = "SDD", n = 0, s = 0, stringsAsFactors = FALSE))
  }
  if (!"I" %in% data$x) {
    data <- rbind_AMR(data, data.frame(x = "I", n = 0, s = 0, stringsAsFactors = FALSE))
  }
  if (!"R" %in% data$x) {
    data <- rbind_AMR(data, data.frame(x = "R", n = 0, s = 0, stringsAsFactors = FALSE))
  }
  if (!"N" %in% data$x) {
    data <- rbind_AMR(data, data.frame(x = "N", n = 0, s = 0, stringsAsFactors = FALSE))
  }
  
  data <- data[!(data$n == 0 & data$x %in% c("SDD", "I", "N")), , drop = FALSE]
  data$x <- factor(data$x, levels = intersect(unique(data$x), c("S", "SDD", "I", "R", "N")), ordered = TRUE)

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


#' @method barplot sir
#' @importFrom graphics barplot axis
#' @export
#' @noRd
barplot.sir <- function(height,
                        main = deparse(substitute(height)),
                        xlab = translate_AMR("Antimicrobial Interpretation", language = language),
                        ylab = translate_AMR("Frequency", language = language),
                        colours_SIR = c("#3CAEA3", "#F6D55C", "#ED553B"),
                        language = get_AMR_locale(),
                        expand = TRUE,
                        ...) {
  meet_criteria(xlab, allow_class = "character", has_length = 1)
  meet_criteria(main, allow_class = "character", has_length = 1, allow_NULL = TRUE)
  meet_criteria(ylab, allow_class = "character", has_length = 1)
  meet_criteria(colours_SIR, allow_class = "character", has_length = c(1, 3))
  language <- validate_language(language)
  meet_criteria(expand, allow_class = "logical", has_length = 1)

  if (length(colours_SIR) == 1) {
    colours_SIR <- rep(colours_SIR, 3)
  }
  # add SDD and N to colours
  colours_SIR <- c(colours_SIR[1:2], colours_SIR[2], colours_SIR[3], "#888888")
  main <- gsub(" +", " ", paste0(main, collapse = " "))

  x <- table(height)
  # remove missing I, SDD, and N
  colours_SIR <- colours_SIR[!(names(x) %in% c("SDD", "I", "N") & x == 0)]
  x <- x[!(names(x) %in% c("SDD", "I", "N") & x == 0)]
  # plot it
  barplot(x,
    col = colours_SIR,
    xlab = xlab,
    main = main,
    ylab = ylab,
    axes = FALSE
  )
  axis(2, seq(0, max(x)))
}

#' @method autoplot sir
#' @rdname plot
# will be exported using s3_register() in R/zzz.R
autoplot.sir <- function(object,
                         title = deparse(substitute(object)),
                         xlab = translate_AMR("Antimicrobial Interpretation", language = language),
                         ylab = translate_AMR("Frequency", language = language),
                         colours_SIR = c("#3CAEA3", "#F6D55C", "#ED553B"),
                         language = get_AMR_locale(),
                         ...) {
  stop_ifnot_installed("ggplot2")
  meet_criteria(title, allow_class = "character", allow_NULL = TRUE)
  meet_criteria(ylab, allow_class = "character", has_length = 1)
  meet_criteria(xlab, allow_class = "character", has_length = 1)
  meet_criteria(colours_SIR, allow_class = "character", has_length = c(1, 3))

  if ("main" %in% names(list(...))) {
    title <- list(...)$main
  }
  if (!is.null(title)) {
    title <- gsub(" +", " ", paste0(title, collapse = " "))
  }

  if (length(colours_SIR) == 1) {
    colours_SIR <- rep(colours_SIR, 3)
  }

  df <- as.data.frame(table(object), stringsAsFactors = TRUE)
  colnames(df) <- c("x", "n")
  df <- df[!(df$n == 0 & df$x %in% c("SDD", "I", "N")), , drop = FALSE]
  ggplot2::ggplot(df) +
    ggplot2::geom_col(ggplot2::aes(x = x, y = n, fill = x)) +
    # limits = force is needed because of a ggplot2 >= 3.3.4 bug (#4511)
    ggplot2::scale_fill_manual(
      values = c(
        "S" = colours_SIR[1],
        "SDD" = colours_SIR[2],
        "I" = colours_SIR[2],
        "R" = colours_SIR[3],
        "N" = "#888888"
      ),
      limits = force
    ) +
    ggplot2::labs(title = title, x = xlab, y = ylab) +
    ggplot2::theme(legend.position = "none")
}

#' @method fortify sir
#' @rdname plot
# will be exported using s3_register() in R/zzz.R
fortify.sir <- function(object, ...) {
  stats::setNames(
    as.data.frame(table(object)),
    c("x", "y")
  )
}

plotrange_as_table <- function(x, expand, keep_operators = "all", mic_range = NULL) {
  x <- x[!is.na(x)]
  if (is.mic(x)) {
    x <- as.mic(x, keep_operators = keep_operators)
    if (expand == TRUE) {
      # expand range for MIC by adding common intermediate factors levels
      extra_range <- COMMON_MIC_VALUES[COMMON_MIC_VALUES > min(x, na.rm = TRUE) & COMMON_MIC_VALUES < max(x, na.rm = TRUE)]
      # remove the ones that are in 25% range of user values
      extra_range <- extra_range[!vapply(FUN.VALUE = logical(1), extra_range, function(r) any(abs(r - x) / x < 0.25, na.rm = TRUE))]
      nms <- extra_range
      extra_range <- rep(0, length(extra_range))
      names(extra_range) <- nms
      x <- table(droplevels(x, as.mic = FALSE))
      extra_range <- extra_range[!names(extra_range) %in% names(x) & names(extra_range) %in% VALID_MIC_LEVELS]
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

ggplot2_get_from_dots <- function(arg, default, ...) {
  dots <- list(...)
  if (!arg %in% names(dots)) {
    default
  } else {
    dots[[arg]]
  }
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

plot_colours_subtitle_guideline <- function(x, mo, ab, guideline, colours_SIR, fn, language, method, breakpoint_type, include_PKPD, ...) {
  stop_if(length(x) == 0, "no observations to plot", call = FALSE)
  
  guideline <- get_guideline(guideline, AMR::clinical_breakpoints)
  
  # store previous interpretations to backup
  sir_history <- AMR_env$sir_interpretation_history
  # and clear previous interpretations
  AMR_env$sir_interpretation_history <- AMR_env$sir_interpretation_history[0, , drop = FALSE]
  
  if (!is.null(mo) && !is.null(ab)) {
    # interpret and give colour based on MIC values
    mo <- as.mo(mo)
    moname <- mo_name(mo, language = language)
    ab <- as.ab(ab)
    abname <- ab_name(ab, language = language)
    
    sir <- suppressWarnings(suppressMessages(as.sir(fn(names(x)), mo = mo, ab = ab, guideline = guideline, include_screening = FALSE, include_PKPD = include_PKPD, breakpoint_type = breakpoint_type, ...)))
    guideline_txt <- guideline
    if (all(is.na(sir))) {
      sir_screening <- suppressWarnings(suppressMessages(as.sir(fn(names(x)), mo = mo, ab = ab, guideline = guideline, include_screening = TRUE, include_PKPD = include_PKPD, breakpoint_type = breakpoint_type, ...)))
      if (!all(is.na(sir_screening))) {
        message_(
          "Only ", guideline, " ", method, " interpretations found for ",
          ab_name(ab, language = NULL, tolower = TRUE), " in ", italicise(moname), " for screening"
        )
        sir <- sir_screening
        guideline_txt <- paste0("(Screen, ", guideline_txt, ")")
      } else {
        message_(
          "No ", guideline, " ", method, " interpretations found for ",
          ab_name(ab, language = NULL, tolower = TRUE), " in ", italicise(moname)
        )
        guideline_txt <- paste0("(", guideline_txt, ")")
      }
    } else {
      if (isTRUE(list(...)$uti)) {
        guideline_txt <- paste("UTIs,", guideline_txt)
      }
      ref_tbl <- paste0('"', unique(AMR_env$sir_interpretation_history$ref_table), '"', collapse = "/")
      guideline_txt <- paste0("(", guideline_txt, ": ", ref_tbl, ")")
    }
    cols <- character(length = length(sir))
    cols[is.na(sir)] <- "#BEBEBE"
    cols[sir == "S"] <- colours_SIR[1]
    cols[sir == "SDD"] <- colours_SIR[2]
    cols[sir == "I"] <- colours_SIR[2]
    cols[sir == "R"] <- colours_SIR[3]
    cols[sir == "N"] <- "#888888"
    sub <- bquote(.(abname) ~ "-" ~ italic(.(moname)) ~ .(guideline_txt))
  } else {
    cols <- "#BEBEBE"
    sub <- NULL
  }
  
  # restore previous interpretations to backup
  AMR_env$sir_interpretation_history <- sir_history
  
  list(cols = cols, count = as.double(x), sub = sub, guideline = guideline)
}
