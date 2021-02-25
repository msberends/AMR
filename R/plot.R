# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Data Analysis for R                   #
#                                                                      #
# SOURCE                                                               #
# https://github.com/msberends/AMR                                     #
#                                                                      #
# LICENCE                                                              #
# (c) 2018-2021 Berends MS, Luz CF et al.                              #
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
#' Functions to plot classes `rsi`, `mic` and `disk`, with support for base R and `ggplot2`.
#' @inheritSection lifecycle Stable Lifecycle
#' @inheritSection AMR Read more on Our Website!
#' @param x MIC values created with [as.mic()] or disk diffusion values created with [as.disk()]
#' @param mapping aesthetic mappings to use for [`ggplot()`][ggplot2::ggplot()]
#' @param main,title title of the plot
#' @param xlab,ylab axis title
#' @param mo any (vector of) text that can be coerced to a valid microorganism code with [as.mo()]
#' @param ab any (vector of) text that can be coerced to a valid antimicrobial code with [as.ab()]
#' @param guideline interpretation guideline to use, defaults to the latest included EUCAST guideline, see *Details*
#' @param colours_RSI colours to use for filling in the bars, must be a vector of three values (in the order R, S and I). The default colours are colour-blind friendly.
#' @param expand logical to indicate whether the range on the x axis should be expanded between the lowest and highest value. For MIC values, intermediate values will be factors of 2 starting from the highest MIC value. For disk diameters, the whole diameter range will be filled.
#' @details For interpreting MIC values as well as disk diffusion diameters, supported guidelines to be used as input for the `guideline` argument are: `r vector_and(AMR::rsi_translation$guideline, quotes = TRUE, reverse = TRUE)`.
#' 
#' Simply using `"CLSI"` or `"EUCAST"` as input will automatically select the latest version of that guideline.
#' @name plot
#' @rdname plot
#' @return The `ggplot` functions return a [`ggplot`][ggplot2::ggplot()] model that is extendible with any `ggplot2` function.
#' @param ... arguments passed on to [as.rsi()]
#' @examples 
#' some_mic_values <- random_mic(size = 100)
#' some_disk_values <- random_disk(size = 100, mo = "Escherichia coli", ab = "cipro")
#' 
#' plot(some_mic_values)
#' plot(some_disk_values)
#' 
#' # when providing the microorganism and antibiotic, colours will show interpretations:
#' plot(some_mic_values, mo = "S. aureus", ab = "ampicillin")
#' plot(some_disk_values, mo = "Escherichia coli", ab = "cipro")
#' 
#' if (require("ggplot2")) {
#'   ggplot(some_mic_values)
#'   ggplot(some_disk_values, mo = "Escherichia coli", ab = "cipro")
#' }
NULL

#' @method plot mic
#' @importFrom graphics barplot axis mtext
#' @export
#' @rdname plot
plot.mic <- function(x,
                     main = paste("MIC values of", deparse(substitute(x))),
                     ylab = "Frequency",
                     xlab = "Minimum Inhibitory Concentration (mg/L)",
                     mo = NULL,
                     ab = NULL,
                     guideline = "EUCAST",
                     colours_RSI = c("#ED553B", "#3CAEA3", "#F6D55C"),
                     expand = TRUE,
                     ...) {
  meet_criteria(main, allow_class = "character")
  meet_criteria(ylab, allow_class = "character", has_length = 1)
  meet_criteria(xlab, allow_class = "character", has_length = 1)
  meet_criteria(mo, allow_class = c("mo", "character"), allow_NULL = TRUE)
  meet_criteria(ab, allow_class = c("ab", "character"), allow_NULL = TRUE)
  meet_criteria(guideline, allow_class = "character", has_length = 1)
  meet_criteria(colours_RSI, allow_class = "character", has_length = c(1, 3))
  if (length(colours_RSI) == 1) {
    colours_RSI <- rep(colours_RSI, 3)
  }
  main <- gsub(" +", " ", paste0(main, collapse = " "))
  
  x <- plot_prepare_table(x, expand = expand)
  
  cols_sub <- plot_colours_and_sub(x = x, 
                                   mo = mo, 
                                   ab = ab,
                                   guideline = guideline,
                                   colours_RSI = colours_RSI, 
                                   fn = as.mic,
                                   ...)
  
  barplot(x,
          col = cols_sub$cols,
          main = main,
          ylim = c(0, max(x) * ifelse(any(colours_RSI %in% cols_sub$cols), 1.1, 1)),
          ylab = ylab,
          xlab = xlab,
          axes = FALSE)
  axis(2, seq(0, max(as.double(x))))
  if (!is.null(cols_sub$sub)) {
    mtext(side = 3, line = 0.5, adj = 0.5, cex = 0.75, cols_sub$sub)
  }
  
  if (any(colours_RSI %in% cols_sub$cols)) {
    legend_txt <- character(0)
    legend_col <- character(0)
    if (colours_RSI[2] %in% cols_sub$cols) {
      legend_txt <- "Susceptible"
      legend_col <- colours_RSI[2]
    }
    if (colours_RSI[3] %in% cols_sub$cols) {
      legend_txt <- c(legend_txt, "Incr. exposure")
      legend_col <- c(legend_col, colours_RSI[3])
    }
    if (colours_RSI[1] %in% cols_sub$cols) {
      legend_txt <- c(legend_txt, "Resistant")
      legend_col <- c(legend_col, colours_RSI[1])
    }
    legend("top", 
           x.intersp = 0.5,
           legend = legend_txt,
           fill = legend_col,
           horiz = TRUE,
           cex = 0.75, 
           box.lwd = 0, 
           bg = "#FFFFFF55")
  }
}

#' @method barplot mic
#' @export
#' @noRd
barplot.mic <- function(height,
                        main = paste("MIC values of", deparse(substitute(height))),
                        ylab = "Frequency",
                        xlab = "Minimum Inhibitory Concentration (mg/L)",
                        mo = NULL,
                        ab = NULL,
                        guideline = "EUCAST",
                        colours_RSI = c("#ED553B", "#3CAEA3", "#F6D55C"),
                        expand = TRUE,
                        ...) {
  meet_criteria(main, allow_class = "character")
  meet_criteria(ylab, allow_class = "character", has_length = 1)
  meet_criteria(xlab, allow_class = "character", has_length = 1)
  meet_criteria(mo, allow_class = c("mo", "character"), allow_NULL = TRUE)
  meet_criteria(ab, allow_class = c("ab", "character"), allow_NULL = TRUE)
  meet_criteria(guideline, allow_class = "character", has_length = 1)
  meet_criteria(colours_RSI, allow_class = "character", has_length = c(1, 3))
  main <- gsub(" +", " ", paste0(main, collapse = " "))
  
  plot(x = height,
       main = main,
       ylab = ylab,
       xlab = xlab,
       mo = mo,
       ab = ab,
       guideline = guideline,
       colours_RSI = colours_RSI,
       ...)
}

#' @method ggplot mic
#' @rdname plot
# will be exported using s3_register() in R/zzz.R
ggplot.mic <- function(data,
                       mapping = NULL,
                       title = paste("MIC values of", deparse(substitute(data))),
                       ylab = "Frequency",
                       xlab = "Minimum Inhibitory Concentration (mg/L)",
                       mo = NULL,
                       ab = NULL,
                       guideline = "EUCAST",
                       colours_RSI = c("#ED553B", "#3CAEA3", "#F6D55C"),
                       expand = TRUE,
                       ...) {
  stop_ifnot_installed("ggplot2")
  meet_criteria(title, allow_class = "character")
  meet_criteria(ylab, allow_class = "character", has_length = 1)
  meet_criteria(xlab, allow_class = "character", has_length = 1)
  meet_criteria(mo, allow_class = c("mo", "character"), allow_NULL = TRUE)
  meet_criteria(ab, allow_class = c("ab", "character"), allow_NULL = TRUE)
  meet_criteria(guideline, allow_class = "character", has_length = 1)
  meet_criteria(colours_RSI, allow_class = "character", has_length = c(1, 3))
  
  title <- gsub(" +", " ", paste0(title, collapse = " "))
  
  x <- plot_prepare_table(data, expand = expand)
  cols_sub <- plot_colours_and_sub(x = x, 
                                   mo = mo, 
                                   ab = ab,
                                   guideline = guideline,
                                   colours_RSI = colours_RSI, 
                                   fn = as.mic,
                                   ...)
  df <- as.data.frame(x, stringsAsFactors = TRUE)
  colnames(df) <- c("mic", "count")
  df$cols <- cols_sub$cols
  df$cols[df$cols == colours_RSI[1]] <- "Resistant"
  df$cols[df$cols == colours_RSI[2]] <- "Susceptible"
  df$cols[df$cols == colours_RSI[3]] <- "Incr. exposure"
  df$cols <- factor(df$cols, 
                    levels = c("Susceptible", "Incr. exposure", "Resistant"),
                    ordered = TRUE)
  if (!is.null(mapping)) {
    p <- ggplot2::ggplot(df, mapping = mapping)
  } else {
    p <- ggplot2::ggplot(df)
  }
  
  if (any(colours_RSI %in% cols_sub$cols)) {
    p <- p +
      ggplot2::geom_col(aes(x = mic, y = count, fill = cols)) + 
      ggplot2::scale_fill_manual(values = c("Resistant" = colours_RSI[1],
                                            "Susceptible" = colours_RSI[2],
                                            "Incr. exposure" = colours_RSI[3]),,
                                 name = NULL)
  } else {
    p <- p +
      ggplot2::geom_col(aes(x = mic, y = count))
  }
  
  p +
    ggplot2::labs(title = title, x = xlab, y = ylab, subtitle = cols_sub$sub)
}


#' @method plot disk
#' @export
#' @importFrom graphics barplot axis mtext
#' @rdname plot
plot.disk <- function(x,
                      main = paste("Disk zones values of", deparse(substitute(x))),
                      ylab = "Frequency",
                      xlab = "Disk diffusion diameter (mm)",
                      mo = NULL,
                      ab = NULL,
                      guideline = "EUCAST",
                      colours_RSI = c("#ED553B", "#3CAEA3", "#F6D55C"),
                      expand = TRUE,
                      ...) {
  meet_criteria(main, allow_class = "character")
  meet_criteria(ylab, allow_class = "character", has_length = 1)
  meet_criteria(xlab, allow_class = "character", has_length = 1)
  meet_criteria(mo, allow_class = c("mo", "character"), allow_NULL = TRUE)
  meet_criteria(ab, allow_class = c("ab", "character"), allow_NULL = TRUE)
  meet_criteria(guideline, allow_class = "character", has_length = 1)
  meet_criteria(colours_RSI, allow_class = "character", has_length = c(1, 3))
  if (length(colours_RSI) == 1) {
    colours_RSI <- rep(colours_RSI, 3)
  }
  main <- gsub(" +", " ", paste0(main, collapse = " "))
  
  x <- plot_prepare_table(x, expand = expand)
  
  cols_sub <- plot_colours_and_sub(x = x, 
                                   mo = mo, 
                                   ab = ab,
                                   guideline = guideline,
                                   colours_RSI = colours_RSI, 
                                   fn = as.disk,
                                   ...)
  
  barplot(x,
          col = cols_sub$cols,
          main = main,
          ylim = c(0, max(x) * ifelse(any(colours_RSI %in% cols_sub$cols), 1.1, 1)),
          ylab = ylab,
          xlab = xlab,
          axes = FALSE)
  axis(2, seq(0, max(x)))
  if (!is.null(cols_sub$sub)) {
    mtext(side = 3, line = 0.5, adj = 0.5, cex = 0.75, cols_sub$sub)
  }
  
  if (any(colours_RSI %in% cols_sub$cols)) {
    legend_txt <- character(0)
    legend_col <- character(0)
    if (colours_RSI[1] %in% cols_sub$cols) {
      legend_txt <- "Resistant"
      legend_col <- colours_RSI[1]
    }
    if (colours_RSI[3] %in% cols_sub$cols) {
      legend_txt <- c(legend_txt, "Incr. exposure")
      legend_col <- c(legend_col, colours_RSI[3])
    }
    if (colours_RSI[2] %in% cols_sub$cols) {
      legend_txt <- c(legend_txt, "Susceptible")
      legend_col <- c(legend_col, colours_RSI[2])
    }
    legend("top", 
           x.intersp = 0.5,
           legend = legend_txt,
           fill = legend_col,
           horiz = TRUE,
           cex = 0.75, 
           box.lwd = 0, 
           bg = "#FFFFFF55")
  }
}

#' @method barplot disk
#' @export
#' @noRd
barplot.disk <- function(height,
                         main = paste("Disk zones values of", deparse(substitute(height))),
                         ylab = "Frequency",
                         xlab = "Disk diffusion diameter (mm)",
                         mo = NULL,
                         ab = NULL,
                         guideline = "EUCAST",
                         colours_RSI = c("#ED553B", "#3CAEA3", "#F6D55C"),
                         expand = TRUE,
                         ...) {
  meet_criteria(main, allow_class = "character")
  meet_criteria(ylab, allow_class = "character", has_length = 1)
  meet_criteria(xlab, allow_class = "character", has_length = 1)
  meet_criteria(mo, allow_class = c("mo", "character"), allow_NULL = TRUE)
  meet_criteria(ab, allow_class = c("ab", "character"), allow_NULL = TRUE)
  meet_criteria(guideline, allow_class = "character", has_length = 1)
  meet_criteria(colours_RSI, allow_class = "character", has_length = c(1, 3))
  
  main <- gsub(" +", " ", paste0(main, collapse = " "))
  
  plot(x = height,
       main = main,
       ylab = ylab,
       xlab = xlab,
       mo = mo,
       ab = ab,
       guideline = guideline,
       colours_RSI = colours_RSI,
       ...)
}

#' @method ggplot disk
#' @rdname plot
# will be exported using s3_register() in R/zzz.R
ggplot.disk <- function(data,
                        mapping = NULL,
                        title = paste("Disk zones values of", deparse(substitute(data))),
                        ylab = "Frequency",
                        xlab = "Disk diffusion diameter (mm)",
                        mo = NULL,
                        ab = NULL,
                        guideline = "EUCAST",
                        colours_RSI = c("#ED553B", "#3CAEA3", "#F6D55C"),
                        expand = TRUE,
                        ...) {
  stop_ifnot_installed("ggplot2")
  meet_criteria(title, allow_class = "character")
  meet_criteria(ylab, allow_class = "character", has_length = 1)
  meet_criteria(xlab, allow_class = "character", has_length = 1)
  meet_criteria(mo, allow_class = c("mo", "character"), allow_NULL = TRUE)
  meet_criteria(ab, allow_class = c("ab", "character"), allow_NULL = TRUE)
  meet_criteria(guideline, allow_class = "character", has_length = 1)
  meet_criteria(colours_RSI, allow_class = "character", has_length = c(1, 3))
  
  title <- gsub(" +", " ", paste0(title, collapse = " "))
  
  x <- plot_prepare_table(data, expand = expand)
  cols_sub <- plot_colours_and_sub(x = x,
                                   mo = mo,
                                   ab = ab,
                                   guideline = guideline,
                                   colours_RSI = colours_RSI,
                                   fn = as.disk,
                                   ...)
  df <- as.data.frame(x, stringsAsFactors = TRUE)
  colnames(df) <- c("disk", "count")
  df$cols <- cols_sub$cols
  df$cols[df$cols == colours_RSI[1]] <- "Resistant"
  df$cols[df$cols == colours_RSI[2]] <- "Susceptible"
  df$cols[df$cols == colours_RSI[3]] <- "Incr. exposure"
  df$cols <- factor(df$cols, 
                    levels = c("Resistant", "Incr. exposure", "Susceptible"),
                    ordered = TRUE)
  if (!is.null(mapping)) {
    p <- ggplot2::ggplot(df, mapping = mapping)
  } else {
    p <- ggplot2::ggplot(df)
  }
  
  if (any(colours_RSI %in% cols_sub$cols)) {
    p <- p +
      ggplot2::geom_col(aes(x = disk, y = count, fill = cols)) + 
      ggplot2::scale_fill_manual(values = c("Resistant" = colours_RSI[1],
                                            "Susceptible" = colours_RSI[2],
                                            "Incr. exposure" = colours_RSI[3]),
                                 name = NULL)
  } else {
    p <- p +
      ggplot2::geom_col(aes(x = disk, y = count))
  }
  
  p +
    ggplot2::labs(title = title, x = xlab, y = ylab, sub = cols_sub$sub)
}

plot_prepare_table <- function(x, expand) {
  if (is.mic(x)) {
    if (expand == TRUE) {
      # expand range for MIC by adding factors of 2 from lowest to highest so all MICs in between also print
      extra_range <- max(as.double(x)) / 2
      while (min(extra_range) / 2 > min(as.double(x))) {
        extra_range <- c(min(extra_range) / 2, extra_range)
      }
      extra_range <- setNames(rep(0, length(extra_range)), extra_range)
      x <- table(droplevels(x, as.mic = FALSE))
      extra_range <- extra_range[!names(extra_range) %in% names(x)]
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

plot_colours_and_sub <- function(x, mo, ab, guideline, colours_RSI, fn, ...) {
  if (!is.null(mo) && !is.null(ab)) {
    # interpret and give colour based on MIC values
    mo <- as.mo(mo)
    ab <- as.ab(ab)
    guideline <- get_guideline(guideline, AMR::rsi_translation)
    rsi <- suppressWarnings(suppressMessages(as.rsi(fn(names(x)), mo = mo, ab = ab, guideline = guideline, ...)))
    cols <- character(length = length(rsi))
    cols[is.na(rsi)] <- "#BEBEBE"
    cols[rsi == "R"] <- colours_RSI[1]
    cols[rsi == "S"] <- colours_RSI[2]
    cols[rsi == "I"] <- colours_RSI[3]
    moname <- mo_name(mo, language = NULL)
    abname <- ab_name(ab, language = NULL)
    if (all(cols == "#BEBEBE")) {
      message_("No ", guideline, " interpretations found for ", 
               ab_name(ab, language = NULL, tolower = TRUE), " in ", moname)
      guideline <- ""
    } else {
      guideline <- paste0("(following ", guideline, ")")
    }
    sub <- bquote(.(abname)~"in"~italic(.(moname))~.(guideline))
  } else {
    cols <- "#BEBEBE"
    sub <- NULL
  }
  list(cols = cols, sub = sub)
}


#' @method plot rsi
#' @export
#' @importFrom graphics plot text axis
#' @rdname plot
plot.rsi <- function(x,
                     ylab = "Percentage",
                     xlab = "Antimicrobial Interpretation",
                     main = paste("Resistance Overview of", deparse(substitute(x))),
                     ...) {
  meet_criteria(ylab, allow_class = "character", has_length = 1)
  meet_criteria(xlab, allow_class = "character", has_length = 1)
  meet_criteria(main, allow_class = "character", has_length = 1)
  
  data <- as.data.frame(table(x), stringsAsFactors = FALSE)
  colnames(data) <- c("x", "n")
  data$s <- round((data$n / sum(data$n)) * 100, 1)
  
  if (!"S" %in% data$x) {
    data <- rbind(data, data.frame(x = "S", n = 0, s = 0, stringsAsFactors = FALSE),
                  stringsAsFactors = FALSE)
  }
  if (!"I" %in% data$x) {
    data <- rbind(data, data.frame(x = "I", n = 0, s = 0, stringsAsFactors = FALSE),
                  stringsAsFactors = FALSE)
  }
  if (!"R" %in% data$x) {
    data <- rbind(data, data.frame(x = "R", n = 0, s = 0, stringsAsFactors = FALSE),
                  stringsAsFactors = FALSE)
  }
  
  data$x <- factor(data$x, levels = c("R", "S", "I"), ordered = TRUE)
  
  ymax <- pm_if_else(max(data$s) > 95, 105, 100)
  
  plot(x = data$x,
       y = data$s,
       lwd = 2,
       ylim = c(0, ymax),
       ylab = ylab,
       xlab = xlab,
       main = main,
       axes = FALSE)
  # x axis
  axis(side = 1, at = 1:pm_n_distinct(data$x), labels = levels(data$x), lwd = 0)
  # y axis, 0-100%
  axis(side = 2, at = seq(0, 100, 5))
  
  text(x = data$x,
       y = data$s + 4,
       labels = paste0(data$s, "% (n = ", data$n, ")"))
}


#' @method barplot rsi
#' @importFrom graphics barplot axis
#' @export
#' @noRd
barplot.rsi <- function(height,
                        main = paste("Resistance Overview of", deparse(substitute(height))),
                        xlab = "Antimicrobial Interpretation",
                        ylab = "Frequency",
                        colours_RSI = c("#ED553B", "#3CAEA3", "#F6D55C"),
                        expand = TRUE,
                        ...) {
  meet_criteria(xlab, allow_class = "character", has_length = 1)
  meet_criteria(main, allow_class = "character", has_length = 1)
  meet_criteria(ylab, allow_class = "character", has_length = 1)
  meet_criteria(colours_RSI, allow_class = "character", has_length = c(1, 3))
  if (length(colours_RSI) == 1) {
    colours_RSI <- rep(colours_RSI, 3)
  }
  main <- gsub(" +", " ", paste0(main, collapse = " "))
  
  x <- table(height)
  x <- x[c(3, 1, 2)]
  barplot(x,
          col = colours_RSI,
          xlab = xlab,
          main = main,
          ylab = ylab,
          axes = FALSE)
  axis(2, seq(0, max(x)))
}
