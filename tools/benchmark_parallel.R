# Benchmark: sequential vs parallel as.sir() across data-set sizes
#
# Run from the repo root with:
#   Rscript tools/benchmark_parallel.R
# or from inside an R session:
#   source("tools/benchmark_parallel.R")
#
# Requires ggplot2 for the output plot; uses devtools::load_all() so the
# package does not need to be installed.

devtools::load_all(".", quiet = TRUE)

sizes <- c(20, 200, 2000, 20000)
n_ab  <- 6   # number of antibiotic columns

make_df <- function(n) {
  set.seed(42)
  mics <- lapply(seq_len(n_ab), function(j) {
    as.mic(sample(c("0.25", "0.5", "1", "2", "4", "8", "16", "32"), n, TRUE))
  })
  names(mics) <- c("AMC", "GEN", "CIP", "TZP", "IPM", "MEM")
  data.frame(mo = "B_ESCHR_COLI", mics, stringsAsFactors = FALSE)
}

results <- do.call(rbind, lapply(sizes, function(n) {
  df <- make_df(n)

  t_seq <- system.time(
    suppressMessages(as.sir(df, col_mo = "mo", info = FALSE, parallel = FALSE))
  )[["elapsed"]]

  t_par <- system.time(
    suppressMessages(as.sir(df, col_mo = "mo", info = FALSE, parallel = TRUE))
  )[["elapsed"]]

  message(sprintf("n = %6d  seq = %.3fs  par = %.3fs  speedup = %.1fx",
                  n, t_seq, t_par, t_seq / t_par))

  data.frame(n = n, mode = c("sequential", "parallel"),
             seconds = c(t_seq, t_par))
}))

if (requireNamespace("ggplot2", quietly = TRUE)) {
  p <- ggplot2::ggplot(results, ggplot2::aes(x = n, y = seconds,
                                             colour = mode, group = mode)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point(size = 3) +
    ggplot2::scale_x_log10(
      breaks = sizes,
      labels = format(sizes, big.mark = ",", scientific = FALSE)
    ) +
    ggplot2::scale_colour_manual(
      values = c(sequential = "#E05C5C", parallel = "#2E86AB")
    ) +
    ggplot2::labs(
      title    = "as.sir() throughput: sequential vs parallel",
      subtitle = sprintf("%d antibiotic columns, E. coli, EUCAST 2025", n_ab),
      x        = "Number of rows (log scale)",
      y        = "Wall-clock time (seconds)",
      colour   = NULL
    ) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(legend.position = "top")

  out_file <- "tools/benchmark_parallel.png"
  ggplot2::ggsave(out_file, p, width = 7, height = 5, dpi = 150)
  message("Plot saved to ", out_file)
} else {
  message("Install ggplot2 to get a plot; raw results:")
  print(results)
}
