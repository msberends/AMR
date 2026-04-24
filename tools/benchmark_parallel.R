# Benchmark: sequential vs parallel as.sir() across data-set shapes
#
# Run from the repo root:
#   Rscript tools/benchmark_parallel.R
# or inside an R session:
#   source("tools/benchmark_parallel.R")
#
# Two panels:
#   Left  – fixed columns (n_ab_fixed), varying rows.
#           Parallel wins at small n; sequential catches up at large n due to
#           memory-bandwidth saturation (all workers compete for the same
#           clinical_breakpoints lookup table in L3 cache / RAM).
#   Right – fixed rows (n_rows_fixed), varying column count.
#           This is the shape that actually benefits: each additional column
#           keeps another core busy.  The "real world" gain for a 2854×65
#           dataset lives here.
#
# Requires ggplot2; uses devtools::load_all() so the package need not be
# installed.

devtools::load_all(".", quiet = TRUE)

# ── configuration ─────────────────────────────────────────────────────────────
row_sizes   <- c(200, 1000, 5000, 20000)
col_sizes   <- c(4, 8, 16, 32, 48)
n_rows_fixed <- 1000
n_ab_fixed   <- 16
n_cores_avail <- AMR:::get_n_cores(Inf)

all_abs <- c("AMC", "GEN", "CIP", "TZP", "IPM", "MEM",
             "AMP", "TMP", "SXT", "NIT", "FOX", "CRO",
             "FEP", "CAZ", "CTX", "TOB", "AMK", "ERY",
             "AZM", "CLI", "VAN", "TEC", "RIF", "MTR",
             "MFX", "LNZ", "TGC", "DOX", "FLC", "OXA",
             "PEN", "CXM", "CZO", "KAN", "COL", "FOS",
             "MUP", "TCY", "TEC", "IPM", "CHL", "FEP",
             "MEM", "TZP", "GEN", "AMC", "AMX", "AMP")
all_abs <- unique(all_abs)

mic_vals <- c("0.25", "0.5", "1", "2", "4", "8", "16", "32")

make_df <- function(n_rows, n_ab) {
  set.seed(42)
  ab_sel <- all_abs[seq_len(min(n_ab, length(all_abs)))]
  mics <- lapply(ab_sel, function(a) as.mic(sample(mic_vals, n_rows, TRUE)))
  names(mics) <- ab_sel
  data.frame(mo = "B_ESCHR_COLI", mics, stringsAsFactors = FALSE)
}

time_both <- function(n_rows, n_ab, label) {
  df <- make_df(n_rows, n_ab)
  t_seq <- system.time(
    suppressMessages(as.sir(df, col_mo = "mo", info = FALSE, parallel = FALSE))
  )[["elapsed"]]
  t_par <- system.time(
    suppressMessages(as.sir(df, col_mo = "mo", info = FALSE, parallel = TRUE))
  )[["elapsed"]]
  message(sprintf("%-28s  seq=%5.2fs  par=%5.2fs  speedup=%.1fx",
                  label, t_seq, t_par, t_seq / t_par))
  data.frame(group = label, mode = c("sequential", "parallel"),
             seconds = c(t_seq, t_par), stringsAsFactors = FALSE)
}

# ── warm-up (avoid first-call overhead biasing results) ───────────────────────
message("Warming up cache ...")
invisible(suppressMessages(as.sir(make_df(100, 6), col_mo = "mo", info = FALSE)))
invisible(suppressMessages(as.sir(make_df(100, 6), col_mo = "mo", info = FALSE, parallel = TRUE)))
sir_interpretation_history(clean = TRUE)

# ── panel 1: vary rows, fixed columns ─────────────────────────────────────────
message(sprintf("\nPanel 1 – varying rows, %d fixed columns:", n_ab_fixed))
res_rows <- do.call(rbind, lapply(row_sizes, function(n) {
  time_both(n, n_ab_fixed, sprintf("rows=%d", n))
}))
res_rows$x     <- rep(row_sizes, each = 2)
res_rows$panel <- "Vary rows (16 fixed AB columns)"

# ── panel 2: vary columns, fixed rows ─────────────────────────────────────────
message(sprintf("\nPanel 2 – varying columns, %d fixed rows:", n_rows_fixed))
res_cols <- do.call(rbind, lapply(col_sizes, function(n_ab) {
  time_both(n_rows_fixed, n_ab, sprintf("cols=%d", n_ab))
}))
res_cols$x     <- rep(col_sizes, each = 2)
res_cols$panel <- sprintf("Vary columns (%d fixed rows)", n_rows_fixed)

results <- rbind(res_rows, res_cols)

if (requireNamespace("ggplot2", quietly = TRUE)) {
  p <- ggplot2::ggplot(
    results,
    ggplot2::aes(x = x, y = seconds, colour = mode, group = mode)
  ) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::facet_wrap(~panel, scales = "free_x") +
    ggplot2::scale_colour_manual(
      values = c(sequential = "#E05C5C", parallel = "#2E86AB")
    ) +
    ggplot2::labs(
      title    = "as.sir() throughput: sequential vs parallel",
      subtitle = sprintf("E. coli, EUCAST 2026, %d cores available", n_cores_avail),
      x        = "Dataset dimension (rows  ·left·  or columns  ·right·)",
      y        = "Wall-clock time (seconds)",
      colour   = NULL
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "top")

  out_file <- "tools/benchmark_parallel.png"
  ggplot2::ggsave(out_file, p, width = 10, height = 5, dpi = 150)
  message("\nPlot saved to ", out_file)
} else {
  message("Install ggplot2 to get a plot; raw results:")
  print(results[, c("panel", "group", "mode", "seconds")])
}
