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
# Center Groningen in the Netherlands, in collaboration with many      #
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
# how to conduct AMR data analysis: https://amr-for-r.org              #
# ==================================================================== #

# Runs reproduction_of_microorganisms.R non-interactively, in chunks, in one R process.
#
# The script itself stays one file (to be run line by line by a human as well). This runner cuts it into chunks
# at its checkpoints (the saveRDS() lines of data-raw/taxonomy*.rds), runs these one by one, logs the start,
# end and duration of each chunk in a progress file, and stops with the name of the failing chunk.
# It can resume from any checkpoint: it then runs the setup chunk, loads the saved results and the parts of
# earlier chunks that later chunks need (marked with '# @resume-block' in the script), and continues.
#
# Usage, from the root of the repository:
#   Rscript data-raw/_reproduction_scripts/run_microorganisms_build.R [options]
# Options:
#   --review-dir=DIR        folder for the CSV files of all reviews (git-tracked, for the human review)
#   --log-dir=DIR           folder for the progress file (default: a temporary folder)
#   --resume-from=NAME      continue after this checkpoint, e.g. --resume-from=taxonomy1b
#   --list                  only show the chunks and checkpoints, and stop

script_file <- "data-raw/_reproduction_scripts/reproduction_of_microorganisms.R"
if (!file.exists(script_file)) {
  stop("Run this from the root of the AMR repository", call. = FALSE)
}

# arguments
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  hit <- args[startsWith(args, paste0("--", name, "="))]
  if (length(hit) == 0) {
    return(default)
  }
  sub("^[^=]+=", "", hit[1])
}
review_dir <- get_arg("review-dir")
log_dir <- get_arg("log-dir", tempfile("microorganisms_build_"))
resume_from <- get_arg("resume-from")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
progress_file <- file.path(log_dir, "progress.tsv")

# cut the script into chunks ----

lines <- readLines(script_file, encoding = "UTF-8")
# the checkpoints in their order; taxonomy_lpsn_missing.rds is no checkpoint to resume from, as the rest of the
# LPSN chunk is needed anyway (and is fast thanks to the LPSN cache)
checkpoints <- c(
  "taxonomy_lpsn", "taxonomy_mycobank", "taxonomy_gbif",
  "taxonomy0", "taxonomy1", "taxonomy1b", "taxonomy1c",
  "taxonomy2", "taxonomy2b", "taxonomy2c", "taxonomy3", "taxonomy3b"
)
checkpoint_line <- vapply(checkpoints, function(cp) {
  hit <- which(startsWith(lines, "saveRDS(") & grepl(paste0("\"data-raw/", cp, ".rds\""), lines, fixed = TRUE))
  if (length(hit) != 1) {
    stop("Checkpoint ", cp, " must be saved exactly once in ", script_file, " (found ", length(hit), ")", call. = FALSE)
  }
  hit
}, integer(1))
if (is.unsorted(checkpoint_line)) {
  stop("The checkpoints in the script are not in the expected order", call. = FALSE)
}
setup_end <- which(startsWith(lines, "# Read LPSN data")) - 1
stopifnot(length(setup_end) == 1)
chunks <- data.frame(
  name = c("setup", checkpoints, "save_to_package"),
  from = c(1L, setup_end + 1L, checkpoint_line + 1L),
  to = c(setup_end, checkpoint_line, length(lines)),
  stringsAsFactors = FALSE
)

# the blocks that are needed again when resuming
block_start <- which(startsWith(lines, "# @resume-block "))
block_end <- which(lines == "# @end-resume-block")
stopifnot(length(block_start) == length(block_end), all(block_start < block_end))
blocks <- data.frame(
  name = sub("^# @resume-block ([^ ]+) .*", "\\1", lines[block_start]),
  after = sub(".* after=([^ ]+).*", "\\1", lines[block_start]),
  from = block_start,
  to = block_end,
  stringsAsFactors = FALSE
)
stopifnot(all(blocks$after %in% checkpoints))

message("Chunks:")
for (i in seq_len(nrow(chunks))) {
  message(sprintf("  %-17s lines %5d-%5d", chunks$name[i], chunks$from[i], chunks$to[i]))
}
message("Blocks needed when resuming:")
for (i in seq_len(nrow(blocks))) {
  message(sprintf("  %-17s lines %5d-%5d, when resuming from %s or later", blocks$name[i], blocks$from[i], blocks$to[i], blocks$after[i]))
}
if ("--list" %in% args) {
  quit(save = "no", status = 0)
}

# run ----

options(
  AMR_build_view = FALSE,
  AMR_build_print_rows = 200,
  warn = 1,
  width = 200
)
if (!is.null(review_dir)) {
  options(AMR_build_review_dir = review_dir)
}
if (!file.exists(progress_file)) {
  writeLines("chunk\tstarted\tfinished\tminutes\tstatus", progress_file)
}
log_progress <- function(chunk, started, status) {
  cat(
    paste(chunk, format(started, "%Y-%m-%d %H:%M:%S"), format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      round(as.numeric(difftime(Sys.time(), started, units = "mins")), 1), status,
      sep = "\t"
    ), "\n",
    file = progress_file, append = TRUE, sep = ""
  )
}
run_lines <- function(name, from, to) {
  message("\n", strrep("=", 100), "\n=== CHUNK ", name, " (lines ", from, "-", to, ") started at ", format(Sys.time()), "\n", strrep("=", 100))
  started <- Sys.time()
  code <- paste(lines[from:to], collapse = "\n")
  ok <- tryCatch(
    {
      source(exprs = parse(text = code, keep.source = FALSE), local = globalenv(), echo = TRUE, max.deparse.length = 300)
      TRUE
    },
    error = function(e) {
      message("\n!!! CHUNK ", name, " FAILED: ", conditionMessage(e))
      log_progress(name, started, paste("FAILED:", gsub("[\t\n]", " ", conditionMessage(e))))
      FALSE
    }
  )
  if (!ok) {
    last_done <- which(chunks$name == name) - 2
    if (length(last_done) == 1 && last_done >= 1) {
      message("\nAfter fixing the cause, resume with: --resume-from=", checkpoints[last_done])
    } else {
      message("\nAfter fixing the cause, run again from the start")
    }
    quit(save = "no", status = 1)
  }
  log_progress(name, started, "OK")
  message("=== CHUNK ", name, " finished in ", round(as.numeric(difftime(Sys.time(), started, units = "mins")), 1), " minutes")
}

run_lines("setup", chunks$from[1], chunks$to[1])

if (is.null(resume_from)) {
  first_chunk <- 2
} else {
  if (!resume_from %in% checkpoints) {
    stop("--resume-from must be one of: ", paste(checkpoints, collapse = ", "), call. = FALSE)
  }
  k <- which(checkpoints == resume_from)
  message("\nResuming after checkpoint ", resume_from)
  started <- Sys.time()
  # the source taxonomies that were saved up to this checkpoint
  for (src in intersect(checkpoints[seq_len(k)], c("taxonomy_lpsn", "taxonomy_mycobank", "taxonomy_gbif"))) {
    assign(src, readRDS(file.path("data-raw", paste0(src, ".rds"))), envir = globalenv())
  }
  # the parts of earlier chunks that later chunks need
  for (i in which(match(blocks$after, checkpoints) <= k)) {
    run_lines(paste0("resume block ", blocks$name[i]), blocks$from[i], blocks$to[i])
  }
  if (!resume_from %in% c("taxonomy_lpsn", "taxonomy_mycobank", "taxonomy_gbif")) {
    assign("taxonomy", readRDS(file.path("data-raw", paste0(resume_from, ".rds"))), envir = globalenv())
  }
  log_progress(paste("resume from", resume_from), started, "OK")
  first_chunk <- which(chunks$name == resume_from) + 1
}

for (i in seq(first_chunk, nrow(chunks))) {
  run_lines(chunks$name[i], chunks$from[i], chunks$to[i])
}
message("\nAll chunks finished. Progress file: ", progress_file)
