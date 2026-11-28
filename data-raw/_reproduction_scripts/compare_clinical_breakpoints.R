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

# Compares the results of as.sir() between two builds of clinical_breakpoints, so that every changed interpretation is seen and explained, instead of only the changed rows.
#
# Each step runs from the command line, and a snapshot is made with the code and the data of the package folder that it is given, so that changes to as.sir() are compared as well:
#
#   Rscript compare_clinical_breakpoints.R grid     <grid.rds> <pkg_folder_old> <pkg_folder_new>
#   Rscript compare_clinical_breakpoints.R snapshot <grid.rds> <pkg_folder> <snapshot.rds> [cores]
#   Rscript compare_clinical_breakpoints.R diff     <old.rds> <new.rds> <changes.csv> [explain.R]
#
# The file explain.R must define explain_change(), which takes the data.frame of changes and returns one reason per row (NA if unexplained). The diff step exits with status 1 if any change remains unexplained.

suppressPackageStartupMessages(library(dplyr))

# Species that test the fallback to a broader breakpoint, on top of the organisms that are in the table for an agent and method. These cover the exceptions and extensions that EUCAST and CLSI state in their tables.
sentinels <- c(
  "Escherichia coli", "Klebsiella pneumoniae", "Klebsiella aerogenes", "Proteus mirabilis", "Morganella", "Morganella morganii", "Salmonella enterica", "Shigella flexneri", "Enterobacter cloacae", "Serratia marcescens", "Yersinia enterocolitica",
  "Pseudomonas aeruginosa", "Pseudomonas putida", "Acinetobacter baumannii", "Stenotrophomonas maltophilia", "Burkholderia pseudomallei", "Burkholderia cepacia",
  "Staphylococcus aureus", "Staphylococcus epidermidis", "Staphylococcus saprophyticus", "Staphylococcus lugdunensis", "Enterococcus faecalis", "Enterococcus faecium",
  "Streptococcus pneumoniae", "Streptococcus pyogenes", "Streptococcus agalactiae", "Streptococcus mitis", "Haemophilus influenzae", "Haemophilus parainfluenzae", "Moraxella catarrhalis", "Neisseria gonorrhoeae", "Neisseria meningitidis",
  "Listeria monocytogenes", "Campylobacter jejuni", "Pasteurella multocida", "Vibrio cholerae", "Aeromonas hydrophila", "Kingella kingae", "Corynebacterium diphtheriae", "Nocardia asteroides", "Bacillus anthracis", "Bacillus cereus", "Brucella melitensis",
  "Bacteroides fragilis", "Bacteroides", "Parabacteroides distasonis", "Phocaeicola vulgatus", "Prevotella bivia", "Clostridioides difficile", "Clostridium perfringens", "Cutibacterium acnes",
  "Candida albicans", "Aspergillus fumigatus"
)

load_breakpoints <- function(pkg) {
  env <- new.env()
  load(file.path(pkg, "data", "clinical_breakpoints.rda"), envir = env)
  as.data.frame(env$clinical_breakpoints, stringsAsFactors = FALSE)
}

# One line per guideline, agent and method, with the organisms and the values to test: every breakpoint, one step to either side of it, and both ends of the scale
make_grid <- function(pkg_old, pkg_new) {
  suppressMessages(devtools::load_all(pkg_new, quiet = TRUE))
  sentinel_mo <- as.character(as.mo(sentinels, info = FALSE))
  bind_rows(load_breakpoints(pkg_old), load_breakpoints(pkg_new)) |>
    filter(type == "human") |>
    mutate(across(c(guideline, method, ab, mo, site), as.character)) |>
    group_by(guideline, ab, method) |>
    summarise(
      has_screening = any(grepl("screen", site, ignore.case = TRUE)),
      mo = list(unique(c(setdiff(mo, "UNKNOWN"), sentinel_mo))),
      breakpoints = list(unique(c(as.double(breakpoint_S), as.double(breakpoint_R)))),
      .groups = "drop"
    ) |>
    mutate(value = Map(test_values, breakpoints, method)) |>
    select(-breakpoints)
}

test_values <- function(breakpoints, method) {
  breakpoints <- breakpoints[!is.na(breakpoints)]
  if (method == "MIC") {
    sort(unique(c(0.001, 512, breakpoints / 2, breakpoints, breakpoints * 2)))
  } else {
    sort(unique(c(6, 50, breakpoints - 1, breakpoints, breakpoints + 1)))
  }
}

snapshot_key <- function(key) {
  screening <- if (key$has_screening) c(FALSE, TRUE) else FALSE
  combinations <- expand.grid(
    mo = key$mo[[1]], value = key$value[[1]], uti = c(NA, FALSE, TRUE), include_screening = screening,
    stringsAsFactors = FALSE
  )
  out <- lapply(split(combinations, paste(combinations$uti, combinations$include_screening)), function(part) {
    uti <- part$uti[1]
    # keep the interpretation log per row, it shows which breakpoint was used
    invisible(sir_interpretation_history(clean = TRUE))
    x <- if (key$method == "MIC") as.mic(part$value) else as.disk(part$value)
    part$sir <- as.character(suppressWarnings(suppressMessages(as.sir(
      x,
      mo = part$mo, ab = key$ab, guideline = key$guideline, uti = uti, include_screening = part$include_screening[1], info = FALSE
    ))))
    used <- as.data.frame(sir_interpretation_history(clean = TRUE), stringsAsFactors = FALSE)
    part$breakpoint_used <- used$breakpoint_S_R[match(seq_len(nrow(part)), used$index)]
    part$ref_table <- used$ref_table[match(seq_len(nrow(part)), used$index)]
    part
  })
  bind_rows(out) |> mutate(guideline = key$guideline, ab = key$ab, method = key$method)
}

make_snapshot <- function(grid, pkg, cores = 1L) {
  suppressMessages(devtools::load_all(pkg, quiet = TRUE))
  keys <- split(grid, seq_len(nrow(grid)))
  results <- parallel::mclapply(keys, function(key) {
    tryCatch(snapshot_key(key), error = function(e) data.frame(guideline = key$guideline, ab = key$ab, method = key$method, error = conditionMessage(e)))
  }, mc.cores = cores, mc.preschedule = FALSE)
  bind_rows(results)
}

compare_snapshots <- function(old, new) {
  by <- c("guideline", "ab", "method", "mo", "value", "uti", "include_screening")
  full_join(old, new, by = by, suffix = c("_old", "_new")) |>
    filter(
      paste(sir_old) != paste(sir_new) | paste(breakpoint_used_old) != paste(breakpoint_used_new)
    )
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  if (args[1] == "grid") {
    saveRDS(make_grid(args[3], args[4]), args[2])
  } else if (args[1] == "snapshot") {
    cores <- if (length(args) >= 5) as.integer(args[5]) else 1L
    snapshot <- make_snapshot(readRDS(args[2]), args[3], cores)
    if ("error" %in% colnames(snapshot) && any(!is.na(snapshot$error))) {
      print(snapshot |> filter(!is.na(error)) |> distinct(guideline, ab, method, error) |> head(50))
      stop("the snapshot has errors, see above")
    }
    saveRDS(snapshot, args[4])
    message("Snapshot of ", nrow(snapshot), " interpretations saved")
  } else if (args[1] == "diff") {
    old <- readRDS(args[2])
    new <- readRDS(args[3])
    message("Compared: ", nrow(old), " and ", nrow(new), " interpretations")
    changes <- compare_snapshots(old, new)
    changes$reason <- NA_character_
    if (length(args) >= 5) {
      source(args[5], local = TRUE)
      changes$reason <- explain_change(changes)
    }
    write.csv(changes, args[4], row.names = FALSE, na = "")
    print(as.data.frame(changes |> count(reason, guideline)), row.names = FALSE)
    message(nrow(changes), " changes, of which ", sum(is.na(changes$reason)), " unexplained")
    if (anyNA(changes$reason)) {
      quit(save = "no", status = 1)
    }
  } else {
    stop("unknown step: ", args[1])
  }
}
