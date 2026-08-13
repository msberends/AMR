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
# how to conduct AMR data analysis: https://amr-for-r.org              #
# ==================================================================== #

# Generic parser for EUCAST Clinical Breakpoint Table xlsx files.
# Works for any organism sheet by auto-detecting:
#   - Header rows (via "MIC breakpoint" in col 2)
#   - Sub-header rows (via "S " pattern in col 2, one row below header)
#   - Organism names (row 1 for single-organism sheets, or discovered above
#     each header row for multi-organism sheets like Anaerobic bacteria)
#   - Data rows (between sub-header and next header/end)
#   - Notes (col 9, resolved per-row via superscript references)
#
# Uses tidyxl to parse rich text cells, separating base values from
# superscript footnote references that EUCAST uses for notes.
#
# Output mimics the AMR::clinical_breakpoints structure with added columns:
#   version, is_screening, note


library(tidyxl)
library(dplyr, warn.conflicts = FALSE)
library(purrr)
devtools::load_all()

# ==============================================================================
# 1. Rich text parser
# ==============================================================================
split_rich_text <- function(fmt_list) {
  map_dfr(fmt_list, function(fmt) {
    if (is.null(fmt) || nrow(fmt) == 0)
      return(tibble(base_text = NA_character_, note_super = NA_character_))
    is_super <- !is.na(fmt$vertAlign) & fmt$vertAlign == "superscript"
    tibble(
      base_text  = paste0(fmt$character[!is_super], collapse = ""),
      note_super = paste0(fmt$character[is_super],  collapse = "")
    )
  })
}

# ==============================================================================
# 2. Notes block parser
# ==============================================================================
parse_notes_block <- function(txt) {
  if (is.na(txt) || txt == "") return(list())
  
  txt <- gsub("\r\n", "\n", txt)
  txt <- gsub("\r", "\n", txt)
  txt <- gsub("\u00A0", " ", txt)
  
  find_keys <- function(pattern, txt) {
    m <- gregexpr(pattern, txt, perl = TRUE)[[1]]
    if (m[1] == -1)
      return(data.frame(pos = integer(0), key = character(0),
                        mlen = integer(0), stringsAsFactors = FALSE))
    lens <- attr(m, "match.length")
    keys <- gsub("[.\\s]+$", "", trimws(substring(txt, m, m + lens - 1)))
    data.frame(pos = as.integer(m), key = keys, mlen = as.integer(lens),
               stringsAsFactors = FALSE)
  }
  
  # Combined key: 1/A. or 5/A.
  r1 <- find_keys("(?:^|(?<=\\n)|(?<=\\.))\\s*\\d+/[A-Z]\\.\\s*", txt)
  # Numbered with dot: "1." "2."
  r2 <- find_keys("(?:^|(?<=\\n)|(?<=\\.))\\s*\\d+\\.\\s+", txt)
  # Lettered: "A." followed by uppercase (not species names like "C. difficile")
  r3 <- find_keys("(?:^|(?<=\\n)|(?<=\\.)|(?<=\\s))\\s*[A-Z]\\.\\s+(?=[A-Z])", txt)
  # Numbered without dot (e.g. C. difficile style): digit space uppercase
  r4 <- find_keys("(?:^|(?<=\\n))\\d+\\s+(?=[A-Z])", txt)
  
  all_keys <- rbind(r1, r2, r3, r4)
  if (nrow(all_keys) == 0) return(list())
  all_keys <- all_keys[order(all_keys$pos), ]
  
  # Remove overlapping matches (keep longer/more specific)
  keep <- rep(TRUE, nrow(all_keys))
  for (i in seq_len(nrow(all_keys))) {
    if (!keep[i]) next
    for (j in seq_len(nrow(all_keys))) {
      if (i == j || !keep[j]) next
      if (all_keys$pos[j] >= all_keys$pos[i] &&
          all_keys$pos[j] < all_keys$pos[i] + all_keys$mlen[i]) {
        if (nchar(all_keys$key[i]) >= nchar(all_keys$key[j])) {
          keep[j] <- FALSE
        } else {
          keep[i] <- FALSE
        }
      }
    }
  }
  all_keys <- all_keys[keep, ]
  
  # Build note lookup, allowing multiple bodies per key
  notes <- list()
  add_note <- function(key, body) {
    if (is.null(notes[[key]])) {
      notes[[key]] <<- body
    } else if (!body %in% notes[[key]]) {
      notes[[key]] <<- c(notes[[key]], body)
    }
  }
  
  for (i in seq_len(nrow(all_keys))) {
    body_start <- all_keys$pos[i] + all_keys$mlen[i]
    body_end <- if (i < nrow(all_keys)) all_keys$pos[i + 1] - 1 else nchar(txt)
    body <- trimws(substr(txt, body_start, body_end))
    key <- all_keys$key[i]
    
    if (grepl("/", key)) {
      subkeys <- unlist(strsplit(key, "/"))
      for (sk in subkeys) add_note(sk, body)
      add_note(key, body)
    } else {
      add_note(key, body)
    }
  }
  notes
}

# ==============================================================================
# 3. Note resolver
# ==============================================================================
resolve_notes <- function(notes_list, mic_super, disk_super) {
  nl <- notes_list
  if (is.null(nl) || length(nl) == 0) return(NA_character_)
  
  mic_refs <- character(0)
  disk_refs <- character(0)
  if (!is.na(mic_super) && mic_super != "")
    mic_refs <- trimws(unlist(strsplit(mic_super, ",")))
  if (!is.na(disk_super) && disk_super != "")
    disk_refs <- trimws(unlist(strsplit(disk_super, ",")))
  all_refs <- unique(c(mic_refs, disk_refs))
  if (length(all_refs) == 0) return(NA_character_)
  
  # Detect combined keys (e.g. "1/A")
  combined_keys <- grep("/", names(nl), value = TRUE)
  used_combined <- character(0)
  consumed_refs <- character(0)
  
  for (ck in combined_keys) {
    parts <- unlist(strsplit(ck, "/"))
    if (all(parts %in% all_refs)) {
      used_combined <- c(used_combined, ck)
      consumed_refs <- c(consumed_refs, parts)
    }
  }
  
  remaining_refs <- setdiff(all_refs, consumed_refs)
  
  parts_out <- character(0)
  bodies_seen <- character(0)
  
  # Combined notes
  for (ck in used_combined) {
    bodies <- nl[[ck]]
    if (!is.null(bodies)) {
      for (body in bodies) {
        if (!body %in% bodies_seen) {
          parts_out <- c(parts_out, paste0("[", ck, "] ", body))
          bodies_seen <- c(bodies_seen, body)
        }
      }
    }
  }
  
  # Remaining individual notes: prefer duplicated over missing
  for (ref in remaining_refs) {
    if (!is.null(nl[[ref]])) {
      for (body in nl[[ref]]) {
        # Include all notes, even if body was already seen under a combined key
        entry <- paste0("[", ref, "] ", body)
        if (!entry %in% parts_out) {
          parts_out <- c(parts_out, entry)
        }
      }
    }
  }
  
  if (length(parts_out) == 0) return(NA_character_)
  paste(parts_out, collapse = " | ")
}

# ==============================================================================
# 4. Cell helpers
# ==============================================================================
get_cell <- function(df, r, c) {
  hit <- df[df$row == r & df$col == c, ]
  if (nrow(hit) == 0) return(list(base = NA_character_, super = NA_character_))
  list(base = hit$base_text[1], super = hit$note_super[1])
}

parse_bp <- function(txt) {
  if (is.na(txt) || txt == "") return(NA_real_)
  t <- trimws(txt)
  # if (toupper(t) %in% c("NOTE", "IE", "IP", "NA", "-")) return(NA_real_)
  if (toupper(t) %in% c("NA", "-")) return(NA_real_)
  t <- gsub("^\\((.+)\\)$", "\\1", t)
  # suppressWarnings(as.numeric(t))
  t
}

is_screening_fn <- function(txt) {
  if (is.na(txt)) return(FALSE)
  grepl("^\\(.+\\)$", trimws(txt))
}

format_disk_dose <- function(dose) {
  if (is.na(dose) || dose == "") return(NA_character_)
  d <- trimws(dose)
  if (grepl("mcg|unit", d, ignore.case = TRUE)) return(d)
  d <- gsub("-", "/", d)
  paste0(d, " mcg")
}

# ==============================================================================
# 5. Detect version and guideline from workbook
# ==============================================================================
detect_version <- function(cells) {
  # Version string is typically in col 5+ of row 1 (or sometimes col 9)
  # e.g. "EUCAST Clinical Breakpoint Tables v. 13.1, valid from 2023-06-29"
  version_cells <- cells |>
    filter(row == 1, !is_blank, col >= 5) |>
    mutate(text = coalesce(character, ""))
  version_text <- paste(version_cells$text, collapse = " ")
  
  version <- regmatches(version_text,
                        regexpr("v\\.?\\s*[\\d.]+", version_text, perl = TRUE))
  version <- gsub("v\\.?\\s*", "", version)
  if (length(version) == 0) version <- NA_character_
  
  major <- suppressWarnings(as.integer(sub("\\..*", "", version)))
  guideline <- if (!is.na(major)) paste0("EUCAST ", major + 2010) else NA_character_
  
  list(version = version, guideline = guideline)
}

# ==============================================================================
# 6. Main parser: parse a single sheet
# ==============================================================================
parse_sheet <- function(xlsx_path, sheet_name) {
  cells <- xlsx_cells(xlsx_path, sheets = sheet_name)
  
  # --- Detect version ---
  ver <- detect_version(cells)
  
  # --- Parse rich text for cols A:I ---
  target <- cells |> filter(col >= 1, col <= 9, !is_blank)
  rich <- split_rich_text(target$character_formatted)
  
  parsed <- target |>
    bind_cols(rich) |>
    mutate(
      base_text = case_when(
        !is.na(base_text) & base_text != "" ~ trimws(base_text),
        !is.na(character)                   ~ trimws(character),
        !is.na(numeric)                     ~ as.character(numeric),
        TRUE                                ~ NA_character_
      ),
      note_super = if_else(is.na(note_super) | note_super == "",
                           NA_character_, note_super)
    ) |>
    select(row, col, base_text, note_super)
  
  # --- Detect header rows (col 2 contains "MIC breakpoint") ---
  header_rows <- parsed |>
    filter(col == 2, grepl("MIC breakpoint", base_text, fixed = TRUE)) |>
    pull(row) |>
    sort()
  
  if (length(header_rows) == 0) {
    message("  No MIC breakpoint headers found in sheet '", sheet_name, "', skipping")
    return(NULL)
  }
  
  # Sub-header rows: one row below each header (contains "S " pattern)
  sub_header_rows <- header_rows + 1
  
  # --- Detect sheet type: single-organism vs multi-organism ---
  # Multi-organism: col A at header rows contains "Antimicrobial agent"
  header_col_a <- parsed |>
    filter(row %in% header_rows, col == 1) |>
    pull(base_text)
  
  is_multi_organism <- all(grepl("Antimicrobial agent", header_col_a, fixed = TRUE))
  
  # --- Detect organisms and their table ranges ---
  max_data_row <- max(parsed$row)
  
  if (is_multi_organism) {
    # Multi-organism: discover organism names above each header row.
    # Organism name is a standalone row (col A only, no data in B:H)
    # found between the previous table's end and this header row.
    rows_with_b <- parsed |> filter(col == 2) |> pull(row)
    
    tables <- list()
    for (i in seq_along(header_rows)) {
      hr <- header_rows[i]
      shr <- sub_header_rows[i]
      
      # Search backwards from the header row for the organism name
      search_from <- if (i == 1) 1 else sub_header_rows[i - 1] + 1
      candidates <- parsed |>
        filter(col == 1, row >= search_from, row < hr,
               !row %in% rows_with_b) |>
        # Exclude known non-organism text patterns
        filter(!grepl("Antimicrobial|MIC determination|Disk diffusion|Breakpoints for|Expert Rules|For abbreviations|For species|Examples of|haze|Isolated|Ignore haemolysis|Numbered|Lettered|Medium:|Inoculum:|Incubation:|Reading:|Quality control:|See disk",
                      base_text))
      
      # Take the last candidate (closest to the header)
      if (nrow(candidates) > 0) {
        org_row <- candidates |> slice_max(row, n = 1)
        organism <- org_row$base_text
      } else {
        # Fallback: use sheet name
        organism <- sheet_name
      }
      
      # Data range: from sub_header + 1 to the row before the next organism
      # or end of data
      first_data <- shr + 1
      if (i < length(header_rows)) {
        # End before the next organism name (or header)
        last_data <- header_rows[i + 1] - 1
        # Walk back to find actual last data row
        while (last_data >= first_data) {
          has_data <- nrow(parsed |> filter(row == last_data, col %in% 2:8)) > 0
          if (has_data) break
          last_data <- last_data - 1
        }
      } else {
        last_data <- max_data_row
        while (last_data >= first_data) {
          has_data <- nrow(parsed |> filter(row == last_data, col %in% 2:8)) > 0
          if (has_data) break
          last_data <- last_data - 1
        }
      }
      
      # Notes cell: col 9, in the data range (usually at first_data, merged)
      note_cell <- parsed |> filter(col == 9, row >= first_data, row <= last_data) |>
        slice_min(row, n = 1)
      note_text <- if (nrow(note_cell) > 0) note_cell$base_text[1] else NA_character_
      notes_parsed <- parse_notes_block(note_text)
      
      tables[[i]] <- list(organism = organism, first_data = first_data,
                          last_data = last_data, notes = notes_parsed)
    }
  } else {
    # Single-organism: organism name from row 1, col A
    org_cell <- parsed |> filter(row == 1, col == 1)
    organism <- if (nrow(org_cell) > 0) org_cell$base_text[1] else sheet_name
    # Clean: remove trailing * or whitespace
    organism <- gsub("[*]+$", "", trimws(organism))
    
    # Notes: col 9, within each class section
    # For single-organism sheets, there is one notes block per class.
    # Notes cell is at the header row or first data row of each class.
    tables <- list()
    for (i in seq_along(header_rows)) {
      hr <- header_rows[i]
      shr <- sub_header_rows[i]
      first_data <- shr + 1
      
      if (i < length(header_rows)) {
        last_data <- header_rows[i + 1] - 1
        while (last_data >= first_data) {
          has_data <- nrow(parsed |> filter(row == last_data, col %in% 2:8)) > 0
          if (has_data) break
          last_data <- last_data - 1
        }
      } else {
        last_data <- max_data_row
        while (last_data >= first_data) {
          has_data <- nrow(parsed |> filter(row == last_data, col %in% 2:8)) > 0
          if (has_data) break
          last_data <- last_data - 1
        }
      }
      
      # Notes: look in col 9 from header row through end of section
      note_cells <- parsed |>
        filter(col == 9, row >= hr, row <= last_data,
               !grepl("^Notes", base_text))
      note_text <- if (nrow(note_cells) > 0) {
        paste(note_cells$base_text, collapse = "\n")
      } else {
        NA_character_
      }
      notes_parsed <- parse_notes_block(note_text)
      
      tables[[i]] <- list(organism = organism, first_data = first_data,
                          last_data = last_data, notes = notes_parsed)
    }
  }
  
  # --- Extract breakpoint rows ---
  results <- vector("list", 500)
  idx <- 0L
  
  for (tbl in tables) {
    org <- tbl$organism
    nl  <- tbl$notes
    
    for (r in tbl$first_data:tbl$last_data) {
      agent <- get_cell(parsed, r, 1)
      if (is.na(agent$base)) next
      # Skip category headers (rows that have col A text but no data in B:H)
      has_any_data <- any(!is.na(c(
        get_cell(parsed, r, 2)$base,
        get_cell(parsed, r, 3)$base,
        get_cell(parsed, r, 5)$base,
        get_cell(parsed, r, 6)$base,
        get_cell(parsed, r, 7)$base
      )))
      if (!has_any_data) next
      
      mic_s   <- get_cell(parsed, r, 2)
      mic_r   <- get_cell(parsed, r, 3)
      mic_atu <- get_cell(parsed, r, 4)
      disk_dose_cell <- get_cell(parsed, r, 5)
      disk_s  <- get_cell(parsed, r, 6)
      disk_r  <- get_cell(parsed, r, 7)
      disk_atu <- get_cell(parsed, r, 8)
      
      note_text <- resolve_notes(nl, mic_s$super, disk_s$super)
      
      # MIC row
      s_val <- parse_bp(mic_s$base)
      r_val <- parse_bp(mic_r$base)
      if (!is.na(s_val) || !is.na(r_val)) {
        idx <- idx + 1L
        results[[idx]] <- tibble(
          guideline    = ver$guideline,
          version      = ver$version,
          type         = "human",
          host         = "human",
          method       = "MIC",
          site         = NA_character_,
          mo           = org,
          rank_index   = NA_integer_,
          ab           = agent$base,
          ref_tbl      = org,
          disk_dose    = NA_character_,
          breakpoint_S = s_val,
          breakpoint_R = r_val,
          uti          = FALSE,
          is_SDD       = FALSE,
          is_screening = is_screening_fn(mic_s$base),
          note         = note_text
        )
      }
      
      # DISK row
      s_val <- parse_bp(disk_s$base)
      r_val <- parse_bp(disk_r$base)
      if (!is.na(s_val) || !is.na(r_val)) {
        idx <- idx + 1L
        results[[idx]] <- tibble(
          guideline    = ver$guideline,
          version      = ver$version,
          type         = "human",
          host         = "human",
          method       = "DISK",
          site         = NA_character_,
          mo           = org,
          rank_index   = NA_integer_,
          ab           = agent$base,
          ref_tbl      = org,
          disk_dose    = format_disk_dose(disk_dose_cell$base),
          breakpoint_S = s_val,
          breakpoint_R = r_val,
          uti          = FALSE,
          is_SDD       = FALSE,
          is_screening = is_screening_fn(disk_s$base),
          note         = note_text
        )
      }
    }
  }
  
  if (idx == 0) return(NULL)
  bind_rows(results[seq_len(idx)])
}

# ==============================================================================
# 7. Parse "Topical agents" sheet (transposed layout)
# ==============================================================================
parse_topical_sheet <- function(xlsx_path) {
  cells <- xlsx_cells(xlsx_path, sheets = "Topical agents")
  
  ver <- detect_version(cells)
  
  target <- cells |> filter(!is_blank)
  rich <- split_rich_text(target$character_formatted)
  parsed <- target |>
    bind_cols(rich) |>
    mutate(
      base_text = case_when(
        !is.na(base_text) & base_text != "" ~ trimws(base_text),
        !is.na(character)                   ~ trimws(character),
        !is.na(numeric)                     ~ as.character(numeric),
        TRUE                                ~ NA_character_
      ),
      note_super = if_else(is.na(note_super) | note_super == "",
                           NA_character_, note_super)
    ) |>
    select(row, col, base_text, note_super)
  
  # Antimicrobial names from row 6, cols 4:18
  ab_cols <- 4:18
  ab_info <- list()
  for (col_idx in ab_cols) {
    name_cell <- get_cell(parsed, 6, col_idx)
    dose_cell <- get_cell(parsed, 10, col_idx)
    if (is.na(name_cell$base)) next
    ab_info[[as.character(col_idx)]] <- list(
      name = gsub("\r\n", " ", name_cell$base),
      name_super = name_cell$super,
      dose = if (!is.na(dose_cell$base) && !dose_cell$base %in% c("-", "ND"))
        paste0(dose_cell$base, " mcg") else NA_character_
    )
  }
  
  # Organism rows: col 1 from row 11 onwards, excluding "Notes"
  org_cells <- parsed |>
    filter(col == 1, row >= 11, !grepl("^Notes", base_text))
  
  organism_pairs <- list()
  for (i in seq_len(nrow(org_cells))) {
    mic_row <- org_cells$row[i]
    disk_row <- mic_row + 1
    mic_check <- get_cell(parsed, mic_row, 2)
    disk_check <- get_cell(parsed, disk_row, 2)
    if (is.na(mic_check$base) || !grepl("MIC", mic_check$base)) next
    if (is.na(disk_check$base) || !grepl("Zone", disk_check$base)) next
    organism_pairs[[length(organism_pairs) + 1]] <- list(
      organism = gsub("\r\n", " ", org_cells$base_text[i]),
      mic_row = mic_row, disk_row = disk_row
    )
  }
  
  # Notes: find in col 1, rows after the last organism data row
  last_data_row <- if (length(organism_pairs) > 0)
    max(vapply(organism_pairs, function(x) x$disk_row, numeric(1))) else 0L
  notes_cells <- parsed |>
    filter(col == 1, row > last_data_row, !grepl("^Notes$", base_text))
  notes_text <- if (nrow(notes_cells) > 0)
    paste(notes_cells$base_text, collapse = "\n") else ""
  if (is.na(notes_text)) notes_text <- ""
  notes_lookup <- parse_notes_block(notes_text)
  
  results <- vector("list", 200)
  idx <- 0L
  
  for (od in organism_pairs) {
    for (col_idx in ab_cols) {
      cc <- as.character(col_idx)
      if (is.null(ab_info[[cc]])) next
      ab_entry <- ab_info[[cc]]
      ab_name <- ab_entry$name
      ab_super <- ab_entry$name_super
      ab_dose <- ab_entry$dose
      
      # MIC
      mic_cell <- get_cell(parsed, od$mic_row, col_idx)
      mic_val <- parse_bp(mic_cell$base)
      # Combine superscripts from the cell and the agent name
      all_supers <- c(mic_cell$super, ab_super)
      all_supers <- all_supers[!is.na(all_supers) & all_supers != ""]
      mic_note <- resolve_notes(notes_lookup, paste(all_supers, collapse = ","), NA)
      
      if (!is.na(mic_val)) {
        idx <- idx + 1L
        results[[idx]] <- tibble(
          guideline = ver$guideline, version = ver$version,
          type = "human", host = "human", method = "MIC",
          site = "Topical", mo = od$organism, rank_index = NA_integer_,
          ab = ab_name, ref_tbl = "Topical agents",
          disk_dose = NA_character_,
          breakpoint_S = mic_val, breakpoint_R = mic_val,
          uti = FALSE, is_SDD = FALSE, is_screening = TRUE,
          note = mic_note
        )
      }
      
      # DISK
      disk_cell <- get_cell(parsed, od$disk_row, col_idx)
      disk_val <- parse_bp(disk_cell$base)
      all_supers_d <- c(disk_cell$super, ab_super)
      all_supers_d <- all_supers_d[!is.na(all_supers_d) & all_supers_d != ""]
      disk_note <- resolve_notes(notes_lookup, NA, paste(all_supers_d, collapse = ","))
      
      if (!is.na(disk_val)) {
        idx <- idx + 1L
        results[[idx]] <- tibble(
          guideline = ver$guideline, version = ver$version,
          type = "human", host = "human", method = "DISK",
          site = "Topical", mo = od$organism, rank_index = NA_integer_,
          ab = ab_name, ref_tbl = "Topical agents",
          disk_dose = ab_dose,
          breakpoint_S = disk_val, breakpoint_R = disk_val,
          uti = FALSE, is_SDD = FALSE, is_screening = TRUE,
          note = disk_note
        )
      }
    }
  }
  
  if (idx == 0) return(NULL)
  bind_rows(results[seq_len(idx)])
}

# ==============================================================================
# 8. Parse all data sheets
# ==============================================================================
parse_workbook <- function(xlsx_path,
                           skip_sheets = c("Content", "Changes", "Notes",
                                           "Guidance", "Dosages",
                                           "Technical uncertainty",
                                           "Non-species related breakpoints",
                                           "PKPD breakpoints",
                                           "PK PD breakpoints")) {
  all_sheets <- xlsx_sheet_names(xlsx_path)
  data_sheets <- setdiff(all_sheets, skip_sheets)
  
  results <- list()
  for (s in data_sheets) {
    message("Parsing ", basename(xlsx_path), ": ", s)
    res <- tryCatch({
      if (s == "Topical agents") {
        parse_topical_sheet(xlsx_path)
      } else {
        parse_sheet(xlsx_path, s)
      }
    },
    error = function(e) {
      message("  ERROR: ", conditionMessage(e))
      NULL
    })
    if (!is.null(res) && nrow(res) > 0) {
      res$sheet <- s
      results[[s]] <- res
      message("  -> ", nrow(res), " rows")
    }
  }
  
  bind_rows(results)
}

# ==============================================================================
# 8. Run
# ==============================================================================

breakpoint_files <- list.files(path = "data-raw",
                               pattern = "breakpoint_table.*xlsx",
                               full.names = TRUE,
                               recursive = FALSE,
                               ignore.case = TRUE)
breakpoint_files <- breakpoint_files[breakpoint_files %unlike% "dosages"]

breakpoints_eucast <- tibble()

for (xlsx_path in breakpoint_files) {
  
  message("Parsing ", basename(xlsx_path))
  result <- suppressMessages(parse_workbook(xlsx_path))
  
  # cat("\n============================\n")
  cat("Total rows:", nrow(result), "\n\n")
  # glimpse(result)
  
  # cat("\n--- Rows per sheet ---\n")
  # result |> count(sheet) |> print(n = 40)
  
  # cat("\n--- Rows per organism (top 20) ---\n")
  # result |> count(mo) |> arrange(desc(n)) |> head(20) |> print()
  
  breakpoints_eucast <- breakpoints_eucast |>
    bind_rows(result)
  
  # saveRDS(result, "/home/claude/eucast_all_breakpoints.rds")
  # write.csv(result, "/home/claude/eucast_all_breakpoints.csv", row.names = FALSE)
  # cat("\nSaved.\n")
}

breakpoints_eucast_raw <- breakpoints_eucast

saveRDS(breakpoints_eucast_raw, "data-raw/breakpoints_eucast_raw.rds")
write.csv(breakpoints_eucast_raw, "data-raw/breakpoints_eucast_raw.csv", row.names = FALSE)

breakpoints_eucast <- breakpoints_eucast_raw |>
  filter(!breakpoint_S %in% c("NOTE", "IE", "IP", "NA", "-"),
         !breakpoint_R %in% c("NOTE", "IE", "IP", "NA", "-")) |>
  mutate(breakpoint_S = as.numeric(breakpoint_S),
         breakpoint_R = as.numeric(breakpoint_R),
         uti = uti | ab %like_case% "UTI") |>
  filter(!(is.na(breakpoint_S) & is.na(breakpoint_R)))

breakpoints_eucast %>% count(guideline, version)
