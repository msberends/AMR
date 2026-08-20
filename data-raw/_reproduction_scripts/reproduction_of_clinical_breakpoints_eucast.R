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
# Output mimics the AMR::clinical_breakpoints structure with an added
# `version` column (guideline already covers the year; version keeps the
# exact x.y string), plus `sheet` and `note`.


library(tidyxl)
library(dplyr, warn.conflicts = FALSE)
library(purrr)
devtools::load_all()

# ==============================================================================
# 1. Rich text parser
# ==============================================================================
# Struck-through runs (tidyxl: character_formatted$strike == TRUE) are always
# disregarded before splitting base text from superscript note references,
# since EUCAST uses strike-through to mark removed/superseded text within a
# cell (most commonly within Notes blocks) that must not end up in the parsed
# output.
split_rich_text <- function(fmt_list) {
  map_dfr(fmt_list, function(fmt) {
    if (is.null(fmt) || nrow(fmt) == 0)
      return(tibble(base_text = NA_character_, note_super = NA_character_))
    is_struck <- !is.na(fmt$strike) & fmt$strike
    fmt <- fmt[!is_struck, , drop = FALSE]
    if (nrow(fmt) == 0)
      return(tibble(base_text = NA_character_, note_super = NA_character_))
    is_super <- !is.na(fmt$vertAlign) & fmt$vertAlign == "superscript"
    tibble(
      base_text  = paste0(fmt$character[!is_super], collapse = ""),
      note_super = paste0(fmt$character[is_super],  collapse = "")
    )
  })
}

# ==============================================================================
# 1b. Agent-name / organism-restriction splitter
# ==============================================================================
# On some sheets, the agent-name cell (column A) additionally restricts the
# row's breakpoints to a subset of the sheet's overall organism scope, e.g.
# Staphylococcus sheet: "Amikacin, Coagulase-negative staphylococci" (the
# same antibiotic gets different breakpoints for S. aureus vs coagulase-
# negative staphylococci), or Enterobacterales sheet: "Cefuroxime iv,
# E. coli, Klebsiella spp. (except K. aerogenes), Raoultella spp. and
# P. mirabilis". Parsing this from prose would be fragile (parenthetical
# route/indication qualifiers like "(uncomplicated UTI only)" are
# syntactically indistinguishable from organism qualifiers in free text).
#
# EUCAST instead encodes this distinction through formatting: the
# antibiotic name (including any parenthetical route/indication qualifier)
# is always bold, matching the sheet's hyperlink-styled agent-name
# convention, while an appended organism restriction is always set in a
# run with bold explicitly turned off. Splitting on this formatting
# boundary is far more robust than any text-pattern approach, and matches
# how EUCAST's own spreadsheet authors distinguish the two.
#
# The split point is the earliest run after which every remaining
# non-superscript run is non-bold (this correctly keeps a bold parenthetical
# that appears mid-cell, e.g. "Meropenem(indications other than
# meningitis), Pseudomonas other than P. aeruginosa" is not split at the
# first non-bold run (an empty separator) but at the true boundary before
# "Pseudomonas"). Superscript note-reference runs (which are always bold,
# by EUCAST convention, regardless of position) are excluded before this
# check, since they are handled separately by split_rich_text() and are
# not part of either the agent name or the organism restriction.
#
# If a multi-run cell's bold pattern never settles into an all-non-bold
# tail (i.e. no valid split point exists), the format is unexpected for
# this convention and the function errors rather than guessing, per design:
# in practice every EUCAST agent-name cell observed follows this
# convention, so a failure here should surface a genuine anomaly in a
# future version of the table rather than be silently absorbed. This is
# why `valid_rows` restricts the check to the actual data rows of the
# sheet's organism table(s) (computed by the caller once the table
# boundaries are known): column A elsewhere on a sheet can contain prose
# (section introductions, cross-references to another sheet's table) that
# was never an antibiotic name and has no reason to follow this
# convention, so it must not be able to trip the error at all.
split_agent_organism <- function(cells, formats, sheet_name, valid_rows) {
  col_a <- cells |> filter(col == 1, !is_blank, row %in% valid_rows)
  if (nrow(col_a) == 0) {
    return(tibble(row = integer(0), ab_text = character(0), mo_override = character(0)))
  }
  base_bold <- formats$local$font$bold[col_a$local_format_id]
  
  out <- vector("list", nrow(col_a))
  for (i in seq_len(nrow(col_a))) {
    fmt <- col_a$character_formatted[[i]]
    row_i <- col_a$row[i]
    
    if (is.null(fmt) || nrow(fmt) == 0) {
      # Plain (non-rich-text) cell: nothing to split, whole text is the agent name
      out[[i]] <- tibble(row = row_i, ab_text = col_a$character[i], mo_override = NA_character_)
      next
    }
    
    # Drop struck-through and superscript runs first, exactly as split_rich_text() does
    is_struck <- !is.na(fmt$strike) & fmt$strike
    fmt <- fmt[!is_struck, , drop = FALSE]
    is_super <- !is.na(fmt$vertAlign) & fmt$vertAlign == "superscript"
    fmt <- fmt[!is_super, , drop = FALSE]
    
    if (nrow(fmt) <= 1) {
      txt <- if (nrow(fmt) == 1) fmt$character[1] else col_a$character[i]
      out[[i]] <- tibble(row = row_i, ab_text = txt, mo_override = NA_character_)
      next
    }
    
    bold_resolved <- ifelse(is.na(fmt$bold), base_bold[i], fmt$bold)
    n <- length(bold_resolved)
    # Find the earliest k (k < n) such that runs (k+1):n are all non-bold
    split_at <- NA_integer_
    for (k in seq_len(n - 1)) {
      if (all(!bold_resolved[(k + 1):n])) {
        split_at <- k
        break
      }
    }
    
    if (is.na(split_at)) {
      # Either the whole cell is bold (no restriction present -> fine, keep
      # as a single agent name) or the bold pattern never settles into an
      # all-non-bold tail, which is the genuinely unexpected case.
      if (all(bold_resolved) || !any(bold_resolved)) {
        out[[i]] <- tibble(row = row_i,
                           ab_text = paste0(fmt$character, collapse = ""),
                           mo_override = NA_character_)
        next
      }
      stop(
        "split_agent_organism(): unexpected bold pattern in column A, sheet '",
        sheet_name, "', row ", row_i, ": bold never settles to a non-bold tail. ",
        "Text: '", paste0(fmt$character, collapse = ""), "'. ",
        "This breaks the assumption that antibiotic names are bold and any ",
        "trailing organism restriction is not; inspect this cell manually."
      )
    }
    
    ab_text <- paste0(fmt$character[seq_len(split_at)], collapse = "")
    mo_override <- paste0(fmt$character[(split_at + 1):n], collapse = "")
    out[[i]] <- tibble(row = row_i, ab_text = ab_text, mo_override = mo_override)
  }
  
  bind_rows(out) |>
    mutate(
      ab_text = trimws(gsub("[\r\n]+", " ", ab_text)),
      ab_text = gsub(",\\s*$", "", ab_text),
      mo_override = trimws(gsub("[\r\n]+", " ", mo_override)),
      # The comma separating the agent name from the restriction is not
      # always formatted consistently: it is usually part of the bold
      # agent-name run (so ab_text's trailing-comma strip above handles
      # it), but is occasionally fused into the start of the non-bold
      # restriction run instead, e.g. Staphylococcus sheet "Daptomycin,
      # other staphylococci" -- leaving a leading ", " on mo_override that
      # would otherwise reach as.mo() and fail to resolve.
      mo_override = sub("^,\\s*", "", mo_override),
      mo_override = na_if(mo_override, "")
    )
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
  # The version string sits in row 1, in whichever column follows the
  # organism name in col A. Its column position is NOT fixed: it depends on
  # how many data columns the sheet has (e.g. col I for sheets with MIC and
  # disk columns, but col C or E for MIC-only sheets such as N. gonorrhoeae
  # or H. pylori, which have fewer columns). Detecting it by requiring
  # col >= 5 therefore misses those narrower sheets. Instead, take every
  # non-blank cell in row 1 other than col A (the organism name), which is
  # robust regardless of sheet width.
  # e.g. "EUCAST Clinical Breakpoint Tables v. 13.1, valid from 2023-06-29"
  version_cells <- cells |>
    filter(row == 1, !is_blank, col > 1) |>
    mutate(text = coalesce(character, ""))
  version_text <- paste(version_cells$text, collapse = " ")
  
  version <- regmatches(version_text,
                        regexpr("v\\.?\\s*[\\d.]+", version_text, perl = TRUE))
  version <- gsub("v\\.?\\s*", "", version)
  if (length(version) == 0) version <- NA_character_
  
  # The guideline year is read directly from the "valid from YYYY-MM-DD"
  # date in the version string, not derived from the version number: the
  # two are not reliably related. EUCAST's bacterial clinical breakpoint
  # tables happen to have incremented their major version roughly once a
  # year since inception, which coincidentally makes "major version + 2010"
  # look right for those files, but the antifungal (AFST) table keeps its
  # own, independent version sequence -- v12.1 of the AFST table is dated
  # "valid from 2026-04-10", i.e. EUCAST 2026, not the "2022" that
  # "12 + 2010" would wrongly produce. Every version string observed
  # carries exactly one 4-digit year in its "valid from" date, so a direct
  # regex on that is both simpler and correct for both table families.
  year <- regmatches(version_text, regexpr("(?<=valid from )\\d{4}", version_text, perl = TRUE))
  if (length(year) == 0) {
    # Fall back to any 4-digit year anywhere in the string, in case a
    # future version string omits the "valid from" wording but still
    # states a year somewhere.
    year <- regmatches(version_text, regexpr("\\d{4}", version_text))
  }
  guideline <- if (length(year) == 1) paste0("EUCAST ", year) else NA_character_
  
  # A human-readable description of the source file/table, for ref_tbl,
  # e.g. "Clinical Breakpoint Tables v. 16.1" or "Antifungal Clinical
  # Breakpoint Table v. 12.1" -- everything in the version string between
  # the leading "EUCAST " and the trailing ", valid from ..." date.
  file_desc <- sub("^EUCAST\\s+", "", version_text)
  file_desc <- trimws(sub(",?\\s*valid from.*$", "", file_desc, ignore.case = TRUE))
  if (!nzchar(file_desc)) file_desc <- NA_character_
  
  list(version = version, guideline = guideline, file_desc = file_desc)
}

# ==============================================================================
# 6. Main parser: parse a single sheet
# ==============================================================================
parse_sheet <- function(xlsx_path, sheet_name) {
  cells <- xlsx_cells(xlsx_path, sheets = sheet_name)
  
  # --- Detect version ---
  ver <- detect_version(cells)
  
  formats <- xlsx_formats(xlsx_path)
  
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
  
  # --- Detect the notes column ---
  # On sheets that report both MIC and disk breakpoints, columns are laid out
  # as: A agent | B-D MIC (S/R/ATU) | E-H disk (dose/S/R/ATU) | I notes.
  # On sheets that report MIC only (no disk diffusion method exists, e.g.
  # N. gonorrhoeae, N. meningitidis, H. pylori, M. tuberculosis), columns D-H
  # are absent and the layout collapses to: A agent | B-D MIC (S/R/ATU) |
  # E notes. The notes column therefore shifts from I to E. Rather than
  # hard-coding either position, detect it per sheet from the header row
  # itself: it is the right-most populated column on the header row (the
  # cell whose text starts with "Notes").
  header_row_cols <- parsed |>
    filter(row %in% header_rows, col > 3, grepl("^Notes", base_text)) |>
    pull(col)
  
  notes_col <- if (length(header_row_cols) > 0) {
    # Normally consistent across header rows within a sheet; take the modal
    # (most frequent) value defensively in case of stray formatting cells.
    as.integer(names(sort(table(header_row_cols), decreasing = TRUE))[1])
  } else {
    9L  # fallback to the conventional position if no "Notes" label is found
  }
  
  # Data columns run up to (but excluding) the notes column; disk columns
  # are only present when notes_col > 5 (i.e. the sheet is wider than a
  # MIC-only layout of A:E).
  has_disk_columns <- notes_col > 5
  
  # --- Detect sheet type: single-organism vs multi-organism ---
  # Multi-organism: col A at header rows contains "Antimicrobial agent".
  # This label alone is not sufficient, though: some genuinely single-
  # organism sheets (e.g. M. tuberculosis) also happen to use "Antimicrobial
  # agent" at their (single) header row, rather than an antibiotic class
  # name. The reliable distinguishing signal is the number of header rows:
  # a sheet with only one MIC breakpoint header is single-organism
  # regardless of that header's col-A label, since true multi-organism
  # sheets (Anaerobic bacteria, Topical agents) repeat the header once per
  # organism.
  header_col_a <- parsed |>
    filter(row %in% header_rows, col == 1) |>
    pull(base_text)
  
  is_multi_organism <- length(header_rows) > 1 &&
    all(grepl("Antimicrobial agent", header_col_a, fixed = TRUE))
  
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
          has_data <- nrow(parsed |> filter(row == last_data, col %in% setdiff(2:8, notes_col))) > 0
          if (has_data) break
          last_data <- last_data - 1
        }
      } else {
        last_data <- max_data_row
        while (last_data >= first_data) {
          has_data <- nrow(parsed |> filter(row == last_data, col %in% setdiff(2:8, notes_col))) > 0
          if (has_data) break
          last_data <- last_data - 1
        }
      }
      
      # Notes cell: notes_col, in the data range (usually at first_data, merged)
      note_cell <- parsed |> filter(col == notes_col, row >= first_data, row <= last_data) |>
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
    
    # Notes: notes_col, within each class section
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
          has_data <- nrow(parsed |> filter(row == last_data, col %in% setdiff(2:8, notes_col))) > 0
          if (has_data) break
          last_data <- last_data - 1
        }
      } else {
        last_data <- max_data_row
        while (last_data >= first_data) {
          has_data <- nrow(parsed |> filter(row == last_data, col %in% setdiff(2:8, notes_col))) > 0
          if (has_data) break
          last_data <- last_data - 1
        }
      }
      
      # Notes: look in notes_col from header row through end of section
      note_cells <- parsed |>
        filter(col == notes_col, row >= hr, row <= last_data,
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
  # Split agent name from any organism restriction in column A (see
  # split_agent_organism() for the formatting-based rationale). This is
  # scoped to the actual data rows of each organism table, computed just
  # above, rather than run across the whole column: prose cells elsewhere
  # in column A (section intros, cross-references to other sheets) are not
  # antibiotic names and are not expected to follow the bold-name /
  # non-bold-restriction convention, so they must not be able to trip the
  # "unexpected format" error at all.
  data_rows <- unique(unlist(map(tables, function(tbl) tbl$first_data:tbl$last_data)))
  agent_org_split <- split_agent_organism(cells, formats, sheet_name, valid_rows = data_rows)
  
  results <- vector("list", 500)
  idx <- 0L
  
  for (tbl in tables) {
    org <- tbl$organism
    nl  <- tbl$notes
    
    for (r in tbl$first_data:tbl$last_data) {
      agent <- get_cell(parsed, r, 1)
      if (is.na(agent$base)) next
      # Skip category headers (rows that have col A text but no data).
      # Columns 5:7 only hold disk data when the sheet actually has disk
      # columns; on MIC-only sheets col 5 is the notes column and must not
      # be treated as a data indicator.
      data_cols <- if (has_disk_columns) c(2, 3, 5, 6, 7) else c(2, 3)
      has_any_data <- any(!is.na(vapply(data_cols, function(cc) get_cell(parsed, r, cc)$base,
                                        character(1))))
      if (!has_any_data) next
      
      # Agent name and any organism restriction encoded in its formatting
      # (see split_agent_organism()). Falls back to the raw cell text and
      # the table-level organism when this row has no split entry (e.g.
      # plain-text cells with no rich formatting at all).
      split_row <- agent_org_split[agent_org_split$row == r, ]
      agent_text <- if (nrow(split_row) > 0) split_row$ab_text[1] else agent$base
      row_mo <- if (nrow(split_row) > 0 && !is.na(split_row$mo_override[1])) {
        split_row$mo_override[1]
      } else {
        org
      }
      
      mic_s   <- get_cell(parsed, r, 2)
      mic_r   <- get_cell(parsed, r, 3)
      mic_atu <- get_cell(parsed, r, 4)
      # Disk columns (E:H) only exist on sheets that report disk diffusion
      # breakpoints. On MIC-only sheets, col E is the notes column instead
      # of the disk dose, so disk cells must not be read from it there.
      if (has_disk_columns) {
        disk_dose_cell <- get_cell(parsed, r, 5)
        disk_s  <- get_cell(parsed, r, 6)
        disk_r  <- get_cell(parsed, r, 7)
        disk_atu <- get_cell(parsed, r, 8)
      } else {
        disk_dose_cell <- list(base = NA_character_, super = NA_character_)
        disk_s  <- list(base = NA_character_, super = NA_character_)
        disk_r  <- list(base = NA_character_, super = NA_character_)
        disk_atu <- list(base = NA_character_, super = NA_character_)
      }
      
      # Notes must apply only to the method they were referenced from: MIC
      # superscripts resolve to the MIC row's note, disk superscripts to the
      # DISK row's note. Resolving both against the union of MIC and disk
      # references (as before) duplicated every note into both rows.
      # A superscript on the agent name itself (col A) is method-agnostic
      # and must be folded into both sides, e.g. L. monocytogenes
      # "Trimethoprim-sulfamethoxazole (all indications)1".
      combine_super <- function(...) {
        parts <- c(...)
        parts <- parts[!is.na(parts) & parts != ""]
        if (length(parts) == 0) return(NA_character_)
        paste(parts, collapse = ",")
      }
      mic_note  <- resolve_notes(nl, combine_super(mic_s$super, agent$super), NA)
      disk_note <- resolve_notes(nl, NA, combine_super(disk_s$super, agent$super))
      
      # MIC row
      s_val <- parse_bp(mic_s$base)
      r_val <- parse_bp(mic_r$base)
      if (!is.na(s_val) || !is.na(r_val)) {
        idx <- idx + 1L
        results[[idx]] <- tibble(
          guideline    = ver$guideline,
          version      = ver$version,
          file_desc    = ver$file_desc,
          type         = "human",
          host         = "human",
          method       = "MIC",
          site         = NA_character_,
          mo           = row_mo,
          rank_index   = NA_integer_,
          ab           = agent_text,
          disk_dose    = NA_character_,
          breakpoint_S = s_val,
          breakpoint_R = r_val,
          note         = mic_note
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
          file_desc    = ver$file_desc,
          type         = "human",
          host         = "human",
          method       = "DISK",
          site         = NA_character_,
          mo           = row_mo,
          rank_index   = NA_integer_,
          ab           = agent_text,
          disk_dose    = format_disk_dose(disk_dose_cell$base),
          breakpoint_S = s_val,
          breakpoint_R = r_val,
          note         = disk_note
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
          guideline = ver$guideline, version = ver$version, file_desc = ver$file_desc,
          type = "human", host = "human", method = "MIC",
          site = "Topical", mo = od$organism, rank_index = NA_integer_,
          ab = ab_name,
          disk_dose = NA_character_,
          breakpoint_S = mic_val, breakpoint_R = mic_val,
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
          guideline = ver$guideline, version = ver$version, file_desc = ver$file_desc,
          type = "human", host = "human", method = "DISK",
          site = "Topical", mo = od$organism, rank_index = NA_integer_,
          ab = ab_name,
          disk_dose = ab_dose,
          breakpoint_S = disk_val, breakpoint_R = disk_val,
          note = disk_note
        )
      }
    }
  }
  
  if (idx == 0) return(NULL)
  bind_rows(results[seq_len(idx)])
}

# ==============================================================================
# 8. Parse antifungal sheets with a transposed species-in-columns layout
#    (EUCAST AFST tables: "Yeast", "Aspergillus")
# ==============================================================================
# Layout differs from the bacterial sheets in two ways:
#  - Species run across column-blocks rather than down rows. Each block is
#    2 columns wide (S/R) or 3 (S/R/ATU), detected dynamically from the
#    sub-header row rather than assumed fixed-width, since Aspergillus mixes
#    2- and 3-column species blocks on the same sheet.
#  - Notes sit in a block below the table, one full note per row (starting
#    "1. ", "2. " etc.) rather than a single merged cell to the right of the
#    table; these are concatenated and handed to the same parse_notes_block()
#    used elsewhere.
# Antifungal breakpoints on these sheets are MIC-only; there is no disk
# diffusion method, so only method = "MIC" rows are produced.
parse_yeast_sheet <- function(xlsx_path, sheet_name) {
  cells <- xlsx_cells(xlsx_path, sheets = sheet_name)
  
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
  
  # --- Locate the "Antifungal agent" header row ---
  hdr_row <- parsed |>
    filter(col == 1, grepl("^Antifungal agent$", base_text)) |>
    pull(row)
  if (length(hdr_row) == 0) {
    message("  No 'Antifungal agent' header found in sheet '", sheet_name, "', skipping")
    return(NULL)
  }
  hdr_row <- hdr_row[1]
  species_row <- hdr_row + 1L
  sub_row <- hdr_row + 2L
  first_data_row <- hdr_row + 3L
  
  # --- Detect species column-blocks from the sub-header row ---
  # Each block starts at a cell reading "S <=" (S \u2264) and runs through the
  # next "R >" and, where present, an immediately following "ATU" column,
  # stopping before the next "S <=" or the end of populated columns.
  sub_cells <- parsed |> filter(row == sub_row) |> arrange(col)
  s_cols <- sub_cells |> filter(grepl("^S\\b", base_text)) |> pull(col)
  
  if (length(s_cols) == 0) {
    message("  No species S/R columns found in sheet '", sheet_name, "', skipping")
    return(NULL)
  }
  
  species_blocks <- list()
  for (i in seq_along(s_cols)) {
    s_col <- s_cols[i]
    r_col <- s_col + 1L
    # ATU present if the cell two columns after S reads "ATU" (i.e. this
    # block is 3 columns wide rather than 2)
    atu_cell <- get_cell(parsed, sub_row, s_col + 2L)
    has_atu <- !is.na(atu_cell$base) && grepl("^ATU\\b", atu_cell$base)
    
    sp_cell <- get_cell(parsed, species_row, s_col)
    species <- if (!is.na(sp_cell$base)) gsub("\r?\n", " ", sp_cell$base) else NA_character_
    if (is.na(species)) next
    # Expand an abbreviated genus initial (e.g. "A. terreus") to the
    # sheet's own genus ("Aspergillus terreus") before this reaches
    # as.mo(): a bare genus initial is ambiguous across unrelated genera
    # with no context to disambiguate, and this is not merely theoretical
    # here -- as.mo("A. terreus") alone resolves to Acinetobacter
    # terrestris, not Aspergillus terreus, since Acinetobacter has a
    # matching species epithet and nothing here signals which genus was
    # meant. The sheet name is used as this genus reference; it consists
    # of only the bare genus for both AFST sheets this parser handles
    # ("Yeast" is actually genus-mixed and not abbreviated in practice, so
    # this only fires for "Aspergillus").
    sheet_genus <- sub("^.*?([A-Z][a-z]+).*$", "\\1", trimws(sheet_name))
    if (grepl("^[A-Z]\\.\\s", species) && substr(species, 1, 1) == substr(sheet_genus, 1, 1)) {
      species <- sub("^[A-Z]\\.", sheet_genus, species)
    }
    
    species_blocks[[length(species_blocks) + 1]] <- list(
      species = species, s_col = s_col, r_col = r_col,
      atu_col = if (has_atu) s_col + 2L else NA_integer_
    )
  }
  
  # --- Locate the end of the data block (row before the "Notes" label) ---
  notes_label_row <- parsed |>
    filter(col == 1, grepl("^Notes$", base_text), row > first_data_row) |>
    pull(row)
  last_data_row <- if (length(notes_label_row) > 0) {
    min(notes_label_row) - 1L
  } else {
    max(parsed$row)
  }
  
  # --- Notes block: one note per row below the "Notes" label, col A ---
  notes_text <- if (length(notes_label_row) > 0) {
    note_cells <- parsed |>
      filter(col == 1, row > min(notes_label_row)) |>
      arrange(row)
    paste(note_cells$base_text, collapse = "\n")
  } else {
    ""
  }
  notes_lookup <- parse_notes_block(notes_text)
  
  # --- Extract breakpoint rows: one per (agent, species) combination ---
  results <- vector("list", 1000)
  idx <- 0L
  
  for (r in first_data_row:last_data_row) {
    agent_cell <- get_cell(parsed, r, 1)
    if (is.na(agent_cell$base)) next
    agent_name <- gsub("\r?\n", " ", agent_cell$base)
    agent_super <- agent_cell$super
    
    for (sb in species_blocks) {
      s_cell <- get_cell(parsed, r, sb$s_col)
      r_cell <- get_cell(parsed, r, sb$r_col)
      
      s_val <- parse_bp(s_cell$base)
      r_val <- parse_bp(r_cell$base)
      if (is.na(s_val) && is.na(r_val)) next
      
      all_supers <- c(s_cell$super, r_cell$super, agent_super)
      all_supers <- all_supers[!is.na(all_supers) & all_supers != ""]
      note_text <- resolve_notes(notes_lookup, paste(all_supers, collapse = ","), NA)
      
      idx <- idx + 1L
      results[[idx]] <- tibble(
        guideline    = ver$guideline,
        version      = ver$version,
        file_desc    = ver$file_desc,
        type         = "human",
        host         = "human",
        method       = "MIC",
        site         = NA_character_,
        mo           = sb$species,
        rank_index   = NA_integer_,
        ab           = agent_name,
        disk_dose    = NA_character_,
        breakpoint_S = s_val,
        breakpoint_R = r_val,
        note         = note_text
      )
    }
  }
  
  if (idx == 0) return(NULL)
  bind_rows(results[seq_len(idx)])
}

# ==============================================================================
# 9. Parse all data sheets
# ==============================================================================
parse_workbook <- function(xlsx_path,
                           skip_sheets = c("Content", "Changes", "Notes",
                                           "Guidance", "Dosages",
                                           "Technical uncertainty",
                                           "Non-species related breakpoints",
                                           "PKPD breakpoints",
                                           "PK PD breakpoints")) {
  all_sheets <- xlsx_sheet_names(xlsx_path)
  
  # Some workbooks (e.g. the antifungal AFST tables) number-prefix their
  # sheet names, e.g. "1. Notes", "7. Dosages", "6. Aspergillus ", and are
  # not always consistent in capitalisation, e.g. "3. Technical Uncertainty"
  # vs "Technical uncertainty". Matching skip_sheets by exact equality would
  # fail to skip these, so match case-insensitively by whether a skip name
  # occurs anywhere in the (trimmed) sheet name instead.
  is_skipped <- vapply(trimws(all_sheets), function(s) {
    any(vapply(skip_sheets, function(skip) grepl(tolower(skip), tolower(s), fixed = TRUE),
               logical(1)))
  }, logical(1))
  data_sheets <- all_sheets[!is_skipped]
  
  results <- list()
  for (s in data_sheets) {
    message("Parsing ", basename(xlsx_path), ": ", s)
    res <- tryCatch({
      if (trimws(s) == "Topical agents") {
        parse_topical_sheet(xlsx_path)
      } else if (trimws(s) %in% c("5. Yeast", "6. Aspergillus")) {
        parse_yeast_sheet(xlsx_path, s)
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
# 10. Run
# ==============================================================================

breakpoint_files <- list.files(path = "data-raw",
                               pattern = "breakpoint.*table.*[.]xlsx$",
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


# Cleanup for `clinical_breakpoints` table ----
# This section derives the AMR::clinical_breakpoints-compatible
# `breakpoints_eucast` table from `breakpoints_eucast_raw`, following the
# same conventions as the WHONET-based `breakpoints_new` build (see
# data-raw/reproduction_of_clinical_breakpoints.R): guideline/type/host/
# method/mo/ab/rank_index/ref_tbl/disk_dose/breakpoint_S/breakpoint_R/uti/
# is_SDD, arranged and de-duplicated the same way.
#
# mo resolution needs pre-processing the WHONET pipeline doesn't, because
# EUCAST expresses organism scope as free text rather than as a single
# coded taxon per row (both at sheet level, e.g. "Bacillus spp. except
# B. anthracis", and at row level via split_agent_organism()'s
# mo_override, e.g. "other staphylococci" or "E. coli, Klebsiella spp.
# (except K. aerogenes), Raoultella spp. and P. mirabilis"):
#
#  1. An "except ..." / "other than ..." clause attached to a taxon does
#     NOT need expanding into every other member of that taxon: the
#     exception is already handled by rank_index precedence, since the
#     named exception (e.g. B. anthracis, K. aerogenes) gets its own more
#     specific row from its own sheet or its own list entry, which takes
#     priority over the broader row at lookup time. This mirrors how the
#     currently published clinical_breakpoints handles it (confirmed:
#     "Bacillus" stays at genus rank alongside separate species-level
#     "Bacillus anthracis" rows; same for Corynebacterium). The clause is
#     therefore just stripped before mo resolution, wherever it occurs.
#  2. Elliptical phrasing that refers back to the sheet's own organism
#     rather than naming a taxon ("other staphylococci", "other
#     enterococci", "other Enterobacterales") is resolved to the sheet
#     name itself, which already resolves correctly via as.mo().
#  3. Cells naming two or more taxa jointly, e.g. "Aerococcus sanguinicola
#     and A. urinae" or "E. coli, Klebsiella spp. (except K. aerogenes),
#     Raoultella spp. and P. mirabilis", are split into one row per named
#     taxon before resolution: as.mo() cannot reliably resolve a joint
#     string (it either mismatches to a single subspecies or fails
#     outright), whereas each individual member resolves correctly on its
#     own, matching how the currently published dataset represents such
#     cases (as separate rows, not a joint entry).
strip_mo_exceptions <- function(x) {
  trimws(sub("\\s*\\(?(except|other than)\\b.*$", "", x, ignore.case = TRUE, perl = TRUE))
}

# Whole-cell variant of the above, used before any list-splitting is
# attempted. Only strips a trailing except/other-than clause when it is NOT
# wrapped in parentheses attached to a single list member -- e.g.
# "Bacillus spp. except B. anthracis" and "Corynebacterium spp. other than
# C. diphtheriae and C. ulcerans" both run genuinely to the end of the cell
# with no further list content, so the whole tail is dropped. By contrast,
# "E. coli, Klebsiella spp. (except K. aerogenes), Raoultella spp. and
# P. mirabilis" has its exception parenthesized and attached to only one
# list member, with more members following afterward; stripping from
# "except" to the end of the whole string there would wrongly discard
# "Raoultella spp. and P. mirabilis" too. That case is intentionally left
# untouched here and handled per-member by strip_mo_exceptions() inside
# split_mo_list() instead, after the list has already been split apart.
strip_mo_exceptions_wholecell <- function(x) {
  trimws(sub("\\s*(except|other than)\\b(?![^(]*\\)).*$", "", x, ignore.case = TRUE, perl = TRUE))
}

# Extracts a single trailing parenthetical qualifier from an antibiotic-name
# string, e.g. "Cefuroxime oral (uncomplicated UTI only)" -> "uncomplicated
# UTI only". Returns NA where no such qualifier is present. This is captured
# into `site` before the qualifier is stripped from `ab`, so that route- or
# indication-specific variants of the same drug (e.g. "... iv" vs "... oral
# (uncomplicated UTI only)") remain distinguishable by the distinct() step
# below, rather than colliding when their S breakpoint happens to match.
extract_trailing_parenthetical <- function(x) {
  m <- regexpr("(?<=\\()[^)]*(?=\\)\\s*$)", x, perl = TRUE)
  out <- rep(NA_character_, length(x))
  found <- m != -1
  out[found] <- regmatches(x, m)
  out
}

resolve_other_x <- function(x, sheet_name) {
  # A few agent-name-embedded organism restrictions use elliptical phrasing
  # that refers back to the sheet's own organism rather than naming a taxon
  # directly, e.g. Staphylococcus sheet: "other staphylococci"; Enterococcus
  # sheet: "other enterococci"; Enterobacterales sheet: "other
  # Enterobacterales". In every observed case this means "the sheet's own
  # genus/order", which the sheet name itself already resolves to directly
  # (as.mo("Staphylococcus") -> genus Staphylococcus, etc.), so the fallback
  # is simply the sheet name, not a species enumeration.
  if_else(grepl("^other\\b", trimws(x), ignore.case = TRUE), sheet_name, x)
}

expand_abbreviated_genus <- function(pieces, fallback_genus = NA_character_) {
  # Tracks the most recently seen fully-spelled genus within a split list
  # and expands any subsequent "X. species" abbreviation sharing its first
  # letter to use it in full, e.g. in "Corynebacterium diphtheriae and
  # C. ulcerans", the "C." in the second piece is expanded to
  # "Corynebacterium" before as.mo() sees it. This matters because a bare
  # genus initial is often genuinely ambiguous across many genera -- e.g.
  # as.mo("C. ulcerans") alone resolves to Campylobacter upsaliensis, not
  # Corynebacterium ulcerans, since there is no context to prefer one
  # genus over another -- whereas the fully-qualified name is unambiguous.
  # A piece is only expanded if its abbreviation letter actually matches
  # the tracked genus's initial, so e.g. the leading "E. coli" in an
  # Enterobacterales list is correctly left alone (nothing precedes it to
  # expand from), and an abbreviation that happens to not match the most
  # recent genus is also left alone rather than expanded incorrectly.
  #
  # fallback_genus seeds this tracking before the list is scanned, for
  # sheets whose organism scope is itself a single genus (e.g. the
  # Staphylococcus sheet's "S. pseudintermedius, S. intermedius,
  # S. schleiferi and S. coagulans" never spells out "Staphylococcus" even
  # once, so without a seed nothing would ever be expanded, and
  # as.mo("S. intermedius") alone resolves to the wrong genus entirely --
  # Streptococcus intermedius, not Staphylococcus intermedius). Only pass
  # a fallback_genus when the sheet's own organism is confirmed to resolve
  # to genus rank (see split_mo_list()); an order-level sheet name like
  # "Enterobacterales" must never be used this way, since it shares a
  # first letter with genera it has nothing to do with (e.g. "E. coli").
  #
  # A bare lowercase species epithet with no genus marker at all (e.g.
  # "urinae" in "Aerococcus sanguinicola and urinae", an older EUCAST
  # phrasing -- see split_mo_list_one()) is likewise expanded against the
  # tracked genus, prefixed rather than substituted since there is no
  # abbreviation to replace.
  current_genus <- fallback_genus
  vapply(pieces, function(p) {
    full_match <- regmatches(p, regexpr("^[A-Z][a-z]+(?=\\s)", p, perl = TRUE))
    if (length(full_match) == 1 && nzchar(full_match)) {
      current_genus <<- full_match
      return(p)
    }
    if (!is.na(current_genus) && grepl("^[A-Z]\\.\\s", p) &&
        substr(p, 1, 1) == substr(current_genus, 1, 1)) {
      return(sub("^[A-Z]\\.", current_genus, p))
    }
    if (!is.na(current_genus) && grepl("^[a-z]+$", p)) {
      return(paste(current_genus, p))
    }
    p
  }, character(1), USE.NAMES = FALSE)
}

split_mo_list <- function(x, sheet_name = NA_character_) {
  # Splits a free-text organism cell into one entry per named taxon. Handles:
  #  - a single taxon (returned unchanged)
  #  - two taxa joined by "and", e.g. "Aerococcus sanguinicola and A. urinae"
  #  - N taxa in a comma/"and" list, e.g. "E. coli, Klebsiella spp. (except
  #    K. aerogenes), Raoultella spp. and P. mirabilis"
  # An "(except ...)" or "except ..." clause attached to one list member is
  # stripped from that member only (not the whole cell), since it scopes to
  # that member alone, e.g. "Klebsiella spp. (except K. aerogenes)" becomes
  # "Klebsiella spp." -- the exclusion itself needs no further handling, by
  # the same rank_index-precedence reasoning as the sheet-level case.
  # Named species groups such as "Streptococcus groups A, B, C and G" or
  # "Viridans group streptococci" must NOT be split; they are recognised by
  # not matching the "Genus species" shape on at least one piece, and are
  # returned unchanged so as.mo() resolves them directly to their
  # species-group code.
  
  # Resolve each distinct sheet's own organism once (not per row -- x and
  # sheet_name are typically thousands of rows spanning ~40 distinct
  # sheets), to use as an abbreviation seed (see expand_abbreviated_genus())
  # only when the sheet's organism is itself a single genus. An order-level
  # sheet (e.g. Enterobacterales) or a sheet whose name doesn't resolve at
  # all correctly gets no fallback, since the sheet-genus seed would either
  # be meaningless (an order is not a genus) or misleading (matching a
  # first letter shared with an unrelated genus, e.g. "Enterobacterales"
  # sharing "E" with "Escherichia" but not being that genus).
  distinct_sheets <- unique(sheet_name[!is.na(sheet_name)])
  genus_lookup <- setNames(rep(NA_character_, length(distinct_sheets)), distinct_sheets)
  for (sn in distinct_sheets) {
    sheet_mo <- suppressWarnings(suppressMessages(as.mo(sn, info = FALSE)))
    if (!is.na(sheet_mo) && isTRUE(mo_rank(sheet_mo) == "genus")) {
      genus_lookup[sn] <- mo_genus(sheet_mo)
    }
  }
  
  map2(x, sheet_name, function(s, sn) {
    sheet_genus_fallback <- if (!is.na(sn) && sn %in% names(genus_lookup)) genus_lookup[[sn]] else NA_character_
    split_mo_list_one(s, sheet_genus_fallback)
  })
}

split_mo_list_one <- function(s, sheet_genus_fallback) {
  s_clean <- trimws(gsub("[\r\n]+", " ", s))
  s_clean <- gsub("\u00A0", " ", s_clean, fixed = TRUE) # normalise non-breaking spaces (seen in some source cells)
  # Normalise a missing space after a genus-abbreviation period, e.g. a
  # v12.0 source typo "S.lugdunensis" -> "S. lugdunensis". Safe as a
  # blanket fix: it only matches an uppercase letter immediately followed
  # by "." and a lowercase letter with no space, which never occurs
  # correctly formed any other way in these cells.
  s_clean <- gsub("([A-Z])\\.(?=[a-z])", "\\1. ", s_clean, perl = TRUE)
  # Split on commas and top-level " and " (not inside parentheses)
  pieces <- unlist(strsplit(s_clean, ",(?![^(]*\\))|\\s+and\\s+(?![^(]*\\))", perl = TRUE))
  pieces <- trimws(pieces)
  pieces <- pieces[pieces != ""]
  pieces <- strip_mo_exceptions(pieces)
  pieces <- trimws(gsub("\\s*\\(\\s*\\)\\s*$", "", pieces)) # drop now-empty "()" leftovers
  
  if (length(pieces) <= 1) return(s_clean)
  
  # Only treat as a genuine multi-taxon list if every piece looks like a
  # taxon on its own: either an abbreviated genus initial ("E.", "P."
  # followed by a species epithet), a capitalised word of 2+ letters
  # (a full genus/organism name, "spp.", or the start of a longer
  # descriptive name such as "Streptococcus groups A"), a bare lowercase
  # species epithet with no genus marker at all, or the specific known
  # descriptive phrase "coagulase-negative staphylococci" (e.g.
  # "S. aureus and coagulase-negative staphylococci except S. epidermidis
  # and S. lugdunensis" on the Staphylococcus sheet, after the trailing
  # "except" clause is stripped by strip_mo_exceptions_wholecell() further
  # up, splits into "S. aureus" and this phrase). Other free-text
  # descriptive phrases are deliberately NOT accepted here even though
  # as.mo() can often resolve them: at its default matching tolerance
  # as.mo() will also resolve genuinely unrelated text to *something*
  # rather than fail, so a broad allowance here would risk silently
  # mapping real garbage to a wrong organism instead of correctly falling
  # back to the unsplit string. The bare lowercase epithet case is a
  # genuine, if unusual, EUCAST phrasing seen in older breakpoint tables
  # (v9.0-v12.0): "Aerococcus sanguinicola and urinae" and "Campylobacter
  # jejuni and coli" both drop the second species' genus abbreviation
  # entirely (newer versions write "and A. urinae" / "and C. coli"). Left
  # unhandled, the whole cell would either fail to resolve (Aerococcus
  # case) or silently resolve to the wrong thing -- as.mo("Campylobacter
  # jejuni and coli") matches a subspecies, not the two intended species
  # -- so it is accepted here and expanded against the preceding piece's
  # genus in expand_abbreviated_genus(), the same way an "X." abbreviation
  # is. A bare single uppercase letter with no following text or period --
  # e.g. the "B", "C", "G" that a naive comma/"and" split produces from
  # "Streptococcus groups A, B, C and G" -- still fails every one of these
  # shapes, correctly signalling that the comma/"and" split misfired on a
  # named group description rather than a genuine taxon list, so the
  # whole cell falls back to the original, unsplit string.
  looks_like_taxon <- grepl("^[A-Z]\\.\\s", pieces) | grepl("^[A-Z][a-z]+(\\.|\\s|$)", pieces) |
    grepl("^[a-z]+$", pieces) | grepl("^coagulase-negative staphylococci$", pieces, ignore.case = TRUE)
  if (!all(looks_like_taxon)) return(s_clean)
  
  # Reject a piece containing more than one capitalised, 3+ letter word
  # (a rough genus/proper-name detector): this signals two taxa were
  # fused into one piece because a separating comma was missing in the
  # source cell, e.g. a v12.0 typo "Citrobacter spp. Klebsiella spp."
  # (every other version correctly has a comma there). as.mo() does not
  # reliably fail on such a fused string -- it can silently match an
  # unrelated species with a similar combined shape (observed:
  # as.mo("Citrobacter spp. Klebsiella spp.") resolves to Citrobacter
  # enshiensis) -- so this must be caught before as.mo() ever sees it,
  # by falling back to the unsplit original string.
  genus_like_count <- vapply(pieces, function(p) length(gregexpr("\\b[A-Z][a-z]{2,}\\b", p, perl = TRUE)[[1]]),
                             integer(1))
  if (any(genus_like_count > 1)) return(s_clean)
  
  expand_abbreviated_genus(pieces, fallback_genus = sheet_genus_fallback)
}

breakpoints_eucast <- breakpoints_eucast_raw |>
  filter(!breakpoint_S %in% c("NOTE", "IE", "IP", "NA", "-"),
         !breakpoint_R %in% c("NOTE", "IE", "IP", "NA", "-")) |>
  mutate(
    mo_text = strip_mo_exceptions_wholecell(mo),
    mo_text = resolve_other_x(mo_text, sheet),
    mo_text = split_mo_list(mo_text, sheet_name = sheet)
  ) |>
  tidyr::unnest(mo_text) |>
  mutate(
    # ref_tbl is pure metadata for a human reviewer to trace a row back to
    # its exact source: the sheet name plus which file/version it came
    # from, e.g. "6. Aspergillus (file: Antifungal Clinical Breakpoint
    # Table v. 12.1)". file_desc is the version string's own description
    # of the table (see detect_version()), so this stays correct
    # automatically as new EUCAST versions are added, without needing any
    # per-file naming convention here.
    ref_tbl = if_else(!is.na(file_desc), paste0(trimws(sheet), " (file: ", file_desc, ")"), sheet),
    mo = as.mo(mo_text)
  ) |>
  select(-file_desc)

# AMR's internal diagnostic caches (AMR:::AMR_env$mo_uncertainties,
# AMR:::AMR_env$ab_previously_coerced) are reset by the *next* top-level
# as.mo()/as.ab() call, not preserved across a pipeline -- by the time the
# whole mutate() chain below finishes (four further mo_rank() calls, then
# as.ab()), only whatever the very last such call happened to leave behind
# would still be there, which is not a reliable record of every uncertain
# match made while building this table. Capture both immediately after the
# call that produced them, into script-owned variables, so they survive
# for inspection regardless of what AMR functions run afterward or how the
# script is subsequently sourced.
mo_uncertainties_eucast <- AMR:::AMR_env$mo_uncertainties
if (nrow(mo_uncertainties_eucast) > 0) {
  message("NOTE: ", nrow(mo_uncertainties_eucast), " organism string(s) were resolved with uncertainty ",
          "(matched below the default confidence, or against multiple candidates). ",
          "Inspect `mo_uncertainties_eucast` and cross-check against the source sheet before trusting these rows.")
}

breakpoints_eucast <- breakpoints_eucast |>
  mutate(
    rank_index = case_when(
      is.na(mo_rank(mo, keep_synonyms = TRUE)) ~ 6, # for UNKNOWN, B_GRAMN, B_ANAER, etc.
      mo_rank(mo, keep_synonyms = TRUE) %like% "(infra|sub)" ~ 1,
      mo_rank(mo, keep_synonyms = TRUE) == "species" ~ 2,
      mo_rank(mo, keep_synonyms = TRUE) == "species group" ~ 2.5,
      mo_rank(mo, keep_synonyms = TRUE) == "genus" ~ 3,
      mo_rank(mo, keep_synonyms = TRUE) == "family" ~ 4,
      mo_rank(mo, keep_synonyms = TRUE) == "order" ~ 5,
      TRUE ~ 6
    ),
    # The route/indication qualifier (if any) is captured into `site` before
    # it is stripped from `ab`, e.g. "Cefuroxime oral (uncomplicated UTI
    # only)" -> ab = "Cefuroxime oral", site = "uncomplicated UTI only".
    # This keeps the iv/oral or indication-specific variants of the same
    # drug+organism+method distinguishable in the distinct() step below
    # (site is part of its key); without it, e.g. "Cefuroxime iv" and
    # "Cefuroxime oral (uncomplicated UTI only)" against the same organism
    # would otherwise collide on identical (ab, mo, method, breakpoint_S)
    # if their S breakpoints happen to match, silently losing one of them.
    # uti is derived from that same qualifier text: every "UTI" mention in
    # the raw agent-name text is confirmed to sit inside a trailing
    # parenthetical qualifier, never in the base drug name itself, so no
    # separate check against the (by this point already-coded) ab is needed.
    qualifier = extract_trailing_parenthetical(ab),
    site = coalesce(site, qualifier),
    ab_input = gsub("\\s*\\(.*$", "", ab), # also fixes a source cell missing a space before "(", e.g. "Meropenem(indications..."
    ab = as.ab(ab_input)
  )

# ab_previously_coerced also logs every intermediate word-by-word fallback
# as.ab() tries internally while resolving a route-suffixed name (e.g. for
# "Ampicillin iv", it separately tries "IV", "AMPICILLIN", and "ORAL"-like
# fragments as candidate matches before landing on the correct whole-string
# result) -- these are internal matching noise, not genuine coercions of
# anything actually in this table, and would otherwise dominate the log.
# Restrict to rows whose logged input is one of the actual, full ab_input
# strings this pipeline passed in, which is what a reviewer needs to see.
ab_previously_coerced_eucast <- AMR:::AMR_env$ab_previously_coerced |>
  filter(tolower(trimws(x_bak)) %in% tolower(trimws(unique(breakpoints_eucast$ab_input))))
breakpoints_eucast <- breakpoints_eucast |> select(-ab_input)
if (nrow(ab_previously_coerced_eucast) > 0) {
  message("NOTE: ", nrow(ab_previously_coerced_eucast), " antibiotic string(s) were resolved via approximate/previous-coercion ",
          "matching rather than an exact code or name match. ",
          "Inspect `ab_previously_coerced_eucast` and cross-check against the source sheet before trusting these rows.")
}

breakpoints_eucast <- breakpoints_eucast |>
  mutate(
    breakpoint_S = as.numeric(breakpoint_S),
    breakpoint_R = as.numeric(breakpoint_R),
    uti = qualifier %like_case% "UTI",
    is_SDD = FALSE # EUCAST has no "susceptible, dose dependent" category (CLSI-only concept)
  ) |>
  select(-mo_text, -qualifier) |>
  # Greek symbols and EM dash symbols are not allowed by CRAN, so replace them with ASCII:
  mutate(disk_dose = disk_dose %>%
           gsub("\u03bc", "mc", ., fixed = TRUE) %>% # this is 'mu', \u03bc
           gsub("\u00b5", "mc", ., fixed = TRUE) %>% # this is 'micro', \u00b5 (yes, they look the same)
           gsub("\u2013", "-", ., fixed = TRUE) %>%
           gsub("(?<=\\d)(?=[a-zA-Z])", " ", ., perl = TRUE)) # keep a space after a number, e.g. "1mcg" to "1 mcg"

# Surface any rows whose organism or antibiotic text failed to resolve,
# rather than letting the final filter() drop them silently -- a growing
# list here across EUCAST versions usually means a new free-text pattern
# needs handling above, not a data error.
#
# as.mo() returns the literal code "UNKNOWN" rather than NA when it finds
# no match at all (unlike as.ab(), which returns NA), so a plain
# is.na(mo) check misses genuine resolution failures entirely. UNKNOWN is
# not automatically wrong to keep -- it is a real, meaningful category
# already used throughout this table (e.g. the rank_index case_when
# above), and the currently published clinical_breakpoints legitimately
# contains hundreds of such rows for cases EUCAST itself left organism-
# unspecific. It is reported here for review rather than dropped
# outright, since only a source-text defect that split_mo_list_one()
# doesn't yet handle (checked by cross-referencing ref_tbl/ab against the
# source sheet) would make a specific UNKNOWN row here spurious rather
# than a genuine EUCAST ambiguity.
unresolved_mo <- breakpoints_eucast |> filter(is.na(mo) | mo == "UNKNOWN") |> distinct(sheet, ref_tbl, ab)
unresolved_ab <- breakpoints_eucast |> filter(is.na(ab)) |> distinct(sheet, ab)
if (nrow(unresolved_mo) > 0) {
  message("NOTE: ", nrow(unresolved_mo), " row(s) have an organism that either could not be resolved (NA, will be dropped) ",
          "or resolved to UNKNOWN (kept, but review whether this is a genuine EUCAST ambiguity or a parsing gap):")
  print(unresolved_mo)
}
if (nrow(unresolved_ab) > 0) {
  message("NOTE: ", nrow(unresolved_ab), " distinct antibiotic string(s) could not be resolved by as.ab() and will be dropped:")
  print(unresolved_ab)
}

breakpoints_eucast <- breakpoints_eucast |>
  filter(!(is.na(breakpoint_S) & is.na(breakpoint_R)), !is.na(mo), !is.na(ab)) |>
  distinct(guideline, type, host, ab, mo, method, site, breakpoint_S, .keep_all = TRUE) |>
  select(guideline, type, host, method, site, mo, rank_index, ab, ref_tbl,
         disk_dose, breakpoint_S, breakpoint_R, uti, is_SDD, note) |>
  arrange(desc(guideline), mo, ab, type, host, method)

breakpoints_eucast %>% count(guideline)
glimpse(breakpoints_eucast)
