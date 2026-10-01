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

# Log of deliberate deviations from the literal cell contents, reported after
# parsing so that every one of them can be reviewed
parse_log <- new.env()
parse_log$disk_dose_inherited <- NULL
parse_log$mic_scope_extensions <- NULL
parse_log$no_breakpoints_determined <- NULL

# ==============================================================================
# 1. Rich text parser
# ==============================================================================
# Struck-through runs (tidyxl: character_formatted$strike == TRUE) are always
# disregarded before splitting base text from superscript note references,
# since EUCAST uses strike-through to mark removed/superseded text within a
# cell (most commonly within Notes blocks) that must not end up in the parsed
# output.
split_rich_text <- function(fmt_list, cell_strike = rep(FALSE, length(fmt_list))) {
  map2_dfr(fmt_list, cell_strike, function(fmt, cell_struck) {
    if (is.null(fmt) || nrow(fmt) == 0)
      return(tibble(is_rich = FALSE, base_text = NA_character_, note_super = NA_character_, full_text = NA_character_))
    # A run without its own strike setting inherits the cell's font
    is_struck <- coalesce(fmt$strike, isTRUE(cell_struck))
    fmt <- fmt[!is_struck, , drop = FALSE]
    if (nrow(fmt) == 0)
      return(tibble(is_rich = TRUE, base_text = NA_character_, note_super = NA_character_, full_text = NA_character_))
    is_super <- !is.na(fmt$vertAlign) & fmt$vertAlign == "superscript"
    tibble(
      is_rich    = TRUE,
      base_text  = paste0(fmt$character[!is_super], collapse = ""),
      note_super = paste0(fmt$character[is_super],  collapse = ""),
      # Text including superscripts, for note bodies, where a superscript is
      # part of the text rather than a reference, e.g. "Ca2+" or "Mg2+"
      full_text  = paste0(fmt$character, collapse = "")
    )
  })
}

# Parses all cells of a sheet into row, col, base_text (without superscripts
# and struck-through text), note_super (superscript note references) and
# full_text (with superscripts, for note bodies). For rich-text cells, the
# result of split_rich_text() is authoritative, also when it is empty because
# the whole cell is struck through: falling back to the raw cell text there
# would bring deleted text back. Plain (non-rich-text) cells with a struck-
# through cell font are disregarded entirely, for the same reason.
parse_cells <- function(cells, formats) {
  cell_strike <- coalesce(formats$local$font$strike[cells$local_format_id], FALSE)
  target <- cells[!cells$is_blank, , drop = FALSE]
  target_strike <- cell_strike[!cells$is_blank]
  rich <- split_rich_text(target$character_formatted, target_strike)
  target |>
    select(-any_of(c("full_text", "base_text", "note_super"))) |>
    bind_cols(rich) |>
    mutate(
      plain_struck = !is_rich & target_strike,
      base_text = case_when(
        plain_struck                        ~ NA_character_,
        is_rich                             ~ na_if(trimws(base_text), ""),
        !is.na(character)                   ~ trimws(character),
        !is.na(numeric)                     ~ as.character(numeric),
        TRUE                                ~ NA_character_
      ),
      full_text = case_when(
        plain_struck                        ~ NA_character_,
        is_rich                             ~ na_if(trimws(full_text), ""),
        TRUE                                ~ base_text
      ),
      note_super = if_else(is.na(note_super) | note_super == "",
                           NA_character_, note_super)
    ) |>
    filter(!is.na(base_text) | !is.na(note_super) | !is.na(full_text)) |>
    select(row, col, base_text, note_super, full_text)
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
  base_strike <- coalesce(formats$local$font$strike[col_a$local_format_id], FALSE)
  
  out <- vector("list", nrow(col_a))
  for (i in seq_len(nrow(col_a))) {
    fmt <- col_a$character_formatted[[i]]
    row_i <- col_a$row[i]
    
    if (is.null(fmt) || nrow(fmt) == 0) {
      # Plain (non-rich-text) cell: nothing to split, whole text is the agent
      # name, unless the whole cell is struck through
      if (!base_strike[i]) {
        out[[i]] <- tibble(row = row_i, ab_text = col_a$character[i], mo_override = NA_character_)
      }
      next
    }
    
    # Drop struck-through and superscript runs first, exactly as split_rich_text() does
    is_struck <- coalesce(fmt$strike, base_strike[i])
    fmt <- fmt[!is_struck, , drop = FALSE]
    is_super <- !is.na(fmt$vertAlign) & fmt$vertAlign == "superscript"
    fmt <- fmt[!is_super, , drop = FALSE]
    
    if (nrow(fmt) <= 1) {
      txt <- if (nrow(fmt) == 1) fmt$character[1] else col_a$character[i]
      out[[i]] <- tibble(row = row_i, ab_text = txt, mo_override = NA_character_)
      next
    }
    
    bold_resolved <- ifelse(is.na(fmt$bold), base_bold[i], fmt$bold)

    # An indication qualifier may also follow the organism restriction as a
    # trailing bold parenthetical, e.g. Staphylococcus sheet v9.0/v10.0:
    # "Ceftaroline" (bold) ", S. aureus" (not bold) "(pneumonia)" (bold),
    # which v11.0 onwards writes as "Ceftaroline (pneumonia), S. aureus".
    # Such a trailing bold run belongs to the agent name, not to the
    # organism, so it is set aside here and appended to the agent name.
    trailing_qualifier <- NA_character_
    n <- length(bold_resolved)
    if (n > 2 && bold_resolved[n] && !bold_resolved[n - 1] &&
        grepl("^\\s*\\([^()]+\\)\\s*$", fmt$character[n])) {
      trailing_qualifier <- trimws(fmt$character[n])
      fmt <- fmt[-n, , drop = FALSE]
      bold_resolved <- bold_resolved[-n]
    }
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
                           ab_text = paste(c(paste0(fmt$character, collapse = ""),
                                             trailing_qualifier[!is.na(trailing_qualifier)]),
                                           collapse = " "),
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
    if (!is.na(trailing_qualifier)) {
      ab_text <- paste(gsub(",\\s*$", "", trimws(ab_text)), trailing_qualifier)
    }
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
      mo_override = na_if(mo_override, ""),
      # a non-bold part that is entirely a parenthetical is an indication, not
      # an organism, e.g. v12.0-v16.1 streptococci: "Mecillinam oral
      # (pivmecillinam)" (bold) "(uncomplicated UTI only)" (not bold)
      is_indication = !is.na(mo_override) & grepl("^\\([^()]+\\)$", mo_override),
      ab_text = if_else(is_indication, paste(ab_text, mo_override), ab_text),
      mo_override = if_else(is_indication, NA_character_, mo_override)
    ) |>
    select(-is_indication)
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
    keys <- gsub("[.\\s]+$", "", trimws(substring(txt, m, m + lens - 1)), perl = TRUE)
    data.frame(pos = as.integer(m), key = keys, mlen = as.integer(lens),
               stringsAsFactors = FALSE)
  }
  
  # A note key only starts at the beginning of a line or directly after the
  # end of a sentence, since EUCAST does not always put each note on its own
  # line (e.g. "...meropenem only. 3.The addition of..." or "...urinary
  # tract infections only.2/A. Susceptibility to..."). A sentence end is a
  # period followed by whitespace, or a period directly preceded by a
  # non-digit: this keeps decimals such as "MIC <=0.5. Isolates" from being
  # read as a key "5". Cross-references within a note, such as "See Note B."
  # or "Notes 5/D and 7/F", are preceded by a word and therefore never match.
  # Keys are not always followed by a space either, e.g. "2.Gentamicin can
  # be used..." or "C.Isolates categorised...", so a directly following
  # uppercase letter (or opening parenthesis) is accepted too.
  key_start <- "(?:^|(?<=\\n)|(?<=\\.\\s)|(?<=[^0-9\\s]\\.))[ \\t]*"
  # Combined key: "1/A." or "5/A."
  r1 <- find_keys(paste0(key_start, "\\d+/[A-Z]\\.(?:\\s+|(?=[A-Z(]))"), txt)
  # Numbered with dot: "1." "2."
  r2 <- find_keys(paste0(key_start, "\\d+\\.(?:\\s+|(?=[A-Z(]))"), txt)
  # Lettered: "A." followed by an uppercase letter, which is not a species
  # name such as "C. difficile"
  r3 <- find_keys(paste0(key_start, "[A-Z]\\.(?:\\s+(?=[A-Z(])|(?=[A-Z][a-z]))"), txt)
  # Numbered without dot at the start of a line (e.g. C. difficile sheet):
  # digit, space, uppercase letter
  r4 <- find_keys("(?:^|(?<=\\n))[ \\t]*\\d+[ \\t]+(?=[A-Z])", txt)
  
  all_keys <- rbind(r1, r2, r3, r4)
  if (nrow(all_keys) == 0) return(list())
  # At the same position, keep the longest (most specific) match, e.g. "1/A."
  # over "1.", and drop any match that starts within an earlier kept match
  all_keys <- all_keys[order(all_keys$pos, -all_keys$mlen), ]
  all_keys <- all_keys[!duplicated(all_keys$pos), ]
  keep <- rep(TRUE, nrow(all_keys))
  for (i in seq_len(nrow(all_keys))) {
    if (!keep[i]) next
    keep[seq_len(nrow(all_keys)) > i & all_keys$pos < all_keys$pos[i] + all_keys$mlen[i]] <- FALSE
  }
  all_keys <- all_keys[keep, ]
  
  notes <- list()
  for (i in seq_len(nrow(all_keys))) {
    body_start <- all_keys$pos[i] + all_keys$mlen[i]
    body_end <- if (i < nrow(all_keys)) all_keys$pos[i + 1] - 1 else nchar(txt)
    body <- trimws(substr(txt, body_start, body_end))
    key <- all_keys$key[i]
    # A combined key "1/A" also serves the MIC reference "1" and the disk
    # reference "A" individually
    for (k in unique(c(key, if (grepl("/", key)) strsplit(key, "/")[[1]]))) {
      if (!nzchar(body)) {
        stop("parse_notes_block(): note key '", key, "' has an empty body, which means a ",
             "cross-reference or other text was taken for a key. Text: '", substr(txt, 1, 300), "'")
      }
      if (!is.null(notes[[k]]) && !identical(notes[[k]], body)) {
        stop("parse_notes_block(): note key '", k, "' occurs more than once with different ",
             "bodies. Text: '", substr(txt, 1, 300), "'")
      }
      notes[[k]] <- body
    }
  }
  notes
}

# ==============================================================================
# 3. Note resolver
# ==============================================================================
resolve_notes <- function(notes_list, mic_super, disk_super, abbreviations = list()) {
  # References are normally separated by commas, but occasionally by a
  # period or a space instead, e.g. "0.125" with superscript "3.4" (v9.0,
  # Staphylococcus) or "1" with superscript "4.5" (v16.1, H. influenzae)
  mic_refs <- character(0)
  disk_refs <- character(0)
  if (!is.na(mic_super) && mic_super != "")
    mic_refs <- trimws(unlist(strsplit(mic_super, "[,.[:space:]]+")))
  if (!is.na(disk_super) && disk_super != "")
    disk_refs <- trimws(unlist(strsplit(disk_super, "[,.[:space:]]+")))
  all_refs <- unique(c(mic_refs, disk_refs))
  all_refs <- all_refs[nzchar(all_refs)]
  if (length(all_refs) == 0) return(NA_character_)

  # A reference that is not a key of the table's own notes block can be an
  # abbreviation from the workbook's general "Notes" sheet, e.g. superscript
  # "HE" in v9.0 ("HE = High exposure for agent (see table of dosages...)").
  # Anything else would silently lose a note, so it is an error.
  nl <- if (is.null(notes_list)) list() else notes_list
  from_abbr <- setdiff(intersect(all_refs, names(abbreviations)), names(nl))
  nl <- c(nl, abbreviations[from_abbr])
  unresolved <- setdiff(all_refs, names(nl))
  if (length(unresolved) > 0) {
    stop("resolve_notes(): note reference(s) '", paste(unresolved, collapse = "', '"),
         "' not found in the notes block (keys found: ", paste(names(notes_list), collapse = " "),
         ") nor in the workbook abbreviations (", paste(names(abbreviations), collapse = " "), ")")
  }
  
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
  base <- hit$base_text[1]
  super <- hit$note_super[1]
  # EUCAST occasionally omits the superscript formatting of a note reference
  # on a breakpoint value, e.g. "33A" (C. acnes, ampicillin-sulbactam disk,
  # v16.1). Move such a trailing letter from the base text to the references.
  m <- regmatches(base, regexec("^(\\(?[0-9.]+\\)?)([A-Z])$", base))[[1]]
  if (length(m) == 3) {
    base <- m[2]
    super <- if (is.na(super) || super == "") m[3] else paste(super, m[3], sep = ",")
  }
  list(base = base, super = super)
}

parse_bp <- function(txt) {
  if (is.na(txt) || txt == "") return(NA_character_)
  t <- trimws(txt)
  # placeholders such as "-" (susceptibility testing not recommended), "IE"
  # or "Note" are kept as they are: the cleanup turns them into blocking rows
  # Parenthesised values like "(2)" are EUCAST's "breakpoints in brackets"
  # (see is_bracketed_bp()). They are preserved as-is in the raw file so
  # they remain distinguishable from clinical breakpoints. The cleanup
  # section later excludes them from the clinical_breakpoints table.
  t
}

# Cells without a breakpoint hold a placeholder: "Note" (see the notes),
# "IE" (insufficient evidence), "IP" (in preparation), "ND" (not determined,
# Topical agents sheet), "NA" (not applicable) or "-"
bp_placeholders <- c("NOTE", "IE", "IP", "ND", "NA", "-")
is_bp_value <- function(txt) {
  !is.na(txt) & nzchar(trimws(txt)) & !toupper(trimws(txt)) %in% bp_placeholders
}

is_bracketed_bp <- function(txt) {
  # A breakpoint cell value wrapped in parentheses, e.g. "(2)" or "(0.5)",
  # is what EUCAST calls a "breakpoint in brackets" (used since v10.0): a
  # value based on ECOFFs that distinguishes isolates without and with
  # phenotypically detectable resistance mechanisms. EUCAST states that
  # clinical evidence as monotherapy is usually lacking and that reporting
  # S or I should be avoided. These are not screening breakpoints: EUCAST
  # marks screening tests in the agent name instead, e.g. "Cefoxitin
  # (screen only)". The definition is read per workbook by
  # parse_workbook_notes() and added to the note of every bracketed row.
  if (is.na(txt) || txt == "") return(FALSE)
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
parse_sheet <- function(xlsx_path, sheet_name, abbreviations = list()) {
  cells <- xlsx_cells(xlsx_path, sheets = sheet_name)
  
  # --- Detect version ---
  ver <- detect_version(cells)
  
  formats <- xlsx_formats(xlsx_path)
  
  # --- Parse rich text for cols A:I ---
  parsed <- parse_cells(cells |> filter(col >= 1, col <= 9), formats)
  
  # --- Statements extending MIC breakpoints to another species ---
  # The H. influenzae sheet states (v9.0-v16.1): "In the absence of specific
  # breakpoints, the H. influenzae MIC breakpoints can be applied to H.
  # parainfluenzae." This applies to MIC breakpoints only, not to zone
  # diameters. Such statements are logged here (with the sheet's own text)
  # and applied in the cleanup section, see `parse_log$mic_scope_extensions`.
  # Any other "breakpoints can be applied to <species>" statement that this
  # pattern does not recognise stops the script, so that it is not missed.
  sheet_text <- unique(parsed$full_text[!is.na(parsed$full_text) & parsed$col == 1])
  scope_pattern <- "In the absence of specific breakpoints, the ([A-Z])\\. ([a-z]+) MIC ?breakpoints can be applied to ([A-Z])\\. ([a-z]+)\\."
  for (txt in sheet_text) {
    txt_squished <- gsub("\\s+", " ", txt)
    m <- regmatches(txt_squished, regexec(scope_pattern, txt_squished))[[1]]
    if (length(m) == 5) {
      sheet_title <- trimws(get_cell(parsed, 1, 1)$base)
      genus <- sub("\\s.*$", "", sheet_title)
      if (substr(genus, 1, 1) != m[2] || m[2] != m[4]) {
        stop("parse_sheet(): MIC scope statement on sheet '", sheet_name, "' does not match the sheet's genus '", genus, "': ", m[1])
      }
      parse_log$mic_scope_extensions <- rbind(
        parse_log$mic_scope_extensions,
        data.frame(guideline = ver$guideline, sheet = sheet_name,
                   from_text = paste(genus, m[3]), to_text = paste(genus, m[5]), statement = m[1])
      )
    } else if (grepl("breakpoints can (also )?be applied to [A-Z]\\.", txt_squished)) {
      stop("parse_sheet(): unrecognised statement about applying breakpoints to another species on sheet '",
           sheet_name, "': ", substr(txt_squished, 1, 300))
    }
  }
  
  # --- Detect header rows (col 2 contains "MIC breakpoint") ---
  header_rows <- parsed |>
    filter(col == 2, grepl("MIC breakpoint", base_text, fixed = TRUE)) |>
    pull(row) |>
    sort()
  
  if (length(header_rows) == 0) {
    # Only acceptable when EUCAST states on the sheet itself that it has not
    # determined breakpoints, e.g. "EUCAST has not determined breakpoints for
    # Burkholderia cepacia complex organisms..." (likewise L. pneumophila).
    # Any other sheet without a breakpoint table means its layout was not
    # recognised, which must not silently drop the sheet.
    # These organisms are logged, so that the non-species related PK-PD
    # breakpoints (v9.0-v13.1) can be blocked for them in the cleanup.
    no_bp_statement <- parsed$full_text[!is.na(parsed$full_text) & grepl("has not determined breakpoints", parsed$full_text, fixed = TRUE)]
    if (length(no_bp_statement) > 0) {
      parse_log$no_breakpoints_determined <- rbind(
        parse_log$no_breakpoints_determined,
        data.frame(guideline = ver$guideline, sheet = sheet_name,
                   organism = trimws(get_cell(parsed, 1, 1)$base),
                   statement = gsub("\\s+", " ", trimws(no_bp_statement[1])))
      )
      message("  EUCAST has not determined breakpoints in sheet '", sheet_name, "', skipping")
      return(NULL)
    }
    # From v14.0 on, the PK-PD sheet holds no table anymore but an explanation
    # why EUCAST withdrew the non-species related breakpoints ("A common
    # misunderstanding is that PK/PD breakpoints are overarching ... This is
    # not the intention."), and from v15.0 on it is called "PK/PD cut-off
    # values". Such a sheet has nothing to extract.
    if (grepl("^PK.?PD", trimws(coalesce(get_cell(parsed, 1, 1)$base, "")))) {
      message("  PK-PD sheet '", sheet_name, "' holds no breakpoint table, skipping")
      return(NULL)
    }
    stop("parse_sheet(): no 'MIC breakpoint' header row found in sheet '", sheet_name, "'")
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
      note_text <- if (nrow(note_cell) > 0) note_cell$full_text[1] else NA_character_
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
        paste(note_cells$full_text, collapse = "\n")
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
    prev_agent_base <- NA_character_
    prev_disk_dose <- NA_character_
    
    for (r in tbl$first_data:tbl$last_data) {
      agent <- get_cell(parsed, r, 1)
      # Skip category headers (rows that have col A text but no data).
      # Columns 5:7 only hold disk data when the sheet actually has disk
      # columns; on MIC-only sheets col 5 is the notes column and must not
      # be treated as a data indicator.
      data_cols <- if (has_disk_columns) c(2, 3, 5, 6, 7) else c(2, 3)
      has_any_data <- any(!is.na(vapply(data_cols, function(cc) get_cell(parsed, r, cc)$base,
                                        character(1))))
      if (is.na(agent$base)) {
        # A row with values but without an agent name would otherwise be
        # skipped silently, e.g. if a future table merged agent cells. Rows
        # holding only "-" (such as an empty spacer row) carry nothing.
        bp_cols <- if (has_disk_columns) c(2, 3, 6, 7) else c(2, 3)
        bp_cells <- vapply(bp_cols, function(cc) coalesce(get_cell(parsed, r, cc)$base, ""), character(1))
        if (any(is_bp_value(bp_cells) | toupper(trimws(bp_cells)) %in% c("IE", "NOTE", "IP"))) {
          stop("parse_sheet(): row ", r, " of sheet '", sheet_name, "' has breakpoint data but no agent name")
        }
        next
      }
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
      
      # EUCAST occasionally leaves the disk content empty in a row that
      # continues the previous row's agent, e.g. v16.1 "Amoxicillin-
      # clavulanic acid, F. nucleatum" directly below "Amoxicillin-
      # clavulanic acid, F. necrophorum" (2-1 mcg), or v15.0 VGS
      # "Benzylpenicillin (endocarditis)" directly below "Benzylpenicillin
      # (indications other than endocarditis)" (1 unit; filled in in
      # v16.1). Only then is the disk content of that previous row used.
      # Every such case is logged; any other disk breakpoint without a disk
      # content stops the script in the integrity checks of the cleanup.
      agent_base <- tolower(trimws(sub("\\s+(iv|oral)$", "", sub("\\s*[(,].*$", "", agent_text), ignore.case = TRUE)))
      disk_dose <- format_disk_dose(disk_dose_cell$base)
      has_disk_bp <- is_bp_value(disk_s$base) || is_bp_value(disk_r$base)
      if (is.na(disk_dose) && has_disk_bp && identical(agent_base, prev_agent_base) && !is.na(prev_disk_dose)) {
        disk_dose <- prev_disk_dose
        parse_log$disk_dose_inherited <- rbind(
          parse_log$disk_dose_inherited,
          data.frame(file = basename(xlsx_path), sheet = sheet_name, row = r, agent = agent_text, disk_dose = disk_dose)
        )
      }
      prev_agent_base <- agent_base
      prev_disk_dose <- disk_dose
      
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
      mic_note  <- resolve_notes(nl, combine_super(mic_s$super, agent$super), NA, abbreviations = abbreviations)
      disk_note <- resolve_notes(nl, NA, combine_super(disk_s$super, agent$super), abbreviations = abbreviations)
      
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
          table_organism = org,
          rank_index   = NA_integer_,
          ab           = agent_text,
          disk_dose    = NA_character_,
          breakpoint_S = s_val,
          breakpoint_R = r_val,
          is_bracketed = is_bracketed_bp(mic_s$base) || is_bracketed_bp(mic_r$base),
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
          table_organism = org,
          rank_index   = NA_integer_,
          ab           = agent_text,
          disk_dose    = disk_dose,
          breakpoint_S = s_val,
          breakpoint_R = r_val,
          is_bracketed = is_bracketed_bp(disk_s$base) || is_bracketed_bp(disk_r$base),
          note         = disk_note
        )
      }
    }
  }
  
  if (idx == 0) return(NULL)
  out <- bind_rows(results[seq_len(idx)])
  
  # The non-species related PK-PD breakpoints (v9.0-v13.1) come with EUCAST's
  # condition for their use, e.g. "These breakpoints are used only when there
  # are no species-specific breakpoints or other recommendations (a dash or a
  # note) in the species-specific tables. ...". This statement is added to
  # the note of every PK-PD row.
  if (grepl("Non-species related", organism, ignore.case = TRUE)) {
    usage <- parsed$full_text[!is.na(parsed$full_text) & parsed$col == 1 &
                                grepl("^These breakpoints are used only when", trimws(parsed$full_text))]
    if (length(usage) != 1) {
      stop("parse_sheet(): the condition for using the PK-PD breakpoints was not found on sheet '", sheet_name, "'")
    }
    usage <- paste0("[PK-PD] ", gsub("\\s+", " ", trimws(usage)))
    out$note <- if_else(is.na(out$note), usage, paste(usage, out$note, sep = " | "))
  }
  out
}

# ==============================================================================
# 7. Parse "Topical agents" sheet (transposed layout)
# ==============================================================================
parse_topical_sheet <- function(xlsx_path, abbreviations = list()) {
  cells <- xlsx_cells(xlsx_path, sheets = "Topical agents")
  
  ver <- detect_version(cells)
  
  parsed <- parse_cells(cells, xlsx_formats(xlsx_path))

  # v9.0 has a different, informational table here ("ECOFFs and systemic
  # clinical breakpoints for antimicrobial agents that are used topically"):
  # EUCAST states it could not reach consensus on topical breakpoints and
  # presents a mix of systemic clinical breakpoints and ECOFFs "for
  # information" only, so it holds no topical breakpoints to extract.
  title <- get_cell(parsed, 1, 1)$base
  if (!is.na(title) && grepl("^ECOFFs and systemic clinical breakpoints", title)) {
    message("  Informational table without topical breakpoints in '", basename(xlsx_path), "', skipping")
    return(NULL)
  }
  # From v10.0 on, the layout is fixed: agent names in row 6 from column D
  # on, disk contents in row 10, and organism rows (MIC row followed by a
  # zone diameter row) from row 11 on. Verify this rather than assume it.
  if (!identical(get_cell(parsed, 6, 1)$base, "Organisms") ||
      !identical(get_cell(parsed, 10, 2)$base, "Disk content") ||
      !grepl("^Screening cut-off values", coalesce(get_cell(parsed, 6, 2)$base, ""))) {
    stop("parse_topical_sheet(): unexpected layout of sheet 'Topical agents' in '", basename(xlsx_path), "'")
  }

  # Antimicrobial names from row 6, column D up to the last populated column
  ab_cols <- seq(4, max(parsed$col[parsed$row == 6]))
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
    paste(notes_cells$full_text, collapse = "\n") else ""
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
      mic_note <- resolve_notes(notes_lookup, paste(all_supers, collapse = ","), NA, abbreviations = abbreviations)
      
      if (!is.na(mic_val)) {
        idx <- idx + 1L
        results[[idx]] <- tibble(
          guideline = ver$guideline, version = ver$version, file_desc = ver$file_desc,
          type = "human", host = "human", method = "MIC",
          site = "Topical", mo = od$organism, table_organism = od$organism, rank_index = NA_integer_,
          ab = ab_name,
          disk_dose = NA_character_,
          breakpoint_S = mic_val, breakpoint_R = mic_val,
          is_bracketed = is_bracketed_bp(mic_cell$base),
          note = mic_note
        )
      }
      
      # DISK
      disk_cell <- get_cell(parsed, od$disk_row, col_idx)
      disk_val <- parse_bp(disk_cell$base)
      all_supers_d <- c(disk_cell$super, ab_super)
      all_supers_d <- all_supers_d[!is.na(all_supers_d) & all_supers_d != ""]
      disk_note <- resolve_notes(notes_lookup, NA, paste(all_supers_d, collapse = ","), abbreviations = abbreviations)
      
      if (!is.na(disk_val)) {
        idx <- idx + 1L
        results[[idx]] <- tibble(
          guideline = ver$guideline, version = ver$version, file_desc = ver$file_desc,
          type = "human", host = "human", method = "DISK",
          site = "Topical", mo = od$organism, table_organism = od$organism, rank_index = NA_integer_,
          ab = ab_name,
          disk_dose = ab_dose,
          breakpoint_S = disk_val, breakpoint_R = disk_val,
          is_bracketed = is_bracketed_bp(disk_cell$base),
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
parse_yeast_sheet <- function(xlsx_path, sheet_name, abbreviations = list()) {
  cells <- xlsx_cells(xlsx_path, sheets = sheet_name)
  
  ver <- detect_version(cells)
  
  parsed <- parse_cells(cells, xlsx_formats(xlsx_path))
  
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
    paste(note_cells$full_text, collapse = "\n")
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
      note_text <- resolve_notes(notes_lookup, paste(all_supers, collapse = ","), NA, abbreviations = abbreviations)
      
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
        table_organism = sb$species,
        rank_index   = NA_integer_,
        ab           = agent_name,
        disk_dose    = NA_character_,
        breakpoint_S = s_val,
        breakpoint_R = r_val,
        is_bracketed = is_bracketed_bp(s_cell$base) || is_bracketed_bp(r_cell$base),
        note         = note_text
      )
    }
  }
  
  if (idx == 0) return(NULL)
  bind_rows(results[seq_len(idx)])
}

# ==============================================================================
# 9. Workbook-level notes: abbreviations and the definition of brackets
# ==============================================================================
# Besides the notes block of each table, every workbook has a general "Notes"
# sheet ("Notes" or e.g. "1. Notes") with two things the tables rely on:
#  - Abbreviations of the form "XX = text", e.g. "HE = High exposure for agent
#    (see table of dosages, last tab in breakpoint table)" in v9.0 and v10.0,
#    where "HE" is used as a superscript on agent names.
#  - EUCAST's definition of "breakpoints in brackets", which differs between
#    versions: within note 7 (v10.0), as "( ) = ..." (v11.0), or as a note of
#    its own starting with "Breakpoints in brackets" (v12.0 on).
# Both are read from the workbook itself, so they always match its version.
# Names of the hidden sheets of a workbook (tidyxl does not report sheet
# visibility, so it is read from the workbook definition itself)
xlsx_hidden_sheets <- function(xlsx_path) {
  tmp <- tempfile()
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE))
  utils::unzip(xlsx_path, files = "xl/workbook.xml", exdir = tmp)
  wb <- paste(readLines(file.path(tmp, "xl", "workbook.xml"), warn = FALSE, encoding = "UTF-8"), collapse = "")
  sheets <- regmatches(wb, gregexpr("<sheet [^>]*/>", wb))[[1]]
  hidden <- sheets[grepl('state="(hidden|veryHidden)"', sheets)]
  nm <- sub('.*name="([^"]*)".*', "\\1", hidden)
  nm <- gsub("&amp;", "&", gsub("&apos;", "'", gsub("&quot;", "\"", nm, fixed = TRUE), fixed = TRUE), fixed = TRUE)
  if (!all(nm %in% xlsx_sheet_names(xlsx_path))) {
    stop("xlsx_hidden_sheets(): could not match hidden sheet names in '", basename(xlsx_path), "'")
  }
  nm
}

parse_workbook_notes <- function(xlsx_path) {
  all_sheets <- xlsx_sheet_names(xlsx_path)
  notes_sheet <- all_sheets[grepl("^\\s*(\\d+\\.\\s*)?Notes\\s*$", all_sheets, ignore.case = TRUE)]
  if (length(notes_sheet) != 1) {
    stop("parse_workbook_notes(): expected exactly one 'Notes' sheet in '", basename(xlsx_path),
         "', found: ", paste(notes_sheet, collapse = ", "))
  }
  txt <- parse_cells(xlsx_cells(xlsx_path, sheets = notes_sheet), xlsx_formats(xlsx_path)) |>
    filter(!is.na(full_text)) |>
    pull(full_text) |>
    gsub(pattern = "\u00A0", replacement = " ", fixed = TRUE) |>
    trimws()
  
  abbr_txt <- txt[grepl("^[A-Z]{2,5}\\s*=\\s*\\S", txt)]
  abbreviations <- as.list(trimws(sub("^[A-Z]{2,5}\\s*=\\s*", "", abbr_txt)))
  names(abbreviations) <- sub("\\s*=.*$", "", abbr_txt)
  
  bracket_definition <- NA_character_
  def_brackets_sign <- txt[grepl("^\\(\\s*\\)\\s*=", txt)]
  def_own_note <- sub("^\\d+\\.\\s*", "", txt[grepl("^(\\d+\\.\\s*)?Breakpoints in brackets", txt)])
  def_in_note <- txt[grepl("Breakpoints in brackets", txt, fixed = TRUE)]
  if (length(def_brackets_sign) > 0) {
    bracket_definition <- sub("^\\(\\s*\\)\\s*=\\s*", "", def_brackets_sign[1])
  } else if (length(def_own_note) > 0) {
    bracket_definition <- def_own_note[1]
  } else if (length(def_in_note) > 0) {
    bracket_definition <- sub("^.*?(Breakpoints in brackets)", "\\1", def_in_note[1])
  }
  bracket_definition <- trimws(gsub("\\s+", " ", bracket_definition))
  
  list(abbreviations = abbreviations, bracket_definition = bracket_definition)
}

# ==============================================================================
# 10. Parse all data sheets
# ==============================================================================
parse_workbook <- function(xlsx_path,
                           skip_sheets = c("Content", "Changes", "Notes",
                                           "Guidance", "Dosages",
                                           "Technical uncertainty",
                                           "Non-species related breakpoints")) {
  all_sheets <- xlsx_sheet_names(xlsx_path)
  
  # Some workbooks (e.g. the antifungal AFST tables) number-prefix their
  # sheet names, e.g. "1. Notes", "7. Dosages", "6. Aspergillus ", and are
  # not always consistent in capitalisation or separators, e.g.
  # "3. Technical Uncertainty" vs "Technical uncertainty", or
  # "PK_PD breakpoints" (v14.0) vs "PK PD breakpoints". Matching skip_sheets
  # by exact equality would fail to skip these, so match case-insensitively
  # (with underscores read as spaces) by whether a skip name occurs anywhere
  # in the (trimmed) sheet name instead.
  is_skipped <- vapply(gsub("_", " ", trimws(all_sheets)), function(s) {
    any(vapply(skip_sheets, function(skip) grepl(tolower(skip), tolower(s), fixed = TRUE),
               logical(1)))
  }, logical(1))
  # Hidden sheets are not part of the published tables: v14.0-v16.1 contain a
  # hidden, outdated copy of the v14.0 PK-PD table with editorial comments
  # ("Remove sheet from BP table"), which must not be read
  hidden <- xlsx_hidden_sheets(xlsx_path)
  if (length(hidden) > 0) {
    message("  skipping hidden sheet(s): ", paste(hidden, collapse = ", "))
  }
  data_sheets <- all_sheets[!is_skipped & !all_sheets %in% hidden]
  
  wb_notes <- parse_workbook_notes(xlsx_path)
  
  results <- list()
  for (s in data_sheets) {
    message("Parsing ", basename(xlsx_path), ": ", s)
    # Errors are deliberately re-thrown, not absorbed: the parsers stop on
    # anything that does not match the expected source structure, and
    # absorbing that here would silently drop the whole sheet instead (as
    # happened before with the Staphylococcus sheet of v9.0 and v10.0).
    res <- tryCatch(
      {
        if (trimws(s) == "Topical agents") {
          parse_topical_sheet(xlsx_path, abbreviations = wb_notes$abbreviations)
        } else if (trimws(s) %in% c("5. Yeast", "6. Aspergillus")) {
          parse_yeast_sheet(xlsx_path, s, abbreviations = wb_notes$abbreviations)
        } else {
          parse_sheet(xlsx_path, s, abbreviations = wb_notes$abbreviations)
        }
      },
      error = function(e) {
        stop("Parsing '", basename(xlsx_path), "', sheet '", s, "' failed: ", conditionMessage(e), call. = FALSE)
      }
    )
    if (!is.null(res) && nrow(res) > 0) {
      res$sheet <- s
      results[[s]] <- res
      message("  -> ", nrow(res), " rows")
    }
  }
  
  out <- bind_rows(results)
  
  # Add EUCAST's own definition of breakpoints in brackets to the note of
  # every bracketed row, so the caveat travels with the values in the raw
  # data (they are excluded from the clinical_breakpoints table).
  if (any(out$is_bracketed)) {
    if (is.na(wb_notes$bracket_definition)) {
      stop("parse_workbook(): '", basename(xlsx_path), "' has bracketed breakpoints, ",
           "but no definition of breakpoints in brackets was found on its Notes sheet")
    }
    bracket_note <- paste0("[( )] ", wb_notes$bracket_definition)
    out$note[out$is_bracketed] <- ifelse(is.na(out$note[out$is_bracketed]),
                                         bracket_note,
                                         paste(bracket_note, out$note[out$is_bracketed], sep = " | "))
  }
  out
}

# ==============================================================================
# 11. Run
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

if (!is.null(parse_log$disk_dose_inherited)) {
  message("NOTE: ", nrow(parse_log$disk_dose_inherited), " disk breakpoint row(s) had no disk content and took it ",
          "from the previous row of the same agent:")
  print(parse_log$disk_dose_inherited)
}

saveRDS(breakpoints_eucast_raw, "data-raw/breakpoints_eucast_raw.rds")
write.csv(breakpoints_eucast_raw, "data-raw/breakpoints_eucast_raw.csv", row.names = FALSE, eol = "\n")


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

# Occasionally, split_agent_organism() fails to separate an organism
# restriction from the antibiotic name in column A -- typically when
# the cell has no rich-text formatting at all (plain text, no bold/non-
# bold boundary to split on), leaving the full string in `ab` while `mo`
# carries only the sheet's default organism. The observed cases all
# follow a pattern of "DrugName, organism-restriction", where the text
# after the comma is recognisably an organism (genus-like capitalised
# word, or "other than" / "except" phrasing) rather than a drug-name
# continuation. For example, EUCAST 2024 Enterococcus sheet:
#   ab = "Vancomycin, enterococci other than E. casseliflavus and E. gallinarum"
# This function detects such cases and returns a two-column tibble
# (ab_clean, mo_from_ab) that the pipeline uses to correct both columns.
split_organism_from_ab <- function(ab, mo) {
  # Only act on ab values that contain a bare comma (i.e. not inside
  # parentheses, which are route/indication qualifiers like "(endocarditis,
  # in combination with...)" that are part of the drug name).
  # Strip parenthetical content first to test for a bare comma.
  ab_no_parens <- gsub("\\([^)]*\\)", "", ab)
  has_bare_comma <- grepl(",", ab_no_parens)
  
  ab_clean <- ab
  mo_from_ab <- rep(NA_character_, length(ab))
  
  for (i in which(has_bare_comma)) {
    # Split at the first bare comma (outside parentheses)
    # Simple approach: find the first comma in ab_no_parens, use its
    # position to split the original ab string
    comma_pos <- regexpr(",", ab_no_parens[i])
    if (comma_pos < 1) next
    before <- trimws(substr(ab[i], 1, comma_pos - 1))
    after  <- trimws(substr(ab[i], comma_pos + 1, nchar(ab[i])))
    if (!nzchar(after)) next
    
    # The part after the comma must look like an organism restriction,
    # not a drug-name continuation. Heuristic: starts with a genus-like
    # capitalised word (3+ letters), or starts with "other" / "except",
    # or contains known taxonomic suffixes (cocci, bacill, etc.)
    is_organism <- grepl(
      "^(other\\b|except\\b|[a-z]*cocc|[a-z]*bacill|[A-Z][a-z]{2,}\\b)",
      after, ignore.case = FALSE, perl = TRUE
    )
    if (!is_organism) next
    
    ab_clean[i] <- before
    mo_from_ab[i] <- after
  }
  
  tibble(ab_clean = ab_clean, mo_from_ab = mo_from_ab)
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

# Route of administration as part of the agent name, e.g. "Cefuroxime iv",
# "Amoxicillin-clavulanic acid oral (uncomplicated UTI only)" or "Ampicillin
# iv (all indications)". Returns "Intravenous", "Oral" or NA.
extract_route <- function(ab) {
  base <- trimws(sub("\\s*\\(.*$", "", ab))
  case_when(
    grepl("\\s(iv|i\\.v\\.)$", base, ignore.case = TRUE) ~ "Intravenous",
    grepl("\\soral$", base, ignore.case = TRUE) ~ "Oral",
    TRUE ~ NA_character_
  )
}

# EUCAST's indication qualifiers mapped to the `site` vocabulary already used
# in clinical_breakpoints (originating from WHONET, e.g. "Non-meningitis",
# "Uncomplicated urinary tract infection" or "Screen"), so that EUCAST rows
# from these tables and from older guidelines are labelled alike. NA means
# the qualifier does not restrict the indication, or is not an indication at
# all: "framycetin" ("Neomycin (framycetin)"), "pivmecillinam" ("Mecillinam
# oral (pivmecillinam)") and "for polymyxin B"
# ("Colistin (for polymyxin B)") name an agent, and "surrogate agent"
# ("Benzylpenicillin (surrogate agent)") describes the agent's use.
# Screening tests become "Screen", which as.sir() leaves out unless
# `include_screening = TRUE`. A qualifier missing from this list stops the
# script, so a new EUCAST wording is never silently taken over as a site.
eucast_site_vocabulary <- c(
  "uncomplicated uti only" = "Uncomplicated urinary tract infection",
  "uti only" = "Urinary tract infection",
  "infections originating from the urinary tract" = "Infections originating from the urinary tract",
  "meningitis" = "Meningitis",
  "indications other than meningitis" = "Non-meningitis",
  "endocarditis" = "Endocarditis",
  "indications other than endocarditis" = "Non-endocarditis",
  "endocarditis and meningitis" = "Meningitis, Endocarditis",
  "indications other than endocarditis and meningitis" = "Non-meningitis, Non-endocarditis",
  "endocarditis, in combination with other antimicrobial treatment" = "Endocarditis with combination treatment",
  "pneumonia" = "Pneumonia",
  "indications other than pneumonia" = "Non-pneumonia",
  "community-acquired pneumonia" = "Community-acquired pneumonia",
  "skin and skin structure infections" = "Skin",
  "systemic infections" = "Systemic infections",
  "other indications" = "Other indications",
  "prophylaxis only" = "Prophylaxis",
  "for prophylaxis only" = "Prophylaxis",
  "screen" = "Screen",
  "screen only" = "Screen",
  "test for high-level aminoglycoside resistance" = "Screen",
  "test for high-level streptomycin resistance" = "Screen",
  "test for acquired aminoglycoside-modifying enzyme" = "Screen",
  "all indications" = NA_character_,
  "all indications including prophylaxis" = NA_character_,
  "all indications, including meningitis and prophylaxis" = NA_character_,
  "surrogate agent" = NA_character_,
  "framycetin" = NA_character_,
  "pivmecillinam" = NA_character_, # "Mecillinam oral (pivmecillinam)", the oral prodrug
  "for polymyxin b" = NA_character_
)

# Builds `site` from the sheet-level site (only "Topical" for the Topical
# agents sheet), the route and the indication qualifier, in the established
# "Route, Indication" form, e.g. "Oral, Uncomplicated urinary tract
# infection". The values on the Topical agents sheet are, in EUCAST's own
# words, "Screening cut-off values" (EUCAST "has not been able to determine
# relevant clinical breakpoints for topical use"), so they get "Topical,
# Screen" and are, like other screening values, only used by as.sir() with
# `include_screening = TRUE`.
build_site <- function(sheet_site, route, qualifier) {
  q <- tolower(gsub("\\s+", " ", trimws(qualifier)))
  unknown <- unique(qualifier[!is.na(q) & !q %in% names(eucast_site_vocabulary)])
  if (length(unknown) > 0) {
    stop("build_site(): unknown indication qualifier(s), add to eucast_site_vocabulary: '",
         paste(unknown, collapse = "', '"), "'")
  }
  indication <- unname(eucast_site_vocabulary[q])
  vapply(seq_along(q), function(i) {
    parts <- c(sheet_site[i], route[i], indication[i])
    if (identical(sheet_site[i], "Topical")) parts <- c(parts, "Screen")
    parts <- unique(parts[!is.na(parts)])
    if (length(parts) == 0) NA_character_ else paste(parts, collapse = ", ")
  }, character(1))
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

# Descriptive organism groups that as.mo() cannot resolve correctly on its
# own, mapped explicitly, following how the WHONET-based breakpoints encode
# them:
#  - "Gram-negative anaerobes" and "Gram-positive anaerobes" (v9.0-v11.0)
#    would resolve to B_GRAMN and B_GRAMP, which as.sir() applies to every
#    Gram-negative or Gram-positive isolate, aerobes included. B_ANAER-NEG
#    and B_ANAER-POS only apply to anaerobes.
#  - "S. anginosus group" would resolve to the subspecies S. anginosus
#    subsp. anginosus. There is no code for the group, so it becomes its
#    three species.
#  - "Streptococcus groups A, C and G" would resolve to groups A, B, C and G,
#    while EUCAST deliberately excludes group B in such rows.
#  - "MRSA" (telavancin) becomes S. aureus, as there is no taxon for MRSA.
#  - "enterococci" is the genus Enterococcus (its exceptions are handled by
#    rank_index, see above).
eucast_organism_groups <- list(
  "gram-negative anaerobes" = "B_ANAER-NEG",
  "gram-positive anaerobes" = "B_ANAER-POS",
  "s. anginosus group" = c("B_STRPT_ANGN", "B_STRPT_CNST", "B_STRPT_INTR"),
  "streptococcus anginosus group" = c("B_STRPT_ANGN", "B_STRPT_CNST", "B_STRPT_INTR"),
  "streptococcus groups a, c and g" = c("B_STRPT_GRPA", "B_STRPT_GRPC", "B_STRPT_GRPG"),
  "mrsa" = "B_STPHY_AURS",
  # "Vancomycin, enterococci other than E. casseliflavus and E. gallinarum"
  "enterococci" = "B_ENTRC",
  # the non-species related PK-PD breakpoints (v9.0-v13.1), which as.sir()
  # uses as the last resort for organism code "UNKNOWN"
  "pk-pd (non-species related) breakpoints" = "UNKNOWN"
)
expand_organism_groups <- function(x) {
  lapply(x, function(txt) {
    key <- tolower(gsub("\\s+", " ", trimws(txt)))
    if (key %in% names(eucast_organism_groups)) eucast_organism_groups[[key]] else txt
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
  # "spp." always ends a taxon, so a capitalised genus directly following it
  # starts the next list member even if the separating comma is missing, e.g.
  # v12.0 Mecillinam: "E. coli, Citrobacter spp.<line break>Klebsiella spp., ..."
  s_clean <- gsub("\\bspp\\.\\s+(?=[A-Z][a-z]{2,})", "spp., ", s_clean, perl = TRUE)
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

# Blocking rows ----
# Where EUCAST gives no breakpoint for an organism, as.sir() must not fall back
# to the breakpoint of a broader taxon (or to the non-species related PK-PD
# breakpoints), e.g. Vibrio fluvialis and cefotaxime ("IE") would otherwise be
# interpreted with the breakpoints for Vibrio spp. Such cases are therefore
# kept as "blocking rows": breakpoint_S and breakpoint_R are NA, so as.sir()
# returns NA, and the note starts with "[No breakpoint]" and the reason. A
# blocking row is only kept where the organism has no breakpoint of its own
# for that agent and method (see "Pruning blocking rows" below), and its
# rank_index is that of its organism + 0.1, so that a breakpoint for the same
# organism (e.g. for another route, or a screening breakpoint when
# as.sir(..., include_screening = TRUE)) always takes precedence. Blocking rows
# arise from:
#  1. cells with "IE", "-", "Note", "IP", "NA" or "ND", and EUCAST's
#     breakpoints in brackets (see is_bracketed_bp());
#  2. organisms excluded with "except" or "other than", e.g. "Bacillus spp.
#     except B. anthracis" or "Klebsiella spp. (except K. aerogenes)";
#  3. restrictions in the notes, see eucast_note_rules;
#  4. the organism scope of a table: an agent listed for only some organisms
#     of a table (e.g. phenoxymethylpenicillin for "Streptococcus groups A, C
#     and G" on the sheet for groups A, B, C and G) has no breakpoint for the
#     other organisms of that table;
#  5. organisms for which EUCAST has not determined breakpoints at all (B.
#     cepacia complex, L. pneumophila), for the years with PK-PD breakpoints.
# The Topical agents sheet never yields blocking rows: that a topical
# screening cut-off value is lacking says nothing about systemic use.
bp_placeholder_meaning <- c(
  "IE" = "IE (insufficient evidence)",
  "NOTE" = "Note (see note)",
  "-" = "- (susceptibility testing is not recommended)",
  "IP" = "IP (in preparation)",
  "NA" = "NA (not applicable)",
  "ND" = "ND (not determined)"
)
blocking_note <- function(reason, note) {
  reason <- paste0("[No breakpoint] ", reason)
  note <- rep_len(as.character(note), length(reason))
  if_else(is.na(note), reason, paste(reason, note, sep = " | "))
}

breakpoints_eucast <- breakpoints_eucast_raw |>
  mutate(
    is_blocker = is_bracketed | (!is_bp_value(breakpoint_S) & !is_bp_value(breakpoint_R)),
    note = case_when(
      is_bracketed ~ blocking_note(paste0("EUCAST gives breakpoints in brackets only (S ", breakpoint_S, ", R ", breakpoint_R, ")"), note),
      is_blocker ~ blocking_note(paste0("EUCAST gives ",
                                        coalesce(bp_placeholder_meaning[toupper(trimws(breakpoint_S))],
                                                 bp_placeholder_meaning[toupper(trimws(breakpoint_R))],
                                                 "no value")), note),
      TRUE ~ note
    ),
    # placeholders become NA; a row keeps a breakpoint as long as one of S and
    # R is a value, e.g. the screening tests of v9.0-v12.0 that give only an S
    # value with "Note" for R
    breakpoint_S = if_else(is_blocker | !is_bp_value(breakpoint_S), NA_character_, breakpoint_S),
    breakpoint_R = if_else(is_blocker | !is_bp_value(breakpoint_R), NA_character_, breakpoint_R)
  ) |>
  filter(!(is_blocker & sheet == "Topical agents"))

# Occasionally, split_agent_organism() could not separate an organism
# restriction from the drug name (see split_organism_from_ab comment).
# Detect and fix those before mo resolution, so the organism restriction
# is processed alongside the regular mo text.
ab_mo_fix <- split_organism_from_ab(breakpoints_eucast$ab, breakpoints_eucast$mo)
fixed_rows <- !is.na(ab_mo_fix$mo_from_ab)
if (any(fixed_rows)) {
  message("NOTE: ", sum(fixed_rows), " row(s) had an organism restriction embedded in the ab column; ",
          "splitting into ab + mo (e.g. '", breakpoints_eucast$ab[which(fixed_rows)[1]], "').")
  breakpoints_eucast$ab[fixed_rows] <- ab_mo_fix$ab_clean[fixed_rows]
  breakpoints_eucast$mo[fixed_rows] <- ab_mo_fix$mo_from_ab[fixed_rows]
}

# "MRSA" (telavancin) can only be coded as S. aureus, so state the scope
breakpoints_eucast <- breakpoints_eucast |>
  mutate(note = case_when(
    trimws(gsub("\\s+", " ", mo)) != "MRSA" ~ note,
    is_blocker ~ trimws(paste0(trimws(sub("^(\\[No breakpoint\\][^|]*).*$", "\\1", note)), " | [Scope] EUCAST gives this for MRSA ",
                               sub("^\\[No breakpoint\\][^|]*", "", note))),
    is.na(note) ~ "[Scope] EUCAST gives this breakpoint for MRSA",
    TRUE ~ paste("[Scope] EUCAST gives this breakpoint for MRSA", note, sep = " | ")
  ))

# Exceptions ("except" / "other than") become blocking rows for the excluded
# organisms, e.g. "Bacillus spp. except B. anthracis" gives a blocking row for
# B. anthracis. An abbreviated genus in the exception is expanded with the
# genus of the organism it is an exception to.
plural_genus <- c("staphylococci" = "Staphylococcus", "enterococci" = "Enterococcus", "streptococci" = "Streptococcus")
extract_mo_exceptions <- function(x) {
  lapply(x, function(txt) {
    txt <- gsub("\\s+", " ", trimws(txt))
    out <- character(0)
    # parenthesised exception attached to one list member, e.g.
    # "Klebsiella spp. (except K. aerogenes)"
    m <- gregexpr("([A-Z][a-z]+)[^,(]*\\(\\s*(?:except|other than)\\s+([^)]+)\\)", txt, perl = TRUE, ignore.case = TRUE)
    for (hit in regmatches(txt, m)[[1]]) {
      genus <- sub("^([A-Z][a-z]+).*$", "\\1", hit)
      exc <- sub("^.*\\(\\s*(?:except|other than)\\s+([^)]+)\\).*$", "\\1", hit, perl = TRUE, ignore.case = TRUE)
      out <- c(out, paste0(genus, "|", exc))
    }
    txt_rest <- gsub("\\([^)]*\\)", "", txt)
    m2 <- regmatches(txt_rest, regexec("^(.*?)\\b(?:except|other than)\\b\\s*(.+)$", txt_rest, perl = TRUE, ignore.case = TRUE))[[1]]
    if (length(m2) == 3) {
      base <- m2[2]
      genus <- regmatches(base, regexpr("[A-Z][a-z]+(?=\\s+spp\\.?|\\s*$|\\s+other|\\s+except)", base, perl = TRUE))
      if (length(genus) == 0) {
        plural <- names(plural_genus)[vapply(names(plural_genus), function(p) grepl(p, base, ignore.case = TRUE), logical(1))]
        genus <- if (length(plural) > 0) unname(plural_genus[plural[1]]) else NA_character_
      }
      out <- c(out, paste0(coalesce(genus[1], ""), "|", m2[3]))
    }
    # split each exception list into taxa and expand abbreviated genera
    unlist(lapply(out, function(o) {
      genus <- sub("\\|.*$", "", o)
      exc <- sub("^[^|]*\\|", "", o)
      pieces <- trimws(unlist(strsplit(exc, ",|\\s+and\\s+")))
      pieces <- pieces[nzchar(pieces)]
      if (nzchar(genus)) {
        pieces <- ifelse(grepl("^[A-Z]\\.\\s*[a-z]", pieces) & substr(pieces, 1, 1) == substr(genus, 1, 1),
                         sub("^[A-Z]\\.\\s*", paste0(genus, " "), pieces), pieces)
      }
      pieces
    }))
  })
}
exceptions <- breakpoints_eucast |>
  filter(sheet != "Topical agents") |>
  mutate(exception = extract_mo_exceptions(mo)) |>
  filter(lengths(exception) > 0) |>
  tidyr::unnest(exception)
if (nrow(exceptions) > 0) {
  exception_blockers <- exceptions |>
    mutate(note = blocking_note(paste0("EUCAST excludes ", exception, " ('", gsub("\\s+", " ", mo), "')"), NA_character_),
           mo = exception, site = NA_character_, is_blocker = TRUE,
           breakpoint_S = NA_character_, breakpoint_R = NA_character_) |>
    select(-exception)
  breakpoints_eucast <- bind_rows(breakpoints_eucast, exception_blockers)
}

breakpoints_eucast <- breakpoints_eucast |>
  mutate(
    mo_text = strip_mo_exceptions_wholecell(mo),
    mo_text = resolve_other_x(mo_text, sheet),
    mo_text = split_mo_list(mo_text, sheet_name = sheet)
  ) |>
  tidyr::unnest(mo_text) |>
  mutate(mo_text = expand_organism_groups(mo_text)) |>
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

# A descriptive organism group (e.g. "... group", "... anaerobes", "MRSA")
# must either be mapped in eucast_organism_groups or resolve to a species
# group code (such as Viridans group streptococci or coagulase-negative
# staphylococci). Otherwise as.mo() has matched it to a single taxon that
# does not cover the group, as happened with "S. anginosus group".
unmapped_groups <- breakpoints_eucast |>
  distinct(mo_text, mo) |>
  filter(grepl("group|anaerob|MRSA|MSSA|-negative|-positive|cocci\\b|bacteria\\b", mo_text, ignore.case = TRUE),
         !grepl("spp\\.$", mo_text),
         # a text that is exactly the name of the taxon it resolved to, such
         # as "Peptostreptococcus anaerobius", is not a descriptive group
         tolower(mo_text) != tolower(mo_name(mo, keep_synonyms = TRUE)),
         !mo_rank(mo, keep_synonyms = TRUE) %in% "species group",
         !mo %in% unlist(eucast_organism_groups))
if (nrow(unmapped_groups) > 0) {
  print(unmapped_groups |> mutate(name = mo_name(mo, keep_synonyms = TRUE)))
  stop(nrow(unmapped_groups), " descriptive organism group(s) resolved to a single taxon, add them to eucast_organism_groups")
}
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
    # The route and indication (if any) are captured into `site` before they
    # are stripped from `ab`, e.g. "Cefuroxime oral (uncomplicated UTI
    # only)" -> ab = "Cefuroxime", site = "Oral, Uncomplicated urinary tract
    # infection" (see build_site()). This keeps the iv/oral or indication-
    # specific variants of the same drug+organism+method distinguishable,
    # both for the distinct() step below (site is part of its key) and for
    # as.sir(), which would otherwise see e.g. "Cefuroxime iv" and
    # "Cefuroxime oral" as two conflicting breakpoints for the same drug.
    qualifier = extract_trailing_parenthetical(ab),
    route = extract_route(ab),
    site = build_site(site, route, qualifier),
    # also fixes a source cell missing a space before "(", e.g. "Meropenem(indications..."
    ab_input = gsub("\\s*\\(.*$", "", ab),
    ab_input = trimws(sub("\\s+(iv|i\\.v\\.|oral)$", "", ab_input, ignore.case = TRUE)),
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

# Restrictions in the notes ----
# Notes can restrict a breakpoint to some organisms, methods or indications,
# e.g. "Zone diameter breakpoints apply to E. coli only. For other
# Enterobacterales, use an MIC method." Every sentence in the notes with such
# wording must match exactly one of these rules, either with an action or as
# informational (e.g. advice on the route of administration, which does not
# change the interpretation of an MIC or zone diameter). A sentence matching
# no rule stops the script, so that a new restriction is never missed.
#  - block:      the governed rows of `method` become blocking rows
#  - site:       `target` is added to the site of the governed rows
#  - restrict:   the governed rows of `method` apply to organism `target`
#                only: a row for a broader taxon containing `target` is set to
#                `target` (with a blocking row for the broader taxon), a row
#                for any other organism becomes a blocking row
#  - exclude:    organism `target` gets a blocking row for `method`
#  - must_block: the governed rows of `method` must already be blocking rows
#                (e.g. "Use an MIC method" for zone diameters given as "Note")
#  - none:       informational
eucast_note_rules <- tibble::tribble(
  ~pattern, ~action, ~method, ~target,
  "zone diameter breakpoint is valid for MSSA only", "block", "DISK", NA,
  "Breakpoints relate to topical use only", "block", NA, NA,
  "Breakpoints apply only to use in the prophylaxis of meningococcal disease", "site", NA, "Prophylaxis",
  "For prophylaxis of meningitis only", "site", NA, "Prophylaxis",
  "Zone diameter breakpoints (apply to|validated for) E\\. coli only", "restrict", "DISK", "B_ESCHR_COLI",
  "For C\\. ?koseri, use an MIC method", "exclude", "DISK", "B_CTRBC_KOSR",
  "For other Enterobacterales, use an MIC method", "none", NA, NA, # see the E. coli restriction
  "^Use an MIC method", "must_block", "DISK", NA,
  "Disk diffusion is unreliable", "must_block", "DISK", NA,
  "^May be tested for epidemiological purposes only \\(ECOFF", "must_block", NA, NA,
  "The .IE. in the table pertains to intravenous therapy only", "none", NA, NA,
  "ATU relevant only if", "none", NA, NA,
  "have been determined for meropenem only", "none", NA, NA,
  "Meropenem is the only carbapenem used for meningitis", "none", NA, NA,
  "Not for meningitis \\(meropenem is the only carbapenem", "none", NA, NA,
  "Breakpointss? apply only to tests performed on Middlebrook", "none", NA, NA,
  "Combinations between penicillins or glycopeptides and aminoglycosides", "none", NA, NA,
  "EUCAST does not recommend systematic screening for BORSA", "none", NA, NA,
  "(For oral administration the breakpoints|Breakpoints for oral administration|Oral administration) (are|is) relevant for (uncomplicated )?urinary tract infections only", "none", NA, NA,
  "For staphylococci other than S\\. aureus, S\\. lugdunensis and S\\. saprophyticus, the cefoxitin MIC is a poorer predictor", "none", NA, NA,
  "only a small number of cases involving these species", "none", NA, NA,
  "Isolates must not be reported susceptible before 24 h incubation", "none", NA, NA,
  "^(Positive|Negative) test:", "none", NA, NA, # see the derived screening breakpoints
  "^The isolate is high-level resistant to gentamicin and other aminoglycosides", "none", NA, NA,
  "Read and interpret the benzylpenicillin disk only for isolates with oxacillin", "none", NA, NA,
  "^Susceptibility inferred from ampicillin", "none", NA, NA,
  "Susceptibility of staphylococci to cephalosporins is inferred from the cefoxitin susceptibility", "none", NA, NA,
  "The breakpoints also apply to meningitis", "none", NA, NA,
  "apply to oral treatment of C\\. difficile infections", "none", NA, NA,
  "The susceptibility of streptococcus groups A, B, C and G to penicillins is inferred", "none", NA, NA,
  "Vancomycin resistance is the expected phenotype for E\\. casseliflavus and E\\. gallinarum", "none", NA, NA,
  "^When the screen is negative", "none", NA, NA,
  "^These breakpoints are used only when there are no species-specific breakpoints", "none", NA, NA, # PK-PD, see parse_sheet()
  "^Include a note that the guidance is based on PK-PD breakpoints only", "none", NA, NA # PK-PD, reporting advice
)
restriction_wording <- "\\bonly\\b|do(es)? not apply|not applicable|not valid|apply to|applies to|\\bexcept\\b|other than|should not be|must not be|cannot be used|not be used|are not reliable|unreliable|not recommended|do not use|use an MIC method|exclusively"
note_sentences <- function(note) {
  bodies <- sub("^\\[[^]]+\\]\\s*", "", unlist(strsplit(note, " \\| ")))
  trimws(unlist(strsplit(gsub("\\s+", " ", bodies), "(?<=[a-z0-9)])\\.\\s+(?=[A-Z])", perl = TRUE)))
}
# (the source notes, i.e. without the "[No breakpoint]" reasons added above)
all_sentences <- unique(note_sentences(breakpoints_eucast_raw$note[!is.na(breakpoints_eucast_raw$note)]))
restriction_sentences <- all_sentences[grepl(restriction_wording, all_sentences, ignore.case = TRUE)]
rule_hits <- vapply(restriction_sentences, function(x) sum(vapply(eucast_note_rules$pattern, function(p) grepl(p, x, perl = TRUE), logical(1))), integer(1))
if (any(rule_hits != 1)) {
  print(data.frame(rules_matched = rule_hits[rule_hits != 1], sentence = substr(restriction_sentences[rule_hits != 1], 1, 200)))
  stop("note sentence(s) with restriction wording that match no rule (or more than one) in eucast_note_rules, see above")
}

taxon_contains <- function(broad, narrow) {
  # TRUE if taxon `narrow` belongs to taxon or species group `broad`
  narrow <- as.mo(narrow, info = FALSE)
  broad <- as.character(broad)
  broad == as.character(narrow) |
    broad == as.character(as.mo(mo_genus(narrow, keep_synonyms = TRUE), info = FALSE)) |
    broad == as.character(as.mo(mo_family(narrow, keep_synonyms = TRUE), info = FALSE)) |
    broad == as.character(as.mo(mo_order(narrow, keep_synonyms = TRUE), info = FALSE)) |
    broad %in% AMR::microorganisms.groups$mo_group[AMR::microorganisms.groups$mo == as.character(narrow)]
}
to_blocker <- function(df, reason) {
  df |> mutate(is_blocker = TRUE, breakpoint_S = NA_character_, breakpoint_R = NA_character_,
               site = NA_character_, note = blocking_note(reason, note))
}
for (i in seq_len(nrow(eucast_note_rules))) {
  rule <- eucast_note_rules[i, ]
  if (rule$action == "none") next
  governed <- !is.na(breakpoints_eucast$note) &
    vapply(breakpoints_eucast$note, function(n) any(grepl(rule$pattern, note_sentences(n), perl = TRUE)), logical(1)) &
    (is.na(rule$method) | breakpoints_eucast$method == rule$method)
  if (!any(governed)) next
  # the reason quotes the matching sentence of each row's own note
  reason_of <- function(rows) {
    vapply(breakpoints_eucast$note[rows], function(n) {
      paste0("EUCAST: '", note_sentences(n)[grepl(rule$pattern, note_sentences(n), perl = TRUE)][1], "'")
    }, character(1), USE.NAMES = FALSE)
  }
  if (rule$action == "must_block") {
    if (any(!breakpoints_eucast$is_blocker[governed])) {
      print(breakpoints_eucast[governed & !breakpoints_eucast$is_blocker, c("guideline", "sheet", "mo", "ab", "method", "breakpoint_S", "breakpoint_R")])
      stop("rows governed by '", rule$pattern, "' have breakpoints, see above")
    }
  } else if (rule$action == "block") {
    rows <- which(governed & !breakpoints_eucast$is_blocker)
    breakpoints_eucast[rows, ] <- to_blocker(breakpoints_eucast[rows, ], reason_of(rows))
  } else if (rule$action == "site") {
    breakpoints_eucast$site[governed] <- if_else(is.na(breakpoints_eucast$site[governed]), rule$target,
                                                 if_else(grepl(rule$target, breakpoints_eucast$site[governed], fixed = TRUE),
                                                         breakpoints_eucast$site[governed],
                                                         paste(breakpoints_eucast$site[governed], rule$target, sep = ", ")))
  } else if (rule$action == "restrict") {
    rows <- which(governed & !breakpoints_eucast$is_blocker)
    for (r in rows) {
      if (as.character(breakpoints_eucast$mo[r]) == rule$target) next
      if (taxon_contains(breakpoints_eucast$mo[r], rule$target)) {
        breakpoints_eucast <- bind_rows(breakpoints_eucast, to_blocker(breakpoints_eucast[r, ], reason_of(r)))
        breakpoints_eucast$mo[r] <- as.mo(rule$target)
      } else {
        breakpoints_eucast[r, ] <- to_blocker(breakpoints_eucast[r, ], reason_of(r))
      }
    }
  } else if (rule$action == "exclude") {
    rows <- which(governed & !breakpoints_eucast$is_blocker)
    for (r in rows) {
      if (as.character(breakpoints_eucast$mo[r]) == rule$target) {
        breakpoints_eucast[r, ] <- to_blocker(breakpoints_eucast[r, ], reason_of(r))
      } else if (taxon_contains(breakpoints_eucast$mo[r], rule$target)) {
        new_row <- to_blocker(breakpoints_eucast[r, ], reason_of(r))
        new_row$mo <- as.mo(rule$target)
        breakpoints_eucast <- bind_rows(breakpoints_eucast, new_row)
      }
    }
  } else {
    stop("unknown action in eucast_note_rules: ", rule$action)
  }
}

# Screening breakpoints stated in the notes ----
# A few screening tests have their breakpoints in the note instead of the
# table cell ("Note"). Only explicit, decisive thresholds are taken over, and
# each is read from the note of the row itself (so it follows the version):
#  1. High-level aminoglycoside resistance in enterococci and viridans group
#     streptococci: "Negative test: Isolates with gentamicin MIC <=128 mg/L or
#     a zone diameter >=8 mm" / "Positive test: ... MIC >128 mg/L or a zone
#     diameter <8 mm", likewise for streptomycin. These are coded as the
#     agents gentamicin-high (GEH) and streptomycin-high (STH), as before.
#  2. Methicillin resistance by cefoxitin MIC in staphylococci: "S. aureus and
#     S. lugdunensis with cefoxitin MIC values >4 mg/L and S. saprophyticus
#     with cefoxitin MIC values >8 mg/L are methicillin resistant". Like the
#     cefoxitin disk ("screen only"), these get site "Screen".
#  3. Cefoxitin zone diameters for coagulase-negative staphylococci that are
#     not identified to species level: "If coagulase-negative staphylococci
#     are not identified to species level, use zone diameter breakpoints
#     S>=25, R<25 mm". The breakpoint for "S. aureus and coagulase-negative
#     staphylococci except S. epidermidis and S. lugdunensis" then applies to
#     the identified CoNS species (one row per species), while the CoNS group
#     itself (an isolate only identified as CoNS) gets 25 mm.
# Statements that are not breakpoints, such as ECOFFs ("There are no clinical
# breakpoints but acquired resistance (indicated by MIC >1 mg/L) should be
# excluded") or hedged ones ("... with oxacillin MIC values >2 mg/L are mostly
# methicillin resistant"), are not taken over.
derived_rows <- list()
# (only the rows of the test itself, e.g. "Gentamicin (test for high-level
# aminoglycoside resistance)"; other aminoglycosides refer to the same note)
hlar_rows <- breakpoints_eucast |>
  filter(grepl("Negative test: Isolates with (gentamicin|streptomycin) MIC \u2264", note),
         ab %in% c("GEN", "STR1"), site %in% "Screen")
for (r in seq_len(nrow(hlar_rows))) {
  row <- hlar_rows[r, ]
  # the note of a row can also hold the test for the other agent, so the
  # agent is taken from the row itself
  agent <- c(GEN = "gentamicin", STR1 = "streptomycin")[as.character(row$ab)]
  if (is.na(agent)) stop("unexpected agent for a high-level aminoglycoside screening test: ", row$ab)
  neg <- regmatches(row$note, regexec(paste0("Negative test: Isolates with ", agent, " MIC \u2264([0-9]+) mg/L(?: or a zone diameter \u2265([0-9]+) mm)?"), row$note, perl = TRUE))[[1]]
  pos <- regmatches(row$note, regexec(paste0("Positive test: Isolates with ", agent, " MIC >([0-9]+) mg/L(?: or a zone diameter <([0-9]+) mm)?"), row$note, perl = TRUE))[[1]]
  if (length(neg) != 3 || length(pos) != 3 || neg[2] != pos[2] || neg[3] != pos[3]) {
    stop("inconsistent high-level aminoglycoside screening test in the note of ", row$guideline, " ", row$sheet, ": ", substr(row$note, 1, 300))
  }
  value <- if (row$method == "MIC") neg[2] else neg[3]
  if (!nzchar(value)) next # e.g. no zone diameter given in this version
  derived_rows[[length(derived_rows) + 1]] <- row |>
    mutate(ab = as.ab(if (agent == "gentamicin") "GEH" else "STH"), site = NA_character_, is_blocker = FALSE,
           breakpoint_S = value, breakpoint_R = value,
           note = sub("^\\[No breakpoint\\][^|]*(\\| )?", "", note))
}
fox_mic <- breakpoints_eucast |>
  filter(ab == "FOX", method == "MIC",
         grepl("S\\. aureus and S\\. lugdunensis with cefoxitin MIC values >[0-9]+ mg/L and S\\. saprophyticus with cefoxitin MIC values >[0-9]+ mg/L are methicillin resistant", note)) |>
  distinct(guideline, .keep_all = TRUE)
for (r in seq_len(nrow(fox_mic))) {
  row <- fox_mic[r, ]
  m <- regmatches(row$note, regexec("S\\. aureus and S\\. lugdunensis with cefoxitin MIC values >([0-9]+) mg/L and S\\. saprophyticus with cefoxitin MIC values >([0-9]+) mg/L", row$note))[[1]]
  for (sp in list(c("B_STPHY_AURS", m[2]), c("B_STPHY_LGDN", m[2]), c("B_STPHY_SPRP", m[3]))) {
    derived_rows[[length(derived_rows) + 1]] <- row |>
      mutate(mo = as.mo(sp[1]), site = "Screen", is_blocker = FALSE, breakpoint_S = sp[2], breakpoint_R = sp[2],
             note = sub("^\\[No breakpoint\\][^|]*(\\| )?", "", note))
  }
}
if (length(derived_rows) > 0) {
  derived <- bind_rows(derived_rows)
  message("NOTE: ", nrow(derived), " screening breakpoint row(s) taken from the notes (",
          paste(unique(paste(derived$guideline, derived$ab, derived$method)), collapse = ", "), ")")
  breakpoints_eucast <- bind_rows(breakpoints_eucast, derived)
}
cons_rule <- "If coagulase-negative staphylococci are not identified to species level, use zone diameter breakpoints S\u2265([0-9]+), R<([0-9]+) mm"
for (g in unique(breakpoints_eucast$guideline[grepl(cons_rule, breakpoints_eucast$note) & breakpoints_eucast$ab == "FOX" & breakpoints_eucast$method == "DISK"])) {
  rows <- which(breakpoints_eucast$guideline == g & breakpoints_eucast$ab == "FOX" & breakpoints_eucast$method == "DISK" &
                  breakpoints_eucast$mo == "B_STPHY_CONS" & !breakpoints_eucast$is_blocker)
  if (length(rows) != 1) stop("expected one cefoxitin zone diameter row for coagulase-negative staphylococci in ", g)
  m <- regmatches(breakpoints_eucast$note[rows], regexec(cons_rule, breakpoints_eucast$note[rows]))[[1]]
  if (length(m) != 3 || m[2] != m[3]) stop("unexpected cefoxitin rule for unidentified CoNS in ", g)
  own_rows <- unique(as.character(breakpoints_eucast$mo[breakpoints_eucast$guideline == g & breakpoints_eucast$ab == "FOX" &
                                                           breakpoints_eucast$method == "DISK" & breakpoints_eucast$mo != "B_STPHY_CONS"]))
  members <- setdiff(as.character(AMR::microorganisms.groups$mo[AMR::microorganisms.groups$mo_group == "B_STPHY_CONS"]), own_rows)
  species_rows <- breakpoints_eucast[rep(rows, length(members)), ]
  species_rows$mo <- as.mo(members)
  breakpoints_eucast$breakpoint_S[rows] <- m[2]
  breakpoints_eucast$breakpoint_R[rows] <- m[3]
  breakpoints_eucast$note[rows] <- paste0("[Scope] EUCAST: '", m[1], "' | ", breakpoints_eucast$note[rows])
  breakpoints_eucast <- bind_rows(breakpoints_eucast, species_rows)
}

# Organism scope of each table ----
# For every agent and method listed in a table, the organism of the table
# (e.g. Enterobacterales, or Streptococcus groups A, B, C and G) gets a
# blocking row if it has no breakpoint of its own for that agent and method,
# so that an organism of the table that is not listed for that agent (e.g.
# group B for phenoxymethylpenicillin, or Salmonella for cefuroxime iv) does
# not fall back to a broader breakpoint. Not for the Topical agents sheet and
# the PK-PD table.
table_scope <- breakpoints_eucast |>
  filter(sheet != "Topical agents", !grepl("^PK", sheet)) |>
  distinct(guideline, sheet, table_organism) |>
  mutate(table_mo = split_mo_list(resolve_other_x(strip_mo_exceptions_wholecell(table_organism), sheet), sheet_name = sheet)) |>
  tidyr::unnest(table_mo) |>
  mutate(table_mo = expand_organism_groups(table_mo)) |>
  tidyr::unnest(table_mo) |>
  mutate(table_mo = as.mo(table_mo, info = FALSE))
if (any(is.na(table_scope$table_mo) | table_scope$table_mo == "UNKNOWN")) {
  print(table_scope |> filter(is.na(table_mo) | table_mo == "UNKNOWN"))
  stop("table organism(s) could not be resolved, see above")
}
scope_blockers <- breakpoints_eucast |>
  filter(sheet != "Topical agents", !grepl("^PK", sheet)) |>
  distinct(guideline, sheet, table_organism, ab, method, .keep_all = TRUE) |>
  inner_join(table_scope, by = c("guideline", "sheet", "table_organism"), relationship = "many-to-many") |>
  mutate(mo = table_mo, note = NA_character_) |>
  select(-table_mo) |>
  to_blocker("EUCAST lists this agent in this table for other organisms only")
breakpoints_eucast <- bind_rows(breakpoints_eucast, scope_blockers)

# Apply the statements that extend MIC breakpoints to another species (see
# parse_sheet(), "Statements extending MIC breakpoints to another species"):
# the MIC rows of the source species on that sheet are copied for the target
# species, with EUCAST's statement in the note. A breakpoint that the source
# gives for the target species itself always takes precedence.
if (!is.null(parse_log$mic_scope_extensions)) {
  scope_ext <- parse_log$mic_scope_extensions |>
    distinct() |>
    mutate(from_mo = as.character(as.mo(from_text, info = FALSE)),
           to_mo = as.character(as.mo(to_text, info = FALSE)))
  if (anyNA(scope_ext$from_mo) || anyNA(scope_ext$to_mo) || any(c(scope_ext$from_mo, scope_ext$to_mo) == "UNKNOWN")) {
    print(scope_ext)
    stop("species in a MIC scope statement could not be resolved, see above")
  }
  scope_rows <- breakpoints_eucast |>
    filter(method == "MIC") |>
    mutate(mo_chr = as.character(mo)) |>
    inner_join(scope_ext |> select(guideline, sheet, mo_chr = from_mo, to_mo, statement),
               by = c("guideline", "sheet", "mo_chr")) |>
    mutate(mo = as.mo(to_mo),
           note = case_when(
             is_blocker ~ trimws(paste0(trimws(sub("^(\\[No breakpoint\\][^|]*).*$", "\\1", note)), " | [Scope] ", statement, " ",
                                        sub("^\\[No breakpoint\\][^|]*", "", note))),
             is.na(note) ~ paste0("[Scope] ", statement),
             TRUE ~ paste0("[Scope] ", statement, " | ", note)
           )) |>
    select(-mo_chr, -to_mo, -statement) |>
    anti_join(breakpoints_eucast |> filter(method == "MIC"), by = c("guideline", "mo", "ab", "method", "site"))
  message("NOTE: ", nrow(scope_rows), " MIC breakpoint row(s) added following EUCAST's statement(s): ",
          paste(unique(paste0(scope_ext$guideline, " ", scope_ext$from_text, " -> ", scope_ext$to_text)), collapse = ", "))
  breakpoints_eucast <- bind_rows(breakpoints_eucast, scope_rows)
}

# Organisms for which EUCAST has not determined breakpoints at all ("EUCAST
# has not determined breakpoints for Burkholderia cepacia complex organisms",
# likewise L. pneumophila) must not be interpreted with the non-species
# related PK-PD breakpoints either (v9.0-v13.1)
pkpd <- breakpoints_eucast |> filter(grepl("^PK", sheet), !is_blocker)
if (nrow(pkpd) > 0 && !is.null(parse_log$no_breakpoints_determined)) {
  no_bp <- parse_log$no_breakpoints_determined |>
    distinct(guideline, organism, .keep_all = TRUE) |>
    mutate(no_bp_mo = as.mo(organism, info = FALSE))
  if (any(is.na(no_bp$no_bp_mo) | no_bp$no_bp_mo == "UNKNOWN")) {
    print(no_bp)
    stop("organism(s) without EUCAST breakpoints could not be resolved, see above")
  }
  no_bp_blockers <- pkpd |>
    distinct(guideline, ab, method, .keep_all = TRUE) |>
    inner_join(no_bp |> select(guideline, no_bp_mo, statement), by = "guideline", relationship = "many-to-many") |>
    mutate(mo = no_bp_mo, note = NA_character_) |>
    select(-no_bp_mo)
  no_bp_blockers <- to_blocker(no_bp_blockers, paste0("EUCAST: '", no_bp_blockers$statement, "'")) |> select(-statement)
  breakpoints_eucast <- bind_rows(breakpoints_eucast, no_bp_blockers)
}

# Any breakpoint that is not purely numeric (placeholders are NA by now)
# would silently become NA in as.numeric() below, so stop here instead.
is_unexpected_bp <- function(x) {
  !is.na(x) & x %unlike_case% "^[0-9.]+$"
}
non_numeric_bp <- breakpoints_eucast |>
  filter(is_unexpected_bp(breakpoint_S) | is_unexpected_bp(breakpoint_R))
if (nrow(non_numeric_bp) > 0) {
  print(non_numeric_bp |> select(guideline, sheet, mo_text, ab, method, breakpoint_S, breakpoint_R))
  stop(nrow(non_numeric_bp), " breakpoint(s) are not numeric and would become NA, see above")
}

breakpoints_eucast <- breakpoints_eucast |>
  mutate(
    # the check above guarantees purely numeric values, so this conversion
    # cannot introduce NAs (blocking rows are NA by design)
    breakpoint_S = as.numeric(breakpoint_S),
    breakpoint_R = as.numeric(breakpoint_R),
    # EUCAST's "arbitrary, off scale" MIC breakpoint S <= 0.001 mg/L (from
    # v10.0 on) is meant to categorise wild-type isolates as "I", not "S".
    # MICs of 0.0005 mg/L and below do occur in vitro, so, as before, it is
    # stored as 0.0001 mg/L to keep such isolates out of "S".
    breakpoint_S = if_else(method == "MIC" & !is.na(breakpoint_S) & breakpoint_S == 0.001, 0.0001, breakpoint_S),
    # same definition as used for the WHONET-based rows
    uti = !is.na(site) & site %like% "(UTI|urinary|urine)",
    is_SDD = FALSE # EUCAST has no "susceptible, dose dependent" category (CLSI-only concept)
  ) |>
  # EUCAST prints the twofold dilutions below 0.125 mg/L rounded, e.g. 0.03 and
  # 0.06 for 0.03125 and 0.0625, while this package (as.mic(), and the WHONET-
  # based breakpoints of all other guidelines) uses 0.032 and 0.064. Since
  # as.sir() compares values literally, an MIC of 0.064 would otherwise fall
  # above an S breakpoint of 0.06. MIC breakpoints are therefore set to the
  # package's dilution level they denote; a value that is not close to any
  # dilution level stops the script. The raw data keep EUCAST's notation.
  mutate(across(c(breakpoint_S, breakpoint_R), function(x) {
    if (!all(method %in% c("MIC", "DISK"))) stop("unexpected method")
    is_mic <- method == "MIC" & !is.na(x)
    levels <- AMR:::COMMON_MIC_VALUES
    nearest <- vapply(x[is_mic], function(z) levels[which.min(abs(log2(z) - log2(levels)))], numeric(1))
    off_level <- abs(log2(x[is_mic]) - log2(nearest)) > 0.15
    if (any(off_level)) {
      stop("MIC breakpoint(s) not on a dilution level: ", paste(unique(x[is_mic][off_level]), collapse = ", "))
    }
    x[is_mic] <- nearest
    x
  })) |>
  # Greek symbols and EM dash symbols are not allowed by CRAN, so replace them with ASCII:
  mutate(disk_dose = disk_dose %>%
           gsub("\u03bc", "mc", ., fixed = TRUE) %>% # this is 'mu', \u03bc
           gsub("\u00b5", "mc", ., fixed = TRUE) %>% # this is 'micro', \u00b5 (yes, they look the same)
           gsub("\u2013", "-", ., fixed = TRUE) %>%
           gsub("(?<=\\d)(?=[a-zA-Z])", " ", ., perl = TRUE)) # keep a space after a number, e.g. "1mcg" to "1 mcg"

# Every row must have resolved to an organism and an antimicrobial. as.mo()
# returns the code "UNKNOWN" rather than NA when it finds no match at all
# (unlike as.ab(), which returns NA), so both are checked. Every EUCAST table
# names its organisms and agents, so an unresolved one always means a gap in
# the parsing above (e.g. a new free-text pattern), never a genuine ambiguity.
# The only exception is the non-species related PK-PD table, which is coded
# as "UNKNOWN" on purpose (as.sir() uses it as the last resort).
unresolved_mo <- breakpoints_eucast |> filter(is.na(mo) | (mo == "UNKNOWN" & !grepl("^PK", sheet))) |> distinct(sheet, ref_tbl, ab, mo_text)
unresolved_ab <- breakpoints_eucast |> filter(is.na(ab)) |> distinct(sheet, ref_tbl)
if (nrow(unresolved_mo) > 0 || nrow(unresolved_ab) > 0) {
  print(unresolved_mo, n = Inf)
  print(unresolved_ab, n = Inf)
  stop(nrow(unresolved_mo), " organism(s) and ", nrow(unresolved_ab), " antimicrobial(s) could not be resolved, see above")
}

# Pruning blocking rows ----
# A blocking row is only needed where its organism has no (non-screening)
# breakpoint of its own for that agent and method, and blocking rows for the
# PK-PD table ("UNKNOWN") have no function. Of several blocking rows for the
# same organism, agent and method, the first is kept.
has_own_bp <- breakpoints_eucast |>
  filter(!is_blocker, !grepl("screen", site, ignore.case = TRUE)) |>
  distinct(guideline, type, host, mo, ab, method)
breakpoints_eucast <- bind_rows(
  breakpoints_eucast |> filter(!is_blocker),
  breakpoints_eucast |>
    filter(is_blocker, mo != "UNKNOWN") |>
    anti_join(has_own_bp, by = c("guideline", "type", "host", "mo", "ab", "method")) |>
    distinct(guideline, type, host, mo, ab, method, .keep_all = TRUE) |>
    mutate(site = NA_character_, uti = FALSE, disk_dose = NA_character_)
)

breakpoints_eucast <- breakpoints_eucast |>
  distinct(guideline, type, host, ab, mo, method, site, breakpoint_S, .keep_all = TRUE) |>
  mutate(
    rank_index = case_when(
      mo == "UNKNOWN" ~ 7, # the PK-PD table, as for the WHONET-based rows
      is.na(mo_rank(mo, keep_synonyms = TRUE)) ~ 6, # for B_GRAMN, B_ANAER, B_ANAER-NEG, etc.
      mo_rank(mo, keep_synonyms = TRUE) %like% "(infra|sub)" ~ 1,
      mo_rank(mo, keep_synonyms = TRUE) == "species" ~ 2,
      mo_rank(mo, keep_synonyms = TRUE) == "species group" ~ 2.5,
      mo_rank(mo, keep_synonyms = TRUE) == "genus" ~ 3,
      mo_rank(mo, keep_synonyms = TRUE) == "family" ~ 4,
      mo_rank(mo, keep_synonyms = TRUE) == "order" ~ 5,
      TRUE ~ 6
    ),
    # a breakpoint of the same organism always precedes its blocking row
    rank_index = rank_index + if_else(is_blocker, 0.1, 0)
  )

# Integrity checks on the final table ----
# 1. One breakpoint per guideline, organism, agent, method and site: two rows
#    sharing these but differing in their values would leave as.sir() to pick
#    one of them arbitrarily.
conflicting_bp <- breakpoints_eucast |>
  group_by(guideline, type, host, ab, mo, method, site) |>
  filter(n() > 1) |>
  ungroup()
if (nrow(conflicting_bp) > 0) {
  print(conflicting_bp |> select(guideline, mo, ab, method, site, disk_dose, breakpoint_S, breakpoint_R, ref_tbl), n = Inf)
  stop(nrow(conflicting_bp), " rows have more than one breakpoint for the same guideline, organism, agent, method and site, see above")
}
# 2. Breakpoints must be in the right order: for MICs S <= R, for zone
#    diameters S >= R.
misordered_bp <- breakpoints_eucast |>
  filter((method == "MIC" & breakpoint_S > breakpoint_R) | (method == "DISK" & breakpoint_S < breakpoint_R))
if (nrow(misordered_bp) > 0) {
  print(misordered_bp |> select(guideline, mo, ab, method, site, breakpoint_S, breakpoint_R, ref_tbl), n = Inf)
  stop(nrow(misordered_bp), " rows have S and R breakpoints in the wrong order, see above")
}
# 3. A disk breakpoint needs a disk content, an MIC breakpoint none.
missing_dose <- breakpoints_eucast |>
  filter(!is_blocker, (method == "DISK" & is.na(disk_dose)) | (method == "MIC" & !is.na(disk_dose)))
if (nrow(missing_dose) > 0) {
  print(missing_dose |> select(guideline, mo, ab, method, disk_dose, ref_tbl), n = Inf)
  stop(nrow(missing_dose), " rows have a disk content that does not match their method, see above")
}
# 4. Blocking rows have no values and state their reason; other rows have an
#    S breakpoint, and only screening tests may lack an R breakpoint.
bad_rows <- breakpoints_eucast |>
  filter((is_blocker & (!is.na(breakpoint_S) | !is.na(breakpoint_R) | !grepl("^\\[No breakpoint\\] ", note))) |
           (!is_blocker & (is.na(breakpoint_S) | (is.na(breakpoint_R) & !grepl("screen", site, ignore.case = TRUE)))))
if (nrow(bad_rows) > 0) {
  print(bad_rows |> select(guideline, mo, ab, method, site, breakpoint_S, breakpoint_R, ref_tbl, note), n = Inf)
  stop(nrow(bad_rows), " rows are neither a complete breakpoint nor a proper blocking row, see above")
}

breakpoints_eucast <- breakpoints_eucast |>
  select(guideline, type, host, method, site, mo, rank_index, ab, ref_tbl,
         disk_dose, breakpoint_S, breakpoint_R, uti, is_SDD, note)

breakpoints_eucast |> count(guideline) |> print()

# Which guideline years these workbooks cover, separately for bacteria and
# fungi: the bacterial tables and the antifungal (AFST) tables are separate
# EUCAST documents, and only these years may replace the WHONET-based rows
# (see merge_eucast_breakpoints()).
eucast_coverage <- breakpoints_eucast_raw |>
  distinct(guideline, file_desc) |>
  mutate(scope = if_else(file_desc %like% "Antifungal", "Fungi", "Bacteria")) |>
  distinct(guideline, scope)

# Merges these EUCAST breakpoints into the (WHONET-based) clinical
# breakpoints: human EUCAST rows are replaced for exactly those guideline
# years and organism groups (bacteria or fungi) that the EUCAST workbooks
# here cover. Everything else is kept as is, e.g. EUCAST antifungal
# breakpoints of guideline years for which no AFST workbook is available,
# EUCAST ECOFFs, EUCAST animal breakpoints and all CLSI breakpoints.
merge_eucast_breakpoints <- function(clinical_breakpoints, breakpoints_eucast, eucast_coverage) {
  # bacteria or fungi; for the non-species related rows ("UNKNOWN"), this
  # follows the agent, e.g. the PK-PD breakpoints of fluconazole belong to
  # the antifungal tables
  scope_of <- function(mo, ab) {
    if_else(mo_kingdom(mo, keep_synonyms = TRUE) %in% "Fungi" |
              (as.character(mo) == "UNKNOWN" & ab_group(ab) %in% "Antifungals"),
            "Fungi", "Bacteria")
  }
  # the organisms in the EUCAST rows must agree with the type of workbook
  # they come from
  new_scope <- breakpoints_eucast |>
    mutate(scope = scope_of(mo, ab), from_afst = ref_tbl %like% "Antifungal")
  mismatch <- new_scope |> filter((scope == "Fungi") != from_afst)
  if (nrow(mismatch) > 0) {
    print(mismatch |> select(guideline, mo, ab, method, ref_tbl), n = Inf)
    stop("merge_eucast_breakpoints(): organisms do not match the type of EUCAST workbook, see above")
  }
  replaced <- clinical_breakpoints |>
    mutate(scope = scope_of(mo, ab)) |>
    semi_join(eucast_coverage, by = c("guideline", "scope")) |>
    filter(type == "human")
  kept <- clinical_breakpoints |>
    anti_join(replaced, by = names(clinical_breakpoints))
  message("Replacing ", nrow(replaced), " WHONET-based human EUCAST rows (",
          paste(sort(unique(paste(replaced$guideline, replaced$scope))), collapse = ", "),
          ") with ", nrow(breakpoints_eucast), " rows from the EUCAST workbooks")
  # rows without any breakpoint are only kept as the blocking rows of the
  # EUCAST workbooks (see "Blocking rows")
  kept |>
    filter(!(is.na(breakpoint_S) & is.na(breakpoint_R)) & !is.na(mo) & !is.na(ab)) |>
    bind_rows(breakpoints_eucast) |>
    arrange(desc(guideline), mo, ab, type, host, method, rank_index) |>
    distinct(guideline, type, host, ab, mo, method, site, breakpoint_S, .keep_all = TRUE) |>
    dataset_UTF8_to_ASCII()
}
