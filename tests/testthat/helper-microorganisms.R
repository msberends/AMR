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

# Integrity rules of the `microorganisms` data set.
# This file is used by the unit tests (test-data-microorganisms.R) AND by the taxonomy build
# (data-raw/_reproduction_scripts/reproduction_of_microorganisms.R), which stops before saving if any rule
# is broken. So the rules exist in one place only.

# the MO code prefix of each domain
MO_DOMAIN_PREFIX <- c(
  "Animalia" = "AN",
  "Archaea" = "A",
  "Bacteria" = "B",
  "Chromista" = "C",
  "Fungi" = "F",
  "Plantae" = "PL",
  "Protozoa" = "P"
)

# read the MO code registry and the approved renames (see data-raw/microorganisms_files/README.md),
# returns NULL if the files are not available (such as in an installed package)
read_mo_registry <- function(root = mo_repository_root()) {
  file <- file.path(root, "data-raw", "microorganisms_files", "mo_code_registry.csv")
  if (is.null(root) || !file.exists(file)) {
    return(NULL)
  }
  utils::read.csv(file, colClasses = "character", na.strings = character(0))
}
read_mo_renames <- function(root = mo_repository_root()) {
  file <- file.path(root, "data-raw", "microorganisms_files", "mo_code_renames.csv")
  if (is.null(root) || !file.exists(file)) {
    return(NULL)
  }
  renames <- utils::read.csv(file, colClasses = "character", na.strings = character(0))
  # only approved renames count
  renames[renames$approved_by != "" & renames$approved_date != "", , drop = FALSE]
}
# MO codes that are retired on purpose, because their taxon was another organism (see the README)
read_mo_retirements <- function(root = mo_repository_root()) {
  file <- file.path(root, "data-raw", "microorganisms_files", "mo_code_retirements.csv")
  if (is.null(root) || !file.exists(file)) {
    return(NULL)
  }
  utils::read.csv(file, colClasses = "character", na.strings = character(0))
}
# known defects of the data set on the main branch (see test-data-microorganisms.R), returns an empty data set
# if the file is not there, as it must be deleted once no defects are left
read_mo_known_defects <- function(root = mo_repository_root()) {
  file <- file.path(root, "data-raw", "microorganisms_files", "microorganisms_known_defects.csv")
  if (is.null(root) || !file.exists(file)) {
    return(data.frame(rule = character(0), record = character(0)))
  }
  utils::read.csv(file, colClasses = "character", na.strings = character(0))
}
# the root of the source repository, found by walking up from the working directory; this works from the root
# itself, from tests/testthat, and from an R CMD check directory inside the repository (as on GitHub Actions).
# Returns NULL elsewhere, such as for an installed package, since data-raw/ never ships with the package.
mo_repository_root <- function() {
  path <- normalizePath(".", mustWork = FALSE)
  for (i in seq_len(6)) {
    if (file.exists(file.path(path, "data-raw", "microorganisms_files", "mo_code_registry.csv")) &&
      file.exists(file.path(path, "DESCRIPTION")) &&
      identical(unname(read.dcf(file.path(path, "DESCRIPTION"), fields = "Package")[1, 1]), "AMR")) {
      return(path)
    }
    parent <- dirname(path)
    if (parent == path) {
      break
    }
    path <- parent
  }
  NULL
}
skip_if_no_mo_repository <- function() {
  testthat::skip_if(is.null(mo_repository_root()), "Source repository not available (data-raw/ does not ship with the package)")
}

# names without a rank suffix, e.g. "Kapabacteria {class}" and "Nitrospira (class)" are both the class
mo_name_without_suffix <- function(x) {
  trimws(gsub(" [{(][a-z ]+[)}]$", "", x))
}

# Returns a named list with one character vector per rule, each describing the records that break the rule.
# All vectors must be empty.
mo_integrity_issues <- function(df, registry = NULL, renames = NULL, retirements = NULL) {
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  df$mo <- as.character(df$mo)
  describe <- function(rows) {
    if (!any(rows, na.rm = TRUE)) {
      return(character(0))
    }
    rows <- which(rows)
    paste0(df$mo[rows], " (", df$fullname[rows], ")")
  }
  issues <- list()

  # identifiers
  issues$duplicated_mo <- describe(duplicated(df$mo) | duplicated(df$mo, fromLast = TRUE))
  issues$duplicated_fullname <- describe(duplicated(df$fullname) | duplicated(df$fullname, fromLast = TRUE))
  issues$empty_mo <- describe(df$mo %in% c("", NA))
  issues$invalid_mo_characters <- describe(df$mo %unlike_case% "^[][A-Z0-9_-]+$")
  issues$repeated_mo_elements <- describe(df$mo %like_case% "^([A-Z]+_[^_]+)_\\1_")

  # the code prefix must match the domain
  prefix <- sub("_.*", "", df$mo)
  expected_prefix <- unname(MO_DOMAIN_PREFIX[df$domain])
  issues$prefix_not_matching_domain <- describe(
    df$mo != "UNKNOWN" & !is.na(expected_prefix) & prefix != expected_prefix
  )
  issues$unknown_domain <- describe(!df$domain %in% c(names(MO_DOMAIN_PREFIX), "(unknown domain)"))

  # the code of a species must start with the code of its genus, of a subspecies with the code of its species
  is_taxon <- df$fullname %unlike% "unknown"
  genus_mo <- df$mo[df$rank == "genus"][match(
    paste(df$domain, df$genus),
    paste(df$domain[df$rank == "genus"], df$genus[df$rank == "genus"])
  )]
  species_mo <- df$mo[df$rank == "species"][match(
    paste(df$domain, df$genus, df$species),
    paste(df$domain[df$rank == "species"], df$genus[df$rank == "species"], df$species[df$rank == "species"])
  )]
  issues$species_without_genus <- describe(df$rank == "species" & is_taxon & is.na(genus_mo))
  issues$species_code_not_under_genus <- describe(
    df$rank == "species" & is_taxon & !is.na(genus_mo) & !startsWith(df$mo, paste0(genus_mo, "_"))
  )
  # (Salmonella serovars are coded as genus + serovar by design, e.g. B_SLMNL_ACHN for Salmonella Aachen)
  issues$subspecies_code_not_under_species <- describe(
    df$rank == "subspecies" & is_taxon & df$genus != "Salmonella" &
      !is.na(species_mo) & !startsWith(df$mo, paste0(species_mo, "_"))
  )
  issues$subspecies_code_not_under_genus <- describe(
    df$rank == "subspecies" & is_taxon & !is.na(genus_mo) & !startsWith(df$mo, paste0(genus_mo, "_"))
  )

  # sources can only contain their own scope (issue #309)
  issues$lpsn_outside_prokaryotes <- describe(df$source == "LPSN" & !df$domain %in% c("Bacteria", "Archaea"))
  issues$mycobank_outside_fungi <- describe(df$source == "MycoBank" & df$domain != "Fungi")

  # allowed values
  issues$invalid_status <- describe(!df$status %in% c("accepted", "synonym", "unknown"))
  issues$invalid_rank <- describe(!df$rank %in% c(
    "domain", "kingdom", "phylum", "class", "order", "family", "genus", "species", "subspecies",
    "species group", "(unknown rank)"
  ))
  issues$invalid_source <- describe(!df$source %in% c("LPSN", "MycoBank", "GBIF", "manually added"))
  issues$invalid_prevalence <- describe(is.na(df$prevalence) | df$prevalence < 1 | df$prevalence > 2)

  # synonyms with a current name must lead to an existing record, without cycles or endless chains; a chain
  # that ends at a synonym without a current name is incomplete data, see mo_synonyms_without_current_name()
  if (all(c("lpsn_renamed_to", "gbif_renamed_to") %in% colnames(df))) {
    chain <- mo_synonym_chain_end(df)
    issues$synonym_with_unknown_target <- describe(chain$has_target & is.na(chain$end))
    issues$synonym_chain_not_ending <- describe(
      chain$has_target & !is.na(chain$end) & chain$end_is_synonym & chain$end_has_target
    )
  }

  # the MO code registry: a registered code must always denote the same taxon
  if (!is.null(registry)) {
    registered_name <- registry$fullname[match(df$mo, registry$mo)]
    if (!is.null(renames) && nrow(renames) > 0) {
      renamed <- match(df$mo, renames$mo)
      registered_name[!is.na(renamed)] <- renames$new_name[renamed[!is.na(renamed)]]
    }
    issues$registered_code_with_other_taxon <- describe(
      !is.na(registered_name) & mo_name_without_suffix(registered_name) != mo_name_without_suffix(df$fullname)
    )
    # a released taxon is never removed, unless its code was retired on purpose
    expected <- registry$fullname
    if (!is.null(renames) && nrow(renames) > 0) {
      renamed <- match(registry$mo, renames$mo)
      expected[!is.na(renamed)] <- renames$new_name[renamed[!is.na(renamed)]]
    }
    missing <- !mo_name_without_suffix(expected) %in% mo_name_without_suffix(df$fullname)
    if (!is.null(retirements)) {
      missing <- missing & !registry$mo %in% retirements$mo
    }
    issues$released_taxon_missing <- if (any(missing)) paste0(registry$mo[missing], " (", expected[missing], ")") else character(0)
  }

  issues
}

# the first available 'renamed to' identifier of each record
mo_renamed_to <- function(df) {
  out <- df$lpsn_renamed_to
  if ("mycobank_renamed_to" %in% colnames(df)) {
    out <- ifelse(is.na(out), df$mycobank_renamed_to, out)
  }
  ifelse(is.na(out), df$gbif_renamed_to, out)
}

# for every synonym with a target: where its chain of renames ends
mo_synonym_chain_end <- function(df) {
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  df$mo <- as.character(df$mo)
  has_target <- df$status == "synonym" & !is.na(mo_renamed_to(df))
  end <- rep(NA_character_, nrow(df))
  end[has_target] <- AMR:::synonym_mo_to_accepted_mo(df$mo[has_target], fill_in_accepted = FALSE, dataset = df)
  end_row <- match(end, df$mo)
  # a target only counts if it exists in the data set
  known_target <- (!is.na(df$lpsn_renamed_to) & df$lpsn_renamed_to %in% df$lpsn) |
    (!is.na(df$gbif_renamed_to) & df$gbif_renamed_to %in% df$gbif)
  if ("mycobank_renamed_to" %in% colnames(df)) {
    known_target <- known_target | (!is.na(df$mycobank_renamed_to) & df$mycobank_renamed_to %in% df$mycobank)
  }
  data.frame(
    has_target = has_target,
    end = end,
    end_is_synonym = df$status[end_row] %in% "synonym",
    end_has_target = known_target[end_row] %in% TRUE
  )
}

# synonyms for which no current name is known (e.g. because the record of the current name was not included),
# directly or at the end of their chain of renames (e.g. Candida glabratus -> Nakaseomyces glabrata, which has no
# current name itself); these are no error, but should be as few as possible
mo_synonyms_without_current_name <- function(df) {
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  chain <- mo_synonym_chain_end(df)
  no_target <- df$status == "synonym" & !chain$has_target
  dead_end_chain <- chain$has_target & !is.na(chain$end) & chain$end_is_synonym & !chain$end_has_target
  out <- df[no_target | dead_end_chain, c("mo", "fullname", "source"), drop = FALSE]
  out$chain_ends_at <- df$fullname[match(chain$end[no_target | dead_end_chain], df$mo)]
  out
}

# a readable summary of all broken rules, or NULL if there are none
mo_integrity_report <- function(issues) {
  issues <- issues[lengths(issues) > 0]
  if (length(issues) == 0) {
    return(NULL)
  }
  paste0(
    names(issues), " (", lengths(issues), "): ",
    vapply(issues, function(x) paste(utils::head(x, 10), collapse = ", "), character(1)),
    ifelse(lengths(issues) > 10, ", ...", ""),
    collapse = "\n"
  )
}
