# Run this file to update the package -------------------------------------
# source("data-raw/internals.R")

# See 'data-raw/eucast_rules.tsv' for the EUCAST reference file
eucast_rules_file <- utils::read.delim(file = "data-raw/eucast_rules.tsv",
                            skip = 10,
                            sep = "\t",
                            stringsAsFactors = FALSE,
                            header = TRUE,
                            strip.white = TRUE,
                            na = c(NA, "", NULL))
# take the order of the reference.rule_group column in the orginal data file
eucast_rules_file$reference.rule_group <- factor(eucast_rules_file$reference.rule_group,
                                                 levels = unique(eucast_rules_file$reference.rule_group),
                                                 ordered = TRUE)
eucast_rules_file <- dplyr::arrange(eucast_rules_file,
  reference.rule_group,
  reference.rule)


# Translations ----
translations_file <- utils::read.delim(file = "data-raw/translations.tsv",
                                       sep = "\t",
                                       stringsAsFactors = FALSE,
                                       header = TRUE,
                                       blank.lines.skip = TRUE,
                                       fill = TRUE,
                                       strip.white = TRUE,
                                       encoding = "UTF-8",
                                       fileEncoding = "UTF-8",
                                       na.strings = c(NA, "", NULL),
                                       allowEscapes = TRUE, # else "\\1" will be imported as "\\\\1"
                                       quote = "")

# Old microorganism codes -------------------------------------------------

microorganisms.translation <- readRDS("data-raw/microorganisms.translation.rds")

# Export to package as internal data ----
usethis::use_data(eucast_rules_file, translations_file, microorganisms.translation,
                  internal = TRUE,
                  overwrite = TRUE,
                  version = 2)

# Remove from global environment ----
rm(eucast_rules_file)
rm(translations_file)
rm(microorganisms.translation)

# Clean mo history ----
usethis::ui_done(paste0("Resetting {usethis::ui_value('mo_history.csv')}"))
tryCatch(
  write.csv(x = data.frame(x = character(0),
                           mo = character(0),
                           uncertainty_level = integer(0),
                           package_version = character(0),
                           stringsAsFactors = FALSE),
            row.names = FALSE,
            file = "inst/mo_history/mo_history.csv"),
  warning = function(w) cat("Warning:", w$message, "\n"),
  error = function(e) cat("Error:", e$message, "\n"))
