# EUCAST rules ----
# For editing the reference file, these values can all be used for target antibiotics:
# "aminoglycosides", "tetracyclines", "polymyxins", "macrolides", "glycopeptides",
# "streptogramins", "cephalosporins", "cephalosporins_without_CAZ", "carbapenems",
# "minopenicillins", "ureidopenicillins", "fluoroquinolones", "all_betalactams",
# and all separate EARS-Net letter codes like "AMC". They can be separated by comma: "AMC, fluoroquinolones".
# The mo_property can be any column name from the AMR::microorganisms data set, or "genus_species" or "gramstain".
# This file contains references to the 'Burkholderia cepacia complex'. The species in this group can be found in:
# LiPuma JJ, 2015 (PMID 16217180).
eucast_rules_file <- dplyr::arrange(
  .data = utils::read.delim(file = "data-raw/eucast_rules.tsv",
                            sep = "\t",
                            stringsAsFactors = FALSE,
                            header = TRUE,
                            strip.white = TRUE,
                            na = c(NA, "", NULL)),
  reference.rule_group,
  reference.rule)

# Translations -----
translations_file <- utils::read.table(file = "data-raw/translations.tsv",
                                       sep = "\t",
                                       stringsAsFactors = FALSE,
                                       header = TRUE,
                                       blank.lines.skip = TRUE,
                                       fill = TRUE,
                                       strip.white = TRUE,
                                       encoding = "UTF-8",
                                       fileEncoding = "UTF-8",
                                       na.strings = c(NA, "", NULL))

# Export to package as internal data ----
usethis::use_data(eucast_rules_file, translations_file,
                  internal = TRUE,
                  overwrite = TRUE,
                  version = 2)

# Remove from global environment ----
rm(eucast_rules_file)
rm(translations_file)
