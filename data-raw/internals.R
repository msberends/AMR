# ---------------------------------------------------------------------------------------------------
# For editing this EUCAST reference file, these values can all be used for target antibiotics:
# all_betalactams, aminoglycosides, carbapenems, cephalosporins, cephalosporins_without_CAZ, fluoroquinolones, 
# glycopeptides, macrolides, minopenicillins, polymyxins, streptogramins, tetracyclines, ureidopenicillins
# and all separate EARS-Net letter codes like AMC. They can be separated by comma: 'AMC, fluoroquinolones'.
# The if_mo_property column can be any column name from the AMR::microorganisms data set, or "genus_species" or "gramstain".
# The EUCAST guideline contains references to the 'Burkholderia cepacia complex'. All species in this group can be found in:
# LiPuma J, Curr Opin Pulm Med. 2005 Nov;11(6):528-33. (PMID 16217180).
# >>>>> IF YOU WANT TO IMPORT THIS FILE INTO YOUR OWN SOFTWARE, HAVE THE FIRST 10 LINES SKIPPED <<<<<
# ---------------------------------------------------------------------------------------------------
eucast_rules_file <- dplyr::arrange(
  .data = utils::read.delim(file = "data-raw/eucast_rules.tsv",
                            skip = 10,
                            sep = "\t",
                            stringsAsFactors = FALSE,
                            header = TRUE,
                            strip.white = TRUE,
                            na = c(NA, "", NULL)),
  reference.rule_group,
  reference.rule)

# Translations -----
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

# Export to package as internal data ----
usethis::use_data(eucast_rules_file, translations_file,
                  internal = TRUE,
                  overwrite = TRUE,
                  version = 2)

# Remove from global environment ----
rm(eucast_rules_file)
rm(translations_file)
