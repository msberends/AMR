# ==================================================================== #
# TITLE                                                                #
# AMR: An R Package for Working with Antimicrobial Resistance Data     #
#                                                                      #
# SOURCE                                                               #
# https://github.com/msberends/AMR                                     #
#                                                                      #
# CITE AS                                                              #
# Berends MS, Luz CF, Friedrich AW, Sinha BNM, Albers CJ, Glasner C    #
# (2022). AMR: An R Package for Working with Antimicrobial Resistance  #
# Data. Journal of Statistical Software, 104(3), 1-31.                 #
# doi:10.18637/jss.v104.i03                                            #
#                                                                      #
# Developed at the University of Groningen, the Netherlands, in        #
# collaboration with non-profit organisations Certe Medical            #
# Diagnostics & Advice, and University Medical Center Groningen.       #
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
# how to conduct AMR data analysis: https://msberends.github.io/AMR/   #
# ==================================================================== #

# add new version numbers here, and add the rules themselves to "data-raw/eucast_rules.tsv" and rsi_translation
# (sourcing "data-raw/_pre_commit_hook.R" will process the TSV file)
EUCAST_VERSION_BREAKPOINTS <- list(
  "12.0" = list(
    version_txt = "v12.0",
    year = 2022,
    title = "'EUCAST Clinical Breakpoint Tables'",
    url = "https://www.eucast.org/clinical_breakpoints/"
  ),
  "11.0" = list(
    version_txt = "v11.0",
    year = 2021,
    title = "'EUCAST Clinical Breakpoint Tables'",
    url = "https://www.eucast.org/clinical_breakpoints/"
  ),
  "10.0" = list(
    version_txt = "v10.0",
    year = 2020,
    title = "'EUCAST Clinical Breakpoint Tables'",
    url = "https://www.eucast.org/ast_of_bacteria/previous_versions_of_documents/"
  )
)
EUCAST_VERSION_EXPERT_RULES <- list(
  "3.1" = list(
    version_txt = "v3.1",
    year = 2016,
    title = "'EUCAST Expert Rules, Intrinsic Resistance and Exceptional Phenotypes'",
    url = "https://www.eucast.org/expert_rules_and_expected_phenotypes/"
  ),
  "3.2" = list(
    version_txt = "v3.2",
    year = 2020,
    title = "'EUCAST Expert Rules' and 'EUCAST Intrinsic Resistance and Unusual Phenotypes'",
    url = "https://www.eucast.org/expert_rules_and_expected_phenotypes/"
  ),
  "3.3" = list(
    version_txt = "v3.3",
    year = 2021,
    title = "'EUCAST Expert Rules' and 'EUCAST Intrinsic Resistance and Unusual Phenotypes'",
    url = "https://www.eucast.org/expert_rules_and_expected_phenotypes/"
  )
)

TAXONOMY_VERSION <- list(
  GBIF = list(
    accessed_date = as.Date("2022-12-11"),
    citation = "GBIF Secretariat (2022). GBIF Backbone Taxonomy. Checklist dataset \\doi{10.15468/39omei}.",
    url = "https://www.gbif.org"
  ),
  LPSN = list(
    accessed_date = as.Date("2022-12-11"),
    citation = "Parte, AC *et al.* (2020). **List of Prokaryotic names with Standing in Nomenclature (LPSN) moves to the DSMZ.** International Journal of Systematic and Evolutionary Microbiology, 70, 5607-5612; \\doi{10.1099/ijsem.0.004332}.",
    url = "https://lpsn.dsmz.de"
  ),
  SNOMED = list(
    accessed_date = as.Date("2021-07-01"),
    citation = "Public Health Information Network Vocabulary Access and Distribution System (PHIN VADS). US Edition of SNOMED CT from 1 September 2020. Value Set Name 'Microoganism', OID 2.16.840.1.114222.4.11.1009 (v12).",
    url = "https://phinvads.cdc.gov"
  ),
  LOINC = list(
    accessed_date = as.Date("2022-10-30"),
    citation = "Logical Observation Identifiers Names and Codes (LOINC), Version 2.73 (8 August, 2022).",
    url = "https://loinc.org"
  )
)

globalVariables(c(
  ".rowid",
  "ab",
  "ab_txt",
  "affect_ab_name",
  "affect_mo_name",
  "angle",
  "antibiotic",
  "antibiotics",
  "atc_group1",
  "atc_group2",
  "base_ab",
  "ci_min",
  "ci_max",
  "code",
  "cols",
  "count",
  "data",
  "disk",
  "dosage",
  "dose",
  "dose_times",
  "fullname",
  "fullname_lower",
  "g_species",
  "genus",
  "gr",
  "group",
  "guideline",
  "hjust",
  "input",
  "intrinsic_resistant",
  "isolates",
  "lang",
  "language",
  "lookup",
  "method",
  "mic",
  "mic ",
  "microorganism",
  "microorganisms",
  "microorganisms.codes",
  "mo",
  "name",
  "new",
  "observations",
  "old",
  "old_name",
  "pattern",
  "R",
  "rank_index",
  "reference.rule",
  "reference.rule_group",
  "reference.version",
  "rowid",
  "rsi",
  "rsi_translation",
  "rule_group",
  "rule_name",
  "se_max",
  "se_min",
  "species",
  "total",
  "txt",
  "type",
  "value",
  "varname",
  "xvar",
  "y",
  "year",
  "yvar"
))
