# ==================================================================== #
# TITLE:                                                               #
# AMR: An R Package for Working with Antimicrobial Resistance Data     #
#                                                                      #
# SOURCE CODE:                                                         #
# https://github.com/msberends/AMR                                     #
#                                                                      #
# PLEASE CITE THIS SOFTWARE AS:                                        #
# Berends MS, Luz CF, Friedrich AW, Sinha BNM, Albers CJ, Glasner C    #
# (2022). AMR: An R Package for Working with Antimicrobial Resistance  #
# Data. Journal of Statistical Software, 104(3), 1-31.                 #
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
# how to conduct AMR data analysis: https://msberends.github.io/AMR/   #
# ==================================================================== #

# add new version numbers here, and add the rules themselves to "data-raw/eucast_rules.tsv" and clinical_breakpoints
# (sourcing "data-raw/_pre_commit_checks.R" will process the TSV file)
EUCAST_VERSION_BREAKPOINTS <- list(
  # "13.0" = list(
  #   version_txt = "v13.0",
  #   year = 2023,
  #   title = "'EUCAST Clinical Breakpoint Tables'",
  #   url = "https://www.eucast.org/clinical_breakpoints/"
  # ),
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
  "3.3" = list(
    version_txt = "v3.3",
    year = 2021,
    title = "'EUCAST Expert Rules' and 'EUCAST Intrinsic Resistance and Unusual Phenotypes'",
    url = "https://www.eucast.org/expert_rules_and_expected_phenotypes"
  ),
  "3.2" = list(
    version_txt = "v3.2",
    year = 2020,
    title = "'EUCAST Expert Rules' and 'EUCAST Intrinsic Resistance and Unusual Phenotypes'",
    url = "https://www.eucast.org/expert_rules_and_expected_phenotypes"
  ),
  "3.1" = list(
    version_txt = "v3.1",
    year = 2016,
    title = "'EUCAST Expert Rules, Intrinsic Resistance and Exceptional Phenotypes'",
    url = "https://www.eucast.org/expert_rules_and_expected_phenotypes"
  )
)
# EUCAST_VERSION_RESISTANTPHENOTYPES <- list(
#   "1.2" = list(
#     version_txt = "v1.2",
#     year = 2023,
#     title = "'Expected Resistant Phenotypes'",
#     url = "https://www.eucast.org/expert_rules_and_expected_phenotypes"
#   )
# )

TAXONOMY_VERSION <- list(
  GBIF = list(
    accessed_date = as.Date("2024-01-08"),
    citation = "GBIF Secretariat (2023). GBIF Backbone Taxonomy. Checklist dataset \\doi{10.15468/39omei}.",
    url = "https://www.gbif.org"
  ),
  LPSN = list(
    accessed_date = as.Date("2022-12-11"),
    citation = "Parte, AC *et al.* (2020). **List of Prokaryotic names with Standing in Nomenclature (LPSN) moves to the DSMZ.** International Journal of Systematic and Evolutionary Microbiology, 70, 5607-5612; \\doi{10.1099/ijsem.0.004332}.",
    url = "https://lpsn.dsmz.de"
  ),
  BacDive = list(
    accessed_date = as.Date("2023-05-12"),
    citation = "Reimer, LC *et al.* (2022). ***BacDive* in 2022: the knowledge base for standardized bacterial and archaeal data.** Nucleic Acids Res., 50(D1):D741-D74; \\doi{10.1093/nar/gkab961}.",
    url = "https://bacdive.dsmz.de"
  ),
  SNOMED = list(
    accessed_date = as.Date("2021-07-01"),
    citation = "Public Health Information Network Vocabulary Access and Distribution System (PHIN VADS). US Edition of SNOMED CT from 1 September 2020. Value Set Name 'Microorganism', OID 2.16.840.1.114222.4.11.1009 (v12).",
    url = "https://phinvads.cdc.gov"
  ),
  LOINC = list(
    accessed_date = as.Date("2023-10-19"),
    citation = "Logical Observation Identifiers Names and Codes (LOINC), Version 2.76 (18 September, 2023).",
    url = "https://loinc.org"
  )
)

globalVariables(c(
  ".GenericCallEnv",
  ".mo",
  ".rowid",
  ".syndromic_group",
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
  "ci_max",
  "ci_min",
  "clinical_breakpoints",
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
  "host_index",
  "host_match",
  "input",
  "intrinsic_resistant",
  "isolates",
  "lang",
  "language",
  "lookup",
  "method",
  "mic ",
  "mic",
  "microorganism",
  "microorganisms",
  "microorganisms.codes",
  "mo",
  "n",
  "name",
  "new",
  "numerator",
  "observations",
  "old",
  "old_name",
  "pattern",
  "R",
  "rank_index",
  "ref_tbl",
  "reference.rule",
  "reference.rule_group",
  "reference.version",
  "rowid",
  "rule_group",
  "rule_name",
  "se_max",
  "se_min",
  "SI",
  "sir",
  "species",
  "syndromic_group",
  "total",
  "txt",
  "type",
  "uti_index",
  "value",
  "varname",
  "x",
  "xvar",
  "y",
  "year",
  "yvar"
))
