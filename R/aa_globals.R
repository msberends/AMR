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
# how to conduct AMR data analysis: https://msberends.github.io/AMR/   #
# ==================================================================== #

# add new version numbers here, and add the rules themselves to "data-raw/eucast_rules.tsv" and clinical_breakpoints
# (sourcing "data-raw/_pre_commit_checks.R" will process the TSV file)
EUCAST_VERSION_BREAKPOINTS <- list(
  "14.0" = list(
    version_txt = "v14.0",
    year = 2024,
    title = "'EUCAST Clinical Breakpoint Tables'",
    url = "https://www.eucast.org/clinical_breakpoints/"
  ),
  "13.1" = list(
    version_txt = "v13.1",
    year = 2023,
    title = "'EUCAST Clinical Breakpoint Tables'",
    url = "https://www.eucast.org/clinical_breakpoints/"
  ),
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
EUCAST_VERSION_EXPECTED_PHENOTYPES <- list(
  "1.2" = list(
    version_txt = "v1.2",
    year = 2023,
    title = "'EUCAST Expected Resistant Phenotypes'",
    url = "https://www.eucast.org/expert_rules_and_expected_phenotypes"
  )
)

TAXONOMY_VERSION <- list(
  GBIF = list(
    name = "Global Biodiversity Information Facility (GBIF)",
    accessed_date = as.Date("2024-06-24"),
    citation = "GBIF Secretariat (2023). GBIF Backbone Taxonomy. Checklist dataset \\doi{10.15468/39omei}.",
    url = "https://www.gbif.org"
  ),
  LPSN = list(
    name = "List of Prokaryotic names with Standing in Nomenclature (LPSN)",
    accessed_date = as.Date("2024-06-24"),
    citation = "Parte, AC *et al.* (2020). **List of Prokaryotic names with Standing in Nomenclature (LPSN) moves to the DSMZ.** International Journal of Systematic and Evolutionary Microbiology, 70, 5607-5612; \\doi{10.1099/ijsem.0.004332}.",
    url = "https://lpsn.dsmz.de"
  ),
  MycoBank = list(
    name = "MycoBank",
    accessed_date = as.Date("2024-06-24"),
    citation = "Vincent, R *et al* (2013). **MycoBank gearing up for new horizons.** IMA Fungus, 4(2), 371-9; \\doi{10.5598/imafungus.2013.04.02.16}.",
    url = "https://www.mycobank.org"
  ),
  BacDive = list(
    name = "BacDive",
    accessed_date = as.Date("2024-07-16"),
    citation = "Reimer, LC *et al.* (2022). ***BacDive* in 2022: the knowledge base for standardized bacterial and archaeal data.** Nucleic Acids Res., 50(D1):D741-D74; \\doi{10.1093/nar/gkab961}.",
    url = "https://bacdive.dsmz.de"
  ),
  SNOMED = list(
    name = "Systematized Nomenclature of Medicine - Clinical Terms (SNOMED-CT)",
    accessed_date = as.Date("2024-07-16"),
    citation = "Public Health Information Network Vocabulary Access and Distribution System (PHIN VADS). US Edition of SNOMED CT from 1 September 2020. Value Set Name 'Microorganism', OID 2.16.840.1.114222.4.11.1009 (v12).",
    url = "https://www.cdc.gov/phin/php/phinvads/"
  ),
  LOINC = list(
    name = "Logical Observation Identifiers Names and Codes (LOINC)",
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
  "antimicrobials",
  "atc_group1",
  "atc_group2",
  "base_ab",
  "beta_posterior_1",
  "beta_posterior_2",
  "ci_max",
  "ci_min",
  "clinical_breakpoints",
  "code",
  "cols",
  "count",
  "coverage",
  "data",
  "disk",
  "dosage",
  "dose",
  "dose_times",
  "fullname",
  "fullname_lower",
  "g_species",
  "gamma_posterior",
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
  "lower",
  "lower_ci",
  "method",
  "mic ",
  "mic",
  "microorganism",
  "microorganisms",
  "microorganisms.codes",
  "mo",
  "n",
  "n_susceptible",
  "n_tested",
  "n_total",
  "name",
  "new",
  "numerator",
  "observations",
  "old",
  "old_name",
  "p_susceptible",
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
  "total_rows",
  "txt",
  "type",
  "upper",
  "upper_ci",
  "uti_index",
  "value",
  "varname",
  "x",
  "xvar",
  "y",
  "year",
  "yvar"
))
