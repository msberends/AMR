# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Data Analysis for R                   #
#                                                                      #
# SOURCE                                                               #
# https://github.com/msberends/AMR                                     #
#                                                                      #
# LICENCE                                                              #
# (c) 2018-2021 Berends MS, Luz CF et al.                              #
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
# (sourcing "data-raw/_internals.R" will process the TSV file)
EUCAST_VERSION_BREAKPOINTS <- list("11.0" = list(version_txt = "v11.0",
                                                 year = 2021, 
                                                 title = "'EUCAST Clinical Breakpoint Tables'",
                                                 url = "https://www.eucast.org/clinical_breakpoints/"),
                                   "10.0" = list(version_txt = "v10.0",
                                                 year = 2020, 
                                                 title = "'EUCAST Clinical Breakpoint Tables'",
                                                 url = "https://www.eucast.org/ast_of_bacteria/previous_versions_of_documents/"))
EUCAST_VERSION_EXPERT_RULES <- list("3.1" = list(version_txt = "v3.1",
                                                 year = 2016, 
                                                 title = "'EUCAST Expert Rules, Intrinsic Resistance and Exceptional Phenotypes'",
                                                 url = "https://www.eucast.org/expert_rules_and_intrinsic_resistance/"),
                                    "3.2" = list(version_txt = "v3.2",
                                                 year = 2020, 
                                                 title = "'EUCAST Expert Rules' and 'EUCAST Intrinsic Resistance and Unusual Phenotypes'",
                                                 url = "https://www.eucast.org/expert_rules_and_intrinsic_resistance/"))

SNOMED_VERSION <- list(title = "Public Health Information Network Vocabulary Access and Distribution System (PHIN VADS)",
                       current_source = "US Edition of SNOMED CT from 1 September 2020",
                       current_version = 12,
                       current_oid = "2.16.840.1.114222.4.11.1009",
                       value_set_name = "Microorganism",
                       url = "https://phinvads.cdc.gov/vads/ViewValueSet.action?oid=2.16.840.1.114222.4.11.1009")

CATALOGUE_OF_LIFE <- list(
  year = 2019,
  version = "Catalogue of Life: {year} Annual Checklist",
  url_CoL = "http://www.catalogueoflife.org/col/",
  url_LPSN = "https://lpsn.dsmz.de",
  yearmonth_LPSN = "5 October 2021"
)

globalVariables(c(".rowid",
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
                  "microorganisms.old",
                  "mo",
                  "name",
                  "new",
                  "observations",
                  "old",
                  "old_name",
                  "pattern",
                  "R",
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
                  "species_id",
                  "total",
                  "txt",
                  "type",
                  "value",
                  "varname",
                  "xvar",
                  "y",
                  "year",
                  "yvar"))
