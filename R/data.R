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

#' Data Sets with `r format(nrow(antibiotics) + nrow(antivirals), big.mark = " ")` Antimicrobial Drugs
#'
#' Two data sets containing all antibiotics/antimycotics and antivirals. Use [as.ab()] or one of the [`ab_*`][ab_property()] functions to retrieve values from the [antibiotics] data set. Three identifiers are included in this data set: an antibiotic ID (`ab`, primarily used in this package) as defined by WHONET/EARS-Net, an ATC code (`atc`) as defined by the WHO, and a Compound ID (`cid`) as found in PubChem. Other properties in this data set are derived from one or more of these codes. Note that some drugs have multiple ATC codes.
#' @format
#' ### For the [antibiotics] data set: a [tibble][tibble::tibble] with `r nrow(antibiotics)` observations and `r ncol(antibiotics)` variables:
#' - `ab`\cr Antibiotic ID as used in this package (such as `AMC`), using the official EARS-Net (European Antimicrobial Resistance Surveillance Network) codes where available
#' - `cid`\cr Compound ID as found in PubChem
#' - `name`\cr Official name as used by WHONET/EARS-Net or the WHO
#' - `group`\cr A short and concise group name, based on WHONET and WHOCC definitions
#' - `atc`\cr ATC codes (Anatomical Therapeutic Chemical) as defined by the WHOCC, like `J01CR02`
#' - `atc_group1`\cr Official pharmacological subgroup (3rd level ATC code) as defined by the WHOCC, like `"Macrolides, lincosamides and streptogramins"`
#' - `atc_group2`\cr Official chemical subgroup (4th level ATC code) as defined by the WHOCC, like `"Macrolides"`
#' - `abbr`\cr List of abbreviations as used in many countries, also for antibiotic susceptibility testing (AST)
#' - `synonyms`\cr Synonyms (often trade names) of a drug, as found in PubChem based on their compound ID
#' - `oral_ddd`\cr Defined Daily Dose (DDD), oral treatment, currently available for `r sum(!is.na(antibiotics$oral_ddd))` drugs
#' - `oral_units`\cr Units of `oral_ddd`
#' - `iv_ddd`\cr Defined Daily Dose (DDD), parenteral (intravenous) treatment, currently available for `r sum(!is.na(antibiotics$iv_ddd))` drugs
#' - `iv_units`\cr Units of `iv_ddd`
#' - `loinc`\cr All LOINC codes (Logical Observation Identifiers Names and Codes) associated with the name of the antimicrobial drug. Use [ab_loinc()] to retrieve them quickly, see [ab_property()].
#'
#' ### For the [antivirals] data set: a [tibble][tibble::tibble] with `r nrow(antivirals)` observations and `r ncol(antivirals)` variables:
#' - `av`\cr Antibiotic ID as used in this package (such as `AMC`), using the official EARS-Net (European Antimicrobial Resistance Surveillance Network) codes where available
#' - `name`\cr Official name as used by WHONET/EARS-Net or the WHO
#' - `atc`\cr ATC codes (Anatomical Therapeutic Chemical) as defined by the WHOCC
#' - `cid`\cr Compound ID as found in PubChem
#' - `atc_group`\cr Official pharmacological subgroup (3rd level ATC code) as defined by the WHOCC
#' - `synonyms`\cr Synonyms (often trade names) of a drug, as found in PubChem based on their compound ID
#' - `oral_ddd`\cr Defined Daily Dose (DDD), oral treatment
#' - `oral_units`\cr Units of `oral_ddd`
#' - `iv_ddd`\cr Defined Daily Dose (DDD), parenteral treatment
#' - `iv_units`\cr Units of `iv_ddd`
#' - `loinc`\cr All LOINC codes (Logical Observation Identifiers Names and Codes) associated with the name of the antimicrobial drug.
#' @details Properties that are based on an ATC code are only available when an ATC is available. These properties are: `atc_group1`, `atc_group2`, `oral_ddd`, `oral_units`, `iv_ddd` and `iv_units`.
#'
#' Synonyms (i.e. trade names) were derived from the PubChem Compound ID (column `cid`) and consequently only available where a CID is available.
#'
#' ### Direct download
#' Like all data sets in this package, these data sets are publicly available for download in the following formats: R, MS Excel, Apache Feather, Apache Parquet, SPSS, SAS, and Stata. Please visit [our website for the download links](https://msberends.github.io/AMR/articles/datasets.html). The actual files are of course available on [our GitHub repository](https://github.com/msberends/AMR/tree/main/data-raw).
#' @source
#'
#' * World Health Organization (WHO) Collaborating Centre for Drug Statistics Methodology (WHOCC): <https://www.whocc.no/atc_ddd_index/>
#'
#' * `r TAXONOMY_VERSION$LOINC$citation` Accessed from <`r TAXONOMY_VERSION$LOINC$url`> on `r documentation_date(TAXONOMY_VERSION$LOINC$accessed_date)`.
#'
#' * European Commission Public Health PHARMACEUTICALS - COMMUNITY REGISTER: <https://ec.europa.eu/health/documents/community-register/html/reg_hum_atc.htm>
#' @inheritSection WHOCC WHOCC
#' @seealso [microorganisms], [intrinsic_resistant]
#' @examples
#' antibiotics
#' antivirals
"antibiotics"

#' @rdname antibiotics
"antivirals"

#' Data Set with `r format(nrow(microorganisms), big.mark = " ")` Microorganisms
#'
#' A data set containing the full microbial taxonomy (**last updated: `r documentation_date(max(TAXONOMY_VERSION$GBIF$accessed_date, TAXONOMY_VERSION$LPSN$accessed_date))`**) of `r nr2char(length(unique(microorganisms$kingdom[!microorganisms$kingdom %like% "unknown"])))` kingdoms from the List of Prokaryotic names with Standing in Nomenclature (LPSN) and the Global Biodiversity Information Facility (GBIF). This data set is the backbone of this `AMR` package. MO codes can be looked up using [as.mo()].
#' @format A [tibble][tibble::tibble] with `r format(nrow(microorganisms), big.mark = " ")` observations and `r ncol(microorganisms)` variables:
#' - `mo`\cr ID of microorganism as used by this package
#' - `fullname`\cr Full name, like `"Escherichia coli"`. For the taxonomic ranks genus, species and subspecies, this is the 'pasted' text of genus, species, and subspecies. For all taxonomic ranks higher than genus, this is the name of the taxon.
#' - `status` \cr Status of the taxon, either `r vector_or(microorganisms$status)`
#' - `kingdom`, `phylum`, `class`, `order`, `family`, `genus`, `species`, `subspecies`\cr Taxonomic rank of the microorganism
#' - `rank`\cr Text of the taxonomic rank of the microorganism, such as `"species"` or `"genus"`
#' - `ref`\cr Author(s) and year of related scientific publication. This contains only the *first surname* and year of the *latest* authors, e.g. "Wallis *et al.* 2006 *emend.* Smith and Jones 2018" becomes "Smith *et al.*, 2018". This field is directly retrieved from the source specified in the column `source`. Moreover, accents were removed to comply with CRAN that only allows ASCII characters, e.g. "V`r "\u00e1\u0148ov\u00e1"`" becomes "Vanova".
#' - `lpsn`\cr Identifier ('Record number') of the List of Prokaryotic names with Standing in Nomenclature (LPSN). This will be the first/highest LPSN identifier to keep one identifier per row. For example, *Acetobacter ascendens* has LPSN Record number 7864 and 11011. Only the first is available in the `microorganisms` data set.
#' - `oxygen_tolerance` \cr Oxygen tolerance, either `r vector_or(microorganisms$oxygen_tolerance)`. These data were retrieved from BacDive (see *Source*). Items that contain "likely" are missing from BacDive and were extrapolated from other species within the same genus to guess the oxygen tolerance. Currently `r round(length(microorganisms$oxygen_tolerance[which(!is.na(microorganisms$oxygen_tolerance))]) / nrow(microorganisms[which(microorganisms$kingdom == "Bacteria"), ]) * 100, 1)`% of all `r format_included_data_number(nrow(microorganisms[which(microorganisms$kingdom == "Bacteria"), ]))` bacteria in the data set contain an oxygen tolerance.
#' - `lpsn_parent`\cr LPSN identifier of the parent taxon
#' - `lpsn_renamed_to`\cr LPSN identifier of the currently valid taxon
#' - `gbif`\cr Identifier ('taxonID') of the Global Biodiversity Information Facility (GBIF)
#' - `gbif_parent`\cr GBIF identifier of the parent taxon
#' - `gbif_renamed_to`\cr GBIF identifier of the currently valid taxon
#' - `source`\cr Either `r vector_or(microorganisms$source)` (see *Source*)
#' - `prevalence`\cr Prevalence of the microorganism according to Bartlett *et al.* (2022, \doi{10.1099/mic.0.001269}), see [mo_matching_score()] for the full explanation
#' - `snomed`\cr Systematized Nomenclature of Medicine (SNOMED) code of the microorganism, version of `r documentation_date(TAXONOMY_VERSION$SNOMED$accessed_date)` (see *Source*). Use [mo_snomed()] to retrieve it quickly, see [mo_property()].
#' @details
#' Please note that entries are only based on the List of Prokaryotic names with Standing in Nomenclature (LPSN) and the Global Biodiversity Information Facility (GBIF) (see below). Since these sources incorporate entries based on (recent) publications in the International Journal of Systematic and Evolutionary Microbiology (IJSEM), it can happen that the year of publication is sometimes later than one might expect.
#'
#' For example, *Staphylococcus pettenkoferi* was described for the first time in Diagnostic Microbiology and Infectious Disease in 2002 (\doi{10.1016/s0732-8893(02)00399-1}), but it was not before 2007 that a publication in IJSEM followed (\doi{10.1099/ijs.0.64381-0}). Consequently, the `AMR` package returns 2007 for `mo_year("S. pettenkoferi")`.
#'
#' @section Included Taxa:
#' Included taxonomic data are:
#' - All `r format_included_data_number(microorganisms[which(microorganisms$kingdom %in% c("Archeae", "Bacteria")), , drop = FALSE])` (sub)species from the kingdoms of Archaea and Bacteria
#' - `r format_included_data_number(microorganisms[which(microorganisms$kingdom == "Fungi"), , drop = FALSE])` (sub)species from the kingdom of Fungi. The kingdom of Fungi is a very large taxon with almost 300,000 different (sub)species, of which most are not microbial (but rather macroscopic, like mushrooms). Because of this, not all fungi fit the scope of this package. Only relevant fungi are covered (such as all species of *Aspergillus*, *Candida*, *Cryptococcus*, *Histoplasma*, *Pneumocystis*, *Saccharomyces* and *Trichophyton*).
#' - `r format_included_data_number(microorganisms[which(microorganisms$kingdom == "Protozoa"), , drop = FALSE])` (sub)species from the kingdom of Protozoa
#' - `r format_included_data_number(microorganisms[which(microorganisms$kingdom == "Animalia"), , drop = FALSE])` (sub)species from `r format_included_data_number(microorganisms[which(microorganisms$kingdom == "Animalia"), "genus", drop = TRUE])` other relevant genera from the kingdom of Animalia (such as *Strongyloides* and *Taenia*)
#' - All `r format_included_data_number(microorganisms[which(microorganisms$status != "accepted"), , drop = FALSE])` previously accepted names of all included (sub)species (these were taxonomically renamed)
#' - The complete taxonomic tree of all included (sub)species: from kingdom to subspecies
#' - The identifier of the parent taxons
#' - The year and first author of the related scientific publication
#'
#' ### Manual additions
#' For convenience, some entries were added manually:
#'
#' - `r format_included_data_number(microorganisms[which(microorganisms$source == "manually added" & microorganisms$genus == "Salmonella"), , drop = FALSE])` entries of *Salmonella*, such as the city-like serovars and groups A to H
#' - `r format_included_data_number(microorganisms[which(microorganisms$source == "manually added" & microorganisms$genus == "Streptococcus"), , drop = FALSE])` entries of *Streptococcus*, such as the beta-haemolytic groups A to K, viridans, and milleri
#' - 2 entries of *Staphylococcus* (coagulase-negative (CoNS) and coagulase-positive (CoPS))
#' - 1 entry of *Blastocystis* (*B.  hominis*), although it officially does not exist (Noel *et al.* 2005, PMID 15634993)
#' - 1 entry of *Moraxella* (*M. catarrhalis*), which was formally named *Branhamella catarrhalis* (Catlin, 1970) though this change was never accepted within the field of clinical microbiology
#' - 8 other 'undefined' entries (unknown, unknown Gram-negatives, unknown Gram-positives, unknown yeast, unknown fungus, and unknown anaerobic Gram-pos/Gram-neg bacteria)
#'
#' The syntax used to transform the original data to a cleansed \R format, can be found here: <https://github.com/msberends/AMR/blob/main/data-raw/reproduction_of_microorganisms.R>.
#'
#' ### Direct download
#' Like all data sets in this package, this data set is publicly available for download in the following formats: R, MS Excel, Apache Feather, Apache Parquet, SPSS, SAS, and Stata. Please visit [our website for the download links](https://msberends.github.io/AMR/articles/datasets.html). The actual files are of course available on [our GitHub repository](https://github.com/msberends/AMR/tree/main/data-raw).
#' @section About the Records from LPSN (see *Source*):
#' LPSN is the main source for bacteriological taxonomy of this `AMR` package.
#'
#' The List of Prokaryotic names with Standing in Nomenclature (LPSN) provides comprehensive information on the nomenclature of prokaryotes. LPSN is a free to use service founded by Jean P. Euzeby in 1997 and later on maintained by Aidan C. Parte.
#' @source
#' * `r TAXONOMY_VERSION$LPSN$citation` Accessed from <`r TAXONOMY_VERSION$LPSN$url`> on `r documentation_date(TAXONOMY_VERSION$LPSN$accessed_date)`.
#'
#' * `r TAXONOMY_VERSION$GBIF$citation` Accessed from <`r TAXONOMY_VERSION$GBIF$url`> on `r documentation_date(TAXONOMY_VERSION$GBIF$accessed_date)`.
#'
#' * `r TAXONOMY_VERSION$SNOMED$citation` URL: <`r TAXONOMY_VERSION$SNOMED$url`>
#'
#' * Grimont *et al.* (2007). Antigenic Formulae of the Salmonella Serovars, 9th Edition. WHO Collaborating Centre for Reference and Research on *Salmonella* (WHOCC-SALM).
#'
#' * Bartlett *et al.* (2022). **A comprehensive list of bacterial pathogens infecting humans** *Microbiology* 168:001269; \doi{10.1099/mic.0.001269}
#' 
#' * Reimer *et al.* (2022). ***BacDive* in 2022: the knowledge base for standardized bacterial and archaeal data.**  *Nucleic Acids Res.* 2022 Jan 7;50(D1):D741-D746; \doi{10.1093/nar/gkab961}
#' @seealso [as.mo()], [mo_property()], [microorganisms.codes], [intrinsic_resistant]
#' @examples
#' microorganisms
"microorganisms"

#' Data Set with `r format(nrow(microorganisms.codes), big.mark = " ")` Common Microorganism Codes
#'
#' A data set containing commonly used codes for microorganisms, from laboratory systems and WHONET. Define your own with [set_mo_source()]. They will all be searched when using [as.mo()] and consequently all the [`mo_*`][mo_property()] functions.
#' @format A [tibble][tibble::tibble] with `r format(nrow(microorganisms.codes), big.mark = " ")` observations and `r ncol(microorganisms.codes)` variables:
#' - `code`\cr Commonly used code of a microorganism
#' - `mo`\cr ID of the microorganism in the [microorganisms] data set
#' @details
#' Like all data sets in this package, this data set is publicly available for download in the following formats: R, MS Excel, Apache Feather, Apache Parquet, SPSS, SAS, and Stata. Please visit [our website for the download links](https://msberends.github.io/AMR/articles/datasets.html). The actual files are of course available on [our GitHub repository](https://github.com/msberends/AMR/tree/main/data-raw).
#' @seealso [as.mo()] [microorganisms]
#' @examples
#' microorganisms.codes
"microorganisms.codes"

#' Data Set with `r format(nrow(example_isolates), big.mark = " ")` Example Isolates
#'
#' A data set containing `r format(nrow(example_isolates), big.mark = " ")` microbial isolates with their full antibiograms. This data set contains randomised fictitious data, but reflects reality and can be used to practise AMR data analysis. For examples, please read [the tutorial on our website](https://msberends.github.io/AMR/articles/AMR.html).
#' @format A [tibble][tibble::tibble] with `r format(nrow(example_isolates), big.mark = " ")` observations and `r ncol(example_isolates)` variables:
#' - `date`\cr Date of receipt at the laboratory
#' - `patient`\cr ID of the patient
#' - `age`\cr Age of the patient
#' - `gender`\cr Gender of the patient, either `r vector_or(example_isolates$gender)`
#' - `ward`\cr Ward type where the patient was admitted, either `r vector_or(example_isolates$ward)`
#' - `mo`\cr ID of microorganism created with [as.mo()], see also the [microorganisms] data set
#' - `PEN:RIF`\cr `r sum(vapply(FUN.VALUE = logical(1), example_isolates, is.sir))` different antibiotics with class [`sir`] (see [as.sir()]); these column names occur in the [antibiotics] data set and can be translated with [set_ab_names()] or [ab_name()]
#' @details
#' Like all data sets in this package, this data set is publicly available for download in the following formats: R, MS Excel, Apache Feather, Apache Parquet, SPSS, SAS, and Stata. Please visit [our website for the download links](https://msberends.github.io/AMR/articles/datasets.html). The actual files are of course available on [our GitHub repository](https://github.com/msberends/AMR/tree/main/data-raw).
#' @examples
#' example_isolates
"example_isolates"

#' Data Set with Unclean Data
#'
#' A data set containing `r format(nrow(example_isolates_unclean), big.mark = " ")` microbial isolates that are not cleaned up and consequently not ready for AMR data analysis. This data set can be used for practice.
#' @format A [tibble][tibble::tibble] with `r format(nrow(example_isolates_unclean), big.mark = " ")` observations and `r ncol(example_isolates_unclean)` variables:
#' - `patient_id`\cr ID of the patient
#' - `date`\cr date of receipt at the laboratory
#' - `hospital`\cr ID of the hospital, from A to C
#' - `bacteria`\cr info about microorganism that can be transformed with [as.mo()], see also [microorganisms]
#' - `AMX:GEN`\cr 4 different antibiotics that have to be transformed with [as.sir()]
#' @details
#' Like all data sets in this package, this data set is publicly available for download in the following formats: R, MS Excel, Apache Feather, Apache Parquet, SPSS, SAS, and Stata. Please visit [our website for the download links](https://msberends.github.io/AMR/articles/datasets.html). The actual files are of course available on [our GitHub repository](https://github.com/msberends/AMR/tree/main/data-raw).
#' @examples
#' example_isolates_unclean
"example_isolates_unclean"

#' Data Set with `r format(nrow(WHONET), big.mark = " ")` Isolates - WHONET Example
#'
#' This example data set has the exact same structure as an export file from WHONET. Such files can be used with this package, as this example data set shows. The antibiotic results are from our [example_isolates] data set. All patient names were created using online surname generators and are only in place for practice purposes.
#' @format A [tibble][tibble::tibble] with `r format(nrow(WHONET), big.mark = " ")` observations and `r ncol(WHONET)` variables:
#' - `Identification number`\cr ID of the sample
#' - `Specimen number`\cr ID of the specimen
#' - `Organism`\cr Name of the microorganism. Before analysis, you should transform this to a valid microbial class, using [as.mo()].
#' - `Country`\cr Country of origin
#' - `Laboratory`\cr Name of laboratory
#' - `Last name`\cr Fictitious last name of patient
#' - `First name`\cr Fictitious initial of patient
#' - `Sex`\cr Fictitious gender of patient
#' - `Age`\cr Fictitious age of patient
#' - `Age category`\cr Age group, can also be looked up using [age_groups()]
#' - `Date of admission`\cr [Date] of hospital admission
#' - `Specimen date`\cr [Date] when specimen was received at laboratory
#' - `Specimen type`\cr Specimen type or group
#' - `Specimen type (Numeric)`\cr Translation of `"Specimen type"`
#' - `Reason`\cr Reason of request with Differential Diagnosis
#' - `Isolate number`\cr ID of isolate
#' - `Organism type`\cr Type of microorganism, can also be looked up using [mo_type()]
#' - `Serotype`\cr Serotype of microorganism
#' - `Beta-lactamase`\cr Microorganism produces beta-lactamase?
#' - `ESBL`\cr Microorganism produces extended spectrum beta-lactamase?
#' - `Carbapenemase`\cr Microorganism produces carbapenemase?
#' - `MRSA screening test`\cr Microorganism is possible MRSA?
#' - `Inducible clindamycin resistance`\cr Clindamycin can be induced?
#' - `Comment`\cr Other comments
#' - `Date of data entry`\cr [Date] this data was entered in WHONET
#' - `AMP_ND10:CIP_EE`\cr `r sum(vapply(FUN.VALUE = logical(1), WHONET, is.sir))` different antibiotics. You can lookup the abbreviations in the [antibiotics] data set, or use e.g. [`ab_name("AMP")`][ab_name()] to get the official name immediately. Before analysis, you should transform this to a valid antibiotic class, using [as.sir()].
#' @details
#' Like all data sets in this package, this data set is publicly available for download in the following formats: R, MS Excel, Apache Feather, Apache Parquet, SPSS, SAS, and Stata. Please visit [our website for the download links](https://msberends.github.io/AMR/articles/datasets.html). The actual files are of course available on [our GitHub repository](https://github.com/msberends/AMR/tree/main/data-raw).
#' @examples
#' WHONET
"WHONET"

#' Data Set with Clinical Breakpoints for SIR Interpretation
#'
#' Data set containing clinical breakpoints to interpret MIC and disk diffusion to SIR values, according to international guidelines. Currently implemented guidelines are EUCAST (`r min(as.integer(gsub("[^0-9]", "", subset(clinical_breakpoints, guideline %like% "EUCAST")$guideline)))`-`r max(as.integer(gsub("[^0-9]", "", subset(clinical_breakpoints, guideline %like% "EUCAST")$guideline)))`) and CLSI (`r min(as.integer(gsub("[^0-9]", "", subset(clinical_breakpoints, guideline %like% "CLSI")$guideline)))`-`r max(as.integer(gsub("[^0-9]", "", subset(clinical_breakpoints, guideline %like% "CLSI")$guideline)))`). Use [as.sir()] to transform MICs or disks measurements to SIR values.
#' @format A [tibble][tibble::tibble] with `r format(nrow(clinical_breakpoints), big.mark = " ")` observations and `r ncol(clinical_breakpoints)` variables:
#' - `guideline`\cr Name of the guideline
#' - `method`\cr Either `r vector_or(clinical_breakpoints$method)`
#' - `site`\cr Body site, e.g. "Oral" or "Respiratory"
#' - `mo`\cr Microbial ID, see [as.mo()]
#' - `rank_index`\cr Taxonomic rank index of `mo` from 1 (subspecies/infraspecies) to 5 (unknown microorganism)
#' - `ab`\cr Antibiotic ID, see [as.ab()]
#' - `ref_tbl`\cr Info about where the guideline rule can be found
#' - `disk_dose`\cr Dose of the used disk diffusion method
#' - `breakpoint_S`\cr Lowest MIC value or highest number of millimetres that leads to "S"
#' - `breakpoint_R`\cr Highest MIC value or lowest number of millimetres that leads to "R"
#' - `uti`\cr A [logical] value (`TRUE`/`FALSE`) to indicate whether the rule applies to a urinary tract infection (UTI)
#' @details
#' Like all data sets in this package, this data set is publicly available for download in the following formats: R, MS Excel, Apache Feather, Apache Parquet, SPSS, SAS, and Stata. Please visit [our website for the download links](https://msberends.github.io/AMR/articles/datasets.html). The actual files are of course available on [our GitHub repository](https://github.com/msberends/AMR/tree/main/data-raw).
#'
#' They **allow for machine reading EUCAST and CLSI guidelines**, which is almost impossible with the MS Excel and PDF files distributed by EUCAST and CLSI.
#' @seealso [intrinsic_resistant]
#' @examples
#' clinical_breakpoints
"clinical_breakpoints"

#' Data Set with Bacterial Intrinsic Resistance
#'
#' Data set containing defined intrinsic resistance by EUCAST of all bug-drug combinations.
#' @format A [tibble][tibble::tibble] with `r format(nrow(intrinsic_resistant), big.mark = " ")` observations and `r ncol(intrinsic_resistant)` variables:
#' - `mo`\cr Microorganism ID
#' - `ab`\cr Antibiotic ID
#' @details
#' This data set is based on `r format_eucast_version_nr(3.3)`.
#'
#' ### Direct download
#' Like all data sets in this package, this data set is publicly available for download in the following formats: R, MS Excel, Apache Feather, Apache Parquet, SPSS, SAS, and Stata. Please visit [our website for the download links](https://msberends.github.io/AMR/articles/datasets.html). The actual files are of course available on [our GitHub repository](https://github.com/msberends/AMR/tree/main/data-raw).
#'
#' They **allow for machine reading EUCAST and CLSI guidelines**, which is almost impossible with the MS Excel and PDF files distributed by EUCAST and CLSI.
#' @examples
#' intrinsic_resistant
"intrinsic_resistant"

#' Data Set with Treatment Dosages as Defined by EUCAST
#'
#' EUCAST breakpoints used in this package are based on the dosages in this data set. They can be retrieved with [eucast_dosage()].
#' @format A [tibble][tibble::tibble] with `r format(nrow(dosage), big.mark = " ")` observations and `r ncol(dosage)` variables:
#' - `ab`\cr Antibiotic ID as used in this package (such as `AMC`), using the official EARS-Net (European Antimicrobial Resistance Surveillance Network) codes where available
#' - `name`\cr Official name of the antimicrobial drug as used by WHONET/EARS-Net or the WHO
#' - `type`\cr Type of the dosage, either `r vector_or(dosage$type)`
#' - `dose`\cr Dose, such as "2 g" or "25 mg/kg"
#' - `dose_times`\cr Number of times a dose must be administered
#' - `administration`\cr Route of administration, either `r vector_or(dosage$administration)`
#' - `notes`\cr Additional dosage notes
#' - `original_txt`\cr Original text in the PDF file of EUCAST
#' - `eucast_version`\cr Version number of the EUCAST Clinical Breakpoints guideline to which these dosages apply
#' @details
#' This data set is based on `r format_eucast_version_nr(12.0)` and `r format_eucast_version_nr(11.0)`.
#'
#' ### Direct download
#' Like all data sets in this package, this data set is publicly available for download in the following formats: R, MS Excel, Apache Feather, Apache Parquet, SPSS, SAS, and Stata. Please visit [our website for the download links](https://msberends.github.io/AMR/articles/datasets.html). The actual files are of course available on [our GitHub repository](https://github.com/msberends/AMR/tree/main/data-raw).
#' @examples
#' dosage
"dosage"
