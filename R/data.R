# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Data Analysis for R                   #
#                                                                      #
# SOURCE                                                               #
# https://github.com/msberends/AMR                                     #
#                                                                      #
# LICENCE                                                              #
# (c) 2018-2022 Berends MS, Luz CF et al.                              #
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

#' Data Sets with `r format(nrow(antibiotics) + nrow(antivirals), big.mark = ",")` Antimicrobial Drugs
#'
#' Two data sets containing all antibiotics/antimycotics and antivirals. Use [as.ab()] or one of the [`ab_*`][ab_property()] functions to retrieve values from the [antibiotics] data set. Three identifiers are included in this data set: an antibiotic ID (`ab`, primarily used in this package) as defined by WHONET/EARS-Net, an ATC code (`atc`) as defined by the WHO, and a Compound ID (`cid`) as found in PubChem. Other properties in this data set are derived from one or more of these codes. Note that some drugs have multiple ATC codes.
#' @format
#' ## For the [antibiotics] data set: a [tibble[tibble::tibble] with `r nrow(antibiotics)` observations and `r ncol(antibiotics)` variables:
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
#' - `loinc`\cr All LOINC codes (Logical Observation Identifiers Names and Codes) associated with the name of the antimicrobial agent. Use [ab_loinc()] to retrieve them quickly, see [ab_property()].
#' 
#' ## For the [antivirals] data set: a [tibble[tibble::tibble] with `r nrow(antivirals)` observations and `r ncol(antivirals)` variables:
#' - `atc`\cr ATC codes (Anatomical Therapeutic Chemical) as defined by the WHOCC
#' - `cid`\cr Compound ID as found in PubChem
#' - `name`\cr Official name as used by WHONET/EARS-Net or the WHO
#' - `atc_group`\cr Official pharmacological subgroup (3rd level ATC code) as defined by the WHOCC
#' - `synonyms`\cr Synonyms (often trade names) of a drug, as found in PubChem based on their compound ID
#' - `oral_ddd`\cr Defined Daily Dose (DDD), oral treatment
#' - `oral_units`\cr Units of `oral_ddd`
#' - `iv_ddd`\cr Defined Daily Dose (DDD), parenteral treatment
#' - `iv_units`\cr Units of `iv_ddd`
#' @details Properties that are based on an ATC code are only available when an ATC is available. These properties are: `atc_group1`, `atc_group2`, `oral_ddd`, `oral_units`, `iv_ddd` and `iv_units`.
#'
#' Synonyms (i.e. trade names) were derived from the Compound ID (`cid`) and consequently only available where a CID is available.
#' 
#' ## Direct download
#' Like all data sets in this package, these data sets are publicly available for download in the following formats: R, MS Excel, Apache Feather, Apache Parquet, SPSS, SAS, and Stata. Please visit [our website for the download links](https://msberends.github.io/AMR/articles/datasets.html). The actual files are of course available on [our GitHub repository](https://github.com/msberends/AMR/tree/main/data-raw).
#' @source World Health Organization (WHO) Collaborating Centre for Drug Statistics Methodology (WHOCC): <https://www.whocc.no/atc_ddd_index/>
#'
#' European Commission Public Health PHARMACEUTICALS - COMMUNITY REGISTER: <https://ec.europa.eu/health/documents/community-register/html/reg_hum_atc.htm>
#' @inheritSection WHOCC WHOCC
#' @seealso [microorganisms], [intrinsic_resistant]
#' @examples
#' antibiotics
#' antivirals
"antibiotics"

#' @rdname antibiotics
"antivirals"

#' Data Set with `r format(nrow(microorganisms), big.mark = ",")` Microorganisms
#'
#' A data set containing the full microbial taxonomy (**last updated: `r CATALOGUE_OF_LIFE$yearmonth_LPSN`**) of `r nr2char(length(unique(microorganisms$kingdom[!microorganisms$kingdom %like% "unknown"])))` kingdoms from the Catalogue of Life (CoL) and the List of Prokaryotic names with Standing in Nomenclature (LPSN). MO codes can be looked up using [as.mo()].
#' @inheritSection catalogue_of_life Catalogue of Life
#' @format A [tibble[tibble::tibble] with `r format(nrow(microorganisms), big.mark = ",")` observations and `r ncol(microorganisms)` variables:
#' - `mo`\cr ID of microorganism as used by this package
#' - `fullname`\cr Full name, like `"Escherichia coli"`
#' - `kingdom`, `phylum`, `class`, `order`, `family`, `genus`, `species`, `subspecies`\cr Taxonomic rank of the microorganism
#' - `rank`\cr Text of the taxonomic rank of the microorganism, like `"species"` or `"genus"`
#' - `ref`\cr Author(s) and year of concerning scientific publication
#' - `species_id`\cr ID of the species as used by the Catalogue of Life
#' - `source`\cr Either `r vector_or(microorganisms$source)` (see *Source*)
#' - `prevalence`\cr Prevalence of the microorganism, see [as.mo()]
#' - `snomed`\cr Systematized Nomenclature of Medicine (SNOMED) code of the microorganism, according to the `r SNOMED_VERSION$current_source` (see *Source*). Use [mo_snomed()] to retrieve it quickly, see [mo_property()].
#' @details 
#' Please note that entries are only based on the Catalogue of Life and the LPSN (see below). Since these sources incorporate entries based on (recent) publications in the International Journal of Systematic and Evolutionary Microbiology (IJSEM), it can happen that the year of publication is sometimes later than one might expect.
#' 
#' For example, *Staphylococcus pettenkoferi* was described for the first time in Diagnostic Microbiology and Infectious Disease in 2002 (\doi{10.1016/s0732-8893(02)00399-1}), but it was not before 2007 that a publication in IJSEM followed (\doi{10.1099/ijs.0.64381-0}). Consequently, the `AMR` package returns 2007 for `mo_year("S. pettenkoferi")`.
#' 
#' ## Manual additions
#' For convenience, some entries were added manually:
#' 
#' - 11 entries of *Streptococcus* (beta-haemolytic: groups A, B, C, D, F, G, H, K and unspecified; other: viridans, milleri)
#' - 2 entries of *Staphylococcus* (coagulase-negative (CoNS) and coagulase-positive (CoPS))
#' - 3 entries of *Trichomonas* (*T. vaginalis*, and its family and genus)
#' - 4 entries of *Toxoplasma* (*T. gondii*, and its order, family and genus)
#' - 1 entry of *Candida* (*C.  krusei*), that is not (yet) in the Catalogue of Life
#' - 1 entry of *Blastocystis* (*B.  hominis*), although it officially does not exist (Noel *et al.* 2005, PMID 15634993)
#' - 1 entry of *Moraxella* (*M. catarrhalis*), which was formally named *Branhamella catarrhalis* (Catlin, 1970) though this change was never accepted within the field of clinical microbiology
#' - 5 other 'undefined' entries (unknown, unknown Gram negatives, unknown Gram positives, unknown yeast and unknown fungus)
#' - 6 families under the Enterobacterales order, according to Adeolu *et al.* (2016, PMID 27620848), that are not (yet) in the Catalogue of Life
#' 
#' ## Direct download
#' Like all data sets in this package, this data set is publicly available for download in the following formats: R, MS Excel, Apache Feather, Apache Parquet, SPSS, SAS, and Stata. Please visit [our website for the download links](https://msberends.github.io/AMR/articles/datasets.html). The actual files are of course available on [our GitHub repository](https://github.com/msberends/AMR/tree/main/data-raw).
#' @section About the Records from LPSN (see *Source*):
#' The List of Prokaryotic names with Standing in Nomenclature (LPSN) provides comprehensive information on the nomenclature of prokaryotes. LPSN is a free to use service founded by Jean P. Euzeby in 1997 and later on maintained by Aidan C. Parte.
#' 
#' As of February 2020, the regularly augmented LPSN database at DSMZ is the basis of the new LPSN service. The new database was implemented for the Type-Strain Genome Server and augmented in 2018 to store all kinds of nomenclatural information. Data from the previous version of LPSN and from the Prokaryotic Nomenclature Up-to-date (PNU) service were imported into the new system. PNU had been established in 1993 as a service of the Leibniz Institute DSMZ, and was curated by Norbert Weiss, Manfred Kracht and Dorothea Gleim.
#' @source 
#' `r gsub("{year}", CATALOGUE_OF_LIFE$year, CATALOGUE_OF_LIFE$version, fixed = TRUE)` as currently implemented in this `AMR` package:
#' 
#' * Annual Checklist (public online taxonomic database), <http://www.catalogueoflife.org>
#' 
#' List of Prokaryotic names with Standing in Nomenclature (`r CATALOGUE_OF_LIFE$yearmonth_LPSN`) as currently implemented in this `AMR` package:
#' 
#' * Parte, A.C., Sarda Carbasse, J., Meier-Kolthoff, J.P., Reimer, L.C. and Goker, M. (2020). List of Prokaryotic names with Standing in Nomenclature (LPSN) moves to the DSMZ. International Journal of Systematic and Evolutionary Microbiology, 70, 5607-5612; \doi{10.1099/ijsem.0.004332}
#' * Parte, A.C. (2018). LPSN - List of Prokaryotic names with Standing in Nomenclature (bacterio.net), 20 years on. International Journal of Systematic and Evolutionary Microbiology, 68, 1825-1829; \doi{10.1099/ijsem.0.002786}
#' * Parte, A.C. (2014). LPSN - List of Prokaryotic names with Standing in Nomenclature. Nucleic Acids Research, 42, Issue D1, D613-D616; \doi{10.1093/nar/gkt1111}
#' * Euzeby, J.P. (1997). List of Bacterial Names with Standing in Nomenclature: a Folder Available on the Internet. International Journal of Systematic Bacteriology, 47, 590-592; \doi{10.1099/00207713-47-2-590}
#' 
#' `r SNOMED_VERSION$current_source` as currently implemented in this `AMR` package:
#' 
#' * Retrieved from the `r SNOMED_VERSION$title`, OID `r SNOMED_VERSION$current_oid`, version `r SNOMED_VERSION$current_version`; url: <`r SNOMED_VERSION$url`>
#' @seealso [as.mo()], [mo_property()], [microorganisms.codes], [intrinsic_resistant]
#' @examples
#' microorganisms
"microorganisms"

#' Data Set with Previously Accepted Taxonomic Names
#'
#' A data set containing old (previously valid or accepted) taxonomic names according to the Catalogue of Life. This data set is used internally by [as.mo()].
#' @inheritSection catalogue_of_life Catalogue of Life
#' @format A [tibble[tibble::tibble] with `r format(nrow(microorganisms.old), big.mark = ",")` observations and `r ncol(microorganisms.old)` variables:
#' - `fullname`\cr Old full taxonomic name of the microorganism
#' - `fullname_new`\cr New full taxonomic name of the microorganism
#' - `ref`\cr Author(s) and year of concerning scientific publication
#' - `prevalence`\cr Prevalence of the microorganism, see [as.mo()]
#' @details 
#' Like all data sets in this package, this data set is publicly available for download in the following formats: R, MS Excel, Apache Feather, Apache Parquet, SPSS, SAS, and Stata. Please visit [our website for the download links](https://msberends.github.io/AMR/articles/datasets.html). The actual files are of course available on [our GitHub repository](https://github.com/msberends/AMR/tree/main/data-raw).
#' @source Catalogue of Life: Annual Checklist (public online taxonomic database), <http://www.catalogueoflife.org> (check included annual version with [catalogue_of_life_version()]).
#' 
#' Parte, A.C. (2018). LPSN - List of Prokaryotic names with Standing in Nomenclature (bacterio.net), 20 years on. International Journal of Systematic and Evolutionary Microbiology, 68, 1825-1829; \doi{10.1099/ijsem.0.002786}
#' @seealso [as.mo()] [mo_property()] [microorganisms]
#' @examples
#' microorganisms.old
"microorganisms.old"

#' Data Set with `r format(nrow(microorganisms.codes), big.mark = ",")` Common Microorganism Codes
#'
#' A data set containing commonly used codes for microorganisms, from laboratory systems and WHONET. Define your own with [set_mo_source()]. They will all be searched when using [as.mo()] and consequently all the [`mo_*`][mo_property()] functions.
#' @format A [tibble[tibble::tibble] with `r format(nrow(microorganisms.codes), big.mark = ",")` observations and `r ncol(microorganisms.codes)` variables:
#' - `code`\cr Commonly used code of a microorganism
#' - `mo`\cr ID of the microorganism in the [microorganisms] data set
#' @details 
#' Like all data sets in this package, this data set is publicly available for download in the following formats: R, MS Excel, Apache Feather, Apache Parquet, SPSS, SAS, and Stata. Please visit [our website for the download links](https://msberends.github.io/AMR/articles/datasets.html). The actual files are of course available on [our GitHub repository](https://github.com/msberends/AMR/tree/main/data-raw).
#' @inheritSection catalogue_of_life Catalogue of Life
#' @seealso [as.mo()] [microorganisms]
#' @examples
#' microorganisms.codes
"microorganisms.codes"

#' Data Set with `r format(nrow(example_isolates), big.mark = ",")` Example Isolates
#'
#' A data set containing `r format(nrow(example_isolates), big.mark = ",")` microbial isolates with their full antibiograms. This data set contains randomised fictitious data, but reflects reality and can be used to practise AMR data analysis. For examples, please read [the tutorial on our website](https://msberends.github.io/AMR/articles/AMR.html).
#' @format A [tibble[tibble::tibble] with `r format(nrow(example_isolates), big.mark = ",")` observations and `r ncol(example_isolates)` variables:
#' - `date`\cr Date of receipt at the laboratory
#' - `patient`\cr ID of the patient
#' - `age`\cr Age of the patient
#' - `gender`\cr Gender of the patient, either `r vector_or(example_isolates$gender)`
#' - `ward`\cr Ward type where the patient was admitted, either `r vector_or(example_isolates$ward)`
#' - `mo`\cr ID of microorganism created with [as.mo()], see also the [microorganisms] data set
#' - `PEN:RIF`\cr `r sum(vapply(FUN.VALUE = logical(1), example_isolates, is.rsi))` different antibiotics with class [`rsi`] (see [as.rsi()]); these column names occur in the [antibiotics] data set and can be translated with [set_ab_names()] or [ab_name()]
#' @details 
#' Like all data sets in this package, this data set is publicly available for download in the following formats: R, MS Excel, Apache Feather, Apache Parquet, SPSS, SAS, and Stata. Please visit [our website for the download links](https://msberends.github.io/AMR/articles/datasets.html). The actual files are of course available on [our GitHub repository](https://github.com/msberends/AMR/tree/main/data-raw).
#' @examples
#' example_isolates
"example_isolates"

#' Data Set with Unclean Data
#'
#' A data set containing `r format(nrow(example_isolates_unclean), big.mark = ",")` microbial isolates that are not cleaned up and consequently not ready for AMR data analysis. This data set can be used for practice.
#' @format A [tibble[tibble::tibble] with `r format(nrow(example_isolates_unclean), big.mark = ",")` observations and `r ncol(example_isolates_unclean)` variables:
#' - `patient_id`\cr ID of the patient
#' - `date`\cr date of receipt at the laboratory
#' - `hospital`\cr ID of the hospital, from A to C
#' - `bacteria`\cr info about microorganism that can be transformed with [as.mo()], see also [microorganisms]
#' - `AMX:GEN`\cr 4 different antibiotics that have to be transformed with [as.rsi()]
#' @details 
#' Like all data sets in this package, this data set is publicly available for download in the following formats: R, MS Excel, Apache Feather, Apache Parquet, SPSS, SAS, and Stata. Please visit [our website for the download links](https://msberends.github.io/AMR/articles/datasets.html). The actual files are of course available on [our GitHub repository](https://github.com/msberends/AMR/tree/main/data-raw).
#' @examples
#' example_isolates_unclean
"example_isolates_unclean"

#' Data Set with `r format(nrow(WHONET), big.mark = ",")` Isolates - WHONET Example
#'
#' This example data set has the exact same structure as an export file from WHONET. Such files can be used with this package, as this example data set shows. The antibiotic results are from our [example_isolates] data set. All patient names are created using online surname generators and are only in place for practice purposes.
#' @format A [tibble[tibble::tibble] with `r format(nrow(WHONET), big.mark = ",")` observations and `r ncol(WHONET)` variables:
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
#' - `AMP_ND10:CIP_EE`\cr `r sum(vapply(FUN.VALUE = logical(1), WHONET, is.rsi))` different antibiotics. You can lookup the abbreviations in the [antibiotics] data set, or use e.g. [`ab_name("AMP")`][ab_name()] to get the official name immediately. Before analysis, you should transform this to a valid antibiotic class, using [as.rsi()].
#' @details 
#' Like all data sets in this package, this data set is publicly available for download in the following formats: R, MS Excel, Apache Feather, Apache Parquet, SPSS, SAS, and Stata. Please visit [our website for the download links](https://msberends.github.io/AMR/articles/datasets.html). The actual files are of course available on [our GitHub repository](https://github.com/msberends/AMR/tree/main/data-raw).
#' @examples
#' WHONET
"WHONET"

#' Data Set for R/SI Interpretation
#'
#' Data set containing reference data to interpret MIC and disk diffusion to R/SI values, according to international guidelines. Currently implemented guidelines are EUCAST (`r min(as.integer(gsub("[^0-9]", "", subset(rsi_translation, guideline %like% "EUCAST")$guideline)))`-`r max(as.integer(gsub("[^0-9]", "", subset(rsi_translation, guideline %like% "EUCAST")$guideline)))`) and CLSI (`r min(as.integer(gsub("[^0-9]", "", subset(rsi_translation, guideline %like% "CLSI")$guideline)))`-`r max(as.integer(gsub("[^0-9]", "", subset(rsi_translation, guideline %like% "CLSI")$guideline)))`). Use [as.rsi()] to transform MICs or disks measurements to R/SI values.
#' @format A [tibble[tibble::tibble] with `r format(nrow(rsi_translation), big.mark = ",")` observations and `r ncol(rsi_translation)` variables:
#' - `guideline`\cr Name of the guideline
#' - `method`\cr Either `r vector_or(rsi_translation$method)`
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
#' rsi_translation
"rsi_translation"

#' Data Set with Bacterial Intrinsic Resistance
#'
#' Data set containing defined intrinsic resistance by EUCAST of all bug-drug combinations.
#' @format A [tibble[tibble::tibble] with `r format(nrow(intrinsic_resistant), big.mark = ",")` observations and `r ncol(intrinsic_resistant)` variables:
#' - `mo`\cr Microorganism ID
#' - `ab`\cr Antibiotic ID
#' @details 
#' This data set is based on `r format_eucast_version_nr(3.3)`.
#' 
#' ## Direct download
#' Like all data sets in this package, this data set is publicly available for download in the following formats: R, MS Excel, Apache Feather, Apache Parquet, SPSS, SAS, and Stata. Please visit [our website for the download links](https://msberends.github.io/AMR/articles/datasets.html). The actual files are of course available on [our GitHub repository](https://github.com/msberends/AMR/tree/main/data-raw).
#' 
#' They **allow for machine reading EUCAST and CLSI guidelines**, which is almost impossible with the MS Excel and PDF files distributed by EUCAST and CLSI.
#' @examples
#' intrinsic_resistant
"intrinsic_resistant"

#' Data Set with Treatment Dosages as Defined by EUCAST
#'
#' EUCAST breakpoints used in this package are based on the dosages in this data set. They can be retrieved with [eucast_dosage()].
#' @format A [tibble[tibble::tibble] with `r format(nrow(dosage), big.mark = ",")` observations and `r ncol(dosage)` variables:
#' - `ab`\cr Antibiotic ID as used in this package (such as `AMC`), using the official EARS-Net (European Antimicrobial Resistance Surveillance Network) codes where available
#' - `name`\cr Official name of the antimicrobial agent as used by WHONET/EARS-Net or the WHO
#' - `type`\cr Type of the dosage, either `r vector_or(dosage$type)`
#' - `dose`\cr Dose, such as "2 g" or "25 mg/kg"
#' - `dose_times`\cr Number of times a dose must be administered
#' - `administration`\cr Route of administration, either `r vector_or(dosage$administration)`
#' - `notes`\cr Additional dosage notes
#' - `original_txt`\cr Original text in the PDF file of EUCAST
#' - `eucast_version`\cr Version number of the EUCAST Clinical Breakpoints guideline to which these dosages apply
#' @details 
#' This data set is based on `r format_eucast_version_nr(11.0)`.
#' 
#' ## Direct download
#' Like all data sets in this package, this data set is publicly available for download in the following formats: R, MS Excel, Apache Feather, Apache Parquet, SPSS, SAS, and Stata. Please visit [our website for the download links](https://msberends.github.io/AMR/articles/datasets.html). The actual files are of course available on [our GitHub repository](https://github.com/msberends/AMR/tree/main/data-raw).
#' @examples
#' dosage
"dosage"
