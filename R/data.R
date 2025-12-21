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
# how to conduct AMR data analysis: https://amr-for-r.org              #
# ==================================================================== #

#' Data Sets with `r format(nrow(antimicrobials) + nrow(antivirals), big.mark = " ")` Antimicrobial Drugs
#'
#' @description
#' Two data sets containing all antimicrobials and antivirals. Use [as.ab()] or one of the [`ab_*`][ab_property()] functions to retrieve values from the [antimicrobials] data set. Three identifiers are included in this data set: an antimicrobial ID (`ab`, primarily used in this package) as defined by WHONET/EARS-Net, an ATC code (`atc`) as defined by the WHO, and a Compound ID (`cid`) as found in PubChem. Other properties in this data set are derived from one or more of these codes. Note that some drugs have multiple ATC codes.
#'
#' **The `antibiotics` data set has been renamed to `antimicrobials`. The old name will be removed in a future version.**
#' @format
#' ### For the [antimicrobials] data set: a [tibble][tibble::tibble] with `r nrow(antimicrobials)` observations and `r ncol(antimicrobials)` variables:
#' - `ab`\cr antimicrobial ID as used in this package (such as `AMC`), using the official EARS-Net (European Antimicrobial Resistance Surveillance Network) codes where available. ***This is a unique identifier.***
#' - `cid`\cr Compound ID as found in PubChem. ***This is a unique identifier.***
#' - `name`\cr Official name as used by WHONET/EARS-Net or the WHO. ***This is a unique identifier.***
#' - `group`\cr A short and concise group name, based on WHONET and WHOCC definitions
#' - `atc`\cr ATC codes (Anatomical Therapeutic Chemical) as defined by the WHOCC, like `J01CR02` (last updated `r documentation_date(TAXONOMY_VERSION$ATC_DDD$accessed_date)`):
#' - `atc_group1`\cr Official pharmacological subgroup (3rd level ATC code) as defined by the WHOCC, like `"Macrolides, lincosamides and streptogramins"`
#' - `atc_group2`\cr Official chemical subgroup (4th level ATC code) as defined by the WHOCC, like `"Macrolides"`
#' - `abbr`\cr List of abbreviations as used in many countries, also for antimicrobial susceptibility testing (AST)
#' - `synonyms`\cr Synonyms (often trade names) of a drug, as found in PubChem based on their compound ID
#'
#' ATC properties (last updated `r documentation_date(TAXONOMY_VERSION$ATC_DDD$accessed_date)`):
#'
#' - `oral_ddd`\cr Defined Daily Dose (DDD), oral treatment, currently available for `r sum(!is.na(AMR::antimicrobials$oral_ddd))` drugs
#' - `oral_units`\cr Units of `oral_ddd`
#' - `iv_ddd`\cr Defined Daily Dose (DDD), parenteral (intravenous) treatment, currently available for `r sum(!is.na(AMR::antimicrobials$iv_ddd))` drugs
#' - `iv_units`\cr Units of `iv_ddd`
#'
#' LOINC:
#'
#' - `loinc`\cr All codes associated with the name of the antimicrobial drug from `r TAXONOMY_VERSION$LOINC$citation` Use [ab_loinc()] to retrieve them quickly, see [ab_property()].
#'
#' ### For the [antivirals] data set: a [tibble][tibble::tibble] with `r nrow(antivirals)` observations and `r ncol(antivirals)` variables:
#' - `av`\cr Antiviral ID as used in this package (such as `ACI`), using the official EARS-Net (European Antimicrobial Resistance Surveillance Network) codes where available. ***This is a unique identifier.*** Combinations are codes that contain a `+` to indicate this, such as `ATA+COBI` for atazanavir/cobicistat.
#' - `name`\cr Official name as used by WHONET/EARS-Net or the WHO. ***This is a unique identifier.***
#' - `atc`\cr ATC codes (Anatomical Therapeutic Chemical) as defined by the WHOCC, see *Details*
#' - `cid`\cr Compound ID as found in PubChem. ***This is a unique identifier.***
#' - `atc_group`\cr Official pharmacological subgroup (3rd level ATC code) as defined by the WHOCC
#' - `synonyms`\cr Synonyms (often trade names) of a drug, as found in PubChem based on their compound ID
#' - `oral_ddd`\cr Defined Daily Dose (DDD), oral treatment
#' - `oral_units`\cr Units of `oral_ddd`
#' - `iv_ddd`\cr Defined Daily Dose (DDD), parenteral treatment
#' - `iv_units`\cr Units of `iv_ddd`
#' - `loinc`\cr All codes associated with the name of the antiviral drug from `r TAXONOMY_VERSION$LOINC$citation` Use [av_loinc()] to retrieve them quickly, see [av_property()].
#' @details Properties that are based on an ATC code are only available when an ATC is available. These properties are: `atc_group1`, `atc_group2`, `oral_ddd`, `oral_units`, `iv_ddd` and `iv_units`. Do note that ATC codes are not unique. For example, J01CR02 is officially the ATC code for "amoxicillin and beta-lactamase inhibitor". Consequently, these two items from the [antimicrobials] data set both return `"J01CR02"`:
#'
#' ```r
#' ab_atc("amoxicillin/clavulanic acid")
#' ab_atc("amoxicillin/sulbactam")
#' ```
#'
#' Synonyms (i.e. trade names) were derived from the PubChem Compound ID (column `cid`) and are consequently only available where a CID is available.
#' @inheritSection AMR Download Our Reference Data
#' @source
#'
#' * `r TAXONOMY_VERSION$ATC_DDD$citation` Accessed from <`r TAXONOMY_VERSION$ATC_DDD$url`> on `r documentation_date(TAXONOMY_VERSION$ATC_DDD$accessed_date)`.
#'
#' * `r TAXONOMY_VERSION$LOINC$citation` Accessed from <`r TAXONOMY_VERSION$LOINC$url`> on `r documentation_date(TAXONOMY_VERSION$LOINC$accessed_date)`.
#'
#' * European Commission Public Health PHARMACEUTICALS - COMMUNITY REGISTER: <https://ec.europa.eu/health/documents/community-register/html/reg_hum_atc.htm>
#' @inheritSection WHOCC WHOCC
#' @seealso [microorganisms], [intrinsic_resistant]
#' @examples
#' antimicrobials
#' antivirals
"antimicrobials"

#' @rdname antimicrobials
"antibiotics"

#' @rdname antimicrobials
"antivirals"

#' Data Set with `r format(nrow(microorganisms), big.mark = " ")` Taxonomic Records of Microorganisms
#'
#' @description
#' A data set containing the full microbial taxonomy (**last updated: `r documentation_date(max(TAXONOMY_VERSION$GBIF$accessed_date, TAXONOMY_VERSION$LPSN$accessed_date, TAXONOMY_VERSION$MycoBank$accessed_date))`**) of `r nr2char(length(unique(microorganisms$kingdom[!microorganisms$kingdom %like% "unknown"])))` kingdoms. This data set is the backbone of this `AMR` package. MO codes can be looked up using [as.mo()] and microorganism properties can be looked up using any of the [`mo_*`][mo_property()] functions.
#'
#' This data set is carefully crafted, yet made 100% reproducible from public and authoritative taxonomic sources (using [this script](https://github.com/msberends/AMR/blob/main/data-raw/_reproduction_scripts/reproduction_of_microorganisms.R)), namely: *`r TAXONOMY_VERSION$LPSN$name`* for bacteria, *`r TAXONOMY_VERSION$MycoBank$name`* for fungi, and *`r TAXONOMY_VERSION$GBIF$name`* for all others taxons.
#' @format A [tibble][tibble::tibble] with `r format(nrow(microorganisms), big.mark = " ")` observations and `r ncol(microorganisms)` variables:
#' - `mo`\cr ID of microorganism as used by this package. ***This is a unique identifier.***
#' - `fullname`\cr Full name, like `"Escherichia coli"`. For the taxonomic ranks genus, species and subspecies, this is the 'pasted' text of genus, species, and subspecies. For all taxonomic ranks higher than genus, this is the name of the taxon. ***This is a unique identifier.***
#' - `status` \cr Status of the taxon, either `r vector_or(microorganisms$status)`
#' - `kingdom`, `phylum`, `class`, `order`, `family`, `genus`, `species`, `subspecies`\cr Taxonomic rank of the microorganism. Note that for fungi, *phylum* is equal to their taxonomic *division*. Also, for fungi, *subkingdom* and *subdivision* were left out since they do not occur in the bacterial taxonomy.
#' - `rank`\cr Text of the taxonomic rank of the microorganism, such as `"species"` or `"genus"`
#' - `ref`\cr Author(s) and year of related scientific publication. This contains only the *first surname* and year of the *latest* authors, e.g. "Wallis *et al.* 2006 *emend.* Smith and Jones 2018" becomes "Smith *et al.*, 2018". This field is directly retrieved from the source specified in the column `source`. Moreover, accents were removed to comply with CRAN that only allows ASCII characters.
#' - `oxygen_tolerance` \cr Oxygen tolerance, either `r vector_or(microorganisms$oxygen_tolerance)`. These data were retrieved from BacDive (see *Source*). Items that contain "likely" are missing from BacDive and were extrapolated from other species within the same genus to guess the oxygen tolerance. Currently `r round(length(microorganisms$oxygen_tolerance[which(!is.na(microorganisms$oxygen_tolerance))]) / nrow(microorganisms[which(microorganisms$kingdom == "Bacteria"), ]) * 100, 1)`% of all `r format_included_data_number(nrow(microorganisms[which(microorganisms$kingdom == "Bacteria"), ]))` bacteria in the data set contain an oxygen tolerance.
#' - `source`\cr Either `r vector_or(microorganisms$source)` (see *Source*)
#' - `lpsn`\cr Identifier ('Record number') of `r TAXONOMY_VERSION$LPSN$name`. This will be the first/highest LPSN identifier to keep one identifier per row. For example, *Acetobacter ascendens* has LPSN Record number 7864 and 11011. Only the first is available in the `microorganisms` data set. ***This is a unique identifier***, though available for only `r format_included_data_number(sum(!is.na(microorganisms$lpsn)))` records.
#' - `lpsn_parent`\cr LPSN identifier of the parent taxon
#' - `lpsn_renamed_to`\cr LPSN identifier of the currently valid taxon
#' - `mycobank`\cr Identifier ('MycoBank #') of `r TAXONOMY_VERSION$MycoBank$name`. ***This is a unique identifier***, though available for only `r format_included_data_number(sum(!is.na(microorganisms$mycobank)))` records.
#' - `mycobank_parent`\cr MycoBank identifier of the parent taxon
#' - `mycobank_renamed_to`\cr MycoBank identifier of the currently valid taxon
#' - `gbif`\cr Identifier ('taxonID') of `r TAXONOMY_VERSION$GBIF$name`. ***This is a unique identifier***, though available for only `r format_included_data_number(sum(!is.na(microorganisms$gbif)))` records.
#' - `gbif_parent`\cr GBIF identifier of the parent taxon
#' - `gbif_renamed_to`\cr GBIF identifier of the currently valid taxon
#' - `prevalence`\cr Prevalence of the microorganism based on Bartlett *et al.* (2022, \doi{10.1099/mic.0.001269}), see [mo_matching_score()] for the full explanation
#' - `snomed`\cr Systematized Nomenclature of Medicine (SNOMED) code of the microorganism, version of `r documentation_date(TAXONOMY_VERSION$SNOMED$accessed_date)` (see *Source*). Use [mo_snomed()] to retrieve it quickly, see [mo_property()].
#' @details
#' Please note that entries are only based on LPSN, MycoBank, and GBIF (see below). Since these sources incorporate entries based on (recent) publications in the International Journal of Systematic and Evolutionary Microbiology (IJSEM), it can happen that the year of publication is sometimes later than one might expect.
#'
#' For example, *Staphylococcus pettenkoferi* was described for the first time in Diagnostic Microbiology and Infectious Disease in 2002 (\doi{10.1016/s0732-8893(02)00399-1}), but it was not until 2007 that a publication in IJSEM followed (\doi{10.1099/ijs.0.64381-0}). Consequently, the `AMR` package returns 2007 for `mo_year("S. pettenkoferi")`.
#'
#' @section Included Taxa:
#' Included taxonomic data from [LPSN](`r TAXONOMY_VERSION$LPSN$url`), [MycoBank](`r TAXONOMY_VERSION$MycoBank$url`), and [GBIF](`r TAXONOMY_VERSION$GBIF$url`) are:
#' - All `r format_included_data_number(microorganisms[which(microorganisms$kingdom %in% c("Archeae", "Bacteria")), , drop = FALSE])` (sub)species from the kingdoms of Archaea and Bacteria
#' - `r format_included_data_number(microorganisms[which(microorganisms$kingdom == "Fungi"), , drop = FALSE])` species from the kingdom of Fungi. The kingdom of Fungi is a very large taxon with almost 300,000 different (sub)species, of which most are not microbial (but rather macroscopic, like mushrooms). Because of this, not all fungi fit the scope of this package. Only relevant fungi are covered (such as all species of *Aspergillus*, *Candida*, *Cryptococcus*, *Histoplasma*, *Pneumocystis*, *Saccharomyces* and *Trichophyton*).
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
#' - `r format_included_data_number(length(which(microorganisms$rank == "species group")))` species groups (such as the beta-haemolytic *Streptococcus* groups A to K, coagulase-negative *Staphylococcus* (CoNS), *Mycobacterium tuberculosis* complex, etc.), of which the group compositions are stored in the [microorganisms.groups] data set
#' - 1 entry of *Blastocystis* (*B.  hominis*), although it officially does not exist (Noel *et al.* 2005, PMID 15634993)
#' - 1 entry of *Moraxella* (*M. catarrhalis*), which was formally named *Branhamella catarrhalis* (Catlin, 1970) though this change was never accepted within the field of clinical microbiology
#' - 8 other 'undefined' entries (unknown, unknown Gram-negatives, unknown Gram-positives, unknown yeast, unknown fungus, and unknown anaerobic Gram-pos/Gram-neg bacteria)
#'
#' The syntax used to transform the original data to a cleansed \R format, can be [found here](https://github.com/msberends/AMR/blob/main/data-raw/_reproduction_scripts/reproduction_of_microorganisms.R).
#' @inheritSection AMR Download Our Reference Data
#' @source
#' Taxonomic entries were imported in this order of importance:
#' 1. `r TAXONOMY_VERSION$LPSN$name`:\cr\cr
#'    `r TAXONOMY_VERSION$LPSN$citation` Accessed from <`r TAXONOMY_VERSION$LPSN$url`> on `r documentation_date(TAXONOMY_VERSION$LPSN$accessed_date)`.
#'
#' 2. `r TAXONOMY_VERSION$MycoBank$name`:\cr\cr
#'    `r TAXONOMY_VERSION$MycoBank$citation` Accessed from <`r TAXONOMY_VERSION$MycoBank$url`> on `r documentation_date(TAXONOMY_VERSION$MycoBank$accessed_date)`.
#'
#' 3. `r TAXONOMY_VERSION$GBIF$name`:\cr\cr
#'    `r TAXONOMY_VERSION$GBIF$citation` Accessed from <`r TAXONOMY_VERSION$GBIF$url`> on `r documentation_date(TAXONOMY_VERSION$GBIF$accessed_date)`.
#'
#' Furthermore, these sources were used for additional details:
#'
#' * `r TAXONOMY_VERSION$BacDive$name`:\cr\cr
#'   `r TAXONOMY_VERSION$BacDive$citation` Accessed from <`r TAXONOMY_VERSION$BacDive$url`> on `r documentation_date(TAXONOMY_VERSION$BacDive$accessed_date)`.
#'
#' * `r TAXONOMY_VERSION$SNOMED$name`:\cr\cr
#'   `r TAXONOMY_VERSION$SNOMED$citation` Accessed from <`r TAXONOMY_VERSION$SNOMED$url`> on `r documentation_date(TAXONOMY_VERSION$SNOMED$accessed_date)`.
#'
#' * Grimont *et al.* (2007). Antigenic Formulae of the Salmonella Serovars, 9th Edition. WHO Collaborating Centre for Reference and Research on *Salmonella* (WHOCC-SALM).
#'
#' * Bartlett *et al.* (2022). **A comprehensive list of bacterial pathogens infecting humans** *Microbiology* 168:001269; \doi{10.1099/mic.0.001269}
#' @seealso [as.mo()], [mo_property()], [microorganisms.groups], [microorganisms.codes], [intrinsic_resistant]
#' @examples
#' microorganisms
"microorganisms"

#' Data Set with `r format(nrow(microorganisms.codes), big.mark = " ")` Common Microorganism Codes
#'
#' A data set containing commonly used codes for microorganisms, from laboratory systems and [WHONET](https://whonet.org). Define your own with [set_mo_source()]. They will all be searched when using [as.mo()] and consequently all the [`mo_*`][mo_property()] functions.
#' @format A [tibble][tibble::tibble] with `r format(nrow(microorganisms.codes), big.mark = " ")` observations and `r ncol(microorganisms.codes)` variables:
#' - `code`\cr Commonly used code of a microorganism. ***This is a unique identifier.***
#' - `mo`\cr ID of the microorganism in the [microorganisms] data set
#' @inheritSection AMR Download Our Reference Data
#' @seealso [as.mo()] [microorganisms]
#' @examples
#' microorganisms.codes
#'
#' # 'ECO' or 'eco' is the WHONET code for E. coli:
#' microorganisms.codes[microorganisms.codes$code == "ECO", ]
#'
#' # and therefore, 'eco' will be understood as E. coli in this package:
#' mo_info("eco")
#'
#' # works for all AMR functions:
#' mo_is_intrinsic_resistant("eco", ab = "vancomycin")
"microorganisms.codes"

#' Data Set with `r format(nrow(microorganisms.groups), big.mark = " ")` Microorganisms In Species Groups
#'
#' A data set containing species groups and microbiological complexes, which are used in [the clinical breakpoints table][clinical_breakpoints].
#' @format A [tibble][tibble::tibble] with `r format(nrow(microorganisms.groups), big.mark = " ")` observations and `r ncol(microorganisms.groups)` variables:
#' - `mo_group`\cr ID of the species group / microbiological complex
#' - `mo`\cr ID of the microorganism belonging in the species group / microbiological complex
#' - `mo_group_name`\cr Name of the species group / microbiological complex, as retrieved with [mo_name()]
#' - `mo_name`\cr Name of the microorganism belonging in the species group / microbiological complex, as retrieved with [mo_name()]
#' @inheritSection AMR Download Our Reference Data
#' @seealso [as.mo()] [microorganisms]
#' @examples
#' microorganisms.groups
#'
#' # these are all species in the Bacteroides fragilis group, as per WHONET:
#' microorganisms.groups[microorganisms.groups$mo_group == "B_BCTRD_FRGL-C", ]
"microorganisms.groups"

#' Data Set with `r format(nrow(example_isolates), big.mark = " ")` Example Isolates
#'
#' A data set containing `r format(nrow(example_isolates), big.mark = " ")` microbial isolates with their full antibiograms. This data set contains randomised fictitious data, but reflects reality and can be used to practise AMR data analysis. For examples, please read [the tutorial on our website](https://amr-for-r.org/articles/AMR.html).
#' @format A [tibble][tibble::tibble] with `r format(nrow(example_isolates), big.mark = " ")` observations and `r ncol(example_isolates)` variables:
#' - `date`\cr Date of receipt at the laboratory
#' - `patient`\cr ID of the patient
#' - `age`\cr Age of the patient
#' - `gender`\cr Gender of the patient, either `r vector_or(example_isolates$gender)`
#' - `ward`\cr Ward type where the patient was admitted, either `r vector_or(example_isolates$ward)`
#' - `mo`\cr ID of microorganism created with [as.mo()], see also the [microorganisms] data set
#' - `PEN:RIF`\cr `r sum(vapply(FUN.VALUE = logical(1), example_isolates, is.sir))` different antimicrobials with class [`sir`] (see [as.sir()]); these column names occur in the [antimicrobials] data set and can be translated with [set_ab_names()] or [ab_name()]
#' @inheritSection AMR Download Our Reference Data
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
#' - `AMX:GEN`\cr 4 different antimicrobials that have to be transformed with [as.sir()]
#' @inheritSection AMR Download Our Reference Data
#' @examples
#' example_isolates_unclean
"example_isolates_unclean"

#' Data Set with `r format(nrow(WHONET), big.mark = " ")` Isolates - WHONET Example
#'
#' This example data set has the exact same structure as an export file from WHONET. Such files can be used with this package, as this example data set shows. The antimicrobial results are from our [example_isolates] data set. All patient names were created using online surname generators and are only in place for practice purposes.
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
#' - `AMP_ND10:CIP_EE`\cr `r sum(vapply(FUN.VALUE = logical(1), WHONET, is.sir))` different antimicrobials. You can lookup the abbreviations in the [antimicrobials] data set, or use e.g. [`ab_name("AMP")`][ab_name()] to get the official name immediately. Before analysis, you should transform this to a valid antimicrobial class, using [as.sir()].
#' @inheritSection AMR Download Our Reference Data
#' @examples
#' WHONET
"WHONET"

#' Data Set with Clinical Breakpoints for SIR Interpretation
#'
#' @description Data set containing clinical breakpoints to interpret MIC and disk diffusion to SIR values, according to international guidelines. This data set contains breakpoints for humans, `r length(unique(clinical_breakpoints$host[!clinical_breakpoints$host %in% clinical_breakpoints$type]))` different animal groups, and ECOFFs.
#'
#' These breakpoints are currently implemented:
#' - For **clinical microbiology**: EUCAST `r min(as.integer(gsub("[^0-9]", "", subset(AMR::clinical_breakpoints, guideline %like% "EUCAST" & type == "human")$guideline)))`-`r max(as.integer(gsub("[^0-9]", "", subset(AMR::clinical_breakpoints, guideline %like% "EUCAST" & type == "human")$guideline)))` and CLSI `r min(as.integer(gsub("[^0-9]", "", subset(AMR::clinical_breakpoints, guideline %like% "CLSI" & type == "human")$guideline)))`-`r max(as.integer(gsub("[^0-9]", "", subset(AMR::clinical_breakpoints, guideline %like% "CLSI" & type == "human")$guideline)))`;
#' - For **veterinary microbiology**: EUCAST `r min(as.integer(gsub("[^0-9]", "", subset(AMR::clinical_breakpoints, guideline %like% "EUCAST" & type == "animal")$guideline)))`-`r max(as.integer(gsub("[^0-9]", "", subset(AMR::clinical_breakpoints, guideline %like% "EUCAST" & type == "animal")$guideline)))` and CLSI `r min(as.integer(gsub("[^0-9]", "", subset(AMR::clinical_breakpoints, guideline %like% "CLSI" & type == "animal")$guideline)))`-`r max(as.integer(gsub("[^0-9]", "", subset(AMR::clinical_breakpoints, guideline %like% "CLSI" & type == "animal")$guideline)))`;
#' - For **ECOFFs** (Epidemiological Cut-off Values): EUCAST `r min(as.integer(gsub("[^0-9]", "", subset(AMR::clinical_breakpoints, guideline %like% "EUCAST" & type == "ECOFF")$guideline)))`-`r max(as.integer(gsub("[^0-9]", "", subset(AMR::clinical_breakpoints, guideline %like% "EUCAST" & type == "ECOFF")$guideline)))` and CLSI `r min(as.integer(gsub("[^0-9]", "", subset(AMR::clinical_breakpoints, guideline %like% "CLSI" & type == "ECOFF")$guideline)))`-`r max(as.integer(gsub("[^0-9]", "", subset(AMR::clinical_breakpoints, guideline %like% "CLSI" & type == "ECOFF")$guideline)))`.
#'
#' Use [as.sir()] to transform MICs or disks measurements to SIR values.
#' @format A [tibble][tibble::tibble] with `r format(nrow(clinical_breakpoints), big.mark = " ")` observations and `r ncol(clinical_breakpoints)` variables:
#' - `guideline`\cr Name of the guideline
#' - `type`\cr Breakpoint type, either `r vector_or(clinical_breakpoints$type)`
#' - `host`\cr Host of infectious agent. This is mostly useful for veterinary breakpoints and is either `r vector_or(clinical_breakpoints$host)`
#' - `method`\cr Testing method, either `r vector_or(clinical_breakpoints$method)`
#' - `site`\cr Body site for which the breakpoint must be applied, e.g. "Oral" or "Respiratory"
#' - `mo`\cr Microbial ID, see [as.mo()]
#' - `rank_index`\cr Taxonomic rank index of `mo` from 1 (subspecies/infraspecies) to 5 (unknown microorganism)
#' - `ab`\cr Antimicrobial code as used by this package, EARS-Net and WHONET, see [as.ab()]
#' - `ref_tbl`\cr Info about where the guideline rule can be found
#' - `disk_dose`\cr Dose of the used disk diffusion method
#' - `breakpoint_S`\cr Lowest MIC value or highest number of millimetres that leads to "S"
#' - `breakpoint_R`\cr Highest MIC value or lowest number of millimetres that leads to "R", can be `NA`
#' - `uti`\cr A [logical] value (`TRUE`/`FALSE`) to indicate whether the rule applies to a urinary tract infection (UTI)
#' - `is_SDD`\cr A [logical] value (`TRUE`/`FALSE`) to indicate whether the intermediate range between "S" and "R" should be interpreted as "SDD", instead of "I". This currently applies to `r sum(clinical_breakpoints$is_SDD)` breakpoints.
#' @details
#' ### Different Types of Breakpoints
#' Supported types of breakpoints are `r vector_and(clinical_breakpoints$type, quote = FALSE)`. ECOFF (Epidemiological cut-off) values are used in antimicrobial susceptibility testing to differentiate between wild-type and non-wild-type strains of bacteria or fungi.
#'
#' The default is `"human"`, which can also be set with the package option [`AMR_breakpoint_type`][AMR-options]. Use [`as.sir(..., breakpoint_type = ...)`][as.sir()] to interpret raw data using a specific breakpoint type, e.g. `as.sir(..., breakpoint_type = "ECOFF")` to use ECOFFs.
#'
#' ### Imported From WHONET
#' Clinical breakpoints in this package were validated through and imported from [WHONET](https://whonet.org), a free desktop Windows application developed and supported by the WHO Collaborating Centre for Surveillance of Antimicrobial Resistance. More can be read on [their website](https://whonet.org). The developers of WHONET and this `AMR` package have been in contact about sharing their work. We highly appreciate their great development on the WHONET software.
#'
#' Our import and reproduction script can be found here: <https://github.com/msberends/AMR/blob/main/data-raw/_reproduction_scripts/reproduction_of_clinical_breakpoints.R>.
#'
#' ### Response From CLSI and EUCAST
#' The CEO of CLSI and the chairman of EUCAST have endorsed the work and public use of this `AMR` package (and consequently the use of their breakpoints) in June 2023, when future development of distributing clinical breakpoints was discussed in a meeting between CLSI, EUCAST, WHO, developers of WHONET software, and developers of this `AMR` package.
#'
#' ### Download Note
#' This `AMR` package (and the WHONET software as well) contains rather complex internal methods to apply the guidelines. For example, some breakpoints must be applied on certain species groups (which are in case of this package available through the [microorganisms.groups] data set). It is important that this is considered when implementing the breakpoints for own use.
#' @inheritSection AMR Download Our Reference Data
#' @seealso [intrinsic_resistant]
#' @examples
#' clinical_breakpoints
"clinical_breakpoints"

#' Data Set Denoting Bacterial Intrinsic Resistance
#'
#' Data set containing `r EUCAST_VERSION_EXPECTED_PHENOTYPES[[1]]$title` of *all* bug-drug combinations between the [microorganisms] and [antimicrobials] data sets.
#' @format A [tibble][tibble::tibble] with `r format(nrow(intrinsic_resistant), big.mark = " ")` observations and `r ncol(intrinsic_resistant)` variables:
#' - `mo`\cr Microorganism ID which occurs in [`microorganisms$mo`][microorganisms]. Names can be retrieved using [mo_name()].
#' - `ab`\cr Antimicrobial ID which occurs in [`antimicrobials$ab`][antimicrobials]. Names can be retrieved using [ab_name()].
#' @details
#' This data set is currently based on `r format_eucast_version_nr(names(EUCAST_VERSION_EXPECTED_PHENOTYPES[1]))`.
#'
#' This data set is internally used by:
#' * [not_intrinsic_resistant()] (an [antimicrobial selector][antimicrobial_selectors])
#' * [mo_is_intrinsic_resistant()]
#' @inheritSection AMR Download Our Reference Data
#' @examples
#' intrinsic_resistant
"intrinsic_resistant"

#' Data Set with Treatment Dosages as Defined by EUCAST
#'
#' EUCAST breakpoints used in this package are based on the dosages in this data set. They can be retrieved with [eucast_dosage()].
#' @format A [tibble][tibble::tibble] with `r format(nrow(dosage), big.mark = " ")` observations and `r ncol(dosage)` variables:
#' - `ab`\cr Antimicrobial ID as used in this package (such as `AMC`), using the official EARS-Net (European Antimicrobial Resistance Surveillance Network) codes where available
#' - `name`\cr Official name of the antimicrobial drug as used by WHONET/EARS-Net or the WHO
#' - `type`\cr Type of the dosage, either `r vector_or(dosage$type)`
#' - `dose`\cr Dose, such as "2 g" or "25 mg/kg"
#' - `dose_times`\cr Number of times a dose must be administered
#' - `administration`\cr Route of administration, either `r vector_or(dosage$administration)`
#' - `notes`\cr Additional dosage notes
#' - `original_txt`\cr Original text in the PDF file of EUCAST
#' - `eucast_version`\cr Version number of the EUCAST Clinical Breakpoints guideline to which these dosages apply, either `r vector_or(dosage$eucast_version, quotes = FALSE, sort = TRUE, reverse = TRUE)`
#' @inheritSection AMR Download Our Reference Data
#' @examples
#' dosage
"dosage"

#' Data Set with `r format(nrow(esbl_isolates), big.mark = " ")` ESBL Isolates
#'
#' A data set containing `r format(nrow(esbl_isolates), big.mark = " ")` microbial isolates with MIC values of common antibiotics and a binary `esbl` column for extended-spectrum beta-lactamase (ESBL) production. This data set contains randomised fictitious data but reflects reality and can be used to practise AMR-related machine learning, e.g., classification modelling with [tidymodels](https://amr-for-r.org/articles/AMR_with_tidymodels.html).
#' @format A [tibble][tibble::tibble] with `r format(nrow(esbl_isolates), big.mark = " ")` observations and `r ncol(esbl_isolates)` variables:
#' - `esbl`\cr Logical indicator if the isolate is ESBL-producing
#' - `genus`\cr Genus of the microorganism
#' - `AMC:COL`\cr MIC values for 17 antimicrobial agents, transformed to class [`mic`] (see [as.mic()])
#' @details See our [tidymodels integration][amr-tidymodels] for an example using this data set.
#' @examples
#' esbl_isolates
"esbl_isolates"
