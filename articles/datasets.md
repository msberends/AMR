# Download data sets for download / own use

All reference data (about microorganisms, antimicrobials, SIR
interpretation, EUCAST rules, etc.) in this `AMR` package are reliable,
up-to-date and freely available. We continually export our data sets to
formats for use in R, MS Excel, Apache Feather, Apache Parquet, SPSS,
and Stata. We also provide tab-separated text files that are
machine-readable and suitable for input in any software program, such as
laboratory information systems.

> If you are working in Python, be sure to use our [AMR for
> Python](https://amr-for-r.org/articles/AMR_for_Python.html) package.
> It allows all relevant AMR data sets to be natively available in
> Python.

## `microorganisms`: Full Microbial Taxonomy

A data set with 78 679 rows and 26 columns, containing the following
column names:  
*mo*, *fullname*, *status*, *kingdom*, *phylum*, *class*, *order*,
*family*, *genus*, *species*, *subspecies*, *rank*, *ref*,
*oxygen_tolerance*, *source*, *lpsn*, *lpsn_parent*, *lpsn_renamed_to*,
*mycobank*, *mycobank_parent*, *mycobank_renamed_to*, *gbif*,
*gbif_parent*, *gbif_renamed_to*, *prevalence*, and *snomed*.

This data set is in R available as `microorganisms`, after you load the
`AMR` package.

It was last updated on 18 September 2025 12:58:34 UTC. Find more info
about the contents, (scientific) source, and structure of this [data set
here](https://amr-for-r.org/reference/microorganisms.html).

**Direct download links:**

- Download as [original R Data Structure (RDS)
  file](https://github.com/msberends/AMR/raw/main/data-raw/datasets/microorganisms.rds)
  (1.8 MB)  
- Download as [tab-separated text
  file](https://github.com/msberends/AMR/raw/main/data-raw/datasets/microorganisms.txt)
  (17.7 MB)  
- Download as [Microsoft Excel
  workbook](https://github.com/msberends/AMR/raw/main/data-raw/datasets/microorganisms.xlsx)
  (8.8 MB)  
- Download as [Apache Feather
  file](https://github.com/msberends/AMR/raw/main/data-raw/datasets/microorganisms.feather)
  (8.4 MB)  
- Download as [Apache Parquet
  file](https://github.com/msberends/AMR/raw/main/data-raw/datasets/microorganisms.parquet)
  (3.8 MB)  
- Download as [IBM SPSS Statistics data
  file](https://github.com/msberends/AMR/raw/main/data-raw/datasets/microorganisms.sav)
  (28.4 MB)  
- Download as [Stata DTA
  file](https://github.com/msberends/AMR/raw/main/data-raw/datasets/microorganisms.dta)
  (89.5 MB)

**NOTE: The exported files for SPSS and Stata contain only the first 50
SNOMED codes per record, as their file size would otherwise exceed 100
MB; the file size limit of GitHub.** Their file structures and
compression techniques are very inefficient. Advice? Use R instead. It’s
free and much better in many ways.

The tab-separated text file and Microsoft Excel workbook both contain
all SNOMED codes as comma separated values.

**Example content**

Included (sub)species per taxonomic kingdom:

|      Kingdom      | Number of (sub)species |
|:-----------------:|:----------------------:|
| (unknown kingdom) |           1            |
|     Animalia      |         1 628          |
|      Archaea      |         1 419          |
|     Bacteria      |         39 249         |
|     Chromista     |          178           |
|       Fungi       |         28 137         |

First 6 rows when filtering on genus *Escherichia*:

|        mo         |          fullname          |  status  | kingdom  |     phylum     |        class        |      order       |       family       |    genus    |    species     | subspecies |    rank    |           ref           |      oxygen_tolerance       | source |  lpsn  | lpsn_parent | lpsn_renamed_to | mycobank | mycobank_parent | mycobank_renamed_to |   gbif   | gbif_parent | gbif_renamed_to | prevalence |                  snomed                   |
|:-----------------:|:--------------------------:|:--------:|:--------:|:--------------:|:-------------------:|:----------------:|:------------------:|:-----------:|:--------------:|:----------:|:----------:|:-----------------------:|:---------------------------:|:------:|:------:|:-----------:|:---------------:|:--------:|:---------------:|:-------------------:|:--------:|:-----------:|:---------------:|:----------:|:-----------------------------------------:|
|      B_ESCHR      |        Escherichia         | accepted | Bacteria | Pseudomonadota | Gammaproteobacteria | Enterobacterales | Enterobacteriaceae | Escherichia |                |            |   genus    | Castellani et al., 1919 |    facultative anaerobe     |  LPSN  | 515602 |     482     |                 |          |                 |                     |          |  11158430   |                 |     1      |    407310004, 407251000, 407281008, …     |
|   B_ESCHR_ADCR    | Escherichia adecarboxylata | synonym  | Bacteria | Pseudomonadota | Gammaproteobacteria | Enterobacterales | Enterobacteriaceae | Escherichia | adecarboxylata |            |  species   |      Leclerc, 1962      | likely facultative anaerobe |  LPSN  | 776052 |   515602    |     777447      |          |                 |                     |          |             |                 |     1      |                                           |
|   B_ESCHR_ALBR    |    Escherichia albertii    | accepted | Bacteria | Pseudomonadota | Gammaproteobacteria | Enterobacterales | Enterobacteriaceae | Escherichia |    albertii    |            |  species   |    Huys et al., 2003    |    facultative anaerobe     |  LPSN  | 776053 |   515602    |                 |          |                 |                     | 5427575  |             |                 |     1      |                 419388003                 |
|   B_ESCHR_BLTT    |    Escherichia blattae     | synonym  | Bacteria | Pseudomonadota | Gammaproteobacteria | Enterobacterales | Enterobacteriaceae | Escherichia |    blattae     |            |  species   |  Burgess et al., 1973   | likely facultative anaerobe |  LPSN  | 776056 |   515602    |     788468      |          |                 |                     |          |             |                 |     1      |                                           |
|   B_ESCHR_COLI    |      Escherichia coli      | accepted | Bacteria | Pseudomonadota | Gammaproteobacteria | Enterobacterales | Enterobacteriaceae | Escherichia |      coli      |            |  species   | Castellani et al., 1919 |    facultative anaerobe     |  LPSN  | 776057 |   515602    |                 |          |                 |                     | 11286021 |             |                 |     1      | 1095001000112106, 715307006, 737528008, … |
| B_ESCHR_COLI_COLI |   Escherichia coli coli    | accepted | Bacteria | Pseudomonadota | Gammaproteobacteria | Enterobacterales | Enterobacteriaceae | Escherichia |      coli      |    coli    | subspecies |                         |                             |  GBIF  |        |   776057    |                 |          |                 |                     | 12233256 |  11286021   |                 |     1      |                                           |

------------------------------------------------------------------------

## `antimicrobials`: Antibiotic and Antifungal Drugs

A data set with 498 rows and 14 columns, containing the following column
names:  
*ab*, *cid*, *name*, *group*, *atc*, *atc_group1*, *atc_group2*,
*abbreviations*, *synonyms*, *oral_ddd*, *oral_units*, *iv_ddd*,
*iv_units*, and *loinc*.

This data set is in R available as `antimicrobials`, after you load the
`AMR` package.

It was last updated on 16 January 2026 09:57:03 UTC. Find more info
about the contents, (scientific) source, and structure of this [data set
here](https://amr-for-r.org/reference/antimicrobials.html).

**Direct download links:**

- Download as [original R Data Structure (RDS)
  file](https://github.com/msberends/AMR/raw/main/data-raw/datasets/antimicrobials.rds)
  (46 kB)  
- Download as [tab-separated text
  file](https://github.com/msberends/AMR/raw/main/data-raw/datasets/antimicrobials.txt)
  (0.1 MB)  
- Download as [Microsoft Excel
  workbook](https://github.com/msberends/AMR/raw/main/data-raw/datasets/antimicrobials.xlsx)
  (78 kB)  
- Download as [Apache Feather
  file](https://github.com/msberends/AMR/raw/main/data-raw/datasets/antimicrobials.feather)
  (0.1 MB)  
- Download as [Apache Parquet
  file](https://github.com/msberends/AMR/raw/main/data-raw/datasets/antimicrobials.parquet)
  (0.1 MB)  
- Download as [IBM SPSS Statistics data
  file](https://github.com/msberends/AMR/raw/main/data-raw/datasets/antimicrobials.sav)
  (0.4 MB)  
- Download as [Stata DTA
  file](https://github.com/msberends/AMR/raw/main/data-raw/datasets/antimicrobials.dta)
  (10 kB)

The tab-separated text, Microsoft Excel, SPSS, and Stata files all
contain the ATC codes, common abbreviations, trade names and LOINC codes
as comma separated values.

**Example content**

| ab  |   cid    |            name             |                     group                      |              atc               |                 atc_group1                  |                          atc_group2                          |    abbreviations    |                        synonyms                         | oral_ddd | oral_units | iv_ddd | iv_units |             loinc              |
|:---:|:--------:|:---------------------------:|:----------------------------------------------:|:------------------------------:|:-------------------------------------------:|:------------------------------------------------------------:|:-------------------:|:-------------------------------------------------------:|:--------:|:----------:|:------:|:--------:|:------------------------------:|
| AMK |  37768   |          Amikacin           |                Aminoglycosides                 | D06AX12, J01GB06, QD06AX12, …  |        Aminoglycoside antibacterials        |                    Other aminoglycosides                     |  ak, ami, amik, …   |          amikacillin, amikacina, amikacine, …           |          |            |  1.0   |    g     |    101493-5, 11-7, 12-5, …     |
| AMX |  33613   |         Amoxicillin         |  Aminopenicillins, Penicillins, Beta-lactams   |  J01CA04, QG51AA03, QJ01CA04   |   Beta-lactam antibacterials, penicillins   |              Penicillins with extended spectrum              | ac, amox, amoxic, … |             acuotricina, alfamox, alfida, …             |   1.5    |     g      |  3.0   |    g     |    101498-4, 15-8, 16-6, …     |
| AMC | 23665637 | Amoxicillin/clavulanic acid | Aminopenicillins, Penicillins, Beta-lactams, … |       J01CR02, QJ01CR02        |   Beta-lactam antibacterials, penicillins   | Combinations of penicillins, incl. beta-lactamase inhibitors |  a/c, amcl, aml, …  |               amocla, amoclan, amoclav, …               |   1.5    |     g      |  3.0   |    g     |                                |
| AMP |   6249   |         Ampicillin          |  Aminopenicillins, Penicillins, Beta-lactams   | J01CA01, QJ01CA01, QJ51CA01, … |   Beta-lactam antibacterials, penicillins   |              Penicillins with extended spectrum              | am, amp, amp100, …  |             adobacillin, alpen, amblosin, …             |   2.0    |     g      |  6.0   |    g     | 101477-8, 101478-6, 18864-9, … |
| AZM |  447043  |        Azithromycin         |                   Macrolides                   | J01FA10, QJ01FA10, QS01AA26, … | Macrolides, lincosamides and streptogramins |                          Macrolides                          |  az, azi, azit, …   |           aritromicina, aruzilina, azasite, …           |   0.3    |     g      |  0.5   |    g     | 100043-9, 16420-2, 16421-0, …  |
| PEN |   5904   |      Benzylpenicillin       |           Penicillins, Beta-lactams            | J01CE01, QJ01CE01, QJ51CE01, … |       Combinations of antibacterials        |                Combinations of antibacterials                | bepe, pen, peni, …  | bencilpenicilina, benzopenicillin, benzylpenicilline, … |          |            |  3.6   |    g     |                                |

------------------------------------------------------------------------

## `clinical_breakpoints`: Interpretation from MIC values & disk diameters to SIR

A data set with 40 217 rows and 14 columns, containing the following
column names:  
*guideline*, *type*, *host*, *method*, *site*, *mo*, *rank_index*, *ab*,
*ref_tbl*, *disk_dose*, *breakpoint_S*, *breakpoint_R*, *uti*, and
*is_SDD*.

This data set is in R available as `clinical_breakpoints`, after you
load the `AMR` package.

It was last updated on 20 April 2025 10:55:31 UTC. Find more info about
the contents, (scientific) source, and structure of this [data set
here](https://amr-for-r.org/reference/clinical_breakpoints.html).

**Direct download links:**

- Download as [original R Data Structure (RDS)
  file](https://github.com/msberends/AMR/raw/main/data-raw/datasets/clinical_breakpoints.rds)
  (88 kB)  
- Download as [tab-separated text
  file](https://github.com/msberends/AMR/raw/main/data-raw/datasets/clinical_breakpoints.txt)
  (3.7 MB)  
- Download as [Microsoft Excel
  workbook](https://github.com/msberends/AMR/raw/main/data-raw/datasets/clinical_breakpoints.xlsx)
  (2.4 MB)  
- Download as [Apache Feather
  file](https://github.com/msberends/AMR/raw/main/data-raw/datasets/clinical_breakpoints.feather)
  (1.8 MB)  
- Download as [Apache Parquet
  file](https://github.com/msberends/AMR/raw/main/data-raw/datasets/clinical_breakpoints.parquet)
  (0.1 MB)  
- Download as [IBM SPSS Statistics data
  file](https://github.com/msberends/AMR/raw/main/data-raw/datasets/clinical_breakpoints.sav)
  (6.6 MB)  
- Download as [Stata DTA
  file](https://github.com/msberends/AMR/raw/main/data-raw/datasets/clinical_breakpoints.dta)
  (11.1 MB)

**Example content**

|  guideline  | type  | host  | method | site |      mo       |          mo_name           | rank_index | ab  |            ab_name            |     ref_tbl     |   disk_dose    | breakpoint_S | breakpoint_R |  uti  | is_SDD |
|:-----------:|:-----:|:-----:|:------:|:----:|:-------------:|:--------------------------:|:----------:|:---:|:-----------------------------:|:---------------:|:--------------:|:------------:|:------------:|:-----:|:------:|
| EUCAST 2025 | human | human |  DISK  |      | B_ACHRMB_XYLS | Achromobacter xylosoxidans |     2      | MEM |           Meropenem           | A. xylosoxidans |     10 mcg     |    26.000    |    20.000    | FALSE | FALSE  |
| EUCAST 2025 | human | human |  MIC   |      | B_ACHRMB_XYLS | Achromobacter xylosoxidans |     2      | MEM |           Meropenem           | A. xylosoxidans |                |    1.000     |    4.000     | FALSE | FALSE  |
| EUCAST 2025 | human | human |  DISK  |      | B_ACHRMB_XYLS | Achromobacter xylosoxidans |     2      | SXT | Trimethoprim/sulfamethoxazole | A. xylosoxidans | 1.25/23.75 mcg |    26.000    |    26.000    | FALSE | FALSE  |
| EUCAST 2025 | human | human |  MIC   |      | B_ACHRMB_XYLS | Achromobacter xylosoxidans |     2      | SXT | Trimethoprim/sulfamethoxazole | A. xylosoxidans |                |    0.125     |    0.125     | FALSE | FALSE  |
| EUCAST 2025 | human | human |  DISK  |      | B_ACHRMB_XYLS | Achromobacter xylosoxidans |     2      | TZP |    Piperacillin/tazobactam    | A. xylosoxidans |    30/6 mcg    |    26.000    |    26.000    | FALSE | FALSE  |
| EUCAST 2025 | human | human |  MIC   |      | B_ACHRMB_XYLS | Achromobacter xylosoxidans |     2      | TZP |    Piperacillin/tazobactam    | A. xylosoxidans |                |    4.000     |    4.000     | FALSE | FALSE  |

------------------------------------------------------------------------

## `microorganisms.groups`: Species Groups and Microbiological Complexes

A data set with 534 rows and 4 columns, containing the following column
names:  
*mo_group*, *mo*, *mo_group_name*, and *mo_name*.

This data set is in R available as `microorganisms.groups`, after you
load the `AMR` package.

It was last updated on 26 March 2025 16:19:17 UTC. Find more info about
the contents, (scientific) source, and structure of this [data set
here](https://amr-for-r.org/reference/microorganisms.groups.html).

**Direct download links:**

- Download as [original R Data Structure (RDS)
  file](https://github.com/msberends/AMR/raw/main/data-raw/datasets/microorganisms.groups.rds)
  (6 kB)  
- Download as [tab-separated text
  file](https://github.com/msberends/AMR/raw/main/data-raw/datasets/microorganisms.groups.txt)
  (50 kB)  
- Download as [Microsoft Excel
  workbook](https://github.com/msberends/AMR/raw/main/data-raw/datasets/microorganisms.groups.xlsx)
  (20 kB)  
- Download as [Apache Feather
  file](https://github.com/msberends/AMR/raw/main/data-raw/datasets/microorganisms.groups.feather)
  (19 kB)  
- Download as [Apache Parquet
  file](https://github.com/msberends/AMR/raw/main/data-raw/datasets/microorganisms.groups.parquet)
  (13 kB)  
- Download as [IBM SPSS Statistics data
  file](https://github.com/msberends/AMR/raw/main/data-raw/datasets/microorganisms.groups.sav)
  (65 kB)  
- Download as [Stata DTA
  file](https://github.com/msberends/AMR/raw/main/data-raw/datasets/microorganisms.groups.dta)
  (83 kB)

**Example content**

|    mo_group    |      mo      |          mo_group_name          |           mo_name           |
|:--------------:|:------------:|:-------------------------------:|:---------------------------:|
| B_ACNTB_BMNN-C | B_ACNTB_BMNN | Acinetobacter baumannii complex |   Acinetobacter baumannii   |
| B_ACNTB_BMNN-C | B_ACNTB_CLCC | Acinetobacter baumannii complex | Acinetobacter calcoaceticus |
| B_ACNTB_BMNN-C | B_ACNTB_LCTC | Acinetobacter baumannii complex | Acinetobacter dijkshoorniae |
| B_ACNTB_BMNN-C | B_ACNTB_NSCM | Acinetobacter baumannii complex | Acinetobacter nosocomialis  |
| B_ACNTB_BMNN-C | B_ACNTB_PITT | Acinetobacter baumannii complex |    Acinetobacter pittii     |
| B_ACNTB_BMNN-C | B_ACNTB_SFRT | Acinetobacter baumannii complex |   Acinetobacter seifertii   |

------------------------------------------------------------------------

## `intrinsic_resistant`: Intrinsic Bacterial Resistance

A data set with 271 905 rows and 2 columns, containing the following
column names:  
*mo* and *ab*.

This data set is in R available as `intrinsic_resistant`, after you load
the `AMR` package.

It was last updated on 28 March 2025 10:17:49 UTC. Find more info about
the contents, (scientific) source, and structure of this [data set
here](https://amr-for-r.org/reference/intrinsic_resistant.html).

**Direct download links:**

- Download as [original R Data Structure (RDS)
  file](https://github.com/msberends/AMR/raw/main/data-raw/datasets/intrinsic_resistant.rds)
  (0.1 MB)  
- Download as [tab-separated text
  file](https://github.com/msberends/AMR/raw/main/data-raw/datasets/intrinsic_resistant.txt)
  (10.1 MB)  
- Download as [Microsoft Excel
  workbook](https://github.com/msberends/AMR/raw/main/data-raw/datasets/intrinsic_resistant.xlsx)
  (2.9 MB)  
- Download as [Apache Feather
  file](https://github.com/msberends/AMR/raw/main/data-raw/datasets/intrinsic_resistant.feather)
  (2.3 MB)  
- Download as [Apache Parquet
  file](https://github.com/msberends/AMR/raw/main/data-raw/datasets/intrinsic_resistant.parquet)
  (0.3 MB)  
- Download as [IBM SPSS Statistics data
  file](https://github.com/msberends/AMR/raw/main/data-raw/datasets/intrinsic_resistant.sav)
  (14.8 MB)  
- Download as [Stata DTA
  file](https://github.com/msberends/AMR/raw/main/data-raw/datasets/intrinsic_resistant.dta)
  (22.6 MB)

**Example content**

Example rows when filtering on *Enterobacter cloacae*:

|    microorganism     |         antibiotic          |
|:--------------------:|:---------------------------:|
| Enterobacter cloacae |      Acetylmidecamycin      |
| Enterobacter cloacae |      Acetylspiramycin       |
| Enterobacter cloacae |         Amoxicillin         |
| Enterobacter cloacae | Amoxicillin/clavulanic acid |
| Enterobacter cloacae |         Ampicillin          |
| Enterobacter cloacae |    Ampicillin/sulbactam     |
| Enterobacter cloacae |          Avoparcin          |
| Enterobacter cloacae |        Azithromycin         |
| Enterobacter cloacae |      Benzylpenicillin       |
| Enterobacter cloacae |          Bleomycin          |
| Enterobacter cloacae |          Cadazolid          |
| Enterobacter cloacae |         Cefadroxil          |
| Enterobacter cloacae |          Cefalexin          |
| Enterobacter cloacae |          Cefalotin          |
| Enterobacter cloacae |          Cefazolin          |
| Enterobacter cloacae |          Cefoxitin          |
| Enterobacter cloacae |       Clarithromycin        |
| Enterobacter cloacae |         Clindamycin         |
| Enterobacter cloacae |         Cycloserine         |
| Enterobacter cloacae |         Dalbavancin         |
| Enterobacter cloacae |        Dirithromycin        |
| Enterobacter cloacae |        Erythromycin         |
| Enterobacter cloacae |       Flurithromycin        |
| Enterobacter cloacae |        Fusidic acid         |
| Enterobacter cloacae |        Gamithromycin        |
| Enterobacter cloacae |          Josamycin          |
| Enterobacter cloacae |         Kitasamycin         |
| Enterobacter cloacae |         Lincomycin          |
| Enterobacter cloacae |          Linezolid          |
| Enterobacter cloacae |         Meleumycin          |
| Enterobacter cloacae |         Midecamycin         |
| Enterobacter cloacae |         Miocamycin          |
| Enterobacter cloacae |        Nafithromycin        |
| Enterobacter cloacae |        Norvancomycin        |
| Enterobacter cloacae |        Oleandomycin         |
| Enterobacter cloacae |         Oritavancin         |
| Enterobacter cloacae |         Pirlimycin          |
| Enterobacter cloacae |        Pristinamycin        |
| Enterobacter cloacae |  Quinupristin/dalfopristin  |
| Enterobacter cloacae |         Ramoplanin          |
| Enterobacter cloacae |         Rifampicin          |
| Enterobacter cloacae |         Rokitamycin         |
| Enterobacter cloacae |        Roxithromycin        |
| Enterobacter cloacae |        Solithromycin        |
| Enterobacter cloacae |         Spiramycin          |
| Enterobacter cloacae |          Tedizolid          |
| Enterobacter cloacae |         Teicoplanin         |
| Enterobacter cloacae |         Telavancin          |
| Enterobacter cloacae |        Telithromycin        |
| Enterobacter cloacae |        Thiacetazone         |
| Enterobacter cloacae |        Tildipirosin         |
| Enterobacter cloacae |         Tilmicosin          |
| Enterobacter cloacae |       Troleandomycin        |
| Enterobacter cloacae |        Tulathromycin        |
| Enterobacter cloacae |           Tylosin           |
| Enterobacter cloacae |         Tylvalosin          |
| Enterobacter cloacae |         Vancomycin          |

------------------------------------------------------------------------

## `dosage`: Dosage Guidelines from EUCAST

A data set with 759 rows and 9 columns, containing the following column
names:  
*ab*, *name*, *type*, *dose*, *dose_times*, *administration*, *notes*,
*original_txt*, and *eucast_version*.

This data set is in R available as `dosage`, after you load the `AMR`
package.

It was last updated on 20 April 2025 10:55:31 UTC. Find more info about
the contents, (scientific) source, and structure of this [data set
here](https://amr-for-r.org/reference/dosage.html).

**Direct download links:**

- Download as [original R Data Structure (RDS)
  file](https://github.com/msberends/AMR/raw/main/data-raw/datasets/dosage.rds)
  (4 kB)  
- Download as [tab-separated text
  file](https://github.com/msberends/AMR/raw/main/data-raw/datasets/dosage.txt)
  (66 kB)  
- Download as [Microsoft Excel
  workbook](https://github.com/msberends/AMR/raw/main/data-raw/datasets/dosage.xlsx)
  (37 kB)  
- Download as [Apache Feather
  file](https://github.com/msberends/AMR/raw/main/data-raw/datasets/dosage.feather)
  (28 kB)  
- Download as [Apache Parquet
  file](https://github.com/msberends/AMR/raw/main/data-raw/datasets/dosage.parquet)
  (9 kB)  
- Download as [IBM SPSS Statistics data
  file](https://github.com/msberends/AMR/raw/main/data-raw/datasets/dosage.sav)
  (97 kB)  
- Download as [Stata DTA
  file](https://github.com/msberends/AMR/raw/main/data-raw/datasets/dosage.dta)
  (0.2 MB)

**Example content**

| ab  |    name     |       type        |    dose     | dose_times | administration | notes |    original_txt    | eucast_version |
|:---:|:-----------:|:-----------------:|:-----------:|:----------:|:--------------:|:-----:|:------------------:|:--------------:|
| AMK |  Amikacin   |  standard_dosage  | 25-30 mg/kg |     1      |       iv       |       | 25-30 mg/kg x 1 iv |       15       |
| AMX | Amoxicillin |    high_dosage    |     2 g     |     6      |       iv       |       |     2 g x 6 iv     |       15       |
| AMX | Amoxicillin |  standard_dosage  |     1 g     |     3      |       iv       |       |    1 g x 3-4 iv    |       15       |
| AMX | Amoxicillin |    high_dosage    |  0.75-1 g   |     3      |      oral      |       | 0.75-1 g x 3 oral  |       15       |
| AMX | Amoxicillin |  standard_dosage  |    0.5 g    |     3      |      oral      |       |   0.5 g x 3 oral   |       15       |
| AMX | Amoxicillin | uncomplicated_uti |    0.5 g    |     3      |      oral      |       |   0.5 g x 3 oral   |       15       |

------------------------------------------------------------------------

## `example_isolates`: Example Data for Practice

A data set with 2 000 rows and 46 columns, containing the following
column names:  
*date*, *patient*, *age*, *gender*, *ward*, *mo*, *PEN*, *OXA*, *FLC*,
*AMX*, *AMC*, *AMP*, *TZP*, *CZO*, *FEP*, *CXM*, *FOX*, *CTX*, *CAZ*,
*CRO*, *GEN*, *TOB*, *AMK*, *KAN*, *TMP*, *SXT*, *NIT*, *FOS*, *LNZ*,
*CIP*, *MFX*, *VAN*, *TEC*, *TCY*, *TGC*, *DOX*, *ERY*, *CLI*, *AZM*,
*IPM*, *MEM*, *MTR*, *CHL*, *COL*, *MUP*, and *RIF*.

This data set is in R available as `example_isolates`, after you load
the `AMR` package.

It was last updated on 15 June 2024 13:33:49 UTC. Find more info about
the contents, (scientific) source, and structure of this [data set
here](https://amr-for-r.org/reference/example_isolates.html).

**Example content**

|    date    | patient | age | gender |   ward   |      mo      | PEN | OXA | FLC | AMX | AMC | AMP | TZP | CZO | FEP | CXM | FOX | CTX | CAZ | CRO | GEN | TOB | AMK | KAN | TMP | SXT | NIT | FOS | LNZ | CIP | MFX | VAN | TEC | TCY | TGC | DOX | ERY | CLI | AZM | IPM | MEM | MTR | CHL | COL | MUP | RIF |
|:----------:|:-------:|:---:|:------:|:--------:|:------------:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| 2002-01-02 | A77334  | 65  |   F    | Clinical | B_ESCHR_COLI |  R  |     |     |     |  I  |     |     |     |     |  I  |     |     |     |     |     |     |     |     |  R  |  R  |     |     |  R  |     |     |  R  |  R  |  R  |     |     |  R  |  R  |  R  |     |     |     |     |     |     |  R  |
| 2002-01-03 | A77334  | 65  |   F    | Clinical | B_ESCHR_COLI |  R  |     |     |     |  I  |     |     |     |     |  I  |     |     |     |     |     |     |     |     |  R  |  R  |     |     |  R  |     |     |  R  |  R  |  R  |     |     |  R  |  R  |  R  |     |     |     |     |     |     |  R  |
| 2002-01-07 | 067927  | 45  |   F    |   ICU    | B_STPHY_EPDR |  R  |     |  R  |     |     |     |     |     |     |  R  |     |     |  R  |     |     |     |     |     |  S  |  S  |     |     |     |     |     |  S  |     |  S  |  S  |  S  |  R  |     |  R  |     |     |     |     |  R  |     |     |
| 2002-01-07 | 067927  | 45  |   F    |   ICU    | B_STPHY_EPDR |  R  |     |  R  |     |     |     |     |     |     |  R  |     |     |  R  |     |     |     |     |     |  S  |  S  |     |     |     |     |     |  S  |     |  S  |  S  |  S  |  R  |     |  R  |     |     |     |     |  R  |     |     |
| 2002-01-13 | 067927  | 45  |   F    |   ICU    | B_STPHY_EPDR |  R  |     |  R  |     |     |     |     |     |     |  R  |     |     |  R  |     |     |     |     |     |  R  |     |     |     |     |     |     |  S  |     |  S  |  S  |  S  |  R  |     |  R  |     |     |     |     |  R  |     |     |
| 2002-01-13 | 067927  | 45  |   F    |   ICU    | B_STPHY_EPDR |  R  |     |  R  |     |     |     |     |     |     |  R  |     |     |  R  |     |     |     |     |     |  R  |     |     |     |     |     |     |  S  |     |  S  |  S  |  S  |  R  |  R  |  R  |     |     |     |     |  R  |     |     |

------------------------------------------------------------------------

## `example_isolates_unclean`: Example Data for Practice

A data set with 3 000 rows and 8 columns, containing the following
column names:  
*patient_id*, *hospital*, *date*, *bacteria*, *AMX*, *AMC*, *CIP*, and
*GEN*.

This data set is in R available as `example_isolates_unclean`, after you
load the `AMR` package.

It was last updated on 27 August 2022 18:49:37 UTC. Find more info about
the contents, (scientific) source, and structure of this [data set
here](https://amr-for-r.org/reference/example_isolates_unclean.html).

**Example content**

| patient_id | hospital |    date    |   bacteria    | AMX | AMC | CIP | GEN |
|:----------:|:--------:|:----------:|:-------------:|:---:|:---:|:---:|:---:|
|     J3     |    A     | 2012-11-21 |    E. coli    |  R  |  I  |  S  |  S  |
|     R7     |    A     | 2018-04-03 | K. pneumoniae |  R  |  I  |  S  |  S  |
|     P3     |    A     | 2014-09-19 |    E. coli    |  R  |  S  |  S  |  S  |
|    P10     |    A     | 2015-12-10 |    E. coli    |  S  |  I  |  S  |  S  |
|     B7     |    A     | 2015-03-02 |    E. coli    |  S  |  S  |  S  |  S  |
|     W3     |    A     | 2018-03-31 |   S. aureus   |  R  |  S  |  R  |  S  |

------------------------------------------------------------------------

## `microorganisms.codes`: Common Laboratory Codes

A data set with 6 036 rows and 2 columns, containing the following
column names:  
*code* and *mo*.

This data set is in R available as `microorganisms.codes`, after you
load the `AMR` package.

It was last updated on 4 May 2025 16:50:25 UTC. Find more info about the
contents, (scientific) source, and structure of this [data set
here](https://amr-for-r.org/reference/microorganisms.codes.html).

**Direct download links:**

- Download as [original R Data Structure (RDS)
  file](https://github.com/msberends/AMR/raw/main/data-raw/datasets/microorganisms.codes.rds)
  (27 kB)  
- Download as [tab-separated text
  file](https://github.com/msberends/AMR/raw/main/data-raw/datasets/microorganisms.codes.txt)
  (0.1 MB)  
- Download as [Microsoft Excel
  workbook](https://github.com/msberends/AMR/raw/main/data-raw/datasets/microorganisms.codes.xlsx)
  (98 kB)  
- Download as [Apache Feather
  file](https://github.com/msberends/AMR/raw/main/data-raw/datasets/microorganisms.codes.feather)
  (0.1 MB)  
- Download as [Apache Parquet
  file](https://github.com/msberends/AMR/raw/main/data-raw/datasets/microorganisms.codes.parquet)
  (68 kB)  
- Download as [IBM SPSS Statistics data
  file](https://github.com/msberends/AMR/raw/main/data-raw/datasets/microorganisms.codes.sav)
  (0.2 MB)  
- Download as [Stata DTA
  file](https://github.com/msberends/AMR/raw/main/data-raw/datasets/microorganisms.codes.dta)
  (0.2 MB)

**Example content**

| code |      mo      |
|:----:|:------------:|
| 1011 |   B_GRAMP    |
| 1012 |   B_GRAMP    |
| 1013 |   B_GRAMN    |
| 1014 |   B_GRAMN    |
| 1015 |   F_YEAST    |
| 103  | B_ESCHR_COLI |

------------------------------------------------------------------------

## `antivirals`: Antiviral Drugs

A data set with 120 rows and 11 columns, containing the following column
names:  
*av*, *name*, *atc*, *cid*, *atc_group*, *synonyms*, *oral_ddd*,
*oral_units*, *iv_ddd*, *iv_units*, and *loinc*.

This data set is in R available as `antivirals`, after you load the
`AMR` package.

It was last updated on 20 October 2023 12:51:48 UTC. Find more info
about the contents, (scientific) source, and structure of this [data set
here](https://amr-for-r.org/reference/antimicrobials.html).

**Direct download links:**

- Download as [original R Data Structure (RDS)
  file](https://github.com/msberends/AMR/raw/main/data-raw/datasets/antivirals.rds)
  (6 kB)  
- Download as [tab-separated text
  file](https://github.com/msberends/AMR/raw/main/data-raw/datasets/antivirals.txt)
  (17 kB)  
- Download as [Microsoft Excel
  workbook](https://github.com/msberends/AMR/raw/main/data-raw/datasets/antivirals.xlsx)
  (16 kB)  
- Download as [Apache Feather
  file](https://github.com/msberends/AMR/raw/main/data-raw/datasets/antivirals.feather)
  (16 kB)  
- Download as [Apache Parquet
  file](https://github.com/msberends/AMR/raw/main/data-raw/datasets/antivirals.parquet)
  (13 kB)  
- Download as [IBM SPSS Statistics data
  file](https://github.com/msberends/AMR/raw/main/data-raw/datasets/antivirals.sav)
  (32 kB)  
- Download as [Stata DTA
  file](https://github.com/msberends/AMR/raw/main/data-raw/datasets/antivirals.dta)
  (78 kB)

The tab-separated text, Microsoft Excel, SPSS, and Stata files all
contain the trade names and LOINC codes as comma separated values.

**Example content**

| av  |        name        |   atc   |    cid    |                             atc_group                              |                       synonyms                        | oral_ddd | oral_units | iv_ddd | iv_units |            loinc             |
|:---:|:------------------:|:-------:|:---------:|:------------------------------------------------------------------:|:-----------------------------------------------------:|:--------:|:----------:|:------:|:--------:|:----------------------------:|
| ABA |      Abacavir      | J05AF06 |  441300   |     Nucleoside and nucleotide reverse transcriptase inhibitors     |          abacavir sulfate, avacavir, ziagen           |   0.6    |     g      |        |          | 29113-8, 30273-7, 30287-7, … |
| ACI |     Aciclovir      | J05AB01 | 135398513 | Nucleosides and nucleotides excl. reverse transcriptase inhibitors |        acicloftal, aciclovier, aciclovirum, …         |   4.0    |     g      |   4    |    g     |                              |
| ADD | Adefovir dipivoxil | J05AF08 |   60871   |     Nucleoside and nucleotide reverse transcriptase inhibitors     | adefovir di, adefovir di ester, adefovir dipivoxyl, … |   10.0   |     mg     |        |          |                              |
| AME |     Amenamevir     | J05AX26 | 11397521  |                          Other antivirals                          |                       amenalief                       |   0.4    |     g      |        |          |                              |
| AMP |     Amprenavir     | J05AE05 |   65016   |                        Protease inhibitors                         |             agenerase, carbamate, prozei              |   1.2    |     g      |        |          | 29114-6, 30296-8, 30297-6, … |
| ASU |    Asunaprevir     | J05AP06 | 16076883  |             Antivirals for treatment of HCV infections             |                sunvepra, sunvepratrade                |   0.2    |     g      |        |          |                              |
