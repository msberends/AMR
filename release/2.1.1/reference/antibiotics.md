# Data Sets with 603 Antimicrobial Drugs

Two data sets containing all antibiotics/antimycotics and antivirals.
Use [`as.ab()`](https://amr-for-r.org/release/2.1.1/reference/as.ab.md)
or one of the
[`ab_*`](https://amr-for-r.org/release/2.1.1/reference/ab_property.md)
functions to retrieve values from the antibiotics data set. Three
identifiers are included in this data set: an antibiotic ID (`ab`,
primarily used in this package) as defined by WHONET/EARS-Net, an ATC
code (`atc`) as defined by the WHO, and a Compound ID (`cid`) as found
in PubChem. Other properties in this data set are derived from one or
more of these codes. Note that some drugs have multiple ATC codes.

## Usage

``` r
antibiotics

antivirals
```

## Format

### For the antibiotics data set: a [tibble](https://tibble.tidyverse.org/reference/tibble.html) with 483 observations and 14 variables:

- `ab`  
  Antibiotic ID as used in this package (such as `AMC`), using the
  official EARS-Net (European Antimicrobial Resistance Surveillance
  Network) codes where available. *This is a unique identifier.*

- `cid`  
  Compound ID as found in PubChem. *This is a unique identifier.*

- `name`  
  Official name as used by WHONET/EARS-Net or the WHO. *This is a unique
  identifier.*

- `group`  
  A short and concise group name, based on WHONET and WHOCC definitions

- `atc`  
  ATC codes (Anatomical Therapeutic Chemical) as defined by the WHOCC,
  like `J01CR02`

- `atc_group1`  
  Official pharmacological subgroup (3rd level ATC code) as defined by
  the WHOCC, like `"Macrolides, lincosamides and streptogramins"`

- `atc_group2`  
  Official chemical subgroup (4th level ATC code) as defined by the
  WHOCC, like `"Macrolides"`

- `abbr`  
  List of abbreviations as used in many countries, also for antibiotic
  susceptibility testing (AST)

- `synonyms`  
  Synonyms (often trade names) of a drug, as found in PubChem based on
  their compound ID

- `oral_ddd`  
  Defined Daily Dose (DDD), oral treatment, currently available for 174
  drugs

- `oral_units`  
  Units of `oral_ddd`

- `iv_ddd`  
  Defined Daily Dose (DDD), parenteral (intravenous) treatment,
  currently available for 146 drugs

- `iv_units`  
  Units of `iv_ddd`

- `loinc`  
  All codes associated with the name of the antimicrobial drug from
  Logical Observation Identifiers Names and Codes (LOINC), Version 2.76
  (18 September, 2023). Use
  [`ab_loinc()`](https://amr-for-r.org/release/2.1.1/reference/ab_property.md)
  to retrieve them quickly, see
  [`ab_property()`](https://amr-for-r.org/release/2.1.1/reference/ab_property.md).

### For the antivirals data set: a [tibble](https://tibble.tidyverse.org/reference/tibble.html) with 120 observations and 11 variables:

- `av`  
  Antiviral ID as used in this package (such as `ACI`), using the
  official EARS-Net (European Antimicrobial Resistance Surveillance
  Network) codes where available. *This is a unique identifier.*
  Combinations are codes that contain a `+` to indicate this, such as
  `ATA+COBI` for atazanavir/cobicistat.

- `name`  
  Official name as used by WHONET/EARS-Net or the WHO. *This is a unique
  identifier.*

- `atc`  
  ATC codes (Anatomical Therapeutic Chemical) as defined by the WHOCC

- `cid`  
  Compound ID as found in PubChem. *This is a unique identifier.*

- `atc_group`  
  Official pharmacological subgroup (3rd level ATC code) as defined by
  the WHOCC

- `synonyms`  
  Synonyms (often trade names) of a drug, as found in PubChem based on
  their compound ID

- `oral_ddd`  
  Defined Daily Dose (DDD), oral treatment

- `oral_units`  
  Units of `oral_ddd`

- `iv_ddd`  
  Defined Daily Dose (DDD), parenteral treatment

- `iv_units`  
  Units of `iv_ddd`

- `loinc`  
  All codes associated with the name of the antiviral drug from Logical
  Observation Identifiers Names and Codes (LOINC), Version 2.76 (18
  September, 2023). Use
  [`av_loinc()`](https://amr-for-r.org/release/2.1.1/reference/av_property.md)
  to retrieve them quickly, see
  [`av_property()`](https://amr-for-r.org/release/2.1.1/reference/av_property.md).

An object of class `tbl_df` (inherits from `tbl`, `data.frame`) with 120
rows and 11 columns.

## Source

- World Health Organization (WHO) Collaborating Centre for Drug
  Statistics Methodology (WHOCC): <https://www.whocc.no/atc_ddd_index/>

- Logical Observation Identifiers Names and Codes (LOINC), Version 2.76
  (18 September, 2023). Accessed from <https://loinc.org> on October
  19th, 2023.

- European Commission Public Health PHARMACEUTICALS - COMMUNITY
  REGISTER:
  <https://ec.europa.eu/health/documents/community-register/html/reg_hum_atc.htm>

## Details

Properties that are based on an ATC code are only available when an ATC
is available. These properties are: `atc_group1`, `atc_group2`,
`oral_ddd`, `oral_units`, `iv_ddd` and `iv_units`.

Synonyms (i.e. trade names) were derived from the PubChem Compound ID
(column `cid`) and consequently only available where a CID is available.

### Direct download

Like all data sets in this package, these data sets are publicly
available for download in the following formats: R, MS Excel, Apache
Feather, Apache Parquet, SPSS, SAS, and Stata. Please visit [our website
for the download
links](https://msberends.github.io/AMR/articles/datasets.html). The
actual files are of course available on [our GitHub
repository](https://github.com/msberends/AMR/tree/main/data-raw).

## WHOCC

This package contains **all ~550 antibiotic, antimycotic and antiviral
drugs** and their Anatomical Therapeutic Chemical (ATC) codes, ATC
groups and Defined Daily Dose (DDD) from the World Health Organization
Collaborating Centre for Drug Statistics Methodology (WHOCC,
<https://www.whocc.no>) and the Pharmaceuticals Community Register of
the European Commission
(<https://ec.europa.eu/health/documents/community-register/html/reg_hum_atc.htm>).

These have become the gold standard for international drug utilisation
monitoring and research.

The WHOCC is located in Oslo at the Norwegian Institute of Public Health
and funded by the Norwegian government. The European Commission is the
executive of the European Union and promotes its general interest.

**NOTE: The WHOCC copyright does not allow use for commercial purposes,
unlike any other info from this package.** See
<https://www.whocc.no/copyright_disclaimer/.>

## See also

[microorganisms](https://amr-for-r.org/release/2.1.1/reference/microorganisms.md),
[intrinsic_resistant](https://amr-for-r.org/release/2.1.1/reference/intrinsic_resistant.md)

## Examples

``` r
antibiotics
antivirals
```
