# Data Sets with 618 Antimicrobial Drugs

Two data sets containing all antimicrobials and antivirals. Use
[`as.ab()`](https://amr-for-r.org/reference/as.ab.md) or one of the
[`ab_*`](https://amr-for-r.org/reference/ab_property.md) functions to
retrieve values from the antimicrobials data set. Three identifiers are
included in this data set: an antimicrobial ID (`ab`, primarily used in
this package) as defined by WHONET/EARS-Net, an ATC code (`atc`) as
defined by the WHO, and a Compound ID (`cid`) as found in PubChem. Other
properties in this data set are derived from one or more of these codes.
Note that some drugs have multiple ATC codes.

**The `antibiotics` data set has been renamed to `antimicrobials`. The
old name will be removed in a future version.**

## Usage

``` r
antimicrobials

antibiotics

antivirals
```

## Format

### For the antimicrobials data set: a [tibble](https://tibble.tidyverse.org/reference/tibble.html) with 498 observations and 14 variables:

- `ab`  
  antimicrobial ID as used in this package (such as `AMC`), using the
  official EARS-Net (European Antimicrobial Resistance Surveillance
  Network) codes where available. ***This is a unique identifier.***

- `cid`  
  Compound ID as found in PubChem. ***This is a unique identifier.***

- `name`  
  Official name as used by WHONET/EARS-Net or the WHO. ***This is a
  unique identifier.***

- `group`  
  One or more short and concise group names, based on WHONET and WHOCC
  definitions

- `atc`  
  ATC codes (Anatomical Therapeutic Chemical) as defined by the WHOCC,
  like `J01CR02` (last updated May 4th, 2025):

- `atc_group1`  
  Official pharmacological subgroup (3rd level ATC code) as defined by
  the WHOCC, like `"Macrolides, lincosamides and streptogramins"`

- `atc_group2`  
  Official chemical subgroup (4th level ATC code) as defined by the
  WHOCC, like `"Macrolides"`

- `abbr`  
  List of abbreviations as used in many countries, also for
  antimicrobial susceptibility testing (AST)

- `synonyms`  
  Synonyms (often trade names) of a drug, as found in PubChem based on
  their compound ID

ATC properties (last updated May 4th, 2025):

- `oral_ddd`  
  Defined Daily Dose (DDD), oral treatment, currently available for 180
  drugs

- `oral_units`  
  Units of `oral_ddd`

- `iv_ddd`  
  Defined Daily Dose (DDD), parenteral (intravenous) treatment,
  currently available for 153 drugs

- `iv_units`  
  Units of `iv_ddd`

LOINC:

- `loinc`  
  All codes associated with the name of the antimicrobial drug from
  Logical Observation Identifiers Names and Codes (LOINC), Version 2.76
  (18 September, 2023). Use
  [`ab_loinc()`](https://amr-for-r.org/reference/ab_property.md) to
  retrieve them quickly, see
  [`ab_property()`](https://amr-for-r.org/reference/ab_property.md).

### For the antivirals data set: a [tibble](https://tibble.tidyverse.org/reference/tibble.html) with 120 observations and 11 variables:

- `av`  
  Antiviral ID as used in this package (such as `ACI`), using the
  official EARS-Net (European Antimicrobial Resistance Surveillance
  Network) codes where available. ***This is a unique identifier.***
  Combinations are codes that contain a `+` to indicate this, such as
  `ATA+COBI` for atazanavir/cobicistat.

- `name`  
  Official name as used by WHONET/EARS-Net or the WHO. ***This is a
  unique identifier.***

- `atc`  
  ATC codes (Anatomical Therapeutic Chemical) as defined by the WHOCC,
  see *Details*

- `cid`  
  Compound ID as found in PubChem. ***This is a unique identifier.***

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
  [`av_loinc()`](https://amr-for-r.org/reference/av_property.md) to
  retrieve them quickly, see
  [`av_property()`](https://amr-for-r.org/reference/av_property.md).

An object of class `deprecated_amr_dataset` (inherits from `tbl_df`,
`tbl`, `data.frame`) with 498 rows and 14 columns.

An object of class `tbl_df` (inherits from `tbl`, `data.frame`) with 120
rows and 11 columns.

## Source

- WHO Collaborating Centre for Drug Statistics Methodology, Guidelines
  for ATC classification and DDD assignment, Oslo Accessed from
  <https://atcddd.fhi.no/atc_ddd_index/> on May 4th, 2025.

- Logical Observation Identifiers Names and Codes (LOINC), Version 2.76
  (18 September, 2023). Accessed from <https://loinc.org> on October
  19th, 2023.

- European Commission Public Health PHARMACEUTICALS - COMMUNITY
  REGISTER:
  <https://ec.europa.eu/health/documents/community-register/html/reg_hum_atc.htm>

## Details

Properties that are based on an ATC code are only available when an ATC
is available. These properties are: `atc_group1`, `atc_group2`,
`oral_ddd`, `oral_units`, `iv_ddd` and `iv_units`. Do note that ATC
codes are not unique. For example, J01CR02 is officially the ATC code
for "amoxicillin and beta-lactamase inhibitor". Consequently, these two
items from the antimicrobials data set both return `"J01CR02"`:

    ab_atc("amoxicillin/clavulanic acid")
    ab_atc("amoxicillin/sulbactam")

Synonyms (i.e. trade names) were derived from the PubChem Compound ID
(column `cid`) and are consequently only available where a CID is
available.

## Download Our Reference Data

All reference data sets in the AMR package - including information on
microorganisms, antimicrobials, and clinical breakpoints - are freely
available for download in multiple formats: R, MS Excel, Apache Feather,
Apache Parquet, SPSS, and Stata.

For maximum compatibility, we also provide machine-readable,
tab-separated plain text files suitable for use in any software,
including laboratory information systems.

Visit [our website for direct download
links](https://amr-for-r.org/articles/datasets.html), or explore the
actual files in [our GitHub
repository](https://github.com/msberends/AMR/tree/main/data-raw/datasets).

## WHOCC

This package contains **all ~550 antibiotic, antimycotic and antiviral
drugs** and their Anatomical Therapeutic Chemical (ATC) codes, ATC
groups and Defined Daily Dose (DDD) from the World Health Organization
Collaborating Centre for Drug Statistics Methodology (WHOCC,
<https://atcddd.fhi.no>) and the Pharmaceuticals Community Register of
the European Commission
(<https://ec.europa.eu/health/documents/community-register/html/reg_hum_atc.htm>).

These have become the gold standard for international drug utilisation
monitoring and research.

The WHOCC is located in Oslo at the Norwegian Institute of Public Health
and funded by the Norwegian government. The European Commission is the
executive of the European Union and promotes its general interest.

**NOTE: The WHOCC copyright does not allow use for commercial purposes,
unlike any other info from this package.** See
<https://atcddd.fhi.no/copyright_disclaimer/.>

## See also

[microorganisms](https://amr-for-r.org/reference/microorganisms.md),
[intrinsic_resistant](https://amr-for-r.org/reference/intrinsic_resistant.md)

## Examples

``` r
antimicrobials
#> # A tibble: 498 × 14
#>    ab        cid name   group atc   atc_group1 atc_group2 abbreviations synonyms
#>    <ab>    <dbl> <chr>  <lis> <lis> <chr>      <chr>      <list>        <named >
#>  1 AMA      4649 4-ami… <chr> <chr> Drugs for… Aminosali… <chr [1]>     <chr>   
#>  2 ACM   6450012 Acety… <chr> <chr> NA         NA         <chr [1]>     <chr>   
#>  3 ASP  49787020 Acety… <chr> <chr> NA         NA         <chr [1]>     <chr>   
#>  4 ALS      8954 Aldes… <chr> <chr> Drugs for… Drugs for… <chr [1]>     <chr>   
#>  5 AMK     37768 Amika… <chr> <chr> Aminoglyc… Other ami… <chr [6]>     <chr>   
#>  6 AKF        NA Amika… <chr> <chr> NA         NA         <chr [1]>     <chr>   
#>  7 AMO     54260 Amoro… <chr> <chr> Antifunga… Other ant… <chr [1]>     <chr>   
#>  8 AMX     33613 Amoxi… <chr> <chr> Beta-lact… Penicilli… <chr [4]>     <chr>   
#>  9 AMC  23665637 Amoxi… <chr> <chr> Beta-lact… Combinati… <chr [6]>     <chr>   
#> 10 AXS    465441 Amoxi… <chr> <chr> NA         NA         <chr [1]>     <chr>   
#> # ℹ 488 more rows
#> # ℹ 5 more variables: oral_ddd <dbl>, oral_units <chr>, iv_ddd <dbl>,
#> #   iv_units <chr>, loinc <list>
antivirals
#> # A tibble: 120 × 11
#>    av       name      atc      cid atc_group synonyms oral_ddd oral_units iv_ddd
#>    <av>     <chr>     <chr>  <dbl> <chr>     <list>      <dbl> <chr>       <dbl>
#>  1 ABA      Abacavir  J05A… 4.41e5 Nucleosi… <chr>         0.6 g              NA
#>  2 ACI      Aciclovir J05A… 1.35e8 Nucleosi… <chr>         4   g               4
#>  3 ADD      Adefovir… J05A… 6.09e4 Nucleosi… <chr>        10   mg             NA
#>  4 AME      Amenamev… J05A… 1.14e7 Other an… <chr>         0.4 g              NA
#>  5 AMP      Amprenav… J05A… 6.50e4 Protease… <chr>         1.2 g              NA
#>  6 ASU      Asunapre… J05A… 1.61e7 Antivira… <chr>         0.2 g              NA
#>  7 ATA      Atazanav… J05A… 1.48e5 Protease… <chr>         0.3 g              NA
#>  8 ATA+COBI Atazanav… J05A… 8.66e7 Antivira… <chr>        NA   NA             NA
#>  9 ATA+RIT  Atazanav… J05A… 2.51e7 Antivira… <chr>         0.3 g              NA
#> 10 BAM      Baloxavi… J05A… 1.24e8 Other an… <chr>        40   mg             NA
#> # ℹ 110 more rows
#> # ℹ 2 more variables: iv_units <chr>, loinc <list>
```
