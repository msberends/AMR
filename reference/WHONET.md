# Data Set with 500 Isolates - WHONET Example

This example data set has the exact same structure as an export file
from WHONET. Such files can be used with this package, as this example
data set shows. The antimicrobial results are from our
[example_isolates](https://amr-for-r.org/reference/example_isolates.md)
data set. All patient names were created using online surname generators
and are only in place for practice purposes.

## Usage

``` r
WHONET
```

## Format

A [tibble](https://tibble.tidyverse.org/reference/tibble.html) with 500
observations and 53 variables:

- `Identification number`  
  ID of the sample

- `Specimen number`  
  ID of the specimen

- `Organism`  
  Name of the microorganism. Before analysis, you should transform this
  to a valid microbial class, using
  [`as.mo()`](https://amr-for-r.org/reference/as.mo.md).

- `Country`  
  Country of origin

- `Laboratory`  
  Name of laboratory

- `Last name`  
  Fictitious last name of patient

- `First name`  
  Fictitious initial of patient

- `Sex`  
  Fictitious gender of patient

- `Age`  
  Fictitious age of patient

- `Age category`  
  Age group, can also be looked up using
  [`age_groups()`](https://amr-for-r.org/reference/age_groups.md)

- `Date of admission`  
  [Date](https://rdrr.io/r/base/Dates.html) of hospital admission

- `Specimen date`  
  [Date](https://rdrr.io/r/base/Dates.html) when specimen was received
  at laboratory

- `Specimen type`  
  Specimen type or group

- `Specimen type (Numeric)`  
  Translation of `"Specimen type"`

- `Reason`  
  Reason of request with Differential Diagnosis

- `Isolate number`  
  ID of isolate

- `Organism type`  
  Type of microorganism, can also be looked up using
  [`mo_type()`](https://amr-for-r.org/reference/mo_property.md)

- `Serotype`  
  Serotype of microorganism

- `Beta-lactamase`  
  Microorganism produces beta-lactamase?

- `ESBL`  
  Microorganism produces extended spectrum beta-lactamase?

- `Carbapenemase`  
  Microorganism produces carbapenemase?

- `MRSA screening test`  
  Microorganism is possible MRSA?

- `Inducible clindamycin resistance`  
  Clindamycin can be induced?

- `Comment`  
  Other comments

- `Date of data entry`  
  [Date](https://rdrr.io/r/base/Dates.html) this data was entered in
  WHONET

- `AMP_ND10:CIP_EE`  
  28 different antimicrobials. You can lookup the abbreviations in the
  [antimicrobials](https://amr-for-r.org/reference/antimicrobials.md)
  data set, or use e.g.
  [`ab_name("AMP")`](https://amr-for-r.org/reference/ab_property.md) to
  get the official name immediately. Before analysis, you should
  transform this to a valid antimicrobial class, using
  [`as.sir()`](https://amr-for-r.org/reference/as.sir.md).

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

## Examples

``` r
WHONET
#> # A tibble: 500 × 53
#>    `Identification number` `Specimen number` Organism Country         Laboratory
#>    <chr>                               <int> <chr>    <chr>           <chr>     
#>  1 fe41d7bafa                           1748 SPN      Belgium         National …
#>  2 91f175ec37                           1767 eco      The Netherlands National …
#>  3 cc4015056e                           1343 eco      The Netherlands National …
#>  4 e864b692f5                           1894 MAP      Denmark         National …
#>  5 3d051fe345                           1739 PVU      Belgium         National …
#>  6 c80762a08d                           1846 103      The Netherlands National …
#>  7 8022d3727c                           1628 103      Denmark         National …
#>  8 f3dc5f553d                           1493 eco      The Netherlands National …
#>  9 15add38f6c                           1847 eco      France          National …
#> 10 fd41248def                           1458 eco      Germany         National …
#> # ℹ 490 more rows
#> # ℹ 48 more variables: `Last name` <chr>, `First name` <chr>, Sex <chr>,
#> #   Age <dbl>, `Age category` <chr>, `Date of admission` <date>,
#> #   `Specimen date` <date>, `Specimen type` <chr>,
#> #   `Specimen type (Numeric)` <dbl>, Reason <chr>, `Isolate number` <int>,
#> #   `Organism type` <chr>, Serotype <chr>, `Beta-lactamase` <lgl>, ESBL <lgl>,
#> #   Carbapenemase <lgl>, `MRSA screening test` <lgl>, …
```
