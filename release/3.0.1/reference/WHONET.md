# Data Set with 500 Isolates - WHONET Example

This example data set has the exact same structure as an export file
from WHONET. Such files can be used with this package, as this example
data set shows. The antimicrobial results are from our
[example_isolates](https://amr-for-r.org/release/3.0.1/reference/example_isolates.md)
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
  [`as.mo()`](https://amr-for-r.org/release/3.0.1/reference/as.mo.md).

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
  [`age_groups()`](https://amr-for-r.org/release/3.0.1/reference/age_groups.md)

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
  [`mo_type()`](https://amr-for-r.org/release/3.0.1/reference/mo_property.md)

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
  [antimicrobials](https://amr-for-r.org/release/3.0.1/reference/antimicrobials.md)
  data set, or use e.g.
  [`ab_name("AMP")`](https://amr-for-r.org/release/3.0.1/reference/ab_property.md)
  to get the official name immediately. Before analysis, you should
  transform this to a valid antimicrobial class, using
  [`as.sir()`](https://amr-for-r.org/release/3.0.1/reference/as.sir.md).

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
```
