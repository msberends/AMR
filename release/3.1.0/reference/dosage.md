# Data Set with Treatment Dosages as Defined by EUCAST

EUCAST breakpoints used in this package are based on the dosages in this
data set. They can be retrieved with
[`eucast_dosage()`](https://amr-for-r.org/release/3.1.0/reference/interpretive_rules.md).

## Usage

``` r
dosage
```

## Format

A [tibble](https://tibble.tidyverse.org/reference/tibble.html) with 759
observations and 9 variables:

- `ab`  
  Antimicrobial ID as used in this package (such as `AMC`), using the
  official EARS-Net (European Antimicrobial Resistance Surveillance
  Network) codes where available

- `name`  
  Official name of the antimicrobial drug as used by WHONET/EARS-Net or
  the WHO

- `type`  
  Type of the dosage, either `"high_dosage"`, `"standard_dosage"`, or
  `"uncomplicated_uti"`

- `dose`  
  Dose, such as "2 g" or "25 mg/kg"

- `dose_times`  
  Number of times a dose must be administered

- `administration`  
  Route of administration, either `""`, `"im"`, `"iv"`, `"oral"`, or NA

- `notes`  
  Additional dosage notes

- `original_txt`  
  Original text in the PDF file of EUCAST

- `eucast_version`  
  Version number of the EUCAST Clinical Breakpoints guideline to which
  these dosages apply, either 15, 14, 13.1, 12, or 11

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
dosage
```
