# Data Set with 2 000 Example Isolates

A data set containing 2 000 microbial isolates with their full
antibiograms. This data set contains randomised fictitious data, but
reflects reality and can be used to practise AMR data analysis. For
examples, please read [the tutorial on our
website](https://amr-for-r.org/articles/AMR.html).

## Usage

``` r
example_isolates
```

## Format

A [tibble](https://tibble.tidyverse.org/reference/tibble.html) with 2
000 observations and 46 variables:

- `date`  
  Date of receipt at the laboratory

- `patient`  
  ID of the patient

- `age`  
  Age of the patient

- `gender`  
  Gender of the patient, either "F" or "M"

- `ward`  
  Ward type where the patient was admitted, either "Clinical", "ICU", or
  "Outpatient"

- `mo`  
  ID of microorganism created with
  [`as.mo()`](https://amr-for-r.org/release/3.0.0/reference/as.mo.md),
  see also the
  [microorganisms](https://amr-for-r.org/release/3.0.0/reference/microorganisms.md)
  data set

- `PEN:RIF`  
  40 different antimicrobials with class
  [`sir`](https://amr-for-r.org/release/3.0.0/reference/as.sir.md) (see
  [`as.sir()`](https://amr-for-r.org/release/3.0.0/reference/as.sir.md));
  these column names occur in the
  [antimicrobials](https://amr-for-r.org/release/3.0.0/reference/antimicrobials.md)
  data set and can be translated with
  [`set_ab_names()`](https://amr-for-r.org/release/3.0.0/reference/ab_property.md)
  or
  [`ab_name()`](https://amr-for-r.org/release/3.0.0/reference/ab_property.md)

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
example_isolates
```
