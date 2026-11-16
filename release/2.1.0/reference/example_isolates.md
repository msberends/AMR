# Data Set with 2 000 Example Isolates

A data set containing 2 000 microbial isolates with their full
antibiograms. This data set contains randomised fictitious data, but
reflects reality and can be used to practise AMR data analysis. For
examples, please read [the tutorial on our
website](https://msberends.github.io/AMR/articles/AMR.html).

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
  [`as.mo()`](https://amr-for-r.org/release/2.1.0/reference/as.mo.md),
  see also the
  [microorganisms](https://amr-for-r.org/release/2.1.0/reference/microorganisms.md)
  data set

- `PEN:RIF`  
  40 different antibiotics with class
  [`sir`](https://amr-for-r.org/release/2.1.0/reference/as.sir.md) (see
  [`as.sir()`](https://amr-for-r.org/release/2.1.0/reference/as.sir.md));
  these column names occur in the
  [antibiotics](https://amr-for-r.org/release/2.1.0/reference/antibiotics.md)
  data set and can be translated with
  [`set_ab_names()`](https://amr-for-r.org/release/2.1.0/reference/ab_property.md)
  or
  [`ab_name()`](https://amr-for-r.org/release/2.1.0/reference/ab_property.md)

## Details

Like all data sets in this package, this data set is publicly available
for download in the following formats: R, MS Excel, Apache Feather,
Apache Parquet, SPSS, SAS, and Stata. Please visit [our website for the
download links](https://msberends.github.io/AMR/articles/datasets.html).
The actual files are of course available on [our GitHub
repository](https://github.com/msberends/AMR/tree/main/data-raw).

## Examples

``` r
example_isolates
```
