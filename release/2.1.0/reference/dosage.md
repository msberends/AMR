# Data Set with Treatment Dosages as Defined by EUCAST

EUCAST breakpoints used in this package are based on the dosages in this
data set. They can be retrieved with
[`eucast_dosage()`](https://amr-for-r.org/release/2.1.0/reference/eucast_rules.md).

## Usage

``` r
dosage
```

## Format

A [tibble](https://tibble.tidyverse.org/reference/tibble.html) with 503
observations and 9 variables:

- `ab`  
  Antibiotic ID as used in this package (such as `AMC`), using the
  official EARS-Net (European Antimicrobial Resistance Surveillance
  Network) codes where available

- `name`  
  Official name of the antimicrobial drug as used by WHONET/EARS-Net or
  the WHO

- `type`  
  Type of the dosage, either "high_dosage", "standard_dosage", or
  "uncomplicated_uti"

- `dose`  
  Dose, such as "2 g" or "25 mg/kg"

- `dose_times`  
  Number of times a dose must be administered

- `administration`  
  Route of administration, either "im", "iv", or "oral"

- `notes`  
  Additional dosage notes

- `original_txt`  
  Original text in the PDF file of EUCAST

- `eucast_version`  
  Version number of the EUCAST Clinical Breakpoints guideline to which
  these dosages apply, either 13, 12, or 11

## Details

Like all data sets in this package, this data set is publicly available
for download in the following formats: R, MS Excel, Apache Feather,
Apache Parquet, SPSS, SAS, and Stata. Please visit [our website for the
download links](https://msberends.github.io/AMR/articles/datasets.html).
The actual files are of course available on [our GitHub
repository](https://github.com/msberends/AMR/tree/main/data-raw).

## Examples

``` r
dosage
```
