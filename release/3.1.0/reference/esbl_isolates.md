# Data Set with 500 ESBL Isolates

A data set containing 500 microbial isolates with MIC values of common
antibiotics and a binary `esbl` column for extended-spectrum
beta-lactamase (ESBL) production. This data set contains randomised
fictitious data but reflects reality and can be used to practise
AMR-related machine learning, e.g., classification modelling with
[tidymodels](https://amr-for-r.org/articles/AMR_with_tidymodels.html).

## Usage

``` r
esbl_isolates
```

## Format

A [tibble](https://tibble.tidyverse.org/reference/tibble.html) with 500
observations and 19 variables:

- `esbl`  
  Logical indicator if the isolate is ESBL-producing

- `genus`  
  Genus of the microorganism

- `AMC:COL`  
  MIC values for 17 antimicrobial drugs, transformed to class
  [`mic`](https://amr-for-r.org/release/3.1.0/reference/as.mic.md) (see
  [`as.mic()`](https://amr-for-r.org/release/3.1.0/reference/as.mic.md))

## Details

See our [tidymodels
integration](https://amr-for-r.org/release/3.1.0/reference/amr-tidymodels.md)
for an example using this data set.

## Examples

``` r
esbl_isolates
```
