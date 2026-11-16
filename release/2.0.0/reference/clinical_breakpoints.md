# Data Set with Clinical Breakpoints for SIR Interpretation

Data set containing clinical breakpoints to interpret MIC and disk
diffusion to SIR values, according to international guidelines.
Currently implemented guidelines are EUCAST (2013-2022) and CLSI
(2013-2022). Use
[`as.sir()`](https://amr-for-r.org/release/2.0.0/reference/as.sir.md) to
transform MICs or disks measurements to SIR values.

## Usage

``` r
clinical_breakpoints
```

## Format

A [tibble](https://tibble.tidyverse.org/reference/tibble.html) with 18
308 observations and 11 variables:

- `guideline`  
  Name of the guideline

- `method`  
  Either "DISK" or "MIC"

- `site`  
  Body site, e.g. "Oral" or "Respiratory"

- `mo`  
  Microbial ID, see
  [`as.mo()`](https://amr-for-r.org/release/2.0.0/reference/as.mo.md)

- `rank_index`  
  Taxonomic rank index of `mo` from 1 (subspecies/infraspecies) to 5
  (unknown microorganism)

- `ab`  
  Antibiotic ID, see
  [`as.ab()`](https://amr-for-r.org/release/2.0.0/reference/as.ab.md)

- `ref_tbl`  
  Info about where the guideline rule can be found

- `disk_dose`  
  Dose of the used disk diffusion method

- `breakpoint_S`  
  Lowest MIC value or highest number of millimetres that leads to "S"

- `breakpoint_R`  
  Highest MIC value or lowest number of millimetres that leads to "R"

- `uti`  
  A [logical](https://rdrr.io/r/base/logical.html) value
  (`TRUE`/`FALSE`) to indicate whether the rule applies to a urinary
  tract infection (UTI)

## Details

Like all data sets in this package, this data set is publicly available
for download in the following formats: R, MS Excel, Apache Feather,
Apache Parquet, SPSS, SAS, and Stata. Please visit [our website for the
download links](https://msberends.github.io/AMR/articles/datasets.html).
The actual files are of course available on [our GitHub
repository](https://github.com/msberends/AMR/tree/main/data-raw).

They **allow for machine reading EUCAST and CLSI guidelines**, which is
almost impossible with the MS Excel and PDF files distributed by EUCAST
and CLSI.

## See also

[intrinsic_resistant](https://amr-for-r.org/release/2.0.0/reference/intrinsic_resistant.md)

## Examples

``` r
clinical_breakpoints
```
