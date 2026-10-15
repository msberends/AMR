# Data Set with 6 293 Common Microorganism Codes

A data set containing commonly used codes for microorganisms, from
laboratory systems and [WHONET](https://whonet.org). Define your own
with [`set_mo_source()`](https://amr-for-r.org/reference/mo_source.md).
They will all be searched when using
[`as.mo()`](https://amr-for-r.org/reference/as.mo.md) and consequently
all the [`mo_*`](https://amr-for-r.org/reference/mo_property.md)
functions.

## Usage

``` r
microorganisms.codes
```

## Format

A [tibble](https://tibble.tidyverse.org/reference/tibble.html) with 6
293 observations and 3 variables:

- `code`  
  Commonly used code of a microorganism. ***This is a unique
  identifier.***

- `mo`  
  ID of the microorganism in the
  [microorganisms](https://amr-for-r.org/reference/microorganisms.md)
  data set. This is the record of the name that the code was named
  after, also if that name is outdated, so that the code and the
  microorganism ID remain relatable: e.g., the WHONET code `EAE` refers
  to *Enterobacter aerogenes* (`B_ENTRBC_AERG`).
  [`as.mo()`](https://amr-for-r.org/reference/as.mo.md) returns the
  current name, here *Klebsiella aerogenes*. Codes for organisms that
  are not in the
  [microorganisms](https://amr-for-r.org/reference/microorganisms.md)
  data set, such as viruses, refer to `UNKNOWN`.

- `source`  
  Source of the code, with the name it refers to, e.g.
  `"WHONET (Enterobacter aerogenes)"`, or `NA` if the source is unknown

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

## See also

[`as.mo()`](https://amr-for-r.org/reference/as.mo.md)
[microorganisms](https://amr-for-r.org/reference/microorganisms.md)

## Examples

``` r
microorganisms.codes
#> # A tibble: 6,293 × 3
#>    code  mo           source                                
#>    <chr> <mo>         <chr>                                 
#>  1 1011  B_GRAMP      NA                                    
#>  2 1012  B_GRAMP      NA                                    
#>  3 1013  B_GRAMN      NA                                    
#>  4 1014  B_GRAMN      NA                                    
#>  5 1015  F_YEAST      NA                                    
#>  6 103   B_ESCHR_COLI WHONET (Escherichia coli O103)        
#>  7 104   B_SLMNL_HMRM WHONET (Salmonella Typhimurium DT 104)
#>  8 1100  B_STRPT      NA                                    
#>  9 1101  B_STRPT_VIRI NA                                    
#> 10 1102  B_STRPT_HAEM NA                                    
#> # ℹ 6,283 more rows

# 'ECO' or 'eco' is the WHONET code for E. coli:
microorganisms.codes[microorganisms.codes$code == "ECO", ]
#> # A tibble: 1 × 3
#>   code  mo           source                   
#>   <chr> <mo>         <chr>                    
#> 1 ECO   B_ESCHR_COLI WHONET (Escherichia coli)

# and therefore, 'eco' will be understood as E. coli in this package:
mo_info("eco")
#> $mo
#> [1] "B_ESCHR_COLI"
#> 
#> $rank
#> [1] "species"
#> 
#> $domain
#> [1] "Bacteria"
#> 
#> $kingdom
#> [1] "Pseudomonadati"
#> 
#> $phylum
#> [1] "Pseudomonadota"
#> 
#> $class
#> [1] "Gammaproteobacteria"
#> 
#> $order
#> [1] "Enterobacterales"
#> 
#> $family
#> [1] "Enterobacteriaceae"
#> 
#> $genus
#> [1] "Escherichia"
#> 
#> $species
#> [1] "coli"
#> 
#> $subspecies
#> [1] ""
#> 
#> $status
#> [1] "accepted"
#> 
#> $synonyms
#>    Migula, 1895 
#> "Bacillus coli" 
#> 
#> $gramstain
#> [1] "Gram-negative"
#> 
#> $morphology
#> [1] "Rods"
#> 
#> $oxygen_tolerance
#> [1] "facultative anaerobe"
#> 
#> $url
#> [1] "https://lpsn.dsmz.de/species/escherichia-coli"
#> 
#> $ref
#> [1] "Castellani et al., 1919"
#> 
#> $snomed
#>  [1] "1095001000112106" "715307006"        "737528008"        "416989002"       
#>  [5] "116397003"        "414097009"        "414098004"        "414099007"       
#>  [9] "414100004"        "116395006"        "735270003"        "116396007"       
#> [13] "83285000"         "116394005"        "112283007"        "710886005"       
#> [17] "710887001"        "710888006"        "710889003"        "414132004"       
#> [21] "721892009"        "416812001"        "416740004"        "417216001"       
#> [25] "457541006"        "710253004"        "416530004"        "417189006"       
#> [29] "409800005"        "713925008"        "444771000124108"  "838549008"       
#> 
#> $lpsn
#> [1] "776057"
#> 
#> $mycobank
#> [1] NA
#> 
#> $gbif
#> [1] "CS33N"
#> 
#> $group_members
#> character(0)
#> 

# works for all AMR functions:
mo_is_intrinsic_resistant("eco", ab = "vancomycin")
#> [1] TRUE
```
