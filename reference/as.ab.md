# Transform Input to an Antibiotic ID

Use this function to determine the antimicrobial drug code of one or
more antimicrobials. The data set
[antimicrobials](https://amr-for-r.org/reference/antimicrobials.md) will
be searched for abbreviations, official names and synonyms (brand
names).

## Usage

``` r
as.ab(x, flag_multiple_results = TRUE, language = get_AMR_locale(),
  info = interactive(), ...)

is.ab(x)

ab_reset_session()
```

## Arguments

- x:

  A [character](https://rdrr.io/r/base/character.html) vector to
  determine to antibiotic ID.

- flag_multiple_results:

  A [logical](https://rdrr.io/r/base/logical.html) to indicate whether a
  note should be printed to the console that probably more than one
  antibiotic drug code or name can be retrieved from a single input
  value.

- language:

  Language to coerce input values from any of the 28 supported
  languages - default to the system language if supported (see
  [`get_AMR_locale()`](https://amr-for-r.org/reference/translate.md)).

- info:

  A [logical](https://rdrr.io/r/base/logical.html) to indicate whether a
  progress bar should be printed - the default is `TRUE` only in
  interactive mode.

- ...:

  Arguments passed on to internal functions.

## Value

A [character](https://rdrr.io/r/base/character.html)
[vector](https://rdrr.io/r/base/vector.html) with additional class `ab`

## Details

All entries in the
[antimicrobials](https://amr-for-r.org/reference/antimicrobials.md) data
set have three different identifiers: a human readable EARS-Net code
(column `ab`, used by ECDC and WHONET), an ATC code (column `atc`, used
by WHO), and a CID code (column `cid`, Compound ID, used by PubChem).
The data set contains more than 5,000 official brand names from many
different countries, as found in PubChem. Not that some drugs contain
multiple ATC codes.

All these properties will be searched for the user input. The `as.ab()`
can correct for different forms of misspelling:

- Wrong spelling of drug names (such as "tobramicin" or "gentamycin"),
  which corrects for most audible similarities such as f/ph, x/ks,
  c/z/s, t/th, etc.

- Too few or too many vowels or consonants

- Switching two characters (such as "mreopenem", often the case in
  clinical data, when doctors typed too fast)

- Digitalised paper records, leaving artefacts like 0/o/O (zero and
  O's), B/8, n/r, etc.

Use the [`ab_*`](https://amr-for-r.org/reference/ab_property.md)
functions to get properties based on the returned antibiotic ID, see
*Examples*.

Note: the `as.ab()` and
[`ab_*`](https://amr-for-r.org/reference/ab_property.md) functions may
use very long regular expression to match brand names of antimicrobial
drugs. This may fail on some systems.

You can add your own manual codes to be considered by `as.ab()` and all
[`ab_*`](https://amr-for-r.org/reference/ab_property.md) functions, see
[`add_custom_antimicrobials()`](https://amr-for-r.org/reference/add_custom_antimicrobials.md).

## Source

World Health Organization (WHO) Collaborating Centre for Drug Statistics
Methodology: <https://atcddd.fhi.no/atc_ddd_index/>

European Commission Public Health PHARMACEUTICALS - COMMUNITY REGISTER:
<https://health.ec.europa.eu/documents/community-register/html/reg_hum_atc.htm>

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

- [antimicrobials](https://amr-for-r.org/reference/antimicrobials.md)
  for the [data.frame](https://rdrr.io/r/base/data.frame.html) that is
  being used to determine ATCs

- [`ab_from_text()`](https://amr-for-r.org/reference/ab_from_text.md)
  for a function to retrieve antimicrobial drugs from clinical text
  (from health care records)

## Examples

``` r
# these examples all return "ERY", the ID of erythromycin:
as.ab("J01FA01")
#> Class 'ab'
#> [1] ERY
as.ab("J 01 FA 01")
#> Class 'ab'
#> [1] ERY
as.ab("Erythromycin")
#> Class 'ab'
#> [1] ERY
as.ab("eryt")
#> Class 'ab'
#> [1] ERY
as.ab("ERYT")
#> Class 'ab'
#> [1] ERY
as.ab("ERY")
#> Class 'ab'
#> [1] ERY
as.ab("eritromicine") # spelled wrong, yet works
#> Class 'ab'
#> [1] ERY
as.ab("Erythrocin") # trade name
#> Class 'ab'
#> [1] ERY

# spelling from different languages and dyslexia are no problem
ab_atc("ceftriaxon")
#> [1] "J01DD04"  "QJ01DD04"
ab_atc("cephtriaxone") # small spelling error
#> [1] "J01DD04"  "QJ01DD04"
ab_atc("cephthriaxone") # or a bit more severe
#> [1] "J01DD04"  "QJ01DD04"
ab_atc("seephthriaaksone") # and even this works
#> [1] "J01DD04"  "QJ01DD04"

# use ab_* functions to get a specific properties (see ?ab_property);
# they use as.ab() internally:
ab_name("J01FA01")
#> [1] "Erythromycin"
ab_name("eryt")
#> [1] "Erythromycin"

# \donttest{
if (require("dplyr")) {
  # you can quickly rename 'sir' columns using set_ab_names() with dplyr:
  example_isolates %>%
    set_ab_names(where(is.sir), property = "atc")
}
#> # A tibble: 2,000 × 46
#>    date       patient   age gender ward     mo           J01CE01 J01CF04 J01CF05
#>    <date>     <chr>   <dbl> <chr>  <chr>    <mo>         <sir>   <sir>   <sir>  
#>  1 2002-01-02 A77334     65 F      Clinical B_ESCHR_COLI   R       NA      NA   
#>  2 2002-01-03 A77334     65 F      Clinical B_ESCHR_COLI   R       NA      NA   
#>  3 2002-01-07 067927     45 F      ICU      B_STPHY_EPDR   R       NA      R    
#>  4 2002-01-07 067927     45 F      ICU      B_STPHY_EPDR   R       NA      R    
#>  5 2002-01-13 067927     45 F      ICU      B_STPHY_EPDR   R       NA      R    
#>  6 2002-01-13 067927     45 F      ICU      B_STPHY_EPDR   R       NA      R    
#>  7 2002-01-14 462729     78 M      Clinical B_STPHY_AURS   R       NA      S    
#>  8 2002-01-14 462729     78 M      Clinical B_STPHY_AURS   R       NA      S    
#>  9 2002-01-16 067927     45 F      ICU      B_STPHY_EPDR   R       NA      R    
#> 10 2002-01-17 858515     79 F      ICU      B_STPHY_EPDR   R       NA      S    
#> # ℹ 1,990 more rows
#> # ℹ 37 more variables: J01CA04 <sir>, J01CR02 <sir>, J01CA01 <sir>,
#> #   J01CR05 <sir>, J01DB04 <sir>, J01DE01 <sir>, J01DC02 <sir>, J01DC01 <sir>,
#> #   J01DD01 <sir>, J01DD02 <sir>, J01DD04 <sir>, J01GB03 <sir>, J01GB01 <sir>,
#> #   J01GB06 <sir>, J01GB04 <sir>, J01EA01 <sir>, J01EE01 <sir>, J01XE01 <sir>,
#> #   J01XX01 <sir>, J01XX08 <sir>, J01MA02 <sir>, J01MA14 <sir>, J01XA01 <sir>,
#> #   J01XA02 <sir>, J01AA07 <sir>, J01AA12 <sir>, J01AA02 <sir>, …
# }
```
