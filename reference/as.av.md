# Transform Input to an Antiviral Drug ID

Use this function to determine the antiviral drug code of one or more
antiviral drugs. The data set
[antivirals](https://amr-for-r.org/reference/antimicrobials.md) will be
searched for abbreviations, official names and synonyms (brand names).

## Usage

``` r
as.av(x, flag_multiple_results = TRUE, info = interactive(), ...)

is.av(x)
```

## Arguments

- x:

  A [character](https://rdrr.io/r/base/character.html) vector to
  determine to antiviral drug ID.

- flag_multiple_results:

  A [logical](https://rdrr.io/r/base/logical.html) to indicate whether a
  note should be printed to the console that probably more than one
  antiviral drug code or name can be retrieved from a single input
  value.

- info:

  A [logical](https://rdrr.io/r/base/logical.html) to indicate whether a
  progress bar should be printed - the default is `TRUE` only in
  interactive mode.

- ...:

  Arguments passed on to internal functions.

## Value

A [character](https://rdrr.io/r/base/character.html)
[vector](https://rdrr.io/r/base/vector.html) with additional class
[`ab`](https://amr-for-r.org/reference/as.ab.md)

## Details

All entries in the
[antivirals](https://amr-for-r.org/reference/antimicrobials.md) data set
have three different identifiers: a human readable EARS-Net code (column
`ab`, used by ECDC and WHONET), an ATC code (column `atc`, used by WHO),
and a CID code (column `cid`, Compound ID, used by PubChem). The data
set contains more than 5,000 official brand names from many different
countries, as found in PubChem. Not that some drugs contain multiple ATC
codes.

All these properties will be searched for the user input. The `as.av()`
can correct for different forms of misspelling:

- Wrong spelling of drug names (such as "acyclovir"), which corrects for
  most audible similarities such as f/ph, x/ks, c/z/s, t/th, etc.

- Too few or too many vowels or consonants

- Switching two characters (such as "aycclovir", often the case in
  clinical data, when doctors typed too fast)

- Digitalised paper records, leaving artefacts like 0/o/O (zero and
  O's), B/8, n/r, etc.

Use the [`av_*`](https://amr-for-r.org/reference/av_property.md)
functions to get properties based on the returned antiviral drug ID, see
*Examples*.

Note: the `as.av()` and
[`av_*`](https://amr-for-r.org/reference/av_property.md) functions may
use very long regular expression to match brand names of antimicrobial
drugs. This may fail on some systems.

## Source

World Health Organization (WHO) Collaborating Centre for Drug Statistics
Methodology: <https://atcddd.fhi.no/atc_ddd_index/>

European Commission Public Health PHARMACEUTICALS - COMMUNITY REGISTER:
<https://ec.europa.eu/health/documents/community-register/html/reg_hum_atc.htm>

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

- [antivirals](https://amr-for-r.org/reference/antimicrobials.md) for
  the [data.frame](https://rdrr.io/r/base/data.frame.html) that is being
  used to determine ATCs

- [`av_from_text()`](https://amr-for-r.org/reference/av_from_text.md)
  for a function to retrieve antimicrobial drugs from clinical text
  (from health care records)

## Examples

``` r
# these examples all return "ACI", the ID of aciclovir:
as.av("J05AB01")
#> Class 'av'
#> [1] ACI
as.av("J 05 AB 01")
#> Class 'av'
#> [1] ACI
as.av("Aciclovir")
#> Class 'av'
#> [1] ACI
as.av("aciclo")
#> Class 'av'
#> [1] ACI
as.av("   aciclo 123")
#> Class 'av'
#> [1] ACI
as.av("ACICL")
#> Class 'av'
#> [1] ACI
as.av("ACI")
#> Class 'av'
#> [1] ACI
as.av("Virorax") # trade name
#> Class 'av'
#> [1] ACI
as.av("Zovirax") # trade name
#> Class 'av'
#> [1] ACI

as.av("acyklofir") # severe spelling error, yet works
#> Class 'av'
#> [1] ACI

# use av_* functions to get a specific properties (see ?av_property);
# they use as.av() internally:
av_name("J05AB01")
#> [1] "Aciclovir"
av_name("acicl")
#> [1] "Aciclovir"
```
