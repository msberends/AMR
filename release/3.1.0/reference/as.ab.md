# Transform Input to an Antibiotic ID

Use this function to determine the antimicrobial drug code of one or
more antimicrobials. The data set
[antimicrobials](https://amr-for-r.org/release/3.1.0/reference/antimicrobials.md)
will be searched for abbreviations, official names and synonyms (brand
names).

## Usage

``` r
as.ab(
  x,
  flag_multiple_results = TRUE,
  language = get_AMR_locale(),
  info = interactive(),
  ...
)

is.ab(x)

ab_reset_session()

NA_ab_
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
  [`get_AMR_locale()`](https://amr-for-r.org/release/3.1.0/reference/translate.md)).

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
[antimicrobials](https://amr-for-r.org/release/3.1.0/reference/antimicrobials.md)
data set have three different identifiers: a human readable EARS-Net
code (column `ab`, used by ECDC and WHONET), an ATC code (column `atc`,
used by WHO), and a CID code (column `cid`, Compound ID, used by
PubChem). The data set contains more than 5,000 official brand names
from many different countries, as found in PubChem. Not that some drugs
contain multiple ATC codes.

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

Use the
[`ab_*`](https://amr-for-r.org/release/3.1.0/reference/ab_property.md)
functions to get properties based on the returned antibiotic ID, see
*Examples*.

Results of non-exact input are remembered during the session, so that
the same input is not matched again. Use `ab_reset_session()` to forget
these results. They are also forgotten after
[`add_custom_antimicrobials()`](https://amr-for-r.org/release/3.1.0/reference/add_custom_antimicrobials.md)
or
[`clear_custom_antimicrobials()`](https://amr-for-r.org/release/3.1.0/reference/add_custom_antimicrobials.md).

Note: the `as.ab()` and
[`ab_*`](https://amr-for-r.org/release/3.1.0/reference/ab_property.md)
functions may use very long regular expression to match brand names of
antimicrobial drugs. This may fail on some systems.

You can add your own manual codes to be considered by `as.ab()` and all
[`ab_*`](https://amr-for-r.org/release/3.1.0/reference/ab_property.md)
functions, see
[`add_custom_antimicrobials()`](https://amr-for-r.org/release/3.1.0/reference/add_custom_antimicrobials.md).

`NA_ab_` is a missing value of the new `ab` class, analogous to e.g.
base R's [`NA_character_`](https://rdrr.io/r/base/NA.html).

## Source

World Health Organization (WHO) Collaborating Centre for Drug Statistics
Methodology: <https://atcddd.fhi.no/atc_ddd_index/>

European Commission Public Health PHARMACEUTICALS - COMMUNITY REGISTER:
<https://ec.europa.eu/health/documents/community-register/html/index_en.htm>

## WHOCC

This package contains **all ~550 antibiotic, antimycotic and antiviral
drugs** and their Anatomical Therapeutic Chemical (ATC) codes, ATC
groups and Defined Daily Dose (DDD) from the World Health Organization
Collaborating Centre for Drug Statistics Methodology (WHOCC,
<https://atcddd.fhi.no>) and the Pharmaceuticals Community Register of
the European Commission
(<https://ec.europa.eu/health/documents/community-register/html/index_en.htm>).

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

- [antimicrobials](https://amr-for-r.org/release/3.1.0/reference/antimicrobials.md)
  for the [data.frame](https://rdrr.io/r/base/data.frame.html) that is
  being used to determine ATCs

- [`ab_from_text()`](https://amr-for-r.org/release/3.1.0/reference/ab_from_text.md)
  for a function to retrieve antimicrobial drugs from clinical text
  (from health care records)

## Examples

``` r
# these examples all return "ERY", the ID of erythromycin:
as.ab("J01FA01")
as.ab("J 01 FA 01")
as.ab("Erythromycin")
as.ab("eryt")
as.ab("ERYT")
as.ab("ERY")
as.ab("eritromicine") # spelled wrong, yet works
as.ab("Erythrocin") # trade name

# spelling from different languages and dyslexia are no problem
ab_atc("ceftriaxon")
ab_atc("cephtriaxone") # small spelling error
ab_atc("cephthriaxone") # or a bit more severe
ab_atc("seephthriaaksone") # and even this works

# use ab_* functions to get a specific properties (see ?ab_property);
# they use as.ab() internally:
ab_name("J01FA01")
ab_name("eryt")

# \donttest{
if (require("dplyr")) {
  # you can quickly rename 'sir' columns using set_ab_names() with dplyr:
  example_isolates %>%
    set_ab_names(where(is.sir), property = "atc")
}
# }
```
