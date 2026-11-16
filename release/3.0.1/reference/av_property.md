# Get Properties of an Antiviral Drug

Use these functions to return a specific property of an antiviral drug
from the
[antivirals](https://amr-for-r.org/release/3.0.1/reference/antimicrobials.md)
data set. All input values will be evaluated internally with
[`as.av()`](https://amr-for-r.org/release/3.0.1/reference/as.av.md).

## Usage

``` r
av_name(x, language = get_AMR_locale(), tolower = FALSE, ...)

av_cid(x, ...)

av_synonyms(x, ...)

av_tradenames(x, ...)

av_group(x, language = get_AMR_locale(), ...)

av_atc(x, ...)

av_loinc(x, ...)

av_ddd(x, administration = "oral", ...)

av_ddd_units(x, administration = "oral", ...)

av_info(x, language = get_AMR_locale(), ...)

av_url(x, open = FALSE, ...)

av_property(x, property = "name", language = get_AMR_locale(), ...)
```

## Arguments

- x:

  Any (vector of) text that can be coerced to a valid antiviral drug
  code with
  [`as.av()`](https://amr-for-r.org/release/3.0.1/reference/as.av.md).

- language:

  Language of the returned text - the default is system language (see
  [`get_AMR_locale()`](https://amr-for-r.org/release/3.0.1/reference/translate.md))
  and can also be set with the package option
  [`AMR_locale`](https://amr-for-r.org/release/3.0.1/reference/AMR-options.md).
  Use `language = NULL` or `language = ""` to prevent translation.

- tolower:

  A [logical](https://rdrr.io/r/base/logical.html) to indicate whether
  the first [character](https://rdrr.io/r/base/character.html) of every
  output should be transformed to a lower case
  [character](https://rdrr.io/r/base/character.html).

- ...:

  Other arguments passed on to
  [`as.av()`](https://amr-for-r.org/release/3.0.1/reference/as.av.md).

- administration:

  Way of administration, either `"oral"` or `"iv"`.

- open:

  Browse the URL using
  [`utils::browseURL()`](https://rdrr.io/r/utils/browseURL.html).

- property:

  One of the column names of one of the
  [antivirals](https://amr-for-r.org/release/3.0.1/reference/antimicrobials.md)
  data set: `vector_or(colnames(antivirals), sort = FALSE)`.

## Value

- An [integer](https://rdrr.io/r/base/integer.html) in case of
  `av_cid()`

- A named [list](https://rdrr.io/r/base/list.html) in case of
  `av_info()` and multiple `av_atc()`/`av_synonyms()`/`av_tradenames()`

- A [double](https://rdrr.io/r/base/double.html) in case of `av_ddd()`

- A [character](https://rdrr.io/r/base/character.html) in all other
  cases

## Details

All output [will be
translated](https://amr-for-r.org/release/3.0.1/reference/translate.md)
where possible.

The function `av_url()` will return the direct URL to the official WHO
website. A warning will be returned if the required ATC code is not
available.

## Source

World Health Organization (WHO) Collaborating Centre for Drug Statistics
Methodology: <https://atcddd.fhi.no/atc_ddd_index/>

European Commission Public Health PHARMACEUTICALS - COMMUNITY REGISTER:
<https://ec.europa.eu/health/documents/community-register/html/reg_hum_atc.htm>

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

[antivirals](https://amr-for-r.org/release/3.0.1/reference/antimicrobials.md)

## Examples

``` r
# all properties:
av_name("ACI")
av_atc("ACI")
av_cid("ACI")
av_synonyms("ACI")
av_tradenames("ACI")
av_group("ACI")
av_url("ACI")

# lowercase transformation
av_name(x = c("ACI", "VALA"))
av_name(x = c("ACI", "VALA"), tolower = TRUE)

# defined daily doses (DDD)
av_ddd("ACI", "oral")
av_ddd_units("ACI", "oral")
av_ddd("ACI", "iv")
av_ddd_units("ACI", "iv")

av_info("ACI") # all properties as a list

# all av_* functions use as.av() internally, so you can go from 'any' to 'any':
av_atc("ACI")
av_group("J05AB01")
av_loinc("abacavir")
av_name("29113-8")
av_name(135398513)
av_name("J05AB01")
```
