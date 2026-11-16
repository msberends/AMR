# Get Properties of an Antibiotic

Use these functions to return a specific property of an antibiotic from
the
[antimicrobials](https://amr-for-r.org/release/3.0.1/reference/antimicrobials.md)
data set. All input values will be evaluated internally with
[`as.ab()`](https://amr-for-r.org/release/3.0.1/reference/as.ab.md).

## Usage

``` r
ab_name(x, language = get_AMR_locale(), tolower = FALSE, ...)

ab_cid(x, ...)

ab_synonyms(x, ...)

ab_tradenames(x, ...)

ab_group(x, language = get_AMR_locale(), ...)

ab_atc(x, only_first = FALSE, ...)

ab_atc_group1(x, language = get_AMR_locale(), ...)

ab_atc_group2(x, language = get_AMR_locale(), ...)

ab_loinc(x, ...)

ab_ddd(x, administration = "oral", ...)

ab_ddd_units(x, administration = "oral", ...)

ab_info(x, language = get_AMR_locale(), ...)

ab_url(x, open = FALSE, ...)

ab_property(x, property = "name", language = get_AMR_locale(), ...)

set_ab_names(data, ..., property = "name", language = get_AMR_locale(),
  snake_case = NULL)
```

## Arguments

- x:

  Any (vector of) text that can be coerced to a valid antibiotic drug
  code with
  [`as.ab()`](https://amr-for-r.org/release/3.0.1/reference/as.ab.md).

- language:

  Language of the returned text - the default is the current system
  language (see
  [`get_AMR_locale()`](https://amr-for-r.org/release/3.0.1/reference/translate.md))
  and can also be set with the package option
  [`AMR_locale`](https://amr-for-r.org/release/3.0.1/reference/AMR-options.md).
  Use `language = NULL` or `language = ""` to prevent translation.

- tolower:

  A [logical](https://rdrr.io/r/base/logical.html) to indicate whether
  the first [character](https://rdrr.io/r/base/character.html) of every
  output should be transformed to a lower case
  [character](https://rdrr.io/r/base/character.html). This will lead to
  e.g. "polymyxin B" and not "polymyxin b".

- ...:

  In case of `set_ab_names()` and `data` is a
  [data.frame](https://rdrr.io/r/base/data.frame.html): columns to
  select (supports tidy selection such as `column1:column4`), otherwise
  other arguments passed on to
  [`as.ab()`](https://amr-for-r.org/release/3.0.1/reference/as.ab.md).

- only_first:

  A [logical](https://rdrr.io/r/base/logical.html) to indicate whether
  only the first ATC code must be returned, with giving preference to
  J0-codes (i.e., the antimicrobial drug group).

- administration:

  Way of administration, either `"oral"` or `"iv"`.

- open:

  Browse the URL using
  [`utils::browseURL()`](https://rdrr.io/r/utils/browseURL.html).

- property:

  One of the column names of one of the
  [antimicrobials](https://amr-for-r.org/release/3.0.1/reference/antimicrobials.md)
  data set: `vector_or(colnames(antimicrobials), sort = FALSE)`.

- data:

  A [data.frame](https://rdrr.io/r/base/data.frame.html) of which the
  columns need to be renamed, or a
  [character](https://rdrr.io/r/base/character.html) vector of column
  names.

- snake_case:

  A [logical](https://rdrr.io/r/base/logical.html) to indicate whether
  the names should be in so-called [snake
  case](https://en.wikipedia.org/wiki/Snake_case): in lower case and all
  spaces/slashes replaced with an underscore (`_`).

## Value

- An [integer](https://rdrr.io/r/base/integer.html) in case of
  `ab_cid()`

- A named [list](https://rdrr.io/r/base/list.html) in case of
  `ab_info()` and multiple `ab_atc()`/`ab_synonyms()`/`ab_tradenames()`

- A [double](https://rdrr.io/r/base/double.html) in case of `ab_ddd()`

- A [data.frame](https://rdrr.io/r/base/data.frame.html) in case of
  `set_ab_names()`

- A [character](https://rdrr.io/r/base/character.html) in all other
  cases

## Details

All output [will be
translated](https://amr-for-r.org/release/3.0.1/reference/translate.md)
where possible.

The function `ab_url()` will return the direct URL to the official WHO
website. A warning will be returned if the required ATC code is not
available.

The function `set_ab_names()` is a special column renaming function for
[data.frame](https://rdrr.io/r/base/data.frame.html)s. It renames
columns names that resemble antimicrobial drugs. It always makes sure
that the new column names are unique. If `property = "atc"` is set,
preference is given to ATC codes from the J-group.

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

[antimicrobials](https://amr-for-r.org/release/3.0.1/reference/antimicrobials.md)

## Examples

``` r
# all properties:
ab_name("AMX")
ab_atc("AMX")
ab_cid("AMX")
ab_synonyms("AMX")
ab_tradenames("AMX")
ab_group("AMX")
ab_atc_group1("AMX")
ab_atc_group2("AMX")
ab_url("AMX")

# smart lowercase transformation
ab_name(x = c("AMC", "PLB"))
ab_name(x = c("AMC", "PLB"), tolower = TRUE)

# defined daily doses (DDD)
ab_ddd("AMX", "oral")
ab_ddd_units("AMX", "oral")
ab_ddd("AMX", "iv")
ab_ddd_units("AMX", "iv")

ab_info("AMX") # all properties as a list

# all ab_* functions use as.ab() internally, so you can go from 'any' to 'any':
ab_atc("AMP")
ab_group("J01CA01")
ab_loinc("ampicillin")
ab_name("21066-6")
ab_name(6249)
ab_name("J01CA01")

# spelling from different languages and dyslexia are no problem
ab_atc("ceftriaxon")
ab_atc("cephtriaxone")
ab_atc("cephthriaxone")
ab_atc("seephthriaaksone")

# use set_ab_names() for renaming columns
colnames(example_isolates)
colnames(set_ab_names(example_isolates))
colnames(set_ab_names(example_isolates, NIT:VAN))
# \donttest{
if (require("dplyr")) {
  example_isolates %>%
    set_ab_names()

  # this does the same:
  example_isolates %>%
    rename_with(set_ab_names)

  # set_ab_names() works with any AB property:
  example_isolates %>%
    set_ab_names(property = "atc")

  example_isolates %>%
    set_ab_names(where(is.sir)) %>%
    colnames()

  example_isolates %>%
    set_ab_names(NIT:VAN) %>%
    colnames()
}
# }
```
