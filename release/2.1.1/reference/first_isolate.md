# Determine First Isolates

Determine first isolates of all microorganisms of every patient per
episode and (if needed) per specimen type. These functions support all
four methods as summarised by Hindler *et al.* in 2007
([doi:10.1086/511864](https://doi.org/10.1086/511864) ). To determine
patient episodes not necessarily based on microorganisms, use
[`is_new_episode()`](https://amr-for-r.org/release/2.1.1/reference/get_episode.md)
that also supports grouping with the `dplyr` package.

## Usage

``` r
first_isolate(
  x = NULL,
  col_date = NULL,
  col_patient_id = NULL,
  col_mo = NULL,
  col_testcode = NULL,
  col_specimen = NULL,
  col_icu = NULL,
  col_keyantimicrobials = NULL,
  episode_days = 365,
  testcodes_exclude = NULL,
  icu_exclude = FALSE,
  specimen_group = NULL,
  type = "points",
  method = c("phenotype-based", "episode-based", "patient-based", "isolate-based"),
  ignore_I = TRUE,
  points_threshold = 2,
  info = interactive(),
  include_unknown = FALSE,
  include_untested_sir = TRUE,
  ...
)

filter_first_isolate(
  x = NULL,
  col_date = NULL,
  col_patient_id = NULL,
  col_mo = NULL,
  episode_days = 365,
  method = c("phenotype-based", "episode-based", "patient-based", "isolate-based"),
  ...
)
```

## Source

Methodology of this function is strictly based on:

- **M39 Analysis and Presentation of Cumulative Antimicrobial
  Susceptibility Test Data, 5th Edition**, 2022, *Clinical and
  Laboratory Standards Institute (CLSI)*.
  <https://clsi.org/standards/products/microbiology/documents/m39/>.

- Hindler JF and Stelling J (2007). **Analysis and Presentation of
  Cumulative Antibiograms: A New Consensus Guideline from the Clinical
  and Laboratory Standards Institute.** Clinical Infectious Diseases,
  44(6), 867-873. [doi:10.1086/511864](https://doi.org/10.1086/511864)

## Arguments

- x:

  a [data.frame](https://rdrr.io/r/base/data.frame.html) containing
  isolates. Can be left blank for automatic determination, see
  *Examples*.

- col_date:

  column name of the result date (or date that is was received on the
  lab) - the default is the first column with a date class

- col_patient_id:

  column name of the unique IDs of the patients - the default is the
  first column that starts with 'patient' or 'patid' (case insensitive)

- col_mo:

  column name of the names or codes of the microorganisms (see
  [`as.mo()`](https://amr-for-r.org/release/2.1.1/reference/as.mo.md)) -
  the default is the first column of class
  [`mo`](https://amr-for-r.org/release/2.1.1/reference/as.mo.md). Values
  will be coerced using
  [`as.mo()`](https://amr-for-r.org/release/2.1.1/reference/as.mo.md).

- col_testcode:

  column name of the test codes. Use `col_testcode = NULL` to **not**
  exclude certain test codes (such as test codes for screening). In that
  case `testcodes_exclude` will be ignored.

- col_specimen:

  column name of the specimen type or group

- col_icu:

  column name of the logicals (`TRUE`/`FALSE`) whether a ward or
  department is an Intensive Care Unit (ICU). This can also be a
  [logical](https://rdrr.io/r/base/logical.html) vector with the same
  length as rows in `x`.

- col_keyantimicrobials:

  (only useful when `method = "phenotype-based"`) column name of the key
  antimicrobials to determine first isolates, see
  [`key_antimicrobials()`](https://amr-for-r.org/release/2.1.1/reference/key_antimicrobials.md).
  The default is the first column that starts with 'key' followed by
  'ab' or 'antibiotics' or 'antimicrobials' (case insensitive). Use
  `col_keyantimicrobials = FALSE` to prevent this. Can also be the
  output of
  [`key_antimicrobials()`](https://amr-for-r.org/release/2.1.1/reference/key_antimicrobials.md).

- episode_days:

  episode in days after which a genus/species combination will be
  determined as 'first isolate' again. The default of 365 days is based
  on the guideline by CLSI, see *Source*.

- testcodes_exclude:

  a [character](https://rdrr.io/r/base/character.html) vector with test
  codes that should be excluded (case-insensitive)

- icu_exclude:

  a [logical](https://rdrr.io/r/base/logical.html) to indicate whether
  ICU isolates should be excluded (rows with value `TRUE` in the column
  set with `col_icu`)

- specimen_group:

  value in the column set with `col_specimen` to filter on

- type:

  type to determine weighed isolates; can be `"keyantimicrobials"` or
  `"points"`, see *Details*

- method:

  the method to apply, either `"phenotype-based"`, `"episode-based"`,
  `"patient-based"` or `"isolate-based"` (can be abbreviated), see
  *Details*. The default is `"phenotype-based"` if antimicrobial test
  results are present in the data, and `"episode-based"` otherwise.

- ignore_I:

  [logical](https://rdrr.io/r/base/logical.html) to indicate whether
  antibiotic interpretations with `"I"` will be ignored when
  `type = "keyantimicrobials"`, see *Details*

- points_threshold:

  minimum number of points to require before differences in the
  antibiogram will lead to inclusion of an isolate when
  `type = "points"`, see *Details*

- info:

  a [logical](https://rdrr.io/r/base/logical.html) to indicate info
  should be printed - the default is `TRUE` only in interactive mode

- include_unknown:

  a [logical](https://rdrr.io/r/base/logical.html) to indicate whether
  'unknown' microorganisms should be included too, i.e. microbial code
  `"UNKNOWN"`, which defaults to `FALSE`. For WHONET users, this means
  that all records with organism code `"con"` (*contamination*) will be
  excluded at default. Isolates with a microbial ID of `NA` will always
  be excluded as first isolate.

- include_untested_sir:

  a [logical](https://rdrr.io/r/base/logical.html) to indicate whether
  also rows without antibiotic results are still eligible for becoming a
  first isolate. Use `include_untested_sir = FALSE` to always return
  `FALSE` for such rows. This checks the data set for columns of class
  `sir` and consequently requires transforming columns with antibiotic
  results using
  [`as.sir()`](https://amr-for-r.org/release/2.1.1/reference/as.sir.md)
  first.

- ...:

  arguments passed on to `first_isolate()` when using
  `filter_first_isolate()`, otherwise arguments passed on to
  [`key_antimicrobials()`](https://amr-for-r.org/release/2.1.1/reference/key_antimicrobials.md)
  (such as `universal`, `gram_negative`, `gram_positive`)

## Value

A [logical](https://rdrr.io/r/base/logical.html) vector

## Details

To conduct epidemiological analyses on antimicrobial resistance data,
only so-called first isolates should be included to prevent
overestimation and underestimation of antimicrobial resistance.
Different methods can be used to do so, see below.

These functions are context-aware. This means that the `x` argument can
be left blank if used inside a
[data.frame](https://rdrr.io/r/base/data.frame.html) call, see
*Examples*.

The `first_isolate()` function is a wrapper around the
[`is_new_episode()`](https://amr-for-r.org/release/2.1.1/reference/get_episode.md)
function, but more efficient for data sets containing microorganism
codes or names.

All isolates with a microbial ID of `NA` will be excluded as first
isolate.

### Different methods

According to Hindler *et al.* (2007,
[doi:10.1086/511864](https://doi.org/10.1086/511864) ), there are
different methods (algorithms) to select first isolates with increasing
reliability: isolate-based, patient-based, episode-based and
phenotype-based. All methods select on a combination of the taxonomic
genus and species (not subspecies).

All mentioned methods are covered in the `first_isolate()` function:

|  |  |
|----|----|
| **Method** | **Function to apply** |
| **Isolate-based** | `first_isolate(x, method = "isolate-based")` |
| *(= all isolates)* |  |
|  |  |
|  |  |
| **Patient-based** | `first_isolate(x, method = "patient-based")` |
| *(= first isolate per patient)* |  |
|  |  |
|  |  |
| **Episode-based** | `first_isolate(x, method = "episode-based")`, or: |
| *(= first isolate per episode)* |  |
| \- 7-Day interval from initial isolate | \- `first_isolate(x, method = "e", episode_days = 7)` |
| \- 30-Day interval from initial isolate | \- `first_isolate(x, method = "e", episode_days = 30)` |
|  |  |
|  |  |
| **Phenotype-based** | `first_isolate(x, method = "phenotype-based")`, or: |
| *(= first isolate per phenotype)* |  |
| \- Major difference in any antimicrobial result | \- `first_isolate(x, type = "points")` |
| \- Any difference in key antimicrobial results | \- `first_isolate(x, type = "keyantimicrobials")` |

### Isolate-based

This method does not require any selection, as all isolates should be
included. It does, however, respect all arguments set in the
`first_isolate()` function. For example, the default setting for
`include_unknown` (`FALSE`) will omit selection of rows without a
microbial ID.

### Patient-based

To include every genus-species combination per patient once, set the
`episode_days` to `Inf`. Although often inappropriate, this method makes
sure that no duplicate isolates are selected from the same patient. In a
large longitudinal data set, this could mean that isolates are
*excluded* that were found years after the initial isolate.

### Episode-based

To include every genus-species combination per patient episode once, set
the `episode_days` to a sensible number of days. Depending on the type
of analysis, this could be 14, 30, 60 or 365. Short episodes are common
for analysing specific hospital or ward data, long episodes are common
for analysing regional and national data.

This is the most common method to correct for duplicate isolates.
Patients are categorised into episodes based on their ID and dates
(e.g., the date of specimen receipt or laboratory result). While this is
a common method, it does not take into account antimicrobial test
results. This means that e.g. a methicillin-resistant *Staphylococcus
aureus* (MRSA) isolate cannot be differentiated from a wildtype
*Staphylococcus aureus* isolate.

### Phenotype-based

This is a more reliable method, since it also *weighs* the antibiogram
(antimicrobial test results) yielding so-called 'first weighted
isolates'. There are two different methods to weigh the antibiogram:

1.  Using `type = "points"` and argument `points_threshold` (default)

    This method weighs *all* antimicrobial drugs available in the data
    set. Any difference from I to S or R (or vice versa) counts as `0.5`
    points, a difference from S to R (or vice versa) counts as `1`
    point. When the sum of points exceeds `points_threshold`, which
    defaults to `2`, an isolate will be selected as a first weighted
    isolate.

    All antimicrobials are internally selected using the
    [`all_antimicrobials()`](https://amr-for-r.org/release/2.1.1/reference/key_antimicrobials.md)
    function. The output of this function does not need to be passed to
    the `first_isolate()` function.

2.  Using `type = "keyantimicrobials"` and argument `ignore_I`

    This method only weighs specific antimicrobial drugs, called *key
    antimicrobials*. Any difference from S to R (or vice versa) in these
    key antimicrobials will select an isolate as a first weighted
    isolate. With `ignore_I = FALSE`, also differences from I to S or R
    (or vice versa) will lead to this.

    Key antimicrobials are internally selected using the
    [`key_antimicrobials()`](https://amr-for-r.org/release/2.1.1/reference/key_antimicrobials.md)
    function, but can also be added manually as a variable to the data
    and set in the `col_keyantimicrobials` argument. Another option is
    to pass the output of the
    [`key_antimicrobials()`](https://amr-for-r.org/release/2.1.1/reference/key_antimicrobials.md)
    function directly to the `col_keyantimicrobials` argument.

The default method is phenotype-based (using `type = "points"`) and
episode-based (using `episode_days = 365`). This makes sure that every
genus-species combination is selected per patient once per year, while
taking into account all antimicrobial test results. If no antimicrobial
test results are available in the data set, only the episode-based
method is applied at default.

## See also

[`key_antimicrobials()`](https://amr-for-r.org/release/2.1.1/reference/key_antimicrobials.md)

## Examples

``` r
# `example_isolates` is a data set available in the AMR package.
# See ?example_isolates.

example_isolates[first_isolate(info = TRUE), ]
# \donttest{
# get all first Gram-negatives
example_isolates[which(first_isolate(info = FALSE) & mo_is_gram_negative()), ]

if (require("dplyr")) {
  # filter on first isolates using dplyr:
  example_isolates %>%
    filter(first_isolate(info = TRUE))
}
if (require("dplyr")) {
  # short-hand version:
  example_isolates %>%
    filter_first_isolate(info = FALSE)
}
if (require("dplyr")) {
  # flag the first isolates per group:
  example_isolates %>%
    group_by(ward) %>%
    mutate(first = first_isolate(info = TRUE)) %>%
    select(ward, date, patient, mo, first)
}
# }
```
