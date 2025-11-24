# Determine First Isolates

Determine first isolates of all microorganisms of every patient per
episode and (if needed) per specimen type. These functions support all
four methods as summarised by Hindler *et al.* in 2007
([doi:10.1086/511864](https://doi.org/10.1086/511864) ). To determine
patient episodes not necessarily based on microorganisms, use
[`is_new_episode()`](https://amr-for-r.org/reference/get_episode.md)
that also supports grouping with the `dplyr` package.

## Usage

``` r
first_isolate(x = NULL, col_date = NULL, col_patient_id = NULL,
  col_mo = NULL, col_testcode = NULL, col_specimen = NULL,
  col_icu = NULL, col_keyantimicrobials = NULL, episode_days = 365,
  testcodes_exclude = NULL, icu_exclude = FALSE, specimen_group = NULL,
  type = "points", method = c("phenotype-based", "episode-based",
  "patient-based", "isolate-based"), ignore_I = TRUE, points_threshold = 2,
  info = interactive(), include_unknown = FALSE,
  include_untested_sir = TRUE, ...)

filter_first_isolate(x = NULL, col_date = NULL, col_patient_id = NULL,
  col_mo = NULL, episode_days = 365, method = c("phenotype-based",
  "episode-based", "patient-based", "isolate-based"), ...)
```

## Source

Methodology of these functions is strictly based on:

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

  A [data.frame](https://rdrr.io/r/base/data.frame.html) containing
  isolates. Can be left blank for automatic determination, see
  *Examples*.

- col_date:

  Column name of the result date (or date that is was received on the
  lab) - the default is the first column with a date class.

- col_patient_id:

  Column name of the unique IDs of the patients - the default is the
  first column that starts with 'patient' or 'patid' (case insensitive).

- col_mo:

  Column name of the names or codes of the microorganisms (see
  [`as.mo()`](https://amr-for-r.org/reference/as.mo.md)) - the default
  is the first column of class
  [`mo`](https://amr-for-r.org/reference/as.mo.md). Values will be
  coerced using [`as.mo()`](https://amr-for-r.org/reference/as.mo.md).

- col_testcode:

  Column name of the test codes. Use `col_testcode = NULL` to **not**
  exclude certain test codes (such as test codes for screening). In that
  case `testcodes_exclude` will be ignored.

- col_specimen:

  Column name of the specimen type or group.

- col_icu:

  Column name of the logicals (`TRUE`/`FALSE`) whether a ward or
  department is an Intensive Care Unit (ICU). This can also be a
  [logical](https://rdrr.io/r/base/logical.html) vector with the same
  length as rows in `x`.

- col_keyantimicrobials:

  (only useful when `method = "phenotype-based"`) column name of the key
  antimicrobials to determine first isolates, see
  [`key_antimicrobials()`](https://amr-for-r.org/reference/key_antimicrobials.md).
  The default is the first column that starts with 'key' followed by
  'ab' or 'antibiotics' or 'antimicrobials' (case insensitive). Use
  `col_keyantimicrobials = FALSE` to prevent this. Can also be the
  output of
  [`key_antimicrobials()`](https://amr-for-r.org/reference/key_antimicrobials.md).

- episode_days:

  Episode in days after which a genus/species combination will be
  determined as 'first isolate' again. The default of 365 days is based
  on the guideline by CLSI, see *Source*.

- testcodes_exclude:

  A [character](https://rdrr.io/r/base/character.html) vector with test
  codes that should be excluded (case-insensitive).

- icu_exclude:

  A [logical](https://rdrr.io/r/base/logical.html) to indicate whether
  ICU isolates should be excluded (rows with value `TRUE` in the column
  set with `col_icu`).

- specimen_group:

  Value in the column set with `col_specimen` to filter on.

- type:

  Type to determine weighed isolates; can be `"keyantimicrobials"` or
  `"points"`, see *Details*.

- method:

  The method to apply, either `"phenotype-based"`, `"episode-based"`,
  `"patient-based"` or `"isolate-based"` (can be abbreviated), see
  *Details*. The default is `"phenotype-based"` if antimicrobial test
  results are present in the data, and `"episode-based"` otherwise.

- ignore_I:

  [logical](https://rdrr.io/r/base/logical.html) to indicate whether
  antibiotic interpretations with `"I"` will be ignored when
  `type = "keyantimicrobials"`, see *Details*.

- points_threshold:

  Minimum number of points to require before differences in the
  antibiogram will lead to inclusion of an isolate when
  `type = "points"`, see *Details*.

- info:

  A [logical](https://rdrr.io/r/base/logical.html) to indicate info
  should be printed - the default is `TRUE` only in interactive mode.

- include_unknown:

  A [logical](https://rdrr.io/r/base/logical.html) to indicate whether
  'unknown' microorganisms should be included too, i.e. microbial code
  `"UNKNOWN"`, which defaults to `FALSE`. For WHONET users, this means
  that all records with organism code `"con"` (*contamination*) will be
  excluded at default. Isolates with a microbial ID of `NA` will always
  be excluded as first isolate.

- include_untested_sir:

  A [logical](https://rdrr.io/r/base/logical.html) to indicate whether
  also rows without antibiotic results are still eligible for becoming a
  first isolate. Use `include_untested_sir = FALSE` to always return
  `FALSE` for such rows. This checks the data set for columns of class
  `sir` and consequently requires transforming columns with antibiotic
  results using [`as.sir()`](https://amr-for-r.org/reference/as.sir.md)
  first.

- ...:

  Arguments passed on to `first_isolate()` when using
  `filter_first_isolate()`, otherwise arguments passed on to
  [`key_antimicrobials()`](https://amr-for-r.org/reference/key_antimicrobials.md)
  (such as `universal`, `gram_negative`, `gram_positive`).

## Value

A [logical](https://rdrr.io/r/base/logical.html) vector

## Details

The methodology implemented in these functions is strictly based on the
recommendations outlined in [CLSI Guideline
M39](https://clsi.org/standards/products/microbiology/documents/m39) and
the research overview by Hindler *et al.* (2007,
[doi:10.1086/511864](https://doi.org/10.1086/511864) ).

To conduct epidemiological analyses on antimicrobial resistance data,
only so-called first isolates should be included to prevent
overestimation and underestimation of antimicrobial resistance.
Different methods can be used to do so, see below.

These functions are context-aware. This means that the `x` argument can
be left blank if used inside a
[data.frame](https://rdrr.io/r/base/data.frame.html) call, see
*Examples*.

The `first_isolate()` function is a wrapper around the
[`is_new_episode()`](https://amr-for-r.org/reference/get_episode.md)
function, but more efficient for data sets containing microorganism
codes or names.

All isolates with a microbial ID of `NA` will be excluded as first
isolate.

### Different methods

According to previously-mentioned sources, there are different methods
(algorithms) to select first isolates with increasing reliability:
isolate-based, patient-based, episode-based and phenotype-based. All
methods select on a combination of the taxonomic genus and species (not
subspecies).

All mentioned methods are covered in the `first_isolate()` function:

|                                                 |                                                        |
|-------------------------------------------------|--------------------------------------------------------|
| **Method**                                      | **Function to apply**                                  |
| **Isolate-based**                               | `first_isolate(x, method = "isolate-based")`           |
| *(= all isolates)*                              |                                                        |
|                                                 |                                                        |
|                                                 |                                                        |
| **Patient-based**                               | `first_isolate(x, method = "patient-based")`           |
| *(= first isolate per patient)*                 |                                                        |
|                                                 |                                                        |
|                                                 |                                                        |
| **Episode-based**                               | `first_isolate(x, method = "episode-based")`, or:      |
| *(= first isolate per episode)*                 |                                                        |
| \- 7-Day interval from initial isolate          | \- `first_isolate(x, method = "e", episode_days = 7)`  |
| \- 30-Day interval from initial isolate         | \- `first_isolate(x, method = "e", episode_days = 30)` |
|                                                 |                                                        |
|                                                 |                                                        |
| **Phenotype-based**                             | `first_isolate(x, method = "phenotype-based")`, or:    |
| *(= first isolate per phenotype)*               |                                                        |
| \- Major difference in any antimicrobial result | \- `first_isolate(x, type = "points")`                 |
| \- Any difference in key antimicrobial results  | \- `first_isolate(x, type = "keyantimicrobials")`      |

**Isolate-based**

*Minimum variables required: Microorganism identifier*

This method does not require any selection, as all isolates should be
included. It does, however, respect all arguments set in the
`first_isolate()` function. For example, the default setting for
`include_unknown` (`FALSE`) will omit selection of rows without a
microbial ID.

**Patient-based**

*Minimum variables required: Microorganism identifier, Patient
identifier*

This method includes every genus-species combination per patient once.
This method makes sure that no duplicate isolates are selected from the
same patient. This method is preferred to e.g. identify the first MRSA
finding of each patient to determine the incidence. Conversely, in a
large longitudinal data set, this could mean that isolates are
*excluded* that were found years after the initial isolate.

**Episode-based**

*Minimum variables required: Microorganism identifier, Patient
identifier, Date*

To include every genus-species combination per patient episode once, set
the `episode_days` to a sensible number of days. Depending on the type
of analysis, this could be e.g., 14, 30, 60 or 365. Short episodes are
common for analysing specific hospital or ward data or ICU cases, long
episodes are common for analysing regional and national data.

This is the most common method to correct for duplicate isolates.
Patients are categorised into episodes based on their ID and dates
(e.g., the date of specimen receipt or laboratory result). While this is
a common method, it does not take into account antimicrobial test
results. This means that e.g. a methicillin-resistant *Staphylococcus
aureus* (MRSA) isolate cannot be differentiated from a wildtype
*Staphylococcus aureus* isolate.

**Phenotype-based**

*Minimum variables required: Microorganism identifier, Patient
identifier, Date, Antimicrobial test results*

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
    [`all_antimicrobials()`](https://amr-for-r.org/reference/key_antimicrobials.md)
    function. The output of this function does not need to be passed to
    the `first_isolate()` function.

2.  Using `type = "keyantimicrobials"` and argument `ignore_I`

    This method only weighs specific antimicrobial drugs, called *key
    antimicrobials*. Any difference from S to R (or vice versa) in these
    key antimicrobials will select an isolate as a first weighted
    isolate. With `ignore_I = FALSE`, also differences from I to S or R
    (or vice versa) will lead to this.

    Key antimicrobials are internally selected using the
    [`key_antimicrobials()`](https://amr-for-r.org/reference/key_antimicrobials.md)
    function, but can also be added manually as a variable to the data
    and set in the `col_keyantimicrobials` argument. Another option is
    to pass the output of the
    [`key_antimicrobials()`](https://amr-for-r.org/reference/key_antimicrobials.md)
    function directly to the `col_keyantimicrobials` argument.

The default method is phenotype-based (using `type = "points"`) and
episode-based (using `episode_days = 365`). This makes sure that every
genus-species combination is selected per patient once per year, while
taking into account all antimicrobial test results. If no antimicrobial
test results are available in the data set, only the episode-based
method is applied at default.

## See also

[`key_antimicrobials()`](https://amr-for-r.org/reference/key_antimicrobials.md)

## Examples

``` r
# `example_isolates` is a data set available in the AMR package.
# See ?example_isolates.

example_isolates[first_isolate(info = TRUE), ]
#> ℹ Determining first isolates using an episode length of 365 days
#> ℹ Using column 'date' as input for `col_date`.
#> ℹ Using column 'patient' as input for `col_patient_id`.
#> ℹ Basing inclusion on all antimicrobial results, using a points threshold
#>   of 2
#> ℹ Excluding 16 isolates with a microbial ID 'UNKNOWN' (in column 'mo')
#> => Found 1,387 'phenotype-based' first isolates (69.4% of total where a
#>    microbial ID was available)
#> # A tibble: 1,387 × 46
#>    date       patient   age gender ward     mo           PEN   OXA   FLC   AMX  
#>    <date>     <chr>   <dbl> <chr>  <chr>    <mo>         <sir> <sir> <sir> <sir>
#>  1 2002-01-02 A77334     65 F      Clinical B_ESCHR_COLI   R     NA    NA    NA 
#>  2 2002-01-07 067927     45 F      ICU      B_STPHY_EPDR   R     NA    R     NA 
#>  3 2002-01-14 462729     78 M      Clinical B_STPHY_AURS   R     NA    S     R  
#>  4 2002-01-16 067927     45 F      ICU      B_STPHY_EPDR   R     NA    R     NA 
#>  5 2002-01-17 858515     79 F      ICU      B_STPHY_EPDR   R     NA    S     NA 
#>  6 2002-01-17 495616     67 M      Clinical B_STPHY_EPDR   R     NA    S     NA 
#>  7 2002-01-19 738003     71 M      Clinical B_ESCHR_COLI   R     NA    NA    NA 
#>  8 2002-01-21 462081     75 F      Clinical B_CTRBC_FRND   R     NA    NA    R  
#>  9 2002-01-22 F35553     50 M      ICU      B_PROTS_MRBL   R     NA    NA    NA 
#> 10 2002-02-03 481442     76 M      ICU      B_STPHY_CONS   R     NA    S     NA 
#> # ℹ 1,377 more rows
#> # ℹ 36 more variables: AMC <sir>, AMP <sir>, TZP <sir>, CZO <sir>, FEP <sir>,
#> #   CXM <sir>, FOX <sir>, CTX <sir>, CAZ <sir>, CRO <sir>, GEN <sir>,
#> #   TOB <sir>, AMK <sir>, KAN <sir>, TMP <sir>, SXT <sir>, NIT <sir>,
#> #   FOS <sir>, LNZ <sir>, CIP <sir>, MFX <sir>, VAN <sir>, TEC <sir>,
#> #   TCY <sir>, TGC <sir>, DOX <sir>, ERY <sir>, CLI <sir>, AZM <sir>,
#> #   IPM <sir>, MEM <sir>, MTR <sir>, CHL <sir>, COL <sir>, MUP <sir>, …
# \donttest{
# get all first Gram-negatives
example_isolates[which(first_isolate(info = FALSE) & mo_is_gram_negative()), ]
#> ℹ Using column 'mo' as input for `mo_is_gram_negative()`
#> # A tibble: 441 × 46
#>    date       patient   age gender ward     mo           PEN   OXA   FLC   AMX  
#>    <date>     <chr>   <dbl> <chr>  <chr>    <mo>         <sir> <sir> <sir> <sir>
#>  1 2002-01-02 A77334     65 F      Clinical B_ESCHR_COLI   R     NA    NA    NA 
#>  2 2002-01-19 738003     71 M      Clinical B_ESCHR_COLI   R     NA    NA    NA 
#>  3 2002-01-21 462081     75 F      Clinical B_CTRBC_FRND   R     NA    NA    R  
#>  4 2002-01-22 F35553     50 M      ICU      B_PROTS_MRBL   R     NA    NA    NA 
#>  5 2002-02-05 067927     45 F      ICU      B_SERRT_MRCS   R     NA    NA    R  
#>  6 2002-02-27 066895     85 F      Clinical B_KLBSL_PNMN   R     NA    NA    R  
#>  7 2002-03-08 4FC193     69 M      Clinical B_ESCHR_COLI   R     NA    NA    R  
#>  8 2002-03-16 4FC193     69 M      Clinical B_PSDMN_AERG   R     NA    NA    R  
#>  9 2002-04-01 496896     46 F      ICU      B_ESCHR_COLI   R     NA    NA    NA 
#> 10 2002-04-23 EE2510     69 F      ICU      B_ESCHR_COLI   R     NA    NA    NA 
#> # ℹ 431 more rows
#> # ℹ 36 more variables: AMC <sir>, AMP <sir>, TZP <sir>, CZO <sir>, FEP <sir>,
#> #   CXM <sir>, FOX <sir>, CTX <sir>, CAZ <sir>, CRO <sir>, GEN <sir>,
#> #   TOB <sir>, AMK <sir>, KAN <sir>, TMP <sir>, SXT <sir>, NIT <sir>,
#> #   FOS <sir>, LNZ <sir>, CIP <sir>, MFX <sir>, VAN <sir>, TEC <sir>,
#> #   TCY <sir>, TGC <sir>, DOX <sir>, ERY <sir>, CLI <sir>, AZM <sir>,
#> #   IPM <sir>, MEM <sir>, MTR <sir>, CHL <sir>, COL <sir>, MUP <sir>, …

if (require("dplyr")) {
  # filter on first isolates using dplyr:
  example_isolates %>%
    filter(first_isolate(info = TRUE))
}
#> ℹ Determining first isolates using an episode length of 365 days
#> ℹ Basing inclusion on all antimicrobial results, using a points threshold
#>   of 2
#> ℹ Excluding 16 isolates with a microbial ID 'UNKNOWN' (in column 'mo')
#> => Found 1,387 'phenotype-based' first isolates (69.4% of total where a
#>    microbial ID was available)
#> # A tibble: 1,387 × 46
#>    date       patient   age gender ward     mo           PEN   OXA   FLC   AMX  
#>    <date>     <chr>   <dbl> <chr>  <chr>    <mo>         <sir> <sir> <sir> <sir>
#>  1 2002-01-02 A77334     65 F      Clinical B_ESCHR_COLI   R     NA    NA    NA 
#>  2 2002-01-07 067927     45 F      ICU      B_STPHY_EPDR   R     NA    R     NA 
#>  3 2002-01-14 462729     78 M      Clinical B_STPHY_AURS   R     NA    S     R  
#>  4 2002-01-16 067927     45 F      ICU      B_STPHY_EPDR   R     NA    R     NA 
#>  5 2002-01-17 858515     79 F      ICU      B_STPHY_EPDR   R     NA    S     NA 
#>  6 2002-01-17 495616     67 M      Clinical B_STPHY_EPDR   R     NA    S     NA 
#>  7 2002-01-19 738003     71 M      Clinical B_ESCHR_COLI   R     NA    NA    NA 
#>  8 2002-01-21 462081     75 F      Clinical B_CTRBC_FRND   R     NA    NA    R  
#>  9 2002-01-22 F35553     50 M      ICU      B_PROTS_MRBL   R     NA    NA    NA 
#> 10 2002-02-03 481442     76 M      ICU      B_STPHY_CONS   R     NA    S     NA 
#> # ℹ 1,377 more rows
#> # ℹ 36 more variables: AMC <sir>, AMP <sir>, TZP <sir>, CZO <sir>, FEP <sir>,
#> #   CXM <sir>, FOX <sir>, CTX <sir>, CAZ <sir>, CRO <sir>, GEN <sir>,
#> #   TOB <sir>, AMK <sir>, KAN <sir>, TMP <sir>, SXT <sir>, NIT <sir>,
#> #   FOS <sir>, LNZ <sir>, CIP <sir>, MFX <sir>, VAN <sir>, TEC <sir>,
#> #   TCY <sir>, TGC <sir>, DOX <sir>, ERY <sir>, CLI <sir>, AZM <sir>,
#> #   IPM <sir>, MEM <sir>, MTR <sir>, CHL <sir>, COL <sir>, MUP <sir>, …
if (require("dplyr")) {
  # short-hand version:
  example_isolates %>%
    filter_first_isolate(info = FALSE)
}
#> # A tibble: 1,387 × 46
#>    date       patient   age gender ward     mo           PEN   OXA   FLC   AMX  
#>    <date>     <chr>   <dbl> <chr>  <chr>    <mo>         <sir> <sir> <sir> <sir>
#>  1 2002-01-02 A77334     65 F      Clinical B_ESCHR_COLI   R     NA    NA    NA 
#>  2 2002-01-07 067927     45 F      ICU      B_STPHY_EPDR   R     NA    R     NA 
#>  3 2002-01-14 462729     78 M      Clinical B_STPHY_AURS   R     NA    S     R  
#>  4 2002-01-16 067927     45 F      ICU      B_STPHY_EPDR   R     NA    R     NA 
#>  5 2002-01-17 858515     79 F      ICU      B_STPHY_EPDR   R     NA    S     NA 
#>  6 2002-01-17 495616     67 M      Clinical B_STPHY_EPDR   R     NA    S     NA 
#>  7 2002-01-19 738003     71 M      Clinical B_ESCHR_COLI   R     NA    NA    NA 
#>  8 2002-01-21 462081     75 F      Clinical B_CTRBC_FRND   R     NA    NA    R  
#>  9 2002-01-22 F35553     50 M      ICU      B_PROTS_MRBL   R     NA    NA    NA 
#> 10 2002-02-03 481442     76 M      ICU      B_STPHY_CONS   R     NA    S     NA 
#> # ℹ 1,377 more rows
#> # ℹ 36 more variables: AMC <sir>, AMP <sir>, TZP <sir>, CZO <sir>, FEP <sir>,
#> #   CXM <sir>, FOX <sir>, CTX <sir>, CAZ <sir>, CRO <sir>, GEN <sir>,
#> #   TOB <sir>, AMK <sir>, KAN <sir>, TMP <sir>, SXT <sir>, NIT <sir>,
#> #   FOS <sir>, LNZ <sir>, CIP <sir>, MFX <sir>, VAN <sir>, TEC <sir>,
#> #   TCY <sir>, TGC <sir>, DOX <sir>, ERY <sir>, CLI <sir>, AZM <sir>,
#> #   IPM <sir>, MEM <sir>, MTR <sir>, CHL <sir>, COL <sir>, MUP <sir>, …
if (require("dplyr")) {
  # flag the first isolates per group:
  example_isolates %>%
    group_by(ward) %>%
    mutate(first = first_isolate(info = TRUE)) %>%
    select(ward, date, patient, mo, first)
}
#> ℹ Determining first isolates using an episode length of 365 days
#> ℹ Basing inclusion on all antimicrobial results, using a points threshold
#>   of 2
#> 
#> Group: ward = "Clinical"
#> ℹ Excluding 9 isolates with a microbial ID 'UNKNOWN' (in column 'mo')
#> => Found 865 'phenotype-based' first isolates (70.1% of total where a
#>    microbial ID was available)
#> 
#> Group: ward = "ICU"
#> ℹ Excluding 6 isolates with a microbial ID 'UNKNOWN' (in column 'mo')
#> => Found 452 'phenotype-based' first isolates (70.0% of total where a
#>    microbial ID was available)
#> 
#> Group: ward = "Outpatient"
#> ℹ Excluding 1 isolates with a microbial ID 'UNKNOWN' (in column 'mo')
#> => Found 99 'phenotype-based' first isolates (82.5% of total where a
#>    microbial ID was available)
#> # A tibble: 2,000 × 5
#> # Groups:   ward [3]
#>    ward     date       patient mo           first
#>    <chr>    <date>     <chr>   <mo>         <lgl>
#>  1 Clinical 2002-01-02 A77334  B_ESCHR_COLI TRUE 
#>  2 Clinical 2002-01-03 A77334  B_ESCHR_COLI FALSE
#>  3 ICU      2002-01-07 067927  B_STPHY_EPDR TRUE 
#>  4 ICU      2002-01-07 067927  B_STPHY_EPDR FALSE
#>  5 ICU      2002-01-13 067927  B_STPHY_EPDR FALSE
#>  6 ICU      2002-01-13 067927  B_STPHY_EPDR FALSE
#>  7 Clinical 2002-01-14 462729  B_STPHY_AURS TRUE 
#>  8 Clinical 2002-01-14 462729  B_STPHY_AURS FALSE
#>  9 ICU      2002-01-16 067927  B_STPHY_EPDR TRUE 
#> 10 ICU      2002-01-17 858515  B_STPHY_EPDR TRUE 
#> # ℹ 1,990 more rows
# }
```
