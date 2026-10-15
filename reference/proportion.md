# Calculate Antimicrobial Resistance

These functions can be used to calculate the (co-)resistance or
susceptibility of microbial isolates (i.e. percentage of S, SI, I, IR or
R). All functions support quasiquotation with pipes, can be used in
[`summarise()`](https://dplyr.tidyverse.org/reference/summarise.html)
from the `dplyr` package and also support grouped variables, see
*Examples*.

`resistance()` should be used to calculate resistance,
`susceptibility()` should be used to calculate susceptibility.  

## Usage

``` r
resistance(
  ...,
  minimum = 30,
  as_percent = FALSE,
  only_all_tested = FALSE,
  guideline = getOption("AMR_guideline", "EUCAST")
)

susceptibility(
  ...,
  minimum = 30,
  as_percent = FALSE,
  only_all_tested = FALSE,
  guideline = getOption("AMR_guideline", "EUCAST")
)

sir_confidence_interval(
  ...,
  ab_result = "R",
  minimum = 30,
  as_percent = FALSE,
  only_all_tested = FALSE,
  confidence_level = 0.95,
  side = "both",
  collapse = FALSE
)

proportion_R(..., minimum = 30, as_percent = FALSE, only_all_tested = FALSE)

proportion_IR(..., minimum = 30, as_percent = FALSE, only_all_tested = FALSE)

proportion_I(..., minimum = 30, as_percent = FALSE, only_all_tested = FALSE)

proportion_SI(..., minimum = 30, as_percent = FALSE, only_all_tested = FALSE)

proportion_S(..., minimum = 30, as_percent = FALSE, only_all_tested = FALSE)

proportion_df(
  data,
  translate_ab = "name",
  language = get_AMR_locale(),
  minimum = 30,
  as_percent = FALSE,
  combine_SI,
  confidence_level = 0.95,
  guideline = getOption("AMR_guideline", "EUCAST")
)

sir_df(
  data,
  translate_ab = "name",
  language = get_AMR_locale(),
  minimum = 30,
  as_percent = FALSE,
  combine_SI,
  confidence_level = 0.95,
  guideline = getOption("AMR_guideline", "EUCAST")
)
```

## Arguments

- ...:

  One or more vectors (or columns) with antibiotic interpretations. They
  will be transformed internally with
  [`as.sir()`](https://amr-for-r.org/reference/as.sir.md) if needed. Use
  multiple columns to calculate (the lack of) co-resistance: the
  probability where one of two drugs have a resistant or susceptible
  result. See *Examples*.

- minimum:

  The minimum allowed number of available (tested) isolates. Any isolate
  count lower than `minimum` will return `NA` with a warning. The
  default number of `30` isolates is advised by the Clinical and
  Laboratory Standards Institute (CLSI) as best practice, see *Source*.

- as_percent:

  A [logical](https://rdrr.io/r/base/logical.html) to indicate whether
  the output must be returned as a hundred fold with % sign (a
  character). A value of `0.123456` will then be returned as `"12.3%"`.

- only_all_tested:

  (for combination therapies, i.e. using more than one variable for
  `...`): a [logical](https://rdrr.io/r/base/logical.html) to indicate
  that isolates must be tested for all antimicrobials, see section
  *Combination Therapy* below. Use
  [`sir_combine()`](https://amr-for-r.org/reference/sir_combine.md) to
  create a new variable for a combination therapy.

- guideline:

  A guideline name, `"EUCAST"` or `"CLSI"` (or any version thereof, e.g.
  `"CLSI 2025"`), to determine whether the 'I' category is considered
  susceptible (EUCAST) or resistant (CLSI). The default is the package
  option
  [`AMR_guideline`](https://amr-for-r.org/reference/AMR-options.md), or
  `"EUCAST"` if it is not set. For functions that return results per
  category, i.e.,
  [`count_df()`](https://amr-for-r.org/reference/count.md),
  `proportion_df()`, `sir_df()`,
  [`ggplot_sir()`](https://amr-for-r.org/reference/ggplot_sir.md),
  [`geom_sir()`](https://amr-for-r.org/reference/ggplot_sir.md), and
  [`labels_sir_count()`](https://amr-for-r.org/reference/ggplot_sir.md),
  this can also be `NULL` or `FALSE` to keep S, I, and R separate. See
  section *Interpretation of SIR* below.

- ab_result:

  Antibiotic results to test against, must be one or more values of "S",
  "SDD", "I", or "R".

- confidence_level:

  The confidence level for the returned confidence interval. For the
  calculation, the number of S or SI isolates, and R isolates are
  compared with the total number of available isolates with R, S, or I
  by using [`binom.test()`](https://rdrr.io/r/stats/binom.test.html),
  i.e., the Clopper-Pearson method.

- side:

  The side of the confidence interval to return. The default is `"both"`
  for a length 2 vector, but can also be (abbreviated as)
  `"min"`/`"left"`/`"lower"`/`"less"` or
  `"max"`/`"right"`/`"higher"`/`"greater"`.

- collapse:

  A [logical](https://rdrr.io/r/base/logical.html) to indicate whether
  the output values should be 'collapsed', i.e. be merged together into
  one value, or a character value to use for collapsing.

- data:

  A [data.frame](https://rdrr.io/r/base/data.frame.html) containing
  columns with class [`sir`](https://amr-for-r.org/reference/as.sir.md)
  (see [`as.sir()`](https://amr-for-r.org/reference/as.sir.md)).

- translate_ab:

  A column name of the
  [antimicrobials](https://amr-for-r.org/reference/antimicrobials.md)
  data set to translate the antibiotic abbreviations to, using
  [`ab_property()`](https://amr-for-r.org/reference/ab_property.md).

- language:

  Language of the returned text - the default is the current system
  language (see
  [`get_AMR_locale()`](https://amr-for-r.org/reference/translate.md))
  and can also be set with the package option
  [`AMR_locale`](https://amr-for-r.org/reference/AMR-options.md). Use
  `language = NULL` or `language = ""` to prevent translation.

- combine_SI:

  *(deprecated in favour of `guideline`)* A
  [logical](https://rdrr.io/r/base/logical.html). If set, it takes
  precedence over `guideline`: `TRUE` merges S, SDD, and I into `"SI"`
  vs. `"R"`, and `FALSE` keeps S, I, and R separate.

## Value

A [double](https://rdrr.io/r/base/double.html) or, when
`as_percent = TRUE`, a
[character](https://rdrr.io/r/base/character.html).

## Details

For a more automated and comprehensive analysis, consider using
[`antibiogram()`](https://amr-for-r.org/reference/antibiogram.md) or
[`wisca()`](https://amr-for-r.org/reference/antibiogram.md), which
streamline many aspects of susceptibility reporting and, importantly,
also support WISCA. The functions described here offer a more hands-on,
manual approach for greater customisation.

**Remember that you should filter your data to let it contain only first
isolates!** This is needed to exclude duplicates and to reduce selection
bias. Use
[`first_isolate()`](https://amr-for-r.org/reference/first_isolate.md) to
determine them in your data set with one of the four available
algorithms.

The function `resistance()` is equal to the function `proportion_R()`.
The function `susceptibility()` is equal to the function
`proportion_SI()`. Since AMR v3.0, `proportion_SI()` and
`proportion_I()` include dose-dependent susceptibility ('SDD').

Use `sir_confidence_interval()` to calculate the confidence interval,
which relies on
[`binom.test()`](https://rdrr.io/r/stats/binom.test.html), i.e., the
Clopper-Pearson method. This function returns a vector of length 2 at
default for antimicrobial *resistance*. Change the `side` argument to
"left"/"min" or "right"/"max" to return a single value, and change the
`ab_result` argument to e.g. `c("S", "I")` to test for antimicrobial
*susceptibility*, see Examples.

These functions are not meant to count isolates, but to calculate the
proportion of resistance/susceptibility. Use the
[`count_*()`](https://amr-for-r.org/reference/count.md) functions to
count isolates. The function `susceptibility()` is essentially equal to
[`count_susceptible()`](https://amr-for-r.org/reference/count.md)`/`[`count_all()`](https://amr-for-r.org/reference/count.md).
*Low counts can influence the outcome - the `proportion_*()` functions
may camouflage this, since they only return the proportion (albeit
dependent on the `minimum` argument).*

The function `proportion_df()` takes any variable from `data` that has
an [`sir`](https://amr-for-r.org/reference/as.sir.md) class (created
with [`as.sir()`](https://amr-for-r.org/reference/as.sir.md)) and
calculates the proportions S, I, and R. It also supports grouped
variables. The function `sir_df()` works exactly like `proportion_df()`,
but adds the number of isolates.

## Interpretation of SIR

In 2019, the European Committee on Antimicrobial Susceptibility Testing
(EUCAST) has decided to change the definitions of susceptibility testing
categories S, I, and R
(<https://www.eucast.org/bacteria/clinical-breakpoints-and-interpretation/definition-of-s-i-and-r/>).

In aggregated analyses, such as proportions, counts, antibiograms, and
plots, the 'I' category must be considered either susceptible or
resistant. In epidemiological AMR data analysis, this follows the
guideline: EUCAST considers 'I' susceptible, while CLSI considers 'I'
resistant. The dose-dependent susceptibility of CLSI ('SDD') is always
considered susceptible, as it equals what EUCAST defines as 'I'. All
functions with a `guideline` argument for this purpose use these
categories:

|             |               |               |
|-------------|---------------|---------------|
| `guideline` | Susceptible   | Resistant     |
| `"EUCAST"`  | S, SDD, I, WT | R, NWT, NS    |
| `"CLSI"`    | S, SDD, WT    | I, R, NWT, NS |

The values 'WT' (wild type) and 'NWT' (non-wild type) result from
epidemiological cut-off values (ECOFFs) instead of clinical breakpoints.
Since they are no clinical categories, 'WT' is always considered
susceptible and 'NWT' always resistant, regardless of the guideline. The
value 'NS' (non-susceptible) is always considered resistant, and 'NI'
(non-interpretable) is always considered a missing result.

The `guideline` argument can be:

- **Left unset**, to use the package option
  [`AMR_guideline`](https://amr-for-r.org/reference/AMR-options.md), or
  `"EUCAST"` if this option is not set. In the latter case, a message is
  shown once per session (but not in R Markdown or Quarto documents);
  set the option to silence it.

- **`"EUCAST"`**, or any EUCAST version such as `"EUCAST 2025"`, to
  consider 'I' susceptible, e.g., to return `"SI"` vs. `"R"` in
  `proportion_df()`.

- **`"CLSI"`**, or any CLSI version such as `"CLSI 2025"`, to consider
  'I' resistant, e.g., to return `"S"` vs. `"IR"` in `proportion_df()`.

- **`NULL` or `FALSE`**, to keep S, I, and R separate. This is only
  possible in functions that return results per category:
  [`count_df()`](https://amr-for-r.org/reference/count.md),
  `proportion_df()`, `sir_df()`,
  [`ggplot_sir()`](https://amr-for-r.org/reference/ggplot_sir.md),
  [`geom_sir()`](https://amr-for-r.org/reference/ggplot_sir.md),
  [`labels_sir_count()`](https://amr-for-r.org/reference/ggplot_sir.md),
  and
  [`bug_drug_combinations()`](https://amr-for-r.org/reference/bug_drug_combinations.md).
  Functions that return a single proportion or count, such as
  `susceptibility()`,
  [`count_resistant()`](https://amr-for-r.org/reference/count.md),
  [`antibiogram()`](https://amr-for-r.org/reference/antibiogram.md),
  [`wisca()`](https://amr-for-r.org/reference/antibiogram.md), and
  [`format()`](https://rdrr.io/r/base/format.html) on the result of
  [`bug_drug_combinations()`](https://amr-for-r.org/reference/bug_drug_combinations.md),
  require 'I' to be considered either susceptible or resistant, and
  return an error.

Use `susceptibility()` and `resistance()` (or
[`count_susceptible()`](https://amr-for-r.org/reference/count.md) and
[`count_resistant()`](https://amr-for-r.org/reference/count.md)) to
follow the guideline, or e.g. `proportion_SI()` and `proportion_R()` to
be explicit regardless of the guideline.

The `combine_SI` argument of the functions that return results per
category,
[`antibiogram()`](https://amr-for-r.org/reference/antibiogram.md), and
[`wisca()`](https://amr-for-r.org/reference/antibiogram.md) is
deprecated in favour of `guideline`. If set, it is followed literally
and takes precedence over `guideline`.

## Combination Therapy

Results are ordered from susceptible to resistant: S (and WT), SDD, I, R
(and NWT, NS). With `method = "most_susceptible"` (in
[`sir_combine()`](https://amr-for-r.org/reference/sir_combine.md)), the
most susceptible result of the antimicrobials is returned, and with
`method = "most_resistant"`, the most resistant result. The value `"NI"`
(non-interpretable) is considered a missing result.

A missing result can only be ignored if it cannot change the conclusion:
a susceptible result for `"most_susceptible"`, and a non-susceptible
result for `"most_resistant"`. Under EUCAST, the 'I' category is
considered susceptible (*susceptible, increased exposure*), and under
CLSI, it is not. With `only_all_tested = TRUE`, any missing result leads
to `NA`.

This is an example for two antimicrobials, Drug A and Drug B, under
EUCAST (the order of the antimicrobials does not matter):

    -------------------------------------------------------------------------
                          Result of method =            Result of method =
                         "most_susceptible"              "most_resistant"
                       -----------------------  -----------------------------
                        only_all_tested =        only_all_tested =
     Drug A    Drug B     FALSE      TRUE          FALSE         TRUE
    --------  --------  ---------  ---------    -----------  -------------
        S        S          S          S             S             S
        S        I          S          S             I             I
        S        R          S          S             R             R
        I        R          I          I             R             R
        R        R          R          R             R             R
        S      <NA>         S        <NA>          <NA>          <NA>
        I      <NA>         I        <NA>          <NA>          <NA>
        R      <NA>       <NA>       <NA>            R           <NA>
      <NA>     <NA>       <NA>       <NA>          <NA>          <NA>
    -------------------------------------------------------------------------

Under CLSI, the combination of I with a missing result is `<NA>` for
`"most_susceptible"` and I for `"most_resistant"`.

The [`count_*()`](https://amr-for-r.org/reference/count.md) and
`proportion_*()` functions apply these rules when more than one variable
is used. With `only_all_tested = FALSE` (default), the functions for
susceptibility, such as `susceptibility()` and
[`count_susceptible()`](https://amr-for-r.org/reference/count.md), count
an isolate as susceptible if any of the antimicrobials is susceptible
(`"most_susceptible"`), and the functions for resistance, such as
`resistance()` and
[`count_resistant()`](https://amr-for-r.org/reference/count.md), count
an isolate as resistant if any of the antimicrobials is resistant
(`"most_resistant"`). An isolate is considered tested if it has such a
result, or if it was tested for all antimicrobials. Therefore, for
combinations with `only_all_tested = FALSE`:

        count_S()    +   count_I()    +   count_R()    >= count_all()
      proportion_S() + proportion_I() + proportion_R() >= 1

With `only_all_tested = TRUE`, all these functions count the result of
the combination therapy (`"most_susceptible"`) of isolates tested for
all antimicrobials, so that:

        count_S()    +   count_I()    +   count_R()    = count_all()
      proportion_S() + proportion_I() + proportion_R() = 1

In [`antibiogram()`](https://amr-for-r.org/reference/antibiogram.md) and
[`wisca()`](https://amr-for-r.org/reference/antibiogram.md), the result
of the combination therapy (`"most_susceptible"`) is used, the 'I'
category is considered susceptible or not following the `guideline`
argument ('SDD' is always considered susceptible), and the result is
then reduced to S or R.

## References

**M39 Analysis and Presentation of Cumulative Antimicrobial
Susceptibility Test Data, 5th Edition**, 2022, *Clinical and Laboratory
Standards Institute (CLSI)*.
<https://clsi.org/standards/products/microbiology/documents/m39/>.

## See also

[`count()`](https://amr-for-r.org/reference/count.md) to count resistant
and susceptible isolates.

## Examples

``` r
# example_isolates is a data set available in the AMR package.
# run ?example_isolates for more info.
example_isolates
#> # A tibble: 2,000 × 46
#>    date       patient   age gender ward     mo           PEN   OXA   FLC   AMX  
#>    <date>     <chr>   <dbl> <chr>  <chr>    <mo>         <sir> <sir> <sir> <sir>
#>  1 2002-01-02 A77334     65 F      Clinical B_ESCHR_COLI   R     NA    NA    NA 
#>  2 2002-01-03 A77334     65 F      Clinical B_ESCHR_COLI   R     NA    NA    NA 
#>  3 2002-01-07 067927     45 F      ICU      B_STPHY_EPDR   R     NA    R     NA 
#>  4 2002-01-07 067927     45 F      ICU      B_STPHY_EPDR   R     NA    R     NA 
#>  5 2002-01-13 067927     45 F      ICU      B_STPHY_EPDR   R     NA    R     NA 
#>  6 2002-01-13 067927     45 F      ICU      B_STPHY_EPDR   R     NA    R     NA 
#>  7 2002-01-14 462729     78 M      Clinical B_STPHY_AURS   R     NA    S     R  
#>  8 2002-01-14 462729     78 M      Clinical B_STPHY_AURS   R     NA    S     R  
#>  9 2002-01-16 067927     45 F      ICU      B_STPHY_EPDR   R     NA    R     NA 
#> 10 2002-01-17 858515     79 F      ICU      B_STPHY_EPDR   R     NA    S     NA 
#> # ℹ 1,990 more rows
#> # ℹ 36 more variables: AMC <sir>, AMP <sir>, TZP <sir>, CZO <sir>, FEP <sir>,
#> #   CXM <sir>, FOX <sir>, CTX <sir>, CAZ <sir>, CRO <sir>, GEN <sir>,
#> #   TOB <sir>, AMK <sir>, KAN <sir>, TMP <sir>, SXT <sir>, NIT <sir>,
#> #   FOS <sir>, LNZ <sir>, CIP <sir>, MFX <sir>, VAN <sir>, TEC <sir>,
#> #   TCY <sir>, TGC <sir>, DOX <sir>, ERY <sir>, CLI <sir>, AZM <sir>,
#> #   IPM <sir>, MEM <sir>, MTR <sir>, CHL <sir>, COL <sir>, MUP <sir>, …


# base R ------------------------------------------------------------
# determines %R
resistance(example_isolates$AMX)
#> [1] 0.5955556
sir_confidence_interval(example_isolates$AMX)
#> [1] 0.5688204 0.6218738
sir_confidence_interval(example_isolates$AMX,
  confidence_level = 0.975
)
#> [1] 0.5650148 0.6255670
sir_confidence_interval(example_isolates$AMX,
  confidence_level = 0.975,
  collapse = ", "
)
#> [1] "0.565, 0.626"

# determines %S+I:
susceptibility(example_isolates$AMX)
#> [1] 0.4044444
sir_confidence_interval(example_isolates$AMX,
  ab_result = c("S", "I")
)
#> [1] 0.3781262 0.4311796

# be more specific
proportion_S(example_isolates$AMX)
#> [1] 0.4022222
proportion_SI(example_isolates$AMX)
#> [1] 0.4044444
proportion_I(example_isolates$AMX)
#> [1] 0.002222222
proportion_IR(example_isolates$AMX)
#> [1] 0.5977778
proportion_R(example_isolates$AMX)
#> [1] 0.5955556

# dplyr -------------------------------------------------------------
# \donttest{
if (require("dplyr")) {
  example_isolates %>%
    group_by(ward) %>%
    summarise(
      r = resistance(CIP),
      n = n_sir(CIP)
    ) # n_sir works like n_distinct in dplyr, see ?n_sir
}
#> # A tibble: 3 × 3
#>   ward           r     n
#>   <chr>      <dbl> <int>
#> 1 Clinical   0.147   869
#> 2 ICU        0.190   447
#> 3 Outpatient 0.161    93
if (require("dplyr")) {
  example_isolates %>%
    group_by(ward) %>%
    summarise(
      cipro_R = resistance(CIP),
      ci_min = sir_confidence_interval(CIP, side = "min"),
      ci_max = sir_confidence_interval(CIP, side = "max"),
    )
}
#> # A tibble: 3 × 4
#>   ward       cipro_R ci_min ci_max
#>   <chr>        <dbl>  <dbl>  <dbl>
#> 1 Clinical     0.147 0.124   0.173
#> 2 ICU          0.190 0.155   0.230
#> 3 Outpatient   0.161 0.0932  0.252
if (require("dplyr")) {
  # scoped dplyr verbs with antimicrobial selectors
  # (you could also use across() of course)
  example_isolates %>%
    group_by(ward) %>%
    summarise_at(
      c(aminoglycosides(), carbapenems()),
      resistance
    )
}
#> ℹ For `aminoglycosides()` using columns GEN (gentamicin), TOB (tobramycin), AMK
#>   (amikacin), and KAN (kanamycin)
#> ℹ For `carbapenems()` using columns IPM (imipenem) and MEM (meropenem)
#> Warning: There was 1 warning in `summarise()`.
#> ℹ In argument: `KAN = (function (..., minimum = 30, as_percent = FALSE,
#>   only_all_tested = FALSE, ...`.
#> ℹ In group 3: `ward = "Outpatient"`.
#> Caused by warning:
#> ! Introducing NA: only 23 results available for KAN in group: ward = "Outpatient"
#> (whilst `minimum = 30`).
#> # A tibble: 3 × 7
#>   ward         GEN   TOB   AMK   KAN    IPM    MEM
#>   <chr>      <dbl> <dbl> <dbl> <dbl>  <dbl>  <dbl>
#> 1 Clinical   0.229 0.315 0.626     1 0.0498 0.0458
#> 2 ICU        0.290 0.400 0.662     1 0.0862 0.0894
#> 3 Outpatient 0.2   0.368 0.605    NA 0.0541 0.0541
if (require("dplyr")) {
  example_isolates %>%
    group_by(ward) %>%
    summarise(
      R = resistance(CIP, as_percent = TRUE),
      SI = susceptibility(CIP, as_percent = TRUE),
      n1 = count_all(CIP), # the actual total; sum of all three
      n2 = n_sir(CIP), # same - analogous to n_distinct
      total = n()
    ) # NOT the number of tested isolates!

  # Calculate co-resistance between amoxicillin/clav acid and gentamicin,
  # so we can see that combination therapy does a lot more than mono therapy:
  example_isolates %>% susceptibility(AMC) # %SI = 76.3%
  example_isolates %>% count_all(AMC) #   n = 1879

  example_isolates %>% susceptibility(GEN) # %SI = 75.4%
  example_isolates %>% count_all(GEN) #   n = 1855

  example_isolates %>% susceptibility(AMC, GEN) # %SI = 94.1%
  example_isolates %>% count_all(AMC, GEN) #   n = 1939


  # See Details on how `only_all_tested` works. Example:
  example_isolates %>%
    summarise(
      numerator = count_susceptible(AMC, GEN),
      denominator = count_all(AMC, GEN),
      proportion = susceptibility(AMC, GEN)
    )

  example_isolates %>%
    summarise(
      numerator = count_susceptible(AMC, GEN, only_all_tested = TRUE),
      denominator = count_all(AMC, GEN, only_all_tested = TRUE),
      proportion = susceptibility(AMC, GEN, only_all_tested = TRUE)
    )


  example_isolates %>%
    group_by(ward) %>%
    summarise(
      cipro_p = susceptibility(CIP, as_percent = TRUE),
      cipro_n = count_all(CIP),
      genta_p = susceptibility(GEN, as_percent = TRUE),
      genta_n = count_all(GEN),
      combination_p = susceptibility(CIP, GEN, as_percent = TRUE),
      combination_n = count_all(CIP, GEN)
    )

  # Get proportions S/I/R immediately of all sir columns
  example_isolates %>%
    select(AMX, CIP) %>%
    proportion_df(translate = FALSE)

  # It also supports grouping variables
  # (use sir_df to also include the count)
  example_isolates %>%
    select(ward, AMX, CIP) %>%
    group_by(ward) %>%
    sir_df(translate = FALSE)
}
#> ℹ `proportion_df()` assumes the EUCAST guideline and thus considers the 'I'
#>   category susceptible (S+SDD+I vs. R). Use `guideline = "CLSI"` for S+SDD vs.
#>   I+R, or `guideline = NULL` to keep S, I, and R separate. Set the `guideline`
#>   argument, or set the `AMR_guideline` option to change the default (see
#>   `?AMR-options`). This message is shown once per session.
#> ℹ `sir_df()` assumes the EUCAST guideline and thus considers the 'I' category
#>   susceptible (S+SDD+I vs. R). Use `guideline = "CLSI"` for S+SDD vs. I+R, or
#>   `guideline = NULL` to keep S, I, and R separate. Set the `guideline`
#>   argument, or set the `AMR_guideline` option to change the default (see
#>   `?AMR-options`). This message is shown once per session.
#> # A tibble: 12 × 7
#>    ward       antibiotic interpretation value ci_min ci_max isolates
#>    <chr>      <chr>      <ord>          <dbl>  <dbl>  <dbl>    <int>
#>  1 Clinical   AMX        SI             0.423 0.389   0.457      357
#>  2 Clinical   AMX        R              0.577 0.543   0.611      487
#>  3 Clinical   CIP        SI             0.853 0.827   0.876      741
#>  4 Clinical   CIP        R              0.147 0.124   0.173      128
#>  5 ICU        AMX        SI             0.369 0.323   0.417      158
#>  6 ICU        AMX        R              0.631 0.583   0.677      270
#>  7 ICU        CIP        SI             0.810 0.770   0.845      362
#>  8 ICU        CIP        R              0.190 0.155   0.230       85
#>  9 Outpatient AMX        SI             0.397 0.288   0.515       31
#> 10 Outpatient AMX        R              0.603 0.485   0.712       47
#> 11 Outpatient CIP        SI             0.839 0.748   0.907       78
#> 12 Outpatient CIP        R              0.161 0.0932  0.252       15
# }
```
