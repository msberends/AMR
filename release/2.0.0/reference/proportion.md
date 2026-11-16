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
resistance(..., minimum = 30, as_percent = FALSE, only_all_tested = FALSE)

susceptibility(..., minimum = 30, as_percent = FALSE, only_all_tested = FALSE)

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
  combine_SI = TRUE,
  confidence_level = 0.95
)

sir_df(
  data,
  translate_ab = "name",
  language = get_AMR_locale(),
  minimum = 30,
  as_percent = FALSE,
  combine_SI = TRUE,
  confidence_level = 0.95
)
```

## Source

**M39 Analysis and Presentation of Cumulative Antimicrobial
Susceptibility Test Data, 5th Edition**, 2022, *Clinical and Laboratory
Standards Institute (CLSI)*.
<https://clsi.org/standards/products/microbiology/documents/m39/>.

## Arguments

- ...:

  one or more vectors (or columns) with antibiotic interpretations. They
  will be transformed internally with
  [`as.sir()`](https://amr-for-r.org/release/2.0.0/reference/as.sir.md)
  if needed. Use multiple columns to calculate (the lack of)
  co-resistance: the probability where one of two drugs have a resistant
  or susceptible result. See *Examples*.

- minimum:

  the minimum allowed number of available (tested) isolates. Any isolate
  count lower than `minimum` will return `NA` with a warning. The
  default number of `30` isolates is advised by the Clinical and
  Laboratory Standards Institute (CLSI) as best practice, see *Source*.

- as_percent:

  a [logical](https://rdrr.io/r/base/logical.html) to indicate whether
  the output must be returned as a hundred fold with % sign (a
  character). A value of `0.123456` will then be returned as `"12.3%"`.

- only_all_tested:

  (for combination therapies, i.e. using more than one variable for
  `...`): a [logical](https://rdrr.io/r/base/logical.html) to indicate
  that isolates must be tested for all antibiotics, see section
  *Combination Therapy* below

- ab_result:

  antibiotic results to test against, must be one or more values of "S",
  "I", or "R"

- confidence_level:

  the confidence level for the returned confidence interval. For the
  calculation, the number of S or SI isolates, and R isolates are
  compared with the total number of available isolates with R, S, or I
  by using [`binom.test()`](https://rdrr.io/r/stats/binom.test.html),
  i.e., the Clopper-Pearson method.

- side:

  the side of the confidence interval to return. The default is `"both"`
  for a length 2 vector, but can also be (abbreviated as)
  `"min"`/`"left"`/`"lower"`/`"less"` or
  `"max"`/`"right"`/`"higher"`/`"greater"`.

- collapse:

  a [logical](https://rdrr.io/r/base/logical.html) to indicate whether
  the output values should be 'collapsed', i.e. be merged together into
  one value, or a character value to use for collapsing

- data:

  a [data.frame](https://rdrr.io/r/base/data.frame.html) containing
  columns with class
  [`sir`](https://amr-for-r.org/release/2.0.0/reference/as.sir.md) (see
  [`as.sir()`](https://amr-for-r.org/release/2.0.0/reference/as.sir.md))

- translate_ab:

  a column name of the
  [antibiotics](https://amr-for-r.org/release/2.0.0/reference/antibiotics.md)
  data set to translate the antibiotic abbreviations to, using
  [`ab_property()`](https://amr-for-r.org/release/2.0.0/reference/ab_property.md)

- language:

  language of the returned text - the default is the current system
  language (see
  [`get_AMR_locale()`](https://amr-for-r.org/release/2.0.0/reference/translate.md))
  and can also be set with the [package
  option](https://amr-for-r.org/release/2.0.0/reference/AMR-options.md)
  [`AMR_locale`](https://amr-for-r.org/release/2.0.0/reference/AMR-options.md).
  Use `language = NULL` or `language = ""` to prevent translation.

- combine_SI:

  a [logical](https://rdrr.io/r/base/logical.html) to indicate whether
  all values of S and I must be merged into one, so the output only
  consists of S+I vs. R (susceptible vs. resistant) - the default is
  `TRUE`

## Value

A [double](https://rdrr.io/r/base/double.html) or, when
`as_percent = TRUE`, a
[character](https://rdrr.io/r/base/character.html).

## Details

The function `resistance()` is equal to the function `proportion_R()`.
The function `susceptibility()` is equal to the function
`proportion_SI()`.

Use `sir_confidence_interval()` to calculate the confidence interval,
which relies on
[`binom.test()`](https://rdrr.io/r/stats/binom.test.html), i.e., the
Clopper-Pearson method. This function returns a vector of length 2 at
default for antimicrobial *resistance*. Change the `side` argument to
"left"/"min" or "right"/"max" to return a single value, and change the
`ab_result` argument to e.g. `c("S", "I")` to test for antimicrobial
*susceptibility*, see Examples.

**Remember that you should filter your data to let it contain only first
isolates!** This is needed to exclude duplicates and to reduce selection
bias. Use
[`first_isolate()`](https://amr-for-r.org/release/2.0.0/reference/first_isolate.md)
to determine them in your data set with one of the four available
algorithms.

These functions are not meant to count isolates, but to calculate the
proportion of resistance/susceptibility. Use the
[`count()`](https://amr-for-r.org/release/2.0.0/reference/count.md)
functions to count isolates. The function `susceptibility()` is
essentially equal to `count_susceptible() / count_all()`. *Low counts
can influence the outcome - the `proportion` functions may camouflage
this, since they only return the proportion (albeit being dependent on
the `minimum` argument).*

The function `proportion_df()` takes any variable from `data` that has
an [`sir`](https://amr-for-r.org/release/2.0.0/reference/as.sir.md)
class (created with
[`as.sir()`](https://amr-for-r.org/release/2.0.0/reference/as.sir.md))
and calculates the proportions S, I, and R. It also supports grouped
variables. The function `sir_df()` works exactly like `proportion_df()`,
but adds the number of isolates.

## Combination Therapy

When using more than one variable for `...` (= combination therapy), use
`only_all_tested` to only count isolates that are tested for all
antibiotics/variables that you test them for. See this example for two
antibiotics, Drug A and Drug B, about how `susceptibility()` works to
calculate the %SI:

    --------------------------------------------------------------------
                        only_all_tested = FALSE  only_all_tested = TRUE
                        -----------------------  -----------------------
     Drug A    Drug B   include as  include as   include as  include as
                        numerator   denominator  numerator   denominator
    --------  --------  ----------  -----------  ----------  -----------
     S or I    S or I       X            X            X            X
       R       S or I       X            X            X            X
      <NA>     S or I       X            X            -            -
     S or I      R          X            X            X            X
       R         R          -            X            -            X
      <NA>       R          -            -            -            -
     S or I     <NA>        X            X            -            -
       R        <NA>        -            -            -            -
      <NA>      <NA>        -            -            -            -
    --------------------------------------------------------------------

Please note that, in combination therapies, for `only_all_tested = TRUE`
applies that:

        count_S()    +   count_I()    +   count_R()    = count_all()
      proportion_S() + proportion_I() + proportion_R() = 1

and that, in combination therapies, for `only_all_tested = FALSE`
applies that:

        count_S()    +   count_I()    +   count_R()    >= count_all()
      proportion_S() + proportion_I() + proportion_R() >= 1

Using `only_all_tested` has no impact when only using one antibiotic as
input.

## Interpretation of SIR

In 2019, the European Committee on Antimicrobial Susceptibility Testing
(EUCAST) has decided to change the definitions of susceptibility testing
categories S, I, and R as shown below
(<https://www.eucast.org/newsiandr/>):

- **S - Susceptible, standard dosing regimen**  
  A microorganism is categorised as "Susceptible, standard dosing
  regimen", when there is a high likelihood of therapeutic success using
  a standard dosing regimen of the agent.

- **I - Susceptible, increased exposure** *  
  A microorganism is categorised as "Susceptible, Increased exposure*"
  when there is a high likelihood of therapeutic success because
  exposure to the agent is increased by adjusting the dosing regimen or
  by its concentration at the site of infection.

- **R = Resistant**  
  A microorganism is categorised as "Resistant" when there is a high
  likelihood of therapeutic failure even when there is increased
  exposure.

  - *Exposure* is a function of how the mode of administration, dose,
    dosing interval, infusion time, as well as distribution and
    excretion of the antimicrobial agent will influence the infecting
    organism at the site of infection.

This AMR package honours this insight. Use `susceptibility()` (equal to
`proportion_SI()`) to determine antimicrobial susceptibility and
[`count_susceptible()`](https://amr-for-r.org/release/2.0.0/reference/count.md)
(equal to
[`count_SI()`](https://amr-for-r.org/release/2.0.0/reference/count.md))
to count susceptible isolates.

## See also

[`count()`](https://amr-for-r.org/release/2.0.0/reference/count.md) to
count resistant and susceptible isolates.

## Examples

``` r
# example_isolates is a data set available in the AMR package.
# run ?example_isolates for more info.
example_isolates


# base R ------------------------------------------------------------
# determines %R
resistance(example_isolates$AMX)
sir_confidence_interval(example_isolates$AMX)
sir_confidence_interval(example_isolates$AMX,
  confidence_level = 0.975
)
sir_confidence_interval(example_isolates$AMX,
  confidence_level = 0.975,
  collapse = ", "
)

# determines %S+I:
susceptibility(example_isolates$AMX)
sir_confidence_interval(example_isolates$AMX,
  ab_result = c("S", "I")
)

# be more specific
proportion_S(example_isolates$AMX)
proportion_SI(example_isolates$AMX)
proportion_I(example_isolates$AMX)
proportion_IR(example_isolates$AMX)
proportion_R(example_isolates$AMX)

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
if (require("dplyr")) {
  example_isolates %>%
    group_by(ward) %>%
    summarise(
      cipro_R = resistance(CIP),
      ci_min = sir_confidence_interval(CIP, side = "min"),
      ci_max = sir_confidence_interval(CIP, side = "max"),
    )
}
if (require("dplyr")) {
  # scoped dplyr verbs with antibiotic selectors
  # (you could also use across() of course)
  example_isolates %>%
    group_by(ward) %>%
    summarise_at(
      c(aminoglycosides(), carbapenems()),
      resistance
    )
}
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
# }
```
