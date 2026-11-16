# Count Available Isolates

These functions can be used to count resistant/susceptible microbial
isolates. All functions support quasiquotation with pipes, can be used
in [`summarise()`](https://dplyr.tidyverse.org/reference/summarise.html)
from the `dplyr` package and also support grouped variables, see
*Examples*.

`count_resistant()` should be used to count resistant isolates,
`count_susceptible()` should be used to count susceptible isolates.

## Usage

``` r
count_resistant(..., only_all_tested = FALSE)

count_susceptible(..., only_all_tested = FALSE)

count_R(..., only_all_tested = FALSE)

count_IR(..., only_all_tested = FALSE)

count_I(..., only_all_tested = FALSE)

count_SI(..., only_all_tested = FALSE)

count_S(..., only_all_tested = FALSE)

count_all(..., only_all_tested = FALSE)

n_sir(..., only_all_tested = FALSE)

count_df(
  data,
  translate_ab = "name",
  language = get_AMR_locale(),
  combine_SI = TRUE
)
```

## Arguments

- ...:

  one or more vectors (or columns) with antibiotic interpretations. They
  will be transformed internally with
  [`as.sir()`](https://amr-for-r.org/release/2.0.0/reference/as.sir.md)
  if needed.

- only_all_tested:

  (for combination therapies, i.e. using more than one variable for
  `...`): a [logical](https://rdrr.io/r/base/logical.html) to indicate
  that isolates must be tested for all antibiotics, see section
  *Combination Therapy* below

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

An [integer](https://rdrr.io/r/base/integer.html)

## Details

These functions are meant to count isolates. Use the
[`resistance()`](https://amr-for-r.org/release/2.0.0/reference/proportion.md)/[`susceptibility()`](https://amr-for-r.org/release/2.0.0/reference/proportion.md)
functions to calculate microbial resistance/susceptibility.

The function `count_resistant()` is equal to the function `count_R()`.
The function `count_susceptible()` is equal to the function
`count_SI()`.

The function `n_sir()` is an alias of `count_all()`. They can be used to
count all available isolates, i.e. where all input antibiotics have an
available result (S, I or R). Their use is equal to
[`n_distinct()`](https://dplyr.tidyverse.org/reference/n_distinct.html).
Their function is equal to
`count_susceptible(...) + count_resistant(...)`.

The function `count_df()` takes any variable from `data` that has an
[`sir`](https://amr-for-r.org/release/2.0.0/reference/as.sir.md) class
(created with
[`as.sir()`](https://amr-for-r.org/release/2.0.0/reference/as.sir.md))
and counts the number of S's, I's and R's. It also supports grouped
variables. The function
[`sir_df()`](https://amr-for-r.org/release/2.0.0/reference/proportion.md)
works exactly like `count_df()`, but adds the percentage of S, I and R.

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

This AMR package honours this insight. Use
[`susceptibility()`](https://amr-for-r.org/release/2.0.0/reference/proportion.md)
(equal to
[`proportion_SI()`](https://amr-for-r.org/release/2.0.0/reference/proportion.md))
to determine antimicrobial susceptibility and `count_susceptible()`
(equal to `count_SI()`) to count susceptible isolates.

## Combination Therapy

When using more than one variable for `...` (= combination therapy), use
`only_all_tested` to only count isolates that are tested for all
antibiotics/variables that you test them for. See this example for two
antibiotics, Drug A and Drug B, about how
[`susceptibility()`](https://amr-for-r.org/release/2.0.0/reference/proportion.md)
works to calculate the %SI:

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

## See also

[`proportion_*`](https://amr-for-r.org/release/2.0.0/reference/proportion.md)
to calculate microbial resistance and susceptibility.

## Examples

``` r
# example_isolates is a data set available in the AMR package.
# run ?example_isolates for more info.

# base R ------------------------------------------------------------
count_resistant(example_isolates$AMX) # counts "R"
count_susceptible(example_isolates$AMX) # counts "S" and "I"
count_all(example_isolates$AMX) # counts "S", "I" and "R"

# be more specific
count_S(example_isolates$AMX)
count_SI(example_isolates$AMX)
count_I(example_isolates$AMX)
count_IR(example_isolates$AMX)
count_R(example_isolates$AMX)

# Count all available isolates
count_all(example_isolates$AMX)
n_sir(example_isolates$AMX)

# n_sir() is an alias of count_all().
# Since it counts all available isolates, you can
# calculate back to count e.g. susceptible isolates.
# These results are the same:
count_susceptible(example_isolates$AMX)
susceptibility(example_isolates$AMX) * n_sir(example_isolates$AMX)

# dplyr -------------------------------------------------------------
# \donttest{
if (require("dplyr")) {
  example_isolates %>%
    group_by(ward) %>%
    summarise(
      R = count_R(CIP),
      I = count_I(CIP),
      S = count_S(CIP),
      n1 = count_all(CIP), # the actual total; sum of all three
      n2 = n_sir(CIP), # same - analogous to n_distinct
      total = n()
    ) # NOT the number of tested isolates!

  # Number of available isolates for a whole antibiotic class
  # (i.e., in this data set columns GEN, TOB, AMK, KAN)
  example_isolates %>%
    group_by(ward) %>%
    summarise(across(aminoglycosides(), n_sir))

  # Count co-resistance between amoxicillin/clav acid and gentamicin,
  # so we can see that combination therapy does a lot more than mono therapy.
  # Please mind that `susceptibility()` calculates percentages right away instead.
  example_isolates %>% count_susceptible(AMC) # 1433
  example_isolates %>% count_all(AMC) # 1879

  example_isolates %>% count_susceptible(GEN) # 1399
  example_isolates %>% count_all(GEN) # 1855

  example_isolates %>% count_susceptible(AMC, GEN) # 1764
  example_isolates %>% count_all(AMC, GEN) # 1936

  # Get number of S+I vs. R immediately of selected columns
  example_isolates %>%
    select(AMX, CIP) %>%
    count_df(translate = FALSE)

  # It also supports grouping variables
  example_isolates %>%
    select(ward, AMX, CIP) %>%
    group_by(ward) %>%
    count_df(translate = FALSE)
}
# }
```
