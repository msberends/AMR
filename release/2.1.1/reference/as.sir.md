# Translate MIC and Disk Diffusion to SIR, or Clean Existing SIR Data

Interpret minimum inhibitory concentration (MIC) values and disk
diffusion diameters according to EUCAST or CLSI, or clean up existing
SIR values. This transforms the input to a new class `sir`, which is an
ordered [factor](https://rdrr.io/pkg/data.table/man/fctr.html) with
levels `S < I < R`.

Currently available **breakpoint guidelines** are EUCAST 2011-2023 and
CLSI 2011-2023, and available **breakpoint types** are "ECOFF",
"animal", and "human".

All breakpoints used for interpretation are publicly available in the
[clinical_breakpoints](https://amr-for-r.org/release/2.1.1/reference/clinical_breakpoints.md)
data set.

## Usage

``` r
as.sir(x, ...)

NA_sir_

is.sir(x)

is_sir_eligible(x, threshold = 0.05)

# S3 method for class 'mic'
as.sir(
  x,
  mo = NULL,
  ab = deparse(substitute(x)),
  guideline = getOption("AMR_guideline", "EUCAST"),
  uti = NULL,
  conserve_capped_values = FALSE,
  add_intrinsic_resistance = FALSE,
  reference_data = AMR::clinical_breakpoints,
  include_screening = getOption("AMR_include_screening", FALSE),
  include_PKPD = getOption("AMR_include_PKPD", TRUE),
  breakpoint_type = getOption("AMR_breakpoint_type", "human"),
  ...
)

# S3 method for class 'disk'
as.sir(
  x,
  mo = NULL,
  ab = deparse(substitute(x)),
  guideline = getOption("AMR_guideline", "EUCAST"),
  uti = NULL,
  add_intrinsic_resistance = FALSE,
  reference_data = AMR::clinical_breakpoints,
  include_screening = getOption("AMR_include_screening", FALSE),
  include_PKPD = getOption("AMR_include_PKPD", TRUE),
  breakpoint_type = getOption("AMR_breakpoint_type", "human"),
  ...
)

# S3 method for class 'data.frame'
as.sir(
  x,
  ...,
  col_mo = NULL,
  guideline = getOption("AMR_guideline", "EUCAST"),
  uti = NULL,
  conserve_capped_values = FALSE,
  add_intrinsic_resistance = FALSE,
  reference_data = AMR::clinical_breakpoints,
  include_screening = getOption("AMR_include_screening", FALSE),
  include_PKPD = getOption("AMR_include_PKPD", TRUE),
  breakpoint_type = getOption("AMR_breakpoint_type", "human")
)

sir_interpretation_history(clean = FALSE)
```

## Format

An object of class `sir` (inherits from `ordered`, `factor`) of length
1.

## Source

For interpretations of minimum inhibitory concentration (MIC) values and
disk diffusion diameters:

- **M39 Analysis and Presentation of Cumulative Antimicrobial
  Susceptibility Test Data**, 2011-2023, *Clinical and Laboratory
  Standards Institute* (CLSI).
  <https://clsi.org/standards/products/microbiology/documents/m39/>.

- **M100 Performance Standard for Antimicrobial Susceptibility
  Testing**, 2011-2023, *Clinical and Laboratory Standards Institute*
  (CLSI).
  <https://clsi.org/standards/products/microbiology/documents/m100/>.

- **Breakpoint tables for interpretation of MICs and zone diameters**,
  2011-2023, *European Committee on Antimicrobial Susceptibility
  Testing* (EUCAST). <https://www.eucast.org/clinical_breakpoints>.

## Arguments

- x:

  vector of values (for class
  [`mic`](https://amr-for-r.org/release/2.1.1/reference/as.mic.md): MIC
  values in mg/L, for class
  [`disk`](https://amr-for-r.org/release/2.1.1/reference/as.disk.md): a
  disk diffusion radius in millimetres)

- ...:

  for using on a [data.frame](https://rdrr.io/r/base/data.frame.html):
  names of columns to apply `as.sir()` on (supports tidy selection such
  as `column1:column4`). Otherwise: arguments passed on to methods.

- threshold:

  maximum fraction of invalid antimicrobial interpretations of `x`, see
  *Examples*

- mo:

  any (vector of) text that can be coerced to valid microorganism codes
  with
  [`as.mo()`](https://amr-for-r.org/release/2.1.1/reference/as.mo.md),
  can be left empty to determine it automatically

- ab:

  any (vector of) text that can be coerced to a valid antimicrobial drug
  code with
  [`as.ab()`](https://amr-for-r.org/release/2.1.1/reference/as.ab.md)

- guideline:

  defaults to EUCAST 2023 (the latest implemented EUCAST guideline in
  the
  [clinical_breakpoints](https://amr-for-r.org/release/2.1.1/reference/clinical_breakpoints.md)
  data set), but can be set with the [package
  option](https://amr-for-r.org/release/2.1.1/reference/AMR-options.md)
  [`AMR_guideline`](https://amr-for-r.org/release/2.1.1/reference/AMR-options.md).
  Currently supports EUCAST (2011-2023) and CLSI (2011-2023), see
  *Details*.

- uti:

  (Urinary Tract Infection) A vector with
  [logical](https://rdrr.io/r/base/logical.html)s (`TRUE` or `FALSE`) to
  specify whether a UTI specific interpretation from the guideline
  should be chosen. For using `as.sir()` on a
  [data.frame](https://rdrr.io/r/base/data.frame.html), this can also be
  a column containing [logical](https://rdrr.io/r/base/logical.html)s or
  when left blank, the data set will be searched for a column
  'specimen', and rows within this column containing 'urin' (such as
  'urine', 'urina') will be regarded isolates from a UTI. See
  *Examples*.

- conserve_capped_values:

  a [logical](https://rdrr.io/r/base/logical.html) to indicate that MIC
  values starting with `">"` (but not `">="`) must always return "R" ,
  and that MIC values starting with `"<"` (but not `"<="`) must always
  return "S"

- add_intrinsic_resistance:

  *(only useful when using a EUCAST guideline)* a
  [logical](https://rdrr.io/r/base/logical.html) to indicate whether
  intrinsic antibiotic resistance must also be considered for applicable
  bug-drug combinations, meaning that e.g. ampicillin will always return
  "R" in *Klebsiella* species. Determination is based on the
  [intrinsic_resistant](https://amr-for-r.org/release/2.1.1/reference/intrinsic_resistant.md)
  data set, that itself is based on ['EUCAST Expert Rules' and 'EUCAST
  Intrinsic Resistance and Unusual Phenotypes'
  v3.3](https://www.eucast.org/expert_rules_and_expected_phenotypes)
  (2021).

- reference_data:

  a [data.frame](https://rdrr.io/r/base/data.frame.html) to be used for
  interpretation, which defaults to the
  [clinical_breakpoints](https://amr-for-r.org/release/2.1.1/reference/clinical_breakpoints.md)
  data set. Changing this argument allows for using own interpretation
  guidelines. This argument must contain a data set that is equal in
  structure to the
  [clinical_breakpoints](https://amr-for-r.org/release/2.1.1/reference/clinical_breakpoints.md)
  data set (same column names and column types). Please note that the
  `guideline` argument will be ignored when `reference_data` is manually
  set.

- include_screening:

  a [logical](https://rdrr.io/r/base/logical.html) to indicate that
  clinical breakpoints for screening are allowed - the default is
  `FALSE`. Can also be set with the [package
  option](https://amr-for-r.org/release/2.1.1/reference/AMR-options.md)
  [`AMR_include_screening`](https://amr-for-r.org/release/2.1.1/reference/AMR-options.md).

- include_PKPD:

  a [logical](https://rdrr.io/r/base/logical.html) to indicate that
  PK/PD clinical breakpoints must be applied as a last resort - the
  default is `TRUE`. Can also be set with the [package
  option](https://amr-for-r.org/release/2.1.1/reference/AMR-options.md)
  [`AMR_include_PKPD`](https://amr-for-r.org/release/2.1.1/reference/AMR-options.md).

- breakpoint_type:

  the type of breakpoints to use, either "ECOFF", "animal", or "human".
  ECOFF stands for Epidemiological Cut-Off values. The default is
  `"human"`, which can also be set with the [package
  option](https://amr-for-r.org/release/2.1.1/reference/AMR-options.md)
  [`AMR_breakpoint_type`](https://amr-for-r.org/release/2.1.1/reference/AMR-options.md).

- col_mo:

  column name of the names or codes of the microorganisms (see
  [`as.mo()`](https://amr-for-r.org/release/2.1.1/reference/as.mo.md)) -
  the default is the first column of class
  [`mo`](https://amr-for-r.org/release/2.1.1/reference/as.mo.md). Values
  will be coerced using
  [`as.mo()`](https://amr-for-r.org/release/2.1.1/reference/as.mo.md).

- clean:

  a [logical](https://rdrr.io/r/base/logical.html) to indicate whether
  previously stored results should be forgotten after returning the
  'logbook' with results

## Value

Ordered [factor](https://rdrr.io/pkg/data.table/man/fctr.html) with new
class `sir`

## Details

*Note: The clinical breakpoints in this package were validated through
and imported from [WHONET](https://whonet.org) and the public use of
this `AMR` package has been endorsed by CLSI and EUCAST, please see
[clinical_breakpoints](https://amr-for-r.org/release/2.1.1/reference/clinical_breakpoints.md)
for more information.*

### How it Works

The `as.sir()` function works in four ways:

1.  For **cleaning raw / untransformed data**. The data will be cleaned
    to only contain values S, I and R and will try its best to determine
    this with some intelligence. For example, mixed values with SIR
    interpretations and MIC values such as `"<0.25; S"` will be coerced
    to `"S"`. Combined interpretations for multiple test methods (as
    seen in laboratory records) such as `"S; S"` will be coerced to
    `"S"`, but a value like `"S; I"` will return `NA` with a warning
    that the input is unclear.

2.  For **interpreting minimum inhibitory concentration (MIC) values**
    according to EUCAST or CLSI. You must clean your MIC values first
    using
    [`as.mic()`](https://amr-for-r.org/release/2.1.1/reference/as.mic.md),
    that also gives your columns the new data class
    [`mic`](https://amr-for-r.org/release/2.1.1/reference/as.mic.md).
    Also, be sure to have a column with microorganism names or codes. It
    will be found automatically, but can be set manually using the `mo`
    argument.

    - Using `dplyr`, SIR interpretation can be done very easily with
      either:

          your_data %>% mutate_if(is.mic, as.sir)
          your_data %>% mutate(across(where(is.mic), as.sir))

    - Operators like "\<=" will be stripped before interpretation. When
      using `conserve_capped_values = TRUE`, an MIC value of e.g. "\>2"
      will always return "R", even if the breakpoint according to the
      chosen guideline is "\>=4". This is to prevent that capped values
      from raw laboratory data would not be treated conservatively. The
      default behaviour (`conserve_capped_values = FALSE`) considers
      "\>2" to be lower than "\>=4" and might in this case return "S" or
      "I".

3.  For **interpreting disk diffusion diameters** according to EUCAST or
    CLSI. You must clean your disk zones first using
    [`as.disk()`](https://amr-for-r.org/release/2.1.1/reference/as.disk.md),
    that also gives your columns the new data class
    [`disk`](https://amr-for-r.org/release/2.1.1/reference/as.disk.md).
    Also, be sure to have a column with microorganism names or codes. It
    will be found automatically, but can be set manually using the `mo`
    argument.

    - Using `dplyr`, SIR interpretation can be done very easily with
      either:

          your_data %>% mutate_if(is.disk, as.sir)
          your_data %>% mutate(across(where(is.disk), as.sir))

4.  For **interpreting a complete data set**, with automatic
    determination of MIC values, disk diffusion diameters, microorganism
    names or codes, and antimicrobial test results. This is done very
    simply by running `as.sir(your_data)`.

**For points 2, 3 and 4: Use `sir_interpretation_history()`** to
retrieve a [data.frame](https://rdrr.io/r/base/data.frame.html) (or
[tibble](https://tibble.tidyverse.org/reference/tibble.html) if the
`tibble` package is installed) with all results of the last `as.sir()`
call.

### Supported Guidelines

For interpreting MIC values as well as disk diffusion diameters,
currently implemented guidelines are EUCAST (2011-2023) and CLSI
(2011-2023).

Thus, the `guideline` argument must be set to e.g., `"EUCAST 2023"` or
`"CLSI 2023"`. By simply using `"EUCAST"` (the default) or `"CLSI"` as
input, the latest included version of that guideline will automatically
be selected. You can set your own data set using the `reference_data`
argument. The `guideline` argument will then be ignored.

You can set the default guideline with the [package
option](https://amr-for-r.org/release/2.1.1/reference/AMR-options.md)
[`AMR_guideline`](https://amr-for-r.org/release/2.1.1/reference/AMR-options.md)
(e.g. in your `.Rprofile` file), such as:

      options(AMR_guideline = "CLSI")
      options(AMR_guideline = "CLSI 2018")
      options(AMR_guideline = "EUCAST 2020")
      # or to reset:
      options(AMR_guideline = NULL)

### After Interpretation

After using `as.sir()`, you can use the
[`eucast_rules()`](https://amr-for-r.org/release/2.1.1/reference/eucast_rules.md)
defined by EUCAST to (1) apply inferred susceptibility and resistance
based on results of other antimicrobials and (2) apply intrinsic
resistance based on taxonomic properties of a microorganism.

### Machine-Readable Clinical Breakpoints

The repository of this package [contains a machine-readable
version](https://github.com/msberends/AMR/blob/main/data-raw/clinical_breakpoints.txt)
of all guidelines. This is a CSV file consisting of 29 747 rows and 12
columns. This file is machine-readable, since it contains one row for
every unique combination of the test method (MIC or disk diffusion), the
antimicrobial drug and the microorganism. **This allows for easy
implementation of these rules in laboratory information systems (LIS)**.
Note that it only contains interpretation guidelines for humans -
interpretation guidelines from CLSI for animals were removed.

### Other

The function `is.sir()` detects if the input contains class `sir`. If
the input is a [data.frame](https://rdrr.io/r/base/data.frame.html), it
iterates over all columns and returns a
[logical](https://rdrr.io/r/base/logical.html) vector.

The function `is_sir_eligible()` returns `TRUE` when a columns contains
at most 5% invalid antimicrobial interpretations (not S and/or I and/or
R), and `FALSE` otherwise. The threshold of 5% can be set with the
`threshold` argument. If the input is a
[data.frame](https://rdrr.io/r/base/data.frame.html), it iterates over
all columns and returns a [logical](https://rdrr.io/r/base/logical.html)
vector.

`NA_sir_` is a missing value of the new `sir` class, analogous to e.g.
base R's [`NA_character_`](https://rdrr.io/r/base/NA.html).

## Interpretation of SIR

In 2019, the European Committee on Antimicrobial Susceptibility Testing
(EUCAST) has decided to change the definitions of susceptibility testing
categories S, I, and R as shown below
(<https://www.eucast.org/newsiandr>):

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
[`susceptibility()`](https://amr-for-r.org/release/2.1.1/reference/proportion.md)
(equal to
[`proportion_SI()`](https://amr-for-r.org/release/2.1.1/reference/proportion.md))
to determine antimicrobial susceptibility and
[`count_susceptible()`](https://amr-for-r.org/release/2.1.1/reference/count.md)
(equal to
[`count_SI()`](https://amr-for-r.org/release/2.1.1/reference/count.md))
to count susceptible isolates.

## Reference Data Publicly Available

All data sets in this `AMR` package (about microorganisms, antibiotics,
SIR interpretation, EUCAST rules, etc.) are publicly and freely
available for download in the following formats: R, MS Excel, Apache
Feather, Apache Parquet, SPSS, SAS, and Stata. We also provide
tab-separated plain text files that are machine-readable and suitable
for input in any software program, such as laboratory information
systems. Please visit [our website for the download
links](https://msberends.github.io/AMR/articles/datasets.html). The
actual files are of course available on [our GitHub
repository](https://github.com/msberends/AMR/tree/main/data-raw).

## See also

[`as.mic()`](https://amr-for-r.org/release/2.1.1/reference/as.mic.md),
[`as.disk()`](https://amr-for-r.org/release/2.1.1/reference/as.disk.md),
[`as.mo()`](https://amr-for-r.org/release/2.1.1/reference/as.mo.md)

## Examples

``` r
example_isolates
summary(example_isolates) # see all SIR results at a glance

# For INTERPRETING disk diffusion and MIC values -----------------------

# a whole data set, even with combined MIC values and disk zones
df <- data.frame(
  microorganism = "Escherichia coli",
  AMP = as.mic(8),
  CIP = as.mic(0.256),
  GEN = as.disk(18),
  TOB = as.disk(16),
  ERY = "R"
)
as.sir(df)

# return a 'logbook' about the results:
sir_interpretation_history()

# for single values
as.sir(
  x = as.mic(2),
  mo = as.mo("S. pneumoniae"),
  ab = "AMP",
  guideline = "EUCAST"
)

as.sir(
  x = as.disk(18),
  mo = "Strep pneu", # `mo` will be coerced with as.mo()
  ab = "ampicillin", # and `ab` with as.ab()
  guideline = "EUCAST"
)

# \donttest{
# the dplyr way
if (require("dplyr")) {
  df %>% mutate_if(is.mic, as.sir)
  df %>% mutate_if(function(x) is.mic(x) | is.disk(x), as.sir)
  df %>% mutate(across(where(is.mic), as.sir))
  df %>% mutate_at(vars(AMP:TOB), as.sir)
  df %>% mutate(across(AMP:TOB, as.sir))

  df %>%
    mutate_at(vars(AMP:TOB), as.sir, mo = .$microorganism)

  # to include information about urinary tract infections (UTI)
  data.frame(
    mo = "E. coli",
    NIT = c("<= 2", 32),
    from_the_bladder = c(TRUE, FALSE)
  ) %>%
    as.sir(uti = "from_the_bladder")

  data.frame(
    mo = "E. coli",
    NIT = c("<= 2", 32),
    specimen = c("urine", "blood")
  ) %>%
    as.sir() # automatically determines urine isolates

  df %>%
    mutate_at(vars(AMP:TOB), as.sir, mo = "E. coli", uti = TRUE)
}

# For CLEANING existing SIR values ------------------------------------

as.sir(c("S", "I", "R", "A", "B", "C"))
as.sir("<= 0.002; S") # will return "S"
sir_data <- as.sir(c(rep("S", 474), rep("I", 36), rep("R", 370)))
is.sir(sir_data)
plot(sir_data) # for percentages
barplot(sir_data) # for frequencies

# the dplyr way
if (require("dplyr")) {
  example_isolates %>%
    mutate_at(vars(PEN:RIF), as.sir)
  # same:
  example_isolates %>%
    as.sir(PEN:RIF)

  # fastest way to transform all columns with already valid AMR results to class `sir`:
  example_isolates %>%
    mutate_if(is_sir_eligible, as.sir)

  # since dplyr 1.0.0, this can also be:
  # example_isolates %>%
  #   mutate(across(where(is_sir_eligible), as.sir))
}
# }
```
