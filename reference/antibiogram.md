# Generate Antibiograms (WISCA, Traditional, Combination, or Syndromic)

Generate antibiograms from antimicrobial susceptibility data, with
support for traditional, combination, syndromic, and WISCA
(Weighted-Incidence Syndromic Combination Antibiogram) methods.

**For empirical therapy guidance, WISCA is the recommended approach.**
When initiating empirical treatment, the causative pathogen is unknown,
and the clinically relevant question is: *"what is the probability that
this regimen will cover whatever pathogen turns out to cause the
infection?"* WISCA answers that question directly by weighting
susceptibility by pathogen incidence within a syndrome and providing
credible intervals via Bayesian Monte Carlo simulation. Traditional
antibiograms remain appropriate for tracking resistance per species for
surveillance purposes. See the section *Explaining WISCA* on this page
and the [WISCA vignette](https://amr-for-r.org/articles/WISCA.html) for
details.

All antibiogram types adhere to previously described approaches (see
*Source*), and the WISCA method implements the Bayesian decision model
by Bielicki *et al.* (2016,
[doi:10.1093/jac/dkv397](https://doi.org/10.1093/jac/dkv397) ). Output
formats include plots and tables, ideal for integration with R Markdown
and Quarto reports.

## Usage

``` r
wisca(
  x,
  antimicrobials = where(is.sir),
  ab_transform = "name",
  syndromic_group = NULL,
  only_all_tested = FALSE,
  digits = 1,
  formatting_type = getOption("AMR_antibiogram_formatting_type", 14),
  col_mo = NULL,
  language = get_AMR_locale(),
  combine_SI,
  guideline = getOption("AMR_guideline", "EUCAST"),
  sep = " + ",
  sort_columns = TRUE,
  simulations = 1000,
  conf_interval = 0.95,
  interval_side = "two-tailed",
  info = interactive(),
  parallel = FALSE,
  ...
)

antibiogram(
  x,
  antimicrobials = where(is.sir),
  mo_transform = "shortname",
  ab_transform = "name",
  syndromic_group = NULL,
  add_total_n = FALSE,
  only_all_tested = FALSE,
  digits = ifelse(wisca, 1, 0),
  formatting_type = getOption("AMR_antibiogram_formatting_type", ifelse(wisca, 14, 18)),
  col_mo = NULL,
  language = get_AMR_locale(),
  minimum = 30,
  combine_SI,
  guideline = getOption("AMR_guideline", "EUCAST"),
  sep = " + ",
  sort_columns = TRUE,
  wisca = FALSE,
  simulations = 1000,
  conf_interval = 0.95,
  interval_side = "two-tailed",
  info = interactive(),
  parallel = FALSE,
  ...
)

retrieve_wisca_parameters(wisca_model, ...)

# S3 method for class 'antibiogram'
plot(x, ...)

# S3 method for class 'antibiogram'
autoplot(
  object,
  geom = c("pointrange", "point", "col", "bar", "errorbar"),
  ci = TRUE,
  sort = TRUE,
  flip = NULL,
  caption = NULL,
  ...
)

wisca_plot(
  wisca_model,
  wisca_plot_type = c("susceptibility_incidence", "posterior_coverage"),
  ...
)

# S3 method for class 'antibiogram'
knit_print(
  x,
  italicise = TRUE,
  na = getOption("knitr.kable.NA", default = ""),
  ...
)
```

## Arguments

- x:

  A [data.frame](https://rdrr.io/r/base/data.frame.html) containing at
  least a column with microorganisms and columns with antimicrobial
  results (class 'sir', see
  [`as.sir()`](https://amr-for-r.org/reference/as.sir.md)).

- antimicrobials:

  A vector specifying the antimicrobials containing SIR values to
  include in the antibiogram (see *Examples*). Will be evaluated using
  [`guess_ab_col()`](https://amr-for-r.org/reference/guess_ab_col.md).
  This can be:

  - Any antimicrobial name or code that could match (see
    [`guess_ab_col()`](https://amr-for-r.org/reference/guess_ab_col.md))
    to any column in `x`

  - Any [antimicrobial
    selector](https://amr-for-r.org/reference/antimicrobial_selectors.md),
    such as
    [`aminoglycosides()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
    or
    [`carbapenems()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)

  - A combination of the above, using
    [`c()`](https://rdrr.io/r/base/c.html), e.g.:

    - `c(aminoglycosides(), "AMP", "AMC")`

    - `c(aminoglycosides(), carbapenems())`

  - Column indices using numbers

  - Combination therapy, indicated by using `"+"`, with or without
    [antimicrobial
    selectors](https://amr-for-r.org/reference/antimicrobial_selectors.md),
    e.g.:

    - `"cipro + genta"`

    - `"TZP+TOB"`

    - `c("TZP", "TZP+GEN", "TZP+TOB")`

    - `carbapenems() + "GEN"`

    - `carbapenems() + c("", "GEN")`

    - `carbapenems() + c("", aminoglycosides())`

- ab_transform:

  A character to transform antimicrobial input - must be one of the
  column names of the
  [antimicrobials](https://amr-for-r.org/reference/antimicrobials.md)
  data set (defaults to `"name"`): `"ab"`, `"cid"`, `"name"`, `"group"`,
  `"atc"`, `"atc_group1"`, `"atc_group2"`, `"abbreviations"`,
  `"synonyms"`, `"oral_ddd"`, `"oral_units"`, `"iv_ddd"`, `"iv_units"`,
  or `"loinc"`. Can also be `NULL` to not transform the input.

- syndromic_group:

  A column name of `x`, or values calculated to split rows of `x`, e.g.
  by using [`ifelse()`](https://rdrr.io/r/base/ifelse.html) or
  [`case_when()`](https://dplyr.tidyverse.org/reference/case-and-replace-when.html).
  See *Examples*.

- only_all_tested:

  (for combination antibiograms): a
  [logical](https://rdrr.io/r/base/logical.html) to indicate that
  isolates must be tested for all antimicrobials, see *Details*.

- digits:

  Number of digits to use for rounding the antimicrobial coverage,
  defaults to 1 for WISCA and 0 otherwise.

- formatting_type:

  Numeric value (1-22 for WISCA, 1-12 for non-WISCA) indicating how the
  'cells' of the antibiogram table should be formatted. See *Details* \>
  *Formatting Type* for a list of options.

- col_mo:

  Column name of the names or codes of the microorganisms (see
  [`as.mo()`](https://amr-for-r.org/reference/as.mo.md)) - the default
  is the first column of class
  [`mo`](https://amr-for-r.org/reference/as.mo.md). Values will be
  coerced using [`as.mo()`](https://amr-for-r.org/reference/as.mo.md).

- language:

  Language to translate text, which defaults to the system language (see
  [`get_AMR_locale()`](https://amr-for-r.org/reference/translate.md)).

- combine_SI:

  *(deprecated in favour of `guideline`)* A
  [logical](https://rdrr.io/r/base/logical.html). If set, it takes
  precedence over `guideline`: `TRUE` considers S, SDD, and I
  susceptible, and `FALSE` only S.

- guideline:

  A guideline name, `"EUCAST"` or `"CLSI"` (or any version thereof, e.g.
  `"CLSI 2025"`), to determine whether the 'I' category is considered
  susceptible (EUCAST) or resistant (CLSI). 'SDD' is always considered
  susceptible. The default is the package option
  [`AMR_guideline`](https://amr-for-r.org/reference/AMR-options.md), or
  `"EUCAST"` if it is not set. See section *Interpretation of SIR*
  below.

- sep:

  A separating character for antimicrobial columns in combination
  antibiograms.

- sort_columns:

  A [logical](https://rdrr.io/r/base/logical.html) to indicate whether
  the antimicrobial columns must be sorted on name.

- simulations:

  (for WISCA) a numerical value to set the number of Monte Carlo
  simulations.

- conf_interval:

  A numerical value to set confidence interval (default is `0.95`).

- interval_side:

  The side of the confidence interval, either `"two-tailed"` (default),
  `"left"` or `"right"`.

- info:

  A [logical](https://rdrr.io/r/base/logical.html) to indicate info
  should be printed - the default is `TRUE` only in interactive mode.

- parallel:

  A [logical](https://rdrr.io/r/base/logical.html) to indicate if
  parallel computing must be used, defaults to `FALSE`. Requires the
  [`future.apply`](https://future.apply.futureverse.org/reference/future_lapply.html)
  package. For WISCA, Monte Carlo simulations are distributed across
  workers; for grouped antibiograms, each group is processed by a
  separate worker. **A non-sequential
  [`future::plan()`](https://future.futureverse.org/reference/plan.html)
  must already be active before setting `parallel = TRUE`** – for
  example, `future::plan(future::multisession)`. An error is thrown if
  `parallel = TRUE` is used without a plan set by the user.

- ...:

  Currently unused.

- mo_transform:

  A character to transform microorganism input - must be `"name"`,
  `"shortname"` (default), `"gramstain"`, or one of the column names of
  the
  [microorganisms](https://amr-for-r.org/reference/microorganisms.md)
  data set: `"mo"`, `"fullname"`, `"status"`, `"domain"`, `"kingdom"`,
  `"phylum"`, `"class"`, `"order"`, `"family"`, `"genus"`, `"species"`,
  `"subspecies"`, `"rank"`, `"ref"`, `"oxygen_tolerance"`,
  `"morphology"`, `"source"`, `"lpsn"`, `"lpsn_parent"`,
  `"lpsn_renamed_to"`, `"mycobank"`, `"mycobank_parent"`,
  `"mycobank_renamed_to"`, `"gbif"`, `"gbif_parent"`,
  `"gbif_renamed_to"`, `"prevalence"`, or `"snomed"`. Can also be `NULL`
  to not transform the input or `NA` to consider all microorganisms
  'unknown'.

- add_total_n:

  *(deprecated in favour of `formatting_type`)* A
  [logical](https://rdrr.io/r/base/logical.html) to indicate whether
  `n_tested` available numbers per pathogen should be added to the table
  (default is `TRUE`). This will add the lowest and highest number of
  available isolates per antimicrobial (e.g., if for *E. coli* 200
  isolates are available for ciprofloxacin and 150 for amoxicillin, the
  returned number will be "150-200"). This option is unavailable when
  `wisca = TRUE`; in that case, use `retrieve_wisca_parameters()` to get
  the parameters used for WISCA.

- minimum:

  The minimum allowed number of available (tested) isolates. Any isolate
  count lower than `minimum` will return `NA` with a warning. The
  default number of `30` isolates is advised by the Clinical and
  Laboratory Standards Institute (CLSI) as best practice, see *Source*.

- wisca:

  A [logical](https://rdrr.io/r/base/logical.html) to indicate whether a
  Weighted-Incidence Syndromic Combination Antibiogram (WISCA) must be
  generated (default is `FALSE`). This will use a Bayesian decision
  model to estimate regimen coverage probabilities using [Monte Carlo
  simulations](https://en.wikipedia.org/wiki/Monte_Carlo_method). Per
  [doi:10.1093/jac/dkv397](https://doi.org/10.1093/jac/dkv397) ,
  susceptibility priors are \\\beta(0.5, 0.5)\\ (Jeffreys) and
  intrinsically resistant pairs (based on
  [intrinsic_resistant](https://amr-for-r.org/reference/intrinsic_resistant.md))
  use \\\beta(1, 9999)\\.

  Set `simulations`, `conf_interval`, and `interval_side` to adjust.

- wisca_model:

  The outcome of `wisca()` or `antibiogram(..., wisca = TRUE)`.

- object:

  An `antibiogram()` object.

- geom:

  The plotting style for the point estimate. One of `"pointrange"`
  (default), `"point"`, `"col"`/`"bar"`, or `"errorbar"`. `"pointrange"`
  is recommended for coverage data: bars imply a meaningful baseline at
  zero, which coverage estimates rarely have.

- ci:

  Logical, whether to draw the credible/confidence interval. Defaults to
  `TRUE`. Ignored (forced `TRUE`) when `geom = "pointrange"` or
  `"errorbar"`, since the interval is intrinsic to those geoms.

- sort:

  Logical, whether to order regimens by coverage. Defaults to `TRUE`.
  When faceted (per pathogen) or grouped (syndromic), ordering is
  applied within each panel/group.

- flip:

  Logical, whether to draw regimens on the y-axis (horizontal). Defaults
  to `NULL`, which flips automatically when any regimen label exceeds 20
  characters (long combination names read poorly on the x-axis). Set
  `TRUE`/`FALSE` to override.

- caption:

  Text to show as caption, will explain non-inferiority for WISCA
  models.

- wisca_plot_type:

  Either `"susceptibility_incidence"` (default) or
  `"posterior_coverage"`.

- italicise:

  A [logical](https://rdrr.io/r/base/logical.html) to indicate whether
  the microorganism names in the
  [knitr](https://rdrr.io/pkg/knitr/man/kable.html) table should be made
  italic, using
  [`italicise_taxonomy()`](https://amr-for-r.org/reference/italicise_taxonomy.md).

- na:

  Character to use for showing `NA` values.

## Details

These functions return a table with values between 0 and 100 for
*susceptibility*, not resistance.

**Remember that you should filter your data to let it contain only first
isolates!** This is needed to exclude duplicates and to reduce selection
bias. Use
[`first_isolate()`](https://amr-for-r.org/reference/first_isolate.md) to
determine them with one of the four available algorithms: isolate-based,
patient-based, episode-based, or phenotype-based.

For estimating antimicrobial coverage, especially when creating a WISCA,
the outcome might become more reliable by only including the top *n*
species encountered in the data. You can filter on this top *n* using
[`top_n_microorganisms()`](https://amr-for-r.org/reference/top_n_microorganisms.md).
For example, use `top_n_microorganisms(your_data, n = 10)` as a
pre-processing step to only include the top 10 species in the data.

The numeric values of an antibiogram are stored in a long format as the
[attribute](https://rdrr.io/r/base/attributes.html) `long_numeric`. You
can retrieve them using `attributes(x)$long_numeric`, where `x` is the
outcome of `antibiogram()` or `wisca()`. This is ideal for e.g. advanced
plotting.

### Formatting Type

The formatting of the 'cells' of the table can be set with the argument
`formatting_type`. In these examples, `5` indicates the antimicrobial
coverage (`4-6` the confidence level), `15` the number of susceptible
isolates, and `300` the number of tested (i.e., available) isolates:

1.  5

2.  15

3.  300

4.  15/300

5.  5 (300)

6.  5% (300)

7.  5 (N=300)

8.  5% (N=300)

9.  5 (15/300)

10. 5% (15/300)

11. 5 (N=15/300)

12. 5% (N=15/300)

13. 5 (4-6)

14. 5% (4-6%) - **default for WISCA**

15. 5 (4-6,300)

16. 5% (4-6%,300)

17. 5 (4-6,N=300)

18. 5% (4-6%,N=300) - **default for non-WISCA**

19. 5 (4-6,15/300)

20. 5% (4-6%,15/300)

21. 5 (4-6,N=15/300)

22. 5% (4-6%,N=15/300)

The default can be set globally with the package option
[`AMR_antibiogram_formatting_type`](https://amr-for-r.org/reference/AMR-options.md),
e.g. `options(AMR_antibiogram_formatting_type = 5)`. Do note that for
WISCA, the total numbers of tested and susceptible isolates are less
useful to report, since these are included in the Bayesian model and
apparent from the susceptibility and its confidence level.

Set `digits` (defaults to `0`) to alter the rounding of the
susceptibility percentages.

### When to Use WISCA vs. Traditional Antibiograms

There are various antibiogram types, as summarised by Klinker *et al.*
(2021,
[doi:10.1177/20499361211011373](https://doi.org/10.1177/20499361211011373)
), and they are all supported by `antibiogram()`: traditional,
combination, syndromic, and WISCA.

**If your goal is to guide empirical therapy, use WISCA.** Traditional
antibiograms fragment susceptibility information by species, but at the
point of prescribing, the clinician does not know which species is
causing the infection. WISCA shifts the unit of analysis from the
isolate to the patient: it estimates the probability that a regimen will
cover the infection, given the local distribution of causative
pathogens. It evaluates combination regimens, weights by pathogen
incidence, and provides credible intervals that honestly communicate
uncertainty. Hebert *et al.* (2012) demonstrated this concretely for the
first time: ciprofloxacin showed 84% susceptibility against *E. coli* in
the traditional antibiogram, but WISCA coverage was only 62% for UTI and
37% for abdominal infections, because other species (including
intrinsically resistant enterococci) contribute substantially to these
syndromes. Note that WISCA is pathogen-agnostic: the outcome is not
stratified by species, but by syndrome.

**Traditional, combination, and syndromic antibiograms remain
appropriate for AMR surveillance**, i.e., tracking resistance trends per
species over time. They are the right tool when the question is *"how
resistant is species X to drug Y in our setting?"* rather than *"what
regimen best covers this syndrome?"*.

All four types are demonstrated in the *Examples* section below.

### Grouped tibbles

For any type of antibiogram, grouped
[tibbles](https://tibble.tidyverse.org/reference/tibble.html) can also
be used to calculate susceptibilities over various groups.

Code example:

    library(dplyr)
    your_data %>%
      group_by(has_sepsis, is_neonate, sex) %>%
      wisca(antimicrobials = c("TZP", "TZP+TOB", "TZP+GEN"))

### Inclusion in Combination Antibiograms

For combination antibiograms, the results of the antimicrobials are
combined using the same rules as
[`sir_combine()`](https://amr-for-r.org/reference/sir_combine.md), see
the section *Combination Therapy* below. Whether isolates must be tested
for all antimicrobials can be set with the `only_all_tested` argument
(default is `FALSE`). Whether the 'I' category is considered susceptible
follows the `guideline` argument, which also determines whether a
missing result can be ignored. To use a combination in other analyses as
well, create a new variable with
[`sir_combine()`](https://amr-for-r.org/reference/sir_combine.md).

### Plotting

All types of antibiograms as listed above can be plotted (using
[`ggplot2::autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
or base R's [`plot()`](https://amr-for-r.org/reference/plot.md) and
[`barplot()`](https://rdrr.io/r/graphics/barplot.html)). As mentioned
above, the numeric values of an antibiogram are stored in a long format
as the [attribute](https://rdrr.io/r/base/attributes.html)
`long_numeric`. You can retrieve them using
`attributes(x)$long_numeric`, where `x` is the outcome of
`antibiogram()` or `wisca()`.

The outcome of `antibiogram()` can also be used directly in R Markdown /
Quarto (i.e., `knitr`) for reports. In this case,
[`knitr::kable()`](https://rdrr.io/pkg/knitr/man/kable.html) will be
applied automatically and microorganism names will even be printed in
italics at default (see argument `italicise`).

You can also use functions from specific 'table reporting' packages to
transform the output of `antibiogram()` to your needs, e.g. with
`flextable::as_flextable()` or `gt::gt()`.

## Explaining WISCA

WISCA (Weighted-Incidence Syndromic Combination Antibiogram) estimates
the probability that an empirical antimicrobial regimen will provide
adequate coverage for a given infection syndrome, before the causative
pathogen has been identified.

It does so by combining two quantities: the relative incidence of each
pathogen within the syndrome (modelled as a Dirichlet distribution) and
the susceptibility of each pathogen to the regimen (modelled as Beta
distributions). These are combined via Monte Carlo simulation to produce
a coverage estimate with a credible interval.

**Prior distributions:** Pathogen incidence uses a non-informative
\\Dirichlet(1, 1, \ldots, 1)\\ prior. Susceptibility proportions use the
Jeffreys prior, \\\beta(0.5, 0.5)\\, except for pathogen-drug
combinations with known intrinsic resistance, which use a strongly
informative \\\beta(1, 9999)\\ prior that forces near-zero
susceptibility regardless of observed data. Intrinsic resistance is
determined using the
[intrinsic_resistant](https://amr-for-r.org/reference/intrinsic_resistant.md)
data set, which is based on ['EUCAST Expected Resistant Phenotypes'
v1.2](https://www.eucast.org/bacteria/important-additional-information/expert-rules/)
(2023).

**Pathogens included in the model:** For each regimen (and group), only
pathogens that occur in the data and have at least one test result for
the regimen are included. Pathogens without any test result for a
regimen do not contribute to its coverage, nor to its `n_total` in
`long_numeric`. Use `retrieve_wisca_parameters()` to get the parameters
that were used.

**Rare pathogens and small groups:** Since every pathogen receives a
prior weight of 1 in the Dirichlet distribution, the prior weighs
heavily when there are few isolates relative to the number of pathogens,
i.e., with pathogens found only once or a few times, or in small groups.
This pulls the incidence towards a uniform distribution, and the very
uncertain susceptibility of rare pathogens can then strongly lower the
estimated coverage, without reflecting the clinical expectation. A
warning is therefore given when this prior weight accounts for more than
5% of the pathogen incidence in any regimen (or group). Consider
limiting the analysis to the most common pathogens using
[`top_n_microorganisms()`](https://amr-for-r.org/reference/top_n_microorganisms.md),
or using fewer and larger groups.

**Interpreting the output:** Overlapping credible intervals between
regimens indicate no significant difference in coverage; if a
narrower-spectrum regimen overlaps with a broader one, the
narrower-spectrum option may be preferred on stewardship grounds.
Non-overlapping intervals indicate a clinically meaningful difference.
For small sample sizes, consider pooling data from multiple sites to
improve precision, provided pathogen distributions are sufficiently
similar (Bielicki *et al.*, 2016).

For the full mathematical derivation and worked examples, see the [WISCA
vignette](https://amr-for-r.org/articles/WISCA.html).

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
[`proportion_*()`](https://amr-for-r.org/reference/proportion.md)
functions apply these rules when more than one variable is used. With
`only_all_tested = FALSE` (default), the functions for susceptibility,
such as
[`susceptibility()`](https://amr-for-r.org/reference/proportion.md) and
[`count_susceptible()`](https://amr-for-r.org/reference/count.md), count
an isolate as susceptible if any of the antimicrobials is susceptible
(`"most_susceptible"`), and the functions for resistance, such as
[`resistance()`](https://amr-for-r.org/reference/proportion.md) and
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

In `antibiogram()` and `wisca()`, the result of the combination therapy
(`"most_susceptible"`) is used, the 'I' category is considered
susceptible or not following the `guideline` argument ('SDD' is always
considered susceptible), and the result is then reduced to S or R.

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
  [`proportion_df()`](https://amr-for-r.org/reference/proportion.md).

- **`"CLSI"`**, or any CLSI version such as `"CLSI 2025"`, to consider
  'I' resistant, e.g., to return `"S"` vs. `"IR"` in
  [`proportion_df()`](https://amr-for-r.org/reference/proportion.md).

- **`NULL` or `FALSE`**, to keep S, I, and R separate. This is only
  possible in functions that return results per category:
  [`count_df()`](https://amr-for-r.org/reference/count.md),
  [`proportion_df()`](https://amr-for-r.org/reference/proportion.md),
  [`sir_df()`](https://amr-for-r.org/reference/proportion.md),
  [`ggplot_sir()`](https://amr-for-r.org/reference/ggplot_sir.md),
  [`geom_sir()`](https://amr-for-r.org/reference/ggplot_sir.md),
  [`labels_sir_count()`](https://amr-for-r.org/reference/ggplot_sir.md),
  and
  [`bug_drug_combinations()`](https://amr-for-r.org/reference/bug_drug_combinations.md).
  Functions that return a single proportion or count, such as
  [`susceptibility()`](https://amr-for-r.org/reference/proportion.md),
  [`count_resistant()`](https://amr-for-r.org/reference/count.md),
  `antibiogram()`, `wisca()`, and
  [`format()`](https://rdrr.io/r/base/format.html) on the result of
  [`bug_drug_combinations()`](https://amr-for-r.org/reference/bug_drug_combinations.md),
  require 'I' to be considered either susceptible or resistant, and
  return an error.

Use [`susceptibility()`](https://amr-for-r.org/reference/proportion.md)
and [`resistance()`](https://amr-for-r.org/reference/proportion.md) (or
[`count_susceptible()`](https://amr-for-r.org/reference/count.md) and
[`count_resistant()`](https://amr-for-r.org/reference/count.md)) to
follow the guideline, or e.g.
[`proportion_SI()`](https://amr-for-r.org/reference/proportion.md) and
[`proportion_R()`](https://amr-for-r.org/reference/proportion.md) to be
explicit regardless of the guideline.

The `combine_SI` argument of the functions that return results per
category, `antibiogram()`, and `wisca()` is deprecated in favour of
`guideline`. If set, it is followed literally and takes precedence over
`guideline`.

## References

- Hebert C *et al.* (2012). **Demonstration of the weighted-incidence
  syndromic combination antibiogram: an empiric prescribing decision
  aid.** *Infection Control & Hospital Epidemiology* 33(4):381-388;
  [doi:10.1086/664768](https://doi.org/10.1086/664768)

- Bielicki JA *et al.* (2016). **Selecting appropriate empirical
  antibiotic regimens for paediatric bloodstream infections: application
  of a Bayesian decision model to local and pooled antimicrobial
  resistance surveillance data.** *Journal of Antimicrobial
  Chemotherapy* 71(3):794-802;
  [doi:10.1093/jac/dkv397](https://doi.org/10.1093/jac/dkv397)

- Cook A *et al.* (2022). **Improving empiric antibiotic prescribing in
  pediatric bloodstream infections: a potential application of
  weighted-incidence syndromic combination antibiograms (WISCA).**
  *Expert Review of Anti-infective Therapy* 20(3):445-456;
  [doi:10.1080/14787210.2021.1967145](https://doi.org/10.1080/14787210.2021.1967145)

- Klinker KP *et al.* (2021). **Antimicrobial stewardship and
  antibiograms: importance of moving beyond traditional antibiograms.**
  *Therapeutic Advances in Infectious Disease*, May
  5;8:20499361211011373;
  [doi:10.1177/20499361211011373](https://doi.org/10.1177/20499361211011373)

- Barbieri E *et al.* (2021). **Development of a Weighted-Incidence
  Syndromic Combination Antibiogram (WISCA) to guide the choice of the
  empiric antibiotic treatment for urinary tract infection in paediatric
  patients: a Bayesian approach.** *Antimicrobial Resistance & Infection
  Control* May 1;10(1):74;
  [doi:10.1186/s13756-021-00939-2](https://doi.org/10.1186/s13756-021-00939-2)

- **M39 Analysis and Presentation of Cumulative Antimicrobial
  Susceptibility Test Data, 5th Edition**, 2022, *Clinical and
  Laboratory Standards Institute (CLSI)*.
  <https://clsi.org/standards/products/microbiology/documents/m39/>.

## Author

Implementation: Dr. Larisse Bolton and Dr. Matthijs Berends

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

# \donttest{
# WISCA antibiogram (recommended for empirical therapy) -----------------

# basic WISCA: empirical coverage per regimen, weighted by pathogen
# incidence, with 95% credible intervals
# (10 most common species, since rare species can distort coverage estimates)
wisca(top_n_microorganisms(example_isolates, n = 10),
  antimicrobials = c("AMC", "AMC+CIP", "AMC+GEN")
)
#> ℹ `wisca()` assumes the EUCAST guideline and thus considers the 'I' category
#>   susceptible, so S, SDD, and I count as covered. Use `guideline = "CLSI"` to
#>   only count S and SDD as covered. Set the `guideline` argument, or set the
#>   `AMR_guideline` option to change the default (see `?AMR-options`). This
#>   message is shown once per session.
#> ℹ Using column mo as input for `col_mo`.
#> # A tibble:       1 × 3
#> # Type:           Weighted-Incidence Syndromic Combination Antibiogram (WISCA)
#> # Cred. interval: 95%
#> # Simulations:    1000 per stratum
#>   `Amoxicillin/clavulanic acid` Amoxicillin/clavulanic …¹ Amoxicillin/clavulan…²
#>   <chr>                         <chr>                     <chr>                 
#> 1 76.9% (74.7-79%)              89.5% (88-91.1%)          93.6% (92.4-94.8%)    
#> # ℹ abbreviated names: ¹​`Amoxicillin/clavulanic acid + Ciprofloxacin`,
#> #   ²​`Amoxicillin/clavulanic acid + Gentamicin`
#> # Use `ggplot2::autoplot()` or base R `plot()` to create a plot of this antibiogram,
#> # and use `wisca_plot()` to assess the simulation outcomes.
#> # Or, use it directly in R Markdown or Quarto, see `antibiogram()`.

# equivalent using antibiogram():
antibiogram(top_n_microorganisms(example_isolates, n = 10),
  antimicrobials = c("AMC", "AMC+CIP", "AMC+GEN"),
  wisca = TRUE
)
#> # A tibble:       1 × 3
#> # Type:           Weighted-Incidence Syndromic Combination Antibiogram (WISCA)
#> # Cred. interval: 95%
#> # Simulations:    1000 per stratum
#>   `Amoxicillin/clavulanic acid` Amoxicillin/clavulanic …¹ Amoxicillin/clavulan…²
#>   <chr>                         <chr>                     <chr>                 
#> 1 76.8% (74.8-78.9%)            89.6% (88-91.1%)          93.7% (92.4-94.9%)    
#> # ℹ abbreviated names: ¹​`Amoxicillin/clavulanic acid + Ciprofloxacin`,
#> #   ²​`Amoxicillin/clavulanic acid + Gentamicin`
#> # Use `ggplot2::autoplot()` or base R `plot()` to create a plot of this antibiogram,
#> # and use `wisca_plot()` to assess the simulation outcomes.
#> # Or, use it directly in R Markdown or Quarto, see `antibiogram()`.

# stratified by syndrome or clinical group (5 most common species, for a readable plot)
out <- wisca(top_n_microorganisms(example_isolates, n = 5),
  antimicrobials = c("TZP", "TZP+TOB", "TZP+GEN"),
  syndromic_group = "ward"
)
#> Warning: In 2 of 9 WISCA regimen(s) or group(s), the prior weight of 1 per pathogen
#> accounts for more than 5% of the pathogen incidence (up to 5.1%), since there
#> are few isolates relative to the number of pathogens, i.e., rare pathogens or a
#> small group. Because WISCA gives every pathogen a prior weight of 1 (Bielicki
#> et al., 2016), this can strongly lower the estimated coverage without
#> reflecting the clinical expectation. Consider limiting the analysis to the most
#> common pathogens with `top_n_microorganisms()`, or using fewer and larger
#> groups.
out
#> # A tibble:       3 × 4
#> # Type:           Weighted-Incidence Syndromic Combination Antibiogram (WISCA)
#> # Cred. interval: 95%
#> # Simulations:    1000 per stratum
#>   `Syndromic Group` `Piperacillin/tazobactam` Piperacillin/tazobactam + Gentam…¹
#>   <chr>             <chr>                     <chr>                             
#> 1 Clinical          76.1% (69-83.5%)          99% (98.2-99.6%)                  
#> 2 ICU               54% (45.3-66.7%)          98.4% (96.9-99.5%)                
#> 3 Outpatient        68.4% (47.5-89.2%)        96.1% (91-99%)                    
#> # ℹ abbreviated name: ¹​`Piperacillin/tazobactam + Gentamicin`
#> # ℹ 1 more variable: `Piperacillin/tazobactam + Tobramycin` <chr>
#> # Use `ggplot2::autoplot()` or base R `plot()` to create a plot of this antibiogram,
#> # and use `wisca_plot()` to assess the simulation outcomes.
#> # Or, use it directly in R Markdown or Quarto, see `antibiogram()`.
if (requireNamespace("ggplot2")) {
  wisca_plot(out)
}


# stratified using grouped tibbles (e.g. by age and gender)
if (requireNamespace("dplyr")) {
  library(dplyr)
  example_isolates %>%
    top_n_microorganisms(n = 10) %>%
    group_by(
      age_group = age_groups(age, c(25, 50, 75)),
      gender
    ) %>%
    wisca(antimicrobials = c("TZP", "TZP+TOB", "TZP+GEN"))
}
#> Warning: Number of tested isolates should exceed 30 for each regimen (and group). WISCA
#> coverage estimates might be inaccurate.
#> Warning: In 12 of 24 WISCA regimen(s) or group(s), the prior weight of 1 per pathogen
#> accounts for more than 5% of the pathogen incidence (up to 33.3%), since there
#> are few isolates relative to the number of pathogens, i.e., rare pathogens or a
#> small group. Because WISCA gives every pathogen a prior weight of 1 (Bielicki
#> et al., 2016), this can strongly lower the estimated coverage without
#> reflecting the clinical expectation. Consider limiting the analysis to the most
#> common pathogens with `top_n_microorganisms()`, or using fewer and larger
#> groups.
#> # A tibble:       8 × 5
#> # Type:           Weighted-Incidence Syndromic Combination Antibiogram (WISCA)
#> # Cred. interval: 95%
#> # Simulations:    1000 per stratum
#>   age_group gender `Piperacillin/tazobactam` Piperacillin/tazobactam + Gentami…¹
#>   <chr>     <chr>  <chr>                     <chr>                              
#> 1 0-24      F      74.7% (16.9-100%)         81.5% (54.9-97.3%)                 
#> 2 0-24      M      83.5% (50.2-99.1%)        86.7% (66.6-98.3%)                 
#> 3 25-49     F      92.1% (80.9-98.6%)        95.2% (88.9-98.7%)                 
#> 4 25-49     M      80.5% (57.5-95.8%)        93.3% (85.3-98.4%)                 
#> 5 50-74     F      75.6% (65.4-88%)          97% (94.7-98.6%)                   
#> 6 50-74     M      68.1% (56.7-78.9%)        96.7% (94.1-98.6%)                 
#> 7 75+       F      75.1% (63.6-85.5%)        97.8% (96-99%)                     
#> 8 75+       M      74.9% (64.4-84.2%)        97.8% (96.1-99.1%)                 
#> # ℹ abbreviated name: ¹​`Piperacillin/tazobactam + Gentamicin`
#> # ℹ 1 more variable: `Piperacillin/tazobactam + Tobramycin` <chr>
#> # Use `ggplot2::autoplot()` or base R `plot()` to create a plot of this antibiogram,
#> # and use `wisca_plot()` to assess the simulation outcomes.
#> # Or, use it directly in R Markdown or Quarto, see `antibiogram()`.


# Traditional antibiogram (for AMR surveillance) ------------------------

antibiogram(example_isolates,
  antimicrobials = c(aminoglycosides(), carbapenems())
)
#> ℹ For `aminoglycosides()` using columns GEN (gentamicin), TOB (tobramycin), AMK
#>   (amikacin), and KAN (kanamycin)
#> ℹ For `carbapenems()` using columns IPM (imipenem) and MEM (meropenem)
#> ℹ `antibiogram()` assumes the EUCAST guideline and thus considers the 'I'
#>   category susceptible, so S, SDD, and I count as susceptible. Use `guideline =
#>   "CLSI"` to only count S and SDD as susceptible. Set the `guideline` argument,
#>   or set the `AMR_guideline` option to change the default (see `?AMR-options`).
#>   This message is shown once per session.
#> # A tibble:       10 × 7
#> # Type:           Traditional Antibiogram
#> # Conf. interval: 95%
#>    Pathogen       Amikacin    Gentamicin Imipenem Kanamycin Meropenem Tobramycin
#>    <chr>          <chr>       <chr>      <chr>    <chr>     <chr>     <chr>     
#>  1 CoNS           0% (0-8%,N… 86% (82-9… 52% (37… 0% (0-8%… 52% (37-… 22% (12-3…
#>  2 E. coli        100% (98-1… 98% (96-9… 100% (9… NA        100% (99… 97% (96-9…
#>  3 E. faecalis    0% (0-9%,N… 0% (0-9%,… 100% (9… 0% (0-9%… NA        0% (0-9%,…
#>  4 K. pneumoniae  NA          90% (79-9… 100% (9… NA        100% (93… 90% (79-9…
#>  5 P. aeruginosa  NA          100% (88-… NA       0% (0-12… NA        100% (88-…
#>  6 P. mirabilis   NA          94% (80-9… 94% (79… NA        NA        94% (80-9…
#>  7 S. aureus      NA          99% (97-1… NA       NA        NA        98% (92-1…
#>  8 S. epidermidis 0% (0-8%,N… 79% (71-8… NA       0% (0-8%… NA        51% (40-6…
#>  9 S. hominis     NA          92% (84-9… NA       NA        NA        85% (74-9…
#> 10 S. pneumoniae  0% (0-3%,N… 0% (0-3%,… NA       0% (0-3%… NA        0% (0-3%,…
#> # Use `ggplot2::autoplot()` or base R `plot()` to create a plot of this antibiogram,
#> # or use it directly in R Markdown or Quarto, see `antibiogram()`.

antibiogram(example_isolates,
  antimicrobials = aminoglycosides(),
  ab_transform = "atc",
  mo_transform = "gramstain"
)
#> ℹ For `aminoglycosides()` using columns GEN (gentamicin), TOB (tobramycin), AMK
#>   (amikacin), and KAN (kanamycin)
#> # A tibble:       2 × 5
#> # Type:           Traditional Antibiogram
#> # Conf. interval: 95%
#>   Pathogen      J01GB01            J01GB03             J01GB04         J01GB06  
#>   <chr>         <chr>              <chr>               <chr>           <chr>    
#> 1 Gram-negative 96% (94-97%,N=686) 96% (95-98%,N=684)  0% (0-10%,N=35) 98% (96-…
#> 2 Gram-positive 34% (31-38%,N=665) 63% (60-66%,N=1170) 0% (0-1%,N=436) 0% (0-1%…
#> # Use `ggplot2::autoplot()` or base R `plot()` to create a plot of this antibiogram,
#> # or use it directly in R Markdown or Quarto, see `antibiogram()`.


# Combination antibiogram (for AMR surveillance) ------------------------

antibiogram(example_isolates,
  antimicrobials = c("TZP", "TZP+TOB", "TZP+GEN"),
  mo_transform = "gramstain"
)
#> # A tibble:       2 × 4
#> # Type:           Combination Antibiogram
#> # Conf. interval: 95%
#>   Pathogen  Piperacillin/tazobac…¹ Piperacillin/tazobac…² Piperacillin/tazobac…³
#>   <chr>     <chr>                  <chr>                  <chr>                 
#> 1 Gram-neg… 88% (85-91%,N=641)     99% (97-99%,N=691)     98% (97-99%,N=693)    
#> 2 Gram-pos… 86% (82-89%,N=345)     98% (96-98%,N=1044)    95% (93-97%,N=550)    
#> # ℹ abbreviated names: ¹​`Piperacillin/tazobactam`,
#> #   ²​`Piperacillin/tazobactam + Gentamicin`,
#> #   ³​`Piperacillin/tazobactam + Tobramycin`
#> # Use `ggplot2::autoplot()` or base R `plot()` to create a plot of this antibiogram,
#> # or use it directly in R Markdown or Quarto, see `antibiogram()`.

# you can use any antimicrobial selector with `+` too:
antibiogram(example_isolates,
  antimicrobials = ureidopenicillins() + c("", "GEN", "tobra"),
  mo_transform = "gramstain"
)
#> ℹ For `ureidopenicillins()` using column TZP (piperacillin/tazobactam)
#> # A tibble:       2 × 4
#> # Type:           Combination Antibiogram
#> # Conf. interval: 95%
#>   Pathogen  Piperacillin/tazobac…¹ Piperacillin/tazobac…² Piperacillin/tazobac…³
#>   <chr>     <chr>                  <chr>                  <chr>                 
#> 1 Gram-neg… 88% (85-91%,N=641)     99% (97-99%,N=691)     98% (97-99%,N=693)    
#> 2 Gram-pos… 86% (82-89%,N=345)     98% (96-98%,N=1044)    95% (93-97%,N=550)    
#> # ℹ abbreviated names: ¹​`Piperacillin/tazobactam`,
#> #   ²​`Piperacillin/tazobactam + Gentamicin`,
#> #   ³​`Piperacillin/tazobactam + Tobramycin`
#> # Use `ggplot2::autoplot()` or base R `plot()` to create a plot of this antibiogram,
#> # or use it directly in R Markdown or Quarto, see `antibiogram()`.

# names of antimicrobials do not need to resemble columns exactly:
antibiogram(example_isolates,
  antimicrobials = c("Cipro", "cipro + genta"),
  mo_transform = "gramstain",
  ab_transform = "name",
  sep = " & "
)
#> # A tibble:       2 × 3
#> # Type:           Traditional Antibiogram
#> # Conf. interval: 95%
#>   Pathogen      Ciprofloxacin      `Ciprofloxacin & Gentamicin`
#>   <chr>         <chr>              <chr>                       
#> 1 Gram-negative 91% (88-93%,N=684) 99% (97-99%,N=694)          
#> 2 Gram-positive 77% (74-80%,N=724) 93% (91-94%,N=847)          
#> # Use `ggplot2::autoplot()` or base R `plot()` to create a plot of this antibiogram,
#> # or use it directly in R Markdown or Quarto, see `antibiogram()`.


# Syndromic antibiogram (for AMR surveillance) --------------------------

antibiogram(example_isolates,
  antimicrobials = c(aminoglycosides(), carbapenems()),
  syndromic_group = "ward"
)
#> ℹ For `aminoglycosides()` using columns GEN (gentamicin), TOB (tobramycin), AMK
#>   (amikacin), and KAN (kanamycin)
#> ℹ For `carbapenems()` using columns IPM (imipenem) and MEM (meropenem)
#> # A tibble:       14 × 8
#> # Type:           Syndromic Antibiogram
#> # Conf. interval: 95%
#>    `Syndromic Group` Pathogen   Amikacin Gentamicin Imipenem Kanamycin Meropenem
#>    <chr>             <chr>      <chr>    <chr>      <chr>    <chr>     <chr>    
#>  1 Clinical          CoNS       NA       89% (84-9… 57% (39… NA        57% (39-…
#>  2 ICU               CoNS       NA       79% (68-8… NA       NA        NA       
#>  3 Outpatient        CoNS       NA       84% (66-9… NA       NA        NA       
#>  4 Clinical          E. coli    100% (9… 98% (96-9… 100% (9… NA        100% (99…
#>  5 ICU               E. coli    100% (9… 99% (95-1… 100% (9… NA        100% (97…
#>  6 Clinical          K. pneumo… NA       92% (81-9… 100% (9… NA        100% (92…
#>  7 Clinical          P. mirabi… NA       100% (88-… NA       NA        NA       
#>  8 Clinical          S. aureus  NA       99% (95-1… NA       NA        NA       
#>  9 ICU               S. aureus  NA       100% (95-… NA       NA        NA       
#> 10 Clinical          S. epider… NA       82% (72-9… NA       NA        NA       
#> 11 ICU               S. epider… NA       72% (60-8… NA       NA        NA       
#> 12 Clinical          S. hominis NA       96% (85-9… NA       NA        NA       
#> 13 Clinical          S. pneumo… 0% (0-5… 0% (0-5%,… NA       0% (0-5%… NA       
#> 14 ICU               S. pneumo… 0% (0-1… 0% (0-12%… NA       0% (0-12… NA       
#> # ℹ 1 more variable: Tobramycin <chr>
#> # Use `ggplot2::autoplot()` or base R `plot()` to create a plot of this antibiogram,
#> # or use it directly in R Markdown or Quarto, see `antibiogram()`.

# with a custom language, though this will be determined automatically
# (i.e., this table will be in Spanish on Spanish systems)
ex1 <- example_isolates[which(mo_genus() == "Escherichia"), ]
#> ℹ Using column mo as input for `mo_genus()`
antibiogram(ex1,
  antimicrobials = aminoglycosides(),
  ab_transform = "name",
  syndromic_group = ifelse(ex1$ward == "ICU",
    "UCI", "No UCI"
  ),
  language = "es"
)
#> ℹ For `aminoglycosides()` using columns GEN (gentamicin), TOB (tobramycin), AMK
#>   (amikacin), and KAN (kanamycin)
#> # A tibble:       2 × 5
#> # Type:           Syndromic Antibiogram
#> # Conf. interval: 95%
#>   `Grupo sindrómico` Patógeno Amikacina            Gentamicina       Tobramicina
#>   <chr>              <chr>    <chr>                <chr>             <chr>      
#> 1 No UCI             E. coli  100% (97-100%,N=119) 98% (96-99%,N=32… 98% (96-99…
#> 2 UCI                E. coli  100% (93-100%,N=52)  99% (95-100%,N=1… 96% (92-99…
#> # Use `ggplot2::autoplot()` or base R `plot()` to create a plot of this antibiogram,
#> # or use it directly in R Markdown or Quarto, see `antibiogram()`.


# Print the output for R Markdown / Quarto -----------------------------

ureido <- wisca(top_n_microorganisms(example_isolates, n = 10),
  antimicrobials = ureidopenicillins(),
  syndromic_group = "ward"
)
#> ℹ For `ureidopenicillins()` using column TZP (piperacillin/tazobactam)
#> Warning: In 1 of 3 WISCA regimen(s) or group(s), the prior weight of 1 per pathogen
#> accounts for more than 5% of the pathogen incidence (up to 5.6%), since there
#> are few isolates relative to the number of pathogens, i.e., rare pathogens or a
#> small group. Because WISCA gives every pathogen a prior weight of 1 (Bielicki
#> et al., 2016), this can strongly lower the estimated coverage without
#> reflecting the clinical expectation. Consider limiting the analysis to the most
#> common pathogens with `top_n_microorganisms()`, or using fewer and larger
#> groups.

# in an Rmd file, you would just need to return `ureido` in a chunk,
# but to be explicit here:
if (requireNamespace("knitr")) {
  cat(knitr::knit_print(ureido))
}
#> 
#> 
#> |Syndromic Group |Piperacillin/tazobactam |
#> |:---------------|:-----------------------|
#> |Clinical        |79.1% (72.2-85.8%)      |
#> |ICU             |57.3% (48.9-69%)        |
#> |Outpatient      |67.8% (49.1-87.9%)      |


# Generate plots with ggplot2 or base R --------------------------------

ab1 <- antibiogram(example_isolates,
  antimicrobials = c("AMC", "CIP", "TZP", "TZP+TOB"),
  mo_transform = "gramstain"
)
ab2 <- wisca(top_n_microorganisms(example_isolates, n = 10),
  antimicrobials = c("AMC", "CIP", "TZP", "TZP+TOB"),
  syndromic_group = "ward"
)
#> Warning: In 4 of 12 WISCA regimen(s) or group(s), the prior weight of 1 per pathogen
#> accounts for more than 5% of the pathogen incidence (up to 6.4%), since there
#> are few isolates relative to the number of pathogens, i.e., rare pathogens or a
#> small group. Because WISCA gives every pathogen a prior weight of 1 (Bielicki
#> et al., 2016), this can strongly lower the estimated coverage without
#> reflecting the clinical expectation. Consider limiting the analysis to the most
#> common pathogens with `top_n_microorganisms()`, or using fewer and larger
#> groups.

if (requireNamespace("ggplot2")) {
  ggplot2::autoplot(ab1)
}

if (requireNamespace("ggplot2")) {
  ggplot2::autoplot(ab2)
}


plot(ab1)

plot(ab2)


# }
```
