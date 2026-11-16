# Plotting for Classes `sir`, `mic` and `disk`

Functions to plot classes `sir`, `mic` and `disk`, with support for base
R and `ggplot2`.

## Usage

``` r
# S3 method for class 'mic'
plot(
  x,
  mo = NULL,
  ab = NULL,
  guideline = "EUCAST",
  main = deparse(substitute(x)),
  ylab = translate_AMR("Frequency", language = language),
  xlab = translate_AMR("Minimum Inhibitory Concentration (mg/L)", language = language),
  colours_SIR = c("#3CAEA3", "#F6D55C", "#ED553B"),
  language = get_AMR_locale(),
  expand = TRUE,
  include_PKPD = getOption("AMR_include_PKPD", TRUE),
  breakpoint_type = getOption("AMR_breakpoint_type", "human"),
  ...
)

# S3 method for class 'mic'
autoplot(
  object,
  mo = NULL,
  ab = NULL,
  guideline = "EUCAST",
  title = deparse(substitute(object)),
  ylab = translate_AMR("Frequency", language = language),
  xlab = translate_AMR("Minimum Inhibitory Concentration (mg/L)", language = language),
  colours_SIR = c("#3CAEA3", "#F6D55C", "#ED553B"),
  language = get_AMR_locale(),
  expand = TRUE,
  include_PKPD = getOption("AMR_include_PKPD", TRUE),
  breakpoint_type = getOption("AMR_breakpoint_type", "human"),
  ...
)

# S3 method for class 'mic'
fortify(object, ...)

# S3 method for class 'disk'
plot(
  x,
  main = deparse(substitute(x)),
  ylab = translate_AMR("Frequency", language = language),
  xlab = translate_AMR("Disk diffusion diameter (mm)", language = language),
  mo = NULL,
  ab = NULL,
  guideline = "EUCAST",
  colours_SIR = c("#3CAEA3", "#F6D55C", "#ED553B"),
  language = get_AMR_locale(),
  expand = TRUE,
  include_PKPD = getOption("AMR_include_PKPD", TRUE),
  breakpoint_type = getOption("AMR_breakpoint_type", "human"),
  ...
)

# S3 method for class 'disk'
autoplot(
  object,
  mo = NULL,
  ab = NULL,
  title = deparse(substitute(object)),
  ylab = translate_AMR("Frequency", language = language),
  xlab = translate_AMR("Disk diffusion diameter (mm)", language = language),
  guideline = "EUCAST",
  colours_SIR = c("#3CAEA3", "#F6D55C", "#ED553B"),
  language = get_AMR_locale(),
  expand = TRUE,
  include_PKPD = getOption("AMR_include_PKPD", TRUE),
  breakpoint_type = getOption("AMR_breakpoint_type", "human"),
  ...
)

# S3 method for class 'disk'
fortify(object, ...)

# S3 method for class 'sir'
plot(
  x,
  ylab = translate_AMR("Percentage", language = language),
  xlab = translate_AMR("Antimicrobial Interpretation", language = language),
  main = deparse(substitute(x)),
  language = get_AMR_locale(),
  ...
)

# S3 method for class 'sir'
autoplot(
  object,
  title = deparse(substitute(object)),
  xlab = translate_AMR("Antimicrobial Interpretation", language = language),
  ylab = translate_AMR("Frequency", language = language),
  colours_SIR = c("#3CAEA3", "#F6D55C", "#ED553B"),
  language = get_AMR_locale(),
  ...
)

# S3 method for class 'sir'
fortify(object, ...)
```

## Arguments

- x, object:

  values created with
  [`as.mic()`](https://amr-for-r.org/release/2.1.1/reference/as.mic.md),
  [`as.disk()`](https://amr-for-r.org/release/2.1.1/reference/as.disk.md)
  or
  [`as.sir()`](https://amr-for-r.org/release/2.1.1/reference/as.sir.md)
  (or their `random_*` variants, such as
  [`random_mic()`](https://amr-for-r.org/release/2.1.1/reference/random.md))

- mo:

  any (vector of) text that can be coerced to a valid microorganism code
  with
  [`as.mo()`](https://amr-for-r.org/release/2.1.1/reference/as.mo.md)

- ab:

  any (vector of) text that can be coerced to a valid antimicrobial drug
  code with
  [`as.ab()`](https://amr-for-r.org/release/2.1.1/reference/as.ab.md)

- guideline:

  interpretation guideline to use - the default is the latest included
  EUCAST guideline, see *Details*

- main, title:

  title of the plot

- xlab, ylab:

  axis title

- colours_SIR:

  colours to use for filling in the bars, must be a vector of three
  values (in the order S, I and R). The default colours are colour-blind
  friendly.

- language:

  language to be used to translate 'Susceptible', 'Increased
  exposure'/'Intermediate' and 'Resistant' - the default is system
  language (see
  [`get_AMR_locale()`](https://amr-for-r.org/release/2.1.1/reference/translate.md))
  and can be overwritten by setting the [package
  option](https://amr-for-r.org/release/2.1.1/reference/AMR-options.md)
  [`AMR_locale`](https://amr-for-r.org/release/2.1.1/reference/AMR-options.md),
  e.g. `options(AMR_locale = "de")`, see
  [translate](https://amr-for-r.org/release/2.1.1/reference/translate.md).
  Use `language = NULL` or `language = ""` to prevent translation.

- expand:

  a [logical](https://rdrr.io/r/base/logical.html) to indicate whether
  the range on the x axis should be expanded between the lowest and
  highest value. For MIC values, intermediate values will be factors of
  2 starting from the highest MIC value. For disk diameters, the whole
  diameter range will be filled.

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

- ...:

  arguments passed on to methods

## Value

The
[`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
functions return a
[`ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html) model
that is extendible with any `ggplot2` function.

The [`fortify()`](https://ggplot2.tidyverse.org/reference/fortify.html)
functions return a [data.frame](https://rdrr.io/r/base/data.frame.html)
as an extension for usage in the
[`ggplot2::ggplot()`](https://ggplot2.tidyverse.org/reference/ggplot.html)
function.

## Details

The interpretation of "I" will be named "Increased exposure" for all
EUCAST guidelines since 2019, and will be named "Intermediate" in all
other cases.

For interpreting MIC values as well as disk diffusion diameters,
supported guidelines to be used as input for the `guideline` argument
are: "EUCAST 2023", "EUCAST 2022", "EUCAST 2021", "EUCAST 2020", "EUCAST
2019", "EUCAST 2018", "EUCAST 2017", "EUCAST 2016", "EUCAST 2015",
"EUCAST 2014", "EUCAST 2013", "EUCAST 2012", "EUCAST 2011", "CLSI 2023",
"CLSI 2022", "CLSI 2021", "CLSI 2020", "CLSI 2019", "CLSI 2018", "CLSI
2017", "CLSI 2016", "CLSI 2015", "CLSI 2014", "CLSI 2013", "CLSI 2012",
and "CLSI 2011".

Simply using `"CLSI"` or `"EUCAST"` as input will automatically select
the latest version of that guideline.

## Examples

``` r
some_mic_values <- random_mic(size = 100)
some_disk_values <- random_disk(size = 100, mo = "Escherichia coli", ab = "cipro")
some_sir_values <- random_sir(50, prob_SIR = c(0.55, 0.05, 0.30))

plot(some_mic_values)
plot(some_disk_values)
plot(some_sir_values)

# when providing the microorganism and antibiotic, colours will show interpretations:
plot(some_mic_values, mo = "S. aureus", ab = "ampicillin")
plot(some_disk_values, mo = "Escherichia coli", ab = "cipro")
plot(some_disk_values, mo = "Escherichia coli", ab = "cipro", language = "nl")

# \donttest{
if (require("ggplot2")) {
  autoplot(some_mic_values)
}
if (require("ggplot2")) {
  autoplot(some_disk_values, mo = "Escherichia coli", ab = "cipro")
}
if (require("ggplot2")) {
  autoplot(some_sir_values)
}
# }
```
