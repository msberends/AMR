# Generate Antibiogram: Traditional, Combined, Syndromic, or Weighted-Incidence Syndromic Combination (WISCA)

Generate an antibiogram, and communicate the results in plots or tables.
These functions follow the logic of Klinker *et al.* and Barbieri *et
al.* (see *Source*), and allow reporting in e.g. R Markdown and Quarto
as well.

## Usage

``` r
antibiogram(
  x,
  antibiotics = where(is.sir),
  mo_transform = "shortname",
  ab_transform = NULL,
  syndromic_group = NULL,
  add_total_n = TRUE,
  only_all_tested = FALSE,
  digits = 0,
  col_mo = NULL,
  language = get_AMR_locale(),
  minimum = 30,
  combine_SI = TRUE,
  sep = " + ",
  info = interactive()
)

# S3 method for class 'antibiogram'
plot(x, ...)

# S3 method for class 'antibiogram'
autoplot(object, ...)

# S3 method for class 'antibiogram'
knit_print(
  x,
  italicise = TRUE,
  na = getOption("knitr.kable.NA", default = ""),
  ...
)
```

## Source

- Klinker KP *et al.* (2021). **Antimicrobial stewardship and
  antibiograms: importance of moving beyond traditional antibiograms**.
  *Therapeutic Advances in Infectious Disease*, May
  5;8:20499361211011373;
  [doi:10.1177/20499361211011373](https://doi.org/10.1177/20499361211011373)

- Barbieri E *et al.* (2021). **Development of a Weighted-Incidence
  Syndromic Combination Antibiogram (WISCA) to guide the choice of the
  empiric antibiotic treatment for urinary tract infection in paediatric
  patients: a Bayesian approach** *Antimicrobial Resistance & Infection
  Control* May 1;10(1):74;
  [doi:10.1186/s13756-021-00939-2](https://doi.org/10.1186/s13756-021-00939-2)

- **M39 Analysis and Presentation of Cumulative Antimicrobial
  Susceptibility Test Data, 5th Edition**, 2022, *Clinical and
  Laboratory Standards Institute (CLSI)*.
  <https://clsi.org/standards/products/microbiology/documents/m39/>.

## Arguments

- x:

  a [data.frame](https://rdrr.io/r/base/data.frame.html) containing at
  least a column with microorganisms and columns with antibiotic results
  (class 'sir', see
  [`as.sir()`](https://amr-for-r.org/release/2.0.0/reference/as.sir.md))

- antibiotics:

  vector of any antibiotic name or code (will be evaluated with
  [`as.ab()`](https://amr-for-r.org/release/2.0.0/reference/as.ab.md),
  column name of `x`, or (any combinations of) [antibiotic
  selectors](https://amr-for-r.org/release/2.0.0/reference/antibiotic_class_selectors.md)
  such as
  [`aminoglycosides()`](https://amr-for-r.org/release/2.0.0/reference/antibiotic_class_selectors.md)
  or
  [`carbapenems()`](https://amr-for-r.org/release/2.0.0/reference/antibiotic_class_selectors.md).
  For combination antibiograms, this can also be set to values separated
  with `"+"`, such as "TZP+TOB" or "cipro + genta", given that columns
  resembling such antibiotics exist in `x`. See *Examples*.

- mo_transform:

  a character to transform microorganism input - must be "name",
  "shortname", "gramstain", or one of the column names of the
  [microorganisms](https://amr-for-r.org/release/2.0.0/reference/microorganisms.md)
  data set: "mo", "fullname", "status", "kingdom", "phylum", "class",
  "order", "family", "genus", "species", "subspecies", "rank", "ref",
  "source", "lpsn", "lpsn_parent", "lpsn_renamed_to", "gbif",
  "gbif_parent", "gbif_renamed_to", "prevalence", or "snomed". Can also
  be `NULL` to not transform the input.

- ab_transform:

  a character to transform antibiotic input - must be one of the column
  names of the
  [antibiotics](https://amr-for-r.org/release/2.0.0/reference/antibiotics.md)
  data set: "ab", "cid", "name", "group", "atc", "atc_group1",
  "atc_group2", "abbreviations", "synonyms", "oral_ddd", "oral_units",
  "iv_ddd", "iv_units", or "loinc". Can also be `NULL` to not transform
  the input.

- syndromic_group:

  a column name of `x`, or values calculated to split rows of `x`, e.g.
  by using [`ifelse()`](https://rdrr.io/r/base/ifelse.html) or
  [`case_when()`](https://dplyr.tidyverse.org/reference/case-and-replace-when.html).
  See *Examples*.

- add_total_n:

  a [logical](https://rdrr.io/r/base/logical.html) to indicate whether
  total available numbers per pathogen should be added to the table
  (default is `TRUE`). This will add the lowest and highest number of
  available isolate per antibiotic (e.g, if for *E. coli* 200 isolates
  are available for ciprofloxacin and 150 for amoxicillin, the returned
  number will be "150-200").

- only_all_tested:

  (for combination antibiograms): a
  [logical](https://rdrr.io/r/base/logical.html) to indicate that
  isolates must be tested for all antibiotics, see *Details*

- digits:

  number of digits to use for rounding

- col_mo:

  column name of the names or codes of the microorganisms (see
  [`as.mo()`](https://amr-for-r.org/release/2.0.0/reference/as.mo.md)) -
  the default is the first column of class
  [`mo`](https://amr-for-r.org/release/2.0.0/reference/as.mo.md). Values
  will be coerced using
  [`as.mo()`](https://amr-for-r.org/release/2.0.0/reference/as.mo.md).

- language:

  language to translate text, which defaults to the system language (see
  [`get_AMR_locale()`](https://amr-for-r.org/release/2.0.0/reference/translate.md))

- minimum:

  the minimum allowed number of available (tested) isolates. Any isolate
  count lower than `minimum` will return `NA` with a warning. The
  default number of `30` isolates is advised by the Clinical and
  Laboratory Standards Institute (CLSI) as best practice, see *Source*.

- combine_SI:

  a [logical](https://rdrr.io/r/base/logical.html) to indicate whether
  all susceptibility should be determined by results of either S or I,
  instead of only S (default is `TRUE`)

- sep:

  a separating character for antibiotic columns in combination
  antibiograms

- info:

  a [logical](https://rdrr.io/r/base/logical.html) to indicate info
  should be printed - the default is `TRUE` only in interactive mode

- ...:

  when used in [R Markdown or
  Quarto](https://rdrr.io/pkg/knitr/man/kable.html): arguments passed on
  to [`knitr::kable()`](https://rdrr.io/pkg/knitr/man/kable.html)
  (otherwise, has no use)

- object:

  an `antibiogram()` object

- italicise:

  a [logical](https://rdrr.io/r/base/logical.html) to indicate whether
  the microorganism names in the
  [knitr](https://rdrr.io/pkg/knitr/man/kable.html) table should be made
  italic, using
  [`italicise_taxonomy()`](https://amr-for-r.org/release/2.0.0/reference/italicise_taxonomy.md).

- na:

  character to use for showing `NA` values

## Details

This function returns a table with values between 0 and 100 for
*susceptibility*, not resistance.

**Remember that you should filter your data to let it contain only first
isolates!** This is needed to exclude duplicates and to reduce selection
bias. Use
[`first_isolate()`](https://amr-for-r.org/release/2.0.0/reference/first_isolate.md)
to determine them in your data set with one of the four available
algorithms.

All types of antibiograms as listed below can be plotted (using
[`ggplot2::autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
or base R
[`plot()`](https://amr-for-r.org/release/2.0.0/reference/plot.md)/[`barplot()`](https://rdrr.io/r/graphics/barplot.html)).
The `antibiogram` object can also be used directly in R Markdown /
Quarto (i.e., `knitr`) for reports. In this case,
[`knitr::kable()`](https://rdrr.io/pkg/knitr/man/kable.html) will be
applied automatically and microorganism names will even be printed in
italics at default (see argument `italicise`). You can also use
functions from specific 'table reporting' packages to transform the
output of `antibiogram()` to your needs, e.g. with
`flextable::as_flextable()` or `gt::gt()`.

### Antibiogram Types

There are four antibiogram types, as proposed by Klinker *et al.* (2021,
[doi:10.1177/20499361211011373](https://doi.org/10.1177/20499361211011373)
), and they are all supported by `antibiogram()`:

1.  **Traditional Antibiogram**

    Case example: Susceptibility of *Pseudomonas aeruginosa* to
    piperacillin/tazobactam (TZP)

    Code example:

        antibiogram(your_data,
                    antibiotics = "TZP")

2.  **Combination Antibiogram**

    Case example: Additional susceptibility of *Pseudomonas aeruginosa*
    to TZP + tobramycin versus TZP alone

    Code example:

        antibiogram(your_data,
                    antibiotics = c("TZP", "TZP+TOB", "TZP+GEN"))

3.  **Syndromic Antibiogram**

    Case example: Susceptibility of *Pseudomonas aeruginosa* to TZP
    among respiratory specimens (obtained among ICU patients only)

    Code example:

        antibiogram(your_data,
                    antibiotics = penicillins(),
                    syndromic_group = "ward")

4.  **Weighted-Incidence Syndromic Combination Antibiogram (WISCA)**

    Case example: Susceptibility of *Pseudomonas aeruginosa* to TZP
    among respiratory specimens (obtained among ICU patients only) for
    male patients age \>=65 years with heart failure

    Code example:

        library(dplyr)
        your_data %>%
          filter(ward == "ICU" & specimen_type == "Respiratory") %>%
          antibiogram(antibiotics = c("TZP", "TZP+TOB", "TZP+GEN"),
                      syndromic_group = ifelse(.$age >= 65 &
                                                 .$gender == "Male" &
                                                 .$condition == "Heart Disease",
                                               "Study Group", "Control Group"))

Note that for combination antibiograms, it is important to realise that
susceptibility can be calculated in two ways, which can be set with the
`only_all_tested` argument (default is `FALSE`). See this example for
two antibiotics, Drug A and Drug B, about how `antibiogram()` works to
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

## Examples

``` r
# example_isolates is a data set available in the AMR package.
# run ?example_isolates for more info.
example_isolates

# \donttest{
# Traditional antibiogram ----------------------------------------------

antibiogram(example_isolates,
  antibiotics = c(aminoglycosides(), carbapenems())
)

antibiogram(example_isolates,
  antibiotics = aminoglycosides(),
  ab_transform = "atc",
  mo_transform = "gramstain"
)

antibiogram(example_isolates,
  antibiotics = carbapenems(),
  ab_transform = "name",
  mo_transform = "name"
)


# Combined antibiogram -------------------------------------------------

# combined antibiotics yield higher empiric coverage
antibiogram(example_isolates,
  antibiotics = c("TZP", "TZP+TOB", "TZP+GEN"),
  mo_transform = "gramstain"
)

# names of antibiotics do not need to resemble columns exactly:
antibiogram(example_isolates,
  antibiotics = c("Cipro", "cipro + genta"),
  mo_transform = "gramstain",
  ab_transform = "name",
  sep = " & "
)


# Syndromic antibiogram ------------------------------------------------

# the data set could contain a filter for e.g. respiratory specimens
antibiogram(example_isolates,
  antibiotics = c(aminoglycosides(), carbapenems()),
  syndromic_group = "ward"
)

# now define a data set with only E. coli
ex1 <- example_isolates[which(mo_genus() == "Escherichia"), ]

# with a custom language, though this will be determined automatically
# (i.e., this table will be in Spanish on Spanish systems)
antibiogram(ex1,
  antibiotics = aminoglycosides(),
  ab_transform = "name",
  syndromic_group = ifelse(ex1$ward == "ICU",
    "UCI", "No UCI"
  ),
  language = "es"
)


# Weighted-incidence syndromic combination antibiogram (WISCA) ---------

# the data set could contain a filter for e.g. respiratory specimens/ICU
antibiogram(example_isolates,
  antibiotics = c("AMC", "AMC+CIP", "TZP", "TZP+TOB"),
  mo_transform = "gramstain",
  minimum = 10, # this should be >=30, but now just as example
  syndromic_group = ifelse(example_isolates$age >= 65 &
    example_isolates$gender == "M",
  "WISCA Group 1", "WISCA Group 2"
  )
)


# Print the output for R Markdown / Quarto -----------------------------

ureido <- antibiogram(example_isolates,
  antibiotics = ureidopenicillins(),
  ab_transform = "name"
)

# in an Rmd file, you would just need to return `ureido` in a chunk,
# but to be explicit here:
if (requireNamespace("knitr")) {
  knitr::knit_print(ureido)
}


# Generate plots with ggplot2 or base R --------------------------------

ab1 <- antibiogram(example_isolates,
  antibiotics = c("AMC", "CIP", "TZP", "TZP+TOB"),
  mo_transform = "gramstain"
)
ab2 <- antibiogram(example_isolates,
  antibiotics = c("AMC", "CIP", "TZP", "TZP+TOB"),
  mo_transform = "gramstain",
  syndromic_group = "ward"
)

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
