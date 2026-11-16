# Calculate the Mean AMR Distance

Calculates a normalised mean for antimicrobial resistance between
multiple observations, to help to identify similar isolates without
comparing antibiograms by hand.

## Usage

``` r
mean_amr_distance(x, ...)

# S3 method for class 'sir'
mean_amr_distance(x, ..., combine_SI = TRUE)

# S3 method for class 'data.frame'
mean_amr_distance(x, ..., combine_SI = TRUE)

amr_distance_from_row(amr_distance, row)
```

## Arguments

- x:

  a vector of class
  [sir](https://amr-for-r.org/release/2.0.0/reference/as.sir.md),
  [mic](https://amr-for-r.org/release/2.0.0/reference/as.mic.md) or
  [disk](https://amr-for-r.org/release/2.0.0/reference/as.disk.md), or a
  [data.frame](https://rdrr.io/r/base/data.frame.html) containing
  columns of any of these classes

- ...:

  variables to select (supports [tidyselect
  language](https://tidyselect.r-lib.org/reference/language.html) such
  as `column1:column4` and `where(is.mic)`, and can thus also be
  [antibiotic
  selectors](https://amr-for-r.org/release/2.0.0/reference/antibiotic_class_selectors.md)

- combine_SI:

  a [logical](https://rdrr.io/r/base/logical.html) to indicate whether
  all values of S and I must be merged into one, so the input only
  consists of S+I vs. R (susceptible vs. resistant) - the default is
  `TRUE`

- amr_distance:

  the outcome of `mean_amr_distance()`

- row:

  an index, such as a row number

## Details

The mean AMR distance is effectively [the
Z-score](https://en.wikipedia.org/wiki/Standard_score); a normalised
numeric value to compare AMR test results which can help to identify
similar isolates, without comparing antibiograms by hand.

MIC values (see
[`as.mic()`](https://amr-for-r.org/release/2.0.0/reference/as.mic.md))
are transformed with [`log2()`](https://rdrr.io/r/base/Log.html) first;
their distance is thus calculated as
`(log2(x) - mean(log2(x))) / sd(log2(x))`.

SIR values (see
[`as.sir()`](https://amr-for-r.org/release/2.0.0/reference/as.sir.md))
are transformed using `"S"` = 1, `"I"` = 2, and `"R"` = 3. If
`combine_SI` is `TRUE` (default), the `"I"` will be considered to be 1.

For data sets, the mean AMR distance will be calculated per column,
after which the mean per row will be returned, see *Examples*.

Use `amr_distance_from_row()` to subtract distances from the distance of
one row, see *Examples*.

## Interpretation

Isolates with distances less than 0.01 difference from each other should
be considered similar. Differences lower than 0.025 should be considered
suspicious.

## Examples

``` r
sir <- random_sir(10)
sir
mean_amr_distance(sir)

mic <- random_mic(10)
mic
mean_amr_distance(mic)
# equal to the Z-score of their log2:
(log2(mic) - mean(log2(mic))) / sd(log2(mic))

disk <- random_disk(10)
disk
mean_amr_distance(disk)

y <- data.frame(
  id = LETTERS[1:10],
  amox = random_sir(10, ab = "amox", mo = "Escherichia coli"),
  cipr = random_disk(10, ab = "cipr", mo = "Escherichia coli"),
  gent = random_mic(10, ab = "gent", mo = "Escherichia coli"),
  tobr = random_mic(10, ab = "tobr", mo = "Escherichia coli")
)
y
mean_amr_distance(y)
y$amr_distance <- mean_amr_distance(y, where(is.mic))
y[order(y$amr_distance), ]

if (require("dplyr")) {
  y %>%
    mutate(
      amr_distance = mean_amr_distance(y),
      check_id_C = amr_distance_from_row(amr_distance, id == "C")
    ) %>%
    arrange(check_id_C)
}
if (require("dplyr")) {
  # support for groups
  example_isolates %>%
    filter(mo_genus() == "Enterococcus" & mo_species() != "") %>%
    select(mo, TCY, carbapenems()) %>%
    group_by(mo) %>%
    mutate(dist = mean_amr_distance(.)) %>%
    arrange(mo, dist)
}
```
