# Random MIC Values/Disk Zones/SIR Generation

These functions can be used for generating random MIC values and disk
diffusion diameters, for AMR data analysis practice. By providing a
microorganism and antimicrobial drug, the generated results will reflect
reality as much as possible.

## Usage

``` r
random_mic(size = NULL, mo = NULL, ab = NULL, ...)

random_disk(size = NULL, mo = NULL, ab = NULL, ...)

random_sir(size = NULL, prob_SIR = c(0.33, 0.33, 0.33), ...)
```

## Arguments

- size:

  desired size of the returned vector. If used in a
  [data.frame](https://rdrr.io/r/base/data.frame.html) call or `dplyr`
  verb, will get the current (group) size if left blank.

- mo:

  any [character](https://rdrr.io/r/base/character.html) that can be
  coerced to a valid microorganism code with
  [`as.mo()`](https://amr-for-r.org/release/2.1.1/reference/as.mo.md)

- ab:

  any [character](https://rdrr.io/r/base/character.html) that can be
  coerced to a valid antimicrobial drug code with
  [`as.ab()`](https://amr-for-r.org/release/2.1.1/reference/as.ab.md)

- ...:

  ignored, only in place to allow future extensions

- prob_SIR:

  a vector of length 3: the probabilities for "S" (1st value), "I" (2nd
  value) and "R" (3rd value)

## Value

class `mic` for `random_mic()` (see
[`as.mic()`](https://amr-for-r.org/release/2.1.1/reference/as.mic.md))
and class `disk` for `random_disk()` (see
[`as.disk()`](https://amr-for-r.org/release/2.1.1/reference/as.disk.md))

## Details

The base R function [`sample()`](https://rdrr.io/r/base/sample.html) is
used for generating values.

Generated values are based on the EUCAST 2023 guideline as implemented
in the
[clinical_breakpoints](https://amr-for-r.org/release/2.1.1/reference/clinical_breakpoints.md)
data set. To create specific generated values per bug or drug, set the
`mo` and/or `ab` argument.

## Examples

``` r
random_mic(25)
random_disk(25)
random_sir(25)

# \donttest{
# make the random generation more realistic by setting a bug and/or drug:
random_mic(25, "Klebsiella pneumoniae") # range 0.0625-64
random_mic(25, "Klebsiella pneumoniae", "meropenem") # range 0.0625-16
random_mic(25, "Streptococcus pneumoniae", "meropenem") # range 0.0625-4

random_disk(25, "Klebsiella pneumoniae") # range 8-50
random_disk(25, "Klebsiella pneumoniae", "ampicillin") # range 11-17
random_disk(25, "Streptococcus pneumoniae", "ampicillin") # range 12-27
# }
```
