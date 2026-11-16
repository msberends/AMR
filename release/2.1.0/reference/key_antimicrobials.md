# (Key) Antimicrobials for First Weighted Isolates

These functions can be used to determine first weighted isolates by
considering the phenotype for isolate selection (see
[`first_isolate()`](https://amr-for-r.org/release/2.1.0/reference/first_isolate.md)).
Using a phenotype-based method to determine first isolates is more
reliable than methods that disregard phenotypes.

## Usage

``` r
key_antimicrobials(
  x = NULL,
  col_mo = NULL,
  universal = c("ampicillin", "amoxicillin/clavulanic acid", "cefuroxime",
    "piperacillin/tazobactam", "ciprofloxacin", "trimethoprim/sulfamethoxazole"),
  gram_negative = c("gentamicin", "tobramycin", "colistin", "cefotaxime", "ceftazidime",
    "meropenem"),
  gram_positive = c("vancomycin", "teicoplanin", "tetracycline", "erythromycin",
    "oxacillin", "rifampin"),
  antifungal = c("anidulafungin", "caspofungin", "fluconazole", "miconazole", "nystatin",
    "voriconazole"),
  only_sir_columns = FALSE,
  ...
)

all_antimicrobials(x = NULL, only_sir_columns = FALSE, ...)

antimicrobials_equal(
  y,
  z,
  type = c("points", "keyantimicrobials"),
  ignore_I = TRUE,
  points_threshold = 2,
  ...
)
```

## Arguments

- x:

  a [data.frame](https://rdrr.io/r/base/data.frame.html) with
  antibiotics columns, like `AMX` or `amox`. Can be left blank to
  determine automatically

- col_mo:

  column name of the names or codes of the microorganisms (see
  [`as.mo()`](https://amr-for-r.org/release/2.1.0/reference/as.mo.md)) -
  the default is the first column of class
  [`mo`](https://amr-for-r.org/release/2.1.0/reference/as.mo.md). Values
  will be coerced using
  [`as.mo()`](https://amr-for-r.org/release/2.1.0/reference/as.mo.md).

- universal:

  names of **broad-spectrum** antimicrobial drugs, case-insensitive. Set
  to `NULL` to ignore. See *Details* for the default antimicrobial drugs

- gram_negative:

  names of antibiotic drugs for **Gram-positives**, case-insensitive.
  Set to `NULL` to ignore. See *Details* for the default antibiotic
  drugs

- gram_positive:

  names of antibiotic drugs for **Gram-negatives**, case-insensitive.
  Set to `NULL` to ignore. See *Details* for the default antibiotic
  drugs

- antifungal:

  names of antifungal drugs for **fungi**, case-insensitive. Set to
  `NULL` to ignore. See *Details* for the default antifungal drugs

- only_sir_columns:

  a [logical](https://rdrr.io/r/base/logical.html) to indicate whether
  only columns must be included that were transformed to class `sir`
  (see
  [`as.sir()`](https://amr-for-r.org/release/2.1.0/reference/as.sir.md))
  on beforehand (default is `FALSE`)

- ...:

  ignored, only in place to allow future extensions

- y, z:

  [character](https://rdrr.io/r/base/character.html) vectors to compare

- type:

  type to determine weighed isolates; can be `"keyantimicrobials"` or
  `"points"`, see *Details*

- ignore_I:

  [logical](https://rdrr.io/r/base/logical.html) to indicate whether
  antibiotic interpretations with `"I"` will be ignored when
  `type = "keyantimicrobials"`, see *Details*

- points_threshold:

  minimum number of points to require before differences in the
  antibiogram will lead to inclusion of an isolate when
  `type = "points"`, see *Details*

## Details

The `key_antimicrobials()` and `all_antimicrobials()` functions are
context-aware. This means that the `x` argument can be left blank if
used inside a [data.frame](https://rdrr.io/r/base/data.frame.html) call,
see *Examples*.

The function `key_antimicrobials()` returns a
[character](https://rdrr.io/r/base/character.html) vector with 12
antimicrobial results for every isolate. The function
`all_antimicrobials()` returns a
[character](https://rdrr.io/r/base/character.html) vector with all
antimicrobial drug results for every isolate. These vectors can then be
compared using `antimicrobials_equal()`, to check if two isolates have
generally the same antibiogram. Missing and invalid values are replaced
with a dot (`"."`) by `key_antimicrobials()` and ignored by
`antimicrobials_equal()`.

Please see the
[`first_isolate()`](https://amr-for-r.org/release/2.1.0/reference/first_isolate.md)
function how these important functions enable the 'phenotype-based'
method for determination of first isolates.

The default antimicrobial drugs used for **all rows** (set in
`universal`) are:

- Ampicillin

- Amoxicillin/clavulanic acid

- Cefuroxime

- Ciprofloxacin

- Piperacillin/tazobactam

- Trimethoprim/sulfamethoxazole

The default antimicrobial drugs used for **Gram-negative bacteria** (set
in `gram_negative`) are:

- Cefotaxime

- Ceftazidime

- Colistin

- Gentamicin

- Meropenem

- Tobramycin

The default antimicrobial drugs used for **Gram-positive bacteria** (set
in `gram_positive`) are:

- Erythromycin

- Oxacillin

- Rifampin

- Teicoplanin

- Tetracycline

- Vancomycin

The default antimicrobial drugs used for **fungi** (set in `antifungal`)
are:

- Anidulafungin

- Caspofungin

- Fluconazole

- Miconazole

- Nystatin

- Voriconazole

## See also

[`first_isolate()`](https://amr-for-r.org/release/2.1.0/reference/first_isolate.md)

## Examples

``` r
# `example_isolates` is a data set available in the AMR package.
# See ?example_isolates.

# output of the `key_antimicrobials()` function could be like this:
strainA <- "SSSRR.S.R..S"
strainB <- "SSSIRSSSRSSS"

# those strings can be compared with:
antimicrobials_equal(strainA, strainB, type = "keyantimicrobials")
# TRUE, because I is ignored (as well as missing values)

antimicrobials_equal(strainA, strainB, type = "keyantimicrobials", ignore_I = FALSE)
# FALSE, because I is not ignored and so the 4th [character] differs

# \donttest{
if (require("dplyr")) {
  # set key antibiotics to a new variable
  my_patients <- example_isolates %>%
    mutate(keyab = key_antimicrobials(antifungal = NULL)) %>% # no need to define `x`
    mutate(
      # now calculate first isolates
      first_regular = first_isolate(col_keyantimicrobials = FALSE),
      # and first WEIGHTED isolates
      first_weighted = first_isolate(col_keyantimicrobials = "keyab")
    )

  # Check the difference in this data set, 'weighted' results in more isolates:
  sum(my_patients$first_regular, na.rm = TRUE)
  sum(my_patients$first_weighted, na.rm = TRUE)
}
# }
```
