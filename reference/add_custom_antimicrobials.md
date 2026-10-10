# Add Custom Antimicrobials

With `add_custom_antimicrobials()` you can add your own custom
antimicrobial drug names and codes.

## Usage

``` r
add_custom_antimicrobials(x)

add_custom_antimicrobial_synonyms(ab, synonyms = NULL)

clear_custom_antimicrobials()
```

## Arguments

- x:

  A [data.frame](https://rdrr.io/r/base/data.frame.html) resembling the
  [antimicrobials](https://amr-for-r.org/reference/antimicrobials.md)
  data set, at least containing columns "ab" and "name". To add synonyms
  to existing antimicrobials, columns "ab" and "synonyms" suffice (see
  *Details*).

- ab:

  For `add_custom_antimicrobial_synonyms()`: a
  [character](https://rdrr.io/r/base/character.html) vector of existing
  antimicrobial codes (see
  [`as.ab()`](https://amr-for-r.org/reference/as.ab.md)), or a
  [data.frame](https://rdrr.io/r/base/data.frame.html) with columns "ab"
  and "synonyms" (one synonym per row).

- synonyms:

  A [character](https://rdrr.io/r/base/character.html) vector of
  synonyms. If `ab` has length 1, all `synonyms` are added to that
  antimicrobial; otherwise `ab` and `synonyms` must have the same
  length.

## Details

**Important:** Due to how R works, the `add_custom_antimicrobials()`
function has to be run in every R session - added antimicrobials are not
stored between sessions and are thus lost when R is exited.

There are two ways to circumvent this and automate the process of adding
antimicrobials:

**Method 1:** Using the package option
[`AMR_custom_ab`](https://amr-for-r.org/reference/AMR-options.md), which
is the preferred method. To use this method:

1.  Create a data set in the structure of the
    [antimicrobials](https://amr-for-r.org/reference/antimicrobials.md)
    data set (containing at the very least columns "ab" and "name", or
    only columns "ab" and "synonyms" for synonyms of existing
    antimicrobials) and save it with
    [`saveRDS()`](https://rdrr.io/r/base/readRDS.html) to a location of
    choice, e.g. `"~/my_custom_ab.rds"`, or any remote location.

2.  Set the file location to the package option
    [`AMR_custom_ab`](https://amr-for-r.org/reference/AMR-options.md):
    `options(AMR_custom_ab = "~/my_custom_ab.rds")`. This can even be a
    remote file location, such as an https URL. Since options are not
    saved between R sessions, it is best to save this option to the
    `.Rprofile` file so that it will be loaded on start-up of R. To do
    this, open the `.Rprofile` file using e.g.
    `utils::file.edit("~/.Rprofile")`, add this text and save the file:

        # Add custom antimicrobial codes:
        options(AMR_custom_ab = "~/my_custom_ab.rds")

    Upon package load, this file will be loaded and run through the
    `add_custom_antimicrobials()` function.

    This file can also hold synonyms for antimicrobials that already
    exist, such as local trade names: rows with an existing code in
    column "ab" and only column "synonyms" filled (one synonym per row,
    codes may repeat) are passed on to
    `add_custom_antimicrobial_synonyms()`. Rows with an existing code
    that have any other column filled still give an error. Synonyms must
    be valid UTF-8 text; convert other encodings first, e.g. with
    `iconv(x, from = "CP949", to = "UTF-8")`.

**Method 2:** Loading the antimicrobial additions directly from your
`.Rprofile` file. Note that the definitions will be stored in a
user-specific R file, which is a suboptimal workflow. To use this
method:

1.  Edit the `.Rprofile` file using e.g.
    `utils::file.edit("~/.Rprofile")`.

2.  Add a text like below and save the file:

         # Add custom antibiotic drug codes:
         AMR::add_custom_antimicrobials(
           data.frame(ab = "TESTAB",
                      name = "Test Antibiotic",
                      group = "Test Group")
         )

Use `add_custom_antimicrobial_synonyms()` to add extra names, such as
local trade names, to antimicrobials that already exist (including
antimicrobials added with `add_custom_antimicrobials()`). These synonyms
are recognised by [`as.ab()`](https://amr-for-r.org/reference/as.ab.md)
and all `ab_*()` functions. Synonyms may contain non-ASCII characters
(e.g. trade names in other scripts): they are matched against the input
as given, before [`as.ab()`](https://amr-for-r.org/reference/as.ab.md)
transliterates the input to ASCII. Matching ignores case, all white
space (including no-break and ideographic spaces) and zero-width
characters. If the input does not match as given, a trailing part in
parentheses or brackets (such as `(piperacillin, tazobactam)`) and a
trailing strength written as a number with a unit symbol (such as
`4.5g`, `4/0.5 g` or `500 mg/mL`) are removed one at a time, and
matching is tried again after each removal, so that the longest
registered synonym wins. Full-width characters are read as their ASCII
equivalents. No words are removed, in any language: to recognise a name
with a dosage form (such as `Inj`), add that variant as a synonym as
well. Only user-added synonyms are matched this way, there is no other
fuzzy matching for them, and a strength on its own is never matched to a
user-added synonym. A synonym cannot be added if it already identifies
another antimicrobial as its code, name, synonym, abbreviation, ATC
code, CID or LOINC code, also when written with other punctuation (such
as `J01-CR05`).

Use `clear_custom_antimicrobials()` to clear the previously added
antimicrobials and synonyms.

## See also

[`add_custom_microorganisms()`](https://amr-for-r.org/reference/add_custom_microorganisms.md)
to add custom microorganisms.

## Examples

``` r
# \donttest{
# returns a wildly guessed result:
as.ab("testab")
#> Class <ab>
#> [1] THA

# now add a custom entry - it will be considered by as.ab() and
# all ab_*() functions
add_custom_antimicrobials(
  data.frame(
    ab = "TESTAB",
    name = "Test Antibiotic",
    # you can add any property present in the
    # 'antimicrobials' data set, such as 'group':
    group = "Test Group"
  )
)
#> ℹ Added one record to the internal `antimicrobials` data set.

# "testab" is now a new antibiotic:
as.ab("testab")
#> Class <ab>
#> [1] TESTAB
ab_name("testab")
#> [1] "Test Antibiotic"
ab_group("testab")
#> [1] "Test Group"

ab_info("testab")
#> $ab
#> [1] "TESTAB"
#> 
#> $cid
#> [1] NA
#> 
#> $name
#> [1] "Test Antibiotic"
#> 
#> $group
#> [1] "Test Group"
#> 
#> $atc
#> [1] NA
#> 
#> $atc_group1
#> [1] NA
#> 
#> $atc_group2
#> [1] NA
#> 
#> $tradenames
#> [1] NA
#> 
#> $loinc
#> [1] NA
#> 
#> $ddd
#> $ddd$oral
#> $ddd$oral$amount
#> [1] NA
#> 
#> $ddd$oral$units
#> [1] NA
#> 
#> 
#> $ddd$iv
#> $ddd$iv$amount
#> [1] NA
#> 
#> $ddd$iv$units
#> [1] NA
#> 
#> 
#> 


# Add Co-fluampicil, which is one of the many J01CR50 codes, see
# https://atcddd.fhi.no/ddd/list_of_ddds_combined_products/
add_custom_antimicrobials(
  data.frame(
    ab = "COFLU",
    name = "Co-fluampicil",
    atc = "J01CR50",
    group = "Beta-lactams/penicillins"
  )
)
#> ℹ Added one record to the internal `antimicrobials` data set.
ab_atc("Co-fluampicil")
#> [1] "J01CR50"
ab_name("J01CR50")
#> [1] "Co-fluampicil"

# even antimicrobial selectors work
# see ?amr_selector
x <- data.frame(
  random_column = "some value",
  coflu = as.sir("S"),
  ampicillin = as.sir("R")
)
x
#>   random_column coflu ampicillin
#> 1    some value     S          R
x[, betalactams()]
#> ℹ For `betalactams()` using columns coflu (co-fluampicil) and ampicillin
#>   coflu ampicillin
#> 1     S          R

# add local trade names as synonyms to existing antimicrobials,
# here a Korean trade name of piperacillin/tazobactam
add_custom_antimicrobial_synonyms("TZP", "\uD0C0\uC870\uC2E0\uC8FC")
#> ℹ Added one synonym to one antimicrobial.
ab_name("\uD0C0\uC870\uC2E0\uC8FC")
#> [1] "Piperacillin/tazobactam"
# a trailing strength, as in many hospital exports, is ignored
as.ab("\uD0C0\uC870\uC2E0\uC8FC 4.5g")
#> Class <ab>
#> [1] TZP
# }
```
