# Get ATC Properties from WHOCC Website

Gets data from the WHOCC website to determine properties of an
Anatomical Therapeutic Chemical (ATC) (e.g. an antimicrobial), such as
the name, defined daily dose (DDD) or standard unit.

## Usage

``` r
atc_online_property(atc_code, property, administration = "O",
  url = "https://atcddd.fhi.no/atc_ddd_index/?code=%s&showdescription=no",
  url_vet = "https://atcddd.fhi.no/atcvet/atcvet_index/?code=%s&showdescription=no")

atc_online_groups(atc_code, ...)

atc_online_ddd(atc_code, ...)

atc_online_ddd_units(atc_code, ...)
```

## Source

<https://atcddd.fhi.no/atc_ddd_alterations__cumulative/ddd_alterations/abbrevations/>

## Arguments

- atc_code:

  A [character](https://rdrr.io/r/base/character.html) (vector) with ATC
  code(s) of antimicrobials, will be coerced with
  [`as.ab()`](https://amr-for-r.org/reference/as.ab.md) and
  [`ab_atc()`](https://amr-for-r.org/reference/ab_property.md)
  internally if not a valid ATC code.

- property:

  Property of an ATC code. Valid values are `"ATC"`, `"Name"`, `"DDD"`,
  `"U"` (`"unit"`), `"Adm.R"`, `"Note"` and `groups`. For this last
  option, all hierarchical groups of an ATC code will be returned, see
  *Examples*.

- administration:

  Type of administration when using `property = "Adm.R"`, see *Details*.

- url:

  URL of website of the WHOCC. The sign `%s` can be used as a
  placeholder for ATC codes.

- url_vet:

  URL of website of the WHOCC for veterinary medicine. The sign `%s` can
  be used as a placeholder for ATC_vet codes (that all start with "Q").

- ...:

  Arguments to pass on to `atc_property`.

## Details

Options for argument `administration`:

- `"Implant"` = Implant

- `"Inhal"` = Inhalation

- `"Instill"` = Instillation

- `"N"` = nasal

- `"O"` = oral

- `"P"` = parenteral

- `"R"` = rectal

- `"SL"` = sublingual/buccal

- `"TD"` = transdermal

- `"V"` = vaginal

Abbreviations of return values when using `property = "U"` (unit):

- `"g"` = gram

- `"mg"` = milligram

- `"mcg"` = microgram

- `"U"` = unit

- `"TU"` = thousand units

- `"MU"` = million units

- `"mmol"` = millimole

- `"ml"` = millilitre (e.g. eyedrops)

**N.B. This function requires an internet connection and only works if
the following packages are installed: `curl`, `rvest`, `xml2`.**

## Examples

``` r
# \donttest{
if (requireNamespace("curl") && requireNamespace("rvest") && requireNamespace("xml2")) {
  # oral DDD (Defined Daily Dose) of amoxicillin
  atc_online_property("J01CA04", "DDD", "O")
  atc_online_ddd(ab_atc("amox"))

  # parenteral DDD (Defined Daily Dose) of amoxicillin
  atc_online_property("J01CA04", "DDD", "P")

  atc_online_property("J01CA04", property = "groups") # search hierarchical groups of amoxicillin
}
#> Loading required namespace: rvest
#> ℹ in `atc_online_property()`: no properties found for ATC QG51AA03. Please
#>   check
#>   https://atcddd.fhi.no/atcvet/atcvet_index/?code=QG51AA03&showdescription=no.
#> ℹ in `atc_online_property()`: no properties found for ATC QJ01CA04. Please
#>   check
#>   https://atcddd.fhi.no/atcvet/atcvet_index/?code=QJ01CA04&showdescription=no.
#> [1] "ANTIINFECTIVES FOR SYSTEMIC USE"        
#> [2] "ANTIBACTERIALS FOR SYSTEMIC USE"        
#> [3] "BETA-LACTAM ANTIBACTERIALS, PENICILLINS"
#> [4] "Penicillins with extended spectrum"     
# }
```
