# Get Properties of a Microorganism

Use these functions to return a specific property of a microorganism
based on the latest accepted taxonomy. All input values will be
evaluated internally with
[`as.mo()`](https://amr-for-r.org/reference/as.mo.md), which makes it
possible to use microbial abbreviations, codes and names as input. See
*Examples*.

## Usage

``` r
mo_name(x, language = get_AMR_locale(),
  keep_synonyms = getOption("AMR_keep_synonyms", FALSE), ...)

mo_fullname(x, language = get_AMR_locale(),
  keep_synonyms = getOption("AMR_keep_synonyms", FALSE), ...)

mo_shortname(x, language = get_AMR_locale(),
  keep_synonyms = getOption("AMR_keep_synonyms", FALSE), ...)

mo_subspecies(x, language = get_AMR_locale(),
  keep_synonyms = getOption("AMR_keep_synonyms", FALSE), ...)

mo_species(x, language = get_AMR_locale(),
  keep_synonyms = getOption("AMR_keep_synonyms", FALSE), ...)

mo_genus(x, language = get_AMR_locale(),
  keep_synonyms = getOption("AMR_keep_synonyms", FALSE), ...)

mo_family(x, language = get_AMR_locale(),
  keep_synonyms = getOption("AMR_keep_synonyms", FALSE), ...)

mo_order(x, language = get_AMR_locale(),
  keep_synonyms = getOption("AMR_keep_synonyms", FALSE), ...)

mo_class(x, language = get_AMR_locale(),
  keep_synonyms = getOption("AMR_keep_synonyms", FALSE), ...)

mo_phylum(x, language = get_AMR_locale(),
  keep_synonyms = getOption("AMR_keep_synonyms", FALSE), ...)

mo_kingdom(x, language = get_AMR_locale(),
  keep_synonyms = getOption("AMR_keep_synonyms", FALSE), ...)

mo_domain(x, language = get_AMR_locale(),
  keep_synonyms = getOption("AMR_keep_synonyms", FALSE), ...)

mo_type(x, language = get_AMR_locale(),
  keep_synonyms = getOption("AMR_keep_synonyms", FALSE), ...)

mo_status(x, language = get_AMR_locale(),
  keep_synonyms = getOption("AMR_keep_synonyms", FALSE), ...)

mo_pathogenicity(x, language = get_AMR_locale(),
  keep_synonyms = getOption("AMR_keep_synonyms", FALSE), ...)

mo_gramstain(x, language = get_AMR_locale(),
  keep_synonyms = getOption("AMR_keep_synonyms", FALSE), ...)

mo_is_gram_negative(x, language = get_AMR_locale(),
  keep_synonyms = getOption("AMR_keep_synonyms", FALSE), ...)

mo_is_gram_positive(x, language = get_AMR_locale(),
  keep_synonyms = getOption("AMR_keep_synonyms", FALSE), ...)

mo_is_yeast(x, language = get_AMR_locale(),
  keep_synonyms = getOption("AMR_keep_synonyms", FALSE), ...)

mo_is_intrinsic_resistant(x, ab, language = get_AMR_locale(),
  keep_synonyms = getOption("AMR_keep_synonyms", FALSE), ...)

mo_oxygen_tolerance(x, language = get_AMR_locale(),
  keep_synonyms = getOption("AMR_keep_synonyms", FALSE), ...)

mo_is_anaerobic(x, language = get_AMR_locale(),
  keep_synonyms = getOption("AMR_keep_synonyms", FALSE), ...)

mo_snomed(x, language = get_AMR_locale(),
  keep_synonyms = getOption("AMR_keep_synonyms", FALSE), ...)

mo_ref(x, language = get_AMR_locale(),
  keep_synonyms = getOption("AMR_keep_synonyms", FALSE), ...)

mo_authors(x, language = get_AMR_locale(),
  keep_synonyms = getOption("AMR_keep_synonyms", FALSE), ...)

mo_year(x, language = get_AMR_locale(),
  keep_synonyms = getOption("AMR_keep_synonyms", FALSE), ...)

mo_lpsn(x, language = get_AMR_locale(),
  keep_synonyms = getOption("AMR_keep_synonyms", FALSE), ...)

mo_mycobank(x, language = get_AMR_locale(),
  keep_synonyms = getOption("AMR_keep_synonyms", FALSE), ...)

mo_gbif(x, language = get_AMR_locale(),
  keep_synonyms = getOption("AMR_keep_synonyms", FALSE), ...)

mo_rank(x, language = get_AMR_locale(),
  keep_synonyms = getOption("AMR_keep_synonyms", FALSE), ...)

mo_taxonomy(x, language = get_AMR_locale(),
  keep_synonyms = getOption("AMR_keep_synonyms", FALSE), ...)

mo_synonyms(x, language = get_AMR_locale(),
  keep_synonyms = getOption("AMR_keep_synonyms", FALSE), ...)

mo_current(x, language = get_AMR_locale(), ...)

mo_group_members(x, language = get_AMR_locale(),
  keep_synonyms = getOption("AMR_keep_synonyms", FALSE), ...)

mo_info(x, language = get_AMR_locale(),
  keep_synonyms = getOption("AMR_keep_synonyms", FALSE), ...)

mo_url(x, open = FALSE, language = get_AMR_locale(),
  keep_synonyms = getOption("AMR_keep_synonyms", FALSE), ...)

mo_property(x, property = "fullname", language = get_AMR_locale(),
  keep_synonyms = getOption("AMR_keep_synonyms", FALSE), ...)
```

## Arguments

- x:

  Any [character](https://rdrr.io/r/base/character.html) (vector) that
  can be coerced to a valid microorganism code with
  [`as.mo()`](https://amr-for-r.org/reference/as.mo.md). Can be left
  blank for auto-guessing the column containing microorganism codes if
  used in a data set, see *Examples*.

- language:

  Language to translate text like "no growth", which defaults to the
  system language (see
  [`get_AMR_locale()`](https://amr-for-r.org/reference/translate.md)).

- keep_synonyms:

  A [logical](https://rdrr.io/r/base/logical.html) to indicate if old,
  previously valid taxonomic names must be preserved and not be
  corrected to currently accepted names. The default is `FALSE`, which
  will return a note if old taxonomic names were processed. The default
  can be set with the package option
  [`AMR_keep_synonyms`](https://amr-for-r.org/reference/AMR-options.md),
  i.e. `options(AMR_keep_synonyms = TRUE)` or
  `options(AMR_keep_synonyms = FALSE)`.

- ...:

  Other arguments passed on to
  [`as.mo()`](https://amr-for-r.org/reference/as.mo.md), such as
  'minimum_matching_score', 'ignore_pattern', and 'remove_from_input'.

- ab:

  Any (vector of) text that can be coerced to a valid antibiotic drug
  code with [`as.ab()`](https://amr-for-r.org/reference/as.ab.md).

- open:

  Browse the URL using
  [`browseURL()`](https://rdrr.io/r/utils/browseURL.html).

- property:

  One of the column names of the
  [microorganisms](https://amr-for-r.org/reference/microorganisms.md)
  data set: "mo", "fullname", "status", "kingdom", "phylum", "class",
  "order", "family", "genus", "species", "subspecies", "rank", "ref",
  "oxygen_tolerance", "source", "lpsn", "lpsn_parent",
  "lpsn_renamed_to", "mycobank", "mycobank_parent",
  "mycobank_renamed_to", "gbif", "gbif_parent", "gbif_renamed_to",
  "prevalence", or "snomed", or must be `"shortname"`.

## Value

- An [integer](https://rdrr.io/r/base/integer.html) in case of
  `mo_year()`

- An [ordered
  factor](https://rdatatable.gitlab.io/data.table/reference/fctr.html)
  in case of `mo_pathogenicity()`

- A [list](https://rdrr.io/r/base/list.html) in case of `mo_taxonomy()`,
  `mo_synonyms()`, `mo_snomed()`, and `mo_info()`

- A [logical](https://rdrr.io/r/base/logical.html) in case of
  `mo_is_anaerobic()`, `mo_is_gram_negative()`, `mo_is_gram_positive()`,
  `mo_is_intrinsic_resistant()`, and `mo_is_yeast()`

- A named [character](https://rdrr.io/r/base/character.html) in case of
  `mo_synonyms()` and `mo_url()`

- A [character](https://rdrr.io/r/base/character.html) in all other
  cases

## Details

All functions will, at default, **not** keep old taxonomic properties,
as synonyms are automatically replaced with the current taxonomy. Take
for example *Enterobacter aerogenes*, which was initially named in 1960
but renamed to *Klebsiella aerogenes* in 2017:

- `mo_genus("Enterobacter aerogenes")` will return `"Klebsiella"` (with
  a note about the renaming)

- `mo_genus("Enterobacter aerogenes", keep_synonyms = TRUE)` will return
  `"Enterobacter"` (with a once-per-session warning that the name is
  outdated)

- `mo_ref("Enterobacter aerogenes")` will return
  `"Tindall et al., 2017"` (with a note about the renaming)

- `mo_ref("Enterobacter aerogenes", keep_synonyms = TRUE)` will return
  `"Hormaeche et al., 1960"` (with a once-per-session warning that the
  name is outdated)

The short name (`mo_shortname()`) returns the first character of the
genus and the full species, such as `"E. coli"`, for species and
subspecies. Exceptions are abbreviations of staphylococci (such as
*"CoNS"*, Coagulase-Negative Staphylococci) and beta-haemolytic
streptococci (such as *"GBS"*, Group B Streptococci). Please bear in
mind that e.g. *E. coli* could mean *Escherichia coli* (kingdom of
Bacteria) as well as *Entamoeba coli* (kingdom of Protozoa). Returning
to the full name will be done using
[`as.mo()`](https://amr-for-r.org/reference/as.mo.md) internally, giving
priority to bacteria and human pathogens, i.e. `"E. coli"` will be
considered *Escherichia coli*. As a result,
`mo_fullname(mo_shortname("Entamoeba coli"))` returns
`"Escherichia coli"`.

Since the top-level of the taxonomy is sometimes referred to as
'kingdom' and sometimes as 'domain', the functions `mo_kingdom()` and
`mo_domain()` return the exact same results.

Determination of human pathogenicity (`mo_pathogenicity()`) is strongly
based on Bartlett *et al.* (2022,
[doi:10.1099/mic.0.001269](https://doi.org/10.1099/mic.0.001269) ). This
function returns a
[factor](https://rdatatable.gitlab.io/data.table/reference/fctr.html)
with the levels *Pathogenic*, *Potentially pathogenic*,
*Non-pathogenic*, and *Unknown*.

Determination of the Gram stain (`mo_gramstain()`) will be based on the
taxonomic kingdom and phylum. Originally, Cavalier-Smith defined the
so-called subkingdoms Negibacteria and Posibacteria (2002, [PMID
11837318](https://pubmed.ncbi.nlm.nih.gov/11837318/)), and only
considered these phyla as Posibacteria: Actinobacteria, Chloroflexi,
Firmicutes, and Tenericutes. These phyla were later renamed to
Actinomycetota, Chloroflexota, Bacillota, and Mycoplasmatota (2021,
[PMID 34694987](https://pubmed.ncbi.nlm.nih.gov/34694987/)). Bacteria in
these phyla are considered Gram-positive in this `AMR` package, except
for members of the class Negativicutes (within phylum Bacillota) which
are Gram-negative. All other bacteria are considered Gram-negative.
Species outside the kingdom of Bacteria will return a value `NA`.
Functions `mo_is_gram_negative()` and `mo_is_gram_positive()` always
return `TRUE` or `FALSE` (or `NA` when the input is `NA` or the MO code
is `UNKNOWN`), thus always return `FALSE` for species outside the
taxonomic kingdom of Bacteria.

Determination of yeasts (`mo_is_yeast()`) will be based on the taxonomic
kingdom and class. *Budding yeasts* are yeasts that reproduce asexually
through a process called budding, where a new cell develops from a small
protrusion on the parent cell. Taxonomically, these are members of the
phylum Ascomycota, class Saccharomycetes (also called Hemiascomycetes)
or Pichiomycetes. *True yeasts* quite specifically refers to yeasts in
the underlying order Saccharomycetales (such as *Saccharomyces
cerevisiae*). Thus, for all microorganisms that are member of the
taxonomic class Saccharomycetes or Pichiomycetes, the function will
return `TRUE`. It returns `FALSE` otherwise (or `NA` when the input is
`NA` or the MO code is `UNKNOWN`).

Determination of intrinsic resistance (`mo_is_intrinsic_resistant()`)
will be based on the
[intrinsic_resistant](https://amr-for-r.org/reference/intrinsic_resistant.md)
data set, which is based on ['EUCAST Expected Resistant Phenotypes'
v1.2](https://www.eucast.org/bacteria/important-additional-information/expert-rules/)
(2023). The `mo_is_intrinsic_resistant()` function can be vectorised
over both argument `x` (input for microorganisms) and `ab` (input for
antimicrobials).

Determination of bacterial oxygen tolerance (`mo_oxygen_tolerance()`)
will be based on BacDive, see *Source*. The function `mo_is_anaerobic()`
only returns `TRUE` if the oxygen tolerance is `"anaerobe"`, indicting
an obligate anaerobic species or genus. It always returns `FALSE` for
species outside the taxonomic kingdom of Bacteria.

The function `mo_url()` will return the direct URL to the online
database entry, which also shows the scientific reference of the
concerned species. [This MycoBank URL](https://www.mycobank.org) will be
used for fungi wherever available , [this LPSN
URL](https://www.mycobank.org) for bacteria wherever available, and
[this GBIF link](https://www.gbif.org) otherwise.

SNOMED codes (`mo_snomed()`) was last updated on July 16th, 2024. See
*Source* and the
[microorganisms](https://amr-for-r.org/reference/microorganisms.md) data
set for more info.

Old taxonomic names (so-called 'synonyms') can be retrieved with
`mo_synonyms()` (which will have the scientific reference as
[name](https://rdrr.io/r/base/names.html)), the current taxonomic name
can be retrieved with `mo_current()`. Both functions return full names.

All output [will be
translated](https://amr-for-r.org/reference/translate.md) where
possible.

## Matching Score for Microorganisms

This function uses [`as.mo()`](https://amr-for-r.org/reference/as.mo.md)
internally, which uses an advanced algorithm to translate arbitrary user
input to valid taxonomy using a so-called matching score. You can read
about this public algorithm on the [MO matching score
page](https://amr-for-r.org/reference/mo_matching_score.md).

## Source

- Berends MS *et al.* (2022). **AMR: An R Package for Working with
  Antimicrobial Resistance Data**. *Journal of Statistical Software*,
  104(3), 1-31;
  [doi:10.18637/jss.v104.i03](https://doi.org/10.18637/jss.v104.i03)

- Parte, AC *et al.* (2020). **List of Prokaryotic names with Standing
  in Nomenclature (LPSN) moves to the DSMZ.** International Journal of
  Systematic and Evolutionary Microbiology, 70, 5607-5612;
  [doi:10.1099/ijsem.0.004332](https://doi.org/10.1099/ijsem.0.004332) .
  Accessed from <https://lpsn.dsmz.de> on June 24th, 2024.

- Vincent, R *et al* (2013). **MycoBank gearing up for new horizons.**
  IMA Fungus, 4(2), 371-9;
  [doi:10.5598/imafungus.2013.04.02.16](https://doi.org/10.5598/imafungus.2013.04.02.16)
  . Accessed from <https://www.mycobank.org> on June 24th, 2024.

- GBIF Secretariat (2023). GBIF Backbone Taxonomy. Checklist dataset
  [doi:10.15468/39omei](https://doi.org/10.15468/39omei) . Accessed from
  <https://www.gbif.org> on June 24th, 2024.

- Reimer, LC *et al.* (2022). ***BacDive* in 2022: the knowledge base
  for standardized bacterial and archaeal data.** Nucleic Acids Res.,
  50(D1):D741-D74;
  [doi:10.1093/nar/gkab961](https://doi.org/10.1093/nar/gkab961) .
  Accessed from <https://bacdive.dsmz.de> on July 16th, 2024.

- Public Health Information Network Vocabulary Access and Distribution
  System (PHIN VADS). US Edition of SNOMED CT from 1 September 2020.
  Value Set Name 'Microorganism', OID 2.16.840.1.114222.4.11.1009 (v12).
  URL: <https://www.cdc.gov/phin/php/phinvads/>

- Bartlett A *et al.* (2022). **A comprehensive list of bacterial
  pathogens infecting humans** *Microbiology* 168:001269;
  [doi:10.1099/mic.0.001269](https://doi.org/10.1099/mic.0.001269)

## Download Our Reference Data

All reference data sets in the AMR package - including information on
microorganisms, antimicrobials, and clinical breakpoints - are freely
available for download in multiple formats: R, MS Excel, Apache Feather,
Apache Parquet, SPSS, and Stata.

For maximum compatibility, we also provide machine-readable,
tab-separated plain text files suitable for use in any software,
including laboratory information systems.

Visit [our website for direct download
links](https://amr-for-r.org/articles/datasets.html), or explore the
actual files in [our GitHub
repository](https://github.com/msberends/AMR/tree/main/data-raw/datasets).

## See also

Data set
[microorganisms](https://amr-for-r.org/reference/microorganisms.md)

## Examples

``` r
# taxonomic tree -----------------------------------------------------------

mo_kingdom("Klebsiella pneumoniae")
#> [1] "Bacteria"
mo_phylum("Klebsiella pneumoniae")
#> [1] "Pseudomonadota"
mo_class("Klebsiella pneumoniae")
#> [1] "Gammaproteobacteria"
mo_order("Klebsiella pneumoniae")
#> [1] "Enterobacterales"
mo_family("Klebsiella pneumoniae")
#> [1] "Enterobacteriaceae"
mo_genus("Klebsiella pneumoniae")
#> [1] "Klebsiella"
mo_species("Klebsiella pneumoniae")
#> [1] "pneumoniae"
mo_subspecies("Klebsiella pneumoniae")
#> [1] ""


# full names and short names -----------------------------------------------

mo_name("Klebsiella pneumoniae")
#> [1] "Klebsiella pneumoniae"
mo_fullname("Klebsiella pneumoniae")
#> [1] "Klebsiella pneumoniae"
mo_shortname("Klebsiella pneumoniae")
#> [1] "K. pneumoniae"


# other properties ---------------------------------------------------------

mo_pathogenicity("Klebsiella pneumoniae")
#> [1] Pathogenic
#> Levels: Pathogenic < Potentially pathogenic < Non-pathogenic < Unknown
mo_gramstain("Klebsiella pneumoniae")
#> [1] "Gram-negative"
mo_snomed("Klebsiella pneumoniae")
#> [[1]]
#> [1] "1098101000112102" "446870005"        "1098201000112108" "409801009"       
#> [5] "56415008"         "714315002"        "713926009"       
#> 
mo_type("Klebsiella pneumoniae")
#> [1] "Bacteria"
mo_rank("Klebsiella pneumoniae")
#> [1] "species"
mo_url("Klebsiella pneumoniae")
#>                                Klebsiella pneumoniae 
#> "https://lpsn.dsmz.de/species/klebsiella-pneumoniae" 
mo_is_yeast(c("Candida", "Trichophyton", "Klebsiella"))
#> [1]  TRUE FALSE FALSE

mo_group_members(c(
  "Streptococcus group A",
  "Streptococcus group C",
  "Streptococcus group G",
  "Streptococcus group L"
))
#> $`Streptococcus Group A`
#> [1] "Streptococcus pyogenes"
#> 
#> $`Streptococcus Group C`
#> [1] "Streptococcus dysgalactiae"             
#> [2] "Streptococcus dysgalactiae dysgalactiae"
#> [3] "Streptococcus dysgalactiae equisimilis" 
#> [4] "Streptococcus equi"                     
#> [5] "Streptococcus equi equi"                
#> [6] "Streptococcus equi ruminatorum"         
#> [7] "Streptococcus equi zooepidemicus"       
#> 
#> $`Streptococcus Group G`
#> [1] "Streptococcus canis"                    
#> [2] "Streptococcus dysgalactiae"             
#> [3] "Streptococcus dysgalactiae dysgalactiae"
#> [4] "Streptococcus dysgalactiae equisimilis" 
#> 
#> $`Streptococcus Group L`
#> [1] "Streptococcus dysgalactiae"             
#> [2] "Streptococcus dysgalactiae dysgalactiae"
#> [3] "Streptococcus dysgalactiae equisimilis" 
#> 


# scientific reference -----------------------------------------------------

mo_ref("Klebsiella aerogenes")
#> [1] "Tindall et al., 2017"
mo_authors("Klebsiella aerogenes")
#> [1] "Tindall et al."
mo_year("Klebsiella aerogenes")
#> [1] 2017
mo_synonyms("Klebsiella aerogenes")
#>   Hormaeche et al., 1960     Bascomb et al., 1971 
#> "Enterobacter aerogenes"     "Klebsiella mobilis" 
mo_lpsn("Klebsiella aerogenes")
#> [1] "777146"
mo_gbif("Klebsiella aerogenes")
#> [1] "9281703"
mo_mycobank("Candida albicans")
#> [1] "256187"
mo_mycobank("Candida krusei")
#> [1] "337013"
mo_mycobank("Candida krusei", keep_synonyms = TRUE)
#> Warning: Function `as.mo()` returned one old taxonomic name. Use `as.mo(...,
#> keep_synonyms = FALSE)` to clean the input to currently accepted taxonomic
#> names, or set the R option `AMR_keep_synonyms` to `FALSE`. This warning
#> will be shown once per session.
#> [1] "268707"


# abbreviations known in the field -----------------------------------------

mo_genus("MRSA")
#> [1] "Staphylococcus"
mo_species("MRSA")
#> [1] "aureus"
mo_shortname("VISA")
#> [1] "S. aureus"
mo_gramstain("VISA")
#> [1] "Gram-positive"

mo_genus("EHEC")
#> [1] "Escherichia"
mo_species("EIEC")
#> [1] "coli"
mo_name("UPEC")
#> [1] "Escherichia coli"


# known subspecies ---------------------------------------------------------

mo_fullname("K. pneu rh")
#> [1] "Klebsiella pneumoniae rhinoscleromatis"
mo_shortname("K. pneu rh")
#> [1] "K. pneumoniae"

# \donttest{
# Becker classification, see ?as.mo ----------------------------------------

mo_fullname("Staph epidermidis")
#> [1] "Staphylococcus epidermidis"
mo_fullname("Staph epidermidis", Becker = TRUE)
#> [1] "Coagulase-negative Staphylococcus (CoNS)"
mo_shortname("Staph epidermidis")
#> [1] "S. epidermidis"
mo_shortname("Staph epidermidis", Becker = TRUE)
#> [1] "CoNS"


# Lancefield classification, see ?as.mo ------------------------------------

mo_fullname("Strep agalactiae")
#> [1] "Streptococcus agalactiae"
mo_fullname("Strep agalactiae", Lancefield = TRUE)
#> [1] "Streptococcus Group B"
mo_shortname("Strep agalactiae")
#> [1] "S. agalactiae"
mo_shortname("Strep agalactiae", Lancefield = TRUE)
#> [1] "GBS"


# language support  --------------------------------------------------------

mo_gramstain("Klebsiella pneumoniae", language = "de") # German
#> [1] "Gramnegativ"
mo_gramstain("Klebsiella pneumoniae", language = "nl") # Dutch
#> [1] "Gram-negatief"
mo_gramstain("Klebsiella pneumoniae", language = "es") # Spanish
#> [1] "Gram negativo"
mo_gramstain("Klebsiella pneumoniae", language = "el") # Greek
#> [1] "Αρνητικό κατά Gram"
mo_gramstain("Klebsiella pneumoniae", language = "uk") # Ukrainian
#> [1] "Грамнегативні"

# mo_type is equal to mo_kingdom, but mo_kingdom will remain untranslated
mo_kingdom("Klebsiella pneumoniae")
#> [1] "Bacteria"
mo_type("Klebsiella pneumoniae")
#> [1] "Bacteria"
mo_kingdom("Klebsiella pneumoniae", language = "zh") # Chinese, no effect
#> [1] "Bacteria"
mo_type("Klebsiella pneumoniae", language = "zh") # Chinese, translated
#> [1] "细菌"

mo_fullname("S. pyogenes", Lancefield = TRUE, language = "de")
#> [1] "Streptococcus Gruppe A"
mo_fullname("S. pyogenes", Lancefield = TRUE, language = "uk")
#> [1] "Streptococcus Група A"


# other --------------------------------------------------------------------

# gram stains and intrinsic resistance can be used as a filter in dplyr verbs
if (require("dplyr")) {
  example_isolates %>%
    filter(mo_is_gram_positive()) %>%
    count(mo_genus(), sort = TRUE)
}
#> ℹ Using column 'mo' as input for `mo_is_gram_positive()`
#> ℹ Using column 'mo' as input for `mo_genus()`
#> # A tibble: 18 × 2
#>    `mo_genus()`        n
#>    <chr>           <int>
#>  1 Staphylococcus    840
#>  2 Streptococcus     275
#>  3 Enterococcus       83
#>  4 Corynebacterium    17
#>  5 Micrococcus         6
#>  6 Gemella             3
#>  7 Aerococcus          2
#>  8 Cutibacterium       1
#>  9 Dermabacter         1
#> 10 Fusibacter          1
#> 11 Globicatella        1
#> 12 Granulicatella      1
#> 13 Lactobacillus       1
#> 14 Leuconostoc         1
#> 15 Listeria            1
#> 16 Paenibacillus       1
#> 17 Rothia              1
#> 18 Schaalia            1
if (require("dplyr")) {
  example_isolates %>%
    filter(mo_is_intrinsic_resistant(ab = "vanco")) %>%
    count(mo_genus(), sort = TRUE)
}
#> ℹ Using column 'mo' as input for `mo_is_intrinsic_resistant()`
#> ℹ Using column 'mo' as input for `mo_genus()`
#> # A tibble: 19 × 2
#>    `mo_genus()`         n
#>    <chr>            <int>
#>  1 Escherichia        467
#>  2 Klebsiella          77
#>  3 Proteus             39
#>  4 Pseudomonas         30
#>  5 Serratia            25
#>  6 Enterobacter        23
#>  7 Citrobacter         11
#>  8 Haemophilus          9
#>  9 Acinetobacter        6
#> 10 Morganella           6
#> 11 Pantoea              4
#> 12 Salmonella           3
#> 13 Neisseria            2
#> 14 Stenotrophomonas     2
#> 15 Campylobacter        1
#> 16 Enterococcus         1
#> 17 Hafnia               1
#> 18 Leuconostoc          1
#> 19 Pseudescherichia     1

# get a list with the complete taxonomy (from kingdom to subspecies)
mo_taxonomy("Klebsiella pneumoniae")
#> $kingdom
#> [1] "Bacteria"
#> 
#> $phylum
#> [1] "Pseudomonadota"
#> 
#> $class
#> [1] "Gammaproteobacteria"
#> 
#> $order
#> [1] "Enterobacterales"
#> 
#> $family
#> [1] "Enterobacteriaceae"
#> 
#> $genus
#> [1] "Klebsiella"
#> 
#> $species
#> [1] "pneumoniae"
#> 
#> $subspecies
#> [1] ""
#> 

# get a list with the taxonomy, the authors, Gram-stain,
# SNOMED codes, and URL to the online database
mo_info("Klebsiella pneumoniae")
#> $mo
#> [1] "B_KLBSL_PNMN"
#> 
#> $rank
#> [1] "species"
#> 
#> $kingdom
#> [1] "Bacteria"
#> 
#> $phylum
#> [1] "Pseudomonadota"
#> 
#> $class
#> [1] "Gammaproteobacteria"
#> 
#> $order
#> [1] "Enterobacterales"
#> 
#> $family
#> [1] "Enterobacteriaceae"
#> 
#> $genus
#> [1] "Klebsiella"
#> 
#> $species
#> [1] "pneumoniae"
#> 
#> $subspecies
#> [1] ""
#> 
#> $status
#> [1] "accepted"
#> 
#> $synonyms
#> NULL
#> 
#> $gramstain
#> [1] "Gram-negative"
#> 
#> $oxygen_tolerance
#> [1] "facultative anaerobe"
#> 
#> $url
#> [1] "https://lpsn.dsmz.de/species/klebsiella-pneumoniae"
#> 
#> $ref
#> [1] "Trevisan, 1887"
#> 
#> $snomed
#> [1] "1098101000112102" "446870005"        "1098201000112108" "409801009"       
#> [5] "56415008"         "714315002"        "713926009"       
#> 
#> $lpsn
#> [1] "777151"
#> 
#> $mycobank
#> [1] NA
#> 
#> $gbif
#> [1] "3221874"
#> 
#> $group_members
#> character(0)
#> 
# }
```
