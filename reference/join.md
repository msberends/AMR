# Join [microorganisms](https://amr-for-r.org/reference/microorganisms.md) to a Data Set

Join the data set
[microorganisms](https://amr-for-r.org/reference/microorganisms.md)
easily to an existing data set or to a
[character](https://rdrr.io/r/base/character.html) vector.

## Usage

``` r
inner_join_microorganisms(x, by = NULL, suffix = c("2", ""), ...)

left_join_microorganisms(x, by = NULL, suffix = c("2", ""), ...)

right_join_microorganisms(x, by = NULL, suffix = c("2", ""), ...)

full_join_microorganisms(x, by = NULL, suffix = c("2", ""), ...)

semi_join_microorganisms(x, by = NULL, ...)

anti_join_microorganisms(x, by = NULL, ...)
```

## Arguments

- x:

  Existing data set to join, or
  [character](https://rdrr.io/r/base/character.html) vector. In case of
  a [character](https://rdrr.io/r/base/character.html) vector, the
  resulting [data.frame](https://rdrr.io/r/base/data.frame.html) will
  contain a column 'x' with these values.

- by:

  A variable to join by - if left empty will search for a column with
  class [`mo`](https://amr-for-r.org/reference/as.mo.md) (created with
  [`as.mo()`](https://amr-for-r.org/reference/as.mo.md)) or will be
  `"mo"` if that column name exists in `x`, could otherwise be a column
  name of `x` with values that exist in `microorganisms$mo` (such as
  `by = "bacteria_id"`), or another column in
  [microorganisms](https://amr-for-r.org/reference/microorganisms.md)
  (but then it should be named, like
  `by = c("bacteria_id" = "fullname")`).

- suffix:

  If there are non-joined duplicate variables in `x` and `y`, these
  suffixes will be added to the output to disambiguate them. Should be a
  [character](https://rdrr.io/r/base/character.html) vector of length 2.

- ...:

  Ignored, only in place to allow future extensions.

## Value

a [data.frame](https://rdrr.io/r/base/data.frame.html)

## Details

**Note:** As opposed to the `join()` functions of `dplyr`,
[character](https://rdrr.io/r/base/character.html) vectors are supported
and at default existing columns will get a suffix `"2"` and the newly
joined columns will not get a suffix.

If the `dplyr` package is installed, their join functions will be used.
Otherwise, the much slower
[`merge()`](https://rdrr.io/pkg/data.table/man/merge.html) and
[`interaction()`](https://rdrr.io/r/base/interaction.html) functions
from base R will be used.

## Examples

``` r
left_join_microorganisms(as.mo("K. pneumoniae"))
#> # A tibble: 1 × 28
#>   mo           fullname    status domain kingdom phylum class order family genus
#>   <mo>         <chr>       <chr>  <chr>  <chr>   <chr>  <chr> <chr> <chr>  <chr>
#> 1 B_KLBSL_PNMN Klebsiella… accep… Bacte… Pseudo… Pseud… Gamm… Ente… Enter… Kleb…
#> # ℹ 18 more variables: species <chr>, subspecies <chr>, rank <chr>, ref <chr>,
#> #   oxygen_tolerance <chr>, morphology <chr>, source <chr>, lpsn <chr>,
#> #   lpsn_parent <chr>, lpsn_renamed_to <chr>, mycobank <chr>,
#> #   mycobank_parent <chr>, mycobank_renamed_to <chr>, gbif <chr>,
#> #   gbif_parent <chr>, gbif_renamed_to <chr>, prevalence <dbl>, snomed <list>
left_join_microorganisms("B_KLBSL_PNMN")
#> # A tibble: 1 × 28
#>   mo           fullname    status domain kingdom phylum class order family genus
#>   <mo>         <chr>       <chr>  <chr>  <chr>   <chr>  <chr> <chr> <chr>  <chr>
#> 1 B_KLBSL_PNMN Klebsiella… accep… Bacte… Pseudo… Pseud… Gamm… Ente… Enter… Kleb…
#> # ℹ 18 more variables: species <chr>, subspecies <chr>, rank <chr>, ref <chr>,
#> #   oxygen_tolerance <chr>, morphology <chr>, source <chr>, lpsn <chr>,
#> #   lpsn_parent <chr>, lpsn_renamed_to <chr>, mycobank <chr>,
#> #   mycobank_parent <chr>, mycobank_renamed_to <chr>, gbif <chr>,
#> #   gbif_parent <chr>, gbif_renamed_to <chr>, prevalence <dbl>, snomed <list>

df <- data.frame(
  date = seq(
    from = as.Date("2018-01-01"),
    to = as.Date("2018-01-07"),
    by = 1
  ),
  bacteria = as.mo(c(
    "S. aureus", "MRSA", "MSSA", "STAAUR",
    "E. coli", "E. coli", "E. coli"
  )),
  stringsAsFactors = FALSE
)
colnames(df)
#> [1] "date"     "bacteria"

df_joined <- left_join_microorganisms(df, "bacteria")
colnames(df_joined)
#>  [1] "date"                "bacteria"            "fullname"           
#>  [4] "status"              "domain"              "kingdom"            
#>  [7] "phylum"              "class"               "order"              
#> [10] "family"              "genus"               "species"            
#> [13] "subspecies"          "rank"                "ref"                
#> [16] "oxygen_tolerance"    "morphology"          "source"             
#> [19] "lpsn"                "lpsn_parent"         "lpsn_renamed_to"    
#> [22] "mycobank"            "mycobank_parent"     "mycobank_renamed_to"
#> [25] "gbif"                "gbif_parent"         "gbif_renamed_to"    
#> [28] "prevalence"          "snomed"             

# \donttest{
if (require("dplyr")) {
  example_isolates %>%
    left_join_microorganisms() %>%
    colnames()
}
#> Joining, by = "mo"
#>  [1] "date"                "patient"             "age"                
#>  [4] "gender"              "ward"                "mo"                 
#>  [7] "PEN"                 "OXA"                 "FLC"                
#> [10] "AMX"                 "AMC"                 "AMP"                
#> [13] "TZP"                 "CZO"                 "FEP"                
#> [16] "CXM"                 "FOX"                 "CTX"                
#> [19] "CAZ"                 "CRO"                 "GEN"                
#> [22] "TOB"                 "AMK"                 "KAN"                
#> [25] "TMP"                 "SXT"                 "NIT"                
#> [28] "FOS"                 "LNZ"                 "CIP"                
#> [31] "MFX"                 "VAN"                 "TEC"                
#> [34] "TCY"                 "TGC"                 "DOX"                
#> [37] "ERY"                 "CLI"                 "AZM"                
#> [40] "IPM"                 "MEM"                 "MTR"                
#> [43] "CHL"                 "COL"                 "MUP"                
#> [46] "RIF"                 "fullname"            "status"             
#> [49] "domain"              "kingdom"             "phylum"             
#> [52] "class"               "order"               "family"             
#> [55] "genus"               "species"             "subspecies"         
#> [58] "rank"                "ref"                 "oxygen_tolerance"   
#> [61] "morphology"          "source"              "lpsn"               
#> [64] "lpsn_parent"         "lpsn_renamed_to"     "mycobank"           
#> [67] "mycobank_parent"     "mycobank_renamed_to" "gbif"               
#> [70] "gbif_parent"         "gbif_renamed_to"     "prevalence"         
#> [73] "snomed"             
# }
```
