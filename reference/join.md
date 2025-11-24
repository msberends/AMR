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
[`merge()`](https://rdatatable.gitlab.io/data.table/reference/merge.html)
and [`interaction()`](https://rdrr.io/r/base/interaction.html) functions
from base R will be used.

## Examples

``` r
left_join_microorganisms(as.mo("K. pneumoniae"))
#> # A tibble: 1 × 26
#>   mo           fullname   status kingdom phylum class order family genus species
#>   <mo>         <chr>      <chr>  <chr>   <chr>  <chr> <chr> <chr>  <chr> <chr>  
#> 1 B_KLBSL_PNMN Klebsiell… accep… Bacter… Pseud… Gamm… Ente… Enter… Kleb… pneumo…
#> # ℹ 16 more variables: subspecies <chr>, rank <chr>, ref <chr>,
#> #   oxygen_tolerance <chr>, source <chr>, lpsn <chr>, lpsn_parent <chr>,
#> #   lpsn_renamed_to <chr>, mycobank <chr>, mycobank_parent <chr>,
#> #   mycobank_renamed_to <chr>, gbif <chr>, gbif_parent <chr>,
#> #   gbif_renamed_to <chr>, prevalence <dbl>, snomed <list>
left_join_microorganisms("B_KLBSL_PNMN")
#> # A tibble: 1 × 26
#>   mo           fullname   status kingdom phylum class order family genus species
#>   <mo>         <chr>      <chr>  <chr>   <chr>  <chr> <chr> <chr>  <chr> <chr>  
#> 1 B_KLBSL_PNMN Klebsiell… accep… Bacter… Pseud… Gamm… Ente… Enter… Kleb… pneumo…
#> # ℹ 16 more variables: subspecies <chr>, rank <chr>, ref <chr>,
#> #   oxygen_tolerance <chr>, source <chr>, lpsn <chr>, lpsn_parent <chr>,
#> #   lpsn_renamed_to <chr>, mycobank <chr>, mycobank_parent <chr>,
#> #   mycobank_renamed_to <chr>, gbif <chr>, gbif_parent <chr>,
#> #   gbif_renamed_to <chr>, prevalence <dbl>, snomed <list>

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
#>  [4] "status"              "kingdom"             "phylum"             
#>  [7] "class"               "order"               "family"             
#> [10] "genus"               "species"             "subspecies"         
#> [13] "rank"                "ref"                 "oxygen_tolerance"   
#> [16] "source"              "lpsn"                "lpsn_parent"        
#> [19] "lpsn_renamed_to"     "mycobank"            "mycobank_parent"    
#> [22] "mycobank_renamed_to" "gbif"                "gbif_parent"        
#> [25] "gbif_renamed_to"     "prevalence"          "snomed"             

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
#> [49] "kingdom"             "phylum"              "class"              
#> [52] "order"               "family"              "genus"              
#> [55] "species"             "subspecies"          "rank"               
#> [58] "ref"                 "oxygen_tolerance"    "source"             
#> [61] "lpsn"                "lpsn_parent"         "lpsn_renamed_to"    
#> [64] "mycobank"            "mycobank_parent"     "mycobank_renamed_to"
#> [67] "gbif"                "gbif_parent"         "gbif_renamed_to"    
#> [70] "prevalence"          "snomed"             
# }
```
