# Determine Clinical or Epidemic Episodes

These functions determine which items in a vector can be considered (the
start of) a new episode. This can be used to determine clinical episodes
for any epidemiological analysis. The `get_episode()` function returns
the index number of the episode per group, while the `is_new_episode()`
function returns `TRUE` for every new `get_episode()` index. Both
absolute and relative episode determination are supported.

## Usage

``` r
get_episode(x, episode_days = NULL, case_free_days = NULL, ...)

is_new_episode(x, episode_days = NULL, case_free_days = NULL, ...)
```

## Arguments

- x:

  Vector of dates (class `Date` or `POSIXt`), will be sorted internally
  to determine episodes.

- episode_days:

  Episode length in days to specify the time period after which a new
  episode begins, can also be less than a day or `Inf`, see *Details*.

- case_free_days:

  (inter-epidemic) interval length in days after which a new episode
  will start, can also be less than a day or `Inf`, see *Details*.

- ...:

  Ignored, only in place to allow future extensions.

## Value

- `get_episode()`: an [integer](https://rdrr.io/r/base/integer.html)
  vector

- `is_new_episode()`: a [logical](https://rdrr.io/r/base/logical.html)
  vector

## Details

Episodes can be determined in two ways: absolute and relative.

1.  Absolute

    This method uses `episode_days` to define an episode length in days,
    after which a new episode will start. A common use case in AMR data
    analysis is microbial epidemiology: episodes of *S. aureus*
    bacteraemia in ICU patients for example. The episode length could
    then be 30 days, so that new *S. aureus* isolates after an ICU
    episode of 30 days will be considered a different (or new) episode.

    Thus, this method counts **since the start of the previous
    episode**.

2.  Relative

    This method uses `case_free_days` to quantify the duration of
    case-free days (the inter-epidemic interval), after which a new
    episode will start. A common use case is infectious disease
    epidemiology: episodes of norovirus outbreaks in a hospital for
    example. The case-free period could then be 14 days, so that new
    norovirus cases after that time will be considered a different (or
    new) episode.

    Thus, this methods counts **since the last case in the previous
    episode**.

In a table:

|            |                          |                            |
|------------|--------------------------|----------------------------|
| Date       | Using `episode_days = 7` | Using `case_free_days = 7` |
| 2023-01-01 | 1                        | 1                          |
| 2023-01-02 | 1                        | 1                          |
| 2023-01-05 | 1                        | 1                          |
| 2023-01-08 | 2\*\*                    | 1                          |
| 2023-02-21 | 3                        | 2\*\*\*                    |
| 2023-02-22 | 3                        | 2                          |
| 2023-02-23 | 3                        | 2                          |
| 2023-02-24 | 3                        | 2                          |
| 2023-03-01 | 4                        | 2                          |

\*\* This marks the start of a new episode, because 8 January 2023 is
more than 7 days since the start of the previous episode (1 January
2023).  
\*\*\* This marks the start of a new episode, because 21 January 2023 is
more than 7 days since the last case in the previous episode (8 January
2023).

Either `episode_days` or `case_free_days` must be provided in the
function.

### Difference between `get_episode()` and `is_new_episode()`

The `get_episode()` function returns the index number of the episode, so
all cases/patients/isolates in the first episode will have the number 1,
all cases/patients/isolates in the second episode will have the number
2, etc.

The `is_new_episode()` function on the other hand, returns `TRUE` for
every new `get_episode()` index.

To specify, when setting `episode_days = 365` (using method 1 as
explained above), this is how the two functions differ:

|         |            |                 |                    |
|---------|------------|-----------------|--------------------|
| patient | date       | `get_episode()` | `is_new_episode()` |
| A       | 2019-01-01 | 1               | TRUE               |
| A       | 2019-03-01 | 1               | FALSE              |
| A       | 2021-01-01 | 2               | TRUE               |
| B       | 2008-01-01 | 1               | TRUE               |
| B       | 2008-01-01 | 1               | FALSE              |
| C       | 2020-01-01 | 1               | TRUE               |

### Other

The
[`first_isolate()`](https://amr-for-r.org/reference/first_isolate.md)
function is a wrapper around the `is_new_episode()` function, but is
more efficient for data sets containing microorganism codes or names and
allows for different isolate selection methods.

The `dplyr` package is not required for these functions to work, but
these episode functions do support [variable
grouping](https://dplyr.tidyverse.org/reference/group_by.html) and work
conveniently inside `dplyr` verbs such as
[`filter()`](https://dplyr.tidyverse.org/reference/filter.html),
[`mutate()`](https://dplyr.tidyverse.org/reference/mutate.html) and
[`summarise()`](https://dplyr.tidyverse.org/reference/summarise.html).

## See also

[`first_isolate()`](https://amr-for-r.org/reference/first_isolate.md)

## Examples

``` r
# difference between absolute and relative determination of episodes:
x <- data.frame(dates = as.Date(c(
  "2021-01-01",
  "2021-01-02",
  "2021-01-05",
  "2021-01-08",
  "2021-02-21",
  "2021-02-22",
  "2021-02-23",
  "2021-02-24",
  "2021-03-01",
  "2021-03-01"
)))
x$absolute <- get_episode(x$dates, episode_days = 7)
x$relative <- get_episode(x$dates, case_free_days = 7)
x
#>         dates absolute relative
#> 1  2021-01-01        1        1
#> 2  2021-01-02        1        1
#> 3  2021-01-05        1        1
#> 4  2021-01-08        2        1
#> 5  2021-02-21        3        2
#> 6  2021-02-22        3        2
#> 7  2021-02-23        3        2
#> 8  2021-02-24        3        2
#> 9  2021-03-01        4        2
#> 10 2021-03-01        4        2


# `example_isolates` is a data set available in the AMR package.
# See ?example_isolates
df <- example_isolates[sample(seq_len(2000), size = 100), ]

get_episode(df$date, episode_days = 60) # indices
#>   [1] 10 47 27 19 35 45 20 18 28  1 20 34 11 31  9 44  6  5 27 42  2 19 33 39  7
#>  [26] 18 21 35  7 45 12  4 16 29 32 23 30 37 25 38 22 41 30 27 22 47 33 10 14 29
#>  [51] 44 40 17 42 35  5  6 13  8 38 45 37  1 43 45 10 46 38 19 43 10 21 44 37 44
#>  [76]  5 33  3 10 43 36 31 27  3 43 13 15 23 32 26 41 13 16 37 16 18 24 26 24 27
is_new_episode(df$date, episode_days = 60) # TRUE/FALSE
#>   [1]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE FALSE  TRUE
#>  [13]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE FALSE  TRUE  TRUE FALSE  TRUE  TRUE
#>  [25]  TRUE FALSE  TRUE FALSE FALSE FALSE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
#>  [37]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE FALSE FALSE FALSE FALSE FALSE FALSE
#>  [49]  TRUE FALSE FALSE  TRUE  TRUE FALSE FALSE FALSE FALSE  TRUE  TRUE FALSE
#>  [61] FALSE FALSE FALSE  TRUE FALSE FALSE  TRUE FALSE FALSE FALSE FALSE FALSE
#>  [73] FALSE FALSE FALSE FALSE FALSE  TRUE FALSE FALSE  TRUE FALSE FALSE FALSE
#>  [85] FALSE FALSE  TRUE FALSE FALSE  TRUE FALSE FALSE FALSE FALSE FALSE FALSE
#>  [97]  TRUE FALSE FALSE FALSE

# filter on results from the third 60-day episode only, using base R
df[which(get_episode(df$date, 60) == 3), ]
#> # A tibble: 2 × 46
#>   date       patient   age gender ward  mo           PEN   OXA   FLC   AMX  
#>   <date>     <chr>   <dbl> <chr>  <chr> <mo>         <sir> <sir> <sir> <sir>
#> 1 2002-08-14 785317     51 F      ICU   B_ESCHR_COLI   R     NA    NA    NA 
#> 2 2002-06-22 FD8039     75 F      ICU   B_ESCHR_COLI   R     NA    NA    NA 
#> # ℹ 36 more variables: AMC <sir>, AMP <sir>, TZP <sir>, CZO <sir>, FEP <sir>,
#> #   CXM <sir>, FOX <sir>, CTX <sir>, CAZ <sir>, CRO <sir>, GEN <sir>,
#> #   TOB <sir>, AMK <sir>, KAN <sir>, TMP <sir>, SXT <sir>, NIT <sir>,
#> #   FOS <sir>, LNZ <sir>, CIP <sir>, MFX <sir>, VAN <sir>, TEC <sir>,
#> #   TCY <sir>, TGC <sir>, DOX <sir>, ERY <sir>, CLI <sir>, AZM <sir>,
#> #   IPM <sir>, MEM <sir>, MTR <sir>, CHL <sir>, COL <sir>, MUP <sir>, RIF <sir>

# the functions also work for less than a day, e.g. to include one per hour:
get_episode(
  c(
    Sys.time(),
    Sys.time() + 60 * 60
  ),
  episode_days = 1 / 24
)
#> [1] 1 2

# \donttest{
if (require("dplyr")) {
  # is_new_episode() can also be used in dplyr verbs to determine patient
  # episodes based on any (combination of) grouping variables:
  df %>%
    mutate(condition = sample(
      x = c("A", "B", "C"),
      size = 100,
      replace = TRUE
    )) %>%
    group_by(patient, condition) %>%
    mutate(new_episode = is_new_episode(date, 365)) %>%
    select(patient, date, condition, new_episode) %>%
    arrange(patient, condition, date)
}
#> # A tibble: 100 × 4
#> # Groups:   patient, condition [100]
#>    patient date       condition new_episode
#>    <chr>   <date>     <chr>     <lgl>      
#>  1 001213  2009-08-03 A         TRUE       
#>  2 067927  2002-01-07 C         TRUE       
#>  3 078381  2014-07-17 C         TRUE       
#>  4 0DBF93  2015-12-03 C         TRUE       
#>  5 105248  2005-06-16 C         TRUE       
#>  6 110918  2016-05-08 C         TRUE       
#>  7 130252  2002-04-08 A         TRUE       
#>  8 16F0F7  2010-01-17 C         TRUE       
#>  9 189795  2008-10-19 B         TRUE       
#> 10 207752  2008-07-11 B         TRUE       
#> # ℹ 90 more rows

if (require("dplyr")) {
  df %>%
    group_by(ward, patient) %>%
    transmute(date,
      patient,
      new_index = get_episode(date, 60),
      new_logical = is_new_episode(date, 60)
    ) %>%
    arrange(patient, ward, date)
}
#> # A tibble: 100 × 5
#> # Groups:   ward, patient [95]
#>    ward     date       patient new_index new_logical
#>    <chr>    <date>     <chr>       <int> <lgl>      
#>  1 Clinical 2009-08-03 001213          1 TRUE       
#>  2 ICU      2002-01-07 067927          1 TRUE       
#>  3 ICU      2014-07-17 078381          1 TRUE       
#>  4 ICU      2015-12-03 0DBF93          1 TRUE       
#>  5 Clinical 2005-06-16 105248          1 TRUE       
#>  6 ICU      2016-05-08 110918          1 TRUE       
#>  7 ICU      2002-04-08 130252          1 TRUE       
#>  8 Clinical 2010-01-17 16F0F7          1 TRUE       
#>  9 Clinical 2008-10-19 189795          1 TRUE       
#> 10 ICU      2008-07-11 207752          1 TRUE       
#> # ℹ 90 more rows

if (require("dplyr")) {
  df %>%
    group_by(ward) %>%
    summarise(
      n_patients = n_distinct(patient),
      n_episodes_365 = sum(is_new_episode(date, episode_days = 365)),
      n_episodes_60 = sum(is_new_episode(date, episode_days = 60)),
      n_episodes_30 = sum(is_new_episode(date, episode_days = 30))
    )
}
#> # A tibble: 3 × 5
#>   ward       n_patients n_episodes_365 n_episodes_60 n_episodes_30
#>   <chr>           <int>          <int>         <int>         <int>
#> 1 Clinical           57             13            36            42
#> 2 ICU                32             11            24            29
#> 3 Outpatient          6              4             7             7

# grouping on patients and microorganisms leads to the same
# results as first_isolate() when using 'episode-based':
if (require("dplyr")) {
  x <- df %>%
    filter_first_isolate(
      include_unknown = TRUE,
      method = "episode-based"
    )

  y <- df %>%
    group_by(patient, mo) %>%
    filter(is_new_episode(date, 365)) %>%
    ungroup()

  identical(x, y)
}
#> [1] FALSE

# but is_new_episode() has a lot more flexibility than first_isolate(),
# since you can now group on anything that seems relevant:
if (require("dplyr")) {
  df %>%
    group_by(patient, mo, ward) %>%
    mutate(flag_episode = is_new_episode(date, 365)) %>%
    select(group_vars(.), flag_episode)
}
#> # A tibble: 100 × 4
#> # Groups:   patient, mo, ward [97]
#>    patient mo           ward       flag_episode
#>    <chr>   <mo>         <chr>      <lgl>       
#>  1 C82046  B_PANTO      ICU        TRUE        
#>  2 C80762  B_ESCHR_COLI ICU        TRUE        
#>  3 A81782  B_STPHY_EPDR Clinical   TRUE        
#>  4 534816  F_CANDD_ALBC Clinical   TRUE        
#>  5 A3C336  B_KLBSL_OXYT Clinical   TRUE        
#>  6 BF4515  B_ESCHR_COLI ICU        TRUE        
#>  7 CD4354  B_ESCHR_COLI ICU        TRUE        
#>  8 C25784  B_STRPT_PNMN Outpatient TRUE        
#>  9 874171  B_STRPT_PNMN Clinical   TRUE        
#> 10 067927  B_STPHY_EPDR ICU        TRUE        
#> # ℹ 90 more rows
# }
```
