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
#>   [1] 21  7  4 46 15  6 11 27 44 25 20  8 30 25 43 33 13 46  2 44 40 15 33 46  2
#>  [26] 32  4 19 22 45 23 45  9 44 39 16  7 46 10 34 37 45  4 20 11 14  5 19  7 40
#>  [51]  6 47 31  4 31 42  7 13 29 38 31 18 28  9 43 17 23 41 12  8 35 13  2 45  4
#>  [76]  3 34 36 15 23 24 15  9 41 12  1 44 17 18 37 32 23 35 38 44 37 38 10 24 26
is_new_episode(df$date, episode_days = 60) # TRUE/FALSE
#>   [1]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
#>  [13]  TRUE FALSE  TRUE  TRUE  TRUE FALSE  TRUE FALSE  TRUE FALSE FALSE FALSE
#>  [25] FALSE  TRUE FALSE  TRUE  TRUE  TRUE  TRUE FALSE  TRUE FALSE  TRUE  TRUE
#>  [37] FALSE FALSE  TRUE  TRUE  TRUE FALSE FALSE FALSE FALSE  TRUE  TRUE FALSE
#>  [49] FALSE FALSE FALSE  TRUE  TRUE FALSE FALSE  TRUE FALSE FALSE  TRUE  TRUE
#>  [61] FALSE  TRUE  TRUE FALSE FALSE  TRUE FALSE  TRUE  TRUE FALSE  TRUE FALSE
#>  [73] FALSE FALSE FALSE  TRUE FALSE  TRUE FALSE FALSE  TRUE FALSE FALSE FALSE
#>  [85] FALSE  TRUE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
#>  [97] FALSE FALSE FALSE  TRUE

# filter on results from the third 60-day episode only, using base R
df[which(get_episode(df$date, 60) == 3), ]
#> # A tibble: 1 × 46
#>   date       patient   age gender ward  mo           PEN   OXA   FLC   AMX  
#>   <date>     <chr>   <dbl> <chr>  <chr> <mo>         <sir> <sir> <sir> <sir>
#> 1 2003-01-06 894506     83 M      ICU   B_STRPT_PNMN   S     NA    NA    S  
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
#> # Groups:   patient, condition [94]
#>    patient date       condition new_episode
#>    <chr>   <date>     <chr>     <lgl>      
#>  1 060287  2007-03-11 A         TRUE       
#>  2 0E2483  2007-11-10 A         TRUE       
#>  3 101305  2006-12-13 A         TRUE       
#>  4 141061  2014-10-22 A         TRUE       
#>  5 151041  2006-02-10 A         TRUE       
#>  6 15D386  2004-08-01 B         TRUE       
#>  7 187841  2008-04-22 A         TRUE       
#>  8 189363  2004-06-22 C         TRUE       
#>  9 208305  2015-09-14 A         TRUE       
#> 10 257844  2011-05-22 A         TRUE       
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
#> # Groups:   ward, patient [89]
#>    ward     date       patient new_index new_logical
#>    <chr>    <date>     <chr>       <int> <lgl>      
#>  1 Clinical 2007-03-11 060287          1 TRUE       
#>  2 Clinical 2007-11-10 0E2483          1 TRUE       
#>  3 Clinical 2006-12-13 101305          1 TRUE       
#>  4 Clinical 2014-10-22 141061          1 TRUE       
#>  5 Clinical 2006-02-10 151041          1 TRUE       
#>  6 ICU      2004-08-01 15D386          1 TRUE       
#>  7 Clinical 2008-04-22 187841          1 TRUE       
#>  8 Clinical 2004-06-22 189363          1 TRUE       
#>  9 Clinical 2015-09-14 208305          1 TRUE       
#> 10 Clinical 2011-05-22 257844          1 TRUE       
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
#> 1 Clinical           55             12            34            39
#> 2 ICU                28             10            25            26
#> 3 Outpatient          6              5             6             6

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
#> # Groups:   patient, mo, ward [93]
#>    patient mo            ward       flag_episode
#>    <chr>   <mo>          <chr>      <lgl>       
#>  1 0E2483  B_ESCHR_COLI  Clinical   TRUE        
#>  2 F24801  B_STRPT_AGLC  ICU        TRUE        
#>  3 A26548  B_STPHY_CONS  ICU        TRUE        
#>  4 422833  B_ENTRC_FCLS  Clinical   TRUE        
#>  5 AFD3B1    UNKNOWN     Clinical   TRUE        
#>  6 F35553  B_ENTRBC_CLOC ICU        TRUE        
#>  7 316893  B_STPHY_CONS  Clinical   TRUE        
#>  8 B54813  B_ENTRC       ICU        TRUE        
#>  9 E06844  B_STPHY_AURS  Outpatient TRUE        
#> 10 B65162  B_STRPT_PNMN  Clinical   TRUE        
#> # ℹ 90 more rows
# }
```
