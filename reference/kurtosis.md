# Kurtosis of the Sample

Kurtosis is a measure of the "tailedness" of the probability
distribution of a real-valued random variable. A normal distribution has
a kurtosis of 3 and a excess kurtosis of 0.

## Usage

``` r
kurtosis(x, na.rm = FALSE, excess = FALSE)

# Default S3 method
kurtosis(x, na.rm = FALSE, excess = FALSE)

# S3 method for class 'matrix'
kurtosis(x, na.rm = FALSE, excess = FALSE)

# S3 method for class 'data.frame'
kurtosis(x, na.rm = FALSE, excess = FALSE)
```

## Arguments

- x:

  A vector of values, a [matrix](https://rdrr.io/r/base/matrix.html) or
  a [data.frame](https://rdrr.io/r/base/data.frame.html).

- na.rm:

  A [logical](https://rdrr.io/r/base/logical.html) to indicate whether
  `NA` values should be stripped before the computation proceeds.

- excess:

  A [logical](https://rdrr.io/r/base/logical.html) to indicate whether
  the *excess kurtosis* should be returned, defined as the kurtosis
  minus 3.

## See also

[`skewness()`](https://amr-for-r.org/reference/skewness.md)

## Examples

``` r
kurtosis(rnorm(10000))
#> [1] 3.071712
kurtosis(rnorm(10000), excess = TRUE)
#> [1] -0.02774835
```
