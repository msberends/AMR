# Skewness of the Sample

Skewness is a measure of the asymmetry of the probability distribution
of a real-valued random variable about its mean.

When negative ('left-skewed'): the left tail is longer; the mass of the
distribution is concentrated on the right of a histogram. When positive
('right-skewed'): the right tail is longer; the mass of the distribution
is concentrated on the left of a histogram. A normal distribution has a
skewness of 0.

## Usage

``` r
skewness(x, na.rm = FALSE)

# Default S3 method
skewness(x, na.rm = FALSE)

# S3 method for class 'matrix'
skewness(x, na.rm = FALSE)

# S3 method for class 'data.frame'
skewness(x, na.rm = FALSE)
```

## Arguments

- x:

  A vector of values, a [matrix](https://rdrr.io/r/base/matrix.html) or
  a [data.frame](https://rdrr.io/r/base/data.frame.html).

- na.rm:

  A [logical](https://rdrr.io/r/base/logical.html) value indicating
  whether `NA` values should be stripped before the computation
  proceeds.

## See also

[`kurtosis()`](https://amr-for-r.org/reference/kurtosis.md)

## Examples

``` r
skewness(runif(1000))
#> [1] -0.03760694
```
