# Transform Input to Minimum Inhibitory Concentrations (MIC)

This transforms vectors to a new class `mic`, which treats the input as
decimal numbers, while maintaining operators (such as "\>=") and only
allowing valid MIC values known to the field of (medical) microbiology.

## Usage

``` r
as.mic(x, na.rm = FALSE)

NA_mic_

is.mic(x)

# S3 method for class 'mic'
droplevels(x, as.mic = FALSE, ...)
```

## Arguments

- x:

  a [character](https://rdrr.io/r/base/character.html) or
  [numeric](https://rdrr.io/r/base/numeric.html) vector

- na.rm:

  a [logical](https://rdrr.io/r/base/logical.html) indicating whether
  missing values should be removed

- as.mic:

  a [logical](https://rdrr.io/r/base/logical.html) to indicate whether
  the `mic` class should be kept - the default is `FALSE`

- ...:

  arguments passed on to methods

## Value

Ordered [factor](https://rdrr.io/pkg/data.table/man/fctr.html) with
additional class `mic`, that in mathematical operations acts as decimal
numbers. Bare in mind that the outcome of any mathematical operation on
MICs will return a [numeric](https://rdrr.io/r/base/numeric.html) value.

## Details

To interpret MIC values as SIR values, use
[`as.sir()`](https://amr-for-r.org/release/2.0.0/reference/as.sir.md) on
MIC values. It supports guidelines from EUCAST (2013-2022) and CLSI
(2013-2022).

This class for MIC values is a quite a special data type: formally it is
an ordered [factor](https://rdrr.io/pkg/data.table/man/fctr.html) with
valid MIC values as
[factor](https://rdrr.io/pkg/data.table/man/fctr.html) levels (to make
sure only valid MIC values are retained), but for any mathematical
operation it acts as decimal numbers:

    x <- random_mic(10)
    x
    #> Class 'mic'
    #>  [1] 16     1      8      8      64     >=128  0.0625 32     32     16

    is.factor(x)
    #> [1] TRUE

    x[1] * 2
    #> [1] 32

    median(x)
    #> [1] 26

This makes it possible to maintain operators that often come with MIC
values, such "\>=" and "\<=", even when filtering using
[numeric](https://rdrr.io/r/base/numeric.html) values in data analysis,
e.g.:

    x[x > 4]
    #> Class 'mic'
    #> [1] 16    8     8     64    >=128 32    32    16

    df <- data.frame(x, hospital = "A")
    subset(df, x > 4) # or with dplyr: df %>% filter(x > 4)
    #>        x hospital
    #> 1     16        A
    #> 5     64        A
    #> 6  >=128        A
    #> 8     32        A
    #> 9     32        A
    #> 10    16        A

The following [generic
functions](https://rdrr.io/r/base/groupGeneric.html) are implemented for
the MIC class: `!`, `!=`, `%%`, `%/%`, `&`, `*`, `+`, `-`, `/`, `<`,
`<=`, `==`, `>`, `>=`, `^`, `|`,
[`abs()`](https://rdrr.io/r/base/MathFun.html),
[`acos()`](https://rdrr.io/r/base/Trig.html),
[`acosh()`](https://rdrr.io/r/base/Hyperbolic.html),
[`all()`](https://rdrr.io/r/base/all.html),
[`any()`](https://rdrr.io/r/base/any.html),
[`asin()`](https://rdrr.io/r/base/Trig.html),
[`asinh()`](https://rdrr.io/r/base/Hyperbolic.html),
[`atan()`](https://rdrr.io/r/base/Trig.html),
[`atanh()`](https://rdrr.io/r/base/Hyperbolic.html),
[`ceiling()`](https://rdrr.io/r/base/Round.html),
[`cos()`](https://rdrr.io/r/base/Trig.html),
[`cosh()`](https://rdrr.io/r/base/Hyperbolic.html),
[`cospi()`](https://rdrr.io/r/base/Trig.html),
[`cummax()`](https://rdrr.io/r/base/cumsum.html),
[`cummin()`](https://rdrr.io/r/base/cumsum.html),
[`cumprod()`](https://rdrr.io/r/base/cumsum.html),
[`cumsum()`](https://rdrr.io/r/base/cumsum.html),
[`digamma()`](https://rdrr.io/r/base/Special.html),
[`exp()`](https://rdrr.io/r/base/Log.html),
[`expm1()`](https://rdrr.io/r/base/Log.html),
[`floor()`](https://rdrr.io/r/base/Round.html),
[`gamma()`](https://rdrr.io/r/base/Special.html),
[`lgamma()`](https://rdrr.io/r/base/Special.html),
[`log()`](https://rdrr.io/r/base/Log.html),
[`log1p()`](https://rdrr.io/r/base/Log.html),
[`log2()`](https://rdrr.io/r/base/Log.html),
[`log10()`](https://rdrr.io/r/base/Log.html),
[`max()`](https://rdrr.io/r/base/Extremes.html),
[`mean()`](https://rdrr.io/r/base/mean.html),
[`min()`](https://rdrr.io/r/base/Extremes.html),
[`prod()`](https://rdrr.io/r/base/prod.html),
[`range()`](https://rdrr.io/r/base/range.html),
[`round()`](https://rdrr.io/r/base/Round.html),
[`sign()`](https://rdrr.io/r/base/sign.html),
[`signif()`](https://rdrr.io/r/base/Round.html),
[`sin()`](https://rdrr.io/r/base/Trig.html),
[`sinh()`](https://rdrr.io/r/base/Hyperbolic.html),
[`sinpi()`](https://rdrr.io/r/base/Trig.html),
[`sqrt()`](https://rdrr.io/r/base/MathFun.html),
[`sum()`](https://rdrr.io/r/base/sum.html),
[`tan()`](https://rdrr.io/r/base/Trig.html),
[`tanh()`](https://rdrr.io/r/base/Hyperbolic.html),
[`tanpi()`](https://rdrr.io/r/base/Trig.html),
[`trigamma()`](https://rdrr.io/r/base/Special.html) and
[`trunc()`](https://rdrr.io/r/base/Round.html). Some functions of the
`stats` package are also implemented:
[`median()`](https://rdrr.io/r/stats/median.html),
[`quantile()`](https://rdrr.io/r/stats/quantile.html),
[`mad()`](https://rdrr.io/r/stats/mad.html),
[`IQR()`](https://rdrr.io/r/stats/IQR.html),
[`fivenum()`](https://rdrr.io/r/stats/fivenum.html). Also,
[`boxplot.stats()`](https://rdrr.io/r/grDevices/boxplot.stats.html) is
supported. Since [`sd()`](https://rdrr.io/r/stats/sd.html) and
[`var()`](https://rdrr.io/r/stats/cor.html) are non-generic functions,
these could not be extended. Use
[`mad()`](https://rdrr.io/r/stats/mad.html) as an alternative, or use
e.g. `sd(as.numeric(x))` where `x` is your vector of MIC values.

Using [`as.double()`](https://rdrr.io/r/base/double.html) or
[`as.numeric()`](https://rdrr.io/r/base/numeric.html) on MIC values will
remove the operators and return a numeric vector. Do **not** use
[`as.integer()`](https://rdrr.io/r/base/integer.html) on MIC values as
by the R convention on
[factor](https://rdrr.io/pkg/data.table/man/fctr.html)s, it will return
the index of the factor levels (which is often useless for regular
users).

Use
[`droplevels()`](https://rdrr.io/pkg/data.table/man/fdroplevels.html) to
drop unused levels. At default, it will return a plain factor. Use
`droplevels(..., as.mic = TRUE)` to maintain the `mic` class.

`NA_mic_` is a missing value of the new `mic` class, analogous to e.g.
base R's [`NA_character_`](https://rdrr.io/r/base/NA.html).

## See also

[`as.sir()`](https://amr-for-r.org/release/2.0.0/reference/as.sir.md)

## Examples

``` r
mic_data <- as.mic(c(">=32", "1.0", "1", "1.00", 8, "<=0.128", "8", "16", "16"))
mic_data
is.mic(mic_data)

# this can also coerce combined MIC/SIR values:
as.mic("<=0.002; S")

# mathematical processing treats MICs as numeric values
fivenum(mic_data)
quantile(mic_data)
all(mic_data < 512)

# interpret MIC values
as.sir(
  x = as.mic(2),
  mo = as.mo("Streptococcus pneumoniae"),
  ab = "AMX",
  guideline = "EUCAST"
)
as.sir(
  x = as.mic(c(0.01, 2, 4, 8)),
  mo = as.mo("Streptococcus pneumoniae"),
  ab = "AMX",
  guideline = "EUCAST"
)

# plot MIC values, see ?plot
plot(mic_data)
plot(mic_data, mo = "E. coli", ab = "cipro")

if (require("ggplot2")) {
  autoplot(mic_data, mo = "E. coli", ab = "cipro")
}
if (require("ggplot2")) {
  autoplot(mic_data, mo = "E. coli", ab = "cipro", language = "nl") # Dutch
}
if (require("ggplot2")) {
  autoplot(mic_data, mo = "E. coli", ab = "cipro", language = "uk") # Ukrainian
}
```
