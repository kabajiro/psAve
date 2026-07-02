# Print a psave object

Prints a one-screen summary of a fitted
[`psave()`](https://kabajiro.github.io/psAve/reference/psave.md) object:
estimand and criterion, the selected mixing weights \\\lambda\\ and
\\\gamma\\ as labeled text bars, the criterion value, a three-row
balance preview (the worst covariates plus the prognostic score), and
then the **literal next call** – echoing the formula and data name from
your own
[`psave()`](https://kabajiro.github.io/psAve/reference/psave.md) call –
that hands the averaged score to
[`MatchIt::matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.html)
or
[`WeightIt::weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.html).

## Usage

``` r
# S3 method for class 'psave'
print(x, digits = 3, ...)
```

## Arguments

- x:

  A `psave` object.

- digits:

  Number of significant digits to print. Default 3.

- ...:

  Ignored.

## Value

`x`, invisibly.

## See also

[`psave()`](https://kabajiro.github.io/psAve/reference/psave.md),
[`summary.psave()`](https://kabajiro.github.io/psAve/reference/summary.psave.md)
