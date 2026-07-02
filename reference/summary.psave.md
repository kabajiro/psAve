# Summarize a psave object

Produces (a) the selected mixing-weight tables \\\lambda\\ and
\\\gamma\\, (b) the `diagnostics` table (all four selection criteria for
every candidate propensity score and for the selected average – the "was
averaging worth it?" comparison), and (c) the full balance table (all
covariates plus the prognostic score; unweighted vs. weighted SMD and
KS, with a `*` marker at weighted SMD \> 0.1).

## Usage

``` r
# S3 method for class 'psave'
summary(object, un = TRUE, candidates = TRUE, ...)

# S3 method for class 'summary.psave'
print(x, digits = 3, ...)
```

## Arguments

- object:

  A `psave` object.

- un:

  If `TRUE` (default), the balance table includes the unweighted
  columns.

- candidates:

  If `TRUE` (default), the per-candidate `diagnostics` table is
  included.

- ...:

  Ignored.

- x:

  A `summary.psave` object.

- digits:

  Number of significant digits to print. Default 3.

## Value

For `summary.psave()`, an object of class `"summary.psave"`: a list with
elements `lambda`, `gamma`, `gamma.mse`, `diagnostics`, `balance`,
`criterion`, `criterion.value`, `prog.target`, `estimand`, `average`,
`nn`, and `call`. `print.summary.psave()` returns `x` invisibly.

## See also

[`psave()`](https://kabajiro.github.io/psAve/reference/psave.md),
[`print.psave()`](https://kabajiro.github.io/psAve/reference/print.psave.md),
[`bal.tab.psave()`](https://kabajiro.github.io/psAve/reference/bal.tab.psave.md)
