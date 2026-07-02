# Predict averaged scores for new data

Computes the model-averaged propensity score (or prognostic score) for
new observations by applying the stored candidate fits to `newdata`,
clipping candidate propensity scores as at fit time, and combining them
with the selected mixing weights. Requires
`psave(..., keep.fits = TRUE)`.

## Usage

``` r
# S3 method for class 'psave'
predict(object, newdata, type = c("ps", "prog"), ...)
```

## Arguments

- object:

  A `psave` object fitted with `keep.fits = TRUE`.

- newdata:

  A data frame containing the variables of the propensity score formula
  (for `type = "ps"`) or of the prognostic specification (for
  `type = "prog"`). Missing values are an error.

- type:

  `"ps"` (default) or `"prog"`.

- ...:

  Ignored.

## Value

A numeric vector with one score per row of `newdata`, named by its
rownames. If `newdata` is missing, the in-sample fitted scores are
returned (equivalent to
[`fitted.psave()`](https://kabajiro.github.io/psAve/reference/fitted.psave.md)).

## See also

[`psave()`](https://kabajiro.github.io/psAve/reference/psave.md),
[`fitted.psave()`](https://kabajiro.github.io/psAve/reference/fitted.psave.md)
