# Extract the averaged propensity or prognostic score

[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) is the
canonical extractor for the model-averaged scores of a
[`psave()`](https://kabajiro.github.io/psAve/reference/psave.md) fit;
[`weights()`](https://rdrr.io/r/stats/weights.html) extracts the
inverse-probability weights implied by the averaged propensity score at
the fitted estimand.

## Usage

``` r
# S3 method for class 'psave'
fitted(object, type = c("ps", "prog"), ...)

# S3 method for class 'psave'
weights(object, ...)
```

## Arguments

- object:

  A `psave` object.

- type:

  `"ps"` (default) for the model-averaged propensity score, or `"prog"`
  for the model-averaged prognostic score.

- ...:

  Ignored.

## Value

For [`fitted()`](https://rdrr.io/r/stats/fitted.values.html), a numeric
vector named by the rownames of the analyzed data. For
[`weights()`](https://rdrr.io/r/stats/weights.html), the numeric vector
`object$weights` (weights at the *fitted* estimand only; for other
estimands use
`WeightIt::get_w_from_ps(fitted(object), object$treat, estimand = ...)`).

## See also

[`psave()`](https://kabajiro.github.io/psAve/reference/psave.md),
[`predict.psave()`](https://kabajiro.github.io/psAve/reference/predict.psave.md)
