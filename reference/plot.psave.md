# Plot a psave object

Three diagnostic displays for a
[`psave()`](https://kabajiro.github.io/psAve/reference/psave.md) fit:

- `"balance"`:

  a Love plot of covariate and prognostic-score balance before/after
  weighting, via
  [`cobalt::love.plot()`](https://ngreifer.github.io/cobalt/reference/love.plot.html)
  (dispatched through
  [`bal.tab.psave()`](https://kabajiro.github.io/psAve/reference/bal.tab.psave.md);
  cobalt is an Import, so always available).

- `"distribution"`:

  the distribution of the propensity scores by treatment group: grey
  curves for the candidate models, a colored curve for the selected
  average (base graphics). Inspect this for extreme candidate scores.

- `"criterion"`:

  the selection criterion over the \\\lambda\\ grid: exact for up to
  three candidates (a curve for M = 2, a colored grid map for M = 3);
  for M \> 3, one profile per candidate (the minimum criterion value
  attainable at each value of \\\lambda_m\\). Requires the stored `path`
  (`keep.path = TRUE` and a grid of at most 100,000 rows).

## Usage

``` r
# S3 method for class 'psave'
plot(x, type = c("balance", "distribution", "criterion"), ...)
```

## Arguments

- x:

  A `psave` object.

- type:

  One of `"balance"` (default), `"distribution"`, `"criterion"`.

- ...:

  For `"balance"`, further arguments to
  [`cobalt::love.plot()`](https://ngreifer.github.io/cobalt/reference/love.plot.html)
  (e.g., `thresholds = 0.1`); otherwise further graphical parameters
  passed to the base plotting calls.

## Value

For `"balance"`, the `ggplot` object from
[`cobalt::love.plot()`](https://ngreifer.github.io/cobalt/reference/love.plot.html)
(invisibly, after printing); otherwise `x`, invisibly.

## See also

[`psave()`](https://kabajiro.github.io/psAve/reference/psave.md),
[`bal.tab.psave()`](https://kabajiro.github.io/psAve/reference/bal.tab.psave.md),
[`cobalt::love.plot()`](https://ngreifer.github.io/cobalt/reference/love.plot.html)
