# Weight by a model-averaged propensity score

Convenience pass-through to
[`WeightIt::weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.html):
constructs balancing weights from the model-averaged propensity score of
a [`psave()`](https://kabajiro.github.io/psAve/reference/psave.md) fit
at the fitted estimand, reusing the stored formula and data. Equivalent
to the canonical explicit call

    WeightIt::weightit(<formula>, data = <data>, ps = fit$ps, estimand = fit$estimand, ...)

All `...` arguments are forwarded verbatim; the return value is an
ordinary `weightit` object.

## Usage

``` r
psave_weight(x, ...)
```

## Arguments

- x:

  A `psave` object.

- ...:

  Arguments forwarded verbatim to
  [`WeightIt::weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.html).

## Value

A `weightit` object; see
[`WeightIt::weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.html).

## See also

[`psave()`](https://kabajiro.github.io/psAve/reference/psave.md),
[`psave_match()`](https://kabajiro.github.io/psAve/reference/psave_match.md),
[`WeightIt::weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.html),
[`WeightIt::get_w_from_ps()`](https://ngreifer.github.io/WeightIt/reference/get_w_from_ps.html)

## Examples

``` r
data("lalonde", package = "MatchIt")
fit <- psave(treat ~ age + educ + married + re74, data = lalonde,
             outcome = ~ re78, ps.methods = "glm", prog.methods = "glm")
w <- psave_weight(fit)
cobalt::bal.tab(w, distance = data.frame(prog = fit$prog))
#> Balance Measures
#>                Type Diff.Adj
#> prog       Distance  -0.0340
#> prop.score Distance   0.0394
#> age         Contin.   0.0402
#> educ        Contin.  -0.0220
#> married      Binary   0.0008
#> re74        Contin.  -0.0271
#> 
#> Effective sample sizes
#>            Control Treated
#> Unadjusted  429.       185
#> Adjusted    278.11     185
```
