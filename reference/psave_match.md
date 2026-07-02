# Match on a model-averaged propensity score

Convenience pass-through to
[`MatchIt::matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.html):
matches on the model-averaged propensity score of a
[`psave()`](https://kabajiro.github.io/psAve/reference/psave.md) fit,
reusing the formula and data stored in the object. Equivalent to the
canonical explicit call

    MatchIt::matchit(<formula>, data = <data>, distance = fit$ps, ...)

but with no opportunity for row misalignment between the two steps. All
`...` arguments are forwarded verbatim; the return value is an ordinary
`matchit` object, so the full MatchIt/cobalt toolkit applies.

## Usage

``` r
psave_match(x, ...)
```

## Arguments

- x:

  A `psave` object.

- ...:

  Arguments forwarded verbatim to
  [`MatchIt::matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.html)
  (e.g., `method`, `caliper`, `ratio`, `replace`).

## Value

A `matchit` object; see
[`MatchIt::matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.html).

## See also

[`psave()`](https://kabajiro.github.io/psAve/reference/psave.md),
[`psave_weight()`](https://kabajiro.github.io/psAve/reference/psave_weight.md),
[`MatchIt::matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.html)

## Examples

``` r
data("lalonde", package = "MatchIt")
fit <- psave(treat ~ age + educ + married + re74, data = lalonde,
             outcome = ~ re78, ps.methods = "glm", prog.methods = "glm")
m <- psave_match(fit, method = "nearest", caliper = 0.2)
cobalt::bal.tab(m, distance = data.frame(prog = fit$prog))
#> Balance Measures
#>              Type Diff.Adj
#> prog     Distance   0.0063
#> distance Distance   0.0752
#> age       Contin.   0.1292
#> educ      Contin.  -0.0349
#> married    Binary  -0.0324
#> re74      Contin.   0.0366
#> 
#> Sample sizes
#>           Control Treated
#> All           429     185
#> Matched       185     185
#> Unmatched     244       0
```
