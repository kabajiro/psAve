# Balance tables for psave objects

A method for
[`cobalt::bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.html):
assesses balance on the covariates of a
[`psave()`](https://kabajiro.github.io/psAve/reference/psave.md) fit
under the implied inverse-probability weights, with the model-averaged
propensity score and (when available) the model-averaged prognostic
score supplied as `distance` measures – the prognostic-score balance
diagnostic of Stuart, Lee and Leacy (2013).

## Usage

``` r
# S3 method for class 'psave'
bal.tab(x, ...)
```

## Arguments

- x:

  A `psave` object.

- ...:

  Further arguments passed on to
  [`cobalt::bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.html)
  (e.g., `un = TRUE`, `thresholds = c(m = 0.1)`).

## Value

A `bal.tab` object; see
[`cobalt::bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.html).

## Details

The call delegates to the default cobalt machinery as
`cobalt::bal.tab(<covariates>, treat = x$treat, weights = x$weights, s.d.denom = x$s.d.denom, distance = data.frame(ps = x$ps, prog = x$prog), ...)`,
so all the usual cobalt arguments (`un`, `stats`, `thresholds`, ...) are
available, and display conventions are cobalt's own (the *selection
criterion* inside
[`psave()`](https://kabajiro.github.io/psAve/reference/psave.md) uses
the paper's uniform sample-SD standardization instead; see
[`psave_criteria()`](https://kabajiro.github.io/psAve/reference/psave_criteria.md)).

## References

Stuart EA, Lee BK, Leacy FP (2013). Prognostic score-based balance
measures can be a useful diagnostic for propensity score methods in
comparative effectiveness research. *Journal of Clinical Epidemiology*,
66(8), S84-S90.
[doi:10.1016/j.jclinepi.2013.01.013](https://doi.org/10.1016/j.jclinepi.2013.01.013)

## See also

[`psave()`](https://kabajiro.github.io/psAve/reference/psave.md),
[`cobalt::bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.html),
[`plot.psave()`](https://kabajiro.github.io/psAve/reference/plot.psave.md)

## Examples

``` r
data("lalonde", package = "MatchIt")
fit <- psave(treat ~ age + educ + married + re74, data = lalonde,
             outcome = ~ re78, ps.methods = "glm", prog.methods = "glm")
cobalt::bal.tab(fit, un = TRUE)
#> Balance Measures
#>             Type Diff.Un Diff.Adj
#> ps      Distance  0.9447   0.0394
#> prog    Distance -0.6949  -0.0340
#> age      Contin. -0.3094   0.0402
#> educ     Contin.  0.0550  -0.0220
#> married   Binary -0.3236   0.0008
#> re74     Contin. -0.7211  -0.0271
#> 
#> Effective sample sizes
#>            Control Treated
#> Unadjusted  429.       185
#> Adjusted    278.11     185
```
