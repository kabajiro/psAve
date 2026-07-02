# Reproducing the published IPW workflow

## Introduction

The worked example in Kabata, Stuart & Shintani (2024) uses the
model-averaged propensity score for *inverse probability weighting*
(IPW): the averaged score $`\bar e(X)`$ is converted into weights, and
the treatment effect is estimated as a weighted regression of the
outcome on treatment using
[`survey::svyglm()`](https://rdrr.io/pkg/survey/man/svyglm.html), which
supplies design-based (sandwich) standard errors. This vignette
reproduces that workflow with `psAve`, for both the ATT (the paper’s
primary estimand) and the ATE (added in the paper’s supplement), and
shows the two equivalent routes to the weights: the ecosystem route
through `WeightIt`, and the direct route through
[`weights()`](https://rdrr.io/r/stats/weights.html).

We use the `lalonde` data from `MatchIt` throughout. See
[`vignette("psAve", package = "psAve")`](https://kabajiro.github.io/psAve/articles/psAve.md)
for an introduction to
[`psave()`](https://kabajiro.github.io/psAve/reference/psave.md) itself.

``` r

library(psAve)
data("lalonde", package = "MatchIt")
```

## ATT weighting

First, fit the model-averaged propensity score with the default settings
— criterion `"prog"` (weighted ASMD of the model-averaged prognostic
score) and `estimand = "ATT"`:

``` r

set.seed(1234)
fit <- psave(treat ~ age + educ + race + married + nodegree + re74 + re75,
             data = lalonde, outcome = ~ re78)
fit
#> A psave object (model-averaged propensity score)
#>  - estimand:  ATT
#>  - criterion: prog (weighted ASMD of the model-averaged prognostic score)
#>  - sample:    614 units (185 treated, 429 control)
#> 
#> lambda (PS mixing weights):
#>   glm      0.000  |                    |
#>   rpart    0.000  |                    |
#>   ranger   0.250  |=====               |
#>   xgboost  0.750  |===============     |
#> 
#> gamma (prognostic mixing weights):
#>   glm      0.000  |                    |
#>   rpart    0.000  |                    |
#>   ranger   0.000  |                    |
#>   xgboost  1.000  |====================|
#> 
#> Criterion value at selected lambda: 0.00506
#> 
#> Balance preview (worst covariates + prognostic score):
#>           smd.un smd.wt ks.un ks.wt
#> racewhite  1.882  0.147 0.558 0.044
#> raceblack  1.762  0.131 0.640 0.048
#> re75       0.290  0.105 0.288 0.121
#> prog       0.315  0.005 0.176 0.142
#> 
#> Next:
#>   MatchIt::matchit(treat ~ age + educ + race + married + nodegree + re74 + re75, data = lalonde, distance = x$ps)
#>     or: psave_match(x)
#>   WeightIt::weightit(treat ~ age + educ + race + married + nodegree + re74 + re75, data = lalonde, ps = x$ps, estimand = "ATT")
#>     or: psave_weight(x)
```

### The ecosystem route: `psave_weight()`

[`psave_weight()`](https://kabajiro.github.io/psAve/reference/psave_weight.md)
hands the averaged score to
[`WeightIt::weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.html)
as a fixed, externally estimated propensity score (`ps = fit$ps`),
reusing the formula and data stored in the `psave` object so no retyping
(and no row-misalignment) can occur. The result is an ordinary
`weightit` object, so the full `WeightIt`/`cobalt` toolkit applies:

``` r

w <- psave_weight(fit)   # a genuine weightit object, estimand = "ATT"
w
#> A weightit object
#>  - method: "glm" (propensity score weighting with GLM)
#>  - number of obs.: 614
#>  - sampling weights: none
#>  - treatment: 2-category
#>  - estimand: ATT (focal: 1)
#>  - covariates: age, educ, race, married, nodegree, re74, re75

cobalt::bal.tab(w, distance = data.frame(prog = fit$prog))
#> Balance Measures
#>                 Type Diff.Adj
#> prog        Distance   0.0051
#> prop.score  Distance   1.1799
#> age          Contin.  -0.0897
#> educ         Contin.  -0.0120
#> race_black    Binary   0.0475
#> race_hispan   Binary  -0.0039
#> race_white    Binary  -0.0437
#> married       Binary  -0.0063
#> nodegree      Binary   0.0339
#> re74         Contin.   0.0294
#> re75         Contin.   0.1050
#> 
#> Effective sample sizes
#>            Control Treated
#> Unadjusted  429.       185
#> Adjusted     36.58     185
```

This is equivalent to the explicit call

``` r

WeightIt::weightit(treat ~ age + educ + race + married + nodegree + re74 + re75,
                   data = lalonde, ps = fit$ps, estimand = "ATT")
```

### The direct route: `weights()`

The `psave` object already contains the IPW weights implied by the
averaged score at the fitted estimand, and `weights(fit)` extracts them.
For the ATT these are

``` math
W_i = \begin{cases} 1 & A_i = 1 \\ \dfrac{\bar e(X_i)}{1 - \bar e(X_i)} & A_i = 0, \end{cases}
```

the standard “weighting by the odds.” They are identical (and
unit-tested to be identical) to
`WeightIt::get_w_from_ps(fit$ps, treat, estimand = "ATT")`.

### The paper’s estimator: `survey::svyglm()`

The published example estimates the ATT by fitting a weighted
generalized linear model of the outcome on treatment within a survey
design, which provides robust (sandwich) standard errors:

``` r

library(survey)
#> Loading required package: grid
#> Loading required package: Matrix
#> Loading required package: survival
#> 
#> Attaching package: 'survey'
#> The following object is masked from 'package:graphics':
#> 
#>     dotchart

des.att <- svydesign(ids = ~1, weights = weights(fit), data = lalonde)
fit.att <- svyglm(re78 ~ treat, design = des.att)
summary(fit.att)
#> 
#> Call:
#> svyglm(formula = re78 ~ treat, design = des.att)
#> 
#> Survey design:
#> svydesign(ids = ~1, weights = weights(fit), data = lalonde)
#> 
#> Coefficients:
#>             Estimate Std. Error t value Pr(>|t|)    
#> (Intercept)   5949.4      702.4    8.47   <2e-16 ***
#> treat          399.8      909.2    0.44     0.66    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> (Dispersion parameter for gaussian family taken to be 52929230)
#> 
#> Number of Fisher Scoring iterations: 2
```

The coefficient on `treat` is the IPW-ATT estimate of the effect of the
treatment on `re78` among the treated.

Two honest caveats, consistent with the paper. First,
[`svyglm()`](https://rdrr.io/pkg/survey/man/svyglm.html)’s sandwich
standard errors treat the weights as *fixed*; they do not account for
the estimation of the propensity score (including the selection of the
mixing weights). This is standard practice in applied IPW analyses and
is what the paper reports; if you want uncertainty that reflects the
whole pipeline, bootstrap the entire procedure from
[`psave()`](https://kabajiro.github.io/psAve/reference/psave.md) onward.
Second, [`svyglm()`](https://rdrr.io/pkg/survey/man/svyglm.html) always
returns a `gaussian`-family fit for a continuous outcome; for binary or
other outcome types, follow the effect-measure guidance in `WeightIt`’s
and `MatchIt`’s estimation vignettes (which use the `marginaleffects`
package) rather than reading off coefficients.

## ATE weighting

For the ATE, set `estimand = "ATE"` in
[`psave()`](https://kabajiro.github.io/psAve/reference/psave.md). This
changes two things internally, both taken from the paper’s supplement:
the weights entering every balance criterion during the grid search, and
the weights returned for estimation. The prognostic model is unchanged —
it is always fit on untreated units only.

``` r

set.seed(1234)
fit.ate <- psave(treat ~ age + educ + race + married + nodegree + re74 + re75,
                 data = lalonde, outcome = ~ re78, estimand = "ATE")
fit.ate
#> A psave object (model-averaged propensity score)
#>  - estimand:  ATE
#>  - criterion: prog (weighted ASMD of the model-averaged prognostic score)
#>  - sample:    614 units (185 treated, 429 control)
#> 
#> lambda (PS mixing weights):
#>   glm      0.100  |==                  |
#>   rpart    0.850  |=================   |
#>   ranger   0.000  |                    |
#>   xgboost  0.050  |=                   |
#> 
#> gamma (prognostic mixing weights):
#>   glm      0.000  |                    |
#>   rpart    0.000  |                    |
#>   ranger   0.000  |                    |
#>   xgboost  1.000  |====================|
#> 
#> Criterion value at selected lambda: 0.00162
#> 
#> Balance preview (worst covariates + prognostic score):
#>            smd.un smd.wt ks.un ks.wt
#> racewhite   1.882  0.741 0.558 0.220
#> racehispan  0.350  0.577 0.083 0.136
#> age         0.309  0.464 0.158 0.184
#> prog        0.315  0.002 0.176 0.264
#> 
#> Next:
#>   MatchIt::matchit(treat ~ age + educ + race + married + nodegree + re74 + re75, data = lalonde, distance = x$ps)
#>     or: psave_match(x)
#>   WeightIt::weightit(treat ~ age + educ + race + married + nodegree + re74 + re75, data = lalonde, ps = x$ps, estimand = "ATE")
#>     or: psave_weight(x)
```

The ATE weights are the usual inverse-probability weights,

``` math
W_i = \begin{cases} \dfrac{1}{\bar e(X_i)} & A_i = 1 \\ \dfrac{1}{1 - \bar e(X_i)} & A_i = 0, \end{cases}
```

and the estimator is the same
[`svyglm()`](https://rdrr.io/pkg/survey/man/svyglm.html) call with the
new weights:

``` r

des.ate <- svydesign(ids = ~1, weights = weights(fit.ate), data = lalonde)
fit.ate.glm <- svyglm(re78 ~ treat, design = des.ate)
summary(fit.ate.glm)
#> 
#> Call:
#> svyglm(formula = re78 ~ treat, design = des.ate)
#> 
#> Survey design:
#> svydesign(ids = ~1, weights = weights(fit.ate), data = lalonde)
#> 
#> Coefficients:
#>             Estimate Std. Error t value Pr(>|t|)    
#> (Intercept)   6414.8      387.0  16.577   <2e-16 ***
#> treat         -356.5      823.5  -0.433    0.665    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> (Dispersion parameter for gaussian family taken to be 46458463)
#> 
#> Number of Fisher Scoring iterations: 2
```

Note that the ATE and ATT fits are *different `psave` objects*: because
the criterion is computed with estimand-specific weights, the selected
mixing weights $`\lambda`$ (and hence the averaged score itself)
generally differ between estimands. Do not reuse an ATT-selected score
for an ATE analysis.

The `weightit` route works the same way — `psave_weight(fit.ate)` passes
`estimand = "ATE"` through automatically.

## Note on estimand-specific formulas

For reference, the conventions used by the selection criteria (details
and derivations in
[`vignette("method-details", package = "psAve")`](https://kabajiro.github.io/psAve/articles/method-details.md)):

- **Weights.** ATT: $`W = 1`$ (treated), $`\bar e/(1-\bar e)`$
  (untreated). ATE: $`1/\bar e`$ (treated), $`1/(1-\bar e)`$
  (untreated). These enter both the balance criteria over the
  $`\lambda`$ grid and the returned `weights`.
- **Standardization.** The weighted absolute standardized mean
  difference (of covariates for `criterion = "smd"`, of the prognostic
  score for `criterion = "prog"`) is standardized by the *unweighted
  standard deviation in the treated group* — for **both** estimands.
  This follows the paper’s supplement exactly; it differs from the
  pooled-SD convention some software uses for the ATE.
- **KS statistic.** `criterion = "ks"` uses the proper weighted
  empirical CDF in each arm, with the estimand-specific weights above.
- **Extreme weights.** Candidate propensity scores are clipped to
  `[0.01, 0.99]` (the `clip` argument) *before* averaging, which bounds
  all weights; the convex average never leaves the clipped range, so no
  re-clipping is applied afterward. Inspect
  `plot(fit.ate, type = "distribution")` and `summary(w)` for remaining
  extreme weights, which are a greater concern for the ATE than for the
  ATT.

Only the ATT and ATE are supported; the supplement’s criterion formulas
were validated for these two estimands, and
[`psave()`](https://kabajiro.github.io/psAve/reference/psave.md) refuses
others with an error rather than guessing.

## References

Kabata, D., Stuart, E. A., & Shintani, A. (2024). Prognostic score-based
model averaging approach for propensity score estimation. *BMC Medical
Research Methodology*, 24, 228. <doi:10.1186/s12874-024-02350-y>

Lumley, T. (2004). Analysis of complex survey samples. *Journal of
Statistical Software*, 9(1), 1–19. <doi:10.18637/jss.v009.i08>
