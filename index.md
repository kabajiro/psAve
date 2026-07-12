# psAve

[![CRAN
status](https://www.r-pkg.org/badges/version/psAve)](https://CRAN.R-project.org/package=psAve)
[![CRAN
downloads](https://cranlogs.r-pkg.org/badges/psAve)](https://CRAN.R-project.org/package=psAve)
[![R-CMD-check](https://github.com/kabajiro/psAve/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/kabajiro/psAve/actions/workflows/R-CMD-check.yaml)

**Model-averaged propensity scores selected by prognostic-score
balance.**

`psAve` implements the propensity score model averaging method of
Kabata, Stuart & Shintani (2024, *BMC Medical Research Methodology*).
Instead of committing to a single propensity score model,
[`psave()`](https://kabajiro.github.io/psAve/reference/psave.md) fits
several candidate models (logistic regression, CART, random forest,
gradient boosting, or any `SuperLearner` wrapper), and combines them as
a convex mixture whose weights are chosen on a simplex grid to optimize
balance on the **prognostic score** — the predicted outcome under the
untreated condition, estimated from untreated units only (Hansen 2008).
Because the prognostic score summarizes the covariates *as they relate
to the outcome*, balancing it targets exactly the covariate directions
that drive bias in the treatment effect; in the paper’s simulations this
“Prog (Ave)” strategy achieved the lowest and most robust bias and RMSE
across 16 scenarios. Covariate-balance (SMD, KS) and prediction-accuracy
(log loss) criteria from the paper are also available.

`psAve` is a companion to the
[MatchIt](https://kosukeimai.github.io/MatchIt/)/[WeightIt](https://ngreifer.github.io/WeightIt/)/[cobalt](https://ngreifer.github.io/cobalt/)
ecosystem, not a replacement for any part of it. The deliverable is a
plain numeric vector of propensity scores, handed to
[`MatchIt::matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.html)
as a distance measure or to
[`WeightIt::weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.html)
as a propensity score (thin
[`psave_match()`](https://kabajiro.github.io/psAve/reference/psave_match.md)
/
[`psave_weight()`](https://kabajiro.github.io/psAve/reference/psave_weight.md)
wrappers do this without retyping the formula); balance assessment —
including prognostic-score balance — works out of the box via
[`cobalt::bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.html);
effect estimation stays where it belongs, in
`MatchIt`/`WeightIt`/`survey`/`marginaleffects`. Estimands: ATT
(default) and ATE.

## Installation

``` r

# From CRAN
install.packages("psAve")

# Development version
# install.packages("remotes")
remotes::install_github("kabajiro/psAve")
```

## Quick example

``` r

library(psAve)
data("lalonde", package = "MatchIt")

set.seed(1234)
fit <- psave(treat ~ age + educ + race + married + nodegree + re74 + re75,
             data = lalonde, outcome = ~ re78)      # defaults: criterion = "prog", ATT
m <- psave_match(fit, method = "nearest")           # a genuine matchit object
cobalt::bal.tab(m, distance = data.frame(prog = fit$prog))
```

## Learn more

- [Getting Started with
  psAve](https://kabajiro.github.io/psAve/articles/psAve.html) —
  motivation, the `lalonde` matching workflow, interpreting output, and
  why using the outcome in the design stage does not bias the analysis.
- [Reproducing the published IPW
  workflow](https://kabajiro.github.io/psAve/articles/weighting.html) —
  ATT and ATE weighting, and the paper’s
  [`survey::svyglm()`](https://rdrr.io/pkg/survey/man/svyglm.html)
  estimator.
- [Method details and design
  decisions](https://kabajiro.github.io/psAve/articles/method-details.html)
  — all formulas, the documented fixes relative to the paper’s reference
  code, relation to other software, and limitations.
- Full documentation: <https://kabajiro.github.io/psAve/> ・ CRAN:
  <https://CRAN.R-project.org/package=psAve>

## Citation

If you use `psAve`, please cite the paper:

> Kabata, D., Stuart, E. A., & Shintani, A. (2024). Prognostic
> score-based model averaging approach for propensity score estimation.
> *BMC Medical Research Methodology*, 24, 228.
> <doi:%5B10.1186/s12874-024-02350-y>\](<https://doi.org/10.1186/s12874-024-02350-y>)

``` bibtex
@article{kabata2024prognostic,
  author  = {Kabata, Daijiro and Stuart, Elizabeth A. and Shintani, Ayumi},
  title   = {Prognostic score-based model averaging approach for propensity score estimation},
  journal = {BMC Medical Research Methodology},
  year    = {2024},
  volume  = {24},
  pages   = {228},
  doi     = {10.1186/s12874-024-02350-y}
}
```

## License

GPL (\>= 2)
