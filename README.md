# psAve

<!-- badges: start -->
[![R-CMD-check](https://github.com/kabajiro/psAve/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/kabajiro/psAve/actions/workflows/R-CMD-check.yaml)
<!-- [![CRAN status](https://www.r-pkg.org/badges/version/psAve)](https://CRAN.R-project.org/package=psAve) -->
<!-- badges: end -->

**Model-averaged propensity scores selected by prognostic-score balance.**

`psAve` implements the propensity score model averaging method of Kabata, Stuart & Shintani (2024, *BMC Medical Research Methodology*). Instead of committing to a single propensity score model, `psave()` fits several candidate models (logistic regression, CART, random forest, gradient boosting, or any `SuperLearner` wrapper), and combines them as a convex mixture whose weights are chosen on a simplex grid to optimize balance on the **prognostic score** — the predicted outcome under the untreated condition, estimated from untreated units only (Hansen 2008). Because the prognostic score summarizes the covariates *as they relate to the outcome*, balancing it targets exactly the covariate directions that drive bias in the treatment effect; in the paper's simulations this "Prog (Ave)" strategy achieved the lowest and most robust bias and RMSE across 16 scenarios. Covariate-balance (SMD, KS) and prediction-accuracy (log loss) criteria from the paper are also available.

`psAve` is a companion to the [MatchIt](https://kosukeimai.github.io/MatchIt/)/[WeightIt](https://ngreifer.github.io/WeightIt/)/[cobalt](https://ngreifer.github.io/cobalt/) ecosystem, not a replacement for any part of it. The deliverable is a plain numeric vector of propensity scores, handed to `MatchIt::matchit()` as a distance measure or to `WeightIt::weightit()` as a propensity score (thin `psave_match()` / `psave_weight()` wrappers do this without retyping the formula); balance assessment — including prognostic-score balance — works out of the box via `cobalt::bal.tab()`; effect estimation stays where it belongs, in `MatchIt`/`WeightIt`/`survey`/`marginaleffects`. Estimands: ATT (default) and ATE.

## Installation

```r
# install.packages("remotes")
remotes::install_github("kabajiro/psAve")
```

## Quick example

```r
library(psAve)
data("lalonde", package = "MatchIt")

set.seed(1234)
fit <- psave(treat ~ age + educ + race + married + nodegree + re74 + re75,
             data = lalonde, outcome = ~ re78)      # defaults: criterion = "prog", ATT
m <- psave_match(fit, method = "nearest")           # a genuine matchit object
cobalt::bal.tab(m, distance = data.frame(prog = fit$prog))
```

## Learn more

* [Getting Started with psAve](vignettes/psAve.Rmd) — motivation, the `lalonde` matching workflow, interpreting output, and why using the outcome in the design stage does not bias the analysis.
* [Reproducing the published IPW workflow](vignettes/weighting.Rmd) — ATT and ATE weighting, and the paper's `survey::svyglm()` estimator.
* [Method details and design decisions](vignettes/method-details.Rmd) — all formulas, the documented fixes relative to the paper's reference code, relation to other software, and limitations.

## Citation

If you use `psAve`, please cite the paper:

> Kabata, D., Stuart, E. A., & Shintani, A. (2024). Prognostic score-based model averaging approach for propensity score estimation. *BMC Medical Research Methodology*, 24, 228. doi:[10.1186/s12874-024-02350-y](https://doi.org/10.1186/s12874-024-02350-y)

```bibtex
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

GPL (>= 2)
