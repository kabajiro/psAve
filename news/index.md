# Changelog

## psAve 1.0.0

Initial release, implementing Kabata, Stuart & Shintani (2024), *BMC
Medical Research Methodology* 24:228, <doi:10.1186/s12874-024-02350-y>.

- [`psave()`](https://kabajiro.github.io/psAve/reference/psave.md):
  model-averaged propensity scores as a convex combination of candidate
  models (`"glm"`, `"rpart"`, `"ranger"`, `"xgboost"` by default; any
  `"SL.*"` SuperLearner wrapper; or user-supplied
  `ps.matrix`/`prog.matrix`), with mixing weights selected on a simplex
  grid.
- `ps.append` / `prog.append`: extra user-supplied candidate score
  columns (a vector of length n, or a matrix/all-numeric data frame with
  unique column names) appended AFTER the candidates from
  `ps.methods`/`ps.matrix` and `prog.methods`/`prog.matrix`. Appended
  propensity columns are validated (strictly in (0, 1)) and clipped like
  every other candidate; grid tie-breaking favors the base candidates.
  Supplying `prog.matrix` or `prog.append` without `outcome` is an
  explicit error (prognostic candidates require the outcome; gamma is
  selected by outcome-prediction MSE among untreated units).
- Four selection criteria from the paper: `"prog"` (weighted ASMD of the
  model-averaged prognostic score, the recommended default;
  per-candidate targets via `prog.target`), `"smd"`, `"ks"`, and
  `"logloss"`. Estimands: ATT (default) and ATE, with the supplement’s
  estimand-specific weight formulas.
- Prognostic scores per Hansen (2008): candidate outcome models fit on
  untreated units only; mixing weights `gamma` selected by unweighted
  untreated-set MSE. [`gaussian()`](https://rdrr.io/r/stats/family.html)
  and [`binomial()`](https://rdrr.io/r/stats/family.html) outcome
  families.
- `average = FALSE` vertex mode selects the single best candidate
  propensity score by the chosen criterion.
- Ecosystem integration: `fit$ps` drops into
  `MatchIt::matchit(distance = )` and `WeightIt::weightit(ps = )`;
  [`psave_match()`](https://kabajiro.github.io/psAve/reference/psave_match.md)
  /
  [`psave_weight()`](https://kabajiro.github.io/psAve/reference/psave_weight.md)
  wrappers reuse the stored formula and data to eliminate
  row-misalignment;
  [`cobalt::bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.html)
  works directly on `psave` objects.
- Full S3 suite: [`print()`](https://rdrr.io/r/base/print.html) (with
  the literal next call),
  [`summary()`](https://rdrr.io/r/base/summary.html) (mixing weights,
  all-criteria diagnostics table, full balance table),
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) (`"balance"`,
  `"distribution"`, `"criterion"`),
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
  [`weights()`](https://rdrr.io/r/stats/weights.html),
  [`predict()`](https://rdrr.io/r/stats/predict.html) (with
  `keep.fits = TRUE`).
- Exported utilities:
  [`simplex_grid()`](https://kabajiro.github.io/psAve/reference/simplex_grid.md)
  (integer-composition simplex enumeration defining the tie-breaking
  order) and
  [`psave_criteria()`](https://kabajiro.github.io/psAve/reference/psave_criteria.md)
  (all four criteria for any propensity score vector).
- Documented fixes relative to the paper’s reference code:
  integer-arithmetic simplex grid (the reference’s floating-point
  `rowSums == 1` filter silently dropped ~10.6% of grid points); proper
  weighted-eCDF KS statistic;
  [`binomial()`](https://rdrr.io/r/stats/family.html) family for binary
  responses; no train/test-inconsistent
  [`scale()`](https://rdrr.io/r/base/scale.html); strict complete-case
  handling (`NA`s error, never dropped). See
  [`vignette("method-details", package = "psAve")`](https://kabajiro.github.io/psAve/articles/method-details.md).
- Three vignettes: Getting Started (matching workflow), Reproducing the
  published IPW workflow (ATT/ATE weighting with
  [`survey::svyglm()`](https://rdrr.io/pkg/survey/man/svyglm.html)), and
  Method details.
