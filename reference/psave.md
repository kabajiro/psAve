# Model-averaged propensity scores selected by prognostic-score balance

`psave()` constructs a model-averaged propensity score \\\bar e(X) =
\sum_m \lambda_m \hat e_m(X)\\: a convex combination of candidate
propensity score models whose mixing weights \\\lambda\\ are selected on
a simplex grid to optimize a balance criterion – by default the weighted
absolute standardized mean difference of a **model-averaged prognostic
score** (the "Prog (Ave)" estimator of Kabata, Stuart and Shintani
2024). The result is deliberately modest: a numeric score vector
designed to be handed to
[`MatchIt::matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.html)
as `distance`, or to
[`WeightIt::weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.html)
as `ps`, with balance assessment via cobalt.

## Usage

``` r
psave(
  formula,
  data,
  outcome = NULL,
  estimand = c("ATT", "ATE"),
  criterion = c("prog", "smd", "ks", "logloss"),
  prog.target = "average",
  ps.methods = c("glm", "rpart", "ranger", "xgboost"),
  prog.methods = c("glm", "rpart", "ranger", "xgboost"),
  ps.matrix = NULL,
  prog.matrix = NULL,
  ps.append = NULL,
  prog.append = NULL,
  average = TRUE,
  family = gaussian(),
  step = 0.05,
  clip = c(0.01, 0.99),
  s.d.denom = "treated",
  cv = 5L,
  control = list(),
  keep.fits = FALSE,
  keep.path = TRUE,
  verbose = FALSE,
  ...
)
```

## Arguments

- formula:

  A two-sided formula `treat ~ x1 + x2 + ...`, exactly as in
  [`MatchIt::matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.html).
  The right-hand side defines both the candidate-PS covariates and the
  balance covariates used by the `smd`/`ks` criteria.

- data:

  A data frame containing the variables in `formula` (and `outcome`).
  Complete cases in all used variables are REQUIRED; any missing value
  is an error, never a silent row drop.

- outcome:

  The outcome specification for the prognostic score: a one-sided
  formula `~ y` (the formula right-hand side is reused as the prognostic
  predictors) or a two-sided formula `y ~ z1 + z2` (a distinct
  prognostic specification). Required when `criterion = "prog"`;
  optional for the outcome-free criteria (`"smd"`, `"ks"`, `"logloss"`),
  where it may still be supplied so that the prognostic score is
  estimated and reported in `diagnostics`, `balance`, and
  [`bal.tab.psave()`](https://kabajiro.github.io/psAve/reference/bal.tab.psave.md).
  Prognostic models are fit on **untreated units only** (Hansen 2008),
  so using the outcome here does not bias effect estimation; see
  Details.

- estimand:

  `"ATT"` (default) or `"ATE"`; determines the inverse-probability
  weights used inside the balance criteria and returned in `weights`.

- criterion:

  The selection criterion for \\\lambda\\: `"prog"` (default) = weighted
  ASMD of the prognostic score (the paper's headline "Prog (Ave)");
  `"smd"` = mean weighted ASMD over the covariates; `"ks"` = mean
  weighted Kolmogorov-Smirnov statistic over the covariates; `"logloss"`
  = negative Bernoulli log-likelihood of treatment assignment (the
  prediction-accuracy criterion in the lineage of Xie et al. 2019).

- prog.target:

  Only used when `criterion = "prog"`: `"average"` (default) targets the
  gamma-mixed prognostic score \\\bar g\\; naming a single element of
  `prog.methods` (or a column of `prog.matrix`) targets that candidate
  prognostic score instead (the paper's "Prog (\\g_k\\)" variants).

- ps.methods:

  Character vector of candidate propensity score learners. The fixed,
  explicit default is `c("glm", "rpart", "ranger", "xgboost")` (no
  auto-detection: a machine-dependent candidate set is irreproducible
  science). Any `"SL.*"` string is accepted verbatim and passed through
  to SuperLearner. The order defines the tie-break preference (see
  Details).

- prog.methods:

  Character vector of candidate learners for the untreated-only
  prognostic models \\g_k(0, X)\\; same menu and default as
  `ps.methods`.

- ps.matrix:

  Optional n x M numeric matrix of user-supplied candidate propensity
  scores (values strictly in (0, 1), column names required). Overrides
  `ps.methods`; the columns are clipped like fitted candidates.

- prog.matrix:

  Optional n x K numeric matrix of user-supplied candidate prognostic
  scores (column names required). Overrides `prog.methods`. Requires
  `outcome` (gamma is selected by outcome-prediction MSE among untreated
  units).

- ps.append:

  Optional extra candidate propensity scores appended AFTER the
  candidates from `ps.methods` or `ps.matrix` (default `NULL`): a
  numeric vector of length n (one candidate, labeled `"append"`) or a
  numeric matrix / all-numeric data frame with n rows (unique, non-empty
  column names required, as for `ps.matrix`). Values must lie strictly
  in (0, 1) and each column is clipped to `clip` before averaging,
  exactly like every other candidate. Because appended columns come
  last, grid tie-breaking (first minimum) favors the base candidates
  (see Details).

- prog.append:

  Optional extra candidate prognostic scores appended AFTER the
  candidates from `prog.methods` or `prog.matrix` (default `NULL`); same
  vector/matrix/data-frame forms and naming rules as `ps.append`, but
  the values are unrestricted finite reals. Requires `outcome`.

- average:

  If `FALSE`, the \\\lambda\\ grid is restricted to the simplex
  VERTICES, i.e., the best single PS candidate by `criterion` is
  selected (the "best single learner" variants of the paper's
  supplement).

- family:

  The prognostic-model family:
  [`gaussian()`](https://rdrr.io/r/stats/family.html) (default) or
  [`binomial()`](https://rdrr.io/r/stats/family.html) only. With
  [`binomial()`](https://rdrr.io/r/stats/family.html), the
  gamma-selection MSE is the Brier score – the formula is unchanged;
  note that the paper's simulations validated continuous outcomes.

- step:

  The simplex-grid increment for BOTH \\\gamma\\ and \\\lambda\\
  (default `0.05`, the paper's value). Must evenly divide 1.

- clip:

  Length-2 numeric: each candidate propensity score column is clipped to
  `[clip[1], clip[2]]` BEFORE averaging (default `c(0.01, 0.99)`, the
  paper's constants). The average is never re-clipped: a convex
  combination of values in the clipping interval stays in it.

- s.d.denom:

  The group whose unweighted standard deviation standardizes mean
  differences in the ASMD-based criteria, passed to
  [`cobalt::col_w_smd()`](https://ngreifer.github.io/cobalt/reference/balance-summary.html):
  `"treated"` (default; the paper's supplement uses the unweighted
  TREATED-group SD for both the ATT and the ATE), `"control"`,
  `"pooled"`, or `"all"`.

- cv:

  Number of cross-validation folds `V` for the
  [`SuperLearner::SuperLearner()`](https://rdrr.io/pkg/SuperLearner/man/SuperLearner.html)
  `cvControl` when `"SL.*"` learners are used (default `5`).

- control:

  A named list of per-learner hyperparameter overrides, e.g.
  `list(ranger = list(num.trees = 1000), xgboost = list(nrounds = 200))`.
  Entries for `"glm"` are passed to
  [`stats::glm()`](https://rdrr.io/r/stats/glm.html), for `"rpart"` to
  [`rpart::rpart.control()`](https://rdrr.io/pkg/rpart/man/rpart.control.html),
  for `"ranger"` to
  [`ranger::ranger()`](http://imbs-hl.github.io/ranger/reference/ranger.md),
  and for `"xgboost"` to the `params` list of
  [`xgboost::xgb.train()`](https://rdrr.io/pkg/xgboost/man/xgb.train.html)
  (with `nrounds` and `verbose` recognized as top-level arguments). The
  resolved values are stored in `info$learners`. Multi-threaded engines
  run single-threaded by default (`ranger` `num.threads = 1`, `xgboost`
  `nthread = 1`) in line with CRAN's at-most-2-cores policy; raise these
  via `control` to speed up real analyses on your own machine.

- keep.fits:

  If `TRUE`, the fitted learner objects are retained in `fits`, enabling
  [`predict.psave()`](https://kabajiro.github.io/psAve/reference/predict.psave.md).
  Default `FALSE`.

- keep.path:

  If `TRUE` (default), the full \\\lambda\\ criterion path is stored in
  `path`. Automatically set to `NULL` with a message if the grid exceeds
  100,000 rows.

- verbose:

  If `TRUE`, progress messages report the learner set, the grid sizes,
  the selected \\\lambda\\/\\\gamma\\, and the criterion value.

- ...:

  Reserved for future use; supplying unused arguments triggers a
  warning.

## Value

An object of class `"psave"`: a list with components

- `ps`:

  numeric(n), named by `rownames(data)`: the model-averaged propensity
  score \\\bar e\\ – the deliverable for `MatchIt::matchit(distance = )`
  / `WeightIt::weightit(ps = )`.

- `prog`:

  numeric(n), named: the model-averaged prognostic score \\\bar g\\;
  `NULL` when `outcome` was not supplied (possible only for the
  outcome-free criteria `"smd"`, `"ks"`, and `"logloss"`).

- `lambda`:

  named numeric(M): the selected PS mixing weights, where M counts the
  base candidates plus any `ps.append` columns.

- `gamma`:

  named numeric(K) (or `NULL`): the selected prognostic mixing weights,
  where K counts any `prog.append` columns.

- `weights`:

  numeric(n): the IPW at `estimand` implied by `ps`.

- `ps.candidates`:

  n x M matrix of clipped candidate propensity scores (fitted or
  user-supplied; `ps.append` columns come last).

- `prog.candidates`:

  n x K matrix of candidate prognostic predictions for all n units (or
  `NULL`); `prog.append` columns come last.

- `criterion`, `criterion.value`:

  the criterion used and its value at the selected \\\lambda\\
  (paper-faithful standardization).

- `diagnostics`:

  (M+1) x 4 data frame: all four criteria
  ([`psave_criteria()`](https://kabajiro.github.io/psAve/reference/psave_criteria.md))
  for each candidate and for the selected average – the "was averaging
  worth it?" table. Its `prog` column always uses the model-averaged
  prognostic score \\\bar g\\, even when `prog.target` names a single
  candidate (see Details); it then need not equal `criterion.value`.

- `path`:

  data frame of the full \\\lambda\\ grid (M columns + `value` + logical
  `selected`) in grid order, or `NULL`.

- `gamma.mse`:

  named numeric(K+1): untreated MSE per prognostic candidate and for the
  selected average (or `NULL`).

- `balance`:

  data frame (covariates + `prog`): unweighted and weighted SMD and KS
  (cobalt display conventions).

- `treat`:

  integer(n) 0/1 treatment as used.

- `covs`:

  numeric n x p balance-covariate matrix with `attr(, "bin.vars")`.

- `estimand`, `s.d.denom`, `prog.target`, `average`:

  scalars, as resolved.

- `outcome.name`:

  name of the outcome variable, or `NA`.

- `formula`, `data`:

  as supplied (they power
  [`psave_match()`](https://kabajiro.github.io/psAve/reference/psave_match.md),
  [`psave_weight()`](https://kabajiro.github.io/psAve/reference/psave_weight.md)
  and
  [`print.psave()`](https://kabajiro.github.io/psAve/reference/print.psave.md);
  note the memory cost of storing `data`).

- `fits`:

  list of fitted learners iff `keep.fits = TRUE` (enables
  [`predict.psave()`](https://kabajiro.github.io/psAve/reference/predict.psave.md));
  otherwise `NULL`.

- `info`:

  list: `step`, `clip`, `cv`, `family`, `grid.size`, `n`, `learners`
  (labels, resolved hyperparameters, engine versions), `psAve.version`.

- `call`:

  the matched call.

## Details

### Algorithm

1.  **Candidate propensity scores.** Each learner in `ps.methods` is fit
    on all \\n\\ units and predicts \\\hat e_m(X_i) = P(A_i = 1 \mid
    X_i)\\ in-sample; each column is clipped to `clip` before averaging.
    Extra user-supplied candidates given via `ps.append` are validated,
    clipped identically, and appended AFTER the base candidates (from
    `ps.methods` or `ps.matrix`), so the first-minimum tie-break favors
    the base set.

2.  **Model-averaged prognostic score.** Each learner in `prog.methods`
    is fit on the **untreated units only** and predicts \\\hat g_k(0,
    X_i)\\ for all \\n\\ units; `prog.append` columns are appended after
    these candidates. The mixing weights \\\gamma\\ minimize the
    unweighted untreated-set mean squared error \$\$\mathrm{MSE}(\gamma)
    = \mathrm{mean}\_{i:A_i=0}\bigl(Y_i - \textstyle\sum_k \gamma_k \hat
    g_k(0, X_i)\bigr)^2\$\$ over `simplex_grid(K, step)`; \\\bar g =
    \sum_k \gamma_k \hat g_k\\.

3.  **Model-averaged propensity score.** The mixing weights \\\lambda\\
    minimize `criterion` over `simplex_grid(M, step)`, where each grid
    row implies \\\bar e\_\lambda = \sum_m \lambda_m \hat e_m\\ and the
    weights \\W_i\\ below.

### Inverse-probability weights (exact supplement formulas)

With \\\bar e_i\\ the averaged propensity score: \$\$\mathrm{ATT:}\quad
W_i = 1 \\ (A_i = 1), \qquad W_i = \bar e_i / (1 - \bar e_i) \\ (A_i =
0);\$\$ \$\$\mathrm{ATE:}\quad W_i = 1 / \bar e_i \\ (A_i = 1), \qquad
W_i = 1 / (1 - \bar e_i) \\ (A_i = 0).\$\$ These are identical to
[`WeightIt::get_w_from_ps()`](https://ngreifer.github.io/WeightIt/reference/get_w_from_ps.html)
at the same estimand.

### Selection criteria

At each grid row \\\lambda\\:

- `logloss`:

  \\-\mathrm{mean}\\A_i \log \bar e_i + (1 - A_i) \log(1 - \bar e_i)\\\\
  (finite by clipping).

- `smd`:

  the mean over covariate columns \\j\\ of \$\$\mathrm{ASMD}\_j =
  \Bigl\| \frac{\sum\_{A_i=1} W_i X\_{ij}}{\sum\_{A_i=1} W_i} -
  \frac{\sum\_{A_i=0} W_i X\_{ij}}{\sum\_{A_i=0} W_i} \Bigr\| / s_j,\$\$
  where \\s_j\\ is the **unweighted sample SD** of \\X_j\\ in the
  `s.d.denom` group (the treated group for both estimands, per the
  paper's supplement).

- `ks`:

  the mean over covariates of the proper weighted-eCDF
  Kolmogorov-Smirnov statistic \\\sup_x \|F^w_1(x) - F^w_0(x)\|\\ with
  \\F^w_a(x) = \sum\_{i:A_i=a} W_i 1(X\_{ij} \le x) / \sum\_{i:A_i=a}
  W_i\\; for binary columns this equals the absolute difference in
  weighted proportions.

- `prog`:

  the \\\mathrm{ASMD}\\ formula applied to the single column \\\bar g\\
  (`prog.target = "average"`) or \\\hat g_k\\ (`prog.target` names a
  learner). THE DEFAULT.

When `prog.target` names a single prognostic candidate,
`criterion.value` is the weighted ASMD of that candidate's prognostic
score, whereas the `prog` column of the `diagnostics` table always
refers to the model-averaged prognostic score \\\bar g\\; the two values
therefore need not coincide, and
[`summary.psave()`](https://kabajiro.github.io/psAve/reference/summary.psave.md)
prints a footnote to this effect. Faithful to the published method, the
`smd` and `prog` criteria standardize **all** columns – including binary
ones – by the plain unweighted sample SD (uniform sample-SD
standardization; `bin.vars = FALSE` is passed to
[`cobalt::col_w_smd()`](https://ngreifer.github.io/cobalt/reference/balance-summary.html)
for every column). The display-oriented `balance` component and
[`bal.tab.psave()`](https://kabajiro.github.io/psAve/reference/bal.tab.psave.md)
use cobalt's native conventions instead; see
[`vignette("method-details", "psAve")`](https://kabajiro.github.io/psAve/articles/method-details.md).

### Simplex grid and tie-breaking

Both grids are enumerated by
[`simplex_grid()`](https://kabajiro.github.io/psAve/reference/simplex_grid.md)
in integer arithmetic (every valid grid point is present by
construction) in a documented order: the first component descends from 1
to 0, then the second on the remainder, and so on – the first row puts
all weight on the first candidate. Ties in any argmin take the FIRST row
attaining the minimum (within a 1e-9 relative numerical tolerance), so
ties favor learners listed earlier in `ps.methods`/`prog.methods`;
candidates appended via `ps.append`/`prog.append` always come last and
therefore lose ties to the base candidates. The tolerance is deliberate:
the criterion values come from floating-point matrix algebra whose
lowest-order bits can differ across BLAS implementations, so an exact
bitwise [`which.min()`](https://rdrr.io/r/base/which.min.html) would not
be reproducible across machines, whereas the tolerant first-minimum rule
is. `gamma` always minimizes the unweighted untreated-set MSE; these
rules are fixed and not arguments.

### Why using the outcome does not bias the design

Prognostic models see the outcomes of **untreated units only**, and the
criterion compares weighted *covariate-like summaries* (the prognostic
score) between arms – it never uses a treated-untreated outcome contrast
(Hansen 2008; Stuart, Lee and Leacy 2013). This is the same argument
under which prognostic-score balance diagnostics are recommended for
propensity score analyses; see
[`vignette("method-details", "psAve")`](https://kabajiro.github.io/psAve/articles/method-details.md).

### Reproducibility

There is no `seed` argument: call
[`set.seed()`](https://rdrr.io/r/base/Random.html) before `psave()` –
stochastic learners are `ranger` and `xgboost`. `info$learners` records
the resolved hyperparameters and engine package versions.

### Relation to other software

[WeightIt::method_super](https://ngreifer.github.io/WeightIt/reference/method_super.html)
with `SL.method = "method.balance"` (Pirracchio and Carone 2018) selects
a SuperLearner combination by *covariate* balance and is weighting-only.
`psave()` targets *prognostic-score* balance on an exhaustive simplex
grid, and returns a score vector equally usable for matching and
weighting. See
[psave-details](https://kabajiro.github.io/psAve/reference/psave-details.md)
and
[`vignette("method-details", "psAve")`](https://kabajiro.github.io/psAve/articles/method-details.md)
for the differences from the paper's reference implementation.

## References

Kabata D, Stuart EA, Shintani A (2024). Prognostic score-based model
averaging approach for propensity score estimation. *BMC Medical
Research Methodology*, 24, 228.
[doi:10.1186/s12874-024-02350-y](https://doi.org/10.1186/s12874-024-02350-y)

Hansen BB (2008). The prognostic analogue of the propensity score.
*Biometrika*, 95(2), 481-488.
[doi:10.1093/biomet/asn004](https://doi.org/10.1093/biomet/asn004)

Stuart EA, Lee BK, Leacy FP (2013). Prognostic score-based balance
measures can be a useful diagnostic for propensity score methods in
comparative effectiveness research. *Journal of Clinical Epidemiology*,
66(8), S84-S90.
[doi:10.1016/j.jclinepi.2013.01.013](https://doi.org/10.1016/j.jclinepi.2013.01.013)

Xie Y, Zhu Y, Cotton CA, Wu P (2019). A model averaging approach for
estimating propensity scores by optimizing balance. *Statistical Methods
in Medical Research*, 28(1), 84-101.
[doi:10.1177/0962280217715487](https://doi.org/10.1177/0962280217715487)

## See also

[`psave_match()`](https://kabajiro.github.io/psAve/reference/psave_match.md),
[`psave_weight()`](https://kabajiro.github.io/psAve/reference/psave_weight.md),
[`simplex_grid()`](https://kabajiro.github.io/psAve/reference/simplex_grid.md),
[`psave_criteria()`](https://kabajiro.github.io/psAve/reference/psave_criteria.md),
[`bal.tab.psave()`](https://kabajiro.github.io/psAve/reference/bal.tab.psave.md),
[`plot.psave()`](https://kabajiro.github.io/psAve/reference/plot.psave.md),
[`summary.psave()`](https://kabajiro.github.io/psAve/reference/summary.psave.md),
[`predict.psave()`](https://kabajiro.github.io/psAve/reference/predict.psave.md)

## Examples

``` r
if (requireNamespace("MatchIt", quietly = TRUE)) {
  data("lalonde", package = "MatchIt")

  # Paper-headline "Prog (Ave)", ATT. For speed, this example restricts the
  # candidate learners; the default is
  # ps.methods = prog.methods = c("glm", "rpart", "ranger", "xgboost").
  meths <- if (requireNamespace("rpart", quietly = TRUE)) c("glm", "rpart") else "glm"
  fit <- psave(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, outcome = ~ re78,
               ps.methods = meths, prog.methods = meths)
  fit
  summary(fit)

  # hand off to MatchIt (canonical explicit call):
  m <- MatchIt::matchit(treat ~ age + educ + race + married + nodegree +
                          re74 + re75,
                        data = lalonde, distance = fit$ps)
  # or, reusing the stored formula and data (no retyping hazard):
  m2 <- psave_match(fit)
}
```
