# Compute all four psAve selection criteria for a propensity score vector

Evaluates, for an arbitrary propensity score vector, the four selection
criteria used by
[`psave()`](https://kabajiro.github.io/psAve/reference/psave.md): the
treatment-assignment log loss, the mean weighted absolute standardized
mean difference (ASMD) of the covariates, the mean weighted
Kolmogorov-Smirnov (KS) statistic of the covariates, and the weighted
ASMD of a prognostic score. This function powers the `diagnostics` table
of a `psave` object (the "was averaging worth it?" comparison) and is
exported as a methods-research utility.

## Usage

``` r
psave_criteria(
  ps,
  treat,
  covs,
  prog = NULL,
  estimand = c("ATT", "ATE"),
  s.d.denom = "treated",
  bin.vars = NULL
)
```

## Arguments

- ps:

  Numeric vector of propensity scores, strictly inside (0, 1).

- treat:

  Treatment vector; coerced to 0/1 like the left-hand side of the
  [`psave()`](https://kabajiro.github.io/psAve/reference/psave.md)
  formula (numeric 0/1, logical, or two-level factor/character with the
  second level treated).

- covs:

  Numeric matrix (or all-numeric data frame) of covariates, one row per
  unit. Factors must already be expanded to dummy columns (as in the
  `covs` component of a `psave` object).

- prog:

  Optional numeric vector: a prognostic score. If `NULL`, the `prog`
  criterion is returned as `NA`.

- estimand:

  `"ATT"` (default) or `"ATE"`; determines the weights (see Details).

- s.d.denom:

  Group whose (unweighted) standard deviation standardizes the mean
  differences: `"treated"` (default; the paper's convention for BOTH
  estimands), `"control"`, `"pooled"`, or `"all"`.

- bin.vars:

  Optional logical vector flagging binary columns of `covs`, used for
  the KS criterion only; if `NULL`, columns with exactly two distinct
  values are detected automatically. The `smd` and `prog` criteria
  always use uniform sample-SD standardization (`bin.vars = FALSE` for
  all columns), which is the convention of the published method (see
  Details).

## Value

A named numeric vector with elements `logloss`, `smd`, `ks`, and `prog`
(the last is `NA` when `prog = NULL`).

## Details

Weights are the inverse-probability weights implied by `ps` at
`estimand`: for the ATT, \\W_i = 1\\ for treated units and
\\e_i/(1-e_i)\\ for untreated units; for the ATE, \\1/e_i\\ and
\\1/(1-e_i)\\. The four criteria are:

- `logloss`:

  \\-\mathrm{mean}\\A_i \log e_i + (1-A_i)\log(1-e_i)\\\\.

- `smd`:

  the mean over covariates \\j\\ of \\\|\bar X\_{1j}^w - \bar
  X\_{0j}^w\| / s_j\\, where \\\bar X\_{aj}^w\\ is the weighted mean of
  \\X_j\\ in arm \\a\\ and \\s_j\\ is the unweighted sample SD of
  \\X_j\\ in the `s.d.denom` group. Computed via
  [`cobalt::col_w_smd()`](https://ngreifer.github.io/cobalt/reference/balance-summary.html).

- `ks`:

  the mean over covariates of the proper weighted-eCDF KS statistic
  \\\sup_x \|F^w_1(x) - F^w_0(x)\|\\, computed via
  [`cobalt::col_w_ks()`](https://ngreifer.github.io/cobalt/reference/balance-summary.html);
  for binary columns this is the absolute difference in weighted
  proportions.

- `prog`:

  the same weighted ASMD formula applied to the single column `prog`.

Faithful to the published method (and its reference implementation), the
`smd` and `prog` criteria standardize **every** column, including binary
ones, by the plain unweighted sample SD
([`sd()`](https://rdrr.io/r/stats/sd.html), the \\n-1\\ formula) of the
`s.d.denom` group – i.e., `bin.vars = FALSE` is passed to
[`cobalt::col_w_smd()`](https://ngreifer.github.io/cobalt/reference/balance-summary.html)
for all columns. cobalt's own display convention (binary columns
standardized by \\\sqrt{p(1-p)}\\) is used only in the display-oriented
`balance` component of a `psave` object and in
[`bal.tab.psave()`](https://kabajiro.github.io/psAve/reference/bal.tab.psave.md).
For the KS criterion the two conventions coincide.

## References

Kabata D, Stuart EA, Shintani A (2024). Prognostic score-based model
averaging approach for propensity score estimation. *BMC Medical
Research Methodology*, 24, 228.
[doi:10.1186/s12874-024-02350-y](https://doi.org/10.1186/s12874-024-02350-y)

Hansen BB (2008). The prognostic analogue of the propensity score.
*Biometrika*, 95(2), 481-488.
[doi:10.1093/biomet/asn004](https://doi.org/10.1093/biomet/asn004)

Xie Y, Zhu Y, Cotton CA, Wu P (2019). A model averaging approach for
estimating propensity scores by optimizing balance. *Statistical Methods
in Medical Research*, 28(1), 84-101.
[doi:10.1177/0962280217715487](https://doi.org/10.1177/0962280217715487)

## See also

[`psave()`](https://kabajiro.github.io/psAve/reference/psave.md),
[`simplex_grid()`](https://kabajiro.github.io/psAve/reference/simplex_grid.md),
[`cobalt::col_w_smd()`](https://ngreifer.github.io/cobalt/reference/balance-summary.html),
[`cobalt::col_w_ks()`](https://ngreifer.github.io/cobalt/reference/balance-summary.html)

## Examples

``` r
set.seed(1)
n <- 200
x1 <- rnorm(n); x2 <- rbinom(n, 1, 0.4)
a <- rbinom(n, 1, plogis(-0.5 + x1 + 0.5 * x2))
ps <- pmin(pmax(fitted(glm(a ~ x1 + x2, family = binomial())), 0.01), 0.99)
g <- 1 + 0.5 * x1 - 0.2 * x2   # a (toy) prognostic score
psave_criteria(ps, a, cbind(x1 = x1, x2 = x2), prog = g, estimand = "ATT")
#>   logloss       smd        ks      prog 
#> 0.5806282 0.1094402 0.1084529 0.2073763 
```
