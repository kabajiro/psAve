# Method details: relation to the paper's reference implementation

psAve implements the published method of Kabata, Stuart and Shintani
(2024) rather than reproducing its reference code line by line. Five
documented defects/quirks of the reference implementation are
deliberately **fixed** (each fix follows the paper's *stated*
definitions):

1.  **Integer simplex grid.** The reference code enumerated the
    mixing-weight grid with
    [`expand.grid()`](https://rdrr.io/r/base/expand.grid.html) and kept
    rows passing an exact floating-point `rowSums(gr) == 1` test,
    silently dropping about 10.6% of the valid grid points for M = 4
    candidates at step 0.05 (1,584 of 1,771 kept).
    [`simplex_grid()`](https://kabajiro.github.io/psAve/reference/simplex_grid.md)
    enumerates integer compositions, so every valid point is present by
    construction, and the enumeration order makes the first-minimum
    tie-break a reproducible rule.

2.  **Proper weighted-eCDF KS statistic.** The reference `Fks` computed
    [`ks.test()`](https://rdrr.io/r/stats/ks.test.html) on covariate
    values *multiplied by* the weights, which is not the paper's
    weighted-eCDF definition. psAve computes \\\sup_x \|F^w_1(x) -
    F^w_0(x)\|\\ with weighted empirical CDFs in each arm (as
    [`cobalt::col_w_ks()`](https://ngreifer.github.io/cobalt/reference/balance-summary.html)
    does).

3.  **Proper binomial family.** The reference code fit binary-response
    SuperLearner models with `gaussian(link = "logit")`; psAve uses
    [`binomial()`](https://rdrr.io/r/stats/family.html) throughout for
    treatment models.

4.  **No per-set [`scale()`](https://rdrr.io/r/base/scale.html).** The
    reference code standardized fitting and prediction sets
    *separately*, an inconsistent transformation. psAve passes raw
    covariates to all engines.

5.  **Strict complete-case alignment.** The reference `Fasmd` applied
    [`na.omit()`](https://rdrr.io/r/stats/na.fail.html) to a covariate
    while using full-length treatment and weight vectors, silently
    misaligning rows in the presence of missing data. psAve refuses
    missing values in any used variable (error, never drop), and names
    all returned score vectors by `rownames(data)`.

## Details

### Clipping

Candidate propensity scores are clipped to `clip` (default \[0.01,
0.99\], the paper's constants) **before** averaging. The average is
never re-clipped: a convex combination of values inside the clipping
interval cannot leave it.

### Criterion vs. display conventions

The `smd`/`prog` *selection criteria* standardize all columns (including
binary ones) by the plain unweighted sample SD of the `s.d.denom` group
– the paper's convention (its reference `Fasmd` uses
[`sd()`](https://rdrr.io/r/stats/sd.html)), implemented by passing
`bin.vars = FALSE` for every column to
[`cobalt::col_w_smd()`](https://ngreifer.github.io/cobalt/reference/balance-summary.html).
The display-oriented `balance` component and
[`bal.tab.psave()`](https://kabajiro.github.io/psAve/reference/bal.tab.psave.md)
follow cobalt's native conventions (binary columns use
\\\sqrt{p(1-p)}\\). For `criterion = "prog"` the denominator is a
positive constant across the \\\lambda\\ grid, so the selected
\\\lambda\\ is invariant to this choice; the *reported*
`criterion.value` uses the paper's convention. For the KS criterion the
two conventions coincide on binary columns.

### Relation to other software

The closest existing functionality is
[WeightIt::method_super](https://ngreifer.github.io/WeightIt/reference/method_super.html)
with `SL.method = "method.balance"` (Pirracchio and Carone 2018): a
covariate-balance-targeted SuperLearner for *weighting only*. psAve
differs in targeting **prognostic-score** balance (Hansen 2008; Stuart,
Lee and Leacy 2013), in searching an exhaustive simplex grid with a
documented tie-break instead of a convex-optimization meta-learner, and
in returning a plain score vector equally usable for **matching**
(`MatchIt::matchit(distance = )`) and weighting. No other package
implements propensity score model averaging in the lineage of Xie et al.
(2019).

### Reproducibility

There is no `seed` argument (setting the global RNG inside a function is
an R anti-pattern): call
[`set.seed()`](https://rdrr.io/r/base/Random.html) before
[`psave()`](https://kabajiro.github.io/psAve/reference/psave.md) when
stochastic learners (`ranger`, `xgboost`) are among the candidates. The
resolved hyperparameters and engine package versions are recorded in
`fit$info$learners`.

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

Pirracchio R, Carone M (2018). The Balance Super Learner: A robust
adaptation of the Super Learner to improve estimation of the average
treatment effect in the treated based on propensity score matching.
*Statistical Methods in Medical Research*, 27(8), 2504-2518.
[doi:10.1177/0962280216682055](https://doi.org/10.1177/0962280216682055)

## See also

[`psave()`](https://kabajiro.github.io/psAve/reference/psave.md),
[`simplex_grid()`](https://kabajiro.github.io/psAve/reference/simplex_grid.md),
[`psave_criteria()`](https://kabajiro.github.io/psAve/reference/psave_criteria.md),
[`vignette("method-details", "psAve")`](https://kabajiro.github.io/psAve/articles/method-details.md)
