# psAve-package.R -- package-level documentation and import directives.

#' psAve: Model-Averaged Propensity Scores Selected by Prognostic-Score Balance
#'
#' Constructs a model-averaged propensity score as a convex combination of
#' candidate propensity score models, with mixing weights selected on a
#' simplex grid to optimize covariate or prognostic-score balance,
#' implementing the method of Kabata, Stuart and Shintani (2024). Prognostic
#' scores follow Hansen (2008): outcome models are fit on untreated units
#' only. The resulting score is designed to be supplied directly to
#' [MatchIt::matchit()] as a distance measure or to [WeightIt::weightit()] as
#' a propensity score, with balance assessment via \pkg{cobalt}.
#'
#' The single estimation function is [psave()]. Its result hands off to the
#' existing ecosystem: [psave_match()] / [psave_weight()] (or the equivalent
#' explicit [MatchIt::matchit()] / [WeightIt::weightit()] calls),
#' [cobalt::bal.tab()] (which has a method for `psave` objects), and effect
#' estimation with the tools you already use (\pkg{survey},
#' \pkg{marginaleffects}, ...). Utilities [simplex_grid()] and
#' [psave_criteria()] expose the grid enumeration and the selection criteria
#' for methods research and testing. See [psave-details] for the differences
#' from the paper's reference implementation.
#'
#' @references
#' Kabata D, Stuart EA, Shintani A (2024). Prognostic score-based model
#' averaging approach for propensity score estimation. *BMC Medical Research
#' Methodology*, 24, 228. \doi{10.1186/s12874-024-02350-y}
#'
#' Hansen BB (2008). The prognostic analogue of the propensity score.
#' *Biometrika*, 95(2), 481-488. \doi{10.1093/biomet/asn004}
#'
#' Stuart EA, Lee BK, Leacy FP (2013). Prognostic score-based balance
#' measures can be a useful diagnostic for propensity score methods in
#' comparative effectiveness research. *Journal of Clinical Epidemiology*,
#' 66(8), S84-S90. \doi{10.1016/j.jclinepi.2013.01.013}
#'
#' Xie Y, Zhu Y, Cotton CA, Wu P (2019). A model averaging approach for
#' estimating propensity scores by optimizing balance. *Statistical Methods in
#' Medical Research*, 28(1), 84-101. \doi{10.1177/0962280217715487}
#'
#' @importFrom stats gaussian binomial weights fitted predict
#' @importFrom graphics plot
#' @keywords internal
"_PACKAGE"

#' Method details: relation to the paper's reference implementation
#'
#' @description
#' \pkg{psAve} implements the published method of Kabata, Stuart and Shintani
#' (2024) rather than reproducing its reference code line by line. Five
#' documented defects/quirks of the reference implementation are deliberately
#' **fixed** (each fix follows the paper's *stated* definitions):
#'
#' 1. **Integer simplex grid.** The reference code enumerated the mixing-weight
#'    grid with `expand.grid()` and kept rows passing an exact floating-point
#'    `rowSums(gr) == 1` test, silently dropping about 10.6% of the valid grid
#'    points for M = 4 candidates at step 0.05 (1,584 of 1,771 kept).
#'    [simplex_grid()] enumerates integer compositions, so every valid point
#'    is present by construction, and the enumeration order makes the
#'    first-minimum tie-break a reproducible rule.
#' 2. **Proper weighted-eCDF KS statistic.** The reference `Fks` computed
#'    `ks.test()` on covariate values *multiplied by* the weights, which is
#'    not the paper's weighted-eCDF definition. \pkg{psAve} computes
#'    \eqn{\sup_x |F^w_1(x) - F^w_0(x)|} with weighted empirical CDFs in each
#'    arm (as [cobalt::col_w_ks()] does).
#' 3. **Proper binomial family.** The reference code fit binary-response
#'    SuperLearner models with `gaussian(link = "logit")`; \pkg{psAve} uses
#'    `binomial()` throughout for treatment models.
#' 4. **No per-set `scale()`.** The reference code standardized fitting and
#'    prediction sets *separately*, an inconsistent transformation. \pkg{psAve}
#'    passes raw covariates to all engines.
#' 5. **Strict complete-case alignment.** The reference `Fasmd` applied
#'    `na.omit()` to a covariate while using full-length treatment and weight
#'    vectors, silently misaligning rows in the presence of missing data.
#'    \pkg{psAve} refuses missing values in any used variable (error, never
#'    drop), and names all returned score vectors by `rownames(data)`.
#'
#' @details
#' ## Clipping
#' Candidate propensity scores are clipped to `clip` (default
#' \[0.01, 0.99\], the paper's constants) **before** averaging. The average is
#' never re-clipped: a convex combination of values inside the clipping
#' interval cannot leave it.
#'
#' ## Criterion vs. display conventions
#' The `smd`/`prog` *selection criteria* standardize all columns (including
#' binary ones) by the plain unweighted sample SD of the `s.d.denom` group --
#' the paper's convention (its reference `Fasmd` uses `sd()`), implemented by
#' passing `bin.vars = FALSE` for every column to [cobalt::col_w_smd()]. The
#' display-oriented `balance` component and [bal.tab.psave()] follow
#' \pkg{cobalt}'s native conventions (binary columns use
#' \eqn{\sqrt{p(1-p)}}). For `criterion = "prog"` the denominator is a
#' positive constant across the \eqn{\lambda} grid, so the selected
#' \eqn{\lambda} is invariant to this choice; the *reported*
#' `criterion.value` uses the paper's convention. For the KS criterion the
#' two conventions coincide on binary columns.
#'
#' ## Relation to other software
#' The closest existing functionality is [WeightIt::method_super] with
#' `SL.method = "method.balance"` (Pirracchio and Carone 2018): a
#' covariate-balance-targeted SuperLearner for *weighting only*. \pkg{psAve}
#' differs in targeting **prognostic-score** balance (Hansen 2008; Stuart,
#' Lee and Leacy 2013), in searching an exhaustive simplex grid with a
#' documented tie-break instead of a convex-optimization meta-learner, and in
#' returning a plain score vector equally usable for **matching**
#' (`MatchIt::matchit(distance = )`) and weighting. No other package
#' implements propensity score model averaging in the lineage of Xie et al.
#' (2019).
#'
#' ## Reproducibility
#' There is no `seed` argument (setting the global RNG inside a function is an
#' R anti-pattern): call `set.seed()` before [psave()] when stochastic
#' learners (`ranger`, `xgboost`) are among the candidates. The resolved
#' hyperparameters and engine package versions are recorded in
#' `fit$info$learners`.
#'
#' @references
#' Kabata D, Stuart EA, Shintani A (2024). Prognostic score-based model
#' averaging approach for propensity score estimation. *BMC Medical Research
#' Methodology*, 24, 228. \doi{10.1186/s12874-024-02350-y}
#'
#' Hansen BB (2008). The prognostic analogue of the propensity score.
#' *Biometrika*, 95(2), 481-488. \doi{10.1093/biomet/asn004}
#'
#' Stuart EA, Lee BK, Leacy FP (2013). Prognostic score-based balance
#' measures can be a useful diagnostic for propensity score methods in
#' comparative effectiveness research. *Journal of Clinical Epidemiology*,
#' 66(8), S84-S90. \doi{10.1016/j.jclinepi.2013.01.013}
#'
#' Xie Y, Zhu Y, Cotton CA, Wu P (2019). A model averaging approach for
#' estimating propensity scores by optimizing balance. *Statistical Methods in
#' Medical Research*, 28(1), 84-101. \doi{10.1177/0962280217715487}
#'
#' Pirracchio R, Carone M (2018). The Balance Super Learner: A robust
#' adaptation of the Super Learner to improve estimation of the average
#' treatment effect in the treated based on propensity score matching.
#' *Statistical Methods in Medical Research*, 27(8), 2504-2518.
#' \doi{10.1177/0962280216682055}
#'
#' @seealso [psave()], [simplex_grid()], [psave_criteria()],
#'   `vignette("method-details", "psAve")`
#' @name psave-details
NULL
