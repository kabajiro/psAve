# psave.R -- the main (and only) estimation function: orchestration only.
# The heavy lifting lives in input.R, learners.R, grid.R, gamma.R, lambda.R,
# criteria.R.

#' Model-averaged propensity scores selected by prognostic-score balance
#'
#' `psave()` constructs a model-averaged propensity score
#' \eqn{\bar e(X) = \sum_m \lambda_m \hat e_m(X)}: a convex combination of
#' candidate propensity score models whose mixing weights \eqn{\lambda} are
#' selected on a simplex grid to optimize a balance criterion -- by default the
#' weighted absolute standardized mean difference of a **model-averaged
#' prognostic score** (the "Prog (Ave)" estimator of Kabata, Stuart and
#' Shintani 2024). The result is deliberately modest: a numeric score vector
#' designed to be handed to [MatchIt::matchit()] as `distance`, or to
#' [WeightIt::weightit()] as `ps`, with balance assessment via \pkg{cobalt}.
#'
#' @param formula A two-sided formula `treat ~ x1 + x2 + ...`, exactly as in
#'   [MatchIt::matchit()]. The right-hand side defines both the candidate-PS
#'   covariates and the balance covariates used by the `smd`/`ks` criteria.
#' @param data A data frame containing the variables in `formula` (and
#'   `outcome`). Complete cases in all used variables are REQUIRED; any
#'   missing value is an error, never a silent row drop.
#' @param outcome The outcome specification for the prognostic score: a
#'   one-sided formula `~ y` (the formula right-hand side is reused as the
#'   prognostic predictors) or a two-sided formula `y ~ z1 + z2` (a distinct
#'   prognostic specification). Required when `criterion = "prog"`; optional
#'   for the outcome-free criteria (`"smd"`, `"ks"`, `"logloss"`), where it may
#'   still be supplied so that the prognostic score is estimated and reported
#'   in `diagnostics`, `balance`, and [bal.tab.psave()].
#'   Prognostic models are fit on **untreated units only** (Hansen 2008), so
#'   using the outcome here does not bias effect estimation; see Details.
#' @param estimand `"ATT"` (default) or `"ATE"`; determines the
#'   inverse-probability weights used inside the balance criteria and returned
#'   in `weights`.
#' @param criterion The selection criterion for \eqn{\lambda}:
#'   `"prog"` (default) = weighted ASMD of the prognostic score (the paper's
#'   headline "Prog (Ave)"); `"smd"` = mean weighted ASMD over the covariates;
#'   `"ks"` = mean weighted Kolmogorov-Smirnov statistic over the covariates;
#'   `"logloss"` = negative Bernoulli log-likelihood of treatment assignment
#'   (the prediction-accuracy criterion in the lineage of Xie et al. 2019).
#' @param prog.target Only used when `criterion = "prog"`: `"average"`
#'   (default) targets the gamma-mixed prognostic score \eqn{\bar g}; naming a
#'   single element of `prog.methods` (or a column of `prog.matrix`) targets
#'   that candidate prognostic score instead (the paper's "Prog
#'   (\eqn{g_k})" variants).
#' @param ps.methods Character vector of candidate propensity score learners.
#'   The fixed, explicit default is `c("glm", "rpart", "ranger", "xgboost")`
#'   (no auto-detection: a machine-dependent candidate set is irreproducible
#'   science). Any `"SL.*"` string is accepted verbatim and passed through to
#'   \pkg{SuperLearner}. The order defines the tie-break preference (see
#'   Details).
#' @param prog.methods Character vector of candidate learners for the
#'   untreated-only prognostic models \eqn{g_k(0, X)}; same menu and default
#'   as `ps.methods`.
#' @param ps.matrix Optional n x M numeric matrix of user-supplied candidate
#'   propensity scores (values strictly in (0, 1), column names required).
#'   Overrides `ps.methods`; the columns are clipped like fitted candidates.
#' @param prog.matrix Optional n x K numeric matrix of user-supplied candidate
#'   prognostic scores (column names required). Overrides `prog.methods`.
#'   Requires `outcome` (gamma is selected by outcome-prediction MSE among
#'   untreated units).
#' @param ps.append Optional extra candidate propensity scores appended
#'   AFTER the candidates from `ps.methods` or `ps.matrix` (default `NULL`):
#'   a numeric vector of length n (one candidate, labeled `"append"`) or a
#'   numeric matrix / all-numeric data frame with n rows (unique, non-empty
#'   column names required, as for `ps.matrix`). Values must lie strictly in
#'   (0, 1) and each column is clipped to `clip` before averaging, exactly
#'   like every other candidate. Because appended columns come last, grid
#'   tie-breaking (first minimum) favors the base candidates (see Details).
#' @param prog.append Optional extra candidate prognostic scores appended
#'   AFTER the candidates from `prog.methods` or `prog.matrix` (default
#'   `NULL`); same vector/matrix/data-frame forms and naming rules as
#'   `ps.append`, but the values are unrestricted finite reals. Requires
#'   `outcome`.
#' @param average If `FALSE`, the \eqn{\lambda} grid is restricted to the
#'   simplex VERTICES, i.e., the best single PS candidate by `criterion` is
#'   selected (the "best single learner" variants of the paper's supplement).
#' @param family The prognostic-model family: `gaussian()` (default) or
#'   `binomial()` only. With `binomial()`, the gamma-selection MSE is the
#'   Brier score -- the formula is unchanged; note that the paper's
#'   simulations validated continuous outcomes.
#' @param step The simplex-grid increment for BOTH \eqn{\gamma} and
#'   \eqn{\lambda} (default `0.05`, the paper's value). Must evenly divide 1.
#' @param clip Length-2 numeric: each candidate propensity score column is
#'   clipped to `[clip[1], clip[2]]` BEFORE averaging (default
#'   `c(0.01, 0.99)`, the paper's constants). The average is never re-clipped:
#'   a convex combination of values in the clipping interval stays in it.
#' @param s.d.denom The group whose unweighted standard deviation standardizes
#'   mean differences in the ASMD-based criteria, passed to
#'   [cobalt::col_w_smd()]: `"treated"` (default; the paper's supplement uses
#'   the unweighted TREATED-group SD for both the ATT and the ATE),
#'   `"control"`, `"pooled"`, or `"all"`.
#' @param cv Number of cross-validation folds `V` for the
#'   `SuperLearner::SuperLearner()` `cvControl` when `"SL.*"` learners are
#'   used (default `5`).
#' @param control A named list of per-learner hyperparameter overrides, e.g.
#'   `list(ranger = list(num.trees = 1000), xgboost = list(nrounds = 200))`.
#'   Entries for `"glm"` are passed to [stats::glm()], for `"rpart"` to
#'   [rpart::rpart.control()], for `"ranger"` to [ranger::ranger()], and for
#'   `"xgboost"` to the `params` list of [xgboost::xgb.train()] (with
#'   `nrounds` and `verbose` recognized as top-level arguments). The resolved
#'   values are stored in `info$learners`.
#' @param keep.fits If `TRUE`, the fitted learner objects are retained in
#'   `fits`, enabling [predict.psave()]. Default `FALSE`.
#' @param keep.path If `TRUE` (default), the full \eqn{\lambda} criterion path
#'   is stored in `path`. Automatically set to `NULL` with a message if the
#'   grid exceeds 100,000 rows.
#' @param verbose If `TRUE`, progress messages report the learner set, the
#'   grid sizes, the selected \eqn{\lambda}/\eqn{\gamma}, and the criterion
#'   value.
#' @param ... Reserved for future use; supplying unused arguments triggers a
#'   warning.
#'
#' @details
#' ## Algorithm
#' 1. **Candidate propensity scores.** Each learner in `ps.methods` is fit on
#'    all \eqn{n} units and predicts \eqn{\hat e_m(X_i) = P(A_i = 1 \mid X_i)}
#'    in-sample; each column is clipped to `clip` before averaging. Extra
#'    user-supplied candidates given via `ps.append` are validated, clipped
#'    identically, and appended AFTER the base candidates (from `ps.methods`
#'    or `ps.matrix`), so the first-minimum tie-break favors the base set.
#' 2. **Model-averaged prognostic score.** Each learner in `prog.methods` is
#'    fit on the **untreated units only** and predicts \eqn{\hat g_k(0, X_i)}
#'    for all \eqn{n} units; `prog.append` columns are appended after these
#'    candidates. The mixing weights \eqn{\gamma} minimize the
#'    unweighted untreated-set mean squared error
#'    \deqn{\mathrm{MSE}(\gamma) = \mathrm{mean}_{i:A_i=0}\bigl(Y_i - \textstyle\sum_k \gamma_k \hat g_k(0, X_i)\bigr)^2}
#'    over `simplex_grid(K, step)`; \eqn{\bar g = \sum_k \gamma_k \hat g_k}.
#' 3. **Model-averaged propensity score.** The mixing weights \eqn{\lambda}
#'    minimize `criterion` over `simplex_grid(M, step)`, where each grid row
#'    implies \eqn{\bar e_\lambda = \sum_m \lambda_m \hat e_m} and the weights
#'    \eqn{W_i} below.
#'
#' ## Inverse-probability weights (exact supplement formulas)
#' With \eqn{\bar e_i} the averaged propensity score:
#' \deqn{\mathrm{ATT:}\quad W_i = 1 \; (A_i = 1), \qquad W_i = \bar e_i / (1 - \bar e_i) \; (A_i = 0);}
#' \deqn{\mathrm{ATE:}\quad W_i = 1 / \bar e_i \; (A_i = 1), \qquad W_i = 1 / (1 - \bar e_i) \; (A_i = 0).}
#' These are identical to [WeightIt::get_w_from_ps()] at the same estimand.
#'
#' ## Selection criteria
#' At each grid row \eqn{\lambda}:
#' \describe{
#'   \item{`logloss`}{\eqn{-\mathrm{mean}\{A_i \log \bar e_i + (1 - A_i) \log(1 - \bar e_i)\}}
#'     (finite by clipping).}
#'   \item{`smd`}{the mean over covariate columns \eqn{j} of
#'     \deqn{\mathrm{ASMD}_j = \Bigl| \frac{\sum_{A_i=1} W_i X_{ij}}{\sum_{A_i=1} W_i} - \frac{\sum_{A_i=0} W_i X_{ij}}{\sum_{A_i=0} W_i} \Bigr| / s_j,}
#'     where \eqn{s_j} is the **unweighted sample SD** of \eqn{X_j} in the
#'     `s.d.denom` group (the treated group for both estimands, per the
#'     paper's supplement).}
#'   \item{`ks`}{the mean over covariates of the proper weighted-eCDF
#'     Kolmogorov-Smirnov statistic
#'     \eqn{\sup_x |F^w_1(x) - F^w_0(x)|} with
#'     \eqn{F^w_a(x) = \sum_{i:A_i=a} W_i 1(X_{ij} \le x) / \sum_{i:A_i=a} W_i};
#'     for binary columns this equals the absolute difference in weighted
#'     proportions.}
#'   \item{`prog`}{the \eqn{\mathrm{ASMD}} formula applied to the single
#'     column \eqn{\bar g} (`prog.target = "average"`) or \eqn{\hat g_k}
#'     (`prog.target` names a learner). THE DEFAULT.}
#' }
#' When `prog.target` names a single prognostic candidate, `criterion.value`
#' is the weighted ASMD of that candidate's prognostic score, whereas the
#' `prog` column of the `diagnostics` table always refers to the
#' model-averaged prognostic score \eqn{\bar g}; the two values therefore
#' need not coincide, and [summary.psave()] prints a footnote to this effect.
#' Faithful to the published method, the `smd` and `prog` criteria standardize
#' **all** columns -- including binary ones -- by the plain unweighted sample
#' SD (uniform sample-SD standardization; `bin.vars = FALSE` is passed to
#' [cobalt::col_w_smd()] for every column). The display-oriented `balance`
#' component and [bal.tab.psave()] use \pkg{cobalt}'s native conventions
#' instead; see `vignette("method-details", "psAve")`.
#'
#' ## Simplex grid and tie-breaking
#' Both grids are enumerated by [simplex_grid()] in integer arithmetic
#' (every valid grid point is present by construction) in a documented order:
#' the first component descends from 1 to 0, then the second on the remainder,
#' and so on -- the first row puts all weight on the first candidate. Ties in
#' any argmin take the FIRST row attaining the minimum (within a 1e-9 relative
#' numerical tolerance), so ties favor learners listed earlier in
#' `ps.methods`/`prog.methods`; candidates appended via
#' `ps.append`/`prog.append` always come last and therefore lose ties to the
#' base candidates. The tolerance is deliberate: the criterion
#' values come from floating-point matrix algebra whose lowest-order bits can
#' differ across BLAS implementations, so an exact bitwise `which.min()`
#' would not be reproducible across machines, whereas the tolerant
#' first-minimum rule is. `gamma` always minimizes the unweighted
#' untreated-set MSE; these rules are fixed and not arguments.
#'
#' ## Why using the outcome does not bias the design
#' Prognostic models see the outcomes of **untreated units only**, and the
#' criterion compares weighted *covariate-like summaries* (the prognostic
#' score) between arms -- it never uses a treated-untreated outcome contrast
#' (Hansen 2008; Stuart, Lee and Leacy 2013). This is the same argument under
#' which prognostic-score balance diagnostics are recommended for propensity
#' score analyses; see `vignette("method-details", "psAve")`.
#'
#' ## Reproducibility
#' There is no `seed` argument: call `set.seed()` before `psave()` --
#' stochastic learners are `ranger` and `xgboost`. `info$learners` records the
#' resolved hyperparameters and engine package versions.
#'
#' ## Relation to other software
#' [WeightIt::method_super] with `SL.method = "method.balance"` (Pirracchio
#' and Carone 2018) selects a SuperLearner combination by *covariate* balance
#' and is weighting-only. `psave()` targets *prognostic-score* balance on an
#' exhaustive simplex grid, and returns a score vector equally usable for
#' matching and weighting. See [psave-details] and
#' `vignette("method-details", "psAve")` for the differences from the paper's
#' reference implementation.
#'
#' @return An object of class `"psave"`: a list with components
#' \describe{
#'   \item{`ps`}{numeric(n), named by `rownames(data)`: the model-averaged
#'     propensity score \eqn{\bar e} -- the deliverable for
#'     `MatchIt::matchit(distance = )` / `WeightIt::weightit(ps = )`.}
#'   \item{`prog`}{numeric(n), named: the model-averaged prognostic score
#'     \eqn{\bar g}; `NULL` when `outcome` was not supplied (possible only for
#'     the outcome-free criteria `"smd"`, `"ks"`, and `"logloss"`).}
#'   \item{`lambda`}{named numeric(M): the selected PS mixing weights, where
#'     M counts the base candidates plus any `ps.append` columns.}
#'   \item{`gamma`}{named numeric(K) (or `NULL`): the selected prognostic
#'     mixing weights, where K counts any `prog.append` columns.}
#'   \item{`weights`}{numeric(n): the IPW at `estimand` implied by `ps`.}
#'   \item{`ps.candidates`}{n x M matrix of clipped candidate propensity
#'     scores (fitted or user-supplied; `ps.append` columns come last).}
#'   \item{`prog.candidates`}{n x K matrix of candidate prognostic predictions
#'     for all n units (or `NULL`); `prog.append` columns come last.}
#'   \item{`criterion`, `criterion.value`}{the criterion used and its value at
#'     the selected \eqn{\lambda} (paper-faithful standardization).}
#'   \item{`diagnostics`}{(M+1) x 4 data frame: all four criteria
#'     ([psave_criteria()]) for each candidate and for the selected average --
#'     the "was averaging worth it?" table. Its `prog` column always uses the
#'     model-averaged prognostic score \eqn{\bar g}, even when `prog.target`
#'     names a single candidate (see Details); it then need not equal
#'     `criterion.value`.}
#'   \item{`path`}{data frame of the full \eqn{\lambda} grid (M columns +
#'     `value` + logical `selected`) in grid order, or `NULL`.}
#'   \item{`gamma.mse`}{named numeric(K+1): untreated MSE per prognostic
#'     candidate and for the selected average (or `NULL`).}
#'   \item{`balance`}{data frame (covariates + `prog`): unweighted and
#'     weighted SMD and KS (cobalt display conventions).}
#'   \item{`treat`}{integer(n) 0/1 treatment as used.}
#'   \item{`covs`}{numeric n x p balance-covariate matrix with
#'     `attr(, "bin.vars")`.}
#'   \item{`estimand`, `s.d.denom`, `prog.target`, `average`}{scalars, as
#'     resolved.}
#'   \item{`outcome.name`}{name of the outcome variable, or `NA`.}
#'   \item{`formula`, `data`}{as supplied (they power [psave_match()],
#'     [psave_weight()] and [print.psave()]; note the memory cost of storing
#'     `data`).}
#'   \item{`fits`}{list of fitted learners iff `keep.fits = TRUE` (enables
#'     [predict.psave()]); otherwise `NULL`.}
#'   \item{`info`}{list: `step`, `clip`, `cv`, `family`, `grid.size`, `n`,
#'     `learners` (labels, resolved hyperparameters, engine versions),
#'     `psAve.version`.}
#'   \item{`call`}{the matched call.}
#' }
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
#' @seealso [psave_match()], [psave_weight()], [simplex_grid()],
#'   [psave_criteria()], [bal.tab.psave()], [plot.psave()], [summary.psave()],
#'   [predict.psave()]
#'
#' @examples
#' if (requireNamespace("MatchIt", quietly = TRUE)) {
#'   data("lalonde", package = "MatchIt")
#'
#'   # Paper-headline "Prog (Ave)", ATT. For speed, this example restricts the
#'   # candidate learners; the default is
#'   # ps.methods = prog.methods = c("glm", "rpart", "ranger", "xgboost").
#'   meths <- if (requireNamespace("rpart", quietly = TRUE)) c("glm", "rpart") else "glm"
#'   fit <- psave(treat ~ age + educ + race + married + nodegree + re74 + re75,
#'                data = lalonde, outcome = ~ re78,
#'                ps.methods = meths, prog.methods = meths)
#'   fit
#'   summary(fit)
#'
#'   # hand off to MatchIt (canonical explicit call):
#'   m <- MatchIt::matchit(treat ~ age + educ + race + married + nodegree +
#'                           re74 + re75,
#'                         data = lalonde, distance = fit$ps)
#'   # or, reusing the stored formula and data (no retyping hazard):
#'   m2 <- psave_match(fit)
#' }
#' @export
psave <- function(formula,
                  data,
                  outcome      = NULL,
                  estimand     = c("ATT", "ATE"),
                  criterion    = c("prog", "smd", "ks", "logloss"),
                  prog.target  = "average",
                  ps.methods   = c("glm", "rpart", "ranger", "xgboost"),
                  prog.methods = c("glm", "rpart", "ranger", "xgboost"),
                  ps.matrix    = NULL,
                  prog.matrix  = NULL,
                  ps.append    = NULL,
                  prog.append  = NULL,
                  average      = TRUE,
                  family       = gaussian(),
                  step         = 0.05,
                  clip         = c(0.01, 0.99),
                  s.d.denom    = "treated",
                  cv           = 5L,
                  control      = list(),
                  keep.fits    = FALSE,
                  keep.path    = TRUE,
                  verbose      = FALSE,
                  ...) {
  cl <- match.call()
  estimand <- match.arg(estimand)
  criterion <- match.arg(criterion)
  s.d.denom <- .validate_sd_denom(s.d.denom)

  dots <- list(...)
  if (length(dots)) {
    warning(sprintf("Ignoring unused argument(s): %s.",
                    paste0("`", names(dots), "`", collapse = ", ")),
            call. = FALSE)
  }

  ## --- scalar argument validation --------------------------------------------
  if (!is.numeric(step) || length(step) != 1L || is.na(step) ||
      step <= 0 || step > 1) {
    stop("`step` must be a single number in (0, 1] that evenly divides 1.",
         call. = FALSE)
  }
  if (!is.numeric(clip) || length(clip) != 2L || anyNA(clip) ||
      clip[1L] <= 0 || clip[2L] >= 1 || clip[1L] >= clip[2L]) {
    stop("`clip` must be a length-2 numeric vector with 0 < clip[1] < clip[2] < 1.",
         call. = FALSE)
  }
  for (nm in c("average", "keep.fits", "keep.path", "verbose")) {
    v <- get(nm, inherits = FALSE)
    if (!is.logical(v) || length(v) != 1L || is.na(v)) {
      stop(sprintf("`%s` must be TRUE or FALSE.", nm), call. = FALSE)
    }
  }
  if (!is.numeric(cv) || length(cv) != 1L || is.na(cv) || cv < 2) {
    stop("`cv` must be a single integer >= 2.", call. = FALSE)
  }
  cv <- as.integer(cv)
  if (is.character(family)) {
    family <- get(family, mode = "function", envir = parent.frame())
  }
  if (is.function(family)) family <- family()
  if (!inherits(family, "family") ||
      !family$family %in% c("gaussian", "binomial")) {
    stop("`family` must be gaussian() or binomial().", call. = FALSE)
  }
  if (!is.list(control) ||
      (length(control) && (is.null(names(control)) || any(!nzchar(names(control)))))) {
    stop("`control` must be a named list of per-learner option lists, e.g. list(ranger = list(num.trees = 1000)).",
         call. = FALSE)
  }

  ## --- inputs (D.1, D.2) ------------------------------------------------------
  inp <- .process_inputs(formula, data, outcome, criterion, family)
  n <- inp$n
  treat <- inp$treat
  .vmsg(verbose, sprintf("n = %d units (%d treated, %d control); %d balance covariate column(s).",
                         n, sum(treat == 1L), sum(treat == 0L), ncol(inp$covs)))

  ## --- candidate propensity scores (D.3, D.4) ---------------------------------
  ps.records <- NULL
  if (!is.null(ps.matrix)) {
    E <- .validate_score_matrix(ps.matrix, n, "ps.matrix", rn = inp$rn, prob = TRUE)
    ps.labels <- colnames(E)
    .vmsg(verbose, sprintf("Using %d user-supplied candidate PS column(s): %s.",
                           ncol(E), paste0('"', ps.labels, '"', collapse = ", ")))
  } else {
    ps.labels <- .validate_methods(ps.methods, "ps.methods")
    .vmsg(verbose, sprintf("PS candidate learners: %s.",
                           paste0('"', ps.labels, '"', collapse = ", ")))
    fit.ps <- .fit_candidates(ps.labels, y = treat, Xnum = inp$Xnum,
                              Xdf = inp$Xdf, family = family, task = "ps",
                              control = control, cv = cv, fit.idx = NULL,
                              verbose = verbose)
    E <- fit.ps$preds
    ps.records <- fit.ps$records
  }
  ## clip each candidate column BEFORE averaging; never re-clip after (D.4)
  E <- pmin(pmax(E, clip[1L]), clip[2L])
  ps.base.labels <- ps.labels
  ps.append.mat <- .normalize_append(ps.append, n, "ps.append", rn = inp$rn,
                                     prob = TRUE)
  if (!is.null(ps.append.mat)) {
    ## appended candidates are clipped exactly like every other candidate,
    ## then go AFTER the base candidates (tie-breaks favor the base set)
    ps.append.mat <- pmin(pmax(ps.append.mat, clip[1L]), clip[2L])
    E <- .append_candidates(E, ps.append.mat, "ps.append")
    ps.labels <- colnames(E)
    .vmsg(verbose, sprintf("Appended %d user-supplied candidate PS column(s): %s.",
                           ncol(ps.append.mat),
                           paste0('"', colnames(ps.append.mat), '"', collapse = ", ")))
  }
  rownames(E) <- inp$rn

  ## --- candidate prognostic scores + gamma (D.6) ------------------------------
  Gmat <- NULL
  prog.records <- NULL
  gsel <- NULL
  prog.base.labels <- NULL
  prog.append.mat <- NULL
  if (is.null(inp$y) && (!is.null(prog.matrix) || !is.null(prog.append))) {
    stop(paste0("`prog.matrix`/`prog.append` were supplied without `outcome`: ",
                "prognostic candidates require the outcome; gamma is selected by ",
                "outcome-prediction MSE among untreated units. Supply `outcome` ",
                "or drop the prognostic candidates."),
         call. = FALSE)
  }
  if (!is.null(inp$y)) {
    if (!is.null(prog.matrix)) {
      Gmat <- .validate_score_matrix(prog.matrix, n, "prog.matrix", rn = inp$rn,
                                     prob = FALSE)
      .vmsg(verbose, sprintf("Using %d user-supplied prognostic candidate column(s): %s.",
                             ncol(Gmat), paste0('"', colnames(Gmat), '"', collapse = ", ")))
    } else {
      prog.labels <- .validate_methods(prog.methods, "prog.methods")
      .vmsg(verbose, sprintf("Prognostic candidate learners (fit on the %d untreated units): %s.",
                             sum(treat == 0L),
                             paste0('"', prog.labels, '"', collapse = ", ")))
      fit.prog <- .fit_candidates(prog.labels, y = inp$y,
                                  Xnum = inp$prog.Xnum, Xdf = inp$prog.Xdf,
                                  family = family, task = "prog",
                                  control = control, cv = cv,
                                  fit.idx = which(treat == 0L),
                                  verbose = verbose)
      Gmat <- fit.prog$preds
      prog.records <- fit.prog$records
    }
    prog.base.labels <- colnames(Gmat)
    prog.append.mat <- .normalize_append(prog.append, n, "prog.append",
                                         rn = inp$rn, prob = FALSE)
    if (!is.null(prog.append.mat)) {
      Gmat <- .append_candidates(Gmat, prog.append.mat, "prog.append")
      .vmsg(verbose, sprintf("Appended %d user-supplied prognostic candidate column(s): %s.",
                             ncol(prog.append.mat),
                             paste0('"', colnames(prog.append.mat), '"', collapse = ", ")))
    }
    rownames(Gmat) <- inp$rn
    gsel <- .select_gamma(inp$y, Gmat, treat, step, verbose = verbose)
  }

  ## unused `control` entries -> warning (append labels are NOT fitted learners)
  used.labels <- unique(c(if (is.null(ps.matrix)) ps.base.labels,
                          if (!is.null(Gmat) && is.null(prog.matrix)) prog.base.labels))
  if (length(control)) {
    unused <- setdiff(names(control), used.labels)
    if (length(unused)) {
      warning(sprintf("`control` entries for %s do not match any fitted learner and were ignored.",
                      paste0('"', unused, '"', collapse = ", ")), call. = FALSE)
    }
  }

  ## --- prog.target resolution (B.1) -------------------------------------------
  target.vec <- NULL
  if (criterion == "prog") {
    if (!is.character(prog.target) || length(prog.target) != 1L || is.na(prog.target)) {
      stop("`prog.target` must be \"average\" or the label of one prognostic candidate.",
           call. = FALSE)
    }
    if (identical(prog.target, "average")) {
      target.vec <- gsel$gbar
    } else if (prog.target %in% colnames(Gmat)) {
      target.vec <- as.numeric(Gmat[, prog.target])
    } else {
      stop(sprintf("`prog.target` must be \"average\" or one of the prognostic candidate labels (%s); got \"%s\".",
                   paste0('"', colnames(Gmat), '"', collapse = ", "), prog.target),
           call. = FALSE)
    }
  } else {
    prog.target <- NA_character_
  }

  ## --- lambda selection (D.5, D.8, D.9) ---------------------------------------
  lsel <- .select_lambda(E, treat, covs = inp$covs, prog.vec = target.vec,
                         criterion = criterion, estimand = estimand,
                         s.d.denom = s.d.denom, step = step,
                         average = average, keep.path = keep.path,
                         verbose = verbose)

  ## --- assemble the averaged score, weights, diagnostics, balance (D.10) ------
  ps <- as.numeric(E %*% lsel$lambda)
  names(ps) <- inp$rn
  gbar <- if (!is.null(gsel)) stats::setNames(gsel$gbar, inp$rn) else NULL
  W <- .ipw_from_ps(ps, treat, estimand)

  M <- ncol(E)
  diag.mat <- t(vapply(seq_len(M), function(m) {
    psave_criteria(E[, m], treat, inp$covs,
                   prog = if (is.null(gsel)) NULL else gsel$gbar,
                   estimand = estimand, s.d.denom = s.d.denom,
                   bin.vars = inp$bin.vars)
  }, numeric(4L)))
  diag.mat <- rbind(diag.mat,
                    psave_criteria(ps, treat, inp$covs,
                                   prog = if (is.null(gsel)) NULL else gsel$gbar,
                                   estimand = estimand, s.d.denom = s.d.denom,
                                   bin.vars = inp$bin.vars))
  diagnostics <- as.data.frame(diag.mat, row.names = c(ps.labels, "average"))

  ## display balance table: cobalt native conventions (H.1)
  balance <- .balance_table(inp$covs, treat, W, s.d.denom, inp$bin.vars,
                            prog = if (is.null(gsel)) NULL else gsel$gbar)

  covs.out <- inp$covs
  attr(covs.out, "bin.vars") <- inp$bin.vars

  fits <- NULL
  if (keep.fits) {
    fits <- list(ps = ps.records,
                 prog = prog.records,
                 design = list(ps = list(terms = inp$ps.tt, xlev = inp$ps.xlev),
                               prog = if (!is.null(inp$prog.tt))
                                 list(terms = inp$prog.tt, xlev = inp$prog.xlev)
                               else NULL))
  }

  ## learner provenance, with appended candidates recorded as user-supplied
  learners.ps <- if (!is.null(ps.records)) .learner_info(ps.records)
                 else "user-supplied ps.matrix"
  if (!is.null(ps.append.mat)) {
    learners.ps <- c(if (is.character(learners.ps)) list(base = learners.ps)
                     else learners.ps,
                     .append_info(colnames(ps.append.mat)))
  }
  learners.prog <- if (!is.null(prog.records)) .learner_info(prog.records)
                   else if (!is.null(Gmat)) "user-supplied prog.matrix"
                   else NULL
  if (!is.null(prog.append.mat)) {
    learners.prog <- c(if (is.character(learners.prog)) list(base = learners.prog)
                       else learners.prog,
                       .append_info(colnames(prog.append.mat)))
  }

  info <- list(step = step,
               clip = clip,
               cv = cv,
               family = family$family,
               grid.size = c(lambda = lsel$grid.size,
                             gamma = if (!is.null(gsel)) gsel$grid.size else NA_integer_),
               n = n,
               learners = list(ps = learners.ps,
                               prog = learners.prog),
               psAve.version = tryCatch(as.character(utils::packageVersion("psAve")),
                                        error = function(e) NA_character_))

  res <- structure(list(ps = ps,
                        prog = gbar,
                        lambda = lsel$lambda,
                        gamma = if (!is.null(gsel)) gsel$gamma else NULL,
                        weights = W,
                        ps.candidates = E,
                        prog.candidates = Gmat,
                        criterion = criterion,
                        criterion.value = lsel$value,
                        diagnostics = diagnostics,
                        path = lsel$path,
                        gamma.mse = if (!is.null(gsel)) gsel$gamma.mse else NULL,
                        balance = balance,
                        treat = treat,
                        covs = covs.out,
                        estimand = estimand,
                        s.d.denom = s.d.denom,
                        prog.target = prog.target,
                        average = average,
                        outcome.name = inp$outcome.name,
                        formula = formula,
                        data = inp$data,
                        fits = fits,
                        info = info,
                        call = cl),
                   class = "psave")
  res
}

# Display-oriented balance table (B.2 `balance`, D.10): cobalt native
# conventions (auto-detected bin.vars; binary columns standardized by
# sqrt(p(1-p))), absolute values. Rows = covariates + "prog".
.balance_table <- function(covs, treat, W, s.d.denom, bin.vars, prog = NULL) {
  smd.un <- cobalt::col_w_smd(covs, treat = treat, weights = NULL, std = TRUE,
                              s.d.denom = s.d.denom, abs = TRUE,
                              bin.vars = bin.vars)
  smd.wt <- cobalt::col_w_smd(covs, treat = treat, weights = W, std = TRUE,
                              s.d.denom = s.d.denom, abs = TRUE,
                              bin.vars = bin.vars)
  ks.un <- cobalt::col_w_ks(covs, treat = treat, weights = NULL,
                            bin.vars = bin.vars)
  ks.wt <- cobalt::col_w_ks(covs, treat = treat, weights = W,
                            bin.vars = bin.vars)
  rows <- colnames(covs)
  if (!is.null(prog)) {
    g <- matrix(as.numeric(prog), ncol = 1L)
    smd.un <- c(smd.un, cobalt::col_w_smd(g, treat = treat, weights = NULL,
                                          std = TRUE, s.d.denom = s.d.denom,
                                          abs = TRUE, bin.vars = FALSE))
    smd.wt <- c(smd.wt, cobalt::col_w_smd(g, treat = treat, weights = W,
                                          std = TRUE, s.d.denom = s.d.denom,
                                          abs = TRUE, bin.vars = FALSE))
    ks.un <- c(ks.un, cobalt::col_w_ks(g, treat = treat, weights = NULL,
                                       bin.vars = FALSE))
    ks.wt <- c(ks.wt, cobalt::col_w_ks(g, treat = treat, weights = W,
                                       bin.vars = FALSE))
    rows <- c(rows, "prog")
  }
  data.frame(smd.un = unname(smd.un), smd.wt = unname(smd.wt),
             ks.un = unname(ks.un), ks.wt = unname(ks.wt),
             row.names = rows)
}
