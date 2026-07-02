# criteria.R -- the exported psave_criteria() plus the internal criterion
# kernels shared with the lambda grid search (D.7, D.8, H.1).

# IPW weights implied by a propensity score (exact supplement formulas, D.7):
#   ATT: W = 1 if A = 1;  e/(1-e)  if A = 0.
#   ATE: W = 1/e if A = 1; 1/(1-e) if A = 0.
.ipw_from_ps <- function(ps, treat, estimand) {
  w <- numeric(length(ps))
  t1 <- treat == 1L
  if (estimand == "ATT") {
    w[t1] <- 1
    w[!t1] <- ps[!t1] / (1 - ps[!t1])
  } else {
    w[t1] <- 1 / ps[t1]
    w[!t1] <- 1 / (1 - ps[!t1])
  }
  w
}

#' Compute all four psAve selection criteria for a propensity score vector
#'
#' Evaluates, for an arbitrary propensity score vector, the four selection
#' criteria used by [psave()]: the treatment-assignment log loss, the mean
#' weighted absolute standardized mean difference (ASMD) of the covariates,
#' the mean weighted Kolmogorov-Smirnov (KS) statistic of the covariates, and
#' the weighted ASMD of a prognostic score. This function powers the
#' `diagnostics` table of a `psave` object (the "was averaging worth it?"
#' comparison) and is exported as a methods-research utility.
#'
#' @param ps Numeric vector of propensity scores, strictly inside (0, 1).
#' @param treat Treatment vector; coerced to 0/1 like the left-hand side of
#'   the [psave()] formula (numeric 0/1, logical, or two-level
#'   factor/character with the second level treated).
#' @param covs Numeric matrix (or all-numeric data frame) of covariates, one
#'   row per unit. Factors must already be expanded to dummy columns (as in
#'   the `covs` component of a `psave` object).
#' @param prog Optional numeric vector: a prognostic score. If `NULL`, the
#'   `prog` criterion is returned as `NA`.
#' @param estimand `"ATT"` (default) or `"ATE"`; determines the weights (see
#'   Details).
#' @param s.d.denom Group whose (unweighted) standard deviation standardizes
#'   the mean differences: `"treated"` (default; the paper's convention for
#'   BOTH estimands), `"control"`, `"pooled"`, or `"all"`.
#' @param bin.vars Optional logical vector flagging binary columns of `covs`,
#'   used for the KS criterion only; if `NULL`, columns with exactly two
#'   distinct values are detected automatically. The `smd` and `prog` criteria
#'   always use uniform sample-SD standardization (`bin.vars = FALSE` for all
#'   columns), which is the convention of the published method (see Details).
#'
#' @details
#' Weights are the inverse-probability weights implied by `ps` at `estimand`:
#' for the ATT, \eqn{W_i = 1} for treated units and \eqn{e_i/(1-e_i)} for
#' untreated units; for the ATE, \eqn{1/e_i} and \eqn{1/(1-e_i)}. The four
#' criteria are:
#' \describe{
#'   \item{`logloss`}{\eqn{-\mathrm{mean}\{A_i \log e_i + (1-A_i)\log(1-e_i)\}}.}
#'   \item{`smd`}{the mean over covariates \eqn{j} of
#'     \eqn{|\bar X_{1j}^w - \bar X_{0j}^w| / s_j}, where \eqn{\bar X_{aj}^w}
#'     is the weighted mean of \eqn{X_j} in arm \eqn{a} and \eqn{s_j} is the
#'     unweighted sample SD of \eqn{X_j} in the `s.d.denom` group. Computed via
#'     [cobalt::col_w_smd()].}
#'   \item{`ks`}{the mean over covariates of the proper weighted-eCDF KS
#'     statistic \eqn{\sup_x |F^w_1(x) - F^w_0(x)|}, computed via
#'     [cobalt::col_w_ks()]; for binary columns this is the absolute
#'     difference in weighted proportions.}
#'   \item{`prog`}{the same weighted ASMD formula applied to the single
#'     column `prog`.}
#' }
#'
#' Faithful to the published method (and its reference implementation), the
#' `smd` and `prog` criteria standardize **every** column, including binary
#' ones, by the plain unweighted sample SD (`sd()`, the \eqn{n-1} formula) of
#' the `s.d.denom` group -- i.e., `bin.vars = FALSE` is passed to
#' [cobalt::col_w_smd()] for all columns. cobalt's own display convention
#' (binary columns standardized by \eqn{\sqrt{p(1-p)}}) is used only in the
#' display-oriented `balance` component of a `psave` object and in
#' [bal.tab.psave()]. For the KS criterion the two conventions coincide.
#'
#' @return A named numeric vector with elements `logloss`, `smd`, `ks`, and
#'   `prog` (the last is `NA` when `prog = NULL`).
#'
#' @references
#' Kabata D, Stuart EA, Shintani A (2024). Prognostic score-based model
#' averaging approach for propensity score estimation. *BMC Medical Research
#' Methodology*, 24, 228. \doi{10.1186/s12874-024-02350-y}
#'
#' Hansen BB (2008). The prognostic analogue of the propensity score.
#' *Biometrika*, 95(2), 481-488. \doi{10.1093/biomet/asn004}
#'
#' Xie Y, Zhu Y, Cotton CA, Wu P (2019). A model averaging approach for
#' estimating propensity scores by optimizing balance. *Statistical Methods in
#' Medical Research*, 28(1), 84-101. \doi{10.1177/0962280217715487}
#'
#' @seealso [psave()], [simplex_grid()], [cobalt::col_w_smd()],
#'   [cobalt::col_w_ks()]
#' @examples
#' set.seed(1)
#' n <- 200
#' x1 <- rnorm(n); x2 <- rbinom(n, 1, 0.4)
#' a <- rbinom(n, 1, plogis(-0.5 + x1 + 0.5 * x2))
#' ps <- pmin(pmax(fitted(glm(a ~ x1 + x2, family = binomial())), 0.01), 0.99)
#' g <- 1 + 0.5 * x1 - 0.2 * x2   # a (toy) prognostic score
#' psave_criteria(ps, a, cbind(x1 = x1, x2 = x2), prog = g, estimand = "ATT")
#' @export
psave_criteria <- function(ps, treat, covs, prog = NULL,
                           estimand = c("ATT", "ATE"),
                           s.d.denom = "treated", bin.vars = NULL) {
  estimand <- match.arg(estimand)
  s.d.denom <- .validate_sd_denom(s.d.denom)
  treat <- .coerce_treat(treat)

  if (is.data.frame(covs)) {
    if (!all(vapply(covs, is.numeric, logical(1L)))) {
      stop(paste0("`covs` must be numeric: expand factors to dummy columns first ",
                  "(e.g., with model.matrix() or cobalt::splitfactor())."),
           call. = FALSE)
    }
    covs <- as.matrix(covs)
  }
  if (!is.matrix(covs) || !is.numeric(covs)) {
    stop("`covs` must be a numeric matrix or an all-numeric data frame.", call. = FALSE)
  }
  n <- length(treat)
  if (nrow(covs) != n) {
    stop("`covs` must have one row per element of `treat`.", call. = FALSE)
  }
  if (!is.numeric(ps) || length(ps) != n) {
    stop("`ps` must be a numeric vector with one value per unit.", call. = FALSE)
  }
  if (anyNA(ps) || any(!is.finite(ps)) || any(ps <= 0) || any(ps >= 1)) {
    stop("All values of `ps` must lie strictly inside (0, 1) (clip extreme scores first).",
         call. = FALSE)
  }
  p <- ncol(covs)
  if (is.null(bin.vars)) {
    bin.vars <- .detect_bin(covs)
  } else {
    if (length(bin.vars) == 1L) bin.vars <- rep(bin.vars, p)
    if (!is.logical(bin.vars) || length(bin.vars) != p || anyNA(bin.vars)) {
      stop("`bin.vars` must be a logical vector with one entry per column of `covs`.",
           call. = FALSE)
    }
  }
  if (!is.null(prog)) {
    prog <- as.numeric(prog)
    if (length(prog) != n || anyNA(prog) || any(!is.finite(prog))) {
      stop("`prog` must be a finite numeric vector with one value per unit.",
           call. = FALSE)
    }
  }

  W <- .ipw_from_ps(ps, treat, estimand)

  ll <- -mean(treat * log(ps) + (1 - treat) * log1p(-ps))
  smd <- mean(cobalt::col_w_smd(covs, treat = treat, weights = W, std = TRUE,
                                s.d.denom = s.d.denom, abs = TRUE,
                                bin.vars = rep(FALSE, p)))
  ks <- mean(cobalt::col_w_ks(covs, treat = treat, weights = W,
                              bin.vars = bin.vars))
  pr <- if (is.null(prog)) NA_real_ else
    unname(cobalt::col_w_smd(matrix(prog, ncol = 1L), treat = treat,
                             weights = W, std = TRUE, s.d.denom = s.d.denom,
                             abs = TRUE, bin.vars = FALSE))

  c(logloss = ll, smd = smd, ks = ks, prog = pr)
}
