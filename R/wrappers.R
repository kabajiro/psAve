# wrappers.R -- thin, guarded pass-throughs to MatchIt and WeightIt (D.11).
# They add NO arguments of their own and reuse the formula and data stored in
# the psave object, which eliminates the row-misalignment hazard of retyping
# `data =` by construction. The explicit two-step call remains canonical.

#' Match on a model-averaged propensity score
#'
#' Convenience pass-through to [MatchIt::matchit()]: matches on the
#' model-averaged propensity score of a [psave()] fit, reusing the formula and
#' data stored in the object. Equivalent to the canonical explicit call
#' \preformatted{MatchIt::matchit(<formula>, data = <data>, distance = fit$ps, ...)}
#' but with no opportunity for row misalignment between the two steps. All
#' `...` arguments are forwarded verbatim; the return value is an ordinary
#' `matchit` object, so the full \pkg{MatchIt}/\pkg{cobalt} toolkit applies.
#'
#' @param x A `psave` object.
#' @param ... Arguments forwarded verbatim to [MatchIt::matchit()] (e.g.,
#'   `method`, `caliper`, `ratio`, `replace`).
#'
#' @return A `matchit` object; see [MatchIt::matchit()].
#' @seealso [psave()], [psave_weight()], [MatchIt::matchit()]
#' @examplesIf requireNamespace("MatchIt", quietly = TRUE)
#' data("lalonde", package = "MatchIt")
#' fit <- psave(treat ~ age + educ + married + re74, data = lalonde,
#'              outcome = ~ re78, ps.methods = "glm", prog.methods = "glm")
#' m <- psave_match(fit, method = "nearest", caliper = 0.2)
#' cobalt::bal.tab(m, distance = data.frame(prog = fit$prog))
#' @export
psave_match <- function(x, ...) {
  if (!inherits(x, "psave")) {
    stop("`x` must be a psave object (the result of psave()).", call. = FALSE)
  }
  .require_pkg("MatchIt", "for psave_match()")
  MatchIt::matchit(formula = x$formula, data = x$data, distance = x$ps, ...)
}

#' Weight by a model-averaged propensity score
#'
#' Convenience pass-through to [WeightIt::weightit()]: constructs balancing
#' weights from the model-averaged propensity score of a [psave()] fit at the
#' fitted estimand, reusing the stored formula and data. Equivalent to the
#' canonical explicit call
#' \preformatted{WeightIt::weightit(<formula>, data = <data>, ps = fit$ps, estimand = fit$estimand, ...)}
#' All `...` arguments are forwarded verbatim; the return value is an ordinary
#' `weightit` object.
#'
#' @param x A `psave` object.
#' @param ... Arguments forwarded verbatim to [WeightIt::weightit()].
#'
#' @return A `weightit` object; see [WeightIt::weightit()].
#' @seealso [psave()], [psave_match()], [WeightIt::weightit()],
#'   [WeightIt::get_w_from_ps()]
#' @examplesIf requireNamespace("MatchIt", quietly = TRUE) && requireNamespace("WeightIt", quietly = TRUE)
#' data("lalonde", package = "MatchIt")
#' fit <- psave(treat ~ age + educ + married + re74, data = lalonde,
#'              outcome = ~ re78, ps.methods = "glm", prog.methods = "glm")
#' w <- psave_weight(fit)
#' cobalt::bal.tab(w, distance = data.frame(prog = fit$prog))
#' @export
psave_weight <- function(x, ...) {
  if (!inherits(x, "psave")) {
    stop("`x` must be a psave object (the result of psave()).", call. = FALSE)
  }
  .require_pkg("WeightIt", "for psave_weight()")
  WeightIt::weightit(formula = x$formula, data = x$data, ps = x$ps,
                     estimand = x$estimand, ...)
}
