# baltab.R -- bal.tab method for psave objects, registered on cobalt's
# generic (D.12) so that cobalt::bal.tab(fit) "just works".

#' Balance tables for psave objects
#'
#' A method for [cobalt::bal.tab()]: assesses balance on the covariates of a
#' [psave()] fit under the implied inverse-probability weights, with the
#' model-averaged propensity score and (when available) the model-averaged
#' prognostic score supplied as `distance` measures -- the prognostic-score
#' balance diagnostic of Stuart, Lee and Leacy (2013).
#'
#' The call delegates to the default \pkg{cobalt} machinery as
#' `cobalt::bal.tab(<covariates>, treat = x$treat, weights = x$weights,
#' s.d.denom = x$s.d.denom, distance = data.frame(ps = x$ps, prog = x$prog),
#' ...)`, so all the usual \pkg{cobalt} arguments (`un`, `stats`,
#' `thresholds`, ...) are available, and display conventions are
#' \pkg{cobalt}'s own (the *selection criterion* inside [psave()] uses the
#' paper's uniform sample-SD standardization instead; see
#' [psave_criteria()]).
#'
#' @param x A `psave` object.
#' @param ... Further arguments passed on to [cobalt::bal.tab()] (e.g.,
#'   `un = TRUE`, `thresholds = c(m = 0.1)`).
#'
#' @return A `bal.tab` object; see [cobalt::bal.tab()].
#'
#' @references
#' Stuart EA, Lee BK, Leacy FP (2013). Prognostic score-based balance
#' measures can be a useful diagnostic for propensity score methods in
#' comparative effectiveness research. *Journal of Clinical Epidemiology*,
#' 66(8), S84-S90. \doi{10.1016/j.jclinepi.2013.01.013}
#'
#' @seealso [psave()], [cobalt::bal.tab()], [plot.psave()]
#' @examplesIf requireNamespace("MatchIt", quietly = TRUE)
#' data("lalonde", package = "MatchIt")
#' fit <- psave(treat ~ age + educ + married + re74, data = lalonde,
#'              outcome = ~ re78, ps.methods = "glm", prog.methods = "glm")
#' cobalt::bal.tab(fit, un = TRUE)
#' @exportS3Method cobalt::bal.tab
bal.tab.psave <- function(x, ...) {
  d <- data.frame(ps = as.numeric(x$ps))
  if (!is.null(x$prog)) d$prog <- as.numeric(x$prog)
  cobalt::bal.tab(as.data.frame(x$covs),
                  treat = x$treat,
                  weights = x$weights,
                  s.d.denom = x$s.d.denom,
                  distance = d,
                  ...)
}
