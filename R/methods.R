# methods.R -- S3 methods for psave objects: print, summary, print.summary,
# fitted, weights, predict (B.3).

#' Print a psave object
#'
#' Prints a one-screen summary of a fitted [psave()] object: estimand and
#' criterion, the selected mixing weights \eqn{\lambda} and \eqn{\gamma} as
#' labeled text bars, the criterion value, a three-row balance preview (the
#' worst covariates plus the prognostic score), and then the **literal next
#' call** -- echoing the formula and data name from your own [psave()] call --
#' that hands the averaged score to [MatchIt::matchit()] or
#' [WeightIt::weightit()].
#'
#' @param x A `psave` object.
#' @param digits Number of significant digits to print. Default 3.
#' @param ... Ignored.
#'
#' @return `x`, invisibly.
#' @seealso [psave()], [summary.psave()]
#' @export
print.psave <- function(x, digits = 3, ...) {
  sub <- substitute(x)
  obj <- if (is.name(sub)) as.character(sub) else "fit"

  n1 <- sum(x$treat == 1L)
  n0 <- sum(x$treat == 0L)
  cat("A psave object (model-averaged propensity score)\n")
  cat(sprintf(" - estimand:  %s\n", x$estimand))
  cat(sprintf(" - criterion: %s (%s)%s\n", x$criterion,
              .criterion_label(x$criterion, x$prog.target),
              if (isTRUE(x$average)) "" else " [average = FALSE: best single candidate]"))
  cat(sprintf(" - sample:    %d units (%d treated, %d control)\n",
              x$info$n, n1, n0))

  cat("\nlambda (PS mixing weights):\n")
  .print_weight_bars(x$lambda, digits)
  if (!is.null(x$gamma)) {
    cat("\ngamma (prognostic mixing weights):\n")
    .print_weight_bars(x$gamma, digits)
  }
  cat(sprintf("\nCriterion value at selected lambda: %s\n",
              format(x$criterion.value, digits = digits)))

  ## 3-row balance preview: worst covariates by weighted SMD + prog
  b <- x$balance
  covrows <- setdiff(rownames(b), "prog")
  worst <- covrows[order(-b[covrows, "smd.wt"])]
  show <- c(utils::head(worst, 3L), intersect("prog", rownames(b)))
  cat("\nBalance preview (worst covariates + prognostic score):\n")
  print(round(b[show, , drop = FALSE], digits))

  ## the literal next call, echoing the user's own formula/data symbols
  f.txt <- tryCatch(deparse1(x$call$formula), error = function(e) NULL)
  if (is.null(f.txt) || !nzchar(f.txt) || f.txt == "NULL") f.txt <- deparse1(x$formula)
  d.txt <- tryCatch(deparse1(x$call$data), error = function(e) "data")
  cat("\nNext:\n")
  cat(sprintf("  MatchIt::matchit(%s, data = %s, distance = %s$ps)\n",
              f.txt, d.txt, obj))
  cat(sprintf("    or: psave_match(%s)\n", obj))
  cat(sprintf("  WeightIt::weightit(%s, data = %s, ps = %s$ps, estimand = \"%s\")\n",
              f.txt, d.txt, obj, x$estimand))
  cat(sprintf("    or: psave_weight(%s)\n", obj))
  invisible(x)
}

# Text-bar rendering of a named weight vector.
.print_weight_bars <- function(w, digits = 3) {
  labs <- format(names(w))
  vals <- format(round(w, digits), nsmall = min(digits, 3L))
  for (i in seq_along(w)) {
    cat(sprintf("  %s  %s  %s\n", labs[i], vals[i], .text_bar(w[i])))
  }
  invisible(NULL)
}

#' Summarize a psave object
#'
#' Produces (a) the selected mixing-weight tables \eqn{\lambda} and
#' \eqn{\gamma}, (b) the `diagnostics` table (all four selection criteria for
#' every candidate propensity score and for the selected average -- the "was
#' averaging worth it?" comparison), and (c) the full balance table (all
#' covariates plus the prognostic score; unweighted vs. weighted SMD and KS,
#' with a `*` marker at weighted SMD > 0.1).
#'
#' @param object A `psave` object.
#' @param un If `TRUE` (default), the balance table includes the unweighted
#'   columns.
#' @param candidates If `TRUE` (default), the per-candidate `diagnostics`
#'   table is included.
#' @param ... Ignored.
#' @param x A `summary.psave` object.
#' @param digits Number of significant digits to print. Default 3.
#'
#' @return For `summary.psave()`, an object of class `"summary.psave"`: a
#'   list with elements `lambda`, `gamma`, `gamma.mse`, `diagnostics`,
#'   `balance`, `criterion`, `criterion.value`, `prog.target`, `estimand`,
#'   `average`, `nn`, and `call`. `print.summary.psave()` returns `x`
#'   invisibly.
#' @seealso [psave()], [print.psave()], [bal.tab.psave()]
#' @export
summary.psave <- function(object, un = TRUE, candidates = TRUE, ...) {
  balance <- object$balance
  if (!isTRUE(un)) {
    balance <- balance[, c("smd.wt", "ks.wt"), drop = FALSE]
  }
  out <- list(call = object$call,
              estimand = object$estimand,
              criterion = object$criterion,
              prog.target = object$prog.target,
              average = object$average,
              criterion.value = object$criterion.value,
              lambda = object$lambda,
              gamma = object$gamma,
              gamma.mse = object$gamma.mse,
              diagnostics = if (isTRUE(candidates)) object$diagnostics else NULL,
              balance = balance,
              nn = c(control = sum(object$treat == 0L),
                     treated = sum(object$treat == 1L)))
  class(out) <- "summary.psave"
  out
}

#' @rdname summary.psave
#' @export
print.summary.psave <- function(x, digits = 3, ...) {
  cat("Summary of a psave fit\n")
  cat(sprintf("Call: %s\n\n", deparse1(x$call)))
  cat(sprintf("Estimand: %s;  criterion: %s (%s)%s\n",
              x$estimand, x$criterion,
              .criterion_label(x$criterion, x$prog.target),
              if (isTRUE(x$average)) "" else " [average = FALSE]"))
  cat(sprintf("Sample: %d treated, %d control\n\n",
              x$nn["treated"], x$nn["control"]))

  cat("Selected mixing weights:\n")
  cat("  lambda (PS):\n")
  print(round(x$lambda, digits))
  if (!is.null(x$gamma)) {
    cat("  gamma (prognostic):\n")
    print(round(x$gamma, digits))
  }
  if (!is.null(x$gamma.mse)) {
    cat("  untreated-set MSE of prognostic candidates:\n")
    print(signif(x$gamma.mse, digits))
  }
  cat(sprintf("\nCriterion value at selected lambda: %s\n",
              format(x$criterion.value, digits = digits)))

  if (!is.null(x$diagnostics)) {
    cat("\nAll criteria, per candidate and for the selected average:\n")
    print(round(x$diagnostics, digits))
    if (identical(x$criterion, "prog") && is.character(x$prog.target) &&
        !is.na(x$prog.target) && !identical(x$prog.target, "average")) {
      cat(sprintf("Note: criterion used ASMD of prognostic candidate '%s'; the 'prog' column of the diagnostics refers to the model-averaged prognostic score.\n",
                  x$prog.target))
    }
  }

  cat("\nBalance (covariates + prognostic score):\n")
  b <- round(x$balance, digits)
  b$` ` <- ifelse(x$balance$smd.wt > 0.1, "*", "")
  print(b)
  cat("---\n'*' = weighted SMD > 0.1\n")
  invisible(x)
}

#' Extract the averaged propensity or prognostic score
#'
#' `fitted()` is the canonical extractor for the model-averaged scores of a
#' [psave()] fit; `weights()` extracts the inverse-probability weights implied
#' by the averaged propensity score at the fitted estimand.
#'
#' @param object A `psave` object.
#' @param type `"ps"` (default) for the model-averaged propensity score, or
#'   `"prog"` for the model-averaged prognostic score.
#' @param ... Ignored.
#'
#' @return For `fitted()`, a numeric vector named by the rownames of the
#'   analyzed data. For `weights()`, the numeric vector `object$weights`
#'   (weights at the *fitted* estimand only; for other estimands use
#'   `WeightIt::get_w_from_ps(fitted(object), object$treat, estimand = ...)`).
#' @seealso [psave()], [predict.psave()]
#' @export
fitted.psave <- function(object, type = c("ps", "prog"), ...) {
  type <- match.arg(type)
  if (type == "prog" && is.null(object$prog)) {
    stop(paste0("No prognostic score is stored: psave() was run with criterion = \"logloss\" ",
                "and no `outcome`."), call. = FALSE)
  }
  object[[type]]
}

#' @rdname fitted.psave
#' @export
weights.psave <- function(object, ...) {
  object$weights
}

#' Predict averaged scores for new data
#'
#' Computes the model-averaged propensity score (or prognostic score) for new
#' observations by applying the stored candidate fits to `newdata`, clipping
#' candidate propensity scores as at fit time, and combining them with the
#' selected mixing weights. Requires `psave(..., keep.fits = TRUE)`.
#'
#' @param object A `psave` object fitted with `keep.fits = TRUE`.
#' @param newdata A data frame containing the variables of the propensity
#'   score formula (for `type = "ps"`) or of the prognostic specification (for
#'   `type = "prog"`). Missing values are an error.
#' @param type `"ps"` (default) or `"prog"`.
#' @param ... Ignored.
#'
#' @return A numeric vector with one score per row of `newdata`, named by its
#'   rownames. If `newdata` is missing, the in-sample fitted scores are
#'   returned (equivalent to [fitted.psave()]).
#' @seealso [psave()], [fitted.psave()]
#' @export
predict.psave <- function(object, newdata, type = c("ps", "prog"), ...) {
  type <- match.arg(type)
  if (missing(newdata) || is.null(newdata)) {
    return(fitted.psave(object, type = type))
  }
  if (is.null(object$fits)) {
    stop(paste0("psave() was called with keep.fits = FALSE, so the fitted learners were ",
                "not retained. Re-run psave() with keep.fits = TRUE to enable predict()."),
         call. = FALSE)
  }
  recs <- object$fits[[type]]
  if (is.null(recs)) {
    if (type == "ps") {
      stop(paste0("Candidate propensity scores were supplied via `ps.matrix`; there are no ",
                  "fitted learners to predict from."), call. = FALSE)
    }
    stop(paste0("No fitted prognostic learners are stored (no `outcome` was used, or the ",
                "candidates came from `prog.matrix`)."), call. = FALSE)
  }
  des <- object$fits$design[[type]]
  newdata <- as.data.frame(newdata)
  bd <- .build_design(des$terms, newdata, xlev = des$xlev,
                      what = sprintf("`newdata` (%s design)", type))
  P <- vapply(recs, function(r) .predict_learner(r, bd$Xnum, bd$Xdf),
              numeric(nrow(newdata)))
  if (!is.matrix(P)) {
    P <- matrix(P, nrow = 1L, dimnames = list(NULL, names(recs)))
  }
  if (type == "ps") {
    clip <- object$info$clip
    P <- pmin(pmax(P, clip[1L]), clip[2L])
    coefs <- object$lambda
  } else {
    coefs <- object$gamma
  }
  ## candidates appended via ps.append/prog.append are unit-specific
  ## user-supplied scores with no fitted model: predictable for new data only
  ## if the selected weights ignore them
  extra <- setdiff(names(coefs), names(recs))
  if (length(extra)) {
    if (any(coefs[extra] > 0)) {
      stop(sprintf(paste0("Cannot predict for new data: the selected %s puts weight on ",
                          "candidate(s) %s appended via `%s`, which are unit-specific ",
                          "user-supplied scores with no fitted model."),
                   if (type == "ps") "lambda" else "gamma",
                   paste0('"', extra, '"', collapse = ", "),
                   if (type == "ps") "ps.append" else "prog.append"),
           call. = FALSE)
    }
    coefs <- coefs[names(recs)]
  }
  out <- as.numeric(P %*% coefs)
  names(out) <- rownames(newdata)
  out
}
