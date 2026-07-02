# input.R -- input validation, treatment coercion, covariate-matrix
# construction (D.1, D.2). Nothing here is exported.

# Coerce a treatment vector to integer 0/1 per MatchIt convention (D.1):
# numeric 0/1 kept; logical -> 0/1; two-level factor/character -> the SECOND
# level (in factor-level order) is treated = 1; anything else errors.
.coerce_treat <- function(a, name = "treatment") {
  if (is.numeric(a)) {
    u <- unique(a[!is.na(a)])
    if (!all(u %in% c(0, 1))) {
      stop(sprintf("The %s variable must be binary: numeric 0/1, logical, or a two-level factor/character.",
                   name), call. = FALSE)
    }
    out <- as.integer(a)
  } else if (is.logical(a)) {
    out <- as.integer(a)
  } else if (is.factor(a) || is.character(a)) {
    f <- droplevels(as.factor(a))
    if (nlevels(f) != 2L) {
      stop(sprintf("The %s variable must have exactly two levels; found %d (%s).",
                   name, nlevels(f), paste0('"', levels(f), '"', collapse = ", ")),
           call. = FALSE)
    }
    out <- as.integer(f == levels(f)[2L])
  } else {
    stop(sprintf("The %s variable must be binary: numeric 0/1, logical, or a two-level factor/character.",
                 name), call. = FALSE)
  }
  out
}

# Error (never drop) on missing values (D.1).
.check_no_na <- function(data, vars, what) {
  vars <- vars[vars %in% names(data)]
  if (!length(vars)) return(invisible(TRUE))
  bad <- vars[vapply(vars, function(v) anyNA(data[[v]]), logical(1L))]
  if (length(bad)) {
    stop(sprintf(paste0("Missing values found in %s variable(s): %s.\n",
                        "psAve requires complete cases in all used variables and never drops rows silently; ",
                        "handle missing data (e.g., by imputation) before calling psave()."),
                 what, paste0('"', bad, '"', collapse = ", ")),
         call. = FALSE)
  }
  invisible(TRUE)
}

# Build the covariate designs from a (response-free) terms object (D.2):
#  - covs: numeric model matrix with NO intercept and FULL dummy expansion
#    (every factor level gets a column; cobalt splitfactor convention), used
#    for the smd/ks criteria and the balance table;
#  - Xnum: the same matrix as a data.frame with syntactic names (for glm,
#    xgboost, and SL.* learners);
#  - Xdf:  the raw model frame of the RHS variables with syntactic names (for
#    rpart and ranger) -- identical rows;
#  - bin.vars: TRUE iff a covs column takes exactly two distinct values;
#  - xlev: factor levels, stored for predict().
.build_design <- function(tt, data, xlev = NULL, what = "formula") {
  mf <- stats::model.frame(tt, data = data, na.action = stats::na.pass,
                           xlev = xlev, drop.unused.levels = is.null(xlev))
  if (anyNA(mf)) {
    stop(sprintf(paste0("Missing values found among the %s variables. ",
                        "psAve requires complete cases and never drops rows silently."),
                 what), call. = FALSE)
  }
  for (j in seq_along(mf)) {
    if (is.character(mf[[j]])) mf[[j]] <- factor(mf[[j]])
  }
  is.fac <- vapply(mf, is.factor, logical(1L))
  contr <- lapply(mf[is.fac], stats::contrasts, contrasts = FALSE)
  X <- stats::model.matrix(tt, data = mf,
                           contrasts.arg = if (length(contr)) contr else NULL)
  X <- X[, setdiff(colnames(X), "(Intercept)"), drop = FALSE]
  if (ncol(X) == 0L) {
    stop("The model formula must contain at least one covariate on the right-hand side.",
         call. = FALSE)
  }
  storage.mode(X) <- "double"

  Xnum <- as.data.frame(X)
  names(Xnum) <- make.names(colnames(X), unique = TRUE)

  Xdf <- as.data.frame(mf)
  names(Xdf) <- make.names(names(Xdf), unique = TRUE)

  list(covs = X,
       Xnum = Xnum,
       Xdf = Xdf,
       bin.vars = .detect_bin(X),
       xlev = stats::.getXlevels(tt, mf))
}

# Validate a user-supplied ps.matrix / prog.matrix (B.1, D.3).
.validate_score_matrix <- function(m, n, arg, rn = NULL, prob = FALSE) {
  if (!is.matrix(m) || !is.numeric(m)) {
    stop(sprintf("`%s` must be a numeric matrix.", arg), call. = FALSE)
  }
  if (nrow(m) != n) {
    stop(sprintf("`%s` must have one row per unit of `data` (%d rows; found %d).",
                 arg, n, nrow(m)), call. = FALSE)
  }
  cn <- colnames(m)
  if (is.null(cn) || anyNA(cn) || any(!nzchar(cn)) || anyDuplicated(cn)) {
    stop(sprintf("`%s` must have unique, non-empty column names (they label the candidates).",
                 arg), call. = FALSE)
  }
  if (anyNA(m) || any(!is.finite(m))) {
    stop(sprintf("`%s` contains missing or non-finite values.", arg), call. = FALSE)
  }
  if (prob && (any(m <= 0) || any(m >= 1))) {
    stop(sprintf("All values of `%s` must lie strictly inside (0, 1).", arg),
         call. = FALSE)
  }
  if (!is.null(rn) && !is.null(rownames(m)) && !identical(rownames(m), rn)) {
    stop(sprintf(paste0("The rownames of `%s` do not match rownames(data). ",
                        "Row alignment is essential; supply the matrix in the row order of `data` ",
                        "(or drop its rownames if alignment is already guaranteed)."), arg),
         call. = FALSE)
  }
  m
}

# Normalize a `ps.append` / `prog.append` argument (B.1, D.3): a numeric
# vector of length n becomes ONE candidate column labeled "append"; a numeric
# matrix or all-numeric data.frame with n rows supplies several candidates
# (unique, non-empty column names required, exactly as for `ps.matrix`).
# Validation is then shared with .validate_score_matrix(). Automatic
# data.frame rownames ("1", "2", ...) are dropped so they cannot trip the
# rownames-alignment guard; explicit rownames are kept and checked.
.normalize_append <- function(x, n, arg, rn = NULL, prob = FALSE) {
  if (is.null(x)) return(NULL)
  if (is.data.frame(x)) {
    if (!all(vapply(x, is.numeric, logical(1L)))) {
      stop(sprintf("All columns of `%s` must be numeric.", arg), call. = FALSE)
    }
    auto.rn <- .row_names_info(x, type = 1L) < 0L
    x <- as.matrix(x)
    if (isTRUE(auto.rn)) rownames(x) <- NULL
  } else if (is.numeric(x) && is.null(dim(x))) {
    if (length(x) != n) {
      stop(sprintf("`%s` as a vector must have one value per unit of `data` (%d values; found %d).",
                   arg, n, length(x)), call. = FALSE)
    }
    x <- matrix(x, ncol = 1L, dimnames = list(names(x), "append"))
  }
  if (!is.matrix(x) || !is.numeric(x)) {
    stop(sprintf("`%s` must be a numeric vector of length %d, or a numeric matrix/all-numeric data.frame with %d rows.",
                 arg, n, n), call. = FALSE)
  }
  .validate_score_matrix(x, n, arg, rn = rn, prob = prob)
}

# Append validated extra candidates AFTER the base candidate columns (D.3):
# simplex_grid() enumerates the first components first, so the first-minimum
# tie-break favors the base candidates over appended ones.
.append_candidates <- function(base, extra, arg) {
  clash <- intersect(colnames(extra), colnames(base))
  if (length(clash)) {
    stop(sprintf(paste0("Column name(s) %s of `%s` collide with existing candidate labels; ",
                        "rename the appended column(s) (for a vector, pass a named ",
                        "1-column matrix to choose a label other than \"append\")."),
                 paste0('"', clash, '"', collapse = ", "), arg),
         call. = FALSE)
  }
  cbind(base, extra)
}

# Full input processing for psave() (D.1, D.2). Returns everything the
# orchestrator needs; errors early with actionable messages.
.process_inputs <- function(formula, data, outcome, criterion, family) {
  if (!inherits(formula, "formula") || length(formula) != 3L) {
    stop("`formula` must be a two-sided formula of the form `treat ~ x1 + x2 + ...`.",
         call. = FALSE)
  }
  if (missing(data) || is.null(data)) {
    stop("`data` must be supplied as a data.frame.", call. = FALSE)
  }
  data <- as.data.frame(data)
  n <- nrow(data)
  if (n == 0L) stop("`data` has zero rows.", call. = FALSE)
  rn <- rownames(data)

  ## --- treatment -----------------------------------------------------------
  tname <- deparse1(formula[[2L]])
  a <- eval(formula[[2L]], data, environment(formula))
  if (length(a) != n) {
    stop(sprintf("The treatment variable `%s` has length %d; expected %d (one per row of `data`).",
                 tname, length(a), n), call. = FALSE)
  }
  if (anyNA(a)) {
    stop(sprintf(paste0("Missing values found in the treatment variable `%s`. ",
                        "psAve requires complete cases and never drops rows silently."),
                 tname), call. = FALSE)
  }
  treat <- .coerce_treat(a, name = sprintf("treatment (`%s`)", tname))
  if (sum(treat == 1L) < 2L || sum(treat == 0L) < 2L) {
    stop(sprintf("At least 2 units are required in each treatment arm; found %d treated and %d control.",
                 sum(treat == 1L), sum(treat == 0L)), call. = FALSE)
  }

  ## --- covariates ----------------------------------------------------------
  ps.tt <- stats::delete.response(stats::terms(formula, data = data))
  if (length(attr(ps.tt, "term.labels")) == 0L) {
    stop("`formula` must contain at least one covariate on the right-hand side.",
         call. = FALSE)
  }
  .check_no_na(data, all.vars(ps.tt), "covariate")
  ps.design <- .build_design(ps.tt, data, what = "covariate (formula right-hand side)")

  ## --- outcome / prognostic design -----------------------------------------
  y <- NULL
  outcome.name <- NA_character_
  prog.design <- NULL
  prog.tt <- NULL
  prog.same <- FALSE
  if (is.null(outcome)) {
    if (criterion == "prog") {
      stop(paste0("criterion = 'prog' requires the outcome variable. ",
                  "Prognostic models are fit on UNTREATED units only and the criterion ",
                  "never uses a treated-untreated outcome contrast, so this does not bias ",
                  "effect estimation (Hansen 2008); see vignette('method-details', 'psAve')."),
           call. = FALSE)
    }
  } else {
    if (!inherits(outcome, "formula")) {
      stop("`outcome` must be a formula: one-sided `~ y` (reusing the formula RHS as prognostic predictors) or two-sided `y ~ z1 + z2`.",
           call. = FALSE)
    }
    if (length(outcome) == 2L) {
      ## one-sided: ~ y -- reuse the PS covariates as prognostic predictors
      yexpr <- outcome[[2L]]
      outcome.name <- deparse1(yexpr)
      .check_no_na(data, all.vars(yexpr), "outcome")
      y <- eval(yexpr, data, environment(outcome))
      prog.tt <- ps.tt
      prog.design <- ps.design
      prog.same <- TRUE
    } else {
      ## two-sided: y ~ z1 + z2 -- distinct prognostic specification
      yexpr <- outcome[[2L]]
      outcome.name <- deparse1(yexpr)
      .check_no_na(data, all.vars(yexpr), "outcome")
      y <- eval(yexpr, data, environment(outcome))
      prog.tt <- stats::delete.response(stats::terms(outcome, data = data))
      if (length(attr(prog.tt, "term.labels")) == 0L) {
        stop("The two-sided `outcome` formula must contain at least one prognostic predictor on its right-hand side.",
             call. = FALSE)
      }
      .check_no_na(data, all.vars(prog.tt), "prognostic predictor")
      if (any(all.vars(yexpr) %in% all.vars(prog.tt))) {
        stop("The outcome variable must not appear among the prognostic predictors in `outcome`.",
             call. = FALSE)
      }
      prog.design <- .build_design(prog.tt, data,
                                   what = "prognostic predictor (outcome right-hand side)")
      prog.same <- FALSE
    }
    if (length(y) != n) {
      stop(sprintf("The outcome variable `%s` has length %d; expected %d (one per row of `data`).",
                   outcome.name, length(y), n), call. = FALSE)
    }
    if (anyNA(y)) {
      stop(sprintf(paste0("Missing values found in the outcome variable `%s`. ",
                          "psAve requires complete cases and never drops rows silently."),
                   outcome.name), call. = FALSE)
    }
    ## outcome leakage guard: Y must not be a PS covariate
    if (any(all.vars(yexpr) %in% all.vars(ps.tt))) {
      stop(sprintf("The outcome variable `%s` must not appear among the propensity score covariates in `formula`.",
                   outcome.name), call. = FALSE)
    }
    ## coerce y per family
    if (identical(family$family, "binomial")) {
      y <- as.numeric(.coerce_treat(y, name = sprintf("binary outcome (`%s`)", outcome.name)))
    } else {
      if (!is.numeric(y)) {
        stop(sprintf("With family = gaussian(), the outcome variable `%s` must be numeric.",
                     outcome.name), call. = FALSE)
      }
      y <- as.numeric(y)
    }
  }

  list(data = data, n = n, rn = rn,
       treat = treat, treat.name = tname,
       covs = ps.design$covs, bin.vars = ps.design$bin.vars,
       Xnum = ps.design$Xnum, Xdf = ps.design$Xdf,
       ps.tt = ps.tt, ps.xlev = ps.design$xlev,
       y = y, outcome.name = outcome.name,
       prog.tt = prog.tt,
       prog.Xnum = if (!is.null(prog.design)) prog.design$Xnum else NULL,
       prog.Xdf = if (!is.null(prog.design)) prog.design$Xdf else NULL,
       prog.xlev = if (!is.null(prog.design)) prog.design$xlev else NULL,
       prog.same = prog.same)
}
