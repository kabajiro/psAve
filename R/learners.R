# learners.R -- the engine registry (D.3, D.6): pinned direct engines
# ("glm", "rpart", "ranger", "xgboost"), the "SL.*" SuperLearner passthrough,
# per-learner `control =` merging, and the matching predict dispatch used by
# predict.psave(). Nothing here is exported.
#
# Conventions:
#  - task = "ps":   binary response (the treatment), family binomial();
#                   fit on ALL n units; in-sample predictions for all n.
#  - task = "prog": outcome model with the user's `family`; fit on the
#                   UNTREATED units only (fit.idx); predictions for ALL n.
#  - Hyperparameters are pinned to the SuperLearner-wrapper defaults used in
#    the paper, so the direct and SL.* routes agree; `control =` entries
#    override them per learner, and the resolved values are stored in
#    `info$learners`.
#  - No scale() anywhere: engines receive raw covariates (fixes the reference
#    code's train/test-inconsistent scaling).

.known_engines <- c("glm", "rpart", "ranger", "xgboost")

.engine_of <- function(label) {
  if (grepl("^SL\\.", label)) "SuperLearner" else label
}

.validate_methods <- function(labels, arg) {
  if (!is.character(labels) || length(labels) == 0L || anyNA(labels) ||
      any(!nzchar(labels))) {
    stop(sprintf("`%s` must be a non-empty character vector of learner labels.", arg),
         call. = FALSE)
  }
  if (anyDuplicated(labels)) {
    stop(sprintf("`%s` contains duplicated labels: %s.", arg,
                 paste0('"', unique(labels[duplicated(labels)]), '"', collapse = ", ")),
         call. = FALSE)
  }
  bad <- labels[!(labels %in% .known_engines | grepl("^SL\\.", labels))]
  if (length(bad)) {
    stop(sprintf(paste0("Unknown learner label(s) in `%s`: %s.\n",
                        "Available: %s, or any \"SL.*\" wrapper name (requires the SuperLearner package)."),
                 arg, paste0('"', bad, '"', collapse = ", "),
                 paste0('"', .known_engines, '"', collapse = ", ")),
         call. = FALSE)
  }
  labels
}

# ---------------------------------------------------------------------------
# glm --------------------------------------------------------------------

# Manual response prediction treating aliased (NA) coefficients as 0. This is
# identical to the fitted values of the pivoted least-squares solution glm
# uses internally, but avoids the rank-deficiency warning that predict.lm()
# emits for the full-dummy-expansion design (which is rank-deficient by
# construction whenever a factor covariate is present).
.predict_glm_manual <- function(fit, Xnum) {
  b <- stats::coef(fit)
  b[is.na(b)] <- 0
  X <- cbind("(Intercept)" = 1, as.matrix(Xnum))
  eta <- drop(X[, names(b), drop = FALSE] %*% b)
  stats::family(fit)$linkinv(eta)
}

.fit_glm <- function(y, Xnum, family, task, control, fit.idx) {
  fam <- if (task == "ps") stats::binomial() else family
  dat <- Xnum
  dat[[".psave_y"]] <- y
  args <- list(formula = .psave_y ~ ., data = dat[fit.idx, , drop = FALSE],
               family = fam)
  if (length(control)) args <- utils::modifyList(args, control)
  fit <- do.call(stats::glm, args)
  pred <- .predict_glm_manual(fit, Xnum)
  params <- args[setdiff(names(args), c("formula", "data"))]
  params$family <- paste0(fam$family, "(link = \"", fam$link, "\")")
  list(engine = "glm", fit = fit, pred = as.numeric(pred), params = params,
       package = "stats", version = as.character(getRversion()))
}

# ---------------------------------------------------------------------------
# rpart -------------------------------------------------------------------

.fit_rpart <- function(y, Xdf, family, task, control, fit.idx, label) {
  .require_pkg("rpart", sprintf("for the \"%s\" learner", label))
  binary <- task == "ps" || identical(family$family, "binomial")
  ctrl <- list(cp = 0.01, minsplit = 20, maxdepth = 30, xval = 0)
  if (length(control)) ctrl <- utils::modifyList(ctrl, control)
  dat <- Xdf
  dat[[".psave_y"]] <- if (binary) factor(y, levels = c(0, 1)) else y
  fit <- rpart::rpart(.psave_y ~ ., data = dat[fit.idx, , drop = FALSE],
                      method = if (binary) "class" else "anova",
                      control = do.call(rpart::rpart.control, ctrl))
  pred <- if (binary) stats::predict(fit, newdata = Xdf, type = "prob")[, "1"]
          else as.numeric(stats::predict(fit, newdata = Xdf, type = "vector"))
  params <- c(list(method = if (binary) "class" else "anova"), ctrl)
  list(engine = "rpart", fit = fit, pred = as.numeric(pred), params = params,
       binary = binary, package = "rpart",
       version = as.character(utils::packageVersion("rpart")))
}

.predict_rpart <- function(rec, Xdf) {
  if (rec$binary) stats::predict(rec$fit, newdata = Xdf, type = "prob")[, "1"]
  else as.numeric(stats::predict(rec$fit, newdata = Xdf, type = "vector"))
}

# ---------------------------------------------------------------------------
# ranger ------------------------------------------------------------------

.fit_ranger <- function(y, Xdf, family, task, control, fit.idx, label) {
  .require_pkg("ranger", sprintf("for the \"%s\" learner", label))
  binary <- task == "ps" || identical(family$family, "binomial")
  args <- list(num.trees = 500,
               mtry = max(1L, floor(sqrt(ncol(Xdf)))),
               min.node.size = if (binary) 1 else 5,
               probability = binary,
               # single-threaded by default (CRAN policy: at most 2 cores);
               # raise via control = list(ranger = list(num.threads = ...))
               num.threads = 1L)
  if (length(control)) args <- utils::modifyList(args, control)
  yfit <- if (binary) factor(y, levels = c(0, 1)) else y
  fit <- do.call(ranger::ranger,
                 c(list(x = Xdf[fit.idx, , drop = FALSE], y = yfit[fit.idx]), args))
  pred <- .predict_ranger_raw(fit, Xdf)
  list(engine = "ranger", fit = fit, pred = as.numeric(pred), params = args,
       binary = binary, package = "ranger",
       version = as.character(utils::packageVersion("ranger")))
}

.predict_ranger_raw <- function(fit, Xdf) {
  pr <- stats::predict(fit, data = Xdf)$predictions
  if (is.matrix(pr)) pr[, "1"]
  else if (is.factor(pr)) as.numeric(as.character(pr))
  else as.numeric(pr)
}

.predict_ranger <- function(rec, Xdf) {
  .predict_ranger_raw(rec$fit, Xdf)
}

# ---------------------------------------------------------------------------
# xgboost -----------------------------------------------------------------

# Pinned to the SL.xgboost defaults (= the paper's GBM): nrounds = 1000,
# max_depth = 4, eta (shrinkage) = 0.1, min_child_weight = 10.
.fit_xgboost <- function(y, Xnum, family, task, control, fit.idx, label) {
  .require_pkg("xgboost", sprintf("for the \"%s\" learner", label))
  binary <- task == "ps" || identical(family$family, "binomial")
  params <- list(objective = if (binary) "binary:logistic" else "reg:squarederror",
                 max_depth = 4, eta = 0.1, min_child_weight = 10,
                 # single-threaded by default (CRAN policy: at most 2 cores);
                 # raise via control = list(xgboost = list(nthread = ...))
                 nthread = 1L)
  nrounds <- 1000
  vrb <- 0
  if (length(control)) {
    if (!is.null(control$nrounds)) { nrounds <- control$nrounds; control$nrounds <- NULL }
    if (!is.null(control$verbose)) { vrb <- control$verbose; control$verbose <- NULL }
    params <- utils::modifyList(params, control)
  }
  X <- as.matrix(Xnum)
  dtrain <- xgboost::xgb.DMatrix(X[fit.idx, , drop = FALSE], label = y[fit.idx])
  fit <- xgboost::xgb.train(params = params, data = dtrain, nrounds = nrounds,
                            verbose = vrb)
  pred <- stats::predict(fit, xgboost::xgb.DMatrix(X))
  list(engine = "xgboost", fit = fit, pred = as.numeric(pred),
       params = c(params, list(nrounds = nrounds)), binary = binary,
       package = "xgboost",
       version = as.character(utils::packageVersion("xgboost")))
}

.predict_xgboost <- function(rec, Xnum) {
  .require_pkg("xgboost", "to predict from a fitted xgboost learner")
  as.numeric(stats::predict(rec$fit, xgboost::xgb.DMatrix(as.matrix(Xnum))))
}

# ---------------------------------------------------------------------------
# SuperLearner passthrough --------------------------------------------------

# Any "SL.*" label is passed verbatim to SuperLearner::SuperLearner() as a
# single-element SL.library; the candidate prediction is that learner's
# full-data-refit column of `library.predict` (the exact-replication path for
# the paper). newX = all n rows, so prognostic candidates fit on untreated
# units still predict for everyone.
.fit_sl <- function(label, y, Xnum, family, task, control, cv, fit.idx) {
  .require_pkg("SuperLearner", sprintf("for the \"%s\" learner", label))
  fam <- if (task == "ps") stats::binomial() else family
  if (length(control)) {
    warning(sprintf(paste0("`control` entries for \"%s\" are ignored: customize SL.* learners ",
                           "by writing a custom wrapper (see SuperLearner::SL.template)."),
                    label), call. = FALSE)
  }
  Xfit <- Xnum[fit.idx, , drop = FALSE]
  Yfit <- y[fit.idx]
  ## env: SuperLearner() looks its learner/screen functions up in `env`;
  ## the package namespace finds the built-in SL.* wrappers even when
  ## SuperLearner is not attached, and its enclosure chain still reaches the
  ## user's global environment for custom wrappers.
  fit <- SuperLearner::SuperLearner(Y = Yfit, X = Xfit, newX = Xnum,
                                    family = fam, SL.library = label,
                                    cvControl = list(V = cv),
                                    env = asNamespace("SuperLearner"))
  pred <- as.numeric(fit$library.predict[, 1L])
  list(engine = "SuperLearner", fit = fit, pred = pred,
       params = list(SL.library = label, cvControl = list(V = cv),
                     family = fam$family),
       package = "SuperLearner",
       version = as.character(utils::packageVersion("SuperLearner")),
       sl.X = Xfit, sl.Y = Yfit)
}

.predict_sl <- function(rec, Xnum) {
  .require_pkg("SuperLearner", "to predict from a fitted SL.* learner")
  pr <- stats::predict(rec$fit, newdata = Xnum, X = rec$sl.X, Y = rec$sl.Y,
                       onlySL = FALSE)
  as.numeric(pr$library.predict[, 1L])
}

# ---------------------------------------------------------------------------
# dispatchers ----------------------------------------------------------------

.fit_learner <- function(label, y, Xnum, Xdf, family, task, control = NULL,
                         cv = 5L, fit.idx = NULL) {
  n <- nrow(Xnum)
  if (is.null(fit.idx)) fit.idx <- seq_len(n)
  control <- if (is.null(control)) list() else control
  rec <- switch(.engine_of(label),
                glm = .fit_glm(y, Xnum, family, task, control, fit.idx),
                rpart = .fit_rpart(y, Xdf, family, task, control, fit.idx, label),
                ranger = .fit_ranger(y, Xdf, family, task, control, fit.idx, label),
                xgboost = .fit_xgboost(y, Xnum, family, task, control, fit.idx, label),
                SuperLearner = .fit_sl(label, y, Xnum, family, task, control, cv, fit.idx),
                stop(sprintf("Unknown learner label \"%s\".", label), call. = FALSE))
  if (anyNA(rec$pred) || any(!is.finite(rec$pred))) {
    stop(sprintf("Learner \"%s\" produced missing or non-finite predictions.", label),
         call. = FALSE)
  }
  rec$label <- label
  rec$task <- task
  rec
}

# Fit all candidates for one task and collect the n x M prediction matrix.
.fit_candidates <- function(labels, y, Xnum, Xdf, family, task, control,
                            cv, fit.idx = NULL, verbose = FALSE) {
  recs <- vector("list", length(labels))
  names(recs) <- labels
  for (i in seq_along(labels)) {
    .vmsg(verbose, sprintf("Fitting %s candidate \"%s\" ...",
                           if (task == "ps") "propensity score" else "prognostic",
                           labels[i]))
    recs[[i]] <- .fit_learner(labels[i], y, Xnum, Xdf, family, task,
                              control[[labels[i]]], cv, fit.idx)
  }
  preds <- vapply(recs, function(r) r$pred, numeric(nrow(Xnum)))
  if (!is.matrix(preds)) preds <- matrix(preds, ncol = length(recs))
  colnames(preds) <- labels
  list(preds = preds, records = recs)
}

# Predict from a stored learner record on a new design (predict.psave).
.predict_learner <- function(rec, Xnum, Xdf) {
  out <- switch(rec$engine,
                glm = .predict_glm_manual(rec$fit, Xnum),
                rpart = .predict_rpart(rec, Xdf),
                ranger = .predict_ranger(rec, Xdf),
                xgboost = .predict_xgboost(rec, Xnum),
                SuperLearner = .predict_sl(rec, Xnum),
                stop(sprintf("Unknown engine \"%s\".", rec$engine), call. = FALSE))
  as.numeric(out)
}

# Summarize learner provenance for info$learners (B.2).
.learner_info <- function(recs) {
  lapply(recs, function(r) r[c("label", "engine", "package", "version", "params")])
}

# info$learners entries for candidates appended via ps.append / prog.append:
# user-supplied score columns, not fitted engines (type "append").
.append_info <- function(labels) {
  stats::setNames(lapply(labels, function(lab) {
    list(label = lab, engine = "append", package = "user-supplied",
         version = NA_character_, params = NULL, type = "append")
  }), labels)
}
