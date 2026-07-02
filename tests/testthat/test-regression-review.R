# Regression tests requested by an independent statistical review of the
# package: (1) KS lambda-grid selection vs. a
# from-scratch weighted-eCDF loop; (2) chunk-boundary integrity of the
# 256-row chunked lambda evaluation (ATT and ATE); (3) direct-engine
# probability orientation P(A = 1 | X); (4) untreated-only prognostic
# fitting; (5) rank-deficient full-dummy GLM prediction; (6) documented
# `prog.target != "average"` behavior of criterion.value vs. diagnostics.
# Everything is deterministic (set.seed) and base-R/psave_criteria checked.

# --------------------------------------------------------------------------
# (1) KS criterion: full lambda-grid path vs. a naive weighted-eCDF KS loop
# --------------------------------------------------------------------------

test_that("KS lambda-grid selection matches a from-scratch weighted-eCDF KS loop", {
  set.seed(101)
  n <- 80
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  A <- rbinom(n, 1L, plogis(0.6 * x1 - 0.4 * x2))
  A[1:2] <- 1L; A[3:4] <- 0L                     # both arms guaranteed
  y <- 1 + x1 + rnorm(n)
  E <- cbind(m1 = pmin(pmax(plogis(0.8 * x1), 0.02), 0.98),
             m2 = pmin(pmax(plogis(0.2 - 0.5 * x2), 0.02), 0.98))
  d <- data.frame(A = A, y = y, x1 = x1, x2 = x2)

  fit <- psave(A ~ x1 + x2, data = d, outcome = ~ y,
               ps.matrix = E, prog.matrix = cbind(g = 1 + x1),
               estimand = "ATT", criterion = "ks", step = 0.2,
               keep.path = TRUE)

  # naive re-derivation: proper weighted-eCDF KS, written from scratch here
  covs <- cbind(x1 = x1, x2 = x2)
  naive_ks_mean <- function(ebar) {
    W <- ifelse(A == 1L, 1, ebar / (1 - ebar))    # ATT weights
    ks_one <- function(x) {
      xs <- sort(unique(x))
      F1 <- vapply(xs, function(z)
        sum(W[A == 1L] * (x[A == 1L] <= z)) / sum(W[A == 1L]), numeric(1L))
      F0 <- vapply(xs, function(z)
        sum(W[A == 0L] * (x[A == 0L] <= z)) / sum(W[A == 0L]), numeric(1L))
      max(abs(F1 - F0))
    }
    mean(apply(covs, 2L, ks_one))
  }

  grid <- simplex_grid(2, step = 0.2)             # 6 rows, package order
  expect_identical(nrow(fit$path), nrow(grid))
  naive <- apply(grid, 1L, function(l) naive_ks_mean(as.numeric(E %*% l)))

  expect_equal(fit$path$value, naive, tolerance = 1e-12)
  expect_identical(which(fit$path$selected), which.min(naive))
  expect_equal(unname(fit$lambda), grid[which.min(naive), ],
               tolerance = 1e-12)
  expect_equal(fit$criterion.value, min(naive), tolerance = 1e-12)
})

# --------------------------------------------------------------------------
# (2) chunk-boundary integrity: 286-row grid crosses the 256-row chunk
#     boundary with a 30-row tail; every path value re-checked against
#     psave_criteria(); ATT and ATE (treated side non-constant under ATE)
# --------------------------------------------------------------------------

test_that("smd path values are exact across the 256-row chunk boundary (ATT and ATE)", {
  set.seed(202)
  n <- 120
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  x3 <- rbinom(n, 1L, 0.5)
  A <- rbinom(n, 1L, plogis(0.5 * x1 - 0.5 * x2 + 0.4 * x3))
  A[1:2] <- 1L; A[3:4] <- 0L
  y <- 1 + x1 - x2 + rnorm(n)
  E <- cbind(m1 = pmin(pmax(plogis(0.9 * x1), 0.02), 0.98),
             m2 = pmin(pmax(plogis(-0.7 * x2), 0.02), 0.98),
             m3 = pmin(pmax(plogis(0.3 + 0.5 * x3), 0.02), 0.98),
             m4 = pmin(pmax(plogis(0.2 * x1 - 0.2 * x2), 0.02), 0.98))
  d <- data.frame(A = A, y = y, x1 = x1, x2 = x2, x3 = x3)

  for (est in c("ATT", "ATE")) {
    fit <- psave(A ~ x1 + x2 + x3, data = d, outcome = ~ y,
                 ps.matrix = E, prog.matrix = cbind(g = 1 + x1 - x2),
                 estimand = est, criterion = "smd", step = 0.1,
                 keep.path = TRUE)
    # M = 4, step = 0.1: choose(13, 3) = 286 rows -> one full 256-row chunk
    # plus a 30-row tail
    expect_identical(nrow(fit$path), 286L)

    L <- as.matrix(fit$path[, colnames(E)])
    ref <- vapply(seq_len(nrow(L)), function(i) {
      ebar <- as.numeric(fit$ps.candidates %*% L[i, ])
      unname(psave_criteria(ebar, d$A, fit$covs, estimand = est)[["smd"]])
    }, numeric(1L))

    expect_lt(max(abs(fit$path$value - ref)), 1e-12)
    expect_equal(fit$criterion.value, min(ref), tolerance = 1e-12)
  }
})

# --------------------------------------------------------------------------
# (3) engine probability orientation: every direct engine (and the SL.*
#     passthrough) must predict P(A = 1 | X), i.e., higher for treated units
# --------------------------------------------------------------------------

# n = 300, ~25% treated, strong treatment-assignment signal
.orient_data <- function() {
  set.seed(303)
  n <- 300
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  x3 <- rbinom(n, 1L, 0.4)
  A <- rbinom(n, 1L, plogis(-1.6 + 1.8 * x1 + 0.8 * x3))
  data.frame(A = A, x1 = x1, x2 = x2, x3 = x3)
}

for (eng in c("glm", "rpart", "ranger", "xgboost", "SL.glm")) {
  test_that(sprintf("PS engine \"%s\" predicts P(A = 1 | X) with the correct orientation", eng), {
    pkg <- switch(eng, rpart = "rpart", ranger = "ranger",
                  xgboost = "xgboost", SL.glm = "SuperLearner", NULL)
    if (!is.null(pkg)) skip_if_not_installed(pkg)

    d <- .orient_data()
    set.seed(404)                                 # ranger/xgboost/SL CV folds
    fit <- psave(A ~ x1 + x2 + x3, data = d, ps.methods = eng,
                 criterion = "logloss", average = FALSE)
    p <- fit$ps.candidates[, eng]
    expect_gt(mean(p[d$A == 1L]), mean(p[d$A == 0L]))
  })
}

# --------------------------------------------------------------------------
# (4) prognostic candidates are fit on the UNTREATED units only
# --------------------------------------------------------------------------

test_that("prognostic glm candidate equals the untreated-only fit, not the all-rows fit", {
  set.seed(505)
  n <- 150
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  A <- rbinom(n, 1L, plogis(0.4 * x1))
  A[1:2] <- 1L; A[3:4] <- 0L
  # treated outcomes shifted by +50: an all-rows fit is grossly different
  y <- 1 + x1 + 0.5 * x2 + rnorm(n) + 50 * A
  d <- data.frame(A = A, y = y, x1 = x1, x2 = x2)

  fit <- psave(A ~ x1 + x2, data = d, outcome = ~ y,
               ps.methods = "glm", prog.methods = "glm")

  g0 <- stats::glm(y ~ x1 + x2, data = d, subset = A == 0L,
                   family = stats::gaussian())
  ref0 <- as.numeric(stats::predict(g0, newdata = d, type = "response"))
  expect_equal(unname(fit$prog.candidates[, "glm"]), ref0, tolerance = 1e-8)

  gall <- stats::glm(y ~ x1 + x2, data = d, family = stats::gaussian())
  refall <- as.numeric(stats::predict(gall, newdata = d, type = "response"))
  expect_gt(max(abs(unname(fit$prog.candidates[, "glm"]) - refall)), 5)
})

# --------------------------------------------------------------------------
# (5) rank-deficient full-dummy GLM design: aliased coefficients handled as
#     the pivoted least-squares solution (predictions, not coefficients)
# --------------------------------------------------------------------------

test_that("rank-deficient full-dummy GLM candidate equals the pivoted-design prediction", {
  set.seed(606)
  n <- 120
  x1 <- rnorm(n)
  x2 <- 2 * x1                                    # perfectly collinear -> aliased
  f <- factor(sample(c("a", "b", "c"), n, replace = TRUE))
  A <- rbinom(n, 1L, plogis(0.8 * x1 + 0.4 * (f == "b")))
  A[1:2] <- 1L; A[3:4] <- 0L
  d <- data.frame(A = A, x1 = x1, x2 = x2, f = f)

  fit <- psave(A ~ x1 + x2 + f, data = d, ps.methods = "glm",
               criterion = "logloss", keep.fits = TRUE)

  # the package design (no intercept + full dummies + glm's own intercept)
  # is rank deficient by construction: the fit must carry aliased (NA) coefs
  expect_true(anyNA(stats::coef(fit$fits$ps$glm$fit)))

  # a pivoted (full-rank) design spanning the same column space:
  # {1, x1, f} == {1, x1, x2, fa, fb, fc} as a linear space
  gref <- stats::glm(A ~ x1 + f, data = d, family = stats::binomial())
  ref <- as.numeric(stats::predict(gref, type = "response"))
  ref <- pmin(pmax(ref, 0.01), 0.99)              # same clipping as candidates
  expect_equal(unname(fit$ps.candidates[, "glm"]), ref, tolerance = 1e-8)
})

# --------------------------------------------------------------------------
# (6) prog.target = named candidate: criterion.value targets that candidate;
#     diagnostics$prog documents the gamma-averaged score gbar
# --------------------------------------------------------------------------

test_that("prog.target names a candidate: criterion.value vs gbar-based diagnostics", {
  fx <- make_fixture()
  # G2: candidate "g2" equals y exactly on the untreated units -> gamma = (0, 1)
  expect_equal(unname(fx$G2[fx$data$A == 0, "g2"]),
               fx$data$y[fx$data$A == 0])

  fit <- psave(A ~ x, data = fx$data, outcome = ~ y,
               ps.matrix = fx$E, prog.matrix = fx$G2,
               estimand = "ATT", criterion = "prog", prog.target = "g1",
               step = 0.5)
  expect_identical(fit$prog.target, "g1")
  expect_equal(fit$gamma, c(g1 = 0, g2 = 1))
  # g1 is control-constant -> criterion identically 0 -> first grid row wins
  expect_equal(fit$lambda, c(e1 = 1, e2 = 0))

  covs <- cbind(x = fx$data$x)
  # criterion.value == psave_criteria at the TARGET candidate (column 1, "g1")
  crit.target <- psave_criteria(fit$ps, fx$data$A, covs,
                                prog = fx$G2[, 1L],
                                estimand = "ATT")[["prog"]]
  expect_equal(fit$criterion.value, unname(crit.target), tolerance = 1e-12)

  # documented behavior: diagnostics always uses the gamma-averaged gbar
  crit.gbar <- psave_criteria(fit$ps, fx$data$A, covs, prog = fit$prog,
                              estimand = "ATT")[["prog"]]
  expect_equal(fit$diagnostics["average", "prog"], unname(crit.gbar),
               tolerance = 1e-12)
  # and here the two deliberately differ (0 vs 31/21)
  expect_gt(abs(fit$diagnostics["average", "prog"] - fit$criterion.value), 0.1)

  # the summary print carries the explanatory footnote
  expect_output(print(summary(fit)),
                "criterion used ASMD of prognostic candidate 'g1'")
})
