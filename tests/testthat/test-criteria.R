# D.14 (ii)-(iii): psave_criteria — exact hand-computed values on the 6-row
# fixture and from-scratch agreement on random data.  All expected numbers
# are derived in helper-fixtures.R comments.

test_that("psave_criteria returns the four named criteria", {
  fx <- make_fixture()
  res <- psave_criteria(ps = fx$E[, "e1"], treat = fx$data$A,
                        covs = cbind(x = fx$data$x), prog = fx$G1[, "g"],
                        estimand = "ATT")
  expect_named(res, c("logloss", "smd", "ks", "prog"))
  expect_true(is.numeric(res))
  # prog slot is NA when no prognostic score is supplied
  res2 <- psave_criteria(ps = fx$E[, "e1"], treat = fx$data$A,
                         covs = cbind(x = fx$data$x), estimand = "ATT")
  expect_true(is.na(res2[["prog"]]))
})

test_that("hand-computed criteria on the 6-row fixture (ATT, ps = e1)", {
  fx <- make_fixture()
  res <- psave_criteria(ps = fx$E[, "e1"], treat = fx$data$A,
                        covs = cbind(x = fx$data$x), prog = fx$G1[, "g"],
                        estimand = "ATT")
  # derivations in helper-fixtures.R:
  # logloss = -(4 log .5 + log .8 + log .2)/6; smd = 12/7; ks = 16/21;
  # prog = 31/21
  expect_equal(res[["logloss"]], -(4 * log(0.5) + log(0.8) + log(0.2)) / 6,
               tolerance = 1e-10)
  expect_equal(res[["smd"]],  12 / 7,  tolerance = 1e-10)
  expect_equal(res[["ks"]],   16 / 21, tolerance = 1e-10)
  expect_equal(res[["prog"]], 31 / 21, tolerance = 1e-10)
})

test_that("hand-computed criteria on the 6-row fixture (ATE, ps = e1)", {
  fx <- make_fixture()
  res <- psave_criteria(ps = fx$E[, "e1"], treat = fx$data$A,
                        covs = cbind(x = fx$data$x), prog = fx$G1[, "g"],
                        estimand = "ATE")
  # derivations in helper-fixtures.R:
  # logloss unchanged (weight-free); smd = 16/11; ks = 20/33; prog = 35/33
  expect_equal(res[["logloss"]], -(4 * log(0.5) + log(0.8) + log(0.2)) / 6,
               tolerance = 1e-10)
  expect_equal(res[["smd"]],  16 / 11, tolerance = 1e-10)
  expect_equal(res[["ks"]],   20 / 33, tolerance = 1e-10)
  expect_equal(res[["prog"]], 35 / 33, tolerance = 1e-10)
})

test_that("weighted ASMD matches a from-scratch computation on random data", {
  # D.14 (ii) + H.1: uniform sample-SD standardization (bin.vars = FALSE
  # convention) for continuous AND binary covariates.
  set.seed(42)
  n <- 60
  treat <- rbinom(n, 1, 0.5)
  ps <- runif(n, 0.05, 0.95)
  xc <- rnorm(n)             # continuous
  xb <- rbinom(n, 1, 0.4)    # binary
  covs <- cbind(xc = xc, xb = xb)
  for (est in c("ATT", "ATE")) {
    w <- ref_weights(ps, treat, est)
    expected <- mean(c(ref_asmd(xc, treat, w), ref_asmd(xb, treat, w)))
    res <- psave_criteria(ps = ps, treat = treat, covs = covs, estimand = est)
    expect_equal(res[["smd"]], expected, tolerance = 1e-10)
  }
})

test_that("weighted KS matches a from-scratch weighted-eCDF computation", {
  # D.14 (ii): proper weighted eCDF in each arm, sup over observed values
  set.seed(7)
  n <- 80
  treat <- rbinom(n, 1, 0.45)
  ps <- runif(n, 0.1, 0.9)
  xc <- rexp(n)
  xb <- rbinom(n, 1, 0.5)    # binary: KS = |diff in weighted proportions|
  covs <- cbind(xc = xc, xb = xb)
  for (est in c("ATT", "ATE")) {
    w <- ref_weights(ps, treat, est)
    expected <- mean(c(ref_ks(xc, treat, w), ref_ks(xb, treat, w)))
    res <- psave_criteria(ps = ps, treat = treat, covs = covs, estimand = est)
    expect_equal(res[["ks"]], expected, tolerance = 1e-10)
  }
})

test_that("logloss matches the from-scratch Bernoulli log-likelihood", {
  set.seed(42)
  n <- 60
  treat <- rbinom(n, 1, 0.5)
  ps <- runif(n, 0.05, 0.95)
  res <- psave_criteria(ps = ps, treat = treat, covs = cbind(x = rnorm(n)),
                        estimand = "ATT")
  expect_equal(res[["logloss"]], ref_logloss(ps, treat), tolerance = 1e-10)
})
