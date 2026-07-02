# D.14 (v): IPW weights — exact supplement formulas and equality with
# WeightIt::get_w_from_ps().

test_that("ATT weights follow the supplement formula", {
  fx <- make_fixture()
  fit <- psave(A ~ x, data = fx$data, outcome = ~ y,
               ps.matrix = fx$E, prog.matrix = fx$G1,
               estimand = "ATT", step = 0.5)
  A <- fx$data$A
  expected <- ifelse(A == 1, 1, fit$ps / (1 - fit$ps))
  expect_equal(unname(fit$weights), unname(expected), tolerance = 1e-12)
  # weights() extractor returns the fitted-estimand weights
  expect_equal(unname(weights(fit)), unname(fit$weights))
  # treated units carry weight exactly 1 under ATT
  expect_true(all(fit$weights[A == 1] == 1))
})

test_that("ATE weights follow the supplement formula", {
  fx <- make_fixture()
  fit <- psave(A ~ x, data = fx$data, outcome = ~ y,
               ps.matrix = fx$E, prog.matrix = fx$G1,
               estimand = "ATE", step = 0.5)
  A <- fx$data$A
  expected <- ifelse(A == 1, 1 / fit$ps, 1 / (1 - fit$ps))
  expect_equal(unname(fit$weights), unname(expected), tolerance = 1e-12)
})

test_that("weights() equals WeightIt::get_w_from_ps for ATT and ATE", {
  skip_if_not_installed("WeightIt")
  fx <- make_fixture()
  for (est in c("ATT", "ATE")) {
    fit <- psave(A ~ x, data = fx$data, outcome = ~ y,
                 ps.matrix = fx$E, prog.matrix = fx$G1,
                 estimand = est, step = 0.5)
    expected <- WeightIt::get_w_from_ps(fit$ps, treat = fx$data$A,
                                        estimand = est)
    expect_equal(unname(weights(fit)), unname(as.numeric(expected)),
                 tolerance = 1e-12)
  }
})

test_that("weights are consistent with the stored averaged ps", {
  # weights must be a deterministic function of fit$ps (no separate re-clip)
  fx <- make_fixture()
  Eclip <- cbind(p1 = c(0.001, 0.5, 0.999, 0.2, 0.5, 0.8), p2 = fx$E[, "e2"])
  fit <- psave(A ~ x, data = fx$data,
               ps.matrix = Eclip, criterion = "logloss", step = 0.5)
  A <- fx$data$A
  expected <- ifelse(A == 1, 1, fit$ps / (1 - fit$ps))
  expect_equal(unname(fit$weights), unname(expected), tolerance = 1e-12)
})
