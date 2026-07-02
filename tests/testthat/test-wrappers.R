# D.14 (xii): psave_match / psave_weight round-trips.  The wrappers reuse
# the stored formula and data, and the distance / ps slot must equal fit$ps.

test_that("psave_match returns a matchit object with distance = fit$ps", {
  skip_if_not_installed("MatchIt")
  sim <- make_sim_data()
  fit <- psave(sim_formula(), data = sim, outcome = ~ y,
               ps.methods = "glm", prog.methods = "glm")
  m <- psave_match(fit)
  expect_s3_class(m, "matchit")
  expect_equal(unname(as.numeric(m$distance)), unname(fit$ps),
               tolerance = 1e-12)
})

test_that("psave_match forwards arguments verbatim", {
  skip_if_not_installed("MatchIt")
  sim <- make_sim_data()
  fit <- psave(sim_formula(), data = sim, outcome = ~ y,
               ps.methods = "glm", prog.methods = "glm")
  m <- psave_match(fit, method = "nearest", caliper = 0.5)
  expect_s3_class(m, "matchit")
  expect_equal(unname(as.numeric(m$distance)), unname(fit$ps),
               tolerance = 1e-12)
  # a caliper was actually applied
  expect_false(is.null(m$caliper))
})

test_that("psave_weight returns a weightit object with ps = fit$ps", {
  skip_if_not_installed("WeightIt")
  sim <- make_sim_data()
  fit <- psave(sim_formula(), data = sim, outcome = ~ y,
               ps.methods = "glm", prog.methods = "glm")
  w <- psave_weight(fit)
  expect_s3_class(w, "weightit")
  expect_equal(unname(as.numeric(w$ps)), unname(fit$ps), tolerance = 1e-12)
  expect_identical(w$estimand, fit$estimand)
  # WeightIt's ATT weights from this ps must equal the stored weights
  expect_equal(unname(as.numeric(w$weights)), unname(fit$weights),
               tolerance = 1e-8)
})
