# D.14 (iii), (viii): gamma selection — hand-checked untreated MSE,
# exact-candidate selection, K = 1 degenerate path, prog.target handling.
# All expected numbers derived in helper-fixtures.R comments.

test_that("gamma selects the candidate that reproduces Y on untreated units", {
  fx <- make_fixture()
  fit <- psave(A ~ x, data = fx$data, outcome = ~ y,
               ps.matrix = fx$E, prog.matrix = fx$G2, step = 0.5)
  # g2 equals y exactly on the untreated arm -> MSE(0,1) = 0, unique minimum
  expect_equal(fit$gamma, c(g1 = 0, g2 = 1))
  # gbar = g2, named by rownames(data)
  expect_equal(unname(fit$prog), unname(fx$G2[, "g2"]))
  expect_identical(names(fit$prog), rownames(fx$data))
})

test_that("gamma.mse reports hand-computed untreated MSEs", {
  fx <- make_fixture()
  fit <- psave(A ~ x, data = fx$data, outcome = ~ y,
               ps.matrix = fx$E, prog.matrix = fx$G2, step = 0.5)
  # y0 = (1,2,4); g1_0 = (2,2,2) -> MSE = (1+0+4)/3 = 5/3; g2_0 = y0 -> 0;
  # selected average (gamma = (0,1)) -> 0
  expect_equal(fit$gamma.mse, c(g1 = 5 / 3, g2 = 0, average = 0),
               tolerance = 1e-12)
})

test_that("K = 1 skips the gamma grid and returns gamma = 1", {
  fx <- make_fixture()
  fit <- psave(A ~ x, data = fx$data, outcome = ~ y,
               ps.matrix = fx$E, prog.matrix = fx$G1, step = 0.5)
  expect_equal(fit$gamma, c(g = 1))
  # sole candidate g equals y on untreated -> MSE = 0 for candidate and average
  expect_equal(fit$gamma.mse, c(g = 0, average = 0), tolerance = 1e-12)
  expect_equal(unname(fit$prog), unname(fx$G1[, "g"]))
})

test_that("prog.target = 'average' and a named candidate select differently", {
  fx <- make_fixture()
  # 'average': gbar = g2 -> criterion path (31/21, 1/3, 2/3) -> lambda (.5,.5)
  fit.a <- psave(A ~ x, data = fx$data, outcome = ~ y,
                 ps.matrix = fx$E, prog.matrix = fx$G2,
                 prog.target = "average", step = 0.5)
  expect_equal(fit.a$lambda, c(e1 = 0.5, e2 = 0.5))
  expect_equal(fit.a$criterion.value, 1 / 3, tolerance = 1e-10)
  # 'g1': control-constant candidate -> criterion identically 0 -> first row
  fit.g <- psave(A ~ x, data = fx$data, outcome = ~ y,
                 ps.matrix = fx$E, prog.matrix = fx$G2,
                 prog.target = "g1", step = 0.5)
  expect_equal(fit.g$lambda, c(e1 = 1, e2 = 0))
  expect_equal(fit.g$criterion.value, 0, tolerance = 1e-12)
  expect_identical(fit.g$prog.target, "g1")
  # gamma selection is independent of prog.target
  expect_equal(fit.g$gamma, c(g1 = 0, g2 = 1))
})

test_that("prog.target must name a supplied prognostic candidate", {
  fx <- make_fixture()
  expect_error(
    psave(A ~ x, data = fx$data, outcome = ~ y,
          ps.matrix = fx$E, prog.matrix = fx$G2,
          prog.target = "nope", step = 0.5),
    "prog.target"
  )
})
