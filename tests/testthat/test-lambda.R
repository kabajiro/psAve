# D.14 (ii), (iv), (x): lambda grid search — hand-computed argmins for every
# criterion, first-row tie-breaking, vertex (average = FALSE) mode, path.
# All expected numbers derived in helper-fixtures.R comments.

test_that("ATT prog criterion selects the hand-computed argmin", {
  fx <- make_fixture()
  fit <- psave(A ~ x, data = fx$data, outcome = ~ y,
               ps.matrix = fx$E, prog.matrix = fx$G1,
               estimand = "ATT", criterion = "prog", step = 0.5)
  expect_s3_class(fit, "psave")
  expect_identical(fit$estimand, "ATT")
  expect_identical(fit$criterion, "prog")
  # grid values: (31/21, 1/3, 2/3) -> argmin at lambda = (.5, .5)
  expect_equal(fit$lambda, c(e1 = 0.5, e2 = 0.5))
  expect_equal(fit$criterion.value, 1 / 3, tolerance = 1e-10)
  # averaged PS = (e1 + e2)/2, named by rownames(data)
  expect_equal(unname(fit$ps), c(0.45, 0.50, 0.55, 0.50, 0.50, 0.50),
               tolerance = 1e-12)
  expect_identical(names(fit$ps), rownames(fx$data))
  # stored path: grid order, hand values, one selected row
  expect_true(is.data.frame(fit$path))
  expect_named(fit$path, c("e1", "e2", "value", "selected"))
  expect_identical(nrow(fit$path), 3L)
  expect_equal(fit$path$e1, c(1, 0.5, 0))
  expect_equal(fit$path$value, c(31 / 21, 1 / 3, 2 / 3), tolerance = 1e-10)
  expect_identical(fit$path$selected, c(FALSE, TRUE, FALSE))
})

test_that("ATE prog criterion selects the hand-computed argmin", {
  fx <- make_fixture()
  fit <- psave(A ~ x, data = fx$data, outcome = ~ y,
               ps.matrix = fx$E, prog.matrix = fx$G1,
               estimand = "ATE", criterion = "prog", step = 0.5)
  # grid values: (35/33, 359/897, 205/1221) -> argmin at lambda = (0, 1)
  expect_equal(fit$lambda, c(e1 = 0, e2 = 1))
  expect_equal(fit$criterion.value, 205 / 1221, tolerance = 1e-10)
  expect_equal(fit$path$value, c(35 / 33, 359 / 897, 205 / 1221),
               tolerance = 1e-10)
  expect_equal(unname(fit$ps), unname(fx$E[, "e2"]))
})

test_that("logloss criterion selects the hand-computed argmin, no outcome", {
  fx <- make_fixture()
  fit <- psave(A ~ x, data = fx$data,
               ps.matrix = fx$E, criterion = "logloss", step = 0.5)
  # grid values 0.76753, 0.69482, 0.77433 -> argmin at lambda = (.5, .5)
  expect_equal(fit$lambda, c(e1 = 0.5, e2 = 0.5))
  expect_equal(fit$criterion.value,
               -(log(0.45) + log(0.55) + 4 * log(0.5)) / 6,
               tolerance = 1e-10)
  # logloss without outcome: no prognostic side
  expect_null(fit$prog)
  expect_null(fit$gamma)
})

test_that("smd criterion selects the hand-computed argmin", {
  fx <- make_fixture()
  fit <- psave(A ~ x, data = fx$data, outcome = ~ y,
               ps.matrix = fx$E, prog.matrix = fx$G1,
               criterion = "smd", step = 0.5)
  # grid values: (12/7, 1, 2/7) -> argmin at lambda = (0, 1)
  expect_equal(fit$lambda, c(e1 = 0, e2 = 1))
  expect_equal(fit$criterion.value, 2 / 7, tolerance = 1e-10)
})

test_that("flat criterion ties select the first grid row", {
  # D.14 (iv): two identical candidate columns -> ebar identical at every
  # grid point -> criterion flat -> lambda = (1, 0) (first row wins)
  fx <- make_fixture()
  Etie <- cbind(p1 = fx$E[, "e1"], p2 = fx$E[, "e1"])
  fit <- psave(A ~ x, data = fx$data,
               ps.matrix = Etie, criterion = "logloss", step = 0.5)
  expect_equal(fit$lambda, c(p1 = 1, p2 = 0))
})

test_that("average = FALSE restricts selection to the simplex vertices", {
  # D.14 (x): vertex values (31/21, 2/3, 89/79) -> best single candidate e2
  fx <- make_fixture()
  fit <- psave(A ~ x, data = fx$data, outcome = ~ y,
               ps.matrix = fx$E3, prog.matrix = fx$G1,
               step = 0.5, average = FALSE)
  expect_equal(fit$lambda, c(e1 = 0, e2 = 1, e3 = 0))
  expect_true(all(fit$lambda %in% c(0, 1)))   # a unit vector
  expect_equal(sum(fit$lambda), 1)
  expect_equal(fit$criterion.value, 2 / 3, tolerance = 1e-10)
  expect_false(fit$average)
  expect_equal(unname(fit$ps), unname(fx$E3[, "e2"]))
})

test_that("keep.path = FALSE drops the path", {
  fx <- make_fixture()
  fit <- psave(A ~ x, data = fx$data, outcome = ~ y,
               ps.matrix = fx$E, prog.matrix = fx$G1,
               step = 0.5, keep.path = FALSE)
  expect_null(fit$path)
  # selection unaffected
  expect_equal(fit$lambda, c(e1 = 0.5, e2 = 0.5))
})

test_that("criterion.value agrees with psave_criteria at the selected ps", {
  fx <- make_fixture()
  fit <- psave(A ~ x, data = fx$data, outcome = ~ y,
               ps.matrix = fx$E, prog.matrix = fx$G1, step = 0.5)
  vals <- psave_criteria(ps = fit$ps, treat = fx$data$A,
                         covs = cbind(x = fx$data$x), prog = fit$prog,
                         estimand = fit$estimand)
  expect_equal(unname(vals[["prog"]]), fit$criterion.value, tolerance = 1e-10)
})
