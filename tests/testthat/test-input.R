# D.14 (vii) and D.1-D.2: input validation, treatment coercion, covariate
# matrix construction.  All tests here run with ps.matrix/prog.matrix or
# glm only (no Suggests needed).

test_that("NA in the treatment errors, never drops", {
  fx <- make_fixture()
  d <- fx$data
  d$A[1] <- NA
  expect_error(
    psave(A ~ x, data = d, outcome = ~ y,
          ps.matrix = fx$E, prog.matrix = fx$G1, step = 0.5),
    "[Mm]issing|NA"
  )
})

test_that("NA in a formula covariate errors, never drops", {
  fx <- make_fixture()
  d <- fx$data
  d$x[2] <- NA
  expect_error(
    psave(A ~ x, data = d, outcome = ~ y,
          ps.matrix = fx$E, prog.matrix = fx$G1, step = 0.5),
    "[Mm]issing|NA"
  )
})

test_that("NA in the outcome errors when the outcome is used", {
  fx <- make_fixture()
  d <- fx$data
  d$y[4] <- NA
  expect_error(
    psave(A ~ x, data = d, outcome = ~ y,
          ps.matrix = fx$E, prog.matrix = fx$G1, step = 0.5),
    "[Mm]issing|NA"
  )
})

test_that("NA in an unused column is harmless", {
  # criterion = "logloss" without outcome: y is never touched
  fx <- make_fixture()
  d <- fx$data
  d$y[4] <- NA
  fit <- psave(A ~ x, data = d, ps.matrix = fx$E,
               criterion = "logloss", step = 0.5)
  expect_s3_class(fit, "psave")
})

test_that("fewer than 2 units per arm errors", {
  fx <- make_fixture()
  d <- fx$data
  d$A <- c(1, 1, 1, 1, 1, 0)   # one control unit only
  expect_error(
    psave(A ~ x, data = d, ps.matrix = fx$E,
          criterion = "logloss", step = 0.5),
    "2"
  )
})

test_that("missing outcome with criterion = 'prog' gives the teaching error", {
  fx <- make_fixture()
  expect_error(
    psave(A ~ x, data = fx$data, ps.matrix = fx$E, prog.matrix = fx$G1,
          criterion = "prog", step = 0.5),
    "requires the outcome"
  )
})

test_that("treatment coercion follows MatchIt conventions", {
  fx <- make_fixture()
  target <- c(1L, 1L, 1L, 0L, 0L, 0L)
  # logical
  d <- fx$data; d$A <- c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE)
  fit <- psave(A ~ x, data = d, ps.matrix = fx$E,
               criterion = "logloss", step = 0.5)
  expect_equal(as.integer(fit$treat), target)
  # two-level factor: SECOND level (factor order) is treated
  d <- fx$data; d$A <- factor(c("t", "t", "t", "c", "c", "c"),
                              levels = c("c", "t"))
  fit <- psave(A ~ x, data = d, ps.matrix = fx$E,
               criterion = "logloss", step = 0.5)
  expect_equal(as.integer(fit$treat), target)
  # two-level character (alphabetical factor order: "b" treated)
  d <- fx$data; d$A <- c("b", "b", "b", "a", "a", "a")
  fit <- psave(A ~ x, data = d, ps.matrix = fx$E,
               criterion = "logloss", step = 0.5)
  expect_equal(as.integer(fit$treat), target)
})

test_that("a treatment with more than two levels errors", {
  fx <- make_fixture()
  d <- fx$data
  d$A <- factor(c("a", "a", "b", "b", "c", "c"))
  expect_error(
    psave(A ~ x, data = d, ps.matrix = fx$E,
          criterion = "logloss", step = 0.5),
    "treat"
  )
})

test_that("estimand and criterion values are validated", {
  fx <- make_fixture()
  expect_error(
    psave(A ~ x, data = fx$data, outcome = ~ y, ps.matrix = fx$E,
          prog.matrix = fx$G1, estimand = "ATO", step = 0.5),
    "ATT|estimand"
  )
  expect_error(
    psave(A ~ x, data = fx$data, outcome = ~ y, ps.matrix = fx$E,
          prog.matrix = fx$G1, criterion = "balance", step = 0.5),
    "criterion|should be one"
  )
})

test_that("family must be gaussian or binomial", {
  fx <- make_fixture()
  expect_error(
    psave(A ~ x, data = fx$data, outcome = ~ y, ps.matrix = fx$E,
          prog.matrix = fx$G1, family = poisson(), step = 0.5),
    "gaussian|binomial"
  )
})

test_that("an invalid step errors", {
  fx <- make_fixture()
  expect_error(
    psave(A ~ x, data = fx$data, outcome = ~ y, ps.matrix = fx$E,
          prog.matrix = fx$G1, step = 0.3),
    "step"
  )
})

test_that("covs uses full dummy expansion with a bin.vars attribute", {
  fx <- make_fixture()
  d <- fx$data
  d$f <- factor(c("a", "b", "c", "a", "b", "c"))
  fit <- psave(A ~ x + f, data = d, outcome = ~ y,
               ps.matrix = fx$E, prog.matrix = fx$G1, step = 0.5)
  # no intercept + EVERY factor level gets a column: x + 3 dummies = 4
  expect_true(is.matrix(fit$covs))
  expect_identical(ncol(fit$covs), 4L)
  bv <- attr(fit$covs, "bin.vars")
  expect_true(is.logical(bv))
  expect_identical(length(bv), 4L)
  # the three dummies are binary; x (4 distinct values) is not
  expect_identical(sum(bv), 3L)
})

test_that("output vectors are named by rownames(data)", {
  fx <- make_fixture()
  fit <- psave(A ~ x, data = fx$data, outcome = ~ y,
               ps.matrix = fx$E, prog.matrix = fx$G1, step = 0.5)
  expect_identical(names(fit$ps), rownames(fx$data))
  expect_identical(names(fit$prog), rownames(fx$data))
  expect_identical(length(fit$weights), nrow(fx$data))
})

test_that("two-sided outcome formula is accepted", {
  # glm-only path: y ~ z uses a distinct prognostic specification
  sim <- make_sim_data()
  fit <- psave(sim_formula(), data = sim, outcome = y ~ X1 + X2,
               ps.methods = "glm", prog.methods = "glm")
  expect_s3_class(fit, "psave")
  expect_identical(fit$outcome.name, "y")
})

test_that("outcome is optional for the outcome-free criteria (smd, ks, logloss)", {
  fx <- make_fixture()
  for (cr in c("smd", "ks", "logloss")) {
    fit <- psave(A ~ x, data = fx$data, criterion = cr,
                 ps.matrix = fx$E, step = 0.5)
    expect_s3_class(fit, "psave")
    expect_null(fit$prog)
    expect_null(fit$gamma)
    expect_true(is.na(fit$diagnostics[["prog"]][1L]))
  }
})

test_that("criterion = 'prog' still requires the outcome", {
  fx <- make_fixture()
  expect_error(
    psave(A ~ x, data = fx$data, criterion = "prog",
          ps.matrix = fx$E, step = 0.5),
    "requires the outcome"
  )
})
