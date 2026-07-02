# D.14 (xi) plus extractors, learner paths, and the glm-only core.
# print/summary are checked with expect_output fragments, NOT snapshots.

test_that("glm-only end-to-end run works (M = 1, K = 1 degenerate path)", {
  sim <- make_sim_data()
  fit <- psave(sim_formula(), data = sim, outcome = ~ y,
               ps.methods = "glm", prog.methods = "glm")
  expect_s3_class(fit, "psave")
  expect_equal(fit$lambda, c(glm = 1))
  expect_equal(fit$gamma,  c(glm = 1))
  expect_true(all(fit$ps >= 0.01 & fit$ps <= 0.99))
  expect_identical(length(fit$ps), nrow(sim))
  expect_identical(length(fit$weights), nrow(sim))
  expect_identical(fit$estimand, "ATT")
  expect_identical(fit$criterion, "prog")
  expect_identical(fit$outcome.name, "y")
})

test_that("print shows estimand, criterion, lambda and the next-step call", {
  fx <- make_fixture()
  fit <- psave(A ~ x, data = fx$data, outcome = ~ y,
               ps.matrix = fx$E, prog.matrix = fx$G1, step = 0.5)
  expect_output(print(fit), "ATT")
  expect_output(print(fit), "prog")
  expect_output(print(fit), "lambda", ignore.case = TRUE)
  expect_output(print(fit), "psave_match|matchit", ignore.case = TRUE)
  # print returns its argument invisibly (ecosystem convention)
  expect_invisible(print(fit))
})

test_that("summary returns a summary.psave with the expected content", {
  fx <- make_fixture()
  fit <- psave(A ~ x, data = fx$data, outcome = ~ y,
               ps.matrix = fx$E, prog.matrix = fx$G1, step = 0.5)
  s <- summary(fit)
  expect_s3_class(s, "summary.psave")
  expect_output(print(s), "lambda", ignore.case = TRUE)
  expect_output(print(s), "smd", ignore.case = TRUE)
  expect_output(print(s), "ATT")
})

test_that("fitted and weights extractors return the stored components", {
  fx <- make_fixture()
  fit <- psave(A ~ x, data = fx$data, outcome = ~ y,
               ps.matrix = fx$E, prog.matrix = fx$G1, step = 0.5)
  expect_equal(fitted(fit), fit$ps)
  expect_equal(fitted(fit, type = "ps"), fit$ps)
  expect_equal(fitted(fit, type = "prog"), fit$prog)
  expect_equal(weights(fit), fit$weights)
})

test_that("diagnostics table covers all candidates plus the average", {
  fx <- make_fixture()
  fit <- psave(A ~ x, data = fx$data, outcome = ~ y,
               ps.matrix = fx$E, prog.matrix = fx$G1, step = 0.5)
  expect_true(is.data.frame(fit$diagnostics))
  expect_identical(dim(fit$diagnostics), c(3L, 4L))   # M + 1 = 3 rows
})

test_that("balance table has the spec columns and a prog row", {
  fx <- make_fixture()
  fit <- psave(A ~ x, data = fx$data, outcome = ~ y,
               ps.matrix = fx$E, prog.matrix = fx$G1, step = 0.5)
  expect_true(is.data.frame(fit$balance))
  expect_named(fit$balance, c("smd.un", "smd.wt", "ks.un", "ks.wt"))
  expect_identical(nrow(fit$balance), ncol(fit$covs) + 1L)  # covariates + prog
})

test_that("info records the fitting configuration", {
  fx <- make_fixture()
  fit <- psave(A ~ x, data = fx$data, outcome = ~ y,
               ps.matrix = fx$E, prog.matrix = fx$G1, step = 0.5)
  expect_true(is.list(fit$info))
  expect_true(fit$info$n == 6)
  expect_true(fit$info$step == 0.5)
  # stored data powers the wrappers verbatim
  expect_identical(fit$data, fx$data)
  expect_true(inherits(fit$formula, "formula"))
})

test_that("predict requires keep.fits = TRUE", {
  sim <- make_sim_data()
  fit <- psave(sim_formula(), data = sim, outcome = ~ y,
               ps.methods = "glm", prog.methods = "glm")
  expect_error(predict(fit, newdata = sim), "keep.fits")
})

test_that("predict reproduces in-sample averaged scores with glm learners", {
  sim <- make_sim_data()
  fit <- psave(sim_formula(), data = sim, outcome = ~ y,
               ps.methods = "glm", prog.methods = "glm", keep.fits = TRUE)
  # glm is deterministic and the in-sample glm PS lie strictly inside the
  # clip bounds for this fixture, so predicting the training data must
  # reproduce the stored averaged scores
  expect_equal(unname(predict(fit, newdata = sim, type = "ps")),
               unname(fit$ps), tolerance = 1e-8)
  expect_equal(unname(predict(fit, newdata = sim, type = "prog")),
               unname(fit$prog), tolerance = 1e-8)
})

test_that("bal.tab works on a psave object", {
  fx <- make_fixture()
  fit <- psave(A ~ x, data = fx$data, outcome = ~ y,
               ps.matrix = fx$E, prog.matrix = fx$G1, step = 0.5)
  bt <- cobalt::bal.tab(fit)
  expect_s3_class(bt, "bal.tab")
})

test_that("verbose = TRUE emits progress messages", {
  fx <- make_fixture()
  expect_message(
    psave(A ~ x, data = fx$data, outcome = ~ y,
          ps.matrix = fx$E, prog.matrix = fx$G1, step = 0.5, verbose = TRUE)
  )
})

test_that("multi-learner fit works with rpart", {
  skip_if_not_installed("rpart")
  sim <- make_sim_data()
  set.seed(101)
  fit <- psave(sim_formula(), data = sim, outcome = ~ y,
               ps.methods = c("glm", "rpart"),
               prog.methods = c("glm", "rpart"), step = 0.25)
  expect_identical(names(fit$lambda), c("glm", "rpart"))
  expect_identical(names(fit$gamma),  c("glm", "rpart"))
  expect_equal(sum(fit$lambda), 1, tolerance = 1e-12)
  expect_equal(sum(fit$gamma), 1, tolerance = 1e-12)
  expect_true(all(fit$lambda >= 0) && all(fit$gamma >= 0))
  expect_true(all(fit$ps >= 0.01 & fit$ps <= 0.99))
  expect_identical(dim(fit$diagnostics), c(3L, 4L))
})

test_that("full default learner set runs", {
  skip_if_not_installed("rpart")
  skip_if_not_installed("ranger")
  skip_if_not_installed("xgboost")
  sim <- make_sim_data()
  set.seed(202)
  fit <- psave(sim_formula(), data = sim, outcome = ~ y, step = 0.25)
  expect_identical(names(fit$lambda),
                   c("glm", "rpart", "ranger", "xgboost"))
  expect_equal(sum(fit$lambda), 1, tolerance = 1e-12)
  expect_true(all(fit$ps >= 0.01 & fit$ps <= 0.99))
})

test_that("SL.* labels are passed through to SuperLearner", {
  skip_if_not_installed("SuperLearner")
  sim <- make_sim_data()
  set.seed(303)
  fit <- psave(sim_formula(), data = sim,
               ps.methods = c("glm", "SL.glm"),
               criterion = "logloss", step = 0.5)
  expect_identical(names(fit$lambda), c("glm", "SL.glm"))
  expect_equal(sum(fit$lambda), 1, tolerance = 1e-12)
})
