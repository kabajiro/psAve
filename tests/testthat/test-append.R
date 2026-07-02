# ps.append / prog.append: extra user-supplied candidates appended AFTER the
# base candidate set (from *.methods or *.matrix).  All tests use the 6-row
# hand-computable fixture (see helper-fixtures.R) with glm-only / matrix-only
# paths, so no Suggests packages are needed and everything is deterministic.
#
# Extra hand derivations used below:
#   ec = (.50, .50, .50, 1/3, .50, .20)
#     ATT control weights e/(1-e) at u4..u6: (1/2, 1, 1/4), sum 7/4
#     weighted control mean of g = (1, 2, 4):
#       (1/2*1 + 1*2 + 1/4*4)/(7/4) = (7/2)/(7/4) = 2 = treated mean of g
#     => prog criterion (ATT) is EXACTLY 0 at full weight on ec.
#   Noise candidates b1 = (.5,.5,.5,.48,.5,.52), b2 = (.5,.5,.5,.52,.5,.48)
#     and their step-0.5 mixtures with ec all have criterion > 0.12, so the
#     argmin is the appended vertex (0, 0, 1).

test_that("vector ps.append becomes one extra candidate labeled \"append\"", {
  fx <- make_fixture()
  e3 <- unname(fx$E3[, "e3"])
  fit <- psave(A ~ x, data = fx$data, ps.matrix = fx$E, ps.append = e3,
               criterion = "logloss", step = 0.5)
  # base M = 2 plus one appended candidate
  expect_identical(names(fit$lambda), c("e1", "e2", "append"))
  expect_identical(colnames(fit$ps.candidates), c("e1", "e2", "append"))
  expect_equal(unname(fit$ps.candidates[, "append"]), e3, tolerance = 1e-15)
  # everything label-driven follows: diagnostics rows and the path columns
  expect_identical(rownames(fit$diagnostics), c("e1", "e2", "append", "average"))
  expect_true("append" %in% names(fit$path))
})

test_that("matrix ps.append appends several named candidates end-to-end", {
  fx <- make_fixture()
  App <- cbind(a1 = unname(fx$E3[, "e3"]),
               a2 = c(0.60, 0.50, 0.40, 0.30, 0.50, 0.70))
  fit <- psave(A ~ x, data = fx$data, ps.matrix = fx$E, ps.append = App,
               criterion = "logloss", step = 0.5)
  expect_identical(names(fit$lambda), c("e1", "e2", "a1", "a2"))
  expect_equal(unname(fit$ps.candidates[, c("a1", "a2")]), unname(App),
               tolerance = 1e-15)
  expect_equal(unname(fit$info$grid.size["lambda"]),
               nrow(simplex_grid(4, 0.5)))
  expect_true(all(fit$ps >= 0.01 & fit$ps <= 0.99))
  # an all-numeric data.frame is accepted and equivalent
  fit2 <- psave(A ~ x, data = fx$data, ps.matrix = fx$E,
                ps.append = as.data.frame(App),
                criterion = "logloss", step = 0.5)
  expect_equal(fit2$lambda, fit$lambda, tolerance = 1e-15)
})

test_that("ps.append is validated like ps.matrix", {
  fx <- make_fixture()
  ec <- c(0.50, 0.50, 0.50, 1 / 3, 0.50, 0.20)
  # unnamed matrix columns
  expect_error(
    psave(A ~ x, data = fx$data, ps.matrix = fx$E,
          ps.append = matrix(ec, ncol = 1L),
          criterion = "logloss", step = 0.5),
    "column names|colnames"
  )
  # name collision with a base candidate label -> suggest renaming
  expect_error(
    psave(A ~ x, data = fx$data, ps.matrix = fx$E,
          ps.append = cbind(e1 = ec),
          criterion = "logloss", step = 0.5),
    "rename"
  )
  # values must lie strictly in (0, 1)
  bad <- ec; bad[1L] <- 1.2
  expect_error(
    psave(A ~ x, data = fx$data, ps.matrix = fx$E, ps.append = bad,
          criterion = "logloss", step = 0.5),
    "strictly inside"
  )
  # no missing values
  bad2 <- ec; bad2[2L] <- NA
  expect_error(
    psave(A ~ x, data = fx$data, ps.matrix = fx$E, ps.append = bad2,
          criterion = "logloss", step = 0.5),
    "missing or non-finite"
  )
  # rownames, when present, must match rownames(data)
  Am <- cbind(a1 = ec)
  rownames(Am) <- paste0("v", 1:6)
  expect_error(
    psave(A ~ x, data = fx$data, ps.matrix = fx$E, ps.append = Am,
          criterion = "logloss", step = 0.5),
    "rownames"
  )
  # a vector must have one value per unit
  expect_error(
    psave(A ~ x, data = fx$data, ps.matrix = fx$E, ps.append = ec[1:5],
          criterion = "logloss", step = 0.5),
    "one value per unit"
  )
})

test_that("appended candidates are clipped before averaging", {
  fx <- make_fixture()
  v <- c(0.001, 0.5, 0.999, 0.2, 0.5, 0.8)
  fit <- psave(A ~ x, data = fx$data, ps.matrix = fx$E, ps.append = v,
               criterion = "logloss", step = 0.5)
  expect_equal(unname(fit$ps.candidates[, "append"]),
               c(0.01, 0.5, 0.99, 0.2, 0.5, 0.8), tolerance = 1e-15)
  expect_true(all(fit$ps >= 0.01 & fit$ps <= 0.99))
  # the averaged ps is exactly the convex combination of clipped candidates
  expect_equal(unname(fit$ps),
               unname(drop(fit$ps.candidates %*% fit$lambda)),
               tolerance = 1e-15)
})

test_that("an appended candidate minimizing the criterion receives full weight", {
  fx <- make_fixture()
  Enoise <- cbind(b1 = c(0.50, 0.50, 0.50, 0.48, 0.50, 0.52),
                  b2 = c(0.50, 0.50, 0.50, 0.52, 0.50, 0.48))
  ec <- c(0.50, 0.50, 0.50, 1 / 3, 0.50, 0.20)   # prog criterion exactly 0
  fit <- psave(A ~ x, data = fx$data, outcome = ~ y, ps.matrix = Enoise,
               ps.append = ec, prog.matrix = fx$G1, step = 0.5)
  expect_equal(unname(fit$lambda), c(0, 0, 1))
  expect_identical(names(fit$lambda), c("b1", "b2", "append"))
  expect_equal(fit$criterion.value, 0, tolerance = 1e-12)
})

test_that("ties between a base candidate and its appended duplicate go to the base", {
  fx <- make_fixture()
  E1 <- fx$E[, "e1", drop = FALSE]
  dup <- cbind(e1copy = unname(fx$E[, "e1"]))
  fit <- psave(A ~ x, data = fx$data, ps.matrix = E1, ps.append = dup,
               criterion = "logloss", step = 0.5)
  # ebar is identical at every grid row -> flat criterion -> first row wins,
  # i.e. all weight on the BASE candidate (appended columns come last)
  expect_equal(fit$lambda, c(e1 = 1, e1copy = 0))
})

test_that("average = FALSE includes appended candidates among the vertices", {
  fx <- make_fixture()
  ec <- c(0.50, 0.50, 0.50, 1 / 3, 0.50, 0.20)   # best vertex by construction
  fit <- psave(A ~ x, data = fx$data, outcome = ~ y, ps.matrix = fx$E,
               ps.append = ec, prog.matrix = fx$G1, step = 0.5,
               average = FALSE)
  # vertex values: e1 = 31/21, e2 = 2/3 (helper header), append = 0
  expect_equal(unname(fit$info$grid.size["lambda"]), 3)
  expect_identical(names(fit$lambda), c("e1", "e2", "append"))
  expect_equal(unname(fit$lambda), c(0, 0, 1))
  expect_equal(fit$criterion.value, 0, tolerance = 1e-12)
})

test_that("prog.append appends prognostic candidates and gamma follows", {
  fx <- make_fixture()
  G1only <- fx$G2[, "g1", drop = FALSE]
  gv <- c(2, 3, 4, 1, 2, 4)   # equals y on the untreated units -> MSE 0
  fit <- psave(A ~ x, data = fx$data, outcome = ~ y, ps.matrix = fx$E,
               prog.matrix = G1only, prog.append = gv, step = 0.5)
  expect_identical(names(fit$gamma), c("g1", "append"))
  expect_equal(fit$gamma, c(g1 = 0, append = 1))
  expect_identical(names(fit$gamma.mse), c("g1", "append", "average"))
  # g1 untreated MSE = 5/3 (helper header); appended column is exact
  expect_equal(unname(fit$gamma.mse), c(5 / 3, 0, 0), tolerance = 1e-12)
  expect_identical(colnames(fit$prog.candidates), c("g1", "append"))
  expect_equal(unname(fit$prog), gv, tolerance = 1e-15)
})

test_that("prognostic candidates without `outcome` are an error, not ignored", {
  fx <- make_fixture()
  expect_error(
    psave(A ~ x, data = fx$data, ps.matrix = fx$E, prog.matrix = fx$G1,
          criterion = "logloss", step = 0.5),
    "outcome"
  )
  expect_error(
    psave(A ~ x, data = fx$data, ps.matrix = fx$E,
          prog.append = c(2, 3, 4, 1, 2, 4),
          criterion = "logloss", step = 0.5),
    "outcome"
  )
})

test_that("ps.append composes with fitted ps.methods candidates", {
  fx <- make_fixture()
  fit <- psave(A ~ x, data = fx$data, ps.methods = "glm",
               ps.append = unname(fx$E[, "e2"]),
               criterion = "logloss", step = 0.5)
  expect_identical(names(fit$lambda), c("glm", "append"))
  expect_identical(colnames(fit$ps.candidates), c("glm", "append"))
  # provenance: the fitted learner keeps its engine; the appended candidate
  # is recorded as a user-supplied "append" entry
  expect_identical(fit$info$learners$ps$glm$engine, "glm")
  expect_identical(fit$info$learners$ps$append$type, "append")
})

test_that("`control` entries naming append labels still warn as unused", {
  fx <- make_fixture()
  expect_warning(
    psave(A ~ x, data = fx$data, ps.methods = "glm",
          ps.append = cbind(extra = unname(fx$E[, "e2"])),
          control = list(extra = list(maxit = 5)),
          criterion = "logloss", step = 0.5),
    "do not match any fitted learner"
  )
})
