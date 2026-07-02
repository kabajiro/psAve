# D.14 (vi), (viii), (ix): ps.matrix / prog.matrix escape hatches,
# clipping semantics, degenerate M = 1 path.

test_that("ps.matrix bypasses fitting and is stored (already in range)", {
  fx <- make_fixture()
  fit <- psave(A ~ x, data = fx$data, outcome = ~ y,
               ps.matrix = fx$E, prog.matrix = fx$G1, step = 0.5)
  # values already inside [0.01, 0.99]: stored candidates identical to input
  expect_equal(unname(fit$ps.candidates), unname(fx$E), tolerance = 1e-15)
  expect_identical(colnames(fit$ps.candidates), c("e1", "e2"))
  expect_identical(names(fit$lambda), c("e1", "e2"))
})

test_that("candidate ps values are clipped to [0.01, 0.99] before averaging", {
  # D.14 (vi): out-of-clip values 0.001 and 0.999 are pulled to the bounds
  fx <- make_fixture()
  Eclip <- cbind(p1 = c(0.001, 0.5, 0.999, 0.2, 0.5, 0.8), p2 = fx$E[, "e2"])
  fit <- psave(A ~ x, data = fx$data,
               ps.matrix = Eclip, criterion = "logloss", step = 0.5)
  expect_equal(unname(fit$ps.candidates[, "p1"]),
               c(0.01, 0.5, 0.99, 0.2, 0.5, 0.8), tolerance = 1e-15)
  # averaged ps stays inside the clip bounds by convexity...
  expect_true(all(fit$ps >= 0.01 & fit$ps <= 0.99))
  # ...and is EXACTLY the convex combination of the clipped candidates
  # (no distorting re-clip after averaging)
  expect_equal(unname(fit$ps),
               unname(drop(fit$ps.candidates %*% fit$lambda)),
               tolerance = 1e-15)
})

test_that("a custom clip argument is respected", {
  fx <- make_fixture()
  Eclip <- cbind(p1 = c(0.05, 0.5, 0.95, 0.2, 0.5, 0.8), p2 = fx$E[, "e2"])
  fit <- psave(A ~ x, data = fx$data, ps.matrix = Eclip,
               criterion = "logloss", step = 0.5, clip = c(0.1, 0.9))
  expect_equal(unname(fit$ps.candidates[, "p1"]),
               c(0.1, 0.5, 0.9, 0.2, 0.5, 0.8), tolerance = 1e-15)
  expect_true(all(fit$ps >= 0.1 & fit$ps <= 0.9))
})

test_that("ps.matrix is validated", {
  fx <- make_fixture()
  # colnames required
  expect_error(
    psave(A ~ x, data = fx$data, ps.matrix = unname(fx$E),
          criterion = "logloss", step = 0.5),
    "column names|colnames"
  )
  # values must lie strictly in (0, 1)
  Ebad <- fx$E; Ebad[1, 1] <- 1
  expect_error(
    psave(A ~ x, data = fx$data, ps.matrix = Ebad,
          criterion = "logloss", step = 0.5),
    "ps.matrix|0"
  )
  Ebad2 <- fx$E; Ebad2[2, 2] <- 0
  expect_error(
    psave(A ~ x, data = fx$data, ps.matrix = Ebad2,
          criterion = "logloss", step = 0.5),
    "ps.matrix|0"
  )
  # row count must match the data
  expect_error(
    psave(A ~ x, data = fx$data, ps.matrix = fx$E[1:5, , drop = FALSE],
          criterion = "logloss", step = 0.5),
    "ps.matrix|row"
  )
})

test_that("prog.matrix bypasses prognostic fitting", {
  fx <- make_fixture()
  fit <- psave(A ~ x, data = fx$data, outcome = ~ y,
               ps.matrix = fx$E, prog.matrix = fx$G2, step = 0.5)
  expect_equal(unname(fit$prog.candidates), unname(fx$G2), tolerance = 1e-15)
  expect_identical(colnames(fit$prog.candidates), c("g1", "g2"))
  expect_identical(names(fit$gamma), c("g1", "g2"))
})

test_that("M = 1 ps.matrix is a valid degenerate path", {
  fx <- make_fixture()
  E1 <- fx$E[, "e1", drop = FALSE]
  fit <- psave(A ~ x, data = fx$data, outcome = ~ y,
               ps.matrix = E1, prog.matrix = fx$G1, step = 0.5)
  expect_equal(fit$lambda, c(e1 = 1))
  expect_equal(unname(fit$ps), unname(fx$E[, "e1"]))
  # criterion value at the single candidate = hand value 31/21 (ATT, prog)
  expect_equal(fit$criterion.value, 31 / 21, tolerance = 1e-10)
})
