# D.14 (i): simplex_grid — row counts, integer arithmetic, documented order

test_that("simplex_grid has the exact composition row counts", {
  # G = choose(n_steps + M - 1, M - 1)
  expect_identical(nrow(simplex_grid(4, 0.05)), as.integer(choose(23, 3)))  # 1771
  expect_identical(nrow(simplex_grid(4, 0.05)), 1771L)
  expect_identical(nrow(simplex_grid(3, 0.1)),  as.integer(choose(12, 2)))  # 66
  expect_identical(nrow(simplex_grid(2, 0.05)), as.integer(choose(21, 1)))  # 21
  expect_identical(nrow(simplex_grid(5, 0.25)), as.integer(choose(8, 4)))   # 70
  expect_identical(ncol(simplex_grid(4, 0.05)), 4L)
})

test_that("simplex_grid rows sum to 1 exactly in integer arithmetic", {
  g <- simplex_grid(4, 0.05)
  expect_true(is.matrix(g))
  cnt <- round(g * 20)                       # back to integer compositions
  expect_lt(max(abs(g * 20 - cnt)), 1e-9)    # entries are k/20 exactly
  expect_true(all(rowSums(cnt) == 20))       # every composition sums to 20
  expect_true(all(g >= 0))
  # no float-filter bug: no duplicated or dropped compositions
  expect_identical(nrow(unique(cnt)), nrow(g))
})

test_that("simplex_grid(2, 0.5) is exactly the documented 3-row grid", {
  g <- simplex_grid(2, 0.5)
  expect_equal(unname(g), rbind(c(1, 0), c(0.5, 0.5), c(0, 1)),
               ignore_attr = TRUE)
})

test_that("simplex_grid order is descending-lexicographic (tie-break order)", {
  # independent recursive reference enumeration: c1 from n_steps down to 0,
  # then c2 descending on the remainder, etc.
  ref_comps <- function(M, n) {
    if (M == 1L) return(matrix(n, 1L, 1L))
    do.call(rbind, lapply(seq.int(n, 0L),
                          function(k) cbind(k, ref_comps(M - 1L, n - k))))
  }
  for (spec in list(c(3, 4), c(4, 20), c(2, 2))) {
    M <- spec[1]; n <- spec[2]
    expected <- ref_comps(M, n) / n
    dimnames(expected) <- NULL
    got <- simplex_grid(M, 1 / n)
    dimnames(got) <- NULL
    expect_equal(got, expected)
  }
  g <- simplex_grid(4, 0.05)
  expect_equal(unname(g[1, ]),       c(1, 0, 0, 0))  # first row = first vertex
  expect_equal(unname(g[nrow(g), ]), c(0, 0, 0, 1))  # last row = last vertex
})

test_that("simplex_grid(1, step) is the single row 1", {
  g <- simplex_grid(1, 0.05)
  expect_identical(nrow(g), 1L)
  expect_equal(unname(g[1, ]), 1)
})

test_that("simplex_grid errors when 1/step is not an integer", {
  expect_error(simplex_grid(3, 0.3), "step")
  expect_error(simplex_grid(2, 0.07), "step")
})
