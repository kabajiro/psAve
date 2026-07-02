# plot.R -- plot.psave: "balance" (cobalt::love.plot), "distribution" (base
# graphics), "criterion" (base graphics over the lambda path).

#' Plot a psave object
#'
#' Three diagnostic displays for a [psave()] fit:
#' \describe{
#'   \item{`"balance"`}{a Love plot of covariate and prognostic-score balance
#'     before/after weighting, via [cobalt::love.plot()] (dispatched through
#'     [bal.tab.psave()]; \pkg{cobalt} is an Import, so always available).}
#'   \item{`"distribution"`}{the distribution of the propensity scores by
#'     treatment group: grey curves for the candidate models, a colored curve
#'     for the selected average (base graphics). Inspect this for extreme
#'     candidate scores.}
#'   \item{`"criterion"`}{the selection criterion over the \eqn{\lambda} grid:
#'     exact for up to three candidates (a curve for M = 2, a colored grid map
#'     for M = 3); for M > 3, one profile per candidate (the minimum criterion
#'     value attainable at each value of \eqn{\lambda_m}). Requires the stored
#'     `path` (`keep.path = TRUE` and a grid of at most 100,000 rows).}
#' }
#'
#' @param x A `psave` object.
#' @param type One of `"balance"` (default), `"distribution"`, `"criterion"`.
#' @param ... For `"balance"`, further arguments to [cobalt::love.plot()]
#'   (e.g., `thresholds = 0.1`); otherwise further graphical parameters passed
#'   to the base plotting calls.
#'
#' @return For `"balance"`, the `ggplot` object from [cobalt::love.plot()]
#'   (invisibly, after printing); otherwise `x`, invisibly.
#' @seealso [psave()], [bal.tab.psave()], [cobalt::love.plot()]
#' @export
plot.psave <- function(x, type = c("balance", "distribution", "criterion"), ...) {
  type <- match.arg(type)
  switch(type,
         balance = .plot_balance(x, ...),
         distribution = .plot_distribution(x, ...),
         criterion = .plot_criterion(x, ...))
}

.plot_balance <- function(x, ...) {
  p <- cobalt::love.plot(x, ...)
  print(p)
  invisible(p)
}

.plot_distribution <- function(x, ...) {
  E <- x$ps.candidates
  labs <- colnames(E)
  op <- graphics::par(mfrow = c(1L, 2L))
  on.exit(graphics::par(op), add = TRUE)
  for (side in c("Treated (A = 1)", "Control (A = 0)")) {
    rows <- if (startsWith(side, "Treated")) x$treat == 1L else x$treat == 0L
    dens <- lapply(seq_len(ncol(E)),
                   function(m) stats::density(E[rows, m], from = 0, to = 1))
    dave <- stats::density(x$ps[rows], from = 0, to = 1)
    ymax <- max(vapply(c(dens, list(dave)), function(d) max(d$y), numeric(1L)))
    graphics::plot(NA, xlim = c(0, 1), ylim = c(0, ymax),
                   xlab = "Propensity score", ylab = "Density",
                   main = side, ...)
    for (d in dens) graphics::lines(d, col = "grey70")
    graphics::lines(dave, col = "firebrick", lwd = 2)
  }
  graphics::legend("topright",
                   legend = c(sprintf("candidates (%s)", paste(labs, collapse = ", ")),
                              "average"),
                   col = c("grey70", "firebrick"), lwd = c(1, 2), bty = "n",
                   cex = 0.8)
  invisible(x)
}

# Pure-base grey ramp (avoids a grDevices dependency): v in [0, 1], 0 = dark.
.grey_hex <- function(v) {
  g <- as.integer(round(38 + v * (225 - 38)))
  sprintf("#%02X%02X%02X", g, g, g)
}

.plot_criterion <- function(x, ...) {
  if (is.null(x$path)) {
    stop(paste0("No criterion path is stored (keep.path = FALSE, or the lambda grid ",
                "exceeded 100,000 rows). Re-run psave() with keep.path = TRUE."),
         call. = FALSE)
  }
  path <- x$path
  M <- length(x$lambda)
  labs <- names(x$lambda)
  val <- path$value
  ylab <- sprintf("criterion (%s)", x$criterion)
  sel <- which(path$selected)

  if (M == 1L) {
    stop("Only one candidate propensity score model: there is no lambda path to plot.",
         call. = FALSE)
  }
  if (M == 2L) {
    graphics::plot(path[[1L]], val, type = "l",
                   xlab = sprintf("lambda[%s]  (lambda[%s] = 1 - lambda[%s])",
                                  labs[1L], labs[2L], labs[1L]),
                   ylab = ylab, ...)
    graphics::points(path[[1L]][sel], val[sel], pch = 19, col = "firebrick")
  } else if (M == 3L) {
    ## exact: every grid point in the (lambda_1, lambda_2) plane, criterion as
    ## grey level (dark = small); the third coordinate is implied.
    v <- (val - min(val)) / max(max(val) - min(val), .Machine$double.eps)
    graphics::plot(path[[1L]], path[[2L]], pch = 15, col = .grey_hex(v),
                   xlab = sprintf("lambda[%s]", labs[1L]),
                   ylab = sprintf("lambda[%s]", labs[2L]),
                   main = sprintf("%s (dark = small; lambda[%s] implied)",
                                  ylab, labs[3L]), ...)
    graphics::points(path[[1L]][sel], path[[2L]][sel], pch = 0, cex = 1.6,
                     col = "firebrick", lwd = 2)
  } else {
    ## M > 3: one profile per candidate -- at each value of lambda_m, the
    ## minimum criterion over the remaining coordinates.
    prof <- lapply(seq_len(M), function(m) {
      xs <- sort(unique(path[[m]]))
      ys <- vapply(xs, function(v) min(val[path[[m]] == v]), numeric(1L))
      list(x = xs, y = ys)
    })
    ylim <- range(unlist(lapply(prof, `[[`, "y")))
    graphics::plot(NA, xlim = c(0, 1), ylim = ylim,
                   xlab = "lambda[m]", ylab = ylab,
                   main = "Profile criterion per candidate", ...)
    for (m in seq_len(M)) {
      graphics::lines(prof[[m]]$x, prof[[m]]$y, col = m, lty = m)
      graphics::points(x$lambda[m], x$criterion.value, col = m, pch = 19)
    }
    graphics::legend("topright", legend = labs, col = seq_len(M),
                     lty = seq_len(M), bty = "n", cex = 0.8)
  }
  invisible(x)
}
