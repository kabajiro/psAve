# gamma.R -- selection of the prognostic mixing weights gamma (D.6).
#
# gamma minimizes the UNWEIGHTED mean squared prediction error among the
# UNTREATED units:
#   MSE(gamma) = mean_{i: A_i = 0} ( Y_i - sum_k gamma_k g_k(X_i) )^2
# over simplex_grid(K, step), with ties resolved to the first grid row (D.5).
# With family = binomial() the same formula is the Brier score. K = 1 skips
# the grid entirely (gamma = 1).

.select_gamma <- function(y, Gmat, treat, step, verbose = FALSE) {
  K <- ncol(Gmat)
  labels <- colnames(Gmat)
  t0 <- treat == 0L
  G0 <- Gmat[t0, , drop = FALSE]
  y0 <- y[t0]

  ## per-candidate untreated MSE (for the gamma.mse field)
  cand.mse <- colMeans((y0 - G0)^2)
  names(cand.mse) <- labels

  if (K == 1L) {
    gamma <- stats::setNames(1, labels)
    sel.mse <- unname(cand.mse[1L])
    grid.size <- 1L
  } else {
    grid <- simplex_grid(K, step)
    grid.size <- nrow(grid)
    if (grid.size > 1e5) {
      warning(sprintf(paste0("The gamma grid has %s rows (K = %d candidates, step = %s); ",
                             "consider a coarser `step` or fewer prognostic candidates."),
                      format(grid.size, big.mark = ","), K, format(step)),
              call. = FALSE, immediate. = TRUE)
    }
    .vmsg(verbose, sprintf("Selecting gamma over %s grid points (K = %d, step = %s) ...",
                           format(grid.size, big.mark = ","), K, format(step)))
    ## chunked vectorized search: one crossprod-style matrix product per chunk
    vals <- numeric(grid.size)
    for (idx in .chunk_indices(grid.size, 256L)) {
      P <- G0 %*% t(grid[idx, , drop = FALSE])   # n0 x B matrix of gbar values
      vals[idx] <- colMeans((y0 - P)^2)
    }
    sel <- .first_min(vals)
    gamma <- stats::setNames(grid[sel, ], labels)
    sel.mse <- vals[sel]
  }

  gbar <- as.numeric(Gmat %*% gamma)
  gamma.mse <- c(cand.mse, average = sel.mse)

  .vmsg(verbose, sprintf("Selected gamma: %s (untreated MSE = %s)",
                         paste0(labels, " = ", format(gamma, digits = 3),
                                collapse = ", "),
                         .fmt(sel.mse)))

  list(gamma = gamma, gbar = gbar, gamma.mse = gamma.mse, grid.size = grid.size)
}
