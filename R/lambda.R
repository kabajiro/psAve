# lambda.R -- the lambda grid search (D.8, D.9): chunked vectorized
# evaluation of the selection criterion at every grid row, plus the
# average = FALSE vertex mode. Nothing here is exported.
#
# Criterion definitions (evaluated at each grid row lambda, with
# ebar = E %*% lambda and W the D.7 weights):
#   logloss: -mean(A*log(ebar) + (1-A)*log(1-ebar))    (finite by clipping)
#   smd:     mean over covariates j of |wmean1_j - wmean0_j| / s_j with s_j
#            the UNWEIGHTED sample SD of the s.d.denom group (uniform
#            sample-SD standardization for ALL columns, incl. binary -- the
#            paper's convention, H.1)
#   ks:      mean over covariates of the proper weighted-eCDF KS statistic
#   prog:    the smd formula applied to the single prognostic-score column
#
# Vectorization (D.9): logloss/smd/prog use chunked matrix algebra
# (Ebar = E %*% t(grid chunk), ~256 rows per chunk; weighted means are one
# crossprod(X0, W0) per chunk; under the ATT the treated side has W = 1 and
# is precomputed once). ks loops over grid rows with per-covariate sort
# orders precomputed once.

# Unweighted SD denominators per column (sample SD, n - 1 formula), matching
# cobalt::col_w_smd() with bin.vars = FALSE and unweighted standardization.
.sd_denom_vec <- function(mat, treat, s.d.denom) {
  v1 <- apply(mat[treat == 1L, , drop = FALSE], 2L, stats::var)
  v0 <- apply(mat[treat == 0L, , drop = FALSE], 2L, stats::var)
  s2 <- switch(s.d.denom,
               treated = v1,
               control = v0,
               pooled = (v1 + v0) / 2,
               all = apply(mat, 2L, stats::var))
  sqrt(s2)
}

# Precompute per-covariate sort orders and distinct-value step positions for
# the weighted-eCDF KS statistic (D.9).
.ks_prep <- function(covs) {
  lapply(seq_len(ncol(covs)), function(j) {
    x <- covs[, j]
    ord <- order(x)
    xs <- x[ord]
    nx <- length(xs)
    keep <- which(xs[-nx] != xs[-1L])   # eCDF evaluation points (distinct x)
    list(ord = ord, keep = keep)
  })
}

# Mean weighted-eCDF KS over covariates for signed normalized weights sw
# (treated: +W/sum(W_treated); untreated: -W/sum(W_untreated)); the running
# cumulative sum of sw in covariate order is F1 - F0 at each observed value.
.ks_mean <- function(sw, prep) {
  vals <- vapply(prep, function(pj) {
    if (!length(pj$keep)) return(0)
    cs <- cumsum(sw[pj$ord])
    max(abs(cs[pj$keep]))
  }, numeric(1L))
  mean(vals)
}

# The lambda grid search. E is the n x M matrix of CLIPPED candidate PS with
# candidate labels as colnames; prog.vec is the prognostic-score target
# column (only used when criterion == "prog").
.select_lambda <- function(E, treat, covs, prog.vec, criterion, estimand,
                           s.d.denom, step, average = TRUE, keep.path = TRUE,
                           verbose = FALSE) {
  M <- ncol(E)
  n <- nrow(E)
  labels <- colnames(E)

  ## --- grid -----------------------------------------------------------------
  if (average) {
    n_steps <- round(1 / step)
    G <- choose(n_steps + M - 1, M - 1)
    if (G > 1e5) {
      warning(sprintf(paste0("The lambda grid has %s rows (M = %d candidates, step = %s); ",
                             "consider a coarser `step` or fewer candidates. ",
                             "The criterion path will not be stored."),
                      format(G, big.mark = ","), M, format(step)),
              call. = FALSE, immediate. = TRUE)
    }
    grid <- simplex_grid(M, step)
  } else {
    ## vertex mode (D.8): the M simplex vertices in candidate order = the
    ## best SINGLE candidate by `criterion`
    grid <- diag(M)
  }
  G <- nrow(grid)
  .vmsg(verbose, sprintf("Selecting lambda over %s grid %s (M = %d, criterion = \"%s\") ...",
                         format(G, big.mark = ","),
                         if (average) "points" else "vertices (average = FALSE)",
                         M, criterion))

  ## --- precomputations --------------------------------------------------------
  t1 <- treat == 1L
  E1 <- E[t1, , drop = FALSE]
  E0 <- E[!t1, , drop = FALSE]

  if (criterion %in% c("smd", "prog")) {
    if (criterion == "smd") {
      Xall <- covs
    } else {
      Xall <- matrix(prog.vec, ncol = 1L,
                     dimnames = list(NULL, "prog"))
    }
    s <- .sd_denom_vec(Xall, treat, s.d.denom)
    if (any(s == 0) || anyNA(s)) {
      zero <- colnames(Xall)[s == 0 | is.na(s)]
      stop(sprintf(paste0("Zero or undefined standard deviation in the \"%s\" group for: %s. ",
                          "The ASMD-based criterion is undefined; remove constant columns or ",
                          "choose a different `s.d.denom` or `criterion`."),
                   s.d.denom, paste0('"', zero, '"', collapse = ", ")),
           call. = FALSE)
    }
    X1 <- Xall[t1, , drop = FALSE]
    X0 <- Xall[!t1, , drop = FALSE]
    mu1 <- colMeans(X1)   # treated side is constant in lambda under the ATT (W = 1)
  } else if (criterion == "ks") {
    prep <- .ks_prep(covs)
  }

  ## --- chunked evaluation -----------------------------------------------------
  values <- numeric(G)
  for (idx in .chunk_indices(G, 256L)) {
    L <- grid[idx, , drop = FALSE]                    # B x M

    if (criterion == "logloss") {
      Ebar <- E %*% t(L)                              # n x B
      values[idx] <- -(as.numeric(crossprod(treat, log(Ebar))) +
                       as.numeric(crossprod(1 - treat, log1p(-Ebar)))) / n

    } else if (criterion %in% c("smd", "prog")) {
      Ebar0 <- E0 %*% t(L)                            # n0 x B
      W0 <- if (estimand == "ATT") Ebar0 / (1 - Ebar0) else 1 / (1 - Ebar0)
      m0 <- crossprod(X0, W0)                         # p x B
      m0 <- sweep(m0, 2L, colSums(W0), "/")
      if (estimand == "ATT") {
        d <- abs(mu1 - m0)                            # mu1 recycles down columns
      } else {
        Ebar1 <- E1 %*% t(L)
        W1 <- 1 / Ebar1
        m1 <- sweep(crossprod(X1, W1), 2L, colSums(W1), "/")
        d <- abs(m1 - m0)
      }
      values[idx] <- colMeans(d / s)                  # s recycles down columns

    } else {                                          # ks: per-row loop
      Ebar <- E %*% t(L)
      for (b in seq_along(idx)) {
        W <- .ipw_from_ps(Ebar[, b], treat, estimand)
        sw <- numeric(n)
        sw[t1] <- W[t1] / sum(W[t1])
        sw[!t1] <- -W[!t1] / sum(W[!t1])
        values[idx[b]] <- .ks_mean(sw, prep)
      }
    }
  }

  ## --- selection (first minimum = documented tie-break, D.5) -----------------
  sel <- .first_min(values)
  lambda <- stats::setNames(grid[sel, ], labels)
  value <- values[sel]

  .vmsg(verbose, sprintf("Selected lambda: %s (criterion value = %s)",
                         paste0(labels, " = ", format(lambda, digits = 3),
                                collapse = ", "),
                         .fmt(value)))

  ## --- path (D.9) -------------------------------------------------------------
  path <- NULL
  if (G > 1e5) {
    message("Criterion path not stored: the lambda grid exceeds 100,000 rows.")
  } else if (isTRUE(keep.path)) {
    path <- as.data.frame(grid)
    names(path) <- labels
    path$value <- values
    path$selected <- seq_len(G) == sel
  }

  list(lambda = lambda, value = value, path = path, grid.size = G)
}
