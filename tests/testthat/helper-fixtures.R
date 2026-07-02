# ==========================================================================
# Shared fixtures for the psAve test suite
# ==========================================================================
#
# (a) make_sim_data(): deterministic simulated dataset mirroring the layout
#     of the paper's DGP-A (10 covariates, binary treatment, continuous
#     outcome).  Seed is fixed INSIDE the helper, so every test sees the
#     same data.  With seed = 20240228 and n = 200 the realized arms are
#     84 treated / 116 untreated and the in-sample glm propensity scores
#     lie in (0.073, 0.962), i.e. strictly inside the default clipping
#     bounds [0.01, 0.99] (verified independently of the package).
#
# (b) make_fixture(): a 6-row HAND-COMPUTABLE fixture.  Every expected
#     number used in the exact-value tests is derived below; nothing was
#     copied from running the package.
#
# --------------------------------------------------------------------------
# 6-ROW FIXTURE — FULL HAND DERIVATION
# --------------------------------------------------------------------------
# Units u1..u6:
#   A  = (1, 1, 1, 0, 0, 0)           (u1-u3 treated, u4-u6 untreated)
#   x  = (1, 2, 3, 2, 3, 4)           (single balance covariate)
#   y  = (2, 3, 4, 1, 2, 4)           (outcome; untreated y0 = (1, 2, 4))
#
# Candidate propensity scores (all within [0.01, 0.99], no clipping):
#   e1 = (.50, .50, .50, .20, .50, .80)
#   e2 = (.40, .50, .60, .80, .50, .20)
#   e3 = (.50, .50, .50, .30, .50, .70)   (used only in the vertex test)
#
# Candidate prognostic scores:
#   g  = (1, 2, 3, 1, 2, 4)   sole column of G1; also column "g2" of G2
#   g1 = (1, 2, 3, 2, 2, 2)   column "g1" of G2 (control-constant)
#
# Treated-group statistics (denominators; H.1: plain n-1 sample SD):
#   mean(x | A=1) = 2,  sd(x | A=1) = sd(1,2,3) = 1
#   mean(g | A=1) = 2,  sd(g | A=1) = 1;   mean(g1|A=1) = 2, sd = 1
#
# simplex_grid(2, 0.5) rows, IN ORDER:  (1,0), (.5,.5), (0,1)
# Averaged PS at the three rows:
#   row1: ebar = e1
#   row2: ebar = (e1+e2)/2 = (.45, .50, .55, .50, .50, .50)
#   row3: ebar = e2
#
# --- ATT weights (W = 1 treated; ebar/(1-ebar) control, units u4,u5,u6) ---
#   row1: (.2/.8, .5/.5, .8/.2) = (1/4, 1, 4),      sum = 21/4
#   row2: (1, 1, 1),                                 sum = 3
#   row3: (4, 1, 1/4),                               sum = 21/4
#
# prog criterion, ATT  (g control values (1, 2, 4); denominator 1):
#   row1: wm0 = (1/4*1 + 1*2 + 4*4)/(21/4) = (73/4)/(21/4) = 73/21
#         |2 - 73/21| = 31/21 = 1.476190476...
#   row2: wm0 = (1+2+4)/3 = 7/3;   |2 - 7/3| = 1/3
#   row3: wm0 = (4*1 + 1*2 + 1/4*4)/(21/4) = 7/(21/4) = 4/3; |2-4/3| = 2/3
#   => argmin = row 2: lambda = (.5, .5), value = 1/3
#
# smd criterion, ATT  (x control values (2, 3, 4); denominator 1):
#   row1: wm0 = (1/4*2 + 1*3 + 4*4)/(21/4) = (39/2)/(21/4) = 26/7
#         |2 - 26/7| = 12/7
#   row2: wm0 = 3;  |2 - 3| = 1
#   row3: wm0 = (4*2 + 1*3 + 1/4*4)/(21/4) = 12/(21/4) = 16/7; |2-16/7| = 2/7
#   => argmin = row 3: lambda = (0, 1), value = 2/7
#
# weighted KS of x at ps = e1, ATT (weighted eCDFs over observed x):
#   F1: x=1: 1/3;  x=2: 2/3;  x=3: 1
#   F0 (weights (1/4,1,4) at x=(2,3,4), total 21/4):
#       F0(2) = (1/4)/(21/4) = 1/21;  F0(3) = (5/4)/(21/4) = 5/21;  F0(4)=1
#   |F1-F0|: x=1: 1/3;  x=2: 14/21-1/21 = 13/21;  x=3: 1-5/21 = 16/21; x=4: 0
#   => KS = 16/21
#
# logloss (weight-free), the three grid rows:
#   row1: -(3*log(.5) + log(.8) + log(.5) + log(.2))/6
#       = -(4*log(.5) + log(.8) + log(.2))/6            = 0.7675283643...
#   row2: -(log(.45) + log(.5) + log(.55) + 3*log(.5))/6
#       = -(log(.45) + log(.55) + 4*log(.5))/6          = 0.6948222365...
#   row3: -(log(.4) + log(.5) + log(.6) + log(.2) + log(.5) + log(.8))/6
#                                                        = 0.7743320301...
#   => argmin = row 2: lambda = (.5, .5)
#
# --- ATE weights (1/ebar treated; 1/(1-ebar) control) ---
#   row1: treated (2, 2, 2);            control (5/4, 2, 5),  sum 33/4
#   row2: treated (20/9, 2, 20/11), sum 598/99;  control (2, 2, 2)
#   row3: treated (5/2, 2, 5/3),  sum 37/6;      control (5, 2, 5/4), sum 33/4
#
# prog criterion, ATE (denominator still sd(g|A=1) = 1 for BOTH estimands):
#   row1: wm1 = 2 (equal weights);  wm0 = (5/4 + 4 + 20)/(33/4) = 101/33
#         |2 - 101/33| = 35/33 = 1.0606060606...
#   row2: wm1 = (20/9*1 + 2*2 + 20/11*3)/(598/99) = (1156/99)/(598/99)
#             = 578/299;  wm0 = 7/3
#         |578/299 - 7/3| = |(1734-2093)/897| = 359/897 = 0.4002229654...
#   row3: wm1 = (5/2 + 4 + 5)/(37/6) = (23/2)/(37/6) = 69/37
#         wm0 = (5*1 + 2*2 + 5/4*4)/(33/4) = 14/(33/4) = 56/33
#         |69/37 - 56/33| = |(2277-2072)/1221| = 205/1221 = 0.1678951679...
#   => argmin = row 3: lambda = (0, 1), value = 205/1221
#
# smd of x at ps = e1, ATE:
#   wm1 = 2;  wm0 = (5/4*2 + 2*3 + 5*4)/(33/4) = (57/2)/(33/4) = 38/11
#   |2 - 38/11| = 16/11
# weighted KS of x at ps = e1, ATE:
#   F1 (equal weights): 1/3, 2/3, 1 at x = 1, 2, 3
#   F0 (weights (5/4,2,5) at x=(2,3,4), total 33/4):
#       F0(2) = 5/33;  F0(3) = 13/33;  F0(4) = 1
#   diffs: 1/3 = 11/33; 22/33-5/33 = 17/33; 33/33-13/33 = 20/33  => KS = 20/33
#
# gamma selection (G2, step = .5; grid rows (1,0), (.5,.5), (0,1)),
# unweighted MSE among untreated, y0 = (1, 2, 4):
#   g1_0 = (2,2,2):  MSE = ((1-2)^2 + 0 + (4-2)^2)/3 = 5/3
#   mix  = (1.5,2,3): MSE = ((1-1.5)^2 + 0 + (4-3)^2)/3 = (5/4)/3 = 5/12
#   g2_0 = (1,2,4):  MSE = 0
#   => gamma = (0, 1); gamma.mse = (g1 = 5/3, g2 = 0, average = 0)
#
# prog.target = "g1": g1 is constant (= 2) on the control arm, so its
#   weighted control mean is exactly 2 = treated mean for EVERY lambda;
#   the criterion is identically 0 and the FIRST grid row wins:
#   lambda = (1, 0), value = 0.
#
# average = FALSE with (e1, e2, e3), prog criterion, ATT (vertex values):
#   e1: 31/21 (above)
#   e2: 2/3   (above)
#   e3: control weights (.3/.7, .5/.5, .7/.3) = (3/7, 1, 7/3), sum 79/21
#       wm0 = (3/7*1 + 1*2 + 7/3*4)/(79/21) = (247/21)/(79/21) = 247/79
#       |2 - 247/79| = 89/79 = 1.126582...
#   => best vertex = e2: lambda = (0, 1, 0), value = 2/3
# ==========================================================================

make_fixture <- function() {
  d <- data.frame(
    A = c(1, 1, 1, 0, 0, 0),
    x = c(1, 2, 3, 2, 3, 4),
    y = c(2, 3, 4, 1, 2, 4),
    row.names = paste0("u", 1:6)
  )
  E <- cbind(
    e1 = c(0.50, 0.50, 0.50, 0.20, 0.50, 0.80),
    e2 = c(0.40, 0.50, 0.60, 0.80, 0.50, 0.20)
  )
  E3 <- cbind(E, e3 = c(0.50, 0.50, 0.50, 0.30, 0.50, 0.70))
  G1 <- cbind(g = c(1, 2, 3, 1, 2, 4))
  G2 <- cbind(
    g1 = c(1, 2, 3, 2, 2, 2),
    g2 = c(1, 2, 3, 1, 2, 4)
  )
  list(data = d, E = E, E3 = E3, G1 = G1, G2 = G2)
}

make_sim_data <- function(n = 200, seed = 20240228) {
  set.seed(seed)
  X <- matrix(rnorm(n * 7), n, 7)
  B <- matrix(rbinom(n * 3, 1L, 0.5), n, 3)
  lp <- 0.4 * X[, 1] - 0.4 * X[, 2] + 0.3 * X[, 3] + 0.3 * X[, 4] -
    0.2 * X[, 5] + 0.5 * B[, 1]
  A <- rbinom(n, 1L, plogis(lp))
  y <- 1 + 0.5 * A + 0.6 * X[, 1] + 0.6 * X[, 2] + 0.4 * X[, 5] +
    0.4 * X[, 6] + 0.5 * B[, 2] + rnorm(n)
  d <- data.frame(A = A, y = y, X, B)
  names(d) <- c("A", "y", paste0("X", 1:10))
  d
}

sim_formula <- function() {
  stats::as.formula(paste("A ~", paste(paste0("X", 1:10), collapse = " + ")))
}

# --------------------------------------------------------------------------
# From-scratch reference implementations of the criteria (base R only).
# These deliberately re-derive the formulas from the paper/supplement and
# never call psAve or cobalt.
# --------------------------------------------------------------------------

ref_weights <- function(ps, treat, estimand) {
  if (estimand == "ATT") ifelse(treat == 1, 1, ps / (1 - ps))
  else                   ifelse(treat == 1, 1 / ps, 1 / (1 - ps))
}

# |weighted treated mean - weighted control mean| / sd(x | treat == 1)
# (H.1: plain sample SD for ALL covariates, including binary ones)
ref_asmd <- function(x, treat, w) {
  wm1 <- sum(w[treat == 1] * x[treat == 1]) / sum(w[treat == 1])
  wm0 <- sum(w[treat == 0] * x[treat == 0]) / sum(w[treat == 0])
  abs(wm1 - wm0) / stats::sd(x[treat == 1])
}

# max over observed x of |F1_w(x) - F0_w(x)| with proper weighted eCDFs
ref_ks <- function(x, treat, w) {
  xs <- sort(unique(x))
  F1 <- vapply(xs, function(z)
    sum(w[treat == 1] * (x[treat == 1] <= z)) / sum(w[treat == 1]), 0)
  F0 <- vapply(xs, function(z)
    sum(w[treat == 0] * (x[treat == 0] <= z)) / sum(w[treat == 0]), 0)
  max(abs(F1 - F0))
}

ref_logloss <- function(ps, treat) {
  -mean(treat * log(ps) + (1 - treat) * log(1 - ps))
}
