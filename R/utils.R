# utils.R -- internal utilities: messaging, dependency guards, chunking,
# tie-breaking, formatting. Nothing here is exported.

# Verbose messaging helper (D.13). All progress output goes through here so
# that `verbose = FALSE` (the default) is completely silent.
.vmsg <- function(verbose, ...) {
  if (isTRUE(verbose)) message(...)
  invisible(NULL)
}

# Guard for Suggests packages (D.13): every conditional dependency error names
# the missing package and the exact install.packages() call.
.require_pkg <- function(pkg, purpose) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf(paste0("Package \"%s\" is required %s but is not installed.\n",
                        "Install it with: install.packages(\"%s\")"),
                 pkg, purpose, pkg),
         call. = FALSE)
  }
  invisible(TRUE)
}

# Split 1:n into consecutive chunks of at most `size` indices (D.9).
.chunk_indices <- function(n, size = 256L) {
  if (n <= 0L) return(list())
  starts <- seq.int(1L, n, by = as.integer(size))
  lapply(starts, function(s) seq.int(s, min(s + size - 1L, n)))
}

# First minimum with a small numerical tolerance (D.5). The documented rule is
# "the FIRST grid row attaining the minimum"; the tolerance makes the rule
# robust to floating-point dust so that exact ties (e.g., identical candidate
# columns) resolve to the earlier row, favoring learners listed earlier in
# `ps.methods`/`prog.methods`.
.first_min <- function(x, tol = 1e-9) {
  vmin <- min(x)
  which(x <= vmin + tol * max(1, abs(vmin)))[1L]
}

# Simple text bar for print.psave().
.text_bar <- function(w, width = 20L) {
  filled <- max(0L, min(width, as.integer(round(w * width))))
  paste0("|", strrep("=", filled), strrep(" ", width - filled), "|")
}

# TRUE for columns taking exactly two distinct values (D.2).
.detect_bin <- function(mat) {
  apply(mat, 2L, function(x) length(unique(x)) == 2L)
}

# Validate the s.d.denom argument. The lambda grid search implements these
# four denominators in vectorized form; they are also the ones whose meaning
# is identical in cobalt::col_w_smd().
.validate_sd_denom <- function(s.d.denom) {
  choices <- c("treated", "control", "pooled", "all")
  if (!is.character(s.d.denom) || length(s.d.denom) != 1L ||
      is.na(s.d.denom) || !s.d.denom %in% choices) {
    stop(sprintf("`s.d.denom` must be one of %s.",
                 paste0('"', choices, '"', collapse = ", ")),
         call. = FALSE)
  }
  s.d.denom
}

# Human-readable one-liner for each criterion (used by print methods).
.criterion_label <- function(criterion, prog.target = NA_character_) {
  switch(criterion,
         prog = if (identical(prog.target, "average"))
           "weighted ASMD of the model-averaged prognostic score"
         else
           sprintf("weighted ASMD of the \"%s\" prognostic score", prog.target),
         smd = "mean weighted ASMD of the covariates",
         ks = "mean weighted KS statistic of the covariates",
         logloss = "treatment-assignment log loss",
         criterion)
}

# Format a numeric scalar for messages.
.fmt <- function(x, digits = 4) format(x, digits = digits)
