# CRAN comments for psAve 1.0.1

## Resubmission

This is a resubmission. The CRAN incoming pre-test (Debian) noted
"Re-building vignettes had CPU time 7.7 times elapsed time". Cause: the
'ranger' and 'xgboost' learners default to using all available cores.
Fixed: both engines now run single-threaded by default inside psAve
(num.threads = 1 / nthread = 1, user-adjustable via the `control`
argument); vignette CPU time is now ~1x elapsed time. Selection results
are unchanged.

## Test environments

* local: macOS (Apple Silicon), R 4.5.0
* GitHub Actions (all passing): ubuntu-latest (R devel, release, oldrel-1),
  windows-latest (R release), macos-latest (R release)
* win-builder: R Under development (unstable) (2026-06-29 r90199 ucrt)

## R CMD check results

Local and GitHub Actions: 0 errors | 0 warnings | 0 notes.
win-builder (R-devel): 0 errors | 0 warnings | 1 note —
"checking CRAN incoming feasibility": New submission; possibly misspelled
words in DESCRIPTION (Kabata, Shintani).

## Comments

* This is a new submission.
* The flagged words (Kabata, Shintani) are author surnames.
* The package implements the published method of Kabata, Stuart and
  Shintani (2024) <doi:10.1186/s12874-024-02350-y>; the maintainer is the
  first author of that paper.
* All Suggests packages ('MatchIt', 'WeightIt', 'SuperLearner', 'rpart',
  'ranger', 'xgboost', 'survey') are used conditionally with
  requireNamespace() guards; core functionality runs with Imports only.
