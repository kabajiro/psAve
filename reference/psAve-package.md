# psAve: Model-Averaged Propensity Scores Selected by Prognostic-Score Balance

Constructs a model-averaged propensity score as a convex combination of
candidate propensity score models, with mixing weights selected on a
simplex grid to optimize covariate or prognostic-score balance,
implementing the method of Kabata, Stuart and Shintani (2024).
Prognostic scores follow Hansen (2008): outcome models are fit on
untreated units only. The resulting score is designed to be supplied
directly to
[`MatchIt::matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.html)
as a distance measure or to
[`WeightIt::weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.html)
as a propensity score, with balance assessment via cobalt.

## Details

The single estimation function is
[`psave()`](https://kabajiro.github.io/psAve/reference/psave.md). Its
result hands off to the existing ecosystem:
[`psave_match()`](https://kabajiro.github.io/psAve/reference/psave_match.md)
/
[`psave_weight()`](https://kabajiro.github.io/psAve/reference/psave_weight.md)
(or the equivalent explicit
[`MatchIt::matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.html)
/
[`WeightIt::weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.html)
calls),
[`cobalt::bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.html)
(which has a method for `psave` objects), and effect estimation with the
tools you already use (survey, marginaleffects, ...). Utilities
[`simplex_grid()`](https://kabajiro.github.io/psAve/reference/simplex_grid.md)
and
[`psave_criteria()`](https://kabajiro.github.io/psAve/reference/psave_criteria.md)
expose the grid enumeration and the selection criteria for methods
research and testing. See
[psave-details](https://kabajiro.github.io/psAve/reference/psave-details.md)
for the differences from the paper's reference implementation.

## References

Kabata D, Stuart EA, Shintani A (2024). Prognostic score-based model
averaging approach for propensity score estimation. *BMC Medical
Research Methodology*, 24, 228.
[doi:10.1186/s12874-024-02350-y](https://doi.org/10.1186/s12874-024-02350-y)

Hansen BB (2008). The prognostic analogue of the propensity score.
*Biometrika*, 95(2), 481-488.
[doi:10.1093/biomet/asn004](https://doi.org/10.1093/biomet/asn004)

Stuart EA, Lee BK, Leacy FP (2013). Prognostic score-based balance
measures can be a useful diagnostic for propensity score methods in
comparative effectiveness research. *Journal of Clinical Epidemiology*,
66(8), S84-S90.
[doi:10.1016/j.jclinepi.2013.01.013](https://doi.org/10.1016/j.jclinepi.2013.01.013)

Xie Y, Zhu Y, Cotton CA, Wu P (2019). A model averaging approach for
estimating propensity scores by optimizing balance. *Statistical Methods
in Medical Research*, 28(1), 84-101.
[doi:10.1177/0962280217715487](https://doi.org/10.1177/0962280217715487)

## See also

Useful links:

- <https://kabajiro.github.io/psAve/>

- <https://github.com/kabajiro/psAve>

- Report bugs at <https://github.com/kabajiro/psAve/issues>

## Author

**Maintainer**: Daijiro Kabata <daijiro.kabata@port.kobe-u.ac.jp>
\[copyright holder\]

Authors:

- Daijiro Kabata <daijiro.kabata@port.kobe-u.ac.jp> \[copyright holder\]
