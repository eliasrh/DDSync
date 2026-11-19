# Per-event uncertainty and selected exact repair — 12 November 2026

**Historical defaults:** the [19 November reporting update](2026-11-19-uncertainty-reporting.md) raises the default maximum to 4096 probes and reorganizes metrics. The 2048-probe synthetic results below remain the actual measurements from this study. Selected-repair and fallback behavior is unchanged.

The October graph bound was mathematically conservative but a poor default for practical weights. Its whole-component trigger also discarded useful Hutchinson estimates whenever even one entry failed the convergence check. This update replaces that policy in Julia and fixes the same behavior in MATLAB. Post-pruning references and graph membership are unchanged.

## What changed

- **`mode="auto"` is the default.** Components with at most 2048 unpinned events go directly to exact inverse-diagonal calculation. Larger components use Hutchinson, still capped at 2048 probes by default.
- **Convergence is acted on per event.** Every passing entry is retained. Failure elsewhere does not replace it or change its value.
- **Selected exact repairs reuse the Hutch factorization.** Up to 128 failing entries receive full-system unit-vector solves. Invalid estimates are prioritized, then the largest relative Monte Carlo errors. The repair budget counts right-hand sides, independently of component size. This is not a neighborhood inverse: external weak connections remain in the calculation.
- **Remaining positive estimates may be retained.** The default `unconverged="keep"` retains estimates with enough samples and relative Monte Carlo standard error below 1, while labeling them `hutch_unconverged`. They have not passed the requested 0.10 tolerance. Nonpositive/nonfinite estimates and near-zero estimates with SE at least as large as the estimate use fallback. `unconverged="fallback"` instead replaces all unrepaired failing entries.
- **The final default fallback is `pseudo_degree`.** It affects only unresolved entries. Graph bounds are explicitly selectable, but no default route selects them. Legacy missing-value fallback aliases now use degree fallback. An explicit `fallback="exact"` repairs all failures within the full exact cap, then uses the selected-repair/degree policy above that cap.
- **Quality is visible per event.** New `std_quality_STA_PH.txt` files give event ID, reference ID, actual method and pre-repair relative MC standard error. Component metrics report the numbers converged, exactly calculated, retained unconverged and assigned fallback. Existing numeric theta/std column layouts are unchanged.

The solve check now uses componentwise relative backward error, with protection against numerical underflow. The old RHS-only residual threshold could reject accurate solves when large terms cancelled. The definition follows [LAPACK's componentwise backward error](https://www.netlib.org/lapack/explore-html/d0/df7/dla__lin__berr_8f_source.html); a small backward error still does not establish small forward error in an ill-conditioned system.

## Large-component synthetic results

The study uses unit residual-noise scale, a ridge of `1e-10`, 2048 Rademacher probes, a 0.10 relative-diagonal-SE tolerance and at most 128 selected repairs. It runs three seeds in each language. All methods are compared to a full exact diagonal; the complete-graph result also has an analytic oracle. The tests deliberately force Hutchinson to expose its behavior. The graphs are larger than the default full exact cap.

| Synthetic graph | Julia and MATLAB result |
| --- | --- |
| **2,400-event complete graph**, unit edge weights | The tree graph bound was **34.641 times exact std** at every non-reference event. That gives only **0.0833% of the exact inverse-variance weight**, or 2.89% of its inverse-std weight, before caps. |
| **5,000 events, mostly well constrained**, with a 40-event clique connected through a weight-0.001 edge | **4,958 entries converged and were preserved; 41 were repaired exactly**, in every seed and both languages. No proxy entries were needed. Across seeds, maximum relative std error was 3.12% in Julia and 6.21% in MATLAB; the 95th-percentile error was below 1.30%. The old all-node rule rejected the entire component. |
| **5,000 events in two hub-and-neighbor groups**, joined by one edge | **2,499 entries converged**, 128 were repaired, and many remained noisy. Default retention preserved roughly 1,500 additional estimates; roughly 800 still needed degree fallback. The 95th-percentile std error remained about **52%**. The worst relative error reached about **106%** across the six runs. This is a case requiring more work before treating all output as accurate uncertainty. |

For a complete graph, the pinned inverse diagonal is `2/n` without the negligible ridge, while each BFS path has resistance 1. Thus the graph-bound/std ratio is `sqrt(n/2)`. The problem is not limited to a stochastic failure: this bound is intrinsically loose on dense graphs.

On the complete graph, no Hutch entry met the 10% tolerance at 2048 probes. After 128 repairs, the default retained 1,913–1,943 positive noisy estimates and used degree fallback for 328–358 entries. Median std/exact ratios were 1.000–1.018, but the 95th-percentile absolute relative error was **35–42%**, with worst errors up to **73%**. Preserving estimates avoids the massive graph-bound inflation; it does not make an undersampled estimator exact. For a component just above the full exact cap, raising `exact_max_nred` modestly can be more sensible than spending many more random probes.

Degree fallback is not universally accurate either. On the complete graph it gives about **0.707 times exact std**, doubling inverse-variance weight. Across the single weak connection, a high-degree remote hub could receive about **0.02 times exact std**, giving roughly 2,500 times the exact inverse-variance weight if the proxy is used there. Degree is therefore a pragmatic final fallback, not a replacement for graph-aware uncertainty. Counts and per-event labels are necessary for judging the result.

Raw results: [Julia CSV](data/2026-11-12-julia.csv), [MATLAB CSV](data/2026-11-12-matlab.csv). `keep` is the default retention policy; `fallback` replaces all unrepaired failures with degree proxies. Errors exclude the pinned reference. Timing columns are illustrative only: first-call compilation, differing threading and concurrent local workloads prevent a controlled speed comparison.

## Validation

Tests were run with Julia 1.12.3 and native MATLAB R2025b on macOS.

| Check | Result |
| --- | --- |
| Julia regressions | **1,342 assertions passed**, including independent covariance oracles, per-entry retention, selected exact repair, repair budgets, exact-first selection, proxy fallback, old options, input-ID mapping and exports. |
| Native MATLAB regressions | **604 checks passed**, including equivalent new cases. |
| Actual GraphSplit 0.4.0 reader and star builder | **603 assertions passed** with forced degree/exact-option/graph-bound/legacy fallbacks, retained unconverged entries and MATLAB exports. Each fixture retained all 36 expected stars with finite sigma and exactly the intended two component gauges. |
| Full bundled dataset, new defaults | Both implementations retained **608205 of 608577 observations** across 89 station-phase groups. All **70 accepted components used exact calculation directly**, with no probes or fallback. |
| Cross-language comparison | All **14086 finite theta/std rows** in the bundled dataset matched reference IDs, theta and std exactly at output precision. Synchronized differential times also matched. Separate synthetic exact and forced-degree outputs matched numerically. |

Regression tests verify that a selected repair includes a remote weak bridge, rather than using a local conditional inverse. They also check that one bad selected solve does not discard other valid columns, negative estimates do not invalidate neighboring positive entries, tiny underflow terms do not reject valid solves, and deliberately skipped Hutch can still use explicitly requested exact fallback within its cap.

GraphSplit tests check the theta-star interface and weighted incidence rank, not whether arbitrary station geometry constrains every relocation coordinate. Accepted components retain finite uncertainties under the new default policy. Rejected/isolated theta rows remain missing; explicit theta-only modes remain available. Invalid noise scales or unrepresentable uncertainty stop explicitly instead of introducing missing sigma into usable constraints.

## Reproduce and adopt

From the repository root:

```bash
julia --startup-file=no DDSyncJulia/test/runtests.jl
julia --startup-file=no DDSyncJulia/test/large_uncertainty_study.jl > julia-large.csv
julia --startup-file=no DDSyncJulia/test/graphsplit_integration.jl /path/to/GraphSplit.jl
```

In MATLAB:

```matlab
addpath('DDSync/test');
run_tests;
large_uncertainty_study('matlab-large.csv');
```

To include MATLAB outputs in the GraphSplit check, use `run_tests('/absolute/path/to/matlab-results')` and pass that directory as the integration script's second argument. The large studies are separate from the quick CI suites; the dense graph has 2,878,800 edges and needs correspondingly more memory. GraphSplit remains an optional test dependency.

Existing explicit configuration values remain in effect. Set `mode="auto"` and `fallback="pseudo_degree"` to adopt the new selection and final fallback; `repair_max_nodes=128` and `unconverged="keep"` are supplied by defaults. Explicit `mode="hutch"` still samples small components. Explicit `fallback="graph_bound"` remains opt-in. The [uncertainty guide](../uncertainty.md) and both manuals explain all fields and the per-event quality files.

For a few failing entries, raise the repair budget. For broadly noisy components, consider more probes or full exact calculation where affordable. Do not loosen the tolerance just to relabel uncertainty as converged, and do not interpret a degree proxy as calibrated confidence. Previously exported graph-bound files need to be regenerated; they do not preserve the stochastic estimates that were discarded.
