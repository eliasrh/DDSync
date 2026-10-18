# References, uncertainties and probe budgets

DDSync solves each station and phase independently. Its theta values are relative to a reference event, and its statistical uncertainty model is conditional on the retained observations and final robust weights. Check both the reference IDs and the uncertainty diagnostics before using the exports.

## References after pruning

Pruning can split an initially connected group. DDSync rebuilds the surviving graph, reapplies the minimum-edge and cycle-redundancy gates, then fits each accepted component independently. The event with the highest surviving degree is pinned to zero; ties use the smallest event ID. In Julia these are the internal sequential IDs when mapping was needed.

Each theta/std row records its own component's reference. A reference star must contain only nodes with that reference. Subtract theta values only within the same station-phase file and finite reference ID. Offsets between separate components are unconstrained. An event that loses all usable observations, or belongs to a rejected fragment, has `NaN` theta, reference and uncertainty and degree zero.

The reference's theta and uncertainty are zero because of the gauge constraint. This does not mean that its physical arrival time is known exactly. For channel-to-channel differential moveout, both the target and chosen anchor must be in the same surviving component at every channel used. Merely finding both IDs in a file is insufficient.

## What the uncertainty calculation estimates

For each accepted component, DDSync forms the weighted graph Laplacian `L`, removes the reference row/column, and adds the configured numerical ridge. Write the resulting matrix as `Lred`. The model standard deviations are

```text
std_theta = sigma_std * sqrt(diag(inv(Lred)))
```

Here `sigma_std` is that component's final residual-noise scale, optionally floored by `min_sigma`. Exact mode computes this diagonal through batched unit-vector solves. Hutchinson mode estimates it with random vectors `z` and solves `Lred*x = z`, averaging the samples `x .* z`.

“Exact” refers to removing Monte Carlo sampling error from this matrix calculation; it does not remove floating-point error, noise-scale uncertainty, selection effects from pruning, model error, or bias from the ridge. Correlation coefficients used as relative weights are not necessarily calibrated inverse observation variances. These exports are not a complete uncertainty budget for a relocated hypocenter or an absolute pick.

Individual Hutchinson samples can be negative even for a valid positive-definite matrix. The averaged diagonal must be finite, strictly positive and sufficiently stable before DDSync accepts it. A negative average is not a small variance.

## Defaults and choosing a probe budget

The default maximum is now **2048 probes**, replacing 100000. DDSync checks convergence after each batch of 128, starting at 128 samples. For every unpinned node it estimates the Monte Carlo standard error of the mean diagonal and requires

```text
standard_error_of_mean / estimated_diagonal <= 0.10
```

All nodes must pass. The default is roughly a 5% one-standard-error diagnostic for the square-root conversion to std, when errors are small. It is **not** a guarantee that every exported std is within 5% or 10% of the exact value, and is not a simultaneous confidence interval. Checking only a median would hide unstable nodes.

Use the following workflow:

1. Start with the defaults. Inspect every `# std_component` line in `sync_metrics.txt` and the actual method used.
2. A budget of 128 is reasonable for an exploratory run only if the checks pass or the reported fallback is acceptable. Setting a smaller budget than `min_probes` always forces fallback. Do not reduce `min_probes` or loosen `rel_tol` simply to label an unstable result as converged.
3. If a larger component reaches its budget with `not_converged` or `nonpositive_estimate`, try 8192 or more probes. For stable sampling behavior, halving Monte Carlo error usually costs about four times as many probes. Graph structure, reference choice and weak connections matter as well as event count.
4. For final analyses, compare representative components at a larger budget or independent random seed. On components within the exact size limit, use `mode="exact"` to remove sampling noise entirely. Assess the accuracy needed by the downstream application.
5. If the reason is `solve_failed` or `invalid_noise_scale`, increasing the probe count does not repair the numerical problem. Inspect input values, weights, weak connections and scale settings. Raising the ridge can bias both theta and uncertainty and should not be a routine way to hide failure.

A study shipped with this update uses ten random seeds per graph and language. A 100-node graph with a strong hub passed at 128 probes in all runs. A 20-node dense graph generally needed about 640. Dense and ring graphs with 100 nodes needed roughly 4000–5200, and a weakly linked graph sometimes failed even at 8192. These are illustrative synthetic cases, not universal thresholds. See the [update report](updates/2026-10-18-pruning-uncertainty.md) for results and reproduction commands.

## Fallback that preserves GraphSplit constraints

The default `fallback="exact"` applies when Hutchinson is disabled (`probes=0`), exceeds `max_nred`, fails a solve, or exhausts its budget without a usable estimate. This includes nonpositive or nonfinite final estimates. Negative samples or intermediate negative averages can recover with more probes; the decision is made when the diagnostic passes or the budget is exhausted.

Exact calculation is limited to **2048 unpinned nodes** (`exact_max_nred`). If the component is larger or the exact solve fails, DDSync uses **`graph_bound`**, a finite conservative bound. It builds a breadth-first spanning tree of the retained weighted component, rooted at its reference, and computes

```text
bound_variance(i) = sigma_std^2 * sum(1 / edge_weight along tree path from reference to i)
```

The sum is a path resistance. Deleting other edges cannot decrease effective resistance, so the tree value bounds the pinned inverse-Laplacian diagonal from above; a nonnegative ridge only lowers that diagonal further. This follows from [Rayleigh monotonicity for weighted graph Laplacians](https://www.control.utoronto.ca/~jwsimpson/papers/2017k/). The bound is conditional on the same noise scale and weighting model as the exact calculation, not a bound on all physical sources of error. It can be loose, especially on dense graphs, and depends on the selected spanning tree. Its cost is linear in nodes plus edges and it requires no factorization or probes.

**Estimator failure does not introduce NaN uncertainty or zero-weight reference constraints for otherwise accepted theta components.** Each non-reference event receives a finite positive bound, and the reference retains its gauge zero. This is necessary for GraphSplit compatibility: its theta reader can replace missing uncertainties with its minimum sigma, which would turn missing information back into apparent confidence. Isolated events and components rejected by the observation gates still have `NaN` theta/reference/std and are not usable constraints in the first place.

The finite bound is explicitly identified as `method=graph_bound`. Use `fallback="graph_bound"` to skip exact fallback, or `mode="graph_bound"` to request the bound directly. `fallback="nan"` and `fallback="none"` are retained as legacy aliases but now select the finite bound after a Hutchinson failure. A disconnected numerical component, invalid noise scale, or numeric overflow that prevents even a finite bound stops the run with an explanation rather than publishing a usable theta with missing uncertainty. A graph-bound fallback needs a positive finite noise scale; for exactly consistent synthetic data, keep `apply_min_sigma=true` and a positive `min_sigma`.

`pseudo_degree` and `pseudo_weight` remain explicit choices for fast practical weighting, including as fallbacks. They are **not calibrated uncertainties or conservative bounds**, especially across weakly connected graphs. They are not the default failure response.

Fallback applies to the whole component so the output does not mix stochastic and proxy meanings within a reference star. The former `min_diag_rel` / MATLAB `min_diagrel` setting is accepted but ignored. No negative variance is replaced by a tiny positive floor.

Explicit `mode="none"` or `export=false` still allows a theta-only workflow. Do not use those settings when the downstream workflow requires supplied uncertainty. Unknown values in this deliberately disabled mode are excluded from median scaling; `thetaStd` output weights are zero, so choose a different output weight mode.

## Configuration examples

Julia TOML, using the new defaults explicitly:

```toml
[std]
mode = "hutch"
probes = 2048
min_probes = 128
batch = 128
rel_tol = 0.10
fallback = "exact"
exact_max_nred = 2048
```

MATLAB equivalents:

```matlab
cfg = ddsync.config_default();
cfg.std.mode = 'hutch';
cfg.std.hutch.probes = 2048;
cfg.std.hutch.min_probes = 128;
cfg.std.hutch.batch = 128;
cfg.std.hutch.rel_tol = 0.10;
cfg.std.fallback = 'exact';
cfg.std.exact_max_nred = 2048;
```

For deterministic uncertainty on small components, set `mode="exact"` / `cfg.std.mode='exact'`. The exact size limit still applies. For a larger stochastic budget change `probes` to 8192; early stopping still applies. Batch size controls memory use and the spacing of convergence checks.

For reproducible investigations, set Julia's RNG with `using Random; Random.seed!(1234)` before calling `DDSync.run(cfg)`, or MATLAB's with `rng(1234)` before `ddsync.run(cfg)`. Julia and MATLAB use different random streams; equal seeds do not imply equal probe vectors. The exact method is the appropriate basis for cross-language numerical comparisons.

Existing configuration files that explicitly set `probes=100000` or `fallback="pseudo_degree"` retain those choices. Copying an updated template is not required; edit only the relevant fields and retain your input/output paths. The new convergence checks also apply to old explicit probe budgets.

## Reading the diagnostics

Each accepted component produces a comment line such as:

```text
# std_component group=TEST_P ref=3 nodes=10 method=exact reason=nonpositive_estimate probes=128 max_rel_se=Inf
```

This means Hutchinson did not provide a valid positive diagonal at its 128-probe budget, but the bounded exact fallback supplied the exported uncertainties. It does not mean the exact calculation failed.

| Field or value | Meaning |
| --- | --- |
| `ref`, `nodes` | Final component reference ID and number of nodes. |
| `method` | Actual output method: `hutch`, `exact`, `graph_bound`, `pseudo_degree`, `pseudo_weight`, `none`, or `disabled`. |
| `reason=requested` | Directly requested method, or uncertainty export disabled. |
| `reason=converged` | Hutchinson passed all per-node diagnostics. |
| `not_converged` | Positive estimates, but at least one has excessive estimated sampling error. |
| `nonpositive_estimate` / `nonfinite_estimate` | Invalid diagonal estimate; no clipping is applied. |
| `insufficient_probes` | Fewer than the required minimum samples. |
| `disabled_probes` / `size_limit` | Hutchinson was skipped. |
| `solve_failed` / `invalid_noise_scale` | Numerical solve or residual-noise scale was unusable. |
| `+exact_size_limit` / `+exact_failed` | Exact fallback was unavailable; the finite `graph_bound` was used. |
| `+legacy_nan_fallback` | A legacy missing-value fallback was replaced by the finite graph bound. |
| `probes` | Number of successful random probe solves used, excluding exact unit-vector solves. |
| `max_rel_se` | Worst estimated relative diagonal MC standard error; `Inf` when invalid/unassessed and `NaN` when no sampling was attempted. |

The group-table `sigma_hat` is a reporting median across final components, not a common scale used for all components. Ignore lines beginning with `#` when parsing the tabular group records. Retain these comments with your analysis for provenance.

## Floors, precision and downstream use

`robust.min_scale` stabilizes pruning and Huber reweighting. `std.min_sigma` floors the residual-noise scale only when `apply_min_sigma=true`; it does not floor each exported row. Neither setting repairs a failed stochastic variance. A valid graph-scaled std can be smaller than `min_sigma`, and a deliberately disabled noise floor can permit zero uncertainty for exactly consistent synthetic data.

Output precision changes decimal formatting only. Temporary data and pass-2 theta/std storage use single precision. Printing nine or twelve decimal places does not establish nanosecond or picosecond accuracy.

Before combining outputs, verify the event-ID mapping, station-phase identity and component reference. Default uncertainty failures now preserve finite bounds for accepted theta rows. Downstream code still needs to skip rows with missing theta/reference, and handle deliberately disabled uncertainty modes explicitly. A user-configured maximum-sigma filter in a downstream program can still remove conservative bounds; DDSync cannot guarantee relocation geometry or override that filter. Rerun DDSync for previously exported datasets affected by these bugs: old theta files alone do not identify the surviving observation graph, and old clipped std files do not reveal the original negative variance estimates.
