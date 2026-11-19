# References, uncertainties and probe budgets

DDSync fits each station-phase component independently. Theta is relative to a reference event; uncertainty is conditional on the retained observations, final robust weights and residual-noise model. The default now uses exact calculation for small components and makes uncertainty decisions **per event** for larger ones. A few unstable Hutchinson estimates no longer replace a whole component with a graph bound.

## References after pruning

DDSync rebuilds the graph after pruning and reapplies the minimum-edge and cycle-redundancy gates. Each accepted component gets its own highest-degree reference, with smallest event ID breaking ties. In Julia, mapped datasets use internal sequential IDs here.

Subtract theta values only within the same station-phase file and reference ID. Offsets between separate components are unconstrained. Isolated events and rejected fragments have missing theta/reference/std and degree zero. Each accepted reference has theta and std zero because it fixes the gauge; its physical arrival time is not known exactly.

## What is being estimated

For each retained event pair, DDSync compares the measured differential travel time `dt_ij` with the difference predicted by the fitted theta values. The **residual** is the mismatch:

```text
r_ij = (theta_i - theta_j) - dt_ij
```

`sigma_hat` measures how widely these mismatches vary after the final fit, in seconds. DDSync estimates it separately for each connected group at one station and phase and uses it as an estimate of noise in the pairwise measurements. Specifically, it uses 1.4826 times the median absolute deviation from the median residual, using observations whose final positive weight is at least half the median positive weight. If that estimate is zero or nonfinite, it uses the ordinary sample standard deviation of residuals with positive weights.

**`min_sigma = 5e-4` means assuming at least 0.5 ms scatter in these travel-time mismatches when calculating uncertainty.** `sigma_std` is the scatter actually used:

```text
sigma_std = max(sigma_hat, min_sigma)   # when apply_min_sigma=true and min_sigma>0
sigma_std = sigma_hat                   # otherwise
```

DDSync then combines this scatter with the weights and arrangement of retained observations to calculate each theta uncertainty. Let `Lred` be the weighted graph Laplacian after removing the reference row and column and adding the numerical ridge. The conditional model is

```text
std_theta = sigma_std * sqrt(diag(inv(Lred)))
```

**These settings do not impose a minimum `std_theta`.** Combining observations can give smaller theta uncertainties. For example, if the measured mismatch scatter is `1e-5` seconds and the minimum is enabled at `5e-4` seconds, then `sigma_std = 5e-4` seconds. The degree proxy at degree 100 gives `std_theta = 5e-4 / sqrt(100) = 5e-5` seconds. This distinction also applies to exact and Hutchinson calculations, not just fallbacks. The earlier template wording “std_theta floor” was incorrect.

The shorter `sigma` used in descriptions such as `sigma / sqrt(degree)` means `sigma_std`; it is not a different fitted noise parameter. Neither value is an event's standard deviation or the standard error of Hutchinson sampling. The report lists both scales per component. In the observation table, `sigma_hat` is the median of the raw component scales in that station-phase group, solely for reporting. Each component still uses its own scale. The `pseudo_weight` formula does not use either scale.

Exact calculation removes Monte Carlo sampling error from the inverse-diagonal calculation. It does not remove floating-point error, ridge bias, uncertainty in the noise scale, pruning effects or physical model error. Correlation coefficients used as relative weights are not necessarily calibrated inverse observation variances.

Hutchinson uses random vectors `z`, solves `Lred*x = z`, and averages `x .* z`. Individual samples and finite-sample averages can be negative even though the true inverse diagonal is positive. Negative estimates are never replaced with a tiny positive variance.

### Squared Cholesky draws with a shared timing shift

Consider a cluster whose events are well timed relative to one another. If only a few observations connect the cluster to the reference, its times can shift together. DDSync can calculate one such shared contribution exactly, then sample the remaining variance.

For the reduced timing matrix `A`, write `C = inv(A)`. One solve gives `v = A \ ones`, with `s = sum(v)`. The contribution of the shared shift to event `i` is `v[i]^2 / s`. DDSync samples a timing vector `y` with covariance `C` using one triangular solve per draw. It subtracts `v * sum(y) / s`, squares each remaining entry, averages across draws and adds the shared contribution back. The covariance of the remaining draws is `C - v*v'/s`. The calculation uses vectors and a sparse factor. It never needs to form either dense covariance matrix.

The per-event precision check uses the variation of the squared remainder divided by the estimated **total** variance. Taking the square root and multiplying by `sigma_std` gives the reported uncertainty in seconds. See the [worked example, settings and memory limits](updates/2026-11-19-cholesky-uncertainties.md).

## Default decision sequence

1. **Small component:** `mode="auto"` computes the exact diagonal directly when there are at most `exact_max_nred=2048` unpinned events. It spends no random probes first.
2. **Larger component:** sample up to `probes=4096`, using the sampler selected by `auto_probe`. Hutchinson checks after each batch once it has at least 128 samples. Cholesky first checks at 256 draws by default and doubles the budget only while needed. Sampling stops early if every entry passes. At the end, preserve every entry whose estimated Monte Carlo standard error divided by its estimated diagonal is at most `rel_tol=0.10`.
3. **Selected exact repair:** solve for up to `repair_max_nodes=128` failed entries, prioritizing invalid estimates, then the largest relative sampling errors. Reuse the factorization already used by the sampler. This limit counts extra right-hand sides, independently of component size.
4. **Remaining positive estimates:** with `unconverged="keep"`, retain entries with sufficient samples and finite positive estimates whose relative MC standard error is below 1. They are labeled `hutch_unconverged` or `cholesky_unconverged`. The cutoff excludes near-zero estimates whose standard error is as large as the estimate; it is not a confidence guarantee. Use `unconverged="fallback"` to send all unrepaired failing entries to fallback instead.
5. **Remaining unusable entries:** use `fallback="pseudo_degree"`, a finite practical proxy `sigma_std / sqrt(degree)`. This also applies when sampling is disabled, exceeds `max_nred=40000`, or cannot obtain usable solves. Successful entries elsewhere in the component are preserved.

The 10% diagnostic applies to the **variance diagonal**. In the small-error regime it corresponds to about 5% relative standard error of std. Neither this diagnostic nor the 100% cutoff for retained noisy entries is a guaranteed error bound or simultaneous confidence interval. A positive estimate can still be inaccurate. Inspect the per-event method labels and do not treat mixed output as uniformly converged.

A selected exact repair solves `Lred*x = e_i` and takes `x[i]`. This uses the **full component**; it is not an inverse of a local neighborhood. A neighborhood calculation would condition away external connections and can miss uncertainty caused by a weak bridge. Extra work is bounded by the number of selected solves and the batch size, although sparse-factorization memory still depends on graph structure.

Linear solves are checked using componentwise relative backward error, `abs(Lred*x-b) / (abs(Lred)*abs(x)+abs(b))`, with tolerance `1e-12` and protection against underflow in the residual denominator. This is the [LAPACK backward-error definition](https://www.netlib.org/lapack/explore-html/d0/df7/dla__lin__berr_8f_source.html). It avoids rejecting an accurate solve just because large terms cancel. A small backward error does not guarantee small forward error in an ill-conditioned system.

## Fallback choices and their limits

**Degree fallback is a weighting proxy, not a calibrated uncertainty.** It ignores weak connections, edge-weight magnitudes and reference geometry. It can substantially underestimate uncertainty, especially for an event across a weak bridge or a high-degree hub far from the reference. Keeping checked Hutch values and attempting selected exact corrections first reduces the use of this approximation. If many rows use a proxy, check the component’s reported reason before choosing how to obtain statistical weights.

`pseudo_weight` is another explicit practical-weight proxy based on incident edge weights. It is also not calibrated. `fallback="exact"` repairs all failing entries if the whole component fits within `exact_max_nred`; beyond that cap it uses the normal selected-repair budget and degree fallback for entries that remain unusable. If sampling was skipped, exact fallback is attempted only within the full exact size cap. A failed exact calculation also ends in degree fallback unless a different fallback was explicitly selected.

`graph_bound` remains available **only when explicitly requested** as a mode or fallback. It builds a breadth-first spanning tree and uses

```text
bound_variance(i) = sigma_std^2 * sum(1 / edge_weight along tree path from reference to i)
```

This upper-bounds the pinned inverse-Laplacian diagonal by [Rayleigh monotonicity](https://www.control.utoronto.ca/~jwsimpson/papers/2017k/); a nonnegative ridge only decreases the diagonal. However, discarding alternative paths can make the bound much larger than the actual uncertainty. On a uniformly weighted complete graph with `n` events, its std is approximately `sqrt(n/2)` times exact. It is unsuitable as a general default for practical weighting.

DDSync's exported `thetaStd` weight scales as `1/std` before its cap. A downstream least-squares program using inverse-variance weights scales as `1/std^2`: a 20-fold std inflation means 20-fold and 400-fold weight reductions respectively. Caps and minimum-sigma filters can further change these ratios.

Legacy `fallback="nan"` and `"none"` now use degree fallback. `min_diag_rel` / MATLAB `min_diagrel` remains accepted but ignored. No automatic path selects a graph bound. Explicit `mode="graph_bound"` bypasses sampling; an explicit graph-bound fallback affects only entries still needing fallback after repair and the configured retention policy. To reproduce whole-component bounds, use the explicit mode.

## Configuration

Julia TOML:

```toml
[std]
mode = "auto"
auto_probe = "hutch"  # Set to "cholesky" for shared-shift sampling on larger components.
exact_max_nred = 2048
probes = 4096
min_probes = 128
batch = 128
rel_tol = 0.10
repair_max_nodes = 128
unconverged = "keep"
fallback = "pseudo_degree"
```

MATLAB:

```matlab
cfg = ddsync.config_default();
cfg.std.mode = 'auto';
cfg.std.auto_probe = 'hutch'; % Choose 'cholesky' for squared draws.
cfg.std.exact_max_nred = 2048;
cfg.std.hutch.probes = 4096;
cfg.std.hutch.min_probes = 128;
cfg.std.hutch.batch = 128;
cfg.std.hutch.rel_tol = 0.10;
cfg.std.repair_max_nodes = 128;
cfg.std.unconverged = 'keep';
cfg.std.fallback = 'pseudo_degree';
```

Use `mode="hutch"` or `mode="cholesky"` to request sampling even on a small component, `mode="exact"` for deterministic calculation up to the full exact size cap, or `mode="auto"` for the default selection. Explicit exact mode above the cap uses fallback; increase the cap if that calculation is intended and its cost is acceptable. `repair_max_nodes=0` disables selected repairs. A probe budget below `min_probes` can still be repaired exactly; it cannot produce accepted or retained stochastic values.

Existing configurations retain explicit settings. To adopt the new selection and fallback, change `mode` to `auto` and `fallback` to `pseudo_degree`; new repair/retention fields are supplied by defaults. Explicit `fallback="graph_bound"` remains graph-bound fallback. Existing files with `mode="hutch"` still sample small components. Explicit old probe budgets remain in effect, including `probes=2048`. To adopt the new maximum, change it to 4096 or remove the explicit override. The full exact cap stays at 2048 unpinned events and the selected-repair budget stays at 128.

## Reading per-event quality and component diagnostics

The numeric `theta_*.txt` and `std_theta_*.txt` layouts are unchanged. Each std file now has a companion **`std_quality_STA_PH.txt`**, containing supported component members only:

```text
# event_id ref_id method relative_mc_se
```

IDs follow the corresponding theta file, including Julia's sequential mapping. Rows are grouped by component and need not be globally sorted. Join by event ID, not row position. Reference rows are labeled `reference`. Other method values are `exact`, `hutch`, `hutch_unconverged`, `cholesky`, `cholesky_unconverged`, `pseudo_degree`, `pseudo_weight`, `graph_bound` or deliberately disabled `none`.

`relative_mc_se` is the pre-repair relative sampling diagnostic. Thus an `exact` repair can carry a large value or `Inf` from the failed stochastic attempt; that value does not describe the accuracy of the repaired std. `NaN` means no sampling diagnostic was taken. These diagnostic missing values are **not** uncertainty values passed to GraphSplit.

### Read the percentages first

`sync_metrics.txt` reports the following categories globally and by station-phase group. The terminal prints the final global breakdown too, even if every Hutch estimate converged.

| Category | Included values |
| --- | --- |
| Exact | Full exact calculation and successful selected exact repairs. |
| Converged Hutch | Stochastic entries that passed the configured per-entry diagnostic. |
| Converged Cholesky | Squared-draw estimates that met the per-event sampling target. |
| Kept unconverged Cholesky | Usable squared-draw estimates retained above the target. |
| Kept unconverged Hutch | Positive stochastic entries retained without passing that diagnostic. |
| Fallback/proxy | Degree or weight proxies and explicitly selected graph bounds, including when directly requested. |

Percentages use the **number of supported non-reference std entries** as their denominator. One event observed at several station-phase groups contributes once in each file. They are not percentages of unique catalog events or of components, and are not averages of component percentages. Gauge references, rejected/isolated rows and deliberately disabled uncertainty are excluded. Disabled entries are counted separately. With no assessed entries, percentages are `N/A`, not a claim of zero uncertainty problems.

The **requires-review share** is the sum of kept unconverged estimates from both samplers and fallback/proxy percentages. It shows which rows lack an exact or converged stochastic diagonal; it is not a measured fraction of physically incorrect uncertainties. Conversely, exact plus converged estimates from both samplers is the numerically checked share, not proof of calibrated physical uncertainty.

At **10% or more requiring review**, DDSync flags the run or the individual station-phase group. Check which category accounts for the flagged entries. If many Hutch estimates miss the sampling target, try Cholesky before increasing the number of probes. For groups that skipped sampling or used a proxy, check the component reason and size limit first. The 10% review threshold counts entries. The separate `rel_tol=0.10` setting tests the sampling precision of each estimated variance.

The relative-MC-SE cutoff at **1** is still active: an unrepaired entry with SE as large as or larger than its estimate is assigned fallback instead of being kept. Its replacement is now visible in the fallback count and review percentage. The 4096-probe default reduces sampling noise but does not guarantee that this cutoff or the convergence tolerance will be met.

### Component details at the end

The observation table uses aligned whitespace-separated columns. Component details are grouped separately at the end of the report and include successful convergence as well as failures. Each record starts with `# std_component`; its fields continue on following comment lines. Parsers must no longer assume all fields occur on the opening line. Each record includes:

| Field | Meaning |
| --- | --- |
| `method` | Actual method, or `hybrid` when different events use different methods. |
| `reason` | Sampling outcome (`converged`, `not_converged`, `nonpositive_estimate`, `nonfinite_estimate`, `insufficient_probes`, `solve_failed`), a skip reason (`disabled_probes`, `size_limit`, `exact_size_limit`), `exact_failed`, or `requested`. It describes the attempt, before repairs. |
| `probes`, `max_rel_se` | Actual sampling count and worst pre-repair relative MC standard error. |
| `sigma_hat`, `sigma_std` | Measured scatter of fitted-minus-observed differential travel times, and the scatter used in uncertainty formulas after applying the optional minimum, both in seconds. |
| `cholesky_converged` | Events whose Cholesky variance estimate met the sampling target. |
| `cholesky_unconverged` | Events whose Cholesky estimate was kept above the target. |
| `hutch_converged` | Events whose Hutchinson variance estimate met the sampling target. |
| `exact_nodes` | Number supplied by full or selected exact calculation. |
| `repair_attempted`, `repaired_nodes` | Selected solves attempted and successfully used. |
| `unconverged_nodes` | Retained positive stochastic estimates that failed the requested tolerance. |
| `fallback_nodes` | Entries supplied by a degree/weight proxy or an explicitly selected graph bound. |

`repaired_nodes` is a subset of `exact_nodes`, not an additional category. The reference is excluded from these counts. `hutch_converged + cholesky_converged + exact_nodes + unconverged_nodes + fallback_nodes` equals the number of unpinned events when uncertainty is enabled. `pseudo_components_used` counts components using at least one proxy entry. The group-table `sigma_hat` is only a reporting median across final components; each component uses its own scale.

## Choosing budgets and using the exports

Start with `mode = "auto"` and the 4096-probe maximum. Read the global and group review percentages, then inspect `unconverged_nodes`, `fallback_nodes`, component reasons and the quality files. If many Hutch estimates miss the target, set `auto_probe = "cholesky"`. For explicit Hutch mode, change `mode = "hutch"` to `mode = "cholesky"`. Keep the probe limit and `rel_tol` fixed to compare both samplers on the same group. If only a few entries need attention, increasing `repair_max_nodes` lets DDSync calculate those variances exactly.

A larger exact calculation is also an option when its cost is affordable. For an exploratory Hutch run, 128 probes can suffice on favorable graphs. Keep the precision target appropriate for the analysis when interpreting the convergence counts.

For final analyses, compare representative components across independent seeds or against exact diagonals. The [November synthetic study](updates/2026-11-12-per-event-uncertainty.md) includes 2,400–5,000-event components, both language implementations, error distributions and failure modes. That study used 2048 probes, the default at the time; its measurements have not been relabeled as 4096-probe results. It showed that retained positive estimates on a large dense graph can have sizable errors. See the [19 November reporting update](updates/2026-11-19-uncertainty-reporting.md) for the current summary layout and default.

For reproducibility, use Julia's `Random.seed!(1234)` before `DDSync.run(cfg)`, or MATLAB's `rng(1234)` before `ddsync.run(cfg)`. The two languages do not share random streams. Exact calculations provide the basis for cross-language equality checks.

Estimator failure does not deliberately export missing std for accepted theta rows. Defaults preserve finite uncertainties and reference constraints for GraphSplit. Invalid noise scales, underflow/overflow or otherwise unusable fallback values stop the run explicitly. Explicit `mode="none"` or `export=false` remains a theta-only workflow; avoid those choices when downstream uncertainty is required. Rejected/isolated events still have missing theta/reference and must be skipped. A downstream maximum-sigma filter can still remove constraints.

`robust.min_scale` stabilizes pruning/Huber weights. `std.min_sigma` sets the minimum assumed scatter in the travel-time mismatches when calculating uncertainty, not each row's final std. Disabling this floor permits zero std for exactly consistent synthetic data. Output precision controls decimal formatting only; temporary data and pass-2 storage use single precision. Additional printed decimals do not establish equivalent physical accuracy.

Rerun affected datasets to benefit from per-event handling. Old std files cannot reveal which Hutch estimates had converged before a whole-component fallback replaced them.
