# DDSync (MATLAB)

Julia is the primary implementation. MATLAB receives bug fixes. See the [19 November 2026 reporting update](../docs/updates/2026-11-19-uncertainty-reporting.md) and the [uncertainty and probe guide](../docs/uncertainty.md).

**Event IDs in this implementation:** MATLAB does not perform the automatic
mapping available in DDSyncJulia. Reindex the catalog and differential-time
pair headers consistently to `1:N` before running; IDs must not exceed
16,777,216 because the temporary group files store them in single precision.
The helpers in `extras/` can be adapted for this preparation. Keep the ID
correspondence, and use the same reindexed catalog with the resulting
`theta/` and `thetastd/` files in GraphSplit. Reordering catalog rows is safe;
renumbering only the catalog after DDSync is not. The automatic indexing feature remains Julia-only; the pruning and uncertainty bug fixes apply to both implementations.

DDSync reads a HypoDD/GrowClust-style `dt.cc` file, processes one station-phase group at a time, and writes a cleaned and synchronized differential-time file. For each station and phase, it estimates one relative arrival-time value, `theta`, for every event that can be constrained by the data. A kept observation between events `i` and `j` is then represented by `theta(i) - theta(j)`.

The MATLAB and Julia versions use the same synchronization defaults and numerical formats. Automatic handling of external IDs is specific to Julia. This README describes MATLAB configuration and its required input preparation.

---

## Quick start

From inside the `DDSync/` folder:

```matlab
% If you already edited ddsync_params.m, keep your file.
% Otherwise copy the template (this replaces the bundled parameter file):
copyfile('ddsync_params_template.m', 'ddsync_params.m')
edit ddsync_params.m
run_ddsync
```

Or build a configuration directly:

```matlab
cfg = ddsync.config_default();
cfg.io.infile_dt    = 'dt.cc';
cfg.io.catalog_file = 'catalog.txt';
out = ddsync.run(cfg);
```

Relative input and output paths are resolved from MATLAB's working directory. The bundled `dt.cc` and `catalog.txt` can be used immediately from this folder. Use a separate output directory for each run: matching output files are overwritten, and old group files are not automatically removed.

The default `auto` mode uses exact calculation directly for at most 2048 unpinned events. Larger components use up to 4096 probes, retain converged entries, and repair up to 128 failing entries exactly. Remaining estimates are individually retained with a quality label or given degree fallback. For fast exploratory weighting you can explicitly choose a degree proxy; it is not a statistical uncertainty:

```matlab
cfg.std.mode = 'pseudo_degree';
cfg.std.fallback = 'pseudo_degree';
```

---

## Requirements

- MATLAB R2021a or later; native regression tests were run with R2025b.
- No extra toolboxes are required.
- Use compact positive event IDs `1:N` in the catalog and pair headers; the maximum is 16,777,216.

---

## Directory layout

- `+ddsync/` - maintained MATLAB package code.
- `run_ddsync.m` - wrapper that loads `ddsync_params.m` if it exists, otherwise uses defaults.
- `ddsync_params_template.m` - user-editable configuration template.
- `run_ddsync_example.m` - small example showing common overrides.
- `dt.cc`, `catalog.txt` - bundled Spanish Springs example data.
- `theta/`, `thetastd/` - output directories created at runtime.
- `extras/` - dataset-specific preprocessing examples and a legacy entry point. `synchro_polished.m` now calls the maintained package; migrate old uppercase knobs to `cfg` fields as shown in that file. Read `extras/README.md` before using the preprocessing scripts.

---

## Input files

### `dt.cc`

DDSync expects block-structured differential-time data:

```text
# i j 0.0
STA  dt  cc  PH
STA  dt  cc  PH
# i j 0.0
...
```

A header line starts with `#` and gives the event pair `i j`. The third header value is ignored. Each following station line supplies one measurement for that event pair.

| Field | Meaning |
| --- | --- |
| `STA` | Station code. It is treated as text. |
| `dt` | Observed differential travel time in seconds for the current event pair. |
| `cc` | Correlation coefficient or quality value. DDSync converts this to a base weight. |
| `PH` | Phase label, commonly `P` or `S`. |

DDSync groups measurements by `(STA, PH)`. For example, all `DYN P` rows are solved as one graph and all `DYN S` rows are solved as another graph.

Event IDs must be suitable for array-style output: DDSync writes rows `1:maxEventID` in each `theta` and `std_theta` file. If your catalog uses arbitrary IDs, reindex before running. The `extras/reindext_dtcc.m` and `extras/reorder.m` scripts are examples of this kind of preprocessing.

### `catalog.txt`

DDSync uses `catalog.txt` only to find the largest event ID. The event ID must be the last whitespace-separated column on each non-empty line. If your catalog uses a different layout, edit the `readMaxEventID` helper in `+ddsync/run.m` accordingly.

---

## Output files

### `dt_sync.cc`

This is a HypoDD-style differential-time file with synchronized `dt` values. For each kept measurement between events `i` and `j`, DDSync writes:

```text
STA  theta(i)-theta(j)  weight  PH
```

The third column is controlled by `cfg.output.dt_weight_mode`.

Pruned edges are omitted by default. If `cfg.output.write_pruned_edges = true`, DDSync keeps the original station line but writes its weight as `0`.

### `theta/theta_<STA>_<PH>.txt`

Each file has three columns:

```text
EventID  theta  refEventID
```

`theta` is a relative arrival-time value in seconds. After pruning, DDSync rebuilds the connected components and selects the event with the highest surviving degree in each component as its reference (ties use the smallest event ID). It refits each component independently and writes that reference in `refEventID`. The reference has `theta = 0`.

Both the minimum-edge gate and the requirement for at least one cycle are checked again after pruning. Isolated events and rejected fragments have `NaN` theta, reference and uncertainty, with degree zero. Their observations are omitted, or written with zero weight when diagnostic pruned output is enabled.

**Compare or subtract theta values only when they belong to the same station, phase and finite `refEventID`.** Different references identify independently floating components. A theta star must connect only an event and its own reference. For a common anchor across stations/channels, check that the target and anchor share a reference in each file before subtracting. Theta is not an absolute pick time.

### `thetastd/std_theta_<STA>_<PH>.txt`

The base file has four columns:

```text
EventID  std_theta  refEventID  degree
```

Optional columns may be added:

```text
EventID  std_theta  refEventID  degree  thetaWeight  nodeWeight
```

`std_theta` is a converged or explicitly labeled unconverged Hutchinson estimate, an exact inverse-diagonal calculation, an explicitly requested graph bound, a fast pseudo standard deviation, or explicitly disabled (`NaN`) uncertainty, depending on `cfg.std.mode` and fallback behavior. `degree` is the number of kept incident edges in the final graph. `thetaWeight` is a practical inverse-uncertainty weight derived from `std_theta`. `nodeWeight` is an optional diagnostic node-weight summary used by `pseudo_weight` mode.

### `sync_metrics.txt`

A readable report with aligned, whitespace-separated station-phase columns. The report gives uncertainty counts and percentages, observation statistics, and group-level uncertainty summaries. Component details are collected separately at the end; they no longer interrupt the observation rows. The terminal also prints the final uncertainty summary, including successful Hutchinson estimates.

```text
group    n_edges  n_kept  pruned_pct  robust_p50  robust_p05  robust_p95  frac_robust_lt_0p99  sigma_hat
```

`group` is the station-phase key, `n_edges` is the number of input measurements, and `n_kept` is the number retained after component gating and pruning. The four uncertainty categories are **exact**, **converged Hutch**, **kept unconverged Hutch**, and **fallback/proxy**. Exact includes selected exact repairs; fallback includes explicitly selected degree/weight proxies or graph bounds. Counts and percentages describe supported non-reference std entries across station-phase files, not unique catalog events. Gauge references, rejected/isolated rows and deliberately disabled uncertainty are excluded; disabled entries are counted separately. Percentages are `N/A` if there are no assessed entries.

The **requires-review share** is kept unconverged Hutch plus fallback. At 10% or more, the report advises increasing probes and/or exact repairs, or considering a larger full exact cap. Groups receive their own advice so a poorly estimated small group cannot disappear inside a good global average. This reporting threshold is an operational prompt, not a statistical confidence bound. A large exact-plus-converged share is favorable; it does not itself trigger advice to increase sampling.

All scale values are in seconds. After the final fit, DDSync compares each retained measured differential travel time with `theta_i - theta_j`. Component `sigma_hat` measures the scatter of these mismatches within one connected group at one station and phase; DDSync uses it to estimate noise in the pairwise measurements.

`sigma_std` is the scale actually used in uncertainty formulas after the optional `min_sigma` floor: normally `max(sigma_hat, min_sigma)` when the floor is enabled, otherwise `sigma_hat`. The `sigma` in the degree-proxy formula means **`sigma_std`**, not a separate noise estimate. Neither scale is a particular event's std or the Monte Carlo sampling error.

The group-table `sigma_hat` is the median of raw final component scales; each component uses its own scale. `pseudo_weight` does not use either scale.

Component records begin with `# std_component` and continue on subsequent comment lines. They include successful methods as well as failures, both scale values, probes, the pre-repair diagnostic and method counts. `method=hybrid` means events used different methods. Use `std_quality_*.txt` to identify individual entries. Parsers should split the observation table on whitespace, and should not assume all component fields are on one line. See the [diagnostic field guide](../docs/uncertainty.md#reading-per-event-quality-and-component-diagnostics).


---

## Precision for small times

The old fixed `%.6f` formatting could round away information when working at microsecond or sub-microsecond scale. DDSync now exposes separate precision controls:

```matlab
cfg.output.theta_decimals    = 9;  % theta files, seconds
cfg.output.thetastd_decimals = 9;  % std_theta files, seconds
```

These are digits after the decimal point. With seconds as the unit, 6 decimals is 1 microsecond, 9 decimals is 1 nanosecond, and 12 decimals is 1 picosecond. Increase these values for field-lab, acoustic, or other small-time experiments where very small differential times matter.

These settings affect only the `theta` and `std_theta` value columns. The synchronized `dt` column in `dt_sync.cc` is controlled by `cfg.output.dt_decimals`, and the weight column in `dt_sync.cc` is controlled by `cfg.output.dt_weight_decimals`.

Set both values to `6` if you need the older text-file formatting. More printed digits do not increase numerical accuracy: temporary observations and pass-2 theta/std storage currently use Float32/single precision. A value printed to nanoseconds is not evidence of nanosecond accuracy.

---

## Residual scale floors: `min_scale` versus `min_sigma`

For each retained pair, the fitted theta values predict a differential travel time `theta_i - theta_j`. The difference between this prediction and the measured differential travel time is the **residual**. DDSync measures the scatter of these residuals separately for each connected group at one station and phase. It uses that scatter as an estimate of noise in the pairwise measurements.

The two settings below impose a minimum on the scatter used for different purposes.

### `cfg.robust.min_scale`

This protects pruning and IRLS. During outlier detection and Huber reweighting, DDSync divides residuals by a robust scale estimate. If the residuals are exactly zero or nearly zero, the scale can become zero and the division becomes unstable. `cfg.robust.min_scale` prevents that. It can affect which edges are pruned and how robust weights are computed.

### `cfg.std.min_sigma` and `cfg.std.apply_min_sigma`

`min_sigma = 5e-4` means: **when calculating theta uncertainties, assume at least 0.5 ms scatter in the mismatches between measured and fitted differential travel times.** DDSync combines that assumed scatter with the number, weights and arrangement of observations (or just the degree for `pseudo_degree`) to calculate `std_theta`. It applies the same minimum for exact, Hutchinson, graph-bound and degree-proxy calculations.

**This is not a minimum theta uncertainty.** Combining observations can give exported `std_theta` values smaller than `min_sigma`.

```matlab
cfg.std.min_sigma = 5e-4;
cfg.std.apply_min_sigma = true;
```

When `cfg.std.apply_min_sigma = true`, DDSync uses at least `cfg.std.min_sigma` for this mismatch scatter before computing exported `std_theta`. When `cfg.std.apply_min_sigma = false` or `cfg.std.min_sigma = []` or `0`, this minimum assumed mismatch scatter is disabled.

Example: suppose the measured scatter of the final travel-time mismatches is `1e-5` seconds and `cfg.std.min_sigma = 5e-4`.

- With `cfg.std.apply_min_sigma = true`, the mismatch scatter used in the formula is `5e-4` seconds. With `pseudo_degree` and degree 100, the resulting `std_theta` is `5e-5` seconds, below `min_sigma`.
- With `cfg.std.apply_min_sigma = false`, the mismatch scatter used in the formula is `1e-5` seconds. With degree 100, the resulting degree-proxy `std_theta` is `1e-6` seconds.

`min_sigma` does not floor every row in `std_theta_*.txt`. It sets a minimum assumed scatter in the travel-time mismatches before converting that scatter into theta uncertainties. The final row values can be smaller or larger, and the pinned reference event remains `0`.

`pseudo_weight` mode does not use `sigma_hat`, so `min_sigma` does not affect `pseudo_weight` standard deviations.

---

## Configuration reference

All defaults come from `ddsync.config_default()`. You usually override them in `ddsync_params.m`.

### I/O fields

| Field | Default | Plain-language meaning |
| --- | ---: | --- |
| `cfg.io.infile_dt` | `'dt.cc'` | Input differential-time file. |
| `cfg.io.catalog_file` | `'catalog.txt'` | Catalog file used to find `maxEventID`; event ID must be in the last column. |
| `cfg.io.out_dt_sync` | `'dt_sync.cc'` | Output synchronized differential-time file. |
| `cfg.io.metrics_file` | `'sync_metrics.txt'` | Output metrics summary. |
| `cfg.io.tmpdir` | `'ddsync_tmp'` | Scratch directory for per-group spool files and edge-decision files. Put this on fast local storage for large jobs. |
| `cfg.io.thetadir` | `'theta'` | Directory for `theta_<STA>_<PH>.txt`. |
| `cfg.io.thetastd_dir` | `'thetastd'` | Directory for `std_theta_<STA>_<PH>.txt`. |
| `cfg.io.prog_every_lines` | `2e6` | Print progress every this many input lines while streaming. |
| `cfg.io.print_new_groups` | `false` | Print a message whenever a new station-phase group is first seen. Useful for debugging input parsing. |

### Base-weight fields

| Field | Default | Plain-language meaning |
| --- | ---: | --- |
| `cfg.weights.base_fun` | `'cc'` | How to turn the input `cc` value into the starting edge weight. Use `'cc'`, `'cc2'`, or `'ones'`. |
| `cfg.weights.base_fun_handle` | `[]` | Optional MATLAB function handle, for example `@(cc) max(cc,0).^2`. If non-empty, it overrides `base_fun`. |

Examples:

```matlab
cfg.weights.base_fun = 'cc2';       % stronger preference for high-correlation pairs
cfg.weights.base_fun_handle = @(cc) double(cc >= 0.7);
```

DDSync never keeps edges whose base weight is non-positive.

### Robust pruning and IRLS fields

| Field | Default | Plain-language meaning |
| --- | ---: | --- |
| `cfg.robust.k_sigma` | `20` | Prune an edge when its residual is larger than this many robust scale units. Smaller values prune more aggressively. |
| `cfg.robust.min_edges` | `30` | Minimum number of positive-weight edges required both before and after pruning for each component. DDSync also requires at least one loop, meaning more measurements than a bare tree. |
| `cfg.robust.irls_iters` | `10` | Maximum Huber IRLS iterations after initial pruning. Set to `0` to disable robust reweighting. |
| `cfg.robust.huber_c` | `1.345` | Huber tuning constant. Smaller values downweight moderate residuals more strongly. |
| `cfg.robust.irls_rel_tol` | `1e-3` | Stop IRLS early when the objective changes by less than this relative amount. |
| `cfg.robust.min_scale` | `5e-4` | Scale floor used only inside pruning and IRLS to avoid divide-by-zero behavior. |

A practical interpretation: `min_edges` prevents tiny isolated families from producing overconfident `theta` values. `k_sigma` controls how far a measurement can disagree with the internally consistent solution before it is removed.

### Numerical field

| Field | Default | Plain-language meaning |
| --- | ---: | --- |
| `cfg.numeric.ridge_eps` | `1e-10` | Small diagonal ridge added to the reduced Laplacian solve for numerical stability. Increase only if solves are unstable; too large a value can slightly bias `theta`. |

### Output and formatting fields

| Field | Default | Plain-language meaning |
| --- | ---: | --- |
| `cfg.output.write_dt_sync` | `true` | Write `dt_sync.cc`. Set false if you only need `theta`, `std_theta`, and metrics. |
| `cfg.output.dt_weight_mode` | `'thetaStd'` | Weight written in the third column of `dt_sync.cc`: `'base'`, `'robust'`, `'combined'`, or `'thetaStd'`. |
| `cfg.output.dt_decimals` | `5` | Digits after the decimal for synchronized `dt` values in `dt_sync.cc`. |
| `cfg.output.dt_weight_decimals` | `4` | Digits after the decimal for the `dt_sync.cc` weight column. |
| `cfg.output.theta_decimals` | `9` | Digits after the decimal for `theta` values in `theta_*.txt`. |
| `cfg.output.thetastd_decimals` | `9` | Digits after the decimal for `std_theta` values in `std_theta_*.txt`. |
| `cfg.output.station_field_width` | `8` | Width used to align the station field in `dt_sync.cc`. |
| `cfg.output.dt_field_width` | `10` | Width used to align the `dt` field in `dt_sync.cc`. |
| `cfg.output.weight_field_width` | `8` | Width used to align the weight field in `dt_sync.cc`. |
| `cfg.output.write_pruned_edges` | `false` | If false, pruned station lines are omitted. If true, they are written with weight `0`. |
| `cfg.output.theta_weight_scale_mode` | `'fixed'` | For `dt_weight_mode='thetaStd'`: use `'fixed'` scale or compute a global `'median'` scale from kept observations. |
| `cfg.output.theta_weight_scale_fixed` | `500` | Fixed divisor for thetaStd weights. Example: if `std_dt=0.002 s`, then `(1/std_dt)/500 = 1`. |
| `cfg.output.theta_weight_cap` | `1.0` | Maximum written thetaStd-derived weight. Use `Inf` to disable the cap. |

`dt_weight_mode` options:

| Mode | Written weight |
| --- | --- |
| `'base'` | Base weight from `cc`. |
| `'robust'` | Huber robust weight only. |
| `'combined'` | Base weight times robust weight. |
| `'thetaStd'` | Inverse uncertainty proxy for `theta(i)-theta(j)`, scaled and capped; unavailable uncertainty gets weight zero. |

For `thetaStd`, DDSync uses `sqrt(std_theta_i^2 + std_theta_j^2)` and ignores covariance. With valid inverse-Laplacian uncertainties this overestimates pair uncertainty when covariance is positive. The pseudo modes do not provide that guarantee. Unavailable uncertainty produces weight zero and is excluded from median scaling. If you disable std export, choose `base`, `robust`, or `combined` output weights instead of `thetaStd`.

### Standard-deviation export fields

| Field | Default | Plain-language meaning |
| --- | ---: | --- |
| `cfg.std.export` | `true` | Write `std_theta_*.txt` and per-event `std_quality_*.txt` in the std directory. |
| `cfg.std.mode` | `'auto'` | Exact up to the size cap, otherwise per-event Hutchinson with selected exact repair. Explicit `hutch`, `exact`, `graph_bound`, `pseudo_degree`, `pseudo_weight`, and `none` remain available. |
| `cfg.std.fallback` | `'pseudo_degree'` | Only entries still unusable after repair/retention use fallback. Options: `pseudo_degree`, `pseudo_weight`, `exact`, `graph_bound`. Legacy `nan`/`none` aliases now use degree fallback. |
| `cfg.std.exact_max_nred` | `2048` | Full exact size cap in unpinned events. `auto` switches to Hutch above it; explicit exact mode above it uses fallback. |
| `cfg.std.repair_max_nodes` | `128` | Maximum selected full-system exact corrections after Hutch, reusing its factorization. Applies independently of component size. |
| `cfg.std.unconverged` | `'keep'` | Retain unrepaired positive estimates with enough samples and relative MC SE below 1, labeled `hutch_unconverged`. Use `fallback` to replace all failing entries instead. |
| `cfg.std.min_sigma` | `5e-4` | Minimum assumed scatter of measured differential travel times around fitted `theta_i - theta_j`, in seconds, when calculating uncertainty. **Not a minimum `std_theta`; exported values can be smaller.** |
| `cfg.std.apply_min_sigma` | `true` | Switch for `min_sigma`. Set false, or set `min_sigma` to `[]` or `0`, to disable this minimum assumed mismatch scatter. These options do not floor final theta uncertainties. |
| `cfg.std.pseudo_weight_source` | `'combined'` | Edge weight source for `pseudo_weight` and `nodeWeight`: `'combined'`, `'robust'`, or `'base'`. |
| `cfg.std.pseudo_weight_eps` | `1e-6` | Minimum node weight before inverting in `pseudo_weight` mode. Prevents division by zero. |
| `cfg.std.write_weight_column` | `true` | Add the practical `thetaWeight` column to `std_theta_*.txt`. |
| `cfg.std.write_alt_weight_column` | `false` | Add `nodeWeight` to `std_theta_*.txt`. If `write_weight_column` is also true, this is the sixth column; otherwise it is the fifth. |

Mode meanings:

| Mode | Meaning |
| --- | --- |
| `'auto'` | Exact directly for small components, otherwise Hutchinson and per-event repair. |
| `'hutch'` | Sample even small components. Preserve converged entries, repair selected failures, and retain or replace only remaining failing entries. |
| `'exact'` | Unit-vector solves up to the full exact cap; no sampling noise. Beyond that limit use configured fallback (degree by default). |
| `'graph_bound'` | Explicit, potentially very loose upper bound from tree path resistances. Never selected automatically by the defaults. |
| `'pseudo_degree'` | Fast proxy: `std_theta ~= sigma_std/sqrt(degree)`. A weighting proxy, not a calibrated or conservative uncertainty. |
| `'pseudo_weight'` | Fast proxy based on incident edge weights, scaled by `theta_weight_scale_fixed`. Ignores `min_sigma`. |
| `'none'` | Do not compute finite standard deviations. If `std.export` is true, files are still written with `NaN` for non-reference events and `0` for each reference. |

### Hutchinson fields

| Field | Default | Plain-language meaning |
| --- | ---: | --- |
| `cfg.std.hutch.probes` | `4096` | Maximum probe budget; stops early if all nodes pass the diagnostic. Set `0` to use fallback directly. |
| `cfg.std.hutch.min_probes` | `128` | Minimum samples before accepting convergence (at least 2). A smaller budget permits exact repair or fallback, but no accepted/retained stochastic values. |
| `cfg.std.hutch.rel_tol` | `0.10` | Maximum empirical standard error divided by estimated diagonal, tested independently for each unpinned node. Not a guaranteed relative error bound. |
| `cfg.std.hutch.probe_dist` | `'rademacher'` | Random probe distribution: `'rademacher'` (`+1/-1`) or `'gaussian'`. |
| `cfg.std.hutch.batch` | `128` | Number of probes solved at once. Larger can be faster but uses more memory. |
| `cfg.std.hutch.report_every_batch` | `10` | Print progress every this many batches. Set `0` to silence batch progress. |
| `cfg.std.hutch.max_nred` | `40000` | Skip Hutchinson if the reduced system is larger than this many unknowns; use fallback instead. |
| `cfg.std.hutch.min_diagrel` | `1e-12` | Deprecated and ignored. Invalid variance estimates trigger fallback; they are never clipped to a positive floor. |

Each std file has a `std_quality_STA_PH.txt` companion listing supported event ID, reference ID, method and pre-repair relative MC standard error. Join by event ID; rows need not be globally sorted. Global and group metrics give counts and percentages for converged, exact, retained unconverged and fallback entries; detailed component records follow at the end. Existing numeric theta/std columns are unchanged. See the [November update](../docs/updates/2026-11-12-per-event-uncertainty.md).

Hutchinson sampling error usually decreases like `1/sqrt(probes)`, with strong graph dependence. See the [shared uncertainty guide](../docs/uncertainty.md) for choosing a budget, checking convergence, and interpreting fallback.

---

## Common configuration examples

### Preserve microsecond or sub-microsecond values

```matlab
cfg.output.theta_decimals = 12;
cfg.output.thetastd_decimals = 12;
cfg.output.dt_decimals = 9;  % also increase dt_sync precision if needed
```

### Use the measured mismatch scatter without a minimum

```matlab
cfg.std.apply_min_sigma = false;
% or
cfg.std.min_sigma = [];
```

This does not change `cfg.robust.min_scale`, so pruning and IRLS remain protected from zero residual scales.

### Use fast practical standard deviations

```matlab
cfg.std.mode = 'pseudo_degree';
cfg.std.fallback = 'pseudo_degree';
```

### Keep pruned lines for diagnostics

```matlab
cfg.output.write_pruned_edges = true;
```

Pruned station lines are then written with weight `0` instead of being omitted.

---

## Notes on MATLAB/Julia parity

The MATLAB and Julia implementations are intended to produce the same output formats. Small floating-point differences can occur because MATLAB and Julia use different sparse linear algebra libraries and random-number streams, especially for Hutchinson standard deviations. For deterministic comparisons, use the same precision settings and `cfg.std.mode = 'exact'` on components within the exact size limit.


## Tests

From this folder, run `addpath('test'); run_tests` in MATLAB. The suite uses synthetic inputs shared with Julia and checks references after pruning, isolated/tree fragments, exact covariance, stochastic convergence, fallback and exported weights. No downloaded data are required. With `test/` on the path, run `large_uncertainty_study('large.csv')` for the new large-component comparison or `probe_study('probe_study.csv')` for the earlier small-graph study.
