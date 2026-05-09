# DDSync (MATLAB)

DDSync reads a HypoDD/GrowClust-style `dt.cc` file, processes one station-phase group at a time, and writes a cleaned and synchronized differential-time file. For each station and phase, it estimates one relative arrival-time value, `theta`, for every event that can be constrained by the data. A kept observation between events `i` and `j` is then represented by `theta(i) - theta(j)`.

The MATLAB and Julia versions are intended to use the same defaults and produce the same file formats. This README uses MATLAB `cfg` names. The Julia README contains the same guidance using Julia/TOML names.

---

## Quick start

From inside the `DDSync/` folder:

```matlab
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

For a fast first run on a large file, disable Hutchinson standard deviations and use the degree proxy:

```matlab
cfg.std.mode = 'pseudo_degree';
cfg.std.fallback = 'pseudo_degree';
```

---

## Requirements

- MATLAB.
- No extra toolboxes are required.
- Event IDs in the input should be positive integers that fit inside MATLAB array indexing.

---

## Directory layout

- `+ddsync/` - maintained MATLAB package code.
- `run_ddsync.m` - wrapper that loads `ddsync_params.m` if it exists, otherwise uses defaults.
- `ddsync_params_template.m` - user-editable configuration template.
- `run_ddsync_example.m` - small example showing common overrides.
- `dt.cc`, `catalog.txt` - bundled Spanish Springs example data.
- `theta/`, `thetastd/` - output directories created at runtime.
- `extras/` - helper and legacy scripts. The parity target is the maintained `+ddsync/` package code.

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

`theta` is the relative arrival-time value for that event in seconds. `refEventID` is the pinned reference event used to fix the arbitrary additive offset in that connected component. The pinned event has `theta = 0`. Events that cannot be constrained are written as `NaN`.

### `thetastd/std_theta_<STA>_<PH>.txt`

The base file has four columns:

```text
EventID  std_theta  refEventID  degree
```

Optional columns may be added:

```text
EventID  std_theta  refEventID  degree  thetaWeight  nodeWeight
```

`std_theta` is either a Hutchinson estimate, a fast pseudo standard deviation, or `NaN`, depending on `cfg.std.mode` and fallback behavior. `degree` is the number of kept incident edges in the final graph. `thetaWeight` is a practical inverse-uncertainty weight derived from `std_theta`. `nodeWeight` is an optional diagnostic node-weight summary used by `pseudo_weight` mode.

### `sync_metrics.txt`

A tab-delimited summary. Each station-phase row has:

```text
group    n_edges  n_kept  pruned_pct  robust_p50  robust_p05  robust_p95  frac_robust_lt_0p99  sigma_hat
```

`group` is the station-phase key, `n_edges` is the number of input measurements, and `n_kept` is the number retained after component gating and pruning. `sigma_hat` is the final robust residual scale for reporting and for `std_theta` export. It is not the same thing as `cfg.std.min_sigma`.

---

## Precision for small times

The old fixed `%.6f` formatting could round away information when working at microsecond or sub-microsecond scale. DDSync now exposes separate precision controls:

```matlab
cfg.output.theta_decimals    = 9;  % theta files, seconds
cfg.output.thetastd_decimals = 9;  % std_theta files, seconds
```

These are digits after the decimal point. With seconds as the unit, 6 decimals is 1 microsecond, 9 decimals is 1 nanosecond, and 12 decimals is 1 picosecond. Increase these values for field-lab, acoustic, or other small-time experiments where very small differential times matter.

These settings affect only the `theta` and `std_theta` value columns. The synchronized `dt` column in `dt_sync.cc` is controlled by `cfg.output.dt_decimals`, and the weight column in `dt_sync.cc` is controlled by `cfg.output.dt_weight_decimals`.

Set both values to `6` if you need the older text-file formatting.

---

## Residual scale floors: `min_scale` versus `min_sigma`

DDSync has two separate safeguards. They intentionally do different jobs.

### `cfg.robust.min_scale`

This protects pruning and IRLS. During outlier detection and Huber reweighting, DDSync divides residuals by a robust scale estimate. If the residuals are exactly zero or nearly zero, the scale can become zero and the division becomes unstable. `cfg.robust.min_scale` prevents that. It can affect which edges are pruned and how robust weights are computed.

### `cfg.std.min_sigma` and `cfg.std.apply_min_sigma`

This applies only when converting the final residual-noise estimate into exported `std_theta` values for Hutchinson and `pseudo_degree` modes.

```matlab
cfg.std.min_sigma = 5e-4;
cfg.std.apply_min_sigma = true;
```

When `cfg.std.apply_min_sigma = true`, DDSync uses at least `cfg.std.min_sigma` as the residual-noise scale before computing exported `std_theta`. When `cfg.std.apply_min_sigma = false` or `cfg.std.min_sigma = []` or `0`, this export floor is disabled.

Example: suppose the final residual scale is `1e-5` seconds and `cfg.std.min_sigma = 5e-4`.

- With `cfg.std.apply_min_sigma = true`, `std_theta` export uses `5e-4` seconds before graph leverage or degree scaling.
- With `cfg.std.apply_min_sigma = false`, `std_theta` export uses `1e-5` seconds.

`min_sigma` does not floor every row in `std_theta_*.txt`. It floors the residual-noise scale before the graph/degree conversion. The final row values can be smaller or larger, and the pinned reference event remains `0`.

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
| `cfg.robust.min_edges` | `30` | Minimum number of positive-weight edges required before a component is trusted. DDSync also requires at least one loop, meaning more measurements than a bare tree. |
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
| `'thetaStd'` | Conservative inverse uncertainty for `theta(i)-theta(j)`, scaled and capped. |

For `thetaStd`, DDSync uses `sqrt(std_theta_i^2 + std_theta_j^2)` and ignores covariance. This is conservative for many relocation-weighting uses.

### Standard-deviation export fields

| Field | Default | Plain-language meaning |
| --- | ---: | --- |
| `cfg.std.export` | `true` | Create the `thetastd/` directory and write `std_theta_*.txt`. |
| `cfg.std.mode` | `'hutch'` | Main method for `std_theta`: `'hutch'`, `'pseudo_degree'`, `'pseudo_weight'`, or `'none'`. |
| `cfg.std.fallback` | `'pseudo_degree'` | Method used when Hutchinson is skipped. Use `'pseudo_degree'`, `'pseudo_weight'`, or `'nan'`. |
| `cfg.std.min_sigma` | `5e-4` | Optional residual-noise floor used only for exported Hutchinson and `pseudo_degree` `std_theta`. |
| `cfg.std.apply_min_sigma` | `true` | Switch for `min_sigma`. Set false, or set `min_sigma` to `[]` or `0`, to disable the export floor. |
| `cfg.std.pseudo_weight_source` | `'combined'` | Edge weight source for `pseudo_weight` and `nodeWeight`: `'combined'`, `'robust'`, or `'base'`. |
| `cfg.std.pseudo_weight_eps` | `1e-6` | Minimum node weight before inverting in `pseudo_weight` mode. Prevents division by zero. |
| `cfg.std.write_weight_column` | `true` | Add the practical `thetaWeight` column to `std_theta_*.txt`. |
| `cfg.std.write_alt_weight_column` | `false` | Add `nodeWeight` to `std_theta_*.txt`. If `write_weight_column` is also true, this is the sixth column; otherwise it is the fifth. |

Mode meanings:

| Mode | Meaning |
| --- | --- |
| `'hutch'` | Estimate the diagonal of the inverse graph Laplacian using random probes. Most statistical, but potentially expensive. |
| `'pseudo_degree'` | Fast proxy: `std_theta ~= sigma/sqrt(degree)`. Good for large jobs where relative weighting is enough. |
| `'pseudo_weight'` | Fast proxy based on incident edge weights, scaled by `theta_weight_scale_fixed`. Ignores `min_sigma`. |
| `'none'` | Do not compute finite standard deviations. If `std.export` is true, files are still written with `NaN` values. |

### Hutchinson fields

| Field | Default | Plain-language meaning |
| --- | ---: | --- |
| `cfg.std.hutch.probes` | `100000` | Number of random probe vectors. Larger is less noisy but slower. Set `0` to skip Hutchinson and use the fallback. |
| `cfg.std.hutch.probe_dist` | `'rademacher'` | Random probe distribution: `'rademacher'` (`+1/-1`) or `'gaussian'`. |
| `cfg.std.hutch.batch` | `250` | Number of probes solved at once. Larger can be faster but uses more memory. |
| `cfg.std.hutch.report_every_batch` | `10` | Print progress every this many batches. Set `0` to silence batch progress. |
| `cfg.std.hutch.max_nred` | `40000` | Skip Hutchinson if the reduced system is larger than this many unknowns; use fallback instead. |
| `cfg.std.hutch.min_diagrel` | `1e-12` | Small relative floor for the estimated inverse-Laplacian diagonal to avoid zero or negative Monte Carlo artifacts. |

Hutchinson error decreases slowly, approximately like `1/sqrt(probes)`. For very large station-phase components, `pseudo_degree` is often the practical choice.

---

## Common configuration examples

### Preserve microsecond or sub-microsecond values

```matlab
cfg.output.theta_decimals = 12;
cfg.output.thetastd_decimals = 12;
cfg.output.dt_decimals = 9;  % also increase dt_sync precision if needed
```

### Disable the exported `std_theta` residual-noise floor

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

The MATLAB and Julia implementations are intended to produce the same output formats. Small floating-point differences can occur because MATLAB and Julia use different sparse linear algebra libraries and random-number streams, especially for Hutchinson standard deviations. For deterministic format comparisons, use the same precision settings and prefer `std.mode = 'pseudo_degree'`.
