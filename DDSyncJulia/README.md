# DDSync (Julia)

DDSync reads a HypoDD/GrowClust-style `dt.cc` file, processes one station-phase group at a time, and writes a cleaned and synchronized differential-time file. For each station and phase, it estimates one relative arrival-time value, `theta`, for every event that can be constrained by the data. A kept observation between events `i` and `j` is then represented by `theta(i) - theta(j)`.

The MATLAB and Julia versions are intended to use the same defaults and produce the same file formats. This README uses Julia dictionary and TOML names. The MATLAB README contains the same guidance using MATLAB `cfg` names.

---

## Quick start

From inside the `DDSyncJulia/` folder:

```bash
cp ddsync_config_template.toml ddsync_config.toml
julia run_ddsync.jl
```

You can also pass the config file and input paths explicitly:

```bash
julia run_ddsync.jl ddsync_config.toml dt.cc catalog.txt
```

For a fast first run on a large file, disable Hutchinson standard deviations and use the degree proxy in the TOML file:

```toml
[std]
mode = "pseudo_degree"
fallback = "pseudo_degree"
```

---

## Requirements

- Julia.
- Only Julia standard libraries are used.
- Event IDs in the input should be positive integers. The output files are written for rows `1:maxEventID`.

---

## Directory layout

- `src/DDSync.jl` - maintained Julia implementation.
- `run_ddsync.jl` - command-line runner with optional TOML configuration.
- `ddsync_config_template.toml` - user-editable configuration template.
- `theta/`, `thetastd/` - output directories created at runtime.
- `ddsync_tmp/` - default scratch directory for per-group spool files and edge-decision files.

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

Event IDs must be suitable for array-style output: DDSync writes rows `1:maxEventID` in each `theta` and `std_theta` file. If your catalog uses arbitrary IDs, reindex before running.

### `catalog.txt`

DDSync uses `catalog.txt` only to find the largest event ID. The event ID must be the last whitespace-separated column on each non-empty line. If your catalog uses a different layout, edit the `readMaxEventID` helper in `src/DDSync.jl` accordingly.

---

## Output files

### `dt_sync.cc`

This is a HypoDD-style differential-time file with synchronized `dt` values. For each kept measurement between events `i` and `j`, DDSync writes:

```text
STA  theta(i)-theta(j)  weight  PH
```

The third column is controlled by `output.weight_mode`.

Pruned edges are omitted by default. If `output.write_pruned_edges = true`, DDSync keeps the original station line but writes its weight as `0`.

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

`std_theta` is either a Hutchinson estimate, a fast pseudo standard deviation, or `NaN`, depending on `std.mode` and fallback behavior. `degree` is the number of kept incident edges in the final graph. `thetaWeight` is a practical inverse-uncertainty weight derived from `std_theta`. `nodeWeight` is an optional diagnostic node-weight summary used by `pseudo_weight` mode.

### `sync_metrics.txt`

A tab-delimited summary. Each station-phase row has:

```text
group    n_edges  n_kept  pruned_pct  robust_p50  robust_p05  robust_p95  frac_robust_lt_0p99  sigma_hat
```

`group` is the station-phase key, `n_edges` is the number of input measurements, and `n_kept` is the number retained after component gating and pruning. `sigma_hat` is the final robust residual scale for reporting and for `std_theta` export. It is not the same thing as `std.min_sigma`.

---

## Precision for small times

The old fixed `%.6f` formatting could round away information when working at microsecond or sub-microsecond scale. DDSync now exposes separate precision controls:

```toml
[output]
theta_decimals = 9       # theta files, seconds
thetastd_decimals = 9    # std_theta files, seconds
```

These are digits after the decimal point. With seconds as the unit, 6 decimals is 1 microsecond, 9 decimals is 1 nanosecond, and 12 decimals is 1 picosecond. Increase these values for field-lab, acoustic, or other small-time experiments where very small differential times matter.

These settings affect only the `theta` and `std_theta` value columns. The synchronized `dt` column in `dt_sync.cc` is controlled by `output.dt_decimals`, and the weight column in `dt_sync.cc` is controlled by `output.dt_weight_decimals`.

Set both values to `6` if you need the older text-file formatting.

---

## Residual scale floors: `min_scale` versus `min_sigma`

DDSync has two separate safeguards. They intentionally do different jobs.

### `robust.min_scale`

This protects pruning and IRLS. During outlier detection and Huber reweighting, DDSync divides residuals by a robust scale estimate. If the residuals are exactly zero or nearly zero, the scale can become zero and the division becomes unstable. `robust.min_scale` prevents that. It can affect which edges are pruned and how robust weights are computed.

### `std.min_sigma` and `std.apply_min_sigma`

This applies only when converting the final residual-noise estimate into exported `std_theta` values for Hutchinson and `pseudo_degree` modes.

```toml
[std]
min_sigma = 5e-4
apply_min_sigma = true
```

When `std.apply_min_sigma = true`, DDSync uses at least `std.min_sigma` as the residual-noise scale before computing exported `std_theta`. When `std.apply_min_sigma = false` or `std.min_sigma = 0`, this export floor is disabled.

Example: suppose the final residual scale is `1e-5` seconds and `std.min_sigma = 5e-4`.

- With `std.apply_min_sigma = true`, `std_theta` export uses `5e-4` seconds before graph leverage or degree scaling.
- With `std.apply_min_sigma = false`, `std_theta` export uses `1e-5` seconds.

`min_sigma` does not floor every row in `std_theta_*.txt`. It floors the residual-noise scale before the graph/degree conversion. The final row values can be smaller or larger, and the pinned reference event remains `0`.

`pseudo_weight` mode does not use `sigma_hat`, so `min_sigma` does not affect `pseudo_weight` standard deviations.

---

## Configuration reference

All defaults come from `DDSync.config_default()`. You usually override them in `ddsync_config.toml`.

### `[io]` fields

| Field | Default | Plain-language meaning |
| --- | ---: | --- |
| `infile_dt` | `"dt.cc"` | Input differential-time file. |
| `catalog_file` | `"catalog.txt"` | Catalog file used to find `maxEventID`; event ID must be in the last column. |
| `out_dt_sync` | `"dt_sync.cc"` | Output synchronized differential-time file. |
| `metrics_file` | `"sync_metrics.txt"` | Output metrics summary. |
| `tmpdir` | `"ddsync_tmp"` | Scratch directory for per-group spool files and edge-decision files. Put this on fast local storage for large jobs. |
| `thetadir` | `"theta"` | Directory for `theta_<STA>_<PH>.txt`. |
| `thetastd_dir` | `"thetastd"` | Directory for `std_theta_<STA>_<PH>.txt`. |

### `[weights]` fields

| Field | Default | Plain-language meaning |
| --- | ---: | --- |
| `base_mode` | `"cc"` | How to turn the input `cc` value into the starting edge weight. Use `"cc"`, `"cc2"`, or `"ones"`. |

Examples:

```toml
[weights]
base_mode = "cc2"
```

DDSync never keeps edges whose base weight is non-positive.

### `[robust]` fields

| Field | Default | Plain-language meaning |
| --- | ---: | --- |
| `K_SIGMA` | `20.0` | Prune an edge when its residual is larger than this many robust scale units. Smaller values prune more aggressively. |
| `MIN_EDGES` | `30` | Minimum number of positive-weight edges required before a component is trusted. DDSync also requires at least one loop, meaning more measurements than a bare tree. |
| `min_scale` | `5e-4` | Scale floor used only inside pruning and IRLS to avoid divide-by-zero behavior. |

A practical interpretation: `MIN_EDGES` prevents tiny isolated families from producing overconfident `theta` values. `K_SIGMA` controls how far a measurement can disagree with the internally consistent solution before it is removed.

### `[irls]` fields

| Field | Default | Plain-language meaning |
| --- | ---: | --- |
| `iters` | `10` | Maximum Huber IRLS iterations after initial pruning. Set to `0` to disable robust reweighting. |
| `C_HUBER` | `1.345` | Huber tuning constant. Smaller values downweight moderate residuals more strongly. |
| `rel_tol` | `1e-3` | Stop IRLS early when the objective changes by less than this relative amount. |

### `[numeric]` fields

| Field | Default | Plain-language meaning |
| --- | ---: | --- |
| `RIDGE_EPS` | `1e-10` | Small diagonal ridge added to the reduced Laplacian solve for numerical stability. Increase only if solves are unstable; too large a value can slightly bias `theta`. |

### `[output]` fields

| Field | Default | Plain-language meaning |
| --- | ---: | --- |
| `write_dt_sync` | `true` | Write `dt_sync.cc`. Set false if you only need `theta`, `std_theta`, and metrics. |
| `weight_mode` | `"thetaStd"` | Weight written in the third column of `dt_sync.cc`: `"base"`, `"robust"`, `"combined"`, or `"thetaStd"`. |
| `dt_decimals` | `5` | Digits after the decimal for synchronized `dt` values in `dt_sync.cc`. |
| `dt_weight_decimals` | `4` | Digits after the decimal for the `dt_sync.cc` weight column. |
| `theta_decimals` | `9` | Digits after the decimal for `theta` values in `theta_*.txt`. |
| `thetastd_decimals` | `9` | Digits after the decimal for `std_theta` values in `std_theta_*.txt`. |
| `station_field_width` | `8` | Width used to align the station field in `dt_sync.cc`. |
| `dt_field_width` | `10` | Width used to align the `dt` field in `dt_sync.cc`. |
| `weight_field_width` | `8` | Width used to align the weight field in `dt_sync.cc`. |
| `write_pruned_edges` | `false` | If false, pruned station lines are omitted. If true, they are written with weight `0`. |
| `thetastd_scale_mode` | `"fixed"` | For `weight_mode="thetaStd"`: use `"fixed"` scale or compute a global `"median"` scale from kept observations. |
| `thetastd_scale_fixed` | `500.0` | Fixed divisor for thetaStd weights. Example: if `std_dt=0.002 s`, then `(1/std_dt)/500 = 1`. |
| `thetastd_weight_cap` | `1.0` | Maximum written thetaStd-derived weight. Use `Inf` in Julia config code, or a very large TOML value, to effectively disable the cap. |

`weight_mode` options:

| Mode | Written weight |
| --- | --- |
| `"base"` | Base weight from `cc`. |
| `"robust"` | Huber robust weight only. |
| `"combined"` | Base weight times robust weight. |
| `"thetaStd"` | Conservative inverse uncertainty for `theta(i)-theta(j)`, scaled and capped. |

For `thetaStd`, DDSync uses `sqrt(std_theta_i^2 + std_theta_j^2)` and ignores covariance. This is conservative for many relocation-weighting uses.

### `[std]` fields

| Field | Default | Plain-language meaning |
| --- | ---: | --- |
| `export` | `true` | Create the `thetastd/` directory and write `std_theta_*.txt`. |
| `mode` | `"hutch"` | Main method for `std_theta`: `"hutch"`, `"pseudo_degree"`, `"pseudo_weight"`, or `"none"`. |
| `fallback` | `"pseudo_degree"` | Method used when Hutchinson is skipped. Use `"pseudo_degree"`, `"pseudo_weight"`, or `"nan"`. |
| `min_sigma` | `5e-4` | Optional residual-noise floor used only for exported Hutchinson and `pseudo_degree` `std_theta`. |
| `apply_min_sigma` | `true` | Switch for `min_sigma`. Set false, or set `min_sigma` to `0`, to disable the export floor. |
| `pseudo_weight_source` | `"combined"` | Edge weight source for `pseudo_weight` and `nodeWeight`: `"combined"`, `"robust"`, or `"base"`. |
| `pseudo_weight_eps` | `1e-6` | Minimum node weight before inverting in `pseudo_weight` mode. Prevents division by zero. |
| `thetastd_write_weightcol` | `true` | Add the practical `thetaWeight` column to `std_theta_*.txt`. |
| `thetastd_write_alt_node_weightcol` | `false` | Add `nodeWeight` to `std_theta_*.txt`. If `thetastd_write_weightcol` is also true, this is the sixth column; otherwise it is the fifth. |
| `probes` | `100000` | Number of Hutchinson random probe vectors. Larger is less noisy but slower. Set `0` to skip Hutchinson and use the fallback. |
| `probe_dist` | `"rademacher"` | Random probe distribution: `"rademacher"` (`+1/-1`) or `"gaussian"`. |
| `batch` | `250` | Number of probes solved at once. Larger can be faster but uses more memory. |
| `report_every_batch` | `10` | Print progress every this many batches. Set `0` to silence batch progress. |
| `max_nred` | `40000` | Skip Hutchinson if the reduced system is larger than this many unknowns; use fallback instead. |
| `min_diag_rel` | `1e-12` | Small relative floor for the estimated inverse-Laplacian diagonal to avoid zero or negative Monte Carlo artifacts. |

Mode meanings:

| Mode | Meaning |
| --- | --- |
| `"hutch"` | Estimate the diagonal of the inverse graph Laplacian using random probes. Most statistical, but potentially expensive. |
| `"pseudo_degree"` | Fast proxy: `std_theta ~= sigma/sqrt(degree)`. Good for large jobs where relative weighting is enough. |
| `"pseudo_weight"` | Fast proxy based on incident edge weights, scaled by `thetastd_scale_fixed`. Ignores `min_sigma`. |
| `"none"` | Do not compute finite standard deviations. If `std.export` is true, files are still written with `NaN` values. |

Hutchinson error decreases slowly, approximately like `1/sqrt(probes)`. For very large station-phase components, `pseudo_degree` is often the practical choice.

### `[runtime]` fields

| Field | Default | Plain-language meaning |
| --- | ---: | --- |
| `prog_every_lines` | `2000000` | Print progress every this many input lines while streaming. |
| `print_new_groups` | `false` | Print a message whenever a new station-phase group is first seen. Useful for debugging input parsing. |
| `store_dense_group_arrays` | `false` | If true, store per-group theta/std arrays as dense vectors of length `maxEventID`. This is faster but can use more memory. If false, use dictionaries keyed by event ID. |
| `decision_chunk_records` | `200000` | Number of edge-decision records read per chunk during pass 2. Larger chunks can reduce I/O overhead but use more memory. |
| `median_scale_cap` | `2000000` | Maximum number of kept observations sampled when computing the global median thetaStd scale. |

---

## Common configuration examples

### Preserve microsecond or sub-microsecond values

```toml
[output]
theta_decimals = 12
thetastd_decimals = 12
dt_decimals = 9  # also increase dt_sync precision if needed
```

### Disable the exported `std_theta` residual-noise floor

```toml
[std]
apply_min_sigma = false
# or
min_sigma = 0
```

This does not change `robust.min_scale`, so pruning and IRLS remain protected from zero residual scales.

### Use fast practical standard deviations

```toml
[std]
mode = "pseudo_degree"
fallback = "pseudo_degree"
```

### Keep pruned lines for diagnostics

```toml
[output]
write_pruned_edges = true
```

Pruned station lines are then written with weight `0` instead of being omitted.

---

## Notes on MATLAB/Julia parity

The MATLAB and Julia implementations are intended to produce the same output formats. Small floating-point differences can occur because MATLAB and Julia use different sparse linear algebra libraries and random-number streams, especially for Hutchinson standard deviations. For deterministic format comparisons, use the same precision settings and prefer `std.mode = "pseudo_degree"`.
