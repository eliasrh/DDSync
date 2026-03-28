# DDSync (MATLAB)

DDSync synchronizes differential travel times from a large `dt.cc` file by solving a **per-station / per-phase graph synchronization** problem:

- Nodes: events
- Edges: observed differential travel time `dt(i,j)` for a given station and phase
- Output: node potential `theta` such that `theta(i) - theta(j)` best matches the observed `dt(i,j)` (weighted, robust)

DDSync is designed to be **scalable** to very large `dt.cc` files via streaming + per-(station,phase) processing (two-pass I/O).

Tested on MATLAB R2023b (should work on newer versions).

---

## Quick start (minimal)

1. Add this folder to your MATLAB path (add the parent folder that contains `+ddsync/`).
2. Copy the parameter template:
   - `ddsync_params_template.m` → `ddsync_params.m`
3. Edit `ddsync_params.m` with your input/output paths.
4. Run:

```matlab
run_ddsync
```

Or directly in MATLAB:

```matlab
cfg = ddsync.config_default();
cfg.io.infile_dt    = 'dt.cc';
cfg.io.catalog_file = 'catalog.txt';
out = ddsync.run(cfg);
```

---

## Requirements

- MATLAB (tested on R2013b; expected to work on newer versions).
- No extra toolboxes are required. DDSync uses base MATLAB functionality (sparse matrices, `chol`, `pcg`, `containers.Map`, etc.).

---

## Directory layout

- `+ddsync/` — core package code.
- `run_ddsync.m` — wrapper that loads `ddsync_params.m` and calls DDSync.
- `run_ddsync_example.m` — example script with common options (commented).
- `ddsync_params_template.m` — parameter file template.
- `dt.cc`, `catalog.txt` — example data (Spanish Springs sequence; Trugman & Shearer, 2017).
- `sync_metrics.txt` — example metrics output.
- `extras/` — helper scripts (reindexing/reordering, compressed variant).

---

## Input files (explicit formats)

### 1) `dt.cc`

DDSync expects a HypoDD/GrowClust-style `dt.cc` file with blocks like:

```
# i j 0.0
STA  dt  cc  PH
STA  dt  cc  PH
...
# i j 0.0
...
```

- `i`, `j` are **integer event IDs**.
- The `0.0` in the header line is ignored (can be any float).
- Data line fields (space-separated):
  - `STA` — station code (string)
  - `dt` — differential time (seconds)
  - `cc` — correlation coefficient or quality value (float)
  - `PH` — phase label (string, e.g., `P`, `S`)

**Important:** Event IDs must be contiguous `1..N`. If your catalog uses arbitrary IDs, reindex them first (see `extras/reindext_dtcc.m` and `extras/reorder.m`).

### 2) `catalog.txt`

DDSync only uses the catalog to determine `maxEventID`. The catalog must contain event IDs in the **last column**. If your catalog format differs, edit `ddsync.readCatalogIDs` inside `+ddsync/` to match your format.

---

## Output files

### 1) `theta/theta_<STA>_<PH>.txt`

Per station–phase node potentials (one value per event ID). Columns:

1. `EventID`
2. `theta` (seconds) — NaN if event not present in the group
3. `refEventID` (reference event ID for that connected component)

### 2) `thetastd/std_theta_<STA>_<PH>.txt`

Per station–phase uncertainty/weight metadata. Default columns:

1. `EventID`
2. `std_theta` (seconds) — true Hutchinson or pseudo (fallback)
3. `refEventID`
4. `degree` (graph degree within the kept component)

Optional columns (enabled by default):

5. `thetaWeight = min(cap, (1/std_theta) / theta_weight_scale_fixed)`

Optional column (disabled by default):

6. `nodeWeight` — aggregated node weight used by pseudo-weight mode.

### 3) `dt_sync.cc`

A pruned + synchronized `dt.cc` file (same block structure) where:

- outliers are removed (weight written as 0 if `write_pruned_edges = true`)
- `dt` is replaced by `theta(i) - theta(j)` for kept edges
- the third column is a **user-selected weight** via `cfg.output.dt_weight_mode`

### 4) `sync_metrics.txt`

A metrics summary with per-station/phase statistics and overall counts.

Header fields (second line) echo key config values. Each data row is:

```
group    n_edges  n_kept  pruned_pct  robust_p50  robust_p05  robust_p95  frac_robust_lt_0p99  sigma_hat
```

- `group` — station + phase key (e.g., `DYN_P`)
- `n_edges` / `n_kept` — total edges and edges kept after pruning
- `pruned_pct` — percent of edges removed
- `robust_p50/p05/p95` — robust weight percentiles (after IRLS)
- `frac_robust_lt_0p99` — fraction of edges downweighted below 0.99
- `sigma_hat` — robust scale estimate (seconds)

The end of the file includes `overall_*` totals and flags indicating whether pseudo-std was used.

---

## Standard deviation / weight export

True uncertainty estimation of `theta` can be expensive for large components. DDSync supports:

### Hutchinson estimator (`cfg.std.mode = 'hutch'`)

Estimates `diag(inv(L))` stochastically for each connected component of a station–phase graph, then converts to `std(theta)`.

Controlled by:

- `cfg.std.hutch.probes`
- `cfg.std.hutch.max_nred` (skips large components)

### Pseudo-std fallback

If Hutch is skipped (too large / disabled), DDSync can write a **pseudo std** instead of NaNs:

- `cfg.std.fallback = 'pseudo_degree'` (default)
  - `std_theta ~= sigma_hat / sqrt(degree)` (fast, stable)

- `cfg.std.fallback = 'pseudo_weight'`
  - derive a per-node weight from edge weights and map
    `std_theta = 1/(nodeWeight * theta_weight_scale_fixed)`

DDSync prints warnings when pseudo std is used and records counts in `sync_metrics.txt`.

> Note: Pseudo std is not a statistical uncertainty; it is intended as a practical weighting proxy when true std is unavailable.

---

## Key configuration fields (plain-language summary)

- `cfg.io.*` — input/output paths, temp dir, progress logging.
- `cfg.weights.*` — base weights derived from `cc` (`cc`, `cc^2`, or uniform).
- `cfg.robust.*` — pruning threshold (`k_sigma`), minimum edges, IRLS controls.
- `cfg.std.*` — whether to export std, and which estimator (Hutchinson vs. pseudo).
- `cfg.output.*` — formatting + which weight is written to `dt_sync.cc`.

See `+ddsync/config_default.m` for all defaults and descriptions.

---

## Example: richer configuration (all optional)

```matlab
cfg = ddsync.config_default();

% I/O
cfg.io.infile_dt    = 'dt.cc';
cfg.io.catalog_file = 'catalog.txt';
cfg.io.out_dt_sync  = 'dt_sync.cc';
cfg.io.metrics_file = 'sync_metrics.txt';
cfg.io.tmpdir       = 'ddsync_tmp'; % put on fast local disk for large jobs

% Weights
% cfg.weights.base_fun = 'cc2'; % use cc^2 instead of cc

% Robust pruning / IRLS
% cfg.robust.k_sigma    = 8;    % prune if |residual| > K_SIGMA * sigma_hat
% cfg.robust.min_edges  = 30;   % minimum edges per component
% cfg.robust.irls_iters = 10;   % IRLS iterations (0 disables)

% Std/weight export
% cfg.std.mode     = 'pseudo_degree'; % fastest; avoids Hutchinson
% cfg.std.fallback = 'pseudo_degree';

% Output weights in dt_sync.cc
% cfg.output.dt_weight_mode = 'thetaStd'; % 'base'|'robust'|'combined'|'thetaStd'
```

---

## Performance notes

- DDSync performs:
  1) Pass 1: stream `dt.cc` and spool each station–phase group to a binary file
  2) Process each group independently (memory ~ size of one group)
  3) Pass 2: stream `dt.cc` again and write `dt_sync.cc`

- Place `cfg.io.tmpdir` on a fast local SSD for large catalogs.

- If Hutchinson is too slow or too memory intensive:
  - set `cfg.std.mode = 'pseudo_degree'` (fastest)
  - or reduce `cfg.std.hutch.probes`

---

## Example data

The example Spanish Springs dataset (from Trugman & Shearer, 2017 ; https://github.com/dttrugman/GrowClust) is included only in the MATLAB folder (`DDSync/`). This Julia package reads the same `dt.cc` and `catalog.txt` formats but does not bundle the files.

---

## References

- Heimisson, E. R., Yu, Y. (2026, in preparation). *Graph-based denoising of differential travel-time observations with applications to pick reconstruction and path-difference tomography.*
- Trugman, D. T., & Shearer, P. M. (2017). *GrowClust: A hierarchical clustering algorithm for relative earthquake relocation, with application to the Spanish Springs and Sheldon, Nevada, earthquake sequences.* **Seismological Research Letters**, 88(2A), 379–391.

---

## License (non-commercial)

DDSync is licensed for **non-commercial** research and educational use only. Commercial or for-profit use (including internal commercial workflows) requires a separate written licensing agreement with the author (eliasrafn@hi.is). See `LICENSE` for full terms.

