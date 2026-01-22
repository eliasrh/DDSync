# DDSync

DDSync synchronizes differential travel times from a large `dt.cc` file by solving a **per station–phase graph synchronization** problem:

- Nodes: events
- Edges: observed differential travel time `dt(i,j)` for a given station and phase
- Output: a node potential `theta` such that `theta(i) - theta(j)` best matches all observed `dt(i,j)` (weighted, robust)

It is designed to be **scalable** to very large `dt.cc` files via streaming + per-(station,phase) processing.

Tested on Matlab2013b

---

## Quick start

1. Put this folder on your MATLAB path (add the parent folder that contains `+ddsync/`).
2. Create a parameter file:

- Copy `ddsync_params_template.m` to `ddsync_params.m`
- Edit the input/output paths

3. Run:

```matlab
run_ddsync
```

Or directly:

```matlab
cfg = ddsync.config_default();
cfg.io.infile_dt    = 'dt.cc';
cfg.io.catalog_file = 'catalog.txt';
out = ddsync.run(cfg);
```

---

## Inputs

### 1) `dt.cc` format

DDSync expects a hypoDD/GrowClust-like `dt.cc` file with blocks:

```
# i j 0.0
STA  dt  cc  PH
STA  dt  cc  PH
...
# i j 0.0
...
```

- `i`, `j` are integer event IDs (1-based).
- Data lines contain:
  - `STA` station code (string)
  - `dt` differential time (seconds)
  - `cc` correlation coefficient or quality value (float)
  - `PH` phase label (string, e.g. `P`, `S`)

### 2) `catalog.txt`

DDSync only uses the catalog to determine `maxEventID`. The catalog must contain event IDs in the last column.

---

## Outputs

### 1) `theta/theta_<STA>_<PH>.txt`

Per station–phase node potentials (one value per event ID):

Columns:
1. `EventID`
2. `theta` (seconds) — NaN if event not present in the group
3. `refEventID` (the reference event ID for that connected component)

### 2) `thetastd/std_theta_<STA>_<PH>.txt`

Per station–phase uncertainty/weight metadata.

Default columns:
1. `EventID`
2. `std_theta` (seconds) — may be true (Hutch) or pseudo (fallback)
3. `refEventID`
4. `degree` (graph degree within the kept component)

Optional extra columns (enabled by default):
5. `thetaWeight = min(cap, (1/std_theta)/theta_weight_scale_fixed)`

Optional extra column (disabled by default):
6. `nodeWeight` — aggregated node weight used by pseudo-weight mode.

### 3) `dt_sync.cc`

A pruned + synchronized `dt.cc` file (same structure) where:
- outliers are removed (weight written as 0)
- `dt` is replaced by `theta(i)-theta(j)` for kept edges
- the third column is a **user-selected weight** via `cfg.output.dt_weight_mode`

### 4) `sync_metrics.txt`

A small metrics file including total edges/kept and whether pseudo-std was used.

---

## Standard deviation / weight export

True uncertainty estimation of `theta` can be expensive for large components. DDSync supports:

### Hutchinson estimator (`cfg.std.mode = 'hutch'`)
Estimates `diag(inv(L))` stochastically for each connected component of a station–phase graph, then converts to `std(theta)`.

This is controlled by:
- `cfg.std.hutch.probes`
- `cfg.std.hutch.max_nred` (skips large components)

### Pseudo-std fallback (NEW)
If Hutch is skipped (too large / disabled), DDSync can write a **pseudo std** instead of NaNs:

- `cfg.std.fallback = 'pseudo_degree'` (default)
  - `std_theta ~= sigma_hat / sqrt(degree)` (fast, stable)

- `cfg.std.fallback = 'pseudo_weight'`
  - derive a per-node weight from edge weights and map
    `std_theta = 1/(nodeWeight * theta_weight_scale_fixed)`

DDSync prints warnings when pseudo std is used and records counts in `sync_metrics.txt`.

> Note: Pseudo std is not a statistical uncertainty; it is intended as a practical weighting proxy when true std is unavailable.

---

## Key configuration fields (high level)

- `cfg.io.*` : input/output paths, temp dir, progress
- `cfg.weights.*` : base weight from CC (`cc`, `cc^2`, or uniform)
- `cfg.robust.*` : pruning + IRLS (Huber) controls
- `cfg.std.*` : thetastd export mode and Hutch/pseudo settings
- `cfg.output.*` : dt_sync weight mode and scaling

See `+ddsync/config_default.m` for full details.

---

## Performance notes

- DDSync performs:
  1) Pass 1: stream `dt.cc` and spool each station–phase group to a binary file
  2) Process each group independently (memory ~ size of one group)
  3) Pass 2: stream `dt.cc` again and write `dt_sync.cc`

- Place `cfg.io.tmpdir` on a fast local SSD for large catalogs.

- If Hutch is too slow or too memory intensive:
  - set `cfg.std.mode = 'pseudo_degree'` (fastest)
  - or reduce `cfg.std.hutch.probes`

---

## Contact / citation

This repository is intended to accompany a DDSync manuscript in preparation. Please cite the associated publication/preprint when available.

