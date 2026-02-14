# DDSync (Julia)

This is the **Julia** implementation of the DDSync standalone MATLAB package. It performs two-pass streaming over large `dt.cc` files and solves per-(station,phase) graph synchronization problems.

Goals:

- **Easy setup**: no non-stdlib Julia packages.
- **Memory-safe** for huge `dt.cc` files: stream + spool; do not load the whole `dt.cc` into RAM.
- **Output formats** intended to match the MATLAB standalone package:
  - `dt_sync.cc` (HypoDD-style)
  - `theta/theta_STA_PH.txt`
  - `thetastd/std_theta_STA_PH.txt`
  - `sync_metrics.txt`

> **Note:** This Julia version is generally faster than MATLAB, and it provides a dense vs. sparse storage option for pass-2 output writing.

---

## Requirements

- Julia **1.9+** recommended (should work on 1.8+)
- No extra packages (only stdlib)

---

## Quick start (minimal)

Place `dt.cc` and `catalog.txt` in the same directory, then run:

```bash
julia run_ddsync.jl
```

Or provide paths explicitly:

```bash
julia run_ddsync.jl path/to/dt.cc path/to/catalog.txt
```

Outputs:

- `dt_sync.cc`
- `sync_metrics.txt`
- `theta/`
- `thetastd/` (if `cfg[:std][:export] = true`)
- temporary spool/decision files under `synchro_tmp/`

---

## Configuration (fool-proof)

Preferred: use a **TOML config file**. Copy `ddsync_config_template.toml` to `ddsync_config.toml` and edit.
`run_ddsync.jl` will automatically load `ddsync_config.toml` if it exists, or you can pass a TOML file explicitly:

```bash
cp ddsync_config_template.toml ddsync_config.toml
julia run_ddsync.jl myconfig.toml
julia run_ddsync.jl myconfig.toml dt.cc catalog.txt
```

You can also modify the dictionary returned by `DDSync.config_default()` directly in `run_ddsync.jl`.

### Key options (plain language)

These mirror the MATLAB standalone config defaults:

- `cfg[:io]` — input/output paths, temp directory, progress logging
- `cfg[:weights][:base_mode]` — base weights from the `cc` column (`"cc"`, `"cc2"`, or `"ones"`)
- `cfg[:robust]` — pruning threshold (`K_SIGMA`) and minimum edges per component (`MIN_EDGES`)
- `cfg[:irls]` — IRLS iteration count and Huber constant (`C_HUBER`)
- `cfg[:std]` — standard-deviation export: Hutchinson vs. pseudo
- `cfg[:output]` — which weight is written to `dt_sync.cc` and formatting options

### Dense vs. sparse storage (memory policy)

The MATLAB package stores per-group `theta`/`std` as **dense arrays** for fast pass-2 output. This can be memory-heavy for large catalogs.

In Julia you can choose:

- `cfg[:runtime][:store_dense_group_arrays] = true`
  - MATLAB-like: fast pass 2, high RAM usage
- `cfg[:runtime][:store_dense_group_arrays] = false` *(default)*
  - Memory-light: stores `Dict{EventID => value}` per group, slower pass 2

---

## Input files (explicit formats)

### 1) `dt.cc`

HypoDD/GrowClust-style blocks:

```
# i j 0.0
STA  dt  cc  P
STA  dt  cc  S
...
# i j 0.0
...
```

- `i`, `j` are integer event IDs.
- The `0.0` in the header line is ignored.
- Only the first 4 whitespace-separated fields on station lines are used.

**Important:** Event IDs must be contiguous `1..N`. Reindex if needed (see the MATLAB `extras/` scripts for a reference approach).

### 2) `catalog.txt`

DDSync only uses the catalog to determine `maxEventID`. The catalog must contain event IDs in the **last column**.

---

## Output files

### `theta/theta_STA_PH.txt`

```
EventID theta refEventID
```

### `thetastd/std_theta_STA_PH.txt`

Tab-separated (optional columns controlled by config):

```
EventID   std_theta   refEventID   degree   [w_theta]   [node_weight]
```

### `dt_sync.cc`

A synchronized `dt.cc` with pruned outliers removed. The third column is the selected weight mode (`base`, `robust`, `combined`, or `thetaStd`).

### `sync_metrics.txt`

Summary metrics per station–phase group and overall counts.

---

## Notes / limitations

- Hutchinson (`cfg[:std][:mode] = "hutch"`, `cfg[:std][:probes] > 0`) can be **computationally expensive**.
- For very large groups, processing a single station-phase can still require substantial RAM (edges for that group are held in memory).

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
