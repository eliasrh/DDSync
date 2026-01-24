# DDSync (Julia)

This is a **Julia** implementation of the DDSync standalone MATLAB package (two-pass streaming, per-(station,phase) spooling).

Goals:

- **Easy setup**: no non-stdlib Julia packages.
- **Memory-safe** on huge `dt.cc` files: stream + spool; do not load the whole `dt.cc` into RAM.
- **Output formats** intended to match the MATLAB standalone package:
  - `dt_sync.cc` (HypoDD-style)
  - `theta/theta_STA_PH.txt`
  - `thetastd/std_theta_STA_PH.txt`
  - `sync_metrics.txt`

> **Note:** This is intended to be a readable, functional reference implementation rather than the fastest possible code.

## Requirements

- Julia (recommended: 1.9+; should work on 1.8+)
- No extra packages

## Quick start

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

## Configuration

Preferred: use a **TOML config file**. Copy `ddsync_config_template.toml` to `ddsync_config.toml` and edit.
`run_ddsync.jl` will automatically load `ddsync_config.toml` if it exists, or you can pass a TOML file as the first argument:

```bash
cp ddsync_config_template.toml ddsync_config.toml 
```
Edit the configuration file as needed
```bash
julia run_ddsync.jl myconfig.toml
julia run_ddsync.jl myconfig.toml dt.cc catalog.txt
```

You can also modify the dictionary returned by `DDSync.config_default()` directly in `run_ddsync.jl`.

Key options mirror the MATLAB standalone config:

- `cfg[:weights][:base_mode]` = `"cc"` | `"cc2"` | `"ones"`
- `cfg[:robust][:K_SIGMA]`
- `cfg[:robust][:MIN_EDGES]`
- `cfg[:irls][:iters]`, `cfg[:irls][:C_HUBER]`
- `cfg[:output][:weight_mode]` = `"base"` | `"robust"` | `"combined"` | `"thetaStd"`
- `cfg[:std][:mode]` = `"hutch"` | `"pseudo_degree"` | `"pseudo_weight"` | `"none"`
- `cfg[:std][:probes]` (set >0 to enable Hutchinson; can be slow)

### Memory policy

The MATLAB standalone package stores full `theta_full`/`std_full` arrays per group for fast pass-2 writing, but this can be very memory-heavy for large catalogs.

In Julia you can choose:

- `cfg[:runtime][:store_dense_group_arrays] = true`  
  (MATLAB-like: fast in pass 2, potentially large RAM use)

- `cfg[:runtime][:store_dense_group_arrays] = false` *(default)*  
  Store sparse `Dict{EventID => value}` per group (lower RAM, slower pass 2)

## Input format (dt.cc)

HypoDD style:

```
# i j 0.0
STA  dt  cc  P
STA  dt  cc  S
...
# i j 0.0
...
```

Only the first 4 whitespace-separated fields on station lines are used.

## Output formats

`dt_sync.cc` station-line formatting follows the MATLAB standalone package (field widths/decimals controlled by `cfg[:output]`).

Theta file:

```
EventID theta refEventID
```

Std-theta file (tab-separated; optional weight column):

```
EventID   std_theta   refEventID   degree   [w_theta]   [node_weight]
```

## Notes / limitations

- Hutchinson (`cfg[:std][:mode]="hutch"`, `cfg[:std][:probes]>0`) is **computationally expensive**.
- For very large groups, processing a single station-phase can still require substantial RAM (edges must be loaded for that group).
