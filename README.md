# DDSync (Julia and MATLAB)

DDSync is a graph-based denoising and synchronization tool for differential travel-time pair files. It reads a HypoDD/GrowClust-style `dt.cc` file, solves each station-phase group independently, and writes synchronized differential times plus per-event `theta` and optional `std_theta` metadata.

The repository contains two implementations:

- **`DDSyncJulia/` — primary implementation**, with TOML configuration and only Julia standard-library dependencies.
- `DDSync/` — MATLAB implementation, maintained for bug fixes.

The synchronization algorithm and its defaults are shared. Input-ID handling
differs: Julia now maps arbitrary integer IDs automatically, while MATLAB
still requires users to prepare compact positive IDs before running.

## Update (19 November 2026): visible uncertainty quality

`sync_metrics.txt` now gives counts and percentages for exact, converged Hutch, kept unconverged Hutch, and fallback std values, globally and by station-phase group. It advises further sampling or exact repairs when at least 10% require review. Observation columns are aligned, and component details are grouped at the end. The default maximum probe budget is now **4096**; explicit older settings remain in effect.

The report also distinguishes raw residual scale `sigma_hat` from the optionally floored `sigma_std` used in uncertainty formulas. Read the [reporting update](docs/updates/2026-11-19-uncertainty-reporting.md) and [uncertainty guide](docs/uncertainty.md).

## Update (12 November 2026): preserve useful Hutchinson estimates

Uncertainty decisions are now made per event. Small components use exact calculation directly; larger ones retain converged Hutchinson estimates and use selected full-system exact solves to repair troublesome entries. Remaining usable positive estimates can be retained with an explicit quality label. The default final fallback is `pseudo_degree`; graph bounds are now opt-in. Numeric theta/std file layouts are unchanged, and new `std_quality_*.txt` files identify each event's method.

Read the [update and large-component test report](docs/updates/2026-11-12-per-event-uncertainty.md) and the [uncertainty guide](docs/uncertainty.md) for defaults, limits and migration. Degree fallback is a practical weighting proxy and can underestimate uncertainty across weak connections.

The [October update](docs/updates/2026-10-18-pruning-uncertainty.md) fixed references after pruning. That fix remains: never connect theta stars across different reference IDs. November supersedes its whole-component graph-bound fallback policy.

## Try the bundled example in Julia

From the repository root, with Julia 1.10 or later installed:

```bash
mkdir example-run
cd example-run
julia ../DDSyncJulia/run_ddsync.jl ../DDSync/dt.cc ../DDSync/catalog.txt
```

Outputs are written into `example-run/`. The example input is already synchronized, so large changes are not expected. For your own data and configuration, start with the [Julia manual](DDSyncJulia/README.md). Relative paths are resolved from the working directory; use separate output directories for separate runs.

## Recent update (28 September 2026): automatic event-ID mapping in Julia

DDSyncJulia now accepts integer event IDs without requiring them to fit in a
machine integer or form a sequential range. When the catalog IDs do not form
`1:N`, it maps them to that range in catalog row order. The mapping is printed
as a notice, recorded in the metrics, and saved in `event_id_map.csv`.

The usual `dt_sync.cc` retains the original IDs. The additional
`dt_sync_seq.cc` and `catalog_seq.txt` use the new IDs, as do both ID columns
in the `theta` and `std_theta` files. Use **`catalog_seq.txt` with those theta
directories in GraphSplit**. No new TOML setting is required. Inputs already
using `1:N` retain their IDs and do not need the extra files.

The September indexing change did not alter MATLAB input-ID handling. Its IDs must still be reindexed
before processing; automatic mapping is available only in Julia. The
[Julia manual](DDSyncJulia/README.md#event-ids-and-graphsplit) explains file
selection, custom output paths, and how to retain the original correspondence.

---

## What DDSync writes

A typical run writes:

- `dt_sync.cc` - synchronized and optionally pruned differential-time file.
- `theta/theta_<STA>_<PH>.txt` - one relative arrival-time value, `theta`, for each event in each station-phase group.
- `thetastd/std_theta_<STA>_<PH>.txt` - standard-deviation or practical-weight metadata for each event.
- `sync_metrics.txt` - aligned observation statistics, uncertainty percentages and review advice, followed by component details.

For a kept observation between events `i` and `j`, the synchronized differential time is `theta(i) - theta(j)`.

---

## Which implementation should I use?

| Situation | Suggested implementation |
| --- | --- |
| You are starting a new workflow | `DDSyncJulia/` |
| You already work in MATLAB and need the supported bug fixes | `DDSync/` |
| You want a command-line workflow with TOML configuration | `DDSyncJulia/` |
| You want lower memory use during pass 2 on large catalogs | `DDSyncJulia/` with `runtime.store_dense_group_arrays = false` |
| You want to try the bundled Spanish Springs example | Either implementation; Julia command above |

Both implementations read the example data in `DDSync/`; there is no need to duplicate the input files.

---

## Input assumptions

DDSync expects `dt.cc` blocks like:

```text
# i j 0.0
STA  dt  cc  PH
STA  dt  cc  PH
# i j 0.0
...
```

The header gives the event pair `i j`; the third header number is ignored. Each station line gives the station code, differential time in seconds, correlation or quality value, and phase label.

The event ID must be the last whitespace-separated catalog column. Each
catalog ID must be unique, and every pair endpoint in `dt.cc` must occur in
the catalog. In Julia, write IDs as decimal integers, without decimal points
or scientific notation. They are read exactly, even when their magnitude
exceeds the range of a 64-bit integer.

Starting with IDs `1:N` is still recommended: it avoids the mapping and extra
output files. Catalog rows do not have to be sorted by ID. MATLAB requires
this compact positive-ID convention because it writes rows `1:maxEventID`
and stores IDs in single precision during temporary processing. Its IDs
must not exceed 16,777,216. Julia applies this limit to the **number of
catalog events**, rather than to the external IDs.

---

## Implementation-specific documentation

Start here:

- [MATLAB README](DDSync/README.md)
- [Julia README](DDSyncJulia/README.md)

Both include quick starts, input/output formats, precision settings, `min_sigma` behavior, standard-deviation modes, examples, and a full parameter reference.

---

## Example data

The MATLAB folder includes example `dt.cc` and `catalog.txt` from the Spanish Springs sequence based on Trugman and Shearer (2017). The event IDs have been made sequential and lightly filtered for demonstration. Note that the dt.cc file provided here is a synchronized version of the one provided here by Trugman: https://github.com/dttrugman/GrowClust3D.jl/tree/master/examples/data/in as xcordata.txt. The user can thus compare to the original file and check if they are running the code correctly. Resynchronizing will not do much to the file unless input parameters are very different.

---

## Shared `min_scale` versus `min_sigma` behavior

DDSync has two separate small-scale safeguards. They are deliberately separate.

`robust.min_scale` is used during robust pruning and Huber IRLS. It prevents divide-by-zero behavior when residuals are exactly zero or nearly zero. It can affect which edges are pruned and how robust weights are computed.

`std.min_sigma` is used only when exporting `std_theta`, and only when `std.apply_min_sigma = true`. It floors the final residual-noise scale before Hutchinson, exact or `pseudo_degree` standard deviations are written.

Disabling `std.apply_min_sigma` disables the exported `std_theta` residual-noise floor. It does not disable the robust pruning/IRLS floor.

`min_sigma` does not say that every row in `std_theta_*.txt` must be at least `min_sigma`. It floors the residual-noise estimate before graph leverage or degree scaling. Individual rows can be smaller or larger, and the pinned reference event remains `0`.

`sigma_hat` is the raw component residual-noise estimate; `sigma_std` is its optionally floored value used in the uncertainty formulas. The group-table `sigma_hat` is only a reporting median across components.

`pseudo_weight` mode does not use `sigma_hat` or `sigma_std`, so `min_sigma` does not affect `pseudo_weight` standard deviations.

---

## Shared precision behavior

Both implementations expose formatting precision for the `theta` and `std_theta` value columns. The default is 9 digits after the decimal point.

MATLAB:

```matlab
cfg.output.theta_decimals = 9;
cfg.output.thetastd_decimals = 9;
```

Julia/TOML:

```toml
[output]
theta_decimals = 9
thetastd_decimals = 9
```

These settings matter when differential times can be microseconds or smaller. The older fixed `%.6f` output has only 6 digits after the decimal point. With seconds as the unit, 6 decimals is 1 microsecond, 9 decimals is 1 nanosecond, and 12 decimals is 1 picosecond. Set these fields to `12` or higher for sub-microsecond applications, or set them to `6` to reproduce older formatting.

Extra printed digits do not increase numerical accuracy: temporary observations and pass-2 theta/std storage use single precision. The synchronized `dt` column in `dt_sync.cc` has a separate precision setting: `cfg.output.dt_decimals` in MATLAB and `dt_decimals` in Julia/TOML.

---

## Tests

From the repository root:

```bash
julia --startup-file=no DDSyncJulia/test/runtests.jl
```

In MATLAB: `addpath('DDSync/test'); run_tests`. Both suites use shared synthetic data and check pruning, references, uncertainty failure handling and export weights. Julia also tests arbitrary event-ID mapping. The [validation report](docs/updates/2026-10-18-pruning-uncertainty.md) gives the separate probe-study commands and results.

---

## References

- Elías Rafn Heimisson. GraphSplit: Sparse Double-Difference Earthquake Relocation From Synchronized Differential-Time Graphs. ESS Open Archive. 02 October 2026. [doi:10.22541/essoar.15009811/v1](https://doi.org/10.22541/essoar.15009811/v1).
Submitted to *Journal of Geophysical Research: Solid Earth*

- Elías Rafn Heimisson, Yifan Yu; *DDSync: Graph‐Based Denoising of Differential Travel‐Time Observations with Applications to Pick Reconstruction and Path‐Difference Tomography.* **Seismological Research Letters** 2026; doi: https://doi.org/10.1785/0220260086
  
- Trugman, D. T., & Shearer, P. M. (2017). *GrowClust: A hierarchical clustering algorithm for relative earthquake relocation, with application to the Spanish Springs and Sheldon, Nevada, earthquake sequences.* **Seismological Research Letters**, 88(2A), 379-391.

---

## License (non-commercial)

DDSync is licensed for non-commercial research and educational use only. Commercial or for-profit use, including internal commercial workflows, requires a separate written licensing agreement with the author (eliasrafn@hi.is). See `LICENSE` for full terms.
