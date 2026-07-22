# DDSync (MATLAB + Julia)

DDSync is a graph-based denoising and synchronization tool for differential travel-time pair files. It reads a HypoDD/GrowClust-style `dt.cc` file, solves each station-phase group independently, and writes synchronized differential times plus per-event `theta` and optional `std_theta` metadata.

There are two maintained implementations:

- `DDSync/` - MATLAB package implementation.
- `DDSyncJulia/` - Julia implementation with TOML configuration and no non-stdlib package dependencies.

The MATLAB and Julia implementations are intended to match in algorithm, default behavior, and output file formats. The implementation-specific README files document the same behavior using the names used by that language.

---

## What DDSync writes

A typical run writes:

- `dt_sync.cc` - synchronized and optionally pruned differential-time file.
- `theta/theta_<STA>_<PH>.txt` - one relative arrival-time value, `theta`, for each event in each station-phase group.
- `thetastd/std_theta_<STA>_<PH>.txt` - standard-deviation or practical-weight metadata for each event.
- `sync_metrics.txt` - per-group and overall diagnostics.

For a kept observation between events `i` and `j`, the synchronized differential time is `theta(i) - theta(j)`.

---

## Which implementation should I use?

| Situation | Suggested implementation |
| --- | --- |
| You already work in MATLAB or want the original reference package | `DDSync/` |
| You want a command-line workflow with TOML configuration | `DDSyncJulia/` |
| You want lower memory use during pass 2 on large catalogs | `DDSyncJulia/` with `runtime.store_dense_group_arrays = false` |
| You want to run the bundled Spanish Springs example immediately | `DDSync/` |

The Julia folder does not include example `dt.cc` and `catalog.txt` files, but it reads the same file formats as MATLAB.

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

The catalog file is used to find `maxEventID`. The event ID must be the last whitespace-separated column. DDSync writes full `1:maxEventID` rows in `theta` and `std_theta` files, so event IDs should be reindexed to a contiguous positive integer range before running.

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

## Lastest Update: Shared `min_scale` versus `min_sigma` behavior

DDSync has two separate small-scale safeguards. They are deliberately separate.

`robust.min_scale` is used during robust pruning and Huber IRLS. It prevents divide-by-zero behavior when residuals are exactly zero or nearly zero. It can affect which edges are pruned and how robust weights are computed.

`std.min_sigma` is used only when exporting `std_theta`, and only when `std.apply_min_sigma = true`. It floors the final residual-noise scale before Hutchinson or `pseudo_degree` standard deviations are written.

Disabling `std.apply_min_sigma` disables the exported `std_theta` residual-noise floor. It does not disable the robust pruning/IRLS floor.

`min_sigma` does not say that every row in `std_theta_*.txt` must be at least `min_sigma`. It floors the residual-noise estimate before graph leverage or degree scaling. Individual rows can be smaller or larger, and the pinned reference event remains `0`.

`pseudo_weight` mode does not use `sigma_hat`, so `min_sigma` does not affect `pseudo_weight` standard deviations.

---

## Lastest update: Shared precision behavior

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

The synchronized `dt` column in `dt_sync.cc` has a separate precision setting: `cfg.output.dt_decimals` in MATLAB and `dt_decimals` in Julia/TOML.

---

## References

- Elías Rafn Heimisson, Yifan Yu; *DDSync: Graph‐Based Denoising of Differential Travel‐Time Observations with Applications to Pick Reconstruction and Path‐Difference Tomography.* **Seismological Research Letters** 2026; doi: https://doi.org/10.1785/0220260086 
- Trugman, D. T., & Shearer, P. M. (2017). *GrowClust: A hierarchical clustering algorithm for relative earthquake relocation, with application to the Spanish Springs and Sheldon, Nevada, earthquake sequences.* **Seismological Research Letters**, 88(2A), 379-391.

---

## License (non-commercial)

DDSync is licensed for non-commercial research and educational use only. Commercial or for-profit use, including internal commercial workflows, requires a separate written licensing agreement with the author (eliasrafn@hi.is). See `LICENSE` for full terms.
