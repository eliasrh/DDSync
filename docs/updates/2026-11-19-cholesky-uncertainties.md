# A second way to estimate timing uncertainty

DDSync reports how precisely each event is timed relative to a reference event at one station and phase. Events connected by many measurements can move together relative to that reference. The new Cholesky option calculates the contribution from one such shared timing shift directly. It estimates the remaining variance by averaging squared random timing changes.

## Turn it on

In the Julia configuration, edit the existing `[std]` section.

```toml
[std]
mode = "auto"                # Calculate small groups exactly.
auto_probe = "cholesky"      # Use squared draws for larger groups.
probes = 4096                # Maximum draws for one group.
rel_tol = 0.10               # Target sampling error as a fraction of each variance.
```

In MATLAB, set `cfg.std.mode = 'auto'` and `cfg.std.auto_probe = 'cholesky'`. The same sampling settings are under `cfg.std.hutch`, including `probes` and `rel_tol`. Setting `auto_probe` to `hutch` selects the original Hutchinson sampler. Existing configurations default to Hutchinson above the exact size limit.

## What 256 draws means

With these settings DDSync first checks all the events after 256 draws. If any still exceed the sampling target, it uses the same samples and adds more, checking at 512, 1024, 2048 and 4096. A lower maximum is respected, even between checkpoints. A minimum above 256 delays the first check. A maximum below the configured minimum leaves exact repair and fallback available.

For example, if a variance is estimated as 0.0001 seconds squared and its sampling standard error is 0.000008 seconds squared, the ratio is 0.08 and meets a 0.10 target. This tests numerical precision. Errors in measured arrivals or in the assumed timing model can still affect the exported uncertainty. The `min_sigma` setting separately controls the assumed scatter of the measured travel-time differences around the fit.

## Read the result

In `sync_metrics.txt`, look at the counts for **Converged Cholesky**, **Kept unconverged Cholesky** and **Fallback / proxy**. In `std_quality_*.txt`, the labels `cholesky` and `cholesky_unconverged` identify individual events. Selected exact repairs appear under `exact`. The existing 10% review advisory includes retained unconverged estimates from either sampler and fallback values. The numeric timing and uncertainty files use their existing columns.

## Which method fits the problem?

If many Hutch estimates miss the target in a station-phase group, try Cholesky before increasing the number of probes. In automatic mode, change `auto_probe` to `"cholesky"` so small groups still use exact variances. If you selected `mode = "hutch"` explicitly, change `mode` to `"cholesky"`. Keep the probe limit and `rel_tol` the same on the first comparison and check how many estimates meet the target in that group.

Small groups can use exact variances. Cholesky draws can save many probes when events share a strong timing shift. Hutchinson can be excellent when timing errors are nearly independent. In particular, with a diagonal system its random-sign estimate is exact. More probes help sampling precision but do not make a sparse factorization fit in memory.

Both samplers currently factor the same sparse matrix. The arrangement of connections determines how much memory that requires. Lowering `batch` reduces the number of random vectors held at once. `max_nred` limits the number of non-reference events per component that DDSync will sample. Components over the limit follow the configured fallback policy and show `reason=size_limit` in the report. A memory limit imposed by the operating system can stop the process before DDSync writes a report.

The [uncertainty guide](../uncertainty.md) gives the formulas and per-event decisions. The [Julia configuration template](../../DDSyncJulia/ddsync_config_template.toml) explains each setting alongside its value.
