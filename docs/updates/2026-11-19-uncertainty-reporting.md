# Uncertainty reporting and probe defaults — 19 November 2026

The per-event estimator introduced on 12 November could still assign many values to degree fallback when their relative Monte Carlo standard error was at least 1. The detailed records exposed those choices, but users had to count them manually, and component records interrupted the observation table. This update makes the proportions visible and raises the default sampling budget in Julia and MATLAB.

## Changes

- The default maximum probe count is **4096**, doubled from 2048. Sampling still stops early when every entry converges. The full exact cap remains 2048 unpinned events and the selected exact repair budget remains 128.
- The metrics report and final terminal summary show counts and percentages for **exact**, **converged Hutch**, **kept unconverged Hutch**, and **fallback/proxy** values. Converged runs are included, not just failures.
- The report gives global and station-phase summaries, using event counts rather than averaging percentages across differently sized components.
- Observation-table headers and values use aligned columns. Component records are collected at the end in readable multiline records, including the final method, original sampling outcome, repair counts and both residual-noise scales.
- If **kept unconverged + fallback is at least 10%**, the report advises more probes and/or exact repairs, or a larger full exact cap when practical. Groups are assessed individually too, so a poor small group is not hidden by the global percentage.

This update changes reporting and the default budget. It does not remove the relative-error guard, change graph references, change the exact solver, or change the existing numeric theta/std file columns. Graph bounds remain opt-in.

## Reading the percentages

The denominator is the number of supported **non-reference std entries**. An event appearing in multiple station-phase files contributes once per file. Gauge references, rejected or isolated rows and deliberately disabled uncertainties are excluded; disabled entries are counted separately. With no assessed entries, percentages are `N/A`.

“Exact” includes full exact calculations and successful selected repairs; repaired entries must not be counted twice. “Fallback/proxy” includes degree and weight proxies and explicitly selected graph bounds, even if a proxy was requested directly rather than triggered by a failure.

The combined unconverged-plus-fallback share is labeled as requiring review. It is a measure of computational quality, not an estimate of the fraction of physically wrong uncertainties. Exact plus converged Hutch is the favorable share; a large value there does not itself require more probes. Even those calculations depend on the residual-noise model and input weights.

The advisory threshold of 10% is an operational rule, **separate from** the per-entry `rel_tol=0.10` convergence criterion. It does not establish a confidence level. Doubling probes reduces stochastic sampling noise but cannot guarantee convergence on every graph. If the component reason is a size limit, disabled sampling, invalid noise scale or numerical failure, address that cause rather than merely increasing the probe count.

The retained-estimate guard is unchanged: a positive but unrepaired estimate is retained only when it has enough samples and relative MC standard error below 1. Otherwise it uses the configured fallback. That replacement contributes to the visible fallback count and review share.

## What `sigma_hat` means

For each retained pair, DDSync compares the measured differential travel time with the fitted difference `theta_i - theta_j`. The mismatch is the residual. The scatter of these mismatches estimates noise in the pairwise measurements, separately for each connected group at one station and phase. All the following scales are in seconds:

| Name | Meaning |
| --- | --- |
| Component `sigma_hat` | Measured scatter of the travel-time mismatches after the final fit in that component. |
| Component `sigma_std` | Mismatch scatter used in the uncertainty formula, including the minimum `min_sigma` when enabled. |
| Group-table `sigma_hat` | Median of raw final component scales within the station-phase group, for reporting only. |
| `std_theta` | The resulting uncertainty for an individual theta value, after graph or degree scaling. |

With the positive noise floor enabled, `sigma_std = max(sigma_hat, min_sigma)`; with the floor disabled, `sigma_std = sigma_hat`. The `sigma` in a description such as `sigma/sqrt(degree)` means `sigma_std`, not a different fitted noise estimate. `pseudo_weight` uses neither residual scale. The Monte Carlo relative standard error measures stochastic diagonal-estimation error and is a separate quantity.

## Configuration and compatibility

For a Julia TOML configuration that explicitly set the old budget:

```toml
[std]
probes = 4096
```

MATLAB equivalent:

```matlab
cfg.std.hutch.probes = 4096;
```

Explicit budgets in existing configurations remain unchanged until edited. Removing the explicit budget also adopts the new default. Existing 2048-probe synthetic studies and their recorded results retain that budget; they have not been relabeled as 4096-probe experiments.

For code that reads `sync_metrics.txt`, parse the observation table on whitespace rather than requiring tabs. A `# std_component` record now spans multiple comment lines; do not assume every field occurs on its opening line. The per-event `std_quality_*.txt` files retain their existing format and are still the way to identify particular fallback or unconverged entries.

See the [uncertainty guide](../uncertainty.md), [Julia manual](../../DDSyncJulia/README.md) and [MATLAB manual](../../DDSync/README.md) for the reporting layout, scale definitions and budget choices.

## Validation

Julia regression tests passed **1,370 assertions**, including 28 reporting checks. Native MATLAB regression tests passed **629 checks**.

The reporting tests cover unequal component sizes, exact repairs counted once, all-converged output, empty and disabled uncertainty, the 10% advisory boundary, a poor group hidden by good overall coverage, aligned headers including long station names, and component records confined to the final section. An end-to-end mixed-method run checks reported percentages against the actual per-event quality files: 8 exact repairs and 28 degree proxies yield 22.22% exact and 77.78% requiring review across 36 non-reference entries.

The historical November 12 study scripts explicitly retain their 2048-probe budgets so their published experiments remain reproducible. No new accuracy claim is inferred from the doubled default alone.

## Configuration-comment correction

The template previously described `apply_min_sigma` as enabling or disabling a “std_theta floor.” That wording was incorrect: both `min_sigma` and `apply_min_sigma` control the minimum assumed scatter of measured differential travel times around the fitted differences `theta_i - theta_j` when calculating uncertainty. They do not guarantee any minimum exported theta standard deviation. For example, if the measured scatter is `1e-5` seconds and the enabled minimum is `5e-4` seconds, the degree proxy at degree 100 is `5e-5` seconds.

The Julia TOML, source-default comments, MATLAB configuration examples, metrics legend and current manuals now define these mismatches explicitly. Parameter names, values and numerical behavior are unchanged. Verification compared parsed TOML values and executable source text before and after the configuration-comment changes; the metrics changes affect only explanatory text. After this correction, all 1,370 Julia regression assertions and 23 targeted MATLAB reporting checks passed.
