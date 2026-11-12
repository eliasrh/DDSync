# DDSync maintenance update — 18 October 2026

**Historical report:** the [12 November update](2026-11-12-per-event-uncertainty.md) supersedes the uncertainty defaults and whole-component fallback policy described here. The post-pruning reference fix remains.

Julia remains the primary implementation. This update fixes the same numerical bugs in the supported MATLAB implementation and updates both manuals and configuration examples.

## Changes

**References follow the surviving graph.** DDSync now rebuilds connected components after pruning, repeats the minimum-edge and cycle-redundancy checks, selects the highest-degree surviving reference in each accepted component, and refits each component independently. Ties use the smallest event ID. Theta, reference IDs, residual scales, degrees and uncertainties all describe that final component. Unsupported fragments and isolated events retain missing theta/reference metadata and are not exported as usable constraints.

**Failed uncertainty estimates no longer become false precision.** Hutchinson estimates are checked for positive finite diagonals, solve accuracy and per-node Monte Carlo standard error. The code no longer replaces negative estimates with a tiny positive floor. Final IRLS fitting also uses exactly the weights passed to the uncertainty calculation and output.

**Fallback keeps accepted theta constraints usable in GraphSplit.** The default first tries a bounded exact inverse-diagonal calculation. If it is too large or fails, a deterministic spanning-tree resistance bound supplies finite, conservative standard deviations for non-reference nodes. It preserves the component's reference star rather than introducing missing sigma or discarding constraints. GraphSplit's current reader can replace missing sigma with its minimum value, so exporting NaN for a failed estimate would not be an appropriate failure policy. The gauge reference itself still has zero theta and zero std.

The graph bound is conditional on DDSync's residual-noise and weighting model. It can be substantially larger than the exact diagonal and is labeled `graph_bound` in the metrics. It is not a calibrated bound on all physical/model errors. An invalid noise scale, disconnected numerical component or overflow that prevents a finite bound stops the run with an explanation. Explicitly requesting a theta-only mode remains possible.

## Defaults and compatibility

| Setting | Previous | New |
| --- | --- | --- |
| Maximum random probes | 100000 fixed | 2048 maximum, with early stopping |
| Minimum samples before convergence | No check | 128 |
| Probe batch | 250 | 128 |
| Per-node relative MC standard error | No check | At most 0.10 for every unpinned node |
| Default fallback | Degree proxy when skipped | Exact up to 2048 unpinned nodes, then finite graph bound |
| Negative diagonal handling | Tiny positive floor | Continue sampling within budget, then fallback if still invalid |

`probes=128` may be enough on favorable graphs; it is not an accuracy guarantee. The check is an empirical standard-error diagnostic, not a strict relative-error bound or a simultaneous confidence interval. Existing configurations keep explicit probe budgets and proxy choices. Legacy `fallback="nan"` / `"none"` now select the finite graph bound. The old `min_diag_rel` / MATLAB `min_diagrel` field is accepted but ignored.

Per-component metrics now record the reference, node count, actual method, fallback reason, probes used and worst relative standard error. The top-level metrics also record convergence and exact-fallback settings. `pseudo_degree` and `pseudo_weight` remain explicit practical-weight options and are documented as neither calibrated uncertainty nor conservative bounds.

Read the [uncertainty and probe guide](../uncertainty.md) for configuration examples, the graph-bound formula, diagnostics, GraphSplit handling and scale-floor behavior.

## Validation

Tests were run locally with **Julia 1.12.3** and **native MATLAB R2025b Update 4** on macOS. MATLAB was not replaced by an Octave approximation.

| Check | Result |
| --- | --- |
| Julia regression suite | 1242 assertions passed: event IDs, graph splitting, covariance, convergence, fallback and exports. |
| MATLAB regression suite | 531 checks passed on the same synthetic observations and independent matrix oracles. |
| GraphSplit integration | 344 assertions passed using GraphSplit 0.4.0's actual theta reader and Stage-1 star builder, including MATLAB-generated fallback exports. |
| Deterministic Julia/MATLAB synthetic comparison | Theta/std values agreed within `1e-10` seconds; synchronized differential-time files were byte-identical. |
| Full bundled Spanish Springs example | Both retained 608205 of 608577 observations across 89 station-phase groups; 70 accepted components. |
| Bundled-example cross-language comparison | All 14086 finite theta rows had matching reference IDs; exported theta and differential-time values matched exactly at the requested output precision. All accepted non-reference std values were finite and positive. |

The split fixture contains two nearly complete ten-event components joined by contradictory links. Pruning removes the links and isolates an extra event. The references become 3 and 13 within the surviving components, and the isolated, zero-weight-only and absent events cannot become theta-star members. Separate tests reject a post-pruning tree fragment, exercise both dense/dictionary Julia storage and IRLS on/off, and compare exact uncertainty to an independently assembled incidence-matrix inverse.

Failure tests cover nonpositive/nonfinite estimates, insufficient samples, nonconvergence, singular exact systems, invalid noise scales, size-limited exact fallback, deliberately disabled uncertainty, optional diagnostic columns and mixed-method median scaling. The finite graph bound is checked against exact covariance and an analytically solvable weighted tree.

The downstream test deliberately sets one probe and disables exact fallback by size, then reads the resulting files with GraphSplit. All 36 expected stars remain finite, every star stays inside its own component, and the weighted incidence matrix has rank 18 for 20 supported events: exactly the two intended component gauge freedoms. Removing one reference per component gives a positive-definite normal matrix. This tests the uncertainty/graph interface, not whether an arbitrary station geometry fully constrains a 3-D relocation.

## Probe study

The reproducible study uses five graph families, four maximum budgets and ten seeds per language. The table shows the number of runs passing the all-node 10% relative-diagonal-SE check. Counts are Julia / MATLAB, each out of ten; failed stochastic runs use fallback in normal DDSync processing.

| Synthetic graph | 128 | 512 | 2048 | 8192 | Typical probes at the 8192 budget, Julia / MATLAB |
| --- | ---: | ---: | ---: | ---: | ---: |
| Dense, 20 nodes | 0 / 0 | 0 / 0 | 10 / 10 | 10 / 10 | 640 / 640 |
| Dense, 100 nodes | 0 / 0 | 0 / 0 | 0 / 0 | 10 / 10 | 3968 / 3840 |
| Ring, 100 nodes | 0 / 0 | 0 / 0 | 0 / 0 | 10 / 10 | 5184 / 5120 |
| Two dense halves with a weak link, 100 nodes | 0 / 0 | 0 / 0 | 0 / 0 | 9 / 8 | 7360 / 7360 |
| Strong hub plus neighbor links, 100 nodes | 10 / 10 | 10 / 10 | 10 / 10 | 10 / 10 | 128 / 128 |

The largest relative std error among accepted runs in these studies was about 14.2%, consistent with a standard-error diagnostic rather than a guaranteed error ceiling. More probes alone are not a universal solution: graph structure matters, and small/moderate components can be handled deterministically more cheaply than a large blanket probe count. The bundled example used exact fallback in 45 Julia components and 44 MATLAB components; the remaining components passed the stochastic check. Different random streams explain the method-count difference.

These synthetic results motivate a modest initial budget with checked fallback. For a large component that cannot use exact calculation, raise the budget if a tighter estimate is needed instead of a conservative graph bound. Compare representative results across seeds or against exact calculations whenever feasible.

## Reproducing the checks

From the repository root:

```bash
julia --startup-file=no DDSyncJulia/test/runtests.jl
julia --startup-file=no DDSyncJulia/test/probe_study.jl
julia --startup-file=no DDSyncJulia/test/graphsplit_integration.jl /path/to/GraphSplit.jl
```

In MATLAB, from the repository root:

```matlab
addpath('DDSync/test');
run_tests;
probe_study('probe_study.csv');
```

To include MATLAB exports in the GraphSplit check, retain them with `run_tests('/absolute/path/to/matlab-results')` and pass that directory as the integration script's second argument. GraphSplit is an optional test dependency, not a dependency of DDSync itself. The repository's Julia CI remains on Julia 1.10 and current stable; a native MATLAB CI job has been added.

## Documentation and entry-point review

The main README now identifies Julia as the primary implementation, links this report, and gives a runnable example using the bundled data. Both manuals and configuration comments explain component references, convergence, fallback, relative paths and single-precision storage versus decimal formatting. The older MATLAB `synchro_polished.m` now calls the supported package instead of retaining duplicate buggy solver code. Its old uppercase settings should be migrated to the documented `cfg` fields.

Additional cleanup removes MATLAB's `mad` toolbox dependency using the equivalent median formula, fixes the optional node-weight column when the theta-weight column is disabled, clarifies the preprocessing helpers' actual output filenames, labels the checked-in historical metrics, and removes a placeholder citation DOI. Numerical fixes do not change file column layouts or Julia's event-ID mapping convention.

Rerun previously affected datasets. Existing theta files alone do not reconstruct the surviving observation graph, and old clipped std files do not preserve the failed stochastic estimates needed for a reliable repair.
