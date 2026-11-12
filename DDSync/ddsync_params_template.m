function cfg = ddsync_params()
%DDSYNC_PARAMS User-editable parameter file for DDSync.
%
% Copy this file from ddsync_params_template.m to ddsync_params.m and edit.

cfg = ddsync.config_default();

% ---- REQUIRED ----
cfg.io.infile_dt    = 'dt.cc';
cfg.io.catalog_file = 'catalog.txt';

% ---- Outputs ----
cfg.io.out_dt_sync  = 'dt_sync.cc';
cfg.io.metrics_file = 'sync_metrics.txt';

% Optional: scratch space (fast local disk recommended for large catalogs)
% cfg.io.tmpdir = 'ddsync_tmp';

% Example: use CC^2 base weights
% cfg.weights.base_fun = 'cc2';

% Example: tune pruning/robust reweighting
% cfg.robust.k_sigma    = 8;
% cfg.robust.min_edges  = 30;
% cfg.robust.irls_iters = 10;
% cfg.robust.min_scale  = 5e-4; % robust scale floor for pruning/IRLS only

% Optional: a fast weighting proxy (not a statistical uncertainty)
% cfg.std.mode = 'pseudo_degree';

% Default uncertainty policy (customize these only when needed):
% cfg.std.mode = 'auto'; % exact for small components, otherwise per-event Hutch/repair
% cfg.std.fallback = 'pseudo_degree'; % only unresolved entries; graph_bound is opt-in
% cfg.std.exact_max_nred = 2048;
% cfg.std.repair_max_nodes = 128; % selected full-system exact corrections
% cfg.std.unconverged = 'keep'; % retain positive estimates with relative MC SE < 1; otherwise fallback
% cfg.std.min_sigma = 5e-4;       % residual-noise floor for exported std_theta
% cfg.std.apply_min_sigma = true; % set false (or min_sigma=[]/0) to disable

% Example: choose output dt weights
% cfg.output.dt_weight_mode = 'thetaStd'; % 'base'|'robust'|'combined'|'thetaStd'

% Example: preserve microsecond/sub-microsecond values in text outputs
% cfg.output.theta_decimals = 12;
% cfg.output.thetastd_decimals = 12;
% Adaptive Hutchinson: maximum budget, early convergence checks, safe fallback.
% cfg.std.hutch.probes = 2048;
% cfg.std.hutch.min_probes = 128;
% cfg.std.hutch.rel_tol = 0.10;
% Read # std_component diagnostics in sync_metrics.txt before using uncertainties.
% See ../docs/uncertainty.md for probe selection and fallback behavior.


end
