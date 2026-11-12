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

% Example: preserve microsecond/sub-microsecond values in text outputs
% cfg.output.theta_decimals = 12;
% cfg.output.thetastd_decimals = 12;

% Example: use CC^2 base weights
% cfg.weights.base_fun = 'cc2';

% Optional: a fast weighting proxy (not a statistical uncertainty)
% cfg.std.mode = 'pseudo_degree';
% Defaults: exact for small components, otherwise per-event Hutch and repair.
% cfg.std.mode = 'auto';
% cfg.std.exact_max_nred = 2048;
% cfg.std.repair_max_nodes = 128;
% cfg.std.unconverged = 'keep'; % retain positive unresolved estimates with relative MC SE < 1
% cfg.std.fallback = 'pseudo_degree'; % practical proxy, not calibrated uncertainty
% cfg.std.hutch.probes = 2048;
% cfg.std.hutch.min_probes = 128;
% cfg.std.hutch.rel_tol = 0.10;
% Read std_quality_*.txt and # std_component diagnostics before using uncertainties.
% See ../docs/uncertainty.md for probe selection and fallback behavior.


end
