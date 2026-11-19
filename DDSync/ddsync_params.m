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
% Defaults: exact for small components, otherwise the selected sampler and repair.
% cfg.std.auto_probe = 'cholesky'; % Optional. Use squared draws on larger groups. Default is 'hutch'.
% cfg.std.mode = 'auto';
% cfg.std.exact_max_nred = 2048;
% cfg.std.repair_max_nodes = 128;
% cfg.std.unconverged = 'keep'; % retain positive unresolved estimates with relative MC SE < 1
% cfg.std.fallback = 'pseudo_degree'; % practical proxy, not calibrated uncertainty
% Both samplers use these settings. Cholesky checks after 256 draws, then doubles
% up to the maximum if any event still exceeds the sampling target.
% cfg.std.hutch.probes = 4096; % Maximum draws per connected group.
% cfg.std.hutch.min_probes = 128;
% cfg.std.hutch.rel_tol = 0.10;
% Read the method percentages and review advice in sync_metrics.txt.
% Unrepaired relative MC SE >= 1 uses fallback and contributes to the review share.
% Component details are at the end; std_quality_*.txt identifies individual rows.
% See ../docs/uncertainty.md for probe selection and fallback behavior.


end
