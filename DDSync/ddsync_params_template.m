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

% Example: disable expensive Hutch std and use pseudo-degree fallback
% cfg.std.mode = 'pseudo_degree';
% cfg.std.fallback = 'pseudo_degree';

% Example: choose output dt weights
% cfg.output.dt_weight_mode = 'thetaStd'; % 'base'|'robust'|'combined'|'thetaStd'

end
