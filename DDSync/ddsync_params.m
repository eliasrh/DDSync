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

% Example: use CC^2 base weights
% cfg.weights.base_fun = 'cc2';

% Example: disable expensive Hutch std and use pseudo-degree fallback
% cfg.std.mode = 'pseudo_degree';

end
