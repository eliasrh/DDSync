% examples/run_ddsync_example.m
%
% Minimal example. Update file paths then run.

%appath() to add path to file if needed

cfg = ddsync.config_default();
cfg.io.infile_dt    = 'dt.cc';
cfg.io.catalog_file = 'catalog.txt';

cfg.io.out_dt_sync  = 'dt_sync.cc';
cfg.io.metrics_file = 'sync_metrics.txt';

cfg.robust.k_sigma     = 8;     % prune if |residual| > K_SIGMA * sigma_hat

% Faster: pseudo std instead of Hutch
% cfg.std.mode = 'pseudo_degree';

ddsync.run(cfg);
