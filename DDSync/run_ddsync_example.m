% examples/run_ddsync_example.m
%
% Minimal example. Update file paths then run.

%appath() to add path to file if needed

cfg = ddsync.config_default();
cfg.io.infile_dt    = 'dt.cc';
cfg.io.catalog_file = 'catalog.txt';

cfg.io.out_dt_sync  = 'dt_sync.cc';
cfg.io.metrics_file = 'sync_metrics.txt';

% Optional: scratch space (fast local disk recommended for large catalogs)
% cfg.io.tmpdir = 'ddsync_tmp';

% Optional: reweighting + pruning
% cfg.robust.k_sigma     = 8;   % prune if |residual| > K_SIGMA * sigma_hat
% cfg.robust.min_edges   = 30;  % minimum kept edges per component
% cfg.robust.irls_iters  = 10;  % IRLS iterations (0 disables)
% cfg.robust.huber_c     = 1.345;

% Optional: base weights from CC
% cfg.weights.base_fun = 'cc2'; % use cc^2 instead of cc

% Optional: standard deviation / weights export
% cfg.std.mode     = 'pseudo_degree'; % faster than Hutch
% cfg.std.fallback = 'pseudo_degree';

% Optional: dt_sync weight mode
% cfg.output.dt_weight_mode = 'thetaStd'; % 'base'|'robust'|'combined'|'thetaStd'

ddsync.run(cfg);
