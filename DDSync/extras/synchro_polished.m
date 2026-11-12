% SYNCHRO_POLISHED Compatibility entry point for the older DDSync script.
% The numerical implementation is shared with +ddsync/run.m so that pruning
% references and uncertainty failure handling receive the same bug fixes.
% For new workflows use run_ddsync.m or DDSyncJulia/run_ddsync.jl.
% Paths below are relative to MATLAB's current working directory.

addpath(fileparts(fileparts(mfilename('fullpath'))));
cfg = ddsync.config_default();

% Legacy input names. Change these to your own files.
cfg.io.infile_dt = 'dt_long.cc';
cfg.io.catalog_file = 'catalog.txt';
cfg.io.out_dt_sync = 'dt_sync.cc';
cfg.io.metrics_file = 'sync_metrics.txt';

% Common replacements for the old uppercase knobs:
% K_SIGMA / MIN_EDGES
% cfg.robust.k_sigma = 20;
% cfg.robust.min_edges = 30;
% IRLS_ITERS / C_HUBER
% cfg.robust.irls_iters = 10;
% cfg.robust.huber_c = 1.345;
% OUTPUT_WEIGHT_MODE
% cfg.output.dt_weight_mode = 'thetaStd';
% STD_PROBES is now a maximum budget; convergence is checked automatically.
% cfg.std.hutch.probes = 2048;
% cfg.std.hutch.min_probes = 128;
% cfg.std.hutch.rel_tol = 0.10;
% cfg.std.fallback = 'pseudo_degree';
% cfg.std.exact_max_nred = 2048;
% cfg.std.repair_max_nodes = 128; % selected full-system exact corrections
% cfg.std.unconverged = 'keep'; % retain positive estimates with relative MC SE < 1; otherwise fallback
% STD_MIN_DIAGREL has been retired; negative variances are never clipped.

out = ddsync.run(cfg);
