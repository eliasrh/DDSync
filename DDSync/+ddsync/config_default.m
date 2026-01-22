function cfg = config_default()
%DDSYNC.CONFIG_DEFAULT Default configuration for DDSync (standalone).
%
% Edit fields in a copy of this struct (recommended) rather than modifying
% this file. See examples/run_ddsync_example.m.
%
% All times are in seconds; distances are implicit via the dt network.

% ---------------- I/O ----------------
cfg.io.infile_dt    = 'dt.cc';
cfg.io.catalog_file = 'catalog.txt';
cfg.io.out_dt_sync  = 'dt_sync.cc';
cfg.io.metrics_file = 'sync_metrics.txt';

cfg.io.tmpdir       = 'ddsync_tmp';
cfg.io.thetadir     = 'theta';
cfg.io.thetastd_dir = 'thetastd';

cfg.io.prog_every_lines  = 2e6;
cfg.io.print_new_groups  = false;

% ---------------- Base weights ----------------
% Either:
%   cfg.weights.base_fun = 'cc'|'cc2'|'ones'
% or:
%   cfg.weights.base_fun_handle = @(cc) ... (function handle)
cfg.weights.base_fun = 'cc';
cfg.weights.base_fun_handle = [];

% ---------------- Robust pruning / IRLS ----------------
cfg.robust.k_sigma     = 20;     % prune if |residual| > K_SIGMA * sigma_hat
cfg.robust.min_edges   = 30;     % minimum kept edges in a station-phase component
cfg.robust.irls_iters  = 10;     % IRLS iterations (0 disables IRLS)
cfg.robust.huber_c     = 1.345;  % Huber tuning constant
cfg.robust.irls_rel_tol= 1e-3;   % stop IRLS if relative change small

% ---------------- Numerical ----------------
cfg.numeric.ridge_eps  = 1e-10;  % small Laplacian ridge for stability

% ---------------- Outputs ----------------
cfg.output.write_dt_sync          = true;
cfg.output.dt_weight_mode         = 'thetaStd'; % 'base'|'robust'|'combined'|'thetaStd'
cfg.output.dt_weight_decimals     = 4;

cfg.output.dt_decimals            = 5;   % decimals for dt values written to dt_sync
cfg.output.station_field_width    = 8;   % station field width for aligned dt_sync output
cfg.output.dt_field_width         = 10;  % dt numeric field width (includes sign)
cfg.output.weight_field_width     = 8;   % weight numeric field width
cfg.output.write_pruned_edges     = false; % if true: write pruned edges with weight=0 (diagnostic)
% thetaStd scaling to weights for writing dt_sync (only if dt_weight_mode='thetaStd')
cfg.output.theta_weight_scale_mode  = 'fixed';  % 'fixed' or 'median'
cfg.output.theta_weight_scale_fixed = 500;      % w=(1/std_dt)/scale
cfg.output.theta_weight_cap         = 1.0;      % cap written weights (Inf disables)

% ---------------- Theta std / weights export ----------------
cfg.std.export = true;

% mode:
%   'hutch'         : try Hutchinson diag(inv(L)); fallback if skipped
%   'pseudo_degree' : always pseudo std from node degree (fast)
%   'pseudo_weight' : always pseudo std from node weight (fast)
%   'none'          : do not write thetastd files (still can write dt_sync)
cfg.std.mode     = 'hutch';
cfg.std.fallback = 'pseudo_degree';  % when hutch skipped: 'pseudo_degree'|'pseudo_weight'|'nan'

% Pseudo-weight settings (used if std.mode or std.fallback uses pseudo_weight)
cfg.std.pseudo_weight_source = 'combined';  % 'combined'|'robust'|'base'
cfg.std.pseudo_weight_eps    = 1e-6;        % floor on node weight before inversion

% Write extra columns in std_theta_*.txt
% Default file has 4 columns:
%   EventID  std_theta  refEventID  degree
% If write_weight_column=true, add 5th column:
%   thetaWeight = min(cap, (1/std_theta)/theta_weight_scale_fixed)
% If write_alt_weight_column=true, add 6th column:
%   nodeWeight = aggregated node weight used by pseudo_weight (NaN otherwise)
cfg.std.write_weight_column     = true;
cfg.std.write_alt_weight_column = false;

% Hutchinson settings (only used if std.mode='hutch')
cfg.std.hutch.probes            = 100000;    % 0 disables Hutch (forces fallback if export=true)
cfg.std.hutch.probe_dist        = 'rademacher'; % 'rademacher'|'gaussian'
cfg.std.hutch.batch             = 250;
cfg.std.hutch.report_every_batch= 10;

cfg.std.hutch.max_nred          = 40000;     % skip Hutch if reduced system larger than this
cfg.std.hutch.min_sigma         = 5e-4;      % floor on sigma_hat (seconds)
cfg.std.hutch.min_diagrel       = 1e-12;     % floor on diag(inv(L)) relative to median(diag)

end
