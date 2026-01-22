function out = run(cfg)
%DDSYNC.RUN  Synchronize differential travel times (dt.cc) per station-phase.
%
% This is a functional, package-friendly wrapper around the authoritative
% "synchro_polished.m" implementation, redesigned for standalone use.
%
% Usage:
%   cfg = ddsync.config_default();
%   cfg.io.infile_dt    = 'dt.cc';
%   cfg.io.catalog_file = 'catalog.txt';
%   out = ddsync.run(cfg);
%
% See also: ddsync.config_default

if nargin < 1 || isempty(cfg)
    cfg = ddsync.config_default();
else
    cfg = ddsync.merge_struct(ddsync.config_default(), cfg);
end

% ---------- unpack config ----------
infile_dt    = cfg.io.infile_dt;
catalog_file = cfg.io.catalog_file;
out_dt_sync  = cfg.io.out_dt_sync;
metrics_file = cfg.io.metrics_file;

TMPDIR       = cfg.io.tmpdir;
THETADIR     = cfg.io.thetadir;
THETASTD_DIR = cfg.io.thetastd_dir;

PROG_EVERY   = cfg.io.prog_every_lines;

if ~exist(TMPDIR,'dir'), mkdir(TMPDIR); end
if ~exist(THETADIR,'dir'), mkdir(THETADIR); end
if cfg.std.export && ~exist(THETASTD_DIR,'dir'), mkdir(THETASTD_DIR); end

% Base weight from CC (or uniform). Keep nonnegative.
weight_fun_base = cfg.weights.base_fun_handle;
if isempty(weight_fun_base)
    weight_fun_base = ddsync.base_weight_fun(cfg.weights.base_fun);
end

% Pruning / robust
K_SIGMA      = cfg.robust.k_sigma;
MIN_EDGES    = cfg.robust.min_edges;
IRLS_ITERS   = cfg.robust.irls_iters;
C_HUBER      = cfg.robust.huber_c;
IRLS_REL_TOL = cfg.robust.irls_rel_tol;

% Numerical stabilization
RIDGE_EPS    = cfg.numeric.ridge_eps;

% ---- Weight written into dt_sync.cc (third column) ----
OUTPUT_WEIGHT_MODE = cfg.output.dt_weight_mode;   % 'base'|'robust'|'combined'|'thetaStd'

% Theta std estimation / export
EXPORT_THETA_STD = cfg.std.export;

% Hutchinson diag(inv(L)) estimator settings
STD_PROBES      = cfg.std.hutch.probes;
STD_PROBE_DIST  = cfg.std.hutch.probe_dist;
STD_BATCH       = cfg.std.hutch.batch;
STD_REPORT_EVERY_BATCH = cfg.std.hutch.report_every_batch;
STD_MAX_NRED    = cfg.std.hutch.max_nred;
STD_MIN_SIGMA   = cfg.std.hutch.min_sigma;
STD_MIN_DIAGREL = cfg.std.hutch.min_diagrel;

% Pseudo std / fallback behavior
THETASTD_MODE            = cfg.std.mode;      % 'hutch'|'pseudo_degree'|'pseudo_weight'|'none'
THETASTD_FALLBACK        = cfg.std.fallback;  % when hutch skipped: 'pseudo_degree'|'pseudo_weight'|'nan'
PSEUDO_WEIGHT_SOURCE     = cfg.std.pseudo_weight_source; % 'combined'|'robust'|'base'
PSEUDO_WEIGHT_EPS        = cfg.std.pseudo_weight_eps;

% thetastd output extras
WRITE_THETASTD_WEIGHTCOL = cfg.std.write_weight_column; % adds 5th column: thetaWeight
WRITE_ALT_NODE_WEIGHTCOL = cfg.std.write_alt_weight_column; % adds 6th column: nodeWeight (if available)

% thetaStd scaling for dt weights (writing only)
THETASTD_SCALE_MODE  = cfg.output.theta_weight_scale_mode; % 'fixed' or 'median'
THETASTD_SCALE_FIXED = cfg.output.theta_weight_scale_fixed;
THETASTD_WEIGHT_CAP  = cfg.output.theta_weight_cap;

% Streaming / progress
PROG_EVERY = cfg.io.prog_every_lines;

% ---------- PRE-WARN on probes ----------
if EXPORT_THETA_STD && STD_PROBES > 0 && strcmpi(THETASTD_MODE,'hutch')
    fprintf('\n[STD WARNING]\n');
    fprintf('  Hutchinson diag(inv(L)) convergence is slow: estimator error ~ O(1/sqrt(K)).\n');
    fprintf('  If you want ms-level stability in std(theta), K may need to be VERY large.\n');
    fprintf('  For relocation weights, relative scaling often matters more than absolute calibration.\n\n');
end

% ---------- STEP 0: read catalog max ID ----------
fprintf('Reading %s to get max EventID...\n', catalog_file);
maxEventID = readMaxEventID(catalog_file);
fprintf('Max EventID = %d\n', maxEventID);

% ---------- STEP 1: PASS 1 spool dt by (STA,PH) ----------
fprintf('Pass 1: streaming %s and spooling groups...\n', infile_dt);
[groupKeys, groupFiles, groupCounts] = spoolGroups(infile_dt, TMPDIR, PROG_EVERY, cfg.io.print_new_groups);
nGroups = numel(groupKeys);
fprintf('Found %d station+phase groups.\n', nGroups);

key2g = containers.Map('KeyType','char','ValueType','double');
for g=1:nGroups, key2g(groupKeys{g}) = g; end

% ---------- STEP 2: process groups ----------
theta_fulls = cell(nGroups,1);
std_fulls   = cell(nGroups,1);
nodew_fulls = cell(nGroups,1); % optional: per-node alternative weights
deg_fulls   = cell(nGroups,1);
decFiles    = cell(nGroups,1);

overall_edges = 0;
overall_kept  = 0;

% pseudo stats
nPseudoGroups = 0;
nPseudoComps  = 0;

% ---------- metrics file (rich per-group table) ----------
fidm = fopen(metrics_file,'w'); assert(fidm>0);
fprintf(fidm,'# ddsync metrics\n');
fprintf(fidm,'# K_SIGMA=%g MIN_EDGES=%d IRLS_ITERS=%d C_HUBER=%g MODE=%s STD_EXPORT=%d STD_MODE=%s STD_PROBES=%d\n', ...
    K_SIGMA, MIN_EDGES, IRLS_ITERS, C_HUBER, OUTPUT_WEIGHT_MODE, EXPORT_THETA_STD, THETASTD_MODE, STD_PROBES);
fprintf(fidm,'group\tn_edges\tn_kept\tpruned_pct\trobust_p50\trobust_p05\trobust_p95\tfrac_robust_lt_0p99\tsigma_hat\n');

for g=1:nGroups
    key = groupKeys{g};
    fn  = groupFiles{g};
    n_edges = groupCounts(g);
    overall_edges = overall_edges + n_edges;

    fprintf('\n[%d/%d] Processing group %s (edges=%d)\n', g, nGroups, key, n_edges);

    fid = fopen(fn,'r'); assert(fid>0);
    raw = fread(fid, [4, Inf], 'single'); fclose(fid);
    raw = raw.'; % rows: [i j dt cc]
    gi = double(raw(:,1));
    gj = double(raw(:,2));
    gd = double(raw(:,3));
    gc = double(raw(:,4));

    w_base = max(weight_fun_base(gc), 0);

    [theta_full, ref_full, keep_mask, wrob, sigma_hat, std_full, deg_full, nodew_full, pseudo_info] = ...
        processStationPhaseGroup(gi, gj, gd, w_base, maxEventID, ...
            K_SIGMA, MIN_EDGES, IRLS_ITERS, C_HUBER, IRLS_REL_TOL, RIDGE_EPS, ...
            EXPORT_THETA_STD, THETASTD_MODE, THETASTD_FALLBACK, ...
            STD_PROBES, STD_PROBE_DIST, STD_BATCH, STD_REPORT_EVERY_BATCH, ...
            STD_MAX_NRED, STD_MIN_SIGMA, STD_MIN_DIAGREL, ...
            PSEUDO_WEIGHT_SOURCE, PSEUDO_WEIGHT_EPS, ...
            THETASTD_SCALE_FIXED);

    if pseudo_info.used_pseudo
        nPseudoGroups = nPseudoGroups + 1;
        nPseudoComps  = nPseudoComps + pseudo_info.n_pseudo_comps;
    end

    n_kept = sum(keep_mask);
    overall_kept = overall_kept + n_kept;

    theta_fulls{g} = single(theta_full);
    std_fulls{g}   = single(std_full);
    nodew_fulls{g} = single(nodew_full);
    deg_fulls{g}   = uint32(deg_full);

    [sta, ph] = splitKey(key);

    % write theta
    theta_fn = fullfile(THETADIR, sprintf('theta_%s_%s.txt', sta, ph));
    writeThetaTriples(theta_fn, theta_full, ref_full);

    % write std + degree (+ optional weight columns)
    if EXPORT_THETA_STD
        std_fn = fullfile(THETASTD_DIR, sprintf('std_theta_%s_%s.txt', sta, ph));
        writeThetaStdWithDegree(std_fn, std_full, ref_full, deg_full, ...
            THETASTD_SCALE_FIXED, THETASTD_WEIGHT_CAP, ...
            WRITE_THETASTD_WEIGHTCOL, WRITE_ALT_NODE_WEIGHTCOL, nodew_full);
    end

    % decisions file (in spooled-edge order)
    dec_fn = fullfile(TMPDIR, sprintf('dec_%s.bin', sanitizeKey(key)));
    writeDecisionFile(dec_fn, keep_mask, wrob);
    decFiles{g} = dec_fn;

    % metrics (per group)
    pruned_pct = 100*(n_edges-n_kept)/max(n_edges,1);
    if n_kept>0
        rob_kept = wrob(keep_mask);
        r50 = median(rob_kept);
        r05 = prctile(rob_kept,5);
        r95 = prctile(rob_kept,95);
        frac_dw = mean(rob_kept < 0.99);
    else
        r50=NaN; r05=NaN; r95=NaN; frac_dw=NaN;
    end

    % write to metrics file (tab-delimited)
    fprintf(fidm,'%s\t%d\t%d\t%.3f\t%.4f\t%.4f\t%.4f\t%.3f\t%.6g\n', ...
        key, n_edges, n_kept, pruned_pct, r50, r05, r95, frac_dw, sigma_hat);

    fprintf('  kept=%d (%.2f%% pruned), sigma_hat=%.4g s, robust_w p05/p50/p95=%.3g/%.3g/%.3g\n', ...
        n_kept, pruned_pct, sigma_hat, r05, r50, r95);
end

% ---------- overall summary ----------
fprintf(fidm,'\n# OVERALL\n');
fprintf(fidm,'overall_edges=%d\n', overall_edges);
fprintf(fidm,'overall_kept=%d\n', overall_kept);
fprintf(fidm,'overall_pruned_pct=%.3f\n', 100*(overall_edges-overall_kept)/max(overall_edges,1));
fprintf(fidm,'pseudo_groups_used=%d\n', nPseudoGroups);
fprintf(fidm,'pseudo_components_used=%d\n', nPseudoComps);
fclose(fidm);

if nPseudoGroups>0
    fprintf('\n[THETASTD] Pseudo std/weights were used in %d/%d groups (components=%d). See %s.\n', ...
        nPseudoGroups, nGroups, nPseudoComps, metrics_file);
end

% ---------- STEP 3: optional median scaling for thetaStd weights ----------
thetaStdScale = THETASTD_SCALE_FIXED;
if strcmpi(OUTPUT_WEIGHT_MODE,'thetaStd') && strcmpi(THETASTD_SCALE_MODE,'median')
    fprintf('Computing global median scale for thetaStd weights (extra streaming pass)...\n');
    thetaStdScale = computeThetaStdMedianScale(infile_dt, key2g, decFiles, std_fulls);
    fprintf('thetaStdScale (median of 1/std_dt) = %.6g\n', thetaStdScale);
    if thetaStdScale<=0 || ~isfinite(thetaStdScale), thetaStdScale = THETASTD_SCALE_FIXED; end
end

% ---------- STEP 4: PASS 2 write dt_sync ----------
if cfg.output.write_dt_sync
    fprintf('Pass 2: writing %s ...\n', out_dt_sync);
    decReaders = openDecisionReaders(decFiles);

    fin  = fopen(infile_dt,'r');  assert(fin>0);
    fout = fopen(out_dt_sync,'w'); assert(fout>0);

    cur_i = NaN; cur_j = NaN;
    nLine = 0;
    STA_W      = cfg.output.station_field_width;
    DT_W       = cfg.output.dt_field_width;
    W_W        = cfg.output.weight_field_width;
    DT_DECIMALS= cfg.output.dt_decimals;
    W_DECIMALS = cfg.output.dt_weight_decimals;
    lineFmt = sprintf('    %%-%ds %%%d.%df %%%d.%df %%s\n', STA_W, DT_W, DT_DECIMALS, W_W, W_DECIMALS);
    write_pruned = isfield(cfg.output,'write_pruned_edges') && cfg.output.write_pruned_edges;

    pendingHeader = '';
    pendingHasHeader = false;
    pendingPrinted = false;

    while true
        L = fgetl(fin); if ~ischar(L), break; end
        nLine = nLine + 1;
        if mod(nLine, PROG_EVERY)==0
            fprintf('  pass2 lines=%d\n', nLine);
        end

        Lt = strtrim(L);
        if isempty(Lt), continue; end

        if L(1)=='#'
            pendingHeader = L;
            pendingHasHeader = true;
            pendingPrinted = false;
            nums = sscanf(L(2:end),'%d %d %f');
            if numel(nums)>=2, cur_i=nums(1); cur_j=nums(2); else, cur_i=NaN; cur_j=NaN; end
            continue;
        end

        toks = textscan(L,'%s %f %f %s');
        if isempty(toks{1}), continue; end
        sta = toks{1}{1};
        dt_in  = toks{2}(1);
        cc_in  = toks{3}(1);
        ph  = toks{4}{1};
        key = [sta '_' ph];

        if ~isKey(key2g, key)
            % Unknown group: pass through, but keep formatting consistent
            if pendingHasHeader && ~pendingPrinted
                fprintf(fout,'%s\n', pendingHeader);
                pendingPrinted = true;
            end
            fprintf(fout, lineFmt, sta, dt_in, cc_in, ph);
            continue;
        end
        g = key2g(key);

        % read decision + robust weight for this edge (spooled order)
        [keep, wrob, decReaders(g)] = readNextDecision(decReaders(g));
        if ~keep
            % pruned edge: either skip (default) or write with weight=0 (diagnostic)
            if write_pruned
                if pendingHasHeader && ~pendingPrinted
                    fprintf(fout,'%s\n', pendingHeader);
                    pendingPrinted = true;
                end
                fprintf(fout, lineFmt, sta, dt_in, 0.0, ph);
            end
            continue;
        end

        % theta correction
        tf = theta_fulls{g};
        dt_corr = dt_in;
        if ~isnan(cur_i) && ~isnan(cur_j) && cur_i>=1 && cur_j>=1 && cur_i<=numel(tf) && cur_j<=numel(tf)
            ti = double(tf(cur_i));
            tj = double(tf(cur_j));
            if isfinite(ti) && isfinite(tj)
                dt_corr = (ti - tj);
            end
        end

        % choose output weight
        switch lower(OUTPUT_WEIGHT_MODE)
            case 'base'
                w_out = max(weight_fun_base(cc_in),0);
            case 'robust'
                w_out = wrob;
            case 'combined'
                w_out = max(weight_fun_base(cc_in),0) * wrob;
            case 'thetastd'
                sf = std_fulls{g};
                s_i = double(sf(cur_i));
                s_j = double(sf(cur_j));
                std_dt = sqrt(max(0, s_i.^2 + s_j.^2)); % conservative (cov ignored)
                if std_dt<=0 || ~isfinite(std_dt)
                    w_out = 0;
                else
                    w_out = (1/std_dt) / thetaStdScale;
                    if isfinite(THETASTD_WEIGHT_CAP) && THETASTD_WEIGHT_CAP>0
                        w_out = min(w_out, THETASTD_WEIGHT_CAP);
                    end
                end
            otherwise
                w_out = max(weight_fun_base(cc_in),0) * wrob;
        end

        % write output line (only when we actually output something)
        if pendingHasHeader && ~pendingPrinted
            fprintf(fout,'%s\n', pendingHeader);
            pendingPrinted = true;
        end
        fprintf(fout, lineFmt, sta, dt_corr, w_out, ph);
    end

    fclose(fin);
    fclose(fout);
    closeDecisionReaders(decReaders);
end

out = struct();
out.maxEventID = maxEventID;
out.nGroups = nGroups;
out.groupKeys = groupKeys;
out.theta_dir = THETADIR;
out.thetastd_dir = THETASTD_DIR;
out.dt_sync = out_dt_sync;
out.metrics_file = metrics_file;
end


function maxID = readMaxEventID(catalog_file)
    fid = fopen(catalog_file,'r'); assert(fid>0);
    maxID = 0;
    while true
        L = fgetl(fid); if ~ischar(L), break; end
        L = strtrim(L); if isempty(L), continue; end
        toks = strsplit(L);
        id = str2double(toks{end});
        if ~isnan(id), maxID = max(maxID, id); end
    end
    fclose(fid);
end

function [keys, files, counts] = spoolGroups(infile_dt, TMPDIR, PROG_EVERY, PRINT_NEW_GROUPS)
    fid = fopen(infile_dt,'r'); assert(fid>0);
    keys = {}; files = {}; counts = [];
    key2idx = containers.Map('KeyType','char','ValueType','double');
    fids = {};

    cur_i = 0; cur_j = 0;
    nLine = 0;

    while true
        L = fgetl(fid); if ~ischar(L), break; end
        nLine = nLine + 1;
        if mod(nLine, PROG_EVERY)==0
            fprintf('  ... streamed %d lines\n', nLine);
        end
        Lt = strtrim(L); if isempty(Lt), continue; end

        if L(1)=='#'
            nums = sscanf(L(2:end),'%d %d %f');
            if numel(nums)>=2, cur_i=nums(1); cur_j=nums(2); end
            continue;
        end

        toks = textscan(L,'%s %f %f %s');
        if isempty(toks{1}), continue; end
        sta = toks{1}{1};
        dt  = toks{2}(1);
        cc  = toks{3}(1);
        ph  = toks{4}{1};
        key = [sta '_' ph];

        if ~isKey(key2idx, key)
            g = numel(keys) + 1;
            key2idx(key) = g;
            keys{g} = key;
            fn = fullfile(TMPDIR, sprintf('grp_%s.bin', sanitizeKey(key)));
            files{g} = fn;
            counts(g,1) = 0;
            fids{g} = fopen(fn,'w'); assert(fids{g}>0);
        else
            g = key2idx(key);
        end

        rec = single([cur_i, cur_j, dt, cc]);
        fwrite(fids{g}, rec, 'single');
        counts(g) = counts(g) + 1;
    end

    fclose(fid);
    for g=1:numel(fids), fclose(fids{g}); end
end

function X = readGroupFile(fn)
    fid = fopen(fn,'r'); assert(fid>0);
    X = fread(fid, [4, inf], 'single');
    fclose(fid);
end

function s = sanitizeKey(key)
    s = regexprep(key,'[^A-Za-z0-9_]','_');
end

function [sta, ph] = splitKey(key)
    k = strsplit(key,'_');
    sta = k{1}; ph = k{2};
end

function writeThetaTriples(fn, theta_full, ref_full)
    fid = fopen(fn,'w'); assert(fid>0);
    n = numel(theta_full);
    for ev=1:n
        fprintf(fid,'%d %.6f %.0f\n', ev, theta_full(ev), ref_full(ev));
    end
    fclose(fid);
end

function writeThetaStdWithDegree(fn, std_full, ref_full, deg_full, scale_fixed, w_cap, write_wcol, write_alt_wcol, nodew_full)
    fid = fopen(fn,'w'); assert(fid>0);
    n = numel(std_full);
    if nargin < 7, write_wcol = false; end
    if nargin < 8, write_alt_wcol = false; end
    if nargin < 9, nodew_full = NaN(n,1); end

    for ev=1:n
        s = std_full(ev);
        r = ref_full(ev);
        d = deg_full(ev);

        if write_wcol
            if ~isfinite(s) || s<=0
                wtheta = 0;
            else
                wtheta = (1/s) / scale_fixed;
                if isfinite(w_cap) && w_cap>0, wtheta = min(wtheta, w_cap); end
            end
            if write_alt_wcol
                fprintf(fid,'%d\t%.6f\t%.0f\t%u\t%.6f\t%.6f\n', ev, s, r, d, wtheta, nodew_full(ev));
            else
                fprintf(fid,'%d\t%.6f\t%.0f\t%u\t%.6f\n', ev, s, r, d, wtheta);
            end
        else
            fprintf(fid,'%d\t%.6f\t%.0f\t%u\n', ev, s, r, d);
        end
    end
    fclose(fid);
end


function writeDecisionFile(fn, keep_mask, wrob)
    fid = fopen(fn,'w'); assert(fid>0);
    M = numel(keep_mask);
    out = zeros(2,M,'single');
    out(1,:) = single(keep_mask(:).');
    out(2,:) = single(wrob(:).');
    fwrite(fid, out, 'single');
    fclose(fid);
end

function readers = openDecisionReaders(decFiles)
    n = numel(decFiles);
    readers = repmat(struct('fid',[], 'buf',[], 'idx',1, 'n',0, 'chunk',200000), n, 1);
    for g=1:n
        readers(g).fid = fopen(decFiles{g}, 'r');
        assert(readers(g).fid>0);
        readers(g).buf = zeros(2,0,'single');
        readers(g).idx = 1;
        readers(g).n   = 0;
    end
end

function closeDecisionReaders(readers)
    for g=1:numel(readers), fclose(readers(g).fid); end
end

function [keep, wrob, reader] = readNextDecision(reader)
    if reader.idx > reader.n
        x = fread(reader.fid, [2, reader.chunk], 'single');
        reader.buf = x;
        reader.n = size(x,2);
        reader.idx = 1;
        if reader.n==0
            keep=false; wrob=0; return;
        end
    end
    keep = (reader.buf(1, reader.idx) > 0.5);
    wrob = double(reader.buf(2, reader.idx));
    reader.idx = reader.idx + 1;
end

function scale = computeThetaStdMedianScale(infile_dt, key2g, decFiles, std_fulls)
    readers = openDecisionReaders(decFiles);
    fin = fopen(infile_dt,'r'); assert(fin>0);

    vals = [];
    vals_cap = 2e6;

    cur_i = NaN; cur_j = NaN;
    while true
        L = fgetl(fin); if ~ischar(L), break; end
        Lt = strtrim(L); if isempty(Lt), continue; end
        if L(1)=='#'
            nums = sscanf(L(2:end),'%d %d %f');
            if numel(nums)>=2, cur_i=nums(1); cur_j=nums(2); end
            continue;
        end
        toks = textscan(L,'%s %f %f %s');
        if isempty(toks{1}), continue; end
        sta = toks{1}{1}; ph=toks{4}{1};
        key = [sta '_' ph];
        if ~isKey(key2g,key), continue; end
        g = key2g(key);

        [keep,~,readers(g)] = readNextDecision(readers(g));
        if ~keep, continue; end

        sf = std_fulls{g};
        s_i = double(sf(cur_i));
        s_j = double(sf(cur_j));
        std_dt = sqrt(max(0, s_i.^2 + s_j.^2));
        std_dt = max(std_dt, 1e-12);
        vals(end+1,1) = 1/std_dt; %#ok<AGROW>

        if numel(vals) >= vals_cap, break; end
    end

    fclose(fin);
    closeDecisionReaders(readers);

    if isempty(vals), scale = 1; else, scale = median(vals); end
end

function [theta_full, ref_full, keep_mask, wrob_full, sigma_hat_global, std_full, deg_full, nodew_full, pseudo_info] = ...
    processStationPhaseGroup(gi, gj, gd, w_base, maxEventID, ...
        K_SIGMA, MIN_EDGES, IRLS_ITERS, C_HUBER, IRLS_REL_TOL, RIDGE_EPS, ...
        EXPORT_THETA_STD, THETASTD_MODE, THETASTD_FALLBACK, ...
        STD_PROBES, STD_PROBE_DIST, STD_BATCH, STD_REPORT_EVERY_BATCH, ...
        STD_MAX_NRED, STD_MIN_SIGMA, STD_MIN_DIAGREL, ...
        PSEUDO_WEIGHT_SOURCE, PSEUDO_WEIGHT_EPS, THETASTD_SCALE_FIXED)

    ev_all = unique([gi(:); gj(:)]);
    nNodes = numel(ev_all);
    [~, li] = ismember(gi, ev_all);
    [~, lj] = ismember(gj, ev_all);

    theta_local = NaN(nNodes,1);
    ref_local   = NaN(nNodes,1);
    std_local   = NaN(nNodes,1);
    deg_local   = zeros(nNodes,1,'uint32');
    nodew_local = NaN(nNodes,1); % optional per-node weight summary

    pseudo_info = struct('used_pseudo',false,'n_pseudo_comps',0);

    keep_mask = (w_base(:) > 0);     % never keep zero-weight edges
    wrob_full = ones(numel(gd),1);
    wrob_full(~keep_mask) = 0;


    valid_edges = (w_base(:) > 0);
    comps = uf_components(nNodes, li(valid_edges), lj(valid_edges));
    nComp = max(comps);

    sigma_hat_global = NaN;
    comp_sigma = NaN(nComp,1);
    comp_w     = zeros(nComp,1);


    for c=1:nComp
        nodes_c = find(comps==c);
        n_c = numel(nodes_c);

        % edges that lie within this connected component (fast label test)
    % Only consider positive-weight edges when deciding if a component is usable
    mE = (comps(li) == c) & (comps(lj) == c) & (w_base(:) > 0);

    if ~any(mE)
        % Nothing usable in this component: leave NaNs in theta/std, and drop edges (already dropped by keep_mask init)
        continue;
    end

    % Component microfamily gate:
    % Require (i) enough edges and (ii) at least one cycle (redundancy beyond a tree)
    nE = sum(mE);
    redundancy = nE - (n_c - 1);   % = 0 for a tree, >0 means at least one loop

    if (nE < MIN_EDGES) || (redundancy < 1)
        % Too small / no redundancy => do not trust theta/std for this component.
        % Drop its edges from dt output and keep NaNs in theta/std arrays.
        idx_e = find(mE);
        keep_mask(idx_e) = false;
        wrob_full(idx_e) = 0;
        deg_local(nodes_c) = uint32(0);
        continue;
    end


        a_all = li(mE); b_all = lj(mE);
        dt_k  = gd(mE);
        wb_k  = w_base(mE);

        map = zeros(nNodes,1);
        map(nodes_c) = 1:n_c;
        a = map(a_all); b = map(b_all);

        deg0 = accumarray([a(:); b(:)], 1, [n_c,1]);
        [~, pin] = max(deg0);
        keep_idx = setdiff(1:n_c, pin);

        theta_c = solvePinnedLS(n_c, a, b, dt_k, wb_k, pin, keep_idx, RIDGE_EPS);
        r = (theta_c(a) - theta_c(b)) - dt_k;

        % prune
        if numel(r) >= MIN_EDGES
            s0 = 1.4826*mad(r,1);
            if s0==0, s0 = max(STD_MIN_SIGMA, std(r)); end
            bad = abs(r) > K_SIGMA*s0;
        else
            bad = false(numel(r),1);
        end
        keep_c = ~bad;

        % LS pruned
        w_eff = wb_k .* keep_c;
        theta_c = solvePinnedLS(n_c, a, b, dt_k, w_eff, pin, keep_idx, RIDGE_EPS);

        % IRLS
        wrob_c = ones(numel(dt_k),1);
        if IRLS_ITERS>0
            Jprev = inf;
            for it=1:IRLS_ITERS
                r = (theta_c(a) - theta_c(b)) - dt_k;
                % robust scale on a "high-weight" subset to avoid inflation
                ww = w_eff(:);
                medw = median(ww(ww>0));
                use = (ww > 0) & (ww >= 0.5*medw);
                if any(use)
                    s = 1.4826*mad(r(use),1);
                else
                    s = 1.4826*mad(r(ww>0),1);
                end
                if s==0, s = max(STD_MIN_SIGMA, std(r(ww>0))); end
                rs = r./s;

                wrob_c = huberWeight(rs, C_HUBER);
                w_eff  = wb_k .* keep_c .* wrob_c;

                J = sum(w_eff .* (r.^2));  % weighted LS objective
                if isfinite(Jprev) && (Jprev - J) <= IRLS_REL_TOL*max(Jprev,eps)
                    break;
                end
                Jprev = J;

                theta_c = solvePinnedLS(n_c, a, b, dt_k, w_eff, pin, keep_idx, RIDGE_EPS);
            end
        end

        % final sigma_hat (again using higher-weight subset)
        r = (theta_c(a) - theta_c(b)) - dt_k;
        ww = w_eff(:);
        medw = median(ww(ww>0));
        use = (ww > 0) & (ww >= 0.5*medw);
        if any(use)
            sigma_hat = 1.4826*mad(r(use),1);
        else
            sigma_hat = 1.4826*mad(r(ww>0),1);
        end
        if sigma_hat==0, sigma_hat = max(STD_MIN_SIGMA, std(r(ww>0))); end
        sigma_hat = max(sigma_hat, STD_MIN_SIGMA);
        sigma_hat_global = sigma_hat;
        comp_sigma(c) = sigma_hat;

        % degree (kept graph, unweighted count)
        kept = (w_eff > 0);
        comp_w(c)     = sum(kept);
        deg_kept = accumarray([a(kept); b(kept)], 1, [n_c,1]);
        deg_local(nodes_c) = uint32(deg_kept);

        % theta std / pseudo std
        std_c = NaN(n_c,1);
        nodew_c = NaN(n_c,1);

        % Build a per-node weight summary if requested/needed.
        % For pseudo_weight: derive node weight from incident edge weights on kept edges.
        if ismember(lower(string(THETASTD_MODE)), ["pseudo_weight","hutch"]) || ismember(lower(string(THETASTD_FALLBACK)), ["pseudo_weight"])
            kept_for_nodew = (w_eff > 0);
            if any(kept_for_nodew)
                aa = a(kept_for_nodew); bb = b(kept_for_nodew);
                switch lower(string(PSEUDO_WEIGHT_SOURCE))
                    case "base"
                        ww = w_base_c(kept_for_nodew);
                    case "robust"
                        ww = wrob_c(kept_for_nodew);
                    otherwise % 'combined'
                        ww = w_eff(kept_for_nodew);
                end
                sumw = accumarray([aa; bb], [ww; ww], [n_c 1], @sum, 0);
                cntw = accumarray([aa; bb], 1, [n_c 1], @sum, 0);
                nodew_c = sumw ./ max(cntw, 1);
            else
                nodew_c(:) = 0;
            end
        end

        % Decide whether to run Hutch or pseudo.
        run_hutch = EXPORT_THETA_STD && strcmpi(THETASTD_MODE,'hutch') && (STD_PROBES>0) && (numel(keep_idx) <= STD_MAX_NRED);

        if EXPORT_THETA_STD && strcmpi(THETASTD_MODE,'none')
            std_c(:) = NaN;
        elseif run_hutch
            std_c = estimateThetaStdHutch(n_c, a, b, w_eff, pin, keep_idx, RIDGE_EPS, ...
                sigma_hat, STD_PROBES, STD_PROBE_DIST, STD_BATCH, STD_REPORT_EVERY_BATCH, STD_MIN_DIAGREL);
        else
            % Hutch skipped (either disabled, too large, or mode not hutch): use fallback
            use_mode = lower(string(THETASTD_MODE));
            if strcmpi(use_mode,'hutch')
                use_mode = lower(string(THETASTD_FALLBACK));
                if (numel(keep_idx) > STD_MAX_NRED)
                    warning('[THETASTD] Hutch std skipped (Nred=%d > STD_MAX_NRED=%d). Using fallback=%s.', numel(keep_idx), STD_MAX_NRED, use_mode);
                elseif (STD_PROBES<=0)
                    warning('[THETASTD] Hutch std disabled (probes<=0). Using fallback=%s.', use_mode);
                end
            end

            if EXPORT_THETA_STD && strcmpi(use_mode,'pseudo_degree')
                % Pseudo std from degree: sigma_hat / sqrt(deg)
                degv = max(double(deg_kept), 1);
                sig = max(sigma_hat, STD_MIN_SIGMA);
                std_c = sig ./ sqrt(degv);
                std_c(pin) = 0; % keep gauge semantics (theta(pin)=0)
                pseudo_info.used_pseudo = true;
                pseudo_info.n_pseudo_comps = pseudo_info.n_pseudo_comps + 1;
            elseif EXPORT_THETA_STD && strcmpi(use_mode,'pseudo_weight')
                % Pseudo std from node weight: std = 1/(w * THETASTD_SCALE_FIXED)
                wnode = max(nodew_c, PSEUDO_WEIGHT_EPS);
                std_c = 1 ./ (wnode * THETASTD_SCALE_FIXED);
                std_c(pin) = 0; % gauge
                pseudo_info.used_pseudo = true;
                pseudo_info.n_pseudo_comps = pseudo_info.n_pseudo_comps + 1;
            else
                % NaN fallback
                if EXPORT_THETA_STD
                    std_c(:)=NaN; std_c(pin)=0;
                end
            end
        end


        theta_local(nodes_c) = theta_c(:);
        ref_local(nodes_c) = nodes_c(pin);
        if EXPORT_THETA_STD, std_local(nodes_c) = std_c(:); end
        nodew_local(nodes_c) = nodew_c(:);

        idx_e = find(mE);
        keep_mask(idx_e) = keep_c;
        wrob_full(idx_e) = wrob_c;
        wrob_full(idx_e(~keep_c)) = 0;
    end

    % Representative sigma_hat per group (for reporting only).
    good = (comp_w > 0) & isfinite(comp_sigma);
    if any(good)
        sigma_hat_global = median(comp_sigma(good));
    end

    theta_full = NaN(maxEventID,1);
    ref_full   = NaN(maxEventID,1);
    std_full   = NaN(maxEventID,1);
    nodew_full = NaN(maxEventID,1);
    deg_full   = zeros(maxEventID,1,'uint32');

    theta_full(ev_all) = theta_local;
    maskRef = isfinite(ref_local) & (ref_local >= 1) & (ref_local <= numel(ev_all));
    ref_full(ev_all(maskRef)) = ev_all(ref_local(maskRef));
    if EXPORT_THETA_STD, std_full(ev_all) = std_local; end
    nodew_full(ev_all) = nodew_local;
    deg_full(ev_all)   = deg_local;
end

function theta = solvePinnedLS(n_c, a, b, dt_k, w, pin, keep_idx, RIDGE_EPS)
    a=a(:); b=b(:); dt_k=dt_k(:); w=w(:);
    theta = zeros(n_c,1);

    valid = (w>0);
    if ~any(valid), theta(pin)=0; return; end

    aa=a(valid); bb=b(valid); ww=w(valid); dd=dt_k(valid);

    Wadj = sparse([aa; bb],[bb; aa],[ww; ww], n_c, n_c);
    dsum = sum(Wadj,2);
    L    = spdiags(dsum,0,n_c,n_c) - Wadj;

    bvec = accumarray(aa, ww.*dd, [n_c,1]) - accumarray(bb, ww.*dd, [n_c,1]);

    if isempty(keep_idx), theta(pin)=0; return; end
    keep_idx = keep_idx(:);

    Lred = L(keep_idx,keep_idx) + RIDGE_EPS*speye(numel(keep_idx));
    theta_red = Lred \ bvec(keep_idx);

    theta(:)=0;
    theta(pin)=0;
    theta(keep_idx)=theta_red;
end

function w = huberWeight(rs, c)
    ar = abs(rs);
    w = ones(size(ar));
    idx = ar>c;
    w(idx) = c ./ max(ar(idx), eps);
end

function comps = uf_components(n, a, b)
    parent = 1:n;
    rankv  = zeros(1,n,'uint16');
    for k=1:numel(a)
        [parent, rankv] = uf_union(parent, rankv, a(k), b(k));
    end
    for i=1:n
        [parent, ~] = uf_find(parent, i);
    end
    [~,~,comps] = unique(parent);
end

function [parent, rankv] = uf_union(parent, rankv, x, y)
    [parent, rx] = uf_find(parent, x);
    [parent, ry] = uf_find(parent, y);
    if rx==ry, return; end
    if rankv(rx) < rankv(ry)
        parent(rx)=ry;
    elseif rankv(rx) > rankv(ry)
        parent(ry)=rx;
    else
        parent(ry)=rx;
        rankv(rx)=rankv(rx)+1;
    end
end

function [parent, root] = uf_find(parent, x)
    root = x;
    while parent(root)~=root
        root = parent(root);
    end
    while parent(x)~=x
        px = parent(x);
        parent(x) = root;
        x = px;
    end
end

function std_c = estimateThetaStdHutch(n_c, a, b, w_eff, pin, keep_idx, RIDGE_EPS, ...
                                      sigma_hat, K, dist, batch, reportEvery, minDiagRel)
    a=a(:); b=b(:); w_eff=w_eff(:);
    valid = (w_eff>0);
    aa=a(valid); bb=b(valid); ww=w_eff(valid);

    Wadj = sparse([aa; bb],[bb; aa],[ww; ww], n_c, n_c);
    dsum = sum(Wadj,2);
    L = spdiags(dsum,0,n_c,n_c) - Wadj;

    keep_idx = keep_idx(:);
    Lred = L(keep_idx,keep_idx) + RIDGE_EPS*speye(numel(keep_idx));

    % Prefer chol if possible (fast, deterministic)
    useChol = true;
    try
        R = chol(Lred,'lower');
    catch
        useChol = false;
    end

    nred = numel(keep_idx);
    diag_est = zeros(nred,1);
    nDone = 0;

    nBatches = ceil(K/batch);
    for bi=1:nBatches
        thisK = min(batch, K - nDone);
        if thisK<=0, break; end

        % Draw thisK probe vectors (columns). Use rademacher probes by default for lower variance.
        switch lower(dist)
            case 'rademacher'
                Z = sign(randn(nred, thisK));
                Z(Z==0) = 1;
            case 'gaussian'
                Z = randn(nred, thisK);
            otherwise
                error('STD_PROBE_DIST must be rademacher or gaussian');
        end

        % Solve Lred * X = Z for multiple RHS at once when possible (much faster than per-probe loops).
        if useChol
            X = R' \ (R \ Z);
        else
            X = zeros(nred, thisK);
            for t = 1:thisK
                z = Z(:,t);
                [x, flag] = pcg(Lred, z, 1e-8, 300);
                if flag ~= 0
                    % last-resort direct solve (slow but safe)
                    x = Lred \ z;
                end
                X(:,t) = x;
            end
        end

        % Hutchinson diagonal estimator: diag(Lred^{-1}) = E[ X .* Z ] (elementwise).
        diag_est = diag_est + sum(X .* Z, 2);

        nDone = nDone + thisK;

        if reportEvery>0 && mod(bi, reportEvery)==0
            tmp = diag_est / max(nDone,1);
            fprintf('    Hutch diag: %d/%d probes (median diag ~ %.3e)\n', nDone, K, median(abs(tmp)));
        end
    end

    diag_est = diag_est / max(nDone,1);

    % Avoid hard zeros due to finite-K noise:
    %  - don’t clamp to 0 aggressively
    %  - instead floor relative to median magnitude
    medDiag = median(abs(diag_est));
    floorDiag = max(medDiag * minDiagRel, 0);
    diag_est = max(diag_est, floorDiag);

    var_red = (sigma_hat^2) * diag_est;
    var_red = max(var_red, 0);
    std_red = sqrt(var_red);

    std_c = NaN(n_c,1);
    std_c(pin)=0;
    std_c(keep_idx)=std_red;
end
