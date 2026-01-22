%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% synchro.m  (authoritative DDSync script)
%
% Purpose:
%   Synchronize differential travel times (dt.cc) per station-phase by solving
%   for a node potential theta such that theta(i)-theta(j) best matches all
%   observed dt(i,j) in a weighted LS sense. Then robustly prune outliers and
%   optionally apply IRLS (Huber) to downweight remaining outliers.
%
% Memory strategy:
%   - Pass 1: stream dt.cc, spool per (station,phase) into temp binary files
%   - Process each group independently (no need to hold the full dt file)
%   - Pass 2: stream dt.cc again and write output dt_sync.cc using:
%       (i) per-group theta (random access by event id)
%       (ii) per-group edge decisions streamed from temp decision files
%
% Pipeline per (station,phase):
%   1) LS(all edges)
%   2) prune via robust MAD on residuals
%   3) RE-SPLIT into subcomponents after pruning (critical for stability)
%   4) LS(pruned) per subcomponent
%   5) IRLS(Huber, optional) per subcomponent
%   6) (optional) estimate std(theta) via Hutchinson diag(inv(Lred))
%
% Outputs:
%   out_dt_sync              synced dt for kept edges (HypoDD style)
%   theta_<STA>_<PH>.txt     [EventID theta refEventID]
%   thetastd/std_theta_*.txt [EventID std_theta refEventID degree]
%   metrics_file             per-group summary
%
% Output dt weights:
%   'base'     : weight_fun_base(cc)
%   'robust'   : huber factor only
%   'combined' : base * robust
%   'thetaStd' : weight = (1/std_dt)/scale, where std_dt is conservative
%                std_dt = sqrt(std_i^2 + std_j^2)  (cov ignored => conservative)
% Notes:
%   - Event IDs are assumed to be 1..N (or at least <= maxEventID), where
%     maxEventID is read from catalog.txt (last column).
%   - Spool files store [i j dt cc] as single precision for disk efficiency.
%     single stores integers exactly up to ~1.6e7, which is usually sufficient.
%   - If pruning disconnects the graph, we solve each subcomponent separately,
%     pinning each at its max-degree node.
%  NOTE on thetaStd calibration:
%   Because cov(theta_i,theta_j)>0 is common, sqrt(si^2+sj^2) usually OVERestimates
%   std(theta_i-theta_j). Expect conservative sigmas and z narrower than N(0,1).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc;

%% ---------------- USER KNOBS ----------------
infile_dt       = 'dt_long.cc';
catalog_file    = 'catalog.txt';
out_dt_sync      = 'dt_sync.cc';
metrics_file    = 'sync_metrics.txt';

% Base weight from CC (or uniform). Keep nonnegative.
weight_fun_base = @(cc) cc;  % @(cc) ones(size(cc)) or @(cc) cc.^2

% Pruning (robust)
K_SIGMA   = 20;
MIN_EDGES = 30; %sets a minimum on how many edges a group has that is then referenced to a reference event.
                %If small there may be theta value with a very small error estimate
                %if they belong to a microcluster. This causes additional
                %pruning in the dt.cc file

% IRLS (Huber); set 0 to disable
IRLS_ITERS   = 10;
C_HUBER      = 1.345;
IRLS_REL_TOL = 1e-3;

% Numerical stabilization
RIDGE_EPS = 1e-10;

% ---- Weight written into dt_sync.cc (third column) ----
OUTPUT_WEIGHT_MODE = 'thetaStd';   % 'base'|'robust'|'combined'|'thetaStd'

% Theta std estimation
EXPORT_THETA_STD = true;

% Hutchinson diag(inv(L)) estimator settings
STD_PROBES      = 100000;     % 0 disables std even if EXPORT_THETA_STD=true
STD_PROBE_DIST  = 'rademacher';% 'rademacher' or 'gaussian', which random distribution is
                               % used to stochastically estimate variance
STD_BATCH       = 250;       % report per batch, also reduces memory jitter
STD_REPORT_EVERY_BATCH = 10;  % print a progress line every N batches
STD_MAX_NRED    = 40000;     % skip std if reduced system larger than this
STD_MIN_SIGMA   = 5e-4;      % floor on sigma_hat (seconds)
STD_MIN_DIAGREL = 1e-12;     % floor on diag(inv(L)) relative to median(diag)

% thetaStd scaling for dt weights (writing only)
THETASTD_SCALE_MODE  = 'fixed'; % 'fixed' or 'median'
THETASTD_SCALE_FIXED = 500;    % typical 500 - 1000, 1000 means that a std of 0.001 s will give 
                                % a weight of 1, w = (1/std)/THETASTD_SCALE_FIXED
THETASTD_WEIGHT_CAP  = 1.0;     % cap on written weights (Inf to disable) often 
                                % weights are caped at 1 in DD relocations

% Formatting
DT_DECIMALS = 5;
W_DECIMALS  = 4;

% Temp folders
TMPDIR = 'synchro_tmp';
if ~exist(TMPDIR,'dir'), mkdir(TMPDIR); end
THETASTD_DIR = 'thetastd';
if EXPORT_THETA_STD && ~exist(THETASTD_DIR,'dir'), mkdir(THETASTD_DIR); end

% Streaming / progress
PROG_EVERY = 2e6;

%% ---------------- PRE-WARN on probes ----------------
if EXPORT_THETA_STD && STD_PROBES > 0
    fprintf('\n[STD WARNING]\n');
    fprintf('  Hutchinson diag(inv(L)) convergence is slow: estimator error ~ O(1/sqrt(K)).\n');
    fprintf('  If you want ms-level stability in std(theta), K may need to be VERY large.\n');
    fprintf('  For relocation weights, relative scaling often matters more than absolute calibration.\n\n');
end

%% ---------------- STEP 0: read catalog max ID ----------------
fprintf('Reading %s to get max EventID...\n', catalog_file);
maxEventID = readMaxEventID(catalog_file);
fprintf('Max EventID = %d\n', maxEventID);

%% ---------------- STEP 1: PASS 1 spool dt by (STA,PH) ----------------
fprintf('Pass 1: streaming %s and spooling groups...\n', infile_dt);
[groupKeys, groupFiles, groupCounts] = spoolGroups(infile_dt, TMPDIR, PROG_EVERY);
nGroups = numel(groupKeys);
fprintf('Found %d station+phase groups.\n', nGroups);

key2g = containers.Map('KeyType','char','ValueType','double');
for g=1:nGroups, key2g(groupKeys{g}) = g; end

%% ---------------- STEP 2: process groups ----------------
theta_fulls = cell(nGroups,1);
std_fulls   = cell(nGroups,1);
deg_fulls   = cell(nGroups,1);   % degree per event (kept graph)
decFiles    = cell(nGroups,1);

fidm = fopen(metrics_file,'w'); assert(fidm>0);
fprintf(fidm,'# synchro.m metrics\n');
fprintf(fidm,'# K_SIGMA=%g MIN_EDGES=%d IRLS_ITERS=%d C_HUBER=%g MODE=%s EXPORT_THETA_STD=%d STD_PROBES=%d DIST=%s\n', ...
    K_SIGMA, MIN_EDGES, IRLS_ITERS, C_HUBER, OUTPUT_WEIGHT_MODE, EXPORT_THETA_STD, STD_PROBES, STD_PROBE_DIST);
fprintf(fidm,'group,n_edges,n_kept,pruned_pct,robust_p50,robust_p05,robust_p95,frac_robust_lt_0p99,sigma_hat\n');

overall_edges = 0; overall_kept = 0;

for g=1:nGroups
    key = groupKeys{g};
    fprintf('Processing %s (%d/%d)\n', key, g, nGroups);

    X = readGroupFile(groupFiles{g});
    gi = double(X(1,:)).';
    gj = double(X(2,:)).';
    gd = double(X(3,:)).';
    gc = double(X(4,:)).';
    clear X;

    n_edges = numel(gd);
    overall_edges = overall_edges + n_edges;

    w_base = max(weight_fun_base(gc), 0);

    [theta_full, ref_full, keep_mask, wrob, sigma_hat, std_full, deg_full] = ...
        processStationPhaseGroup(gi, gj, gd, w_base, maxEventID, ...
            K_SIGMA, MIN_EDGES, IRLS_ITERS, C_HUBER, IRLS_REL_TOL, RIDGE_EPS, ...
            EXPORT_THETA_STD, STD_PROBES, STD_PROBE_DIST, STD_BATCH, STD_REPORT_EVERY_BATCH, ...
            STD_MAX_NRED, STD_MIN_SIGMA, STD_MIN_DIAGREL);

    n_kept = sum(keep_mask);
    overall_kept = overall_kept + n_kept;

    theta_fulls{g} = single(theta_full);
    std_fulls{g}   = single(std_full);
    deg_fulls{g}   = uint32(deg_full);

    [sta, ph] = splitKey(key);

    % write theta
    theta_fn = sprintf('theta_%s_%s.txt', sta, ph);
    writeThetaTriples(theta_fn, theta_full, ref_full);

    % write std + degree
    if EXPORT_THETA_STD
        std_fn = fullfile(THETASTD_DIR, sprintf('std_theta_%s_%s.txt', sta, ph));
        writeThetaStdWithDegree(std_fn, std_full, ref_full, deg_full);
    end

    % decisions file (in spooled-edge order)
    dec_fn = fullfile(TMPDIR, sprintf('dec_%s.bin', sanitizeKey(key)));
    writeDecisionFile(dec_fn, keep_mask, wrob);
    decFiles{g} = dec_fn;

    % metrics
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
    fprintf(fidm,'%s,%d,%d,%.3f,%.4f,%.4f,%.4f,%.3f,%.6g\n', ...
        key, n_edges, n_kept, pruned_pct, r50, r05, r95, frac_dw, sigma_hat);

    clear gi gj gd gc w_base theta_full ref_full keep_mask wrob sigma_hat std_full deg_full;
end

fprintf(fidm,'\n# OVERALL\n');
fprintf(fidm,'overall_edges=%d\n', overall_edges);
fprintf(fidm,'overall_kept=%d\n', overall_kept);
fprintf(fidm,'overall_pruned_pct=%.3f\n', 100*(overall_edges-overall_kept)/max(overall_edges,1));
fclose(fidm);

%% ---------------- STEP 3: optional median scaling for thetaStd weights ----------------
thetaStdScale = THETASTD_SCALE_FIXED;
if strcmpi(OUTPUT_WEIGHT_MODE,'thetaStd') && strcmpi(THETASTD_SCALE_MODE,'median')
    fprintf('Computing global median scale for thetaStd weights (extra streaming pass)...\n');
    thetaStdScale = computeThetaStdMedianScale(infile_dt, key2g, decFiles, std_fulls);
    fprintf('thetaStdScale (median of 1/std_dt) = %.6g\n', thetaStdScale);
    fprintf('If THETASTD_WEIGHT_CAP = 1, using median will saturate about half of the weights to 1')
    if thetaStdScale<=0 || ~isfinite(thetaStdScale), thetaStdScale = THETASTD_SCALE_FIXED; end
end

%% ---------------- STEP 4: PASS 2 write dt_sync ----------------
fprintf('Pass 2: writing %s ...\n', out_dt_sync);
decReaders = openDecisionReaders(decFiles);

fin  = fopen(infile_dt,'r'); assert(fin>0);
fout = fopen(out_dt_sync,'w'); assert(fout>0);

cur_i = NaN; cur_j = NaN;
fmt_dt = sprintf('%% .%df', DT_DECIMALS);
fmt_w  = sprintf('%%.%df', W_DECIMALS);

while true
    L = fgetl(fin); if ~ischar(L), break; end
    Lt = strtrim(L); if isempty(Lt), continue; end

    if L(1)=='#'
        fprintf(fout,'%s\n',L);
        nums = sscanf(L(2:end),'%d %d %f');
        if numel(nums)>=2, cur_i=nums(1); cur_j=nums(2); else, cur_i=NaN; cur_j=NaN; end
        continue;
    end

    toks = textscan(L,'%s %f %f %s');
    if isempty(toks{1}), continue; end
    sta = toks{1}{1};
    cc_in = toks{3}(1);
    ph  = toks{4}{1};

    key = [sta '_' ph];
    if ~isKey(key2g,key), continue; end
    g = key2g(key);

    [keep, wrob, decReaders(g)] = readNextDecision(decReaders(g));
    if ~keep, continue; end

    tf = theta_fulls{g};
    dt_sync = double(tf(cur_i) - tf(cur_j));
    if ~isfinite(dt_sync)
    continue;
    end
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
            std_dt = max(std_dt, 1e-12);
            w_raw = 1/std_dt;
            w_out = w_raw / max(thetaStdScale, eps);
        otherwise
            error('Bad OUTPUT_WEIGHT_MODE.');
    end

    if isfinite(THETASTD_WEIGHT_CAP)
        w_out = min(w_out, THETASTD_WEIGHT_CAP);
    end

    fprintf(fout, ['    %-4s ' fmt_dt ' ' fmt_w ' %s\n'], sta, dt_sync, w_out, ph);
end

fclose(fin); fclose(fout);
closeDecisionReaders(decReaders);

fprintf('Done.\n  Synced dt: %s\n  Metrics:  %s\n', out_dt_sync, metrics_file);
if EXPORT_THETA_STD, fprintf('  Theta std dir: %s/\n', THETASTD_DIR); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------------------------- LOCAL FUNCTIONS ------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

function [keys, files, counts] = spoolGroups(infile_dt, TMPDIR, PROG_EVERY)
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

function writeThetaStdWithDegree(fn, std_full, ref_full, deg_full)
    fid = fopen(fn,'w'); assert(fid>0);
    n = numel(std_full);
    for ev=1:n
        fprintf(fid,'%d %.6f %.0f %u\n', ev, std_full(ev), ref_full(ev), deg_full(ev));
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

function [theta_full, ref_full, keep_mask, wrob_full, sigma_hat_global, std_full, deg_full] = ...
    processStationPhaseGroup(gi, gj, gd, w_base, maxEventID, ...
        K_SIGMA, MIN_EDGES, IRLS_ITERS, C_HUBER, IRLS_REL_TOL, RIDGE_EPS, ...
        EXPORT_THETA_STD, STD_PROBES, STD_PROBE_DIST, STD_BATCH, STD_REPORT_EVERY_BATCH, ...
        STD_MAX_NRED, STD_MIN_SIGMA, STD_MIN_DIAGREL)

    ev_all = unique([gi(:); gj(:)]);
    nNodes = numel(ev_all);
    [~, li] = ismember(gi, ev_all);
    [~, lj] = ismember(gj, ev_all);

    theta_local = NaN(nNodes,1);
    ref_local   = NaN(nNodes,1);
    std_local   = NaN(nNodes,1);
    deg_local   = zeros(nNodes,1,'uint32');

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

        % theta std
        std_c = NaN(n_c,1);
        if (numel(keep_idx) > STD_MAX_NRED)
        warning('STD_MAX_NRED threshold is met, consider increasing')
        end
        if EXPORT_THETA_STD && (STD_PROBES>0) && (numel(keep_idx) <= STD_MAX_NRED)
            std_c = estimateThetaStdHutch(n_c, a, b, w_eff, pin, keep_idx, RIDGE_EPS, ...
                sigma_hat, STD_PROBES, STD_PROBE_DIST, STD_BATCH, STD_REPORT_EVERY_BATCH, STD_MIN_DIAGREL);
        else
            if EXPORT_THETA_STD
                std_c(:)=NaN; std_c(pin)=0;
            end
        end

        theta_local(nodes_c) = theta_c(:);
        ref_local(nodes_c) = nodes_c(pin);
        if EXPORT_THETA_STD, std_local(nodes_c) = std_c(:); end

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
    deg_full   = zeros(maxEventID,1,'uint32');

    theta_full(ev_all) = theta_local;
    maskRef = isfinite(ref_local) & (ref_local >= 1) & (ref_local <= numel(ev_all));
    ref_full(ev_all(maskRef)) = ev_all(ref_local(maskRef));
    if EXPORT_THETA_STD, std_full(ev_all) = std_local; end
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
