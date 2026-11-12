function large_uncertainty_study(output_path)
%LARGE_UNCERTAINTY_STUDY Three seeds, 2048 probes, unit residual scale.
% Timings are illustrative, not a controlled benchmark.
addpath(fileparts(fileparts(mfilename('fullpath'))));
if nargin<1,output_path='large_uncertainty_matlab.csv';end
fid=fopen(output_path,'w');cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'graph,seed,policy,seconds,converged,exact,kept,fallback,median_std_ratio,p95_abs_std_error,max_abs_std_error,graph_bound_median_ratio,graph_bound_p95_ratio,reason,probes\n');
for graph={'mostly_good5000','weak_bridge5000','dense2400'}
    name=graph{1};[n,a,b,w]=graph_data(name);free=(2:n)';ridge=1e-10;
    W=sparse([a;b],[b;a],[w;w],n,n);L=spdiags(sum(W,2),0,n,n)-W;
    L=L(free,free)+ridge*speye(n-1);R=chol(L,'lower');
    oracle=sqrt(ddsync.selected_inverse_diagonal(L,R,1:n-1,128));assert(all(isfinite(oracle)));
    if strcmp(name,'dense2400'),assert(max(abs(oracle/sqrt(2/n)-1))<1e-6);end
    bound=ddsync.graph_bound_std(n,a,b,w,1,1);br=bound(free)./oracle;
    cfg=ddsync.config_default();opts=cfg.std.hutch;opts.mode='hutch';opts.fallback='pseudo_degree';
    opts.exact_max_nred=2048;opts.repair_max_nodes=128;opts.report_every_batch=0;
    for seed=1:3
        for policy={'keep','fallback'}
            opts.unconverged=policy{1};rng(seed);tic;
            [sd,info]=ddsync.estimate_theta_std(n,a,b,w,1,free,ridge,1,opts);elapsed=toc;
            ratio=sd(free)./oracle;err=abs(ratio-1);
            fprintf(fid,'%s,%d,%s,%.3f,%d,%d,%d,%d,%.6g,%.6g,%.6g,%.6g,%.6g,%s,%d\n', ...
                name,seed,policy{1},elapsed,info.hutch_converged,info.exact_nodes,info.unconverged_nodes,info.fallback_nodes, ...
                median(ratio),percentile95(err),max(err),median(br),percentile95(br),info.reason,info.probes);
        end
    end
end
end
function q=percentile95(x)
x=sort(x);index=1+0.95*(numel(x)-1);lo=floor(index);hi=ceil(index);q=x(lo)+(index-lo)*(x(hi)-x(lo));
end
function [n,a,b,w]=graph_data(name)
if strcmp(name,'dense2400')
    n=2400;[a,b]=find(triu(ones(n),1));w=ones(numel(a),1);return;
end
n=5000;a=[];b=[];w=[];
if strcmp(name,'mostly_good5000')
    split=4960;j=(2:split)';a=[ones(numel(j),1);(2:split-1)'];b=[j;(3:split)'];w=ones(numel(a),1);
    [aa,bb]=find(triu(ones(40),1));a=[a;aa+split;2];b=[b;bb+split;split+1];w=[w;ones(numel(aa),1);0.001];
else
    for start=[1,2501]
        j=(start+1:start+2499)';a=[a;repmat(start,numel(j),1);(start+1:start+2498)'];b=[b;j;(start+2:start+2499)']; %#ok<AGROW>
    end
    a=[a;1];b=[b;2501];w=ones(numel(a),1);
end
end
