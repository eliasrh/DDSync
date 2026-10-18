function run_tests(output_root)
%RUN_TESTS Native MATLAB regression tests on shared synthetic observations.
package_dir=fileparts(fileparts(mfilename('fullpath')));
addpath(package_dir);
if nargin<1
    output_root=tempname;mkdir(output_root);
    cleanup=onCleanup(@()rmdir(output_root,'s')); %#ok<NASGU>
end
checks=0;
for irls=[0,10]
    cfg=split_config(fullfile(output_root,sprintf('split_%d',irls)));
    cfg.robust.irls_iters=irls;
    evalc('ddsync.run(cfg);');
    for phase={'P','S'}
        ph=phase{1};
        theta=readmatrix(fullfile(cfg.io.thetadir,['theta_TEST_' ph '.txt']));
        sd=readmatrix(fullfile(cfg.io.thetastd_dir,['std_theta_TEST_' ph '.txt']));
        check(all(theta(1:10,3)==3));check(all(theta(11:20,3)==13));
        check(all(theta([3,13],2)==0));
        check(all(isnan(theta(21:23,2:3)),'all'));
        check(all(isnan(sd(21:23,2:3)),'all'));
        check(all(sd(21:23,4)==0));
        check(all(sd(setdiff(1:20,[3,13]),2)>1e-5));
        check(all(sd([1,2,11,12],4)==8));check(all(sd([3,13],4)==9));
        factor=1;if strcmp(ph,'S'),factor=1.7;end
        check(max(abs(theta(1:10,2)-0.013*((1:10)'-3)*factor))<0.0002);
    end
    lines=readlines(cfg.io.out_dt_sync); headers=lines(startsWith(lines,'#'));
    check(numel(headers)==88);
    for k=1:numel(headers)
        ids=sscanf(extractAfter(headers(k),1),'%d %d');i=ids(1);j=ids(2);
        check((i<=10 && j<=10)||(i>=11 && i<=20 && j>=11 && j<=20));
    end
    check(contains(fileread(cfg.io.metrics_file),'ref=3 nodes=10 method=exact'));
end
cfg=split_config(fullfile(output_root,'gate'));cfg.robust.min_edges=45;
evalc('ddsync.run(cfg);');
check(isempty(fileread(cfg.io.out_dt_sync)));
theta=readmatrix(fullfile(cfg.io.thetadir,'theta_TEST_P.txt'));
check(all(isnan(theta(:,2:3)),'all'));
cfg=split_config(fullfile(output_root,'fallback'));cfg.std.mode='hutch';cfg.std.hutch.probes=1;
evalc('ddsync.run(cfg);');
metrics=fileread(cfg.io.metrics_file);
check(contains(metrics,'method=exact'));check(contains(metrics,'probes=1'));
sd=readmatrix(fullfile(cfg.io.thetastd_dir,'std_theta_TEST_P.txt'));
check(all(sd(setdiff(1:20,[3,13]),2)>1e-5));
cfg=split_config(fullfile(output_root,'bound'));cfg.std.mode='hutch';cfg.std.hutch.probes=1;cfg.std.exact_max_nred=0;
evalc('ddsync.run(cfg);');
sd=readmatrix(fullfile(cfg.io.thetastd_dir,'std_theta_TEST_P.txt'));
check(all(isfinite(sd(1:20,2))));check(all(sd(setdiff(1:20,[3,13]),2)>0));
check(all(sd(setdiff(1:20,[3,13]),5)>0));check(all(isnan(sd(21:23,2))));
check(contains(fileread(cfg.io.metrics_file),'exact_size_limit'));
check(contains(fileread(cfg.io.metrics_file),'method=graph_bound'));
lines=readlines(cfg.io.out_dt_sync); obs=lines(~startsWith(lines,'#') & strlength(lines)>0);
for k=1:numel(obs), fields=split(strtrim(obs(k)));check(str2double(fields(3))>0);end
cfg=split_config(fullfile(output_root,'split_tree'));cfg.robust.min_edges=3;
cfg.io.infile_dt=fullfile(fileparts(cfg.io.infile_dt),'split_tree.cc');
evalc('ddsync.run(cfg);');
theta=readmatrix(fullfile(cfg.io.thetadir,'theta_TEST_P.txt'));
check(all(theta(1:3,3)==1));check(all(theta(4:6,3)==4));check(all(isnan(theta(7:23,2:3)),'all'));
cfg=split_config(fullfile(output_root,'median'));cfg.robust.min_edges=3;
cfg.io.infile_dt=fullfile(fileparts(cfg.io.infile_dt),'median.cc');
cfg.std.exact_max_nred=2;cfg.output.theta_weight_scale_mode='median';
evalc('ddsync.run(cfg);');
lines=readlines(cfg.io.out_dt_sync);obs=lines(~startsWith(lines,'#') & strlength(lines)>0);
for k=1:numel(obs)
    fields=split(strtrim(obs(k)));weight=str2double(fields(3));
    check(weight>0.5);
end
cfg=split_config(fullfile(output_root,'columns'));cfg.std.mode='pseudo_weight';
cfg.std.write_weight_column=false;cfg.std.write_alt_weight_column=true;
evalc('ddsync.run(cfg);');
sd=readmatrix(fullfile(cfg.io.thetastd_dir,'std_theta_TEST_P.txt'));
check(size(sd,2)==5);check(all(abs(sd(1:20,5)-0.9)<1e-6));
cfg=split_config(fullfile(output_root,'no_export'));cfg.std.export=false;cfg.output.dt_weight_mode='combined';
evalc('ddsync.run(cfg);');
check(~isfolder(cfg.io.thetastd_dir));check(~isempty(fileread(cfg.io.out_dt_sync)));
for k=1:3
    cfg=split_config(fullfile(output_root,sprintf('mode_%d',k)));cfg.output.write_pruned_edges=true;
    if k==1,cfg.std.mode='hutch';cfg.std.fallback='pseudo_degree';cfg.std.hutch.probes=0;end
    if k==2,cfg.std.mode='hutch';cfg.std.fallback='pseudo_weight';cfg.std.hutch.max_nred=0;end
    if k==3,cfg.std.mode='none';end
    evalc('ddsync.run(cfg);');
    lines=readlines(cfg.io.out_dt_sync);check(sum(startsWith(lines,'#'))==95);
    sd=readmatrix(fullfile(cfg.io.thetastd_dir,'std_theta_TEST_P.txt'));
    if k==3,check(all(isnan(sd(setdiff(1:20,[3,13]),2))));else,check(all(isfinite(sd(1:20,2))));end
end
% Independent incidence-matrix oracle for the conditional covariance.
a=[];b=[];
for i=1:9,for j=i+1:10,a(end+1,1)=i;b(end+1,1)=j;end,end %#ok<AGROW>
w=0.5+(1:numel(a))'/numel(a);free=2:10;ridge=1e-10;
A=zeros(numel(a),10);
for k=1:numel(a),A(k,a(k))=1;A(k,b(k))=-1;end
expected=sqrt(diag(inv(A(:,free)'*diag(w)*A(:,free)+ridge*eye(9))));
opts=struct('mode','exact','fallback','exact','probes',8192,'min_probes',128, ...
    'batch',128,'probe_dist','rademacher','rel_tol',0.05,'max_nred',40000, ...
    'exact_max_nred',512,'report_every_batch',0);
[sd,info]=ddsync.estimate_theta_std(10,a,b,w,1,free,ridge,1,opts);
check(strcmp(info.method,'exact'));check(max(abs(sd(free)./expected-1))<1e-12);
bound=ddsync.graph_bound_std(10,a,b,w,1,1);
check(all(bound(free)>=expected-1e-12));check(all(isfinite(bound)) && bound(1)==0);
tree=ddsync.graph_bound_std(4,[1;2;3],[2;3;4],[2;4;0.5],1,1);
check(max(abs(tree-sqrt([0;0.5;0.75;2.75])))<1e-12);
for dist={'rademacher','gaussian'}
    rng(42);opts.mode='hutch';opts.probe_dist=dist{1};
    [sd,info]=ddsync.estimate_theta_std(10,a,b,w,1,free,ridge,1,opts);
    check(strcmp(info.reason,'converged'));check(info.probes>=128 && info.probes<8192);
    check(info.max_rel_se<=0.05);check(max(abs(sd(free)./expected-1))<0.12);
end
[reason,~]=ddsync.hutch_diagnostic([-1;1],[1;1],128,128,0.1);check(strcmp(reason,'nonpositive_estimate'));
[reason,~]=ddsync.hutch_diagnostic([NaN;1],[1;1],128,128,0.1);check(strcmp(reason,'nonfinite_estimate'));
[reason,~]=ddsync.hutch_diagnostic(1,1,1,128,0.1);check(strcmp(reason,'insufficient_probes'));
[reason,~]=ddsync.hutch_diagnostic(1,1e8,128,128,0.1);check(strcmp(reason,'not_converged'));
opts.mode='exact';
failed=false;
try,ddsync.estimate_theta_std(3,1,2,1,1,[2,3],0,1,opts);catch,failed=true;end
check(failed); % cannot invent finite uncertainty for a disconnected system
failed=false;
try,ddsync.estimate_theta_std(10,a,b,w,1,free,ridge,NaN,opts);catch,failed=true;end
check(failed); % stop rather than export a usable theta with missing sigma
cfg=ddsync.config_default();cfg.std.hutch.batch=0;failed=false;
try,ddsync.validate_std_config(cfg.std);catch,failed=true;end
check(failed);
fprintf('MATLAB regression tests: %d checks passed.\n',checks);
    function check(condition)
        assert(all(condition,'all'),'Regression check %d failed',checks+1);
        checks=checks+1;
    end
end

function cfg=split_config(root)
if ~exist(root,'dir'),mkdir(root);end
package_dir=fileparts(fileparts(mfilename('fullpath')));
fixture=fullfile(package_dir,'..','test','fixtures');
cfg=ddsync.config_default();
cfg.io.infile_dt=fullfile(fixture,'split.cc');cfg.io.catalog_file=fullfile(fixture,'catalog.txt');
cfg.io.out_dt_sync=fullfile(root,'dt_sync.cc');cfg.io.metrics_file=fullfile(root,'metrics.txt');
cfg.io.tmpdir=fullfile(root,'tmp');cfg.io.thetadir=fullfile(root,'theta');cfg.io.thetastd_dir=fullfile(root,'thetastd');
cfg.robust.k_sigma=4;cfg.std.mode='exact';cfg.std.hutch.report_every_batch=0;
cfg.output.dt_decimals=9;cfg.output.theta_decimals=12;cfg.output.thetastd_decimals=12;
end
