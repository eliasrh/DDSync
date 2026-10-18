function probe_study(output_file)
%PROBE_STUDY Reproducible sampling study, separate from fast regression tests.
% From repository root: addpath('DDSync/test'); probe_study('study.csv')
addpath(fileparts(fileparts(mfilename('fullpath'))));
fid=fopen(output_file,'w');assert(fid>0);cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'graph,budget,accepted,seeds,median_used,max_relative_std_error_accepted\n');
opts=struct('mode','exact','fallback','nan','probes',2048,'min_probes',128, ...
    'batch',128,'probe_dist','rademacher','rel_tol',0.10,'max_nred',40000, ...
    'exact_max_nred',512,'report_every_batch',0);
for graph={'dense20','dense100','ring100','weak_bridge100','hub100'}
    graph=graph{1};n=100;if strcmp(graph,'dense20'),n=20;end
    a=[];b=[];w=[];
    for i=1:n-1
        for j=i+1:n
            keep=startsWith(graph,'dense') || (strcmp(graph,'ring100') && (j==i+1 || (i==1 && j==n))) || ...
                (strcmp(graph,'hub100') && (i==1 || j==i+1)) || ...
                (strcmp(graph,'weak_bridge100') && ((i<=50 && j<=50) || (i>50 && j>50) || (i==50 && j==51)));
            if ~keep,continue;end
            a(end+1,1)=i;b(end+1,1)=j;weight=1; %#ok<AGROW>
            if strcmp(graph,'weak_bridge100') && i==50 && j==51,weight=1e-3;end
            w(end+1,1)=weight; %#ok<AGROW>
        end
    end
    degree=accumarray([a;b],1,[n,1]);[~,pin]=max(degree);free=setdiff(1:n,pin);
    opts.mode='exact';[oracle,info]=ddsync.estimate_theta_std(n,a,b,w,pin,free,1e-10,1,opts);
    assert(strcmp(info.method,'exact'));
    for budget=[128,512,2048,8192]
        accepted=0;used=[];errors=[];opts.mode='hutch';opts.probes=budget;
        for seed=1:10
            rng(seed);
            [sd,info]=ddsync.estimate_theta_std(n,a,b,w,pin,free,1e-10,1,opts);
            used(end+1)=info.probes; %#ok<AGROW>
            if strcmp(info.reason,'converged')
                accepted=accepted+1;errors(end+1)=max(abs(sd(free)./oracle(free)-1)); %#ok<AGROW>
            end
        end
        maxerr=NaN;if ~isempty(errors),maxerr=max(errors);end
        fprintf(fid,'%s,%d,%d,10,%.0f,%.6f\n',graph,budget,accepted,median(used),maxerr);
    end
end
end
