function sample = hutch_samples(L,R,opts)
%HUTCH_SAMPLES Retain per-entry estimates and diagnostics at the probe budget.
n=size(L,1);mu=zeros(n,1);m2=zeros(n,1);count=0;reason='insufficient_probes';
for bi=1:ceil(opts.probes/opts.batch)
    k=min(opts.batch,opts.probes-count);
    try
        if strcmp(opts.probe_dist,'rademacher'),Z=2*(rand(n,k)>=0.5)-1;else,Z=randn(n,k);end
        X=R'\(R\Z);
        if ~all(ddsync.checked_solve_columns(L,X,Z))
            reason='solve_failed';break;
        end
        samples=X.*Z;bm=mean(samples,2);bm2=sum((samples-bm).^2,2);
        delta=bm-mu;total=count+k;
        m2=m2+bm2+delta.^2*(count*k/total);mu=mu+delta*(k/total);count=total;
    catch
        reason='solve_failed';break;
    end
    [reason,rel,good]=ddsync.hutch_diagnostic(mu,m2,count,opts.min_probes,opts.rel_tol);
    if opts.report_every_batch>0 && mod(bi,opts.report_every_batch)==0
        fprintf('    Hutch diag: %d/%d probes, %d/%d entries converged, max relative MC SE=%.3g\n',count,opts.probes,sum(good),n,rel);
    end
    if strcmp(reason,'converged'),break;end
end
[~,rel,good,relative]=ddsync.hutch_diagnostic(mu,m2,count,opts.min_probes,opts.rel_tol);
sample=struct('mu',mu,'m2',m2,'converged',good,'rel_se',relative,'reason',reason,'probes',count,'max_rel_se',rel);
end
