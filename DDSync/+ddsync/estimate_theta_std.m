function [out,info] = estimate_theta_std(n,a,b,w,pin,free,ridge,sigma,opts)
%ESTIMATE_THETA_STD Checked Hutchinson estimates and bounded exact fallback.
% Pseudo modes are filled by the caller. No failed variance is floored.
out=NaN(n,1); out(pin)=0;
info=struct('method',opts.mode,'reason','requested','probes',0,'max_rel_se',NaN);
valid=isfinite(w) & w>0; aa=a(valid); bb=b(valid); ww=w(valid);
W=sparse([aa;bb],[bb;aa],[ww;ww],n,n);
L=spdiags(sum(W,2),0,n,n)-W;
L=L(free,free)+ridge*speye(numel(free));
nr=numel(free);
if strcmp(info.method,'hutch')
    if opts.probes<=0
        info.reason='disabled_probes';
    elseif nr>opts.max_nred
        info.reason='size_limit';
    else
        try
            R=chol(L,'lower');
            mu=zeros(nr,1); m2=zeros(nr,1); count=0;
            info.reason='insufficient_probes';
            for bi=1:ceil(opts.probes/opts.batch)
                k=min(opts.batch,opts.probes-count);
                if strcmp(opts.probe_dist,'rademacher')
                    Z=2*(rand(nr,k)>=0.5)-1;
                else
                    Z=randn(nr,k);
                end
                X=R'\(R\Z);
                if any(~isfinite(X(:))) || any(max(abs(L*X-Z),[],1)>1e-8*max(1,max(abs(Z),[],1)))
                    info.reason='solve_failed'; break;
                end
                samples=X.*Z;
                bm=mean(samples,2); bm2=sum((samples-bm).^2,2);
                delta=bm-mu; total=count+k;
                m2=m2+bm2+delta.^2*(count*k/total);
                mu=mu+delta*(k/total); count=total;
                info.probes=count;
                [info.reason,info.max_rel_se]=ddsync.hutch_diagnostic(mu,m2,count,opts.min_probes,opts.rel_tol);
                if opts.report_every_batch>0 && mod(bi,opts.report_every_batch)==0
                    fprintf('    Hutch diag: %d/%d probes, max relative MC SE=%.3g (%s)\n',count,opts.probes,info.max_rel_se,info.reason);
                end
                if strcmp(info.reason,'converged'), break; end
            end
            if strcmp(info.reason,'converged')
                if ~isfinite(sigma), info.reason='invalid_noise_scale';
                else, out(free)=sigma*sqrt(mu); end
            end
        catch
            info.reason='solve_failed';
        end
    end
    if ~strcmp(info.reason,'converged'), info.method=opts.fallback; end
end
if strcmp(info.method,'exact')
    if nr>opts.exact_max_nred
        info.method='graph_bound'; info.reason=[info.reason '+exact_size_limit'];
    else
        try
            R=chol(L,'lower'); diagonal=zeros(nr,1);
            for first=1:opts.batch:nr
                last=min(first+opts.batch-1,nr); E=zeros(nr,last-first+1);
                for j=1:size(E,2), E(first+j-1,j)=1; end
                X=R'\(R\E);
                assert(all(isfinite(X(:))) && max(abs(L*X-E),[],'all')<=1e-8,'Exact solve failed');
                for j=1:size(E,2), diagonal(first+j-1)=X(first+j-1,j); end
            end
            assert(isfinite(sigma) && all(isfinite(diagonal) & diagonal>0),'Invalid variance');
            out(free)=sigma*sqrt(diagonal);
        catch
            info.method='graph_bound'; info.reason=[info.reason '+exact_failed'];
        end
    end
end
% Missing sigma can become a minimum-sigma constraint in downstream readers.
if strcmp(opts.mode,'hutch') && ismember(info.method,{'nan','none'})
    info.method='graph_bound'; info.reason=[info.reason '+legacy_nan_fallback'];
end
if strcmp(info.method,'graph_bound')
    out=ddsync.graph_bound_std(n,a,b,w,pin,sigma);
elseif strcmp(info.method,'none')
    out(:)=NaN; out(pin)=0;
end
end
