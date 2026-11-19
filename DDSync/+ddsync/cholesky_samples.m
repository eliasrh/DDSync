function sample = cholesky_samples(A, lower, opts)
%CHOLESKY_SAMPLES Estimate each timing variance using squared random draws.
% A is the precision with the reference event removed. The shared response
% v = A \ ones has covariance contribution v*v'/sum(v). Use its diagonal
% directly, and draw only the remaining uncertainty. No dense inverse is built.
n = size(A,1);
mu = zeros(n,1);
m2 = zeros(n,1);
count = 0;
batches = 0;
reason = 'insufficient_probes';
upper = lower';
rhs = ones(n,1);
v = upper \ (lower \ rhs);
s = sum(v);
if ~ddsync.checked_solve_columns(A,v,rhs) || ~isfinite(s) || s <= 0
    reason = 'solve_failed';
else
    weight = v/s;
    shared_variance = (v/sqrt(s)).^2;
    check_at = min(opts.probes,max(256,opts.min_probes));
    while count < opts.probes
        k = min(opts.batch,check_at-count);
        if strcmp(opts.probe_dist,'rademacher')
            z = 2*(rand(n,k)>0.5)-1;
        else
            z = randn(n,k);
        end
        draws = upper \ z;
        if ~all(ddsync.checked_solve_columns(upper,draws,z))
            reason = 'solve_failed';
            break
        end
        for j = 1:k
            y = draws(:,j);
            variance = (y-weight*sum(y)).^2 + shared_variance;
            count = count + 1;
            delta = variance-mu;
            mu = mu + delta/count;
            m2 = m2 + delta.*(variance-mu);
        end
        batches = batches + 1;
        if count == check_at
            reason = ddsync.hutch_diagnostic(mu,m2,count,opts.min_probes,opts.rel_tol);
            if strcmp(reason,'converged'), break; end
            check_at = count + min(count,opts.probes-count);
        end
        if opts.report_every_batch>0 && mod(batches,opts.report_every_batch)==0
            fprintf('    Cholesky sampling: %d/%d draws completed\n',count,opts.probes);
        end
    end
end
[~,rel,good,relative] = ddsync.hutch_diagnostic(mu,m2,count,opts.min_probes,opts.rel_tol);
sample = struct('mu',mu,'m2',m2,'converged',good,'rel_se',relative, ...
    'reason',reason,'probes',count,'max_rel_se',rel);
end
