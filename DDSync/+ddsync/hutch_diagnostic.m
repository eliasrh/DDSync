function [reason,rel] = hutch_diagnostic(mu,m2,count,min_probes,rel_tol)
%HUTCH_DIAGNOSTIC Per-node empirical MC error, not a confidence bound.
rel=Inf;
if any(~isfinite(mu)) || any(~isfinite(m2))
    reason='nonfinite_estimate';
elseif any(mu<=0)
    reason='nonpositive_estimate';
elseif count<max(2,min_probes)
    reason='insufficient_probes';
else
    se=sqrt(max(m2,0)/((count-1)*count));
    rel=max(se./mu);
    if rel<=rel_tol, reason='converged'; else, reason='not_converged'; end
end
end
