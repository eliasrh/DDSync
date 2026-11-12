function [reason,rel,good,relative] = hutch_diagnostic(mu,m2,count,min_probes,rel_tol)
%HUTCH_DIAGNOSTIC Per-entry empirical MC error, not a confidence bound.
relative=Inf(size(mu));
valid=isfinite(mu) & mu>0 & isfinite(m2) & m2>=0;
if count>=max(2,min_probes)
    relative(valid)=sqrt(m2(valid)/((count-1)*count))./mu(valid);
end
good=valid & relative<=rel_tol;rel=max(relative);
if any(~isfinite(mu)) || any(~isfinite(m2)),reason='nonfinite_estimate';
elseif any(mu<=0),reason='nonpositive_estimate';
elseif count<max(2,min_probes),reason='insufficient_probes';
elseif all(good),reason='converged';
else,reason='not_converged';end
end
