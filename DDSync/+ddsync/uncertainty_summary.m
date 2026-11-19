function summary = uncertainty_summary(components)
%UNCERTAINTY_SUMMARY Count event/station/phase std estimates, excluding gauges.
% exact includes full exact calculations and selected repairs; proxy/bound
% entries count as fallback even when explicitly requested as a mode.
summary=struct('exact',0,'hutch_converged',0,'hutch_unconverged',0, ...
    'fallback',0,'disabled',0,'references',numel(components),'total',0,'review',0, ...
    'review_percent',NaN,'advise',false);
for d=components
    counts=[d.exact_nodes,d.hutch_converged,d.unconverged_nodes,d.fallback_nodes];
    summary.exact=summary.exact+counts(1);
    summary.hutch_converged=summary.hutch_converged+counts(2);
    summary.hutch_unconverged=summary.hutch_unconverged+counts(3);
    summary.fallback=summary.fallback+counts(4);
    summary.disabled=summary.disabled+max(0,d.nodes-1-sum(counts));
end
summary.total=summary.exact+summary.hutch_converged+summary.hutch_unconverged+summary.fallback;
summary.review=summary.hutch_unconverged+summary.fallback;
if summary.total>0
    summary.review_percent=100*summary.review/summary.total;
    summary.advise=summary.review_percent>=10;
end
end
