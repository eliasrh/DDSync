function write_uncertainty_summary(fid,summary,prefix)
%WRITE_UNCERTAINTY_SUMMARY Always report all categories, including success.
if nargin<3,prefix='# ';end
fprintf(fid,'%sUNCERTAINTY SUMMARY\n',prefix);
fprintf(fid,'%sDenominator: %d supported non-reference std entries (station-phase/event rows).\n',prefix,summary.total);
fprintf(fid,'%sExcluded: %d gauge references; %d disabled entries; rejected/absent events.\n',prefix,summary.references,summary.disabled);
fprintf(fid,'%s%-30s %12s %12s\n',prefix,'Category','Count','Percent');
labels={'Exact (including repairs)','Converged Hutch','Kept unconverged Hutch','Converged Cholesky','Kept unconverged Cholesky','Fallback / proxy'};
counts=[summary.exact,summary.hutch_converged,summary.hutch_unconverged,summary.cholesky_converged,summary.cholesky_unconverged,summary.fallback];
for k=1:numel(labels)
    fprintf(fid,'%s%-30s %12d %12s\n',prefix,labels{k},counts(k),pct(counts(k),summary.total));
end
fprintf(fid,'%sA draw is one random sample. Relative MC SE is the estimated sampling error divided by the event timing variance estimate.\n',prefix);
fprintf(fid,'%sFor example, 0.10 means 10%% sampling error on the variance. The uncertainty exported in seconds is its square root after scaling.\n',prefix);
fprintf(fid,'%sRequires review (kept unconverged + fallback/proxy): %d (%s)\n',prefix,summary.review,pct(summary.review,summary.total));
fprintf(fid,'%sNumerically checked (exact + converged estimates): %d (%s)\n',prefix,summary.exact+summary.hutch_converged+summary.cholesky_converged,pct(summary.exact+summary.hutch_converged+summary.cholesky_converged,summary.total));
fprintf(fid,'%sFallback/proxy includes explicitly requested proxies/bounds. Checked does not mean physically calibrated.\n',prefix);
fprintf(fid,'%sSampled entries with relative MC SE >= 1, invalid estimates or too few probes cannot be kept;\n',prefix);
fprintf(fid,'%safter selected repairs these entries use fallback. They are included in the fallback/proxy counts above.\n',prefix);
if summary.total==0
    fprintf(fid,'%sNo assessed std entries; percentages are N/A. Enable uncertainty if required by the application.\n',prefix);
elseif summary.advise
    fprintf(fid,'%sADVICE: At least 10%% require review. Check the categories above and std_quality_*.txt.\n',prefix);
    fprintf(fid,'%s        If many Hutch estimates miss the target, try Cholesky before increasing the number of probes.\n',prefix);
    fprintf(fid,'%s        Set cfg.std.auto_probe = ''cholesky'' in auto mode, or change explicit hutch mode to cfg.std.mode = ''cholesky''.\n',prefix);
    fprintf(fid,'%s        For a few unresolved entries, consider raising cfg.std.repair_max_nodes. For fallbacks, check the component reason.\n',prefix);
else
    fprintf(fid,'%sBelow the 10%% review advisory threshold. Inspect flagged entries and group summaries before weighting.\n',prefix);
end
end
function value=pct(count,total)
if total==0,value='N/A';else,value=sprintf('%.2f%%',100*count/total);end
end
