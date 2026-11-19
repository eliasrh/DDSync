function write_metrics_report(fid,groups,components,overall)
%WRITE_METRICS_REPORT Readable report with uninterrupted observation rows.
summary=ddsync.uncertainty_summary(components);
group_summaries=summarize_groups(groups,components);
fprintf(fid,'\n');ddsync.write_uncertainty_summary(fid,summary);
fprintf(fid,'# sigma_hat: raw final residual-noise scale in seconds; each group reports the median of component scales.\n');
fprintf(fid,'# sigma_std: component sigma_hat after the optional min_sigma floor; this is sigma in the std formulas.\n');
fprintf(fid,'# Neither scale is a theta std or a Monte Carlo error. pseudo_weight does not use sigma_std.\n');
fprintf(fid,'\n# STATION-PHASE OBSERVATIONS\n');
width=16;
for g=groups,width=max(width,numel(g.group));end
head=['%-' num2str(width) 's %10s %10s %11s %12s %12s %12s %20s %13s\n'];
row=['%-' num2str(width) 's %10d %10d %11.3f %12.4f %12.4f %12.4f %20.3f %13.6g\n'];
fprintf(fid,head,'group','n_edges','n_kept','pruned_pct','robust_p50','robust_p05','robust_p95','frac_robust_lt_0p99','sigma_hat');
for g=groups
    fprintf(fid,row,g.group,g.n_edges,g.n_kept,g.pruned_pct,g.r50,g.r05,g.r95,g.frac_dw,g.sigma_hat);
end
fprintf(fid,'\n# OVERALL\n');
fprintf(fid,'overall_edges=%d\noverall_kept=%d\noverall_pruned_pct=%.3f\n',overall.edges,overall.kept,100*(overall.edges-overall.kept)/max(overall.edges,1));
fprintf(fid,'pseudo_groups_used=%d\npseudo_components_used=%d\n',overall.pseudo_groups,overall.pseudo_components);
fprintf(fid,'\n# UNCERTAINTY BY GROUP (counts and percentages use each group''s assessed std entries)\n');
head=['# %-' num2str(width) 's %8s %18s %18s %18s %18s %12s %10s\n'];
row=['# %-' num2str(width) 's %8d %18s %18s %18s %18s %12s %10d\n'];
fprintf(fid,head,'group','n_std','exact','hutch_converged','hutch_unconverged','fallback/proxy','review_pct','disabled');
for k=1:numel(groups)
    g=groups(k);q=group_summaries(k);
    fprintf(fid,row,g.group,q.total,cellstr_count(q.exact,q.total),cellstr_count(q.hutch_converged,q.total), ...
        cellstr_count(q.hutch_unconverged,q.total),cellstr_count(q.fallback,q.total),pct(q.review,q.total),q.disabled);
end
flagged={};
for k=1:numel(groups)
    g=groups(k);q=group_summaries(k);
    if q.advise,flagged{end+1}=g.group;end %#ok<AGROW>
end
if ~isempty(flagged)
    fprintf(fid,'# GROUP ADVICE: At least 10%% require review in each group below, even if the overall share is small.\n');
    fprintf(fid,'# Increase std.hutch.probes and/or std.repair_max_nodes, or std.exact_max_nred where affordable;\n');
    fprintf(fid,'# inspect skip/solve reasons and explicit proxy settings before rerunning.\n');
    for k=1:numel(flagged),fprintf(fid,'#   %s\n',flagged{k});end
end
fprintf(fid,'\n# COMPONENT DETAILS (all methods, including converged Hutch)\n');
fprintf(fid,'# max_rel_se is the pre-repair diagnostic, not the error of exact repairs; scales are seconds.\n');
for d=components
    fprintf(fid,'# std_component group=%s ref=%d nodes=%d method=%s\n',d.group,d.ref,d.nodes,d.method);
    fprintf(fid,'#   reason=%s probes=%d max_rel_se=%.6g\n',d.reason,d.probes,d.max_rel_se);
    fprintf(fid,'#   sigma_hat=%.6g sigma_std=%.6g repair_attempted=%d repaired_nodes=%d\n',d.sigma_hat,d.sigma_std,d.repair_attempted,d.repaired_nodes);
    fprintf(fid,'#   exact_nodes=%d hutch_converged=%d unconverged_nodes=%d fallback_nodes=%d\n\n',d.exact_nodes,d.hutch_converged,d.unconverged_nodes,d.fallback_nodes);
end
end
function value=pct(count,total)
if total==0,value='N/A';else,value=sprintf('%.2f%%',100*count/total);end
end
function value=cellstr_count(count,total)
value=sprintf('%d (%s)',count,pct(count,total));
end

function summaries=summarize_groups(groups,components)
% Pre-index groups once: reporting remains O(groups + components) for DAS.
summaries=repmat(ddsync.uncertainty_summary(struct([])),1,numel(groups));
index=containers.Map('KeyType','char','ValueType','double');
for k=1:numel(groups),index(groups(k).group)=k;end
fields={'exact','hutch_converged','hutch_unconverged','fallback','disabled','references','total','review'};
for d=components
    if ~isKey(index,d.group),continue;end
    k=index(d.group);part=ddsync.uncertainty_summary(d);q=summaries(k);
    for f=1:numel(fields),name=fields{f};q.(name)=q.(name)+part.(name);end
    summaries(k)=q;
end
for k=1:numel(summaries)
    if summaries(k).total>0
        summaries(k).review_percent=100*summaries(k).review/summaries(k).total;
        summaries(k).advise=summaries(k).review_percent>=10;
    end
end
end
