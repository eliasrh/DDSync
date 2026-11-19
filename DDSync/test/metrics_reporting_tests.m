function checks=metrics_reporting_tests()
%METRICS_REPORTING_TESTS Percentages count estimates, not component averages.
addpath(fileparts(fileparts(mfilename('fullpath'))));
checks=0;
base=struct('group','LARGE_STATION_P','ref',1,'nodes',101,'method','hybrid', ...
    'reason','not_converged','probes',4096,'max_rel_se',1.2,'sigma_hat',1e-5,'sigma_std',5e-4, ...
    'repair_attempted',3,'repaired_nodes',3,'exact_nodes',70,'hutch_converged',20,'unconverged_nodes',5,'fallback_nodes',5);
s=ddsync.uncertainty_summary(base);
check(s.total==100 && s.exact==70 && s.hutch_converged==20 && s.hutch_unconverged==5 && s.fallback==5);
check(s.review_percent==10 && s.advise && s.references==1 && s.disabled==0);
text=evalc('ddsync.write_uncertainty_summary(1,s,'''');');
check(contains(text,'Exact (including repairs)') && contains(text,'Converged Hutch') && contains(text,'Kept unconverged Hutch') && contains(text,'Fallback / proxy'));
check(contains(text,'relative MC SE >= 1') && contains(text,'included in the fallback/proxy counts'));
check(contains(text,'70.00%') && contains(text,'20.00%') && contains(text,'5.00%') && contains(text,'ADVICE:'));
justbelow=base;justbelow.exact_nodes=71;justbelow.fallback_nodes=4;
s=ddsync.uncertainty_summary(justbelow);check(~s.advise);
allgood=base;allgood.nodes=91;allgood.exact_nodes=0;allgood.hutch_converged=90;allgood.unconverged_nodes=0;allgood.fallback_nodes=0;
text=evalc('ddsync.write_uncertainty_summary(1,ddsync.uncertainty_summary(allgood),'''');');
check(contains(text,'Converged Hutch') && contains(text,'100.00%') && ~contains(text,'ADVICE:'));
exact=allgood;exact.exact_nodes=90;exact.hutch_converged=0;exact.method='exact';exact.reason='requested';exact.probes=0;
s=ddsync.uncertainty_summary(exact);check(s.exact==90 && ~s.advise && s.review_percent==0);
disabled=exact;disabled.method='disabled';disabled.exact_nodes=0;
s=ddsync.uncertainty_summary(disabled);check(s.total==0 && s.disabled==90 && isnan(s.review_percent) && ~s.advise);
text=evalc('ddsync.write_uncertainty_summary(1,s,'''');');check(contains(text,'N/A') && ~contains(text,'100.00%'));
s=ddsync.uncertainty_summary(struct([]));check(s.total==0 && s.references==0 && ~s.advise);
% Ninety exact std and one fallback = 1/91, not an average of 0% and 100%.
small=base;small.group='VERY_LONG_STATION_NAME_PHASE_S';small.nodes=2;small.exact_nodes=0;small.hutch_converged=0;small.unconverged_nodes=0;small.fallback_nodes=1;
s=ddsync.uncertainty_summary([exact,small]);check(abs(s.review_percent-100/91)<1e-12 && ~s.advise);
group=struct('group',exact.group,'n_edges',1000,'n_kept',999,'pruned_pct',0.1,'r50',1,'r05',0.5,'r95',1,'frac_dw',0.1,'sigma_hat',1e-5);
group2=group;group2.group=small.group;
overall=struct('edges',2000,'kept',1998,'pseudo_groups',1,'pseudo_components',1);
path=tempname;cleanup=onCleanup(@()delete(path)); %#ok<NASGU>
fid=fopen(path,'w');ddsync.write_metrics_report(fid,[group,group2],[exact,small],overall);fclose(fid);
text=fileread(path);lines=splitlines(text);
check(strfind(text,'# UNCERTAINTY SUMMARY')<strfind(text,'# STATION-PHASE OBSERVATIONS'));
check(strfind(text,'# UNCERTAINTY BY GROUP')<strfind(text,'# COMPONENT DETAILS'));
check(strfind(text,'# COMPONENT DETAILS')<strfind(text,'# std_component'));
check(contains(text,'# GROUP ADVICE:') && contains(text,'#   VERY_LONG_STATION_NAME_PHASE_S'));
check(~contains(text,'# ADVICE:')); % low global percentage does not hide the bad group
check(contains(text,'sigma_hat=1e-05 sigma_std=0.0005'));
header=find(startsWith(lines,'group'),1);h=char(lines{header});r=char(lines{header+1});
he=regexp(h,'\S+','end');re=regexp(r,'\S+','end');check(isequal(he(2:end),re(2:end)));
check(~any(contains(lines(header+1:header+2),'#'))); % no component records interleaved
% All-converged details remain present even with no fallback reason.
allgood.method='hutch';allgood.reason='converged';allgood.probes=128;
fid=fopen(path,'w');ddsync.write_metrics_report(fid,group,allgood,overall);fclose(fid);
text=fileread(path);check(contains(text,'method=hutch') && contains(text,'reason=converged') && contains(text,'hutch_converged=90'));
% Empty input remains a legible report with no misleading 100% success.
fid=fopen(path,'w');ddsync.write_metrics_report(fid,struct([]),struct([]),overall);fclose(fid);
text=fileread(path);check(contains(text,'N/A') && contains(text,'# COMPONENT DETAILS') && ~contains(text,'100.00%'));
cfg=ddsync.config_default();check(cfg.std.hutch.probes==4096);
    function check(condition)
        assert(all(condition,'all'),'Metrics reporting regression check %d failed',checks+1);checks=checks+1;
    end
end
