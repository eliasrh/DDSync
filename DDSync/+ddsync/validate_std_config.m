function validate_std_config(cfg)
%VALIDATE_STD_CONFIG Reject invalid estimator options before writing outputs.
assert(ismember(lower(cfg.mode),{'hutch','exact','graph_bound','pseudo_degree','pseudo_weight','none'}),'Invalid std.mode');
assert(ismember(lower(cfg.fallback),{'exact','graph_bound','pseudo_degree','pseudo_weight','nan','none'}),'Invalid std.fallback');
assert(ismember(lower(cfg.hutch.probe_dist),{'rademacher','gaussian'}),'Invalid std.hutch.probe_dist');
v=[cfg.hutch.probes,cfg.hutch.min_probes,cfg.hutch.batch,cfg.hutch.max_nred,cfg.exact_max_nred];
assert(all(isfinite(v) & v>=0 & v==round(v)),'Probe counts and size limits must be nonnegative integers');
assert(cfg.hutch.batch>=1 && cfg.hutch.min_probes>=2,'batch must be positive and min_probes at least 2');
assert(isfinite(cfg.hutch.rel_tol) && cfg.hutch.rel_tol>0,'std.hutch.rel_tol must be finite and positive');
end
