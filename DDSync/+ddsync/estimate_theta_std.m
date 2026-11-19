function [out,info] = estimate_theta_std(n,a,b,w,pin,free,ridge,sigma,opts,pseudo_weight)
%ESTIMATE_THETA_STD Per-event Hutch decisions, selected exact repair and fallback.
if nargin<10,pseudo_weight=NaN(n,1);end
if ~isfield(opts,'repair_max_nodes'),opts.repair_max_nodes=128;end
if ~isfield(opts,'unconverged'),opts.unconverged='keep';end
if ~isfield(opts,'auto_probe'),opts.auto_probe='hutch';end
free=free(:);nr=numel(free);out=NaN(n,1);out(pin)=0;
methods=repmat({'unavailable'},n,1);methods{pin}='reference';relative=NaN(n,1);
reason='requested';probes=0;max_rel_se=NaN;repaired=0;attempted=0;mode=opts.mode;
if strcmp(mode,'auto')
    if nr<=opts.exact_max_nred,mode='exact';else,mode=opts.auto_probe;end
end
if ~strcmp(mode,'none')
    assert(isfinite(sigma) && sigma>=0,'Uncertainty needs a finite nonnegative noise scale; inspect residuals and std.min_sigma.');
end
fallback=opts.fallback;
if ismember(fallback,{'exact','nan','none'}),fallback='pseudo_degree';end
if ismember(mode,{'hutch','cholesky'}) && (opts.probes<=0 || nr>opts.max_nred) && strcmp(opts.fallback,'exact') && nr<=opts.exact_max_nred
    if opts.probes<=0,reason='disabled_probes';else,reason='size_limit';end
    mode='exact';
end
if ismember(mode,{'exact','hutch','cholesky'})
    if strcmp(mode,'exact') && nr>opts.exact_max_nred
        reason='exact_size_limit';
    elseif ismember(mode,{'hutch','cholesky'}) && (opts.probes<=0 || nr>opts.max_nred)
        if opts.probes<=0,reason='disabled_probes';else,reason='size_limit';end
    else
        valid=isfinite(w) & w>0;aa=a(valid);bb=b(valid);ww=w(valid);
        W=sparse([aa;bb],[bb;aa],[ww;ww],n,n);L=spdiags(sum(W,2),0,n,n)-W;
        L=L(free,free)+ridge*speye(nr);
        try,R=chol(L,'lower');catch,R=[];end
        if isempty(R)
            reason='solve_failed';
        elseif strcmp(mode,'exact')
            diagonal=ddsync.selected_inverse_diagonal(L,R,1:nr,opts.batch);ok=isfinite(diagonal);
            out(free(ok))=sigma*sqrt(diagonal(ok));methods(free(ok))={'exact'};
            if ~all(ok),reason='exact_failed';end
        else
            if strcmp(mode,'cholesky')
                sample=ddsync.cholesky_samples(L,R,opts);
            else
                sample=ddsync.hutch_samples(L,R,opts);
            end
            reason=sample.reason;probes=sample.probes;max_rel_se=sample.max_rel_se;
            relative(free)=sample.rel_se;good=sample.converged;
            out(free(good))=sigma*sqrt(sample.mu(good));methods(free(good))={mode};
            bad=find(~good);
            [~,order]=sortrows([-sample.rel_se(bad),bad],[1,2]);bad=bad(order);
            limit=opts.repair_max_nodes;
            if strcmp(opts.fallback,'exact') && nr<=opts.exact_max_nred,limit=numel(bad);end
            targets=bad(1:min(numel(bad),limit));attempted=numel(targets);
            exact=ddsync.selected_inverse_diagonal(L,R,targets,opts.batch);
            for j=1:numel(targets)
                i=targets(j);
                if isfinite(exact(j)),out(free(i))=sigma*sqrt(exact(j));methods{free(i)}='exact';repaired=repaired+1;end
            end
            if strcmp(opts.unconverged,'keep') && probes>=opts.min_probes
                for j=1:numel(bad)
                    i=bad(j);
                    % Near-zero estimates with SE >= estimate are not retained.
                    if strcmp(methods{free(i)},'unavailable') && isfinite(sample.mu(i)) && sample.mu(i)>0 && sample.rel_se(i)<1
                        out(free(i))=sigma*sqrt(sample.mu(i));methods{free(i)}=[mode '_unconverged'];
                    end
                end
            end
        end
    end
elseif strcmp(mode,'none')
    methods(free)={'none'};
else
    fallback=mode;
end
missing=strcmp(methods,'unavailable');
if any(missing)
    switch fallback
        case 'graph_bound',proxy=ddsync.graph_bound_std(n,a,b,w,pin,sigma);
        case 'pseudo_weight',proxy=pseudo_weight;
        otherwise,degree=accumarray([a;b],1,[n,1]);proxy=sigma./sqrt(degree);
    end
    out(missing)=proxy(missing);methods(missing)={fallback};
end
if ~strcmp(mode,'none')
    assert(all(isfinite(out) & out>=0 & isfinite(single(out))),'Uncertainty is outside the supported numeric range or the retained component is disconnected.');
    assert(sigma==0 || all(out(free)>0 & single(out(free))>0),'Uncertainty underflow; check weight/noise scaling.');
end
kinds=unique(methods(free));method='hybrid';if numel(kinds)==1,method=kinds{1};end
info=struct('method',method,'reason',reason,'probes',probes,'max_rel_se',max_rel_se, ...
    'hutch_converged',sum(strcmp(methods,'hutch')),'exact_nodes',sum(strcmp(methods,'exact')), ...
    'cholesky_converged',sum(strcmp(methods,'cholesky')),'cholesky_unconverged',sum(strcmp(methods,'cholesky_unconverged')), ...
    'repaired_nodes',repaired,'repair_attempted',attempted,'unconverged_nodes',sum(ismember(methods,{'hutch_unconverged','cholesky_unconverged'})), ...
    'fallback_nodes',sum(ismember(methods,{'pseudo_degree','pseudo_weight','graph_bound'})), ...
    'node_methods',{methods},'node_rel_se',relative);
end
