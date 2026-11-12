function checks = per_event_tests()
%PER_EVENT_TESTS Independent covariance oracle and mixed-method regressions.
checks=0;n=300;a=[];b=[];w=[];
for j=2:260,edge(1,j,1);if j>2,edge(j-1,j,1);end,end
for i=261:299,for j=i+1:300,edge(i,j,1);end,end
edge(2,261,0.001);free=(2:n)';ridge=1e-10;
A=zeros(numel(a),n);
for k=1:numel(a),A(k,a(k))=1;A(k,b(k))=-1;end
L=A(:,free)'*diag(w)*A(:,free)+ridge*eye(n-1);expected=sqrt(diag(inv(L)));
cfg=ddsync.config_default();opts=cfg.std.hutch;opts.mode='hutch';opts.fallback='pseudo_degree';
opts.exact_max_nred=10;opts.repair_max_nodes=64;opts.unconverged='keep';opts.probes=512;opts.report_every_batch=0;
R=chol(sparse(L),'lower');rng(42);sample=ddsync.hutch_samples(sparse(L),R,opts);
check(sum(sample.converged)>250);check(sum(~sample.converged)>0 && sum(~sample.converged)<=64);
rng(42);[sd,info]=ddsync.estimate_theta_std(n,a,b,w,1,free,ridge,1,opts);
good=free(sample.converged);bad=free(~sample.converged);
check(max(abs(sd(good)-sqrt(sample.mu(sample.converged))))<1e-10);
check(all(strcmp(info.node_methods(good),'hutch')));check(all(strcmp(info.node_methods(bad),'exact')));
check(info.repaired_nodes==numel(bad) && info.fallback_nodes==0);
check(max(abs(sd(bad)./expected(bad-1)-1))<1e-8);check(strcmp(info.method,'hybrid'));
selected=ddsync.selected_inverse_diagonal(sparse(L),R,[299,2],2);
check(max(abs(sqrt(selected)./expected([299,2])-1))<1e-8);check(selected(1)>100);
partial=ddsync.selected_inverse_diagonal(eye(3),diag([1,sqrt(2),1]),[1,2,3],3);
check(all(partial([1,3])==1) && isnan(partial(2)));
chain=spdiags(repmat([-1,3,-1],900,1),[-1,0,1],900,900);
endpoint=ddsync.selected_inverse_diagonal(chain,chol(chain,'lower'),1,1);
check(abs(endpoint-(3-sqrt(5))/2)<1e-12);
[~,~,goodmask]=ddsync.hutch_diagnostic([-1;2],[1;0.5],128,128,0.1);check(isequal(goodmask,[false;true]));
opts.repair_max_nodes=0;rng(42);[kept,ki]=ddsync.estimate_theta_std(n,a,b,w,1,free,ridge,1,opts);
check(max(abs(kept(good)-sd(good)))<1e-12);check(ki.unconverged_nodes>0 && ki.repair_attempted==0);
for i=free'
    if strcmp(ki.node_methods{i},'hutch_unconverged')
        check(abs(kept(i)-sqrt(sample.mu(i-1)))<1e-8 && ki.node_rel_se(i)>0.1 && ki.node_rel_se(i)<1);
    end
end
opts.unconverged='fallback';rng(42);[proxy,pi]=ddsync.estimate_theta_std(n,a,b,w,1,free,ridge,1,opts);
check(max(abs(proxy(good)-sd(good)))<1e-12);check(pi.fallback_nodes==numel(bad) && pi.unconverged_nodes==0);
degree=accumarray([a;b],1,[n,1]);check(max(abs(proxy(bad)-1./sqrt(degree(bad))))<1e-12);
opts.unconverged='keep';opts.probes=1;opts.repair_max_nodes=2;rng(1);
[low,li]=ddsync.estimate_theta_std(n,a,b,w,1,free,ridge,1,opts);
check(li.exact_nodes==2 && li.fallback_nodes==n-3 && li.unconverged_nodes==0);check(all(isfinite(low)) && all(low(free)>0));
opts.mode='auto';opts.exact_max_nred=400;rng(10);next=rand;rng(10);
[exact,ei]=ddsync.estimate_theta_std(n,a,b,w,1,free,ridge,1,opts);
check(rand==next);check(strcmp(ei.method,'exact') && ei.probes==0 && ei.exact_nodes==n-1);
check(max(abs(exact(free)./expected-1))<1e-8);
opts.exact_max_nred=10;opts.probes=512;opts.repair_max_nodes=64;rng(42);
[large,ai]=ddsync.estimate_theta_std(n,a,b,w,1,free,ridge,1,opts);
check(isequal(large,sd) && strcmp(ai.method,'hybrid'));
opts.max_nred=0;[skip,si]=ddsync.estimate_theta_std(n,a,b,w,1,free,ridge,1,opts);
check(si.probes==0 && strcmp(si.method,'pseudo_degree') && all(isfinite(skip)));
opts.mode='hutch';opts.probes=0;opts.fallback='exact';opts.exact_max_nred=400;
[skipped,sk]=ddsync.estimate_theta_std(n,a,b,w,1,free,ridge,1,opts);
check(strcmp(sk.method,'exact') && strcmp(sk.reason,'disabled_probes') && sk.probes==0);
check(max(abs(skipped(free)./expected-1))<1e-8);
    function edge(i,j,v)
        a(end+1,1)=i;b(end+1,1)=j;w(end+1,1)=v;
    end
    function check(condition)
        assert(all(condition,'all'),'Per-event regression check %d failed',checks+1);checks=checks+1;
    end
end
