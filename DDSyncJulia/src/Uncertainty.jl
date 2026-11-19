# Per-event uncertainty decisions. No failed variance is clipped to a floor.
function validateStdConfig(cfg)
    lowercase(String(cfg[:mode])) in ("auto", "hutch", "cholesky", "exact", "graph_bound", "pseudo_degree", "pseudo_weight", "none") || error("Invalid std.mode")
    lowercase(String(get(cfg, :auto_probe, "hutch"))) in ("hutch", "cholesky") || error("std.auto_probe must be hutch or cholesky")
    lowercase(String(cfg[:fallback])) in ("exact", "graph_bound", "pseudo_degree", "pseudo_weight", "nan", "none") || error("Invalid std.fallback")
    lowercase(String(cfg[:unconverged])) in ("keep", "fallback") || error("Invalid std.unconverged")
    lowercase(String(cfg[:probe_dist])) in ("rademacher", "gaussian") || error("Invalid std.probe_dist")
    for key in (:probes, :min_probes, :batch, :max_nred, :exact_max_nred, :repair_max_nodes)
        value = cfg[key]
        (value isa Real && isfinite(value) && value == round(value) && value >= 0) || error("std.$key must be a nonnegative integer")
    end
    cfg[:batch] >= 1 || error("std.batch must be positive")
    cfg[:min_probes] >= 2 || error("std.min_probes must be at least 2")
    isfinite(cfg[:rel_tol]) && cfg[:rel_tol] > 0 || error("std.rel_tol must be finite and positive")
end

function uncertaintyLaplacian(n, a, b, w, free, ridge)
    valid = isfinite.(w) .& (w .> 0)
    aa=a[valid]; bb=b[valid]; ww=w[valid]
    W = sparse(vcat(aa,bb), vcat(bb,aa), vcat(ww,ww), n,n)
    L = spdiagm(0 => vec(sum(W,dims=2))) - W
    return L[free,free] + ridge * speye(length(free))
end

# Empirical errors, not simultaneous confidence bounds. Test each entry separately.
function hutchNodeDiagnostics(mu, m2, count, min_probes, rel_tol)
    rel = fill(Inf, length(mu))
    valid = isfinite.(mu) .& (mu .> 0) .& isfinite.(m2) .& (m2 .>= 0)
    if count >= max(2, min_probes)
        rel[valid] .= sqrt.(m2[valid] ./ ((count-1)*count)) ./ mu[valid]
    end
    return valid .& (rel .<= rel_tol), rel
end

function hutchDiagnostic(mu, m2, count, min_probes, rel_tol)
    good, rel = hutchNodeDiagnostics(mu,m2,count,min_probes,rel_tol)
    reason = any(!isfinite,mu) || any(!isfinite,m2) ? "nonfinite_estimate" :
        any(mu .<= 0) ? "nonpositive_estimate" : count < max(2,min_probes) ? "insufficient_probes" :
        all(good) ? "converged" : "not_converged"
    return reason, maximum(rel)
end

# Componentwise relative backward error (LAPACK BERR definition). Scaling
# only by RHS incorrectly rejects accurate solves with large cancelling terms.
function checkedSolveColumns(L,X,B)
    residual=abs.(L*X-B); scale=abs.(L)*abs.(X)+abs.(B)
    return vec(all(isfinite.(X) .& isfinite.(scale) .& (residual .<= 1e-12 .* max.(scale,(size(L,1)+1)*floatmin(Float64)/eps(Float64))),dims=1))
end

function hutchSamples(L,F,K,dist,batch,reportEvery,min_probes,rel_tol)
    n=size(L,1); mu=zeros(n); m2=zeros(n); count=0; reason="insufficient_probes"
    for bi in 1:cld(K,batch)
        k=min(batch,K-count)
        try
            Z=dist=="rademacher" ? ifelse.(rand(Bool,n,k),1.0,-1.0) : randn(n,k)
            X=F \ Z
            if !all(checkedSolveColumns(L,X,Z))
                reason="solve_failed"; break
            end
            samples=X.*Z; bm=vec(mean(samples,dims=2))
            bm2=vec(sum((samples.-bm).^2,dims=2))
            delta=bm.-mu; total=count+k
            m2 .+= bm2 .+ delta.^2 .* (count*k/total)
            mu .+= delta .* (k/total); count=total
        catch err
            err isa InterruptException && rethrow()
            reason="solve_failed"; break
        end
        reason,rel=hutchDiagnostic(mu,m2,count,min_probes,rel_tol)
        if reportEvery>0 && bi%reportEvery==0
            good,_=hutchNodeDiagnostics(mu,m2,count,min_probes,rel_tol)
            @printf("    Hutch diag: %d/%d probes, %d/%d entries converged, max relative MC SE=%.3g\n",count,K,sum(good),n,rel)
        end
        reason=="converged" && break
    end
    good,rel=hutchNodeDiagnostics(mu,m2,count,min_probes,rel_tol)
    return (mu=mu,m2=m2,converged=good,rel_se=rel,reason=reason,probes=count,max_rel_se=maximum(rel))
end

# Pinning one reference can leave a timing shift shared by many events.
# For C = inv(A), v = C*ones and s = sum(v), its covariance is v*v'/s.
# Draw one sample of C with a triangular solve, remove this shared response,
# and add its variance back exactly. Never form C or the dense outer product.
function choleskySamples(A, F, budget, dist, batch, reportEvery, min_probes, rel_tol)
    n = size(A, 1)
    mean_variance = zeros(n)
    m2 = zeros(n)
    count = 0
    reason = "insufficient_probes"
    try
        ones_rhs = ones(n, 1)
        common = F \ ones_rhs
        all(checkedSolveColumns(A, common, ones_rhs)) || throw(LinearAlgebra.SingularException(0))
        v = vec(common)
        s = sum(v)
        if !isfinite(s) || s <= 0
            throw(LinearAlgebra.PosDefException(0))
        end
        weight = v ./ s
        shared_variance = (v ./ sqrt(s)).^2
        # CHOLMOD stores A[p,p] = L*L'. UP restores the original event order.
        upper = sparse(transpose(sparse(F.L)))
        perm = F.p
        check_at = min(budget, max(256, min_probes))
        batches = 0
        while count < budget
            k = min(batch, check_at - count)
            z = dist == "rademacher" ? ifelse.(rand(Bool,n,k),1.0,-1.0) : randn(n,k)
            draws = F.UP \ z
            if !all(checkedSolveColumns(upper, draws[perm,:], z))
                reason = "solve_failed"
                break
            end
            for y in eachcol(draws)
                shift = sum(y)
                count += 1
                @inbounds for i in eachindex(y)
                    sample = abs2(y[i] - weight[i]*shift) + shared_variance[i]
                    delta = sample - mean_variance[i]
                    mean_variance[i] += delta / count
                    m2[i] += delta * (sample - mean_variance[i])
                end
            end
            batches += 1
            if count == check_at
                reason, _ = hutchDiagnostic(mean_variance, m2, count, min_probes, rel_tol)
                reason == "converged" && break
                check_at = count + min(count, budget-count)
            end
            if reportEvery > 0 && batches % reportEvery == 0
                @printf("    Cholesky sampling: %d/%d draws completed\n", count, budget)
            end
        end
    catch err
        err isa LinearAlgebra.SingularException || err isa LinearAlgebra.PosDefException || rethrow()
        reason = "solve_failed"
    end
    good, relative = hutchNodeDiagnostics(mean_variance, m2, count, min_probes, rel_tol)
    return (mu=mean_variance, m2=m2, converged=good, rel_se=relative,
            reason=reason, probes=count, max_rel_se=maximum(relative))
end

# Selected unit-vector solves against the FULL pinned system. A neighborhood
# inverse would condition away external weak links and is not an exact marginal.
# Reuse the Hutch factorization. Memory is bounded by nred * batch.
function selectedInverseDiagonal(L,F,indices,batch)
    diagonal=fill(NaN,length(indices)); nr=size(L,1)
    for first in 1:batch:length(indices)
        last=min(first+batch-1,length(indices)); ids=indices[first:last]
        try
            E=zeros(nr,length(ids))
            for (j,i) in enumerate(ids); E[i,j]=1.0; end
            X=F \ E; solve_ok=checkedSolveColumns(L,X,E)
            for (j,i) in enumerate(ids)
                if solve_ok[j] && isfinite(X[i,j]) && X[i,j]>0
                    diagonal[first+j-1]=X[i,j]
                end
            end
        catch err
            err isa InterruptException && rethrow()
        end
    end
    return diagonal
end

# Compatibility helper for studies of the raw estimator. Partial positive
# estimates are returned too; callers must inspect the per-entry diagnostics.
function estimateThetaStdHutch(n_c::Int,a::Vector{Int},b::Vector{Int},w::Vector{Float64},pin::Int,free::Vector{Int},ridge::Float64,
    sigma::Float64,minSigma::Float64,applyMinSigma::Bool,K::Int,dist::String,batch::Int,reportEvery::Int,minDiagRel::Float64;
    min_probes::Int=128,rel_tol::Float64=0.10)
    out=fill(NaN,n_c); out[pin]=0.0
    try
        L=uncertaintyLaplacian(n_c,a,b,w,free,ridge); F=cholesky(Symmetric(L))
        info=hutchSamples(L,F,K,dist,batch,reportEvery,min_probes,rel_tol)
        sig=sigmaForThetaStd(sigma,minSigma,applyMinSigma)
        valid=isfinite.(info.mu) .& (info.mu .> 0)
        out[free[valid]].=sig.*sqrt.(info.mu[valid])
        return out,info
    catch err
        err isa InterruptException && rethrow()
        return out,(reason="solve_failed",probes=0,max_rel_se=Inf)
    end
end

function estimateThetaStdExact(n,a,b,w,pin,free,ridge,sigma,minSigma,applyMinSigma,batch)
    out=fill(NaN,n); out[pin]=0.0
    try
        L=uncertaintyLaplacian(n,a,b,w,free,ridge); F=cholesky(Symmetric(L))
        diagonal=selectedInverseDiagonal(L,F,collect(eachindex(free)),batch)
        sig=sigmaForThetaStd(sigma,minSigma,applyMinSigma)
        all(isfinite,diagonal) && isfinite(sig) || return out,false
        values=sig.*sqrt.(diagonal)
        all(isfinite,values) || return out,false
        out[free].=values
        return out,true
    catch err
        err isa InterruptException && rethrow()
        return out,false
    end
end

function estimateThetaStd(n,a,b,w,pin,free,ridge,sigma,opts; pseudo_weight=fill(NaN,n))
    nr=length(free); out=fill(NaN,n); out[pin]=0.0
    methods=fill("unavailable",n); methods[pin]="reference"; rel_se=fill(NaN,n)
    reason="requested"; probes=0; max_rel_se=NaN; repaired=0; attempted=0
    mode=opts.mode=="auto" ? (nr<=opts.exact_max_nred ? "exact" : get(opts, :auto_probe, "hutch")) : opts.mode
    if mode!="none"
        isfinite(sigma) && sigma>=0 || error("Uncertainty needs a finite nonnegative noise scale; inspect residuals and std.min_sigma.")
    end
    fallback=opts.fallback in ("exact","nan","none") ? "pseudo_degree" : opts.fallback
    if mode in ("hutch","cholesky") && (opts.probes<=0 || nr>opts.max_nred) && opts.fallback=="exact" && nr<=opts.exact_max_nred
        reason=opts.probes<=0 ? "disabled_probes" : "size_limit"
        mode="exact"
    end
    if mode in ("exact","hutch","cholesky")
        if mode=="exact" && nr>opts.exact_max_nred
            reason="exact_size_limit"
        elseif mode in ("hutch","cholesky") && (opts.probes<=0 || nr>opts.max_nred)
            reason=opts.probes<=0 ? "disabled_probes" : "size_limit"
        else
            L=uncertaintyLaplacian(n,a,b,w,free,ridge)
            F=try cholesky(Symmetric(L)) catch err
                err isa InterruptException && rethrow()
                nothing
            end
            if F===nothing
                reason="solve_failed"
            elseif mode=="exact"
                diagonal=selectedInverseDiagonal(L,F,collect(1:nr),opts.batch)
                ok=isfinite.(diagonal)
                out[free[ok]].=sigma.*sqrt.(diagonal[ok]); methods[free[ok]].="exact"
                all(ok) || (reason="exact_failed")
            else
                sampler = mode == "cholesky" ? choleskySamples : hutchSamples
                sample=sampler(L,F,opts.probes,opts.probe_dist,opts.batch,opts.report_every_batch,opts.min_probes,opts.rel_tol)
                reason=sample.reason; probes=sample.probes; max_rel_se=sample.max_rel_se
                rel_se[free].=sample.rel_se
                good=sample.converged
                out[free[good]].=sigma.*sqrt.(sample.mu[good]); methods[free[good]].=mode
                bad=findall(.!good)
                # Invalid/near-zero estimates have infinite/large relative error
                # and are repaired first. Ties are deterministic by reduced index.
                sort!(bad,by=i->(-sample.rel_se[i],i))
                limit=opts.fallback=="exact" && nr<=opts.exact_max_nred ? length(bad) : opts.repair_max_nodes
                targets=bad[1:min(length(bad),limit)]; attempted=length(targets)
                exact=selectedInverseDiagonal(L,F,targets,opts.batch)
                for (j,i) in enumerate(targets)
                    if isfinite(exact[j])
                        out[free[i]]=sigma*sqrt(exact[j]); methods[free[i]]="exact"; repaired+=1
                    end
                end
                # Never retain negative/NaN estimates or a near-zero value whose
                # empirical SE is as large as the estimate. Those use fallback.
                if opts.unconverged=="keep" && probes>=opts.min_probes
                    for i in bad
                        if methods[free[i]]=="unavailable" && isfinite(sample.mu[i]) && sample.mu[i]>0 && sample.rel_se[i]<1.0
                            out[free[i]]=sigma*sqrt(sample.mu[i]); methods[free[i]]=mode * "_unconverged"
                        end
                    end
                end
            end
        end
    elseif mode=="none"
        methods[free].="none"
    else
        fallback=mode
    end
    missing=findall(==("unavailable"),methods)
    if !isempty(missing)
        proxy = fallback=="graph_bound" ? estimateThetaStdGraphBound(n,a,b,w,pin,sigma) :
            fallback=="pseudo_weight" ? pseudo_weight : sigma./sqrt.(nodeDegrees(n,a,b))
        out[missing].=proxy[missing]; methods[missing].=fallback
    end
    if mode!="none"
        # Fail explicitly instead of emitting Inf/NaN that a reader might floor.
        all(x->isfinite(x)&&x>=0&&isfinite(Float32(x)),out) || error("Uncertainty is outside the supported numeric range or the retained component is disconnected.")
        sigma==0 || all(x->x>0&&Float32(x)>0,out[free]) || error("Uncertainty underflow; check weight/noise scaling.")
    end
    kinds=unique(methods[free]); method=length(kinds)==1 ? only(kinds) : "hybrid"
    info=(method=method,reason=reason,probes=probes,max_rel_se=max_rel_se,
        hutch_converged=count(==("hutch"),methods),cholesky_converged=count(==("cholesky"),methods),
        cholesky_unconverged=count(==("cholesky_unconverged"),methods),exact_nodes=count(==("exact"),methods),
        repaired_nodes=repaired,repair_attempted=attempted,unconverged_nodes=count(m->m in ("hutch_unconverged","cholesky_unconverged"),methods),
        fallback_nodes=count(x->x in ("pseudo_degree","pseudo_weight","graph_bound"),methods),
        node_methods=methods,node_rel_se=rel_se)
    return out,info
end

"""
A finite upper bound on each pinned inverse-Laplacian diagonal.
A BFS spanning tree supplies a retained path to the reference. Its resistance
sum(1/w) bounds effective resistance from above (Rayleigh monotonicity).
A nonnegative ridge only decreases the inverse diagonal. O(nodes + edges).
"""
function estimateThetaStdGraphBound(n, a, b, w, pin, sigma)
    isfinite(sigma) && sigma > 0 || error("Graph-bound fallback needs a positive finite noise scale; check inputs or enable std.min_sigma.")
    adjacency = [Tuple{Int,Float64}[] for _ in 1:n]
    for k in eachindex(a)
        if isfinite(w[k]) && w[k] > 0
            resistance = 1.0 / w[k]
            push!(adjacency[a[k]], (b[k], resistance))
            push!(adjacency[b[k]], (a[k], resistance))
        end
    end
    distance = fill(Inf, n); distance[pin] = 0.0
    queue = Int[pin]
    head = 1
    while head <= length(queue)
        node = queue[head]; head += 1
        for (other, resistance) in adjacency[node]
            if isinf(distance[other])
                distance[other] = distance[node] + resistance
                isfinite(distance[other]) || error("Graph-bound fallback overflowed; check edge-weight scaling.")
                push!(queue, other)
            end
        end
    end
    all(isfinite, distance) || error("Graph-bound fallback requires a connected retained component.")
    out = sigma .* sqrt.(distance)
    all(x -> isfinite(x) && x <= floatmax(Float32), out) || error("Graph-bound uncertainty is outside the supported numeric range.")
    return out
end

function writeStdQuality(path,components)
    open(path,"w") do io
        println(io,"# event_id ref_id method relative_mc_se")
        for info in components, (i,id) in enumerate(info.event_ids)
            @printf(io,"%d %d %s %.9g\n",id,info.ref,info.node_methods[i],info.node_rel_se[i])
        end
    end
end
