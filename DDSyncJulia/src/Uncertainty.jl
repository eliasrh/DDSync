# Numerical uncertainty diagnostics. Failed estimates are never clipped to a floor.
function validateStdConfig(cfg)
    lowercase(String(cfg[:mode])) in ("hutch", "exact", "graph_bound", "pseudo_degree", "pseudo_weight", "none") || error("Invalid std.mode")
    lowercase(String(cfg[:fallback])) in ("exact", "graph_bound", "pseudo_degree", "pseudo_weight", "nan", "none") || error("Invalid std.fallback")
    lowercase(String(cfg[:probe_dist])) in ("rademacher", "gaussian") || error("Invalid std.probe_dist")
    for key in (:probes, :min_probes, :batch, :max_nred, :exact_max_nred)
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

# Per-node empirical Monte Carlo standard error. This diagnoses sampling noise,
# not model error, and is not a simultaneous confidence bound.
function hutchDiagnostic(mu, m2, count, min_probes, rel_tol)
    if any(x -> !isfinite(x), mu) || any(x -> !isfinite(x), m2)
        return "nonfinite_estimate", Inf
    elseif any(mu .<= 0)
        return "nonpositive_estimate", Inf
    elseif count < max(2, min_probes)
        return "insufficient_probes", Inf
    end
    se = sqrt.(max.(m2, 0.0) ./ ((count-1) * count))
    rel = maximum(se ./ mu)
    return rel <= rel_tol ? "converged" : "not_converged", rel
end

function estimateThetaStdHutch(n_c::Int, a::Vector{Int}, b::Vector{Int}, w_eff::Vector{Float64},
    pin::Int, keep_idx::Vector{Int}, ridge::Float64,
    sigma_hat::Float64, minSigma::Float64, applyMinSigma::Bool, K::Int, dist::String,
    batch::Int, reportEvery::Int, minDiagRel::Float64; min_probes::Int=128, rel_tol::Float64=0.10)
    # minDiagRel is retained for old callers only; no diagonal clipping is performed.
    out = fill(NaN,n_c); out[pin]=0.0
    count = 0; rel = Inf
    try
        L = uncertaintyLaplacian(n_c,a,b,w_eff,keep_idx,ridge)
        F = cholesky(Symmetric(L))
        n = length(keep_idx)
        mu = zeros(n); m2 = zeros(n)
        reason = "insufficient_probes"
        for bi in 1:cld(K,batch)
            k = min(batch,K-count)
            Z = dist == "rademacher" ? ifelse.(rand(Bool,n,k),1.0,-1.0) : randn(n,k)
            X = F \ Z
            # A failed/ill-conditioned solve must not look like a valid random sample.
            residual = L * X - Z
            if !all(isfinite,X) || any(vec(maximum(abs.(residual),dims=1)) .> 1e-8 .* max.(1.0,vec(maximum(abs.(Z),dims=1))))
                return out, (reason="solve_failed", probes=count, max_rel_se=Inf)
            end
            samples = X .* Z
            batch_mu = vec(mean(samples,dims=2))
            batch_m2 = vec(sum((samples .- batch_mu).^2,dims=2))
            delta = batch_mu .- mu
            total = count + k
            m2 .+= batch_m2 .+ delta.^2 .* (count*k/total)
            mu .+= delta .* (k/total)
            count = total
            reason, rel = hutchDiagnostic(mu,m2,count,min_probes,rel_tol)
            if reportEvery > 0 && bi % reportEvery == 0
                @printf("    Hutch diag: %d/%d probes, max relative MC SE=%.3g (%s)\n", count,K,rel,reason)
            end
            reason == "converged" && break
        end
        if reason == "converged"
            sig = sigmaForThetaStd(sigma_hat,minSigma,applyMinSigma)
            if !isfinite(sig)
                reason = "invalid_noise_scale"
            else
                out[keep_idx] .= sig .* sqrt.(mu)
            end
        end
        return out, (reason=reason, probes=count, max_rel_se=rel)
    catch err
        err isa InterruptException && rethrow()
        return out, (reason="solve_failed", probes=count, max_rel_se=rel)
    end
end

function estimateThetaStdExact(n,a,b,w,pin,free,ridge,sigma,minSigma,applyMinSigma,batch)
    out = fill(NaN,n); out[pin]=0.0
    try
        L = uncertaintyLaplacian(n,a,b,w,free,ridge)
        F = cholesky(Symmetric(L))
        nr = length(free); diagonal = zeros(nr)
        # Unit-vector solves in bounded batches; do not form a dense inverse.
        for first in 1:batch:nr
            last = min(first+batch-1,nr)
            E = zeros(nr,last-first+1)
            for (j,i) in enumerate(first:last)
                E[i,j] = 1.0
            end
            X = F \ E
            if !all(isfinite,X) || maximum(abs.(L*X-E)) > 1e-8
                return out,false
            end
            for (j,i) in enumerate(first:last)
                diagonal[i] = X[i,j]
            end
        end
        sig = sigmaForThetaStd(sigma,minSigma,applyMinSigma)
        if !isfinite(sig) || any(x -> !isfinite(x) || x <= 0, diagonal)
            return out,false
        end
        out[free] .= sig .* sqrt.(diagonal)
        return out,true
    catch err
        err isa InterruptException && rethrow()
        return out,false
    end
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
