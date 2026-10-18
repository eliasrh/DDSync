function out = graph_bound_std(n,a,b,w,pin,sigma)
%GRAPH_BOUND_STD Conservative finite uncertainty from a retained spanning tree.
% Path resistance sum(1/w) bounds pinned inverse-Laplacian diagonals above.
% A nonnegative ridge decreases those diagonals. Complexity: O(nodes+edges).
assert(isfinite(sigma) && sigma>0, ...
    'Graph-bound fallback needs a positive finite noise scale; check inputs or enable std.min_sigma.');
% Linked adjacency lists avoid repeatedly growing per-node arrays.
head=zeros(n,1);next=zeros(2*numel(a),1);to=next;resistance=next;count=0;
for k=1:numel(a)
    if ~isfinite(w(k)) || w(k)<=0,continue;end
    count=count+1;to(count)=b(k);resistance(count)=1/w(k);next(count)=head(a(k));head(a(k))=count;
    count=count+1;to(count)=a(k);resistance(count)=1/w(k);next(count)=head(b(k));head(b(k))=count;
end
% Traverse edges in input order, matching Julia's deterministic BFS tree.
for node=1:n
    current=head(node);previous=0;
    while current~=0
        following=next(current);next(current)=previous;previous=current;current=following;
    end
    head(node)=previous;
end
distance=Inf(n,1);distance(pin)=0;queue=zeros(n,1);queue(1)=pin;first=1;last=1;
while first<=last
    node=queue(first);first=first+1;edge=head(node);
    while edge~=0
        other=to(edge);
        if isinf(distance(other))
            distance(other)=distance(node)+resistance(edge);
            assert(isfinite(distance(other)),'Graph-bound fallback overflowed; check edge-weight scaling.');
            last=last+1;queue(last)=other;
        end
        edge=next(edge);
    end
end
assert(all(isfinite(distance)),'Graph-bound fallback requires a connected retained component.');
out=sigma*sqrt(distance);
assert(all(isfinite(out) & out<=double(realmax('single'))),'Graph-bound uncertainty is outside the supported numeric range.');
end
