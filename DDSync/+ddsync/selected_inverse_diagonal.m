function diagonal = selected_inverse_diagonal(L,R,indices,batch)
%SELECTED_INVERSE_DIAGONAL Unit RHS solves against the FULL pinned system.
% Reuses the Hutch factorization. A neighborhood inverse is not an exact marginal.
diagonal=NaN(numel(indices),1);nr=size(L,1);
for first=1:batch:numel(indices)
    last=min(first+batch-1,numel(indices));ids=indices(first:last);
    try
        E=zeros(nr,numel(ids));
        for j=1:numel(ids),E(ids(j),j)=1;end
        X=R'\(R\E);solve_ok=ddsync.checked_solve_columns(L,X,E);
        for j=1:numel(ids)
            if solve_ok(j) && isfinite(X(ids(j),j)) && X(ids(j),j)>0
                diagonal(first+j-1)=X(ids(j),j);
            end
        end
    catch
        % Failed entries remain missing and are handled individually by caller.
    end
end
end
