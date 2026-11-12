function valid = checked_solve_columns(L,X,B)
%CHECKED_SOLVE_COLUMNS Componentwise backward error (LAPACK BERR definition).
% Scaling only by RHS rejects accurate solves with large cancelling terms.
residual=abs(L*X-B);scale=abs(L)*abs(X)+abs(B);
valid=all(isfinite(X) & isfinite(scale) & residual<=1e-12*max(scale,(size(L,1)+1)*realmin/eps),1);
end
