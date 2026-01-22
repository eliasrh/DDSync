function out = merge_struct(a, b)
%MERGE_STRUCT Recursively merge two structs (b overrides a).
out = a;
if isempty(b), return; end
f = fieldnames(b);
for k=1:numel(f)
    key = f{k};
    if isstruct(b.(key))
        if isfield(out, key) && isstruct(out.(key))
            out.(key) = ddsync.merge_struct(out.(key), b.(key));
        else
            out.(key) = b.(key);
        end
    else
        out.(key) = b.(key);
    end
end
end
