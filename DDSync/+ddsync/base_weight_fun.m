function fh = base_weight_fun(mode)
%BASE_WEIGHT_FUN Return a function handle mapping CC -> base weight.
% mode:
%   'cc'   : w = cc
%   'cc2'  : w = cc.^2
%   'ones' : w = 1
if isa(mode,'function_handle')
    fh = mode;
    return;
end
if isempty(mode), mode = 'cc'; end
mode = lower(string(mode));
switch mode
    case "cc"
        fh = @(cc) cc;
    case "cc2"
        fh = @(cc) cc.^2;
    case "ones"
        fh = @(cc) ones(size(cc));
    otherwise
        error('Unknown cfg.weights.base_fun mode: %s', mode);
end
end
