% run_ddsync.m  -- DDSync standalone runner
%
% Edit ddsync_params.m (optional) to customize settings.
% If ddsync_params.m is missing, defaults are used and you must set cfg.io paths.

this_dir = fileparts(mfilename('fullpath'));
addpath(this_dir);  % adds +ddsync package

if exist(fullfile(this_dir,'ddsync_params.m'),'file')
    cfg = ddsync_params();
else
    cfg = ddsync.config_default();
    fprintf('[DDSync] Using config_default(). Tip: copy ddsync_params_template.m -> ddsync_params.m and edit.\n');
end

ddsync.run(cfg);
