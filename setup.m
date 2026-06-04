function setup()
%SETUP  Add Mixture-LDP source folders to the path for this MATLAB session.
%   Run once after starting MATLAB (from the repository root):
%       >> setup
%   Afterwards every function in src/ and every script in experiments/ is
%   callable by name. Clearing the workspace (clear/clear all) does not
%   remove the path, so a single call per session is enough.

    here = fileparts(mfilename('fullpath'));
    addpath(genpath(fullfile(here, 'src')));
    addpath(genpath(fullfile(here, 'experiments')));
    fprintf('Mixture-LDP paths added. You can now run any experiments/*.m script.\n');
end
