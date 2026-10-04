function result = runMainModel(cfg)
%RUNMAINMODEL Fit all six main-model parameters, export results and figures.
root=fileparts(mfilename('fullpath'));
addpath(fullfile(root,'src'));
if nargin<1, cfg=defaultConfig(); end
assert(license('test','Optimization_Toolbox') && license('test','GADS_Toolbox'), ...
    'methane:Toolboxes','Optimization Toolbox and Global Optimization Toolbox are required.');
result=fitMethaneModel(cfg);
exportResults(result);
plotResults(result);
end
