function validateConfig(cfg)
%VALIDATECONFIG Reject invalid units, bounds, or unsupported model variants.
assert(strcmp(cfg.ageScenario,'reservoirCorrected') && ...
    strcmp(cfg.fmTarget,'observations') && strcmp(cfg.productionFamily,'bump'), ...
    'methane:ModelChoice','Use the supplied reservoir-corrected main model.');
assert(strcmp(cfg.outerObjective,'Fm'),'methane:Objective', ...
    'Production minimizes E_C; the outer fit minimizes E_F.');
positive=[cfg.bottomDepth,cfg.meshStep,cfg.scaling,cfg.D0, ...
    cfg.shapeTolerance,cfg.outerTolerance];
assert(all(isfinite(positive) & positive>0),'methane:Config','Invalid positive setting.');
assert(isscalar(cfg.lossDepth) && isfinite(cfg.lossDepth) && ...
    cfg.lossDepth>0 && cfg.lossDepth<cfg.bottomDepth,'methane:LossDepth', ...
    'Loss depth must be within the domain.');
assert(isequal(cfg.deepLossMultiplier,0),'methane:LossDepth', ...
    'This model has no loss below z_loss.');
for name={'qBounds','kLossBounds'}
    v=cfg.(name{1});
    assert(isrow(v) && numel(v)==2 && all(isfinite(v)) && ...
        v(1)>=0 && v(2)>v(1),'methane:Bounds','Invalid parameter bounds.');
end
for name={'qGrid','kLossGrid'}
    v=cfg.(name{1});
    assert(isrow(v) && ~isempty(v) && all(isfinite(v) & v>=0), ...
        'methane:Grid','Invalid optimization grid.');
end
for name={'fitDarcy','fitLoss'}
    assert(islogical(cfg.(name{1})) && isscalar(cfg.(name{1})), ...
        'methane:Config','Fit switches must be scalar logicals.');
end
assert(isscalar(cfg.qMagnitude) && isfinite(cfg.qMagnitude) && cfg.qMagnitude>=0);
assert(isscalar(cfg.kLoss) && isfinite(cfg.kLoss) && cfg.kLoss>=0);
n=[cfg.shapeStarts,cfg.outerStarts,cfg.outerMaxIterations,cfg.concentrationFitPoints];
assert(all(isfinite(n) & n>=1 & n==round(n)) && cfg.shapeStarts<=15 && ...
    cfg.concentrationFitPoints>=2,'methane:Config','Invalid iteration/count setting.');
assert(ismember(cfg.boundaryMode,{'gradient','fixed'}),'methane:BoundaryMode', ...
    'Unknown surface boundary condition.');
if strcmp(cfg.boundaryMode,'fixed')
    assert(isscalar(cfg.fixedKsw) && isfinite(cfg.fixedKsw) && cfg.fixedKsw>0);
end
bumpShapeBounds(cfg);
end
