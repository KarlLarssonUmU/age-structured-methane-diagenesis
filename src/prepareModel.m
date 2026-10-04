function model = prepareModel(cfg,data)
%PREPAREMODEL Assemble fixed P1 matrices and sampling maps once.
validateConfig(cfg);
rate=reactionSettings(cfg);
bumpShapeBounds(cfg);
N = round(cfg.bottomDepth/cfg.meshStep); h = cfg.bottomDepth/N;
assert(abs(h-cfg.meshStep)<1e-10,'Mesh step must divide the domain.');
assert(abs(rate.depth/h-round(rate.depth/h))<1e-8, ...
    'The loss discontinuity must lie on a mesh node.');
coeffData = cfg.scaling*cfg.D0*data.porosity./(1-log(data.porosity.^2));
coefficient = @(z) interp1(data.porosityDepth,coeffData,z,'linear');
zero = @(z) 0;
[D,M,~,z] = assemble1dAdvectionDiffusion(0,cfg.bottomDepth,N,'p1',coefficient,zero,zero,zero);
V = assemble1dAdvectionDiffusion(0,cfg.bottomDepth,N,'p1',zero,zero,zero,@(z) 1);
reaction = @(z) cfg.scaling*rate.perDayFactor/86400* ...
    ((z<=rate.depth)+(z>rate.depth)*rate.deepMultiplier);
R = assemble1dAdvectionDiffusion(0,cfg.bottomDepth,N,'p1',zero,reaction,zero,zero);
model.cfg = cfg; model.data = data; model.z = z;
model.rateSettings=rate;
model.D = D; model.M = M; model.V = V; model.R = R;
model.surfaceMatrix = sparse(1,1,1,N+1,N+1);
model.surfaceSample = sparse(1,1,1,1,N+1);
model.Csample = sampleMatrixP1(data.concentrationDepth,z);
model.Fsample = sampleMatrixP1(data.fmDepth,z);
model.observedFsample = sampleMatrixP1(data.observedFmDepth,z);
model.sourceFm = interp1(data.sedimentDepth,data.sedimentFm,z,'linear');
model.sourceAge = interp1(data.sedimentDepth,data.sedimentAge,z,'linear');
model.coefficient = coefficient(z)/cfg.scaling;
model.porosity = interp1(data.porosityDepth,data.porosity,z,'linear');
model.surfaceK0 = coefficient(0)/cfg.scaling*data.surfaceLine(1)/data.surfaceLine(2);
model.midpointCoefficient = coefficient((z(1:end-1)+z(2:end))/2)/cfg.scaling;
end
