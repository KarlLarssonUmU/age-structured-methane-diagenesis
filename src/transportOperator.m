function [A,kSw] = transportOperator(model,qMagnitude,reactionAmplitude)
%TRANSPORTOPERATOR qMagnitude >= 0 denotes upward Darcy flow in m/s.
% Third argument: the effective loss coefficient kLoss [day^-1].
assert(qMagnitude>=0 && isfinite(qMagnitude) && ...
    reactionAmplitude>=0 && isfinite(reactionAmplitude));
cfg = model.cfg;
switch cfg.boundaryMode
    case 'gradient', kSw = model.surfaceK0+qMagnitude;
    case 'fixed', kSw = cfg.fixedKsw;
    otherwise, error('methane:BoundaryMode','Unknown boundary mode.');
end
assert(isfinite(kSw) && kSw>0);
A = model.D-qMagnitude*cfg.scaling*model.V+reactionAmplitude*model.R ...
    +kSw*cfg.scaling*model.surfaceMatrix;
end
