function solution = solveMethaneModel(model,qMagnitude,reactionAmplitude,productionScaled)
%SOLVEMETHANEMODEL Two linear solves: total methane and Fm-weighted methane.
[A,kSw] = transportOperator(model,qMagnitude,reactionAmplitude);
rhs = model.M*[productionScaled,model.sourceFm.*productionScaled];
u = A\rhs;
solution.qMagnitude = qMagnitude;
solution.qDarcy = -qMagnitude;
rate=reactionSettings(model.cfg);
solution.reactionAmplitude=reactionAmplitude;
solution.kLoss_per_day=reactionAmplitude*rate.perDayFactor;
solution.kSw = kSw;
solution.productionScaled = productionScaled;
solution.concentration = u(:,1);
solution.fmConcentration = u(:,2);
assert(all(u(:,1)>0),'methane:Nonpositive','Nonpositive methane solution.');
solution.fm = u(:,2)./u(:,1);
solution.predictedC = model.Csample*u(:,1);
solution.predictedF = (model.Fsample*u(:,2))./(model.Fsample*u(:,1));
solution.predictedObservedF = (model.observedFsample*u(:,2))./(model.observedFsample*u(:,1));
solution.E_C = sum((solution.predictedC-model.data.concentration).^2);
solution.E_F = sum((solution.predictedF-model.data.fm).^2);
solution.E_F_observations = sum((solution.predictedObservedF-model.data.observedFm).^2);
solution.linearResidual = norm(A*u-rhs,'fro')/max(norm(rhs,'fro'),eps);
end
