function candidate = evaluateMethaneCandidate(model,qMagnitude,reactionAmplitude)
%EVALUATEMETHANECANDIDATE Refit production at fixed Darcy and reaction coefficients.
A = transportOperator(model,qMagnitude,reactionAmplitude);
% Adjoint sampling avoids forming the dense full inverse response.
response = (A'\model.Csample')'*model.M;
productionFit = fitProductionConcentration(model,response);
candidate = solveMethaneModel(model,qMagnitude,reactionAmplitude,productionFit.productionScaled);
candidate.productionFit = productionFit;
end
