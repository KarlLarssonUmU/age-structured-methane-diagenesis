function d = computeDiagnostics(model,solution)
%COMPUTEDIAGNOSTICS Conservative budgets and normalized nodal-source fluxes.
cfg = model.cfg; z = model.z; P = solution.productionScaled;
rate=reactionSettings(cfg);
amplitude=solution.reactionAmplitude;
[A,kSw] = transportOperator(model,solution.qMagnitude,amplitude);
one = ones(size(z));
topResponse = full(model.M'*(A'\model.surfaceSample'));
topContributions = topResponse.*P;
assert(min(topContributions)>-1e-10*sum(topContributions), ...
    'methane:NegativeClassFlux','Mesh does not preserve nonnegative class fluxes.');
topContributions = max(topContributions,0); % only roundoff-level negatives
d.sourceFraction = topContributions/sum(topContributions);
d.sourceFlux = kSw*topContributions;
d.upwardFlux = kSw*solution.concentration(1);
d.grossProduction = one'*model.M*P/cfg.scaling;
d.loss = amplitude*one'*model.R*solution.concentration/cfg.scaling;
d.massBalanceRelative = abs(d.upwardFlux-d.grossProduction+d.loss)/max(d.grossProduction,eps);
d.superpositionRelative = abs(sum(topContributions)-solution.concentration(1))/solution.concentration(1);
d.fmSuperpositionError = abs(sum(d.sourceFraction.*model.sourceFm)-solution.fm(1));
d.advectiveUpwardFlux = solution.qMagnitude*solution.concentration(1);
% Conservative diffusive component obtained from total minus advective flux.
d.diffusiveUpwardFlux = d.upwardFlux-d.advectiveUpwardFlux;
d.gradientDiffusiveUpwardFlux = model.coefficient(1)* ...
    (solution.concentration(2)-solution.concentration(1))/(z(2)-z(1));
d.gradientFluxRelativeError = abs(d.gradientDiffusiveUpwardFlux+d.advectiveUpwardFlux-d.upwardFlux)/d.upwardFlux;
d.lossFraction = d.loss/d.grossProduction;
% Attribute production and loss to each nodal source by conservation.
d.sourceGrossProduction=P.*full(model.M'*one)/cfg.scaling;
d.sourceLoss=d.sourceGrossProduction-d.sourceFlux;
assert(min(d.sourceLoss)>-1e-10*max(d.grossProduction,eps));
d.sourceBudgetRelative=abs(sum(d.sourceLoss)-d.loss)/max(d.grossProduction,eps);
d.cumulativeDepthFraction = cumsum(d.sourceFraction);
[d.uniqueSourceAge,~,group] = unique(model.sourceAge);
d.ageFraction = accumarray(group,d.sourceFraction);
d.cumulativeAgeFraction = cumsum(d.ageFraction);
d.ageProbabilities = [0.10,0.25,0.50,0.75,0.90];
d.ageQuantiles = arrayfun(@(p) weightedQuantile(d.uniqueSourceAge,d.cumulativeAgeFraction,p),d.ageProbabilities);
% Thresholds are conventional apparent radiocarbon ages, not storage times.
d.fractionApparentAgeOver120 = sum(d.sourceFraction(model.sourceAge>120));
d.fractionApparentAgeOver2300 = sum(d.sourceFraction(model.sourceAge>2300));
d.legacyFraction = sum(d.sourceFraction(z>0.10+1e-12));
d.legacyExport = sum(d.sourceFlux(z>0.10+1e-12));
d.modernFraction = sum(d.sourceFraction(model.sourceAge<=0));
d.maxCellPe = max(solution.qMagnitude*diff(z)./(2*model.midpointCoefficient));
d.minimumConcentration = min(solution.concentration);
d.surfaceFm = solution.fm(1);
d.productionPeak = max(P)/cfg.scaling;
[~,i] = max(P); d.productionPeakDepth = z(i);
d.rmseConcentration = sqrt(solution.E_C/numel(model.data.concentration));
d.rmseFm = sqrt(solution.E_F/numel(model.data.fm));
d.rmseObservedFm = sqrt(solution.E_F_observations/numel(model.data.observedFm));
d.lossProfile = amplitude*rate.perDayFactor/86400* ...
    ((z<=rate.depth)+(z>rate.depth)*rate.deepMultiplier).*solution.concentration;
end

function value = weightedQuantile(age,cumulative,p)
i = find(cumulative>=p,1,'first');
if isempty(i), value=age(end); else, value=age(i); end
end
