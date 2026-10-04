function fit = fitProductionConcentration(model,response)
%FITPRODUCTIONCONCENTRATION Optimize all production parameters using E_C.
cfg = model.cfg; y = model.data.concentration;
objective = @(x) shapeError(x);
[lb,ub] = bumpShapeBounds(cfg);
[a,b] = ndgrid([0.03,0.12,0.23,0.40,0.58],[0.12,0.30,0.55]);
starts = min(max([a(:),b(:)],lb),ub);
values = zeros(size(starts,1),1);
for i=1:size(starts,1), values(i)=objective(starts(i,:)); end
[~,order] = sort(values);
options = optimoptions('fmincon','Algorithm','sqp','Display','off', ...
    'OptimalityTolerance',cfg.shapeTolerance,'StepTolerance',1e-10, ...
    'MaxIterations',150,'MaxFunctionEvaluations',1500,'FiniteDifferenceType','central');
trialShape = starts(order(1),:); bestError = values(order(1)); bestFlag = 0;
fits = zeros(min(cfg.shapeStarts,numel(order)),4);
for i=1:size(fits,1)
    [shape,err,flag] = fmincon(objective,starts(order(i),:),[],[],[],[],lb,ub,[],options);
    fits(i,:) = [shape,err,flag];
    if err<bestError || (err<=bestError+1e-12 && bestFlag==0)
        trialShape=shape; bestError=err; bestFlag=flag;
    end
end
fit.shape = trialShape; fit.exitflag = bestFlag; fit.startResults = fits;
fit.effectiveShapeBounds = [lb;ub];
fit.shapeAtBound = any(abs(fit.shape-lb)<1e-5 | abs(fit.shape-ub)<1e-5);
fit.basis = productionBasis(model.z,cfg.productionFamily,fit.shape);
[fit.amplitudesScaled,fit.E_C] = nonnegativeAmplitudes(response*fit.basis,y);
fit.productionScaled = fit.basis*fit.amplitudesScaled;

    function err = shapeError(shape)
        B = productionBasis(model.z,cfg.productionFamily,shape);
        [~,err] = nonnegativeAmplitudes(response*B,y);
    end
end
