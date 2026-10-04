function result = fitMethaneModel(cfg,data)
%FITMETHANEMODEL Concentration-only production inside an outer Fm search.
if nargin<1, cfg=defaultConfig(); end
if nargin<2, data=loadModelData(cfg); end
model=prepareModel(cfg,data); timer=tic;
rate=model.rateSettings;
history=zeros(0,5); cache=containers.Map('KeyType','char','ValueType','any');
if cfg.fitDarcy, qs=unique([cfg.qBounds,cfg.qGrid]); else, qs=cfg.qMagnitude; end
if rate.fit, rs=unique([rate.bounds,rate.grid]); else, rs=rate.value; end
if cfg.fitDarcy, qs=qs(qs>=cfg.qBounds(1) & qs<=cfg.qBounds(2)); end
if rate.fit, rs=rs(rs>=rate.bounds(1) & rs<=rate.bounds(2)); end
scores=zeros(numel(qs),numel(rs));
for j=1:numel(rs)
    for i=1:numel(qs), scores(i,j)=objective([qs(i)/1e-9,rs(j)/rate.scale]); end
end
[bestScore,index]=min(scores(:)); [i,j]=ind2sub(size(scores),index);
bestX=[qs(i)/1e-9,rs(j)/rate.scale]; outerFlags=[];
if cfg.fitDarcy && rate.fit
    [~,order]=sort(scores(:)); seeds=zeros(0,2);
    for k=order'
        [i,j]=ind2sub(size(scores),k); x=[qs(i)/1e-9,rs(j)/rate.scale];
        if isempty(seeds) || all(vecnorm(seeds-x,2,2)>0.4)
            seeds(end+1,:)=x; %#ok<AGROW>
        end
        if size(seeds,1)>=cfg.outerStarts, break; end
    end
    opts=optimoptions('patternsearch','Display','off','UseCompletePoll',true, ...
        'MeshTolerance',cfg.outerTolerance,'MaxIterations',cfg.outerMaxIterations, ...
        'MaxFunctionEvaluations',1500,'FunctionTolerance',1e-10);
    for s=1:size(seeds,1)
        [x,e,flag]=patternsearch(@objective,seeds(s,:),[],[],[],[], ...
            [cfg.qBounds(1)/1e-9,rate.bounds(1)/rate.scale], ...
            [cfg.qBounds(2)/1e-9,rate.bounds(2)/rate.scale],[],opts);
        outerFlags(end+1)=flag; %#ok<AGROW>
        accept(x,e);
    end
    % Optimize every boundary, including exactly zero flow and loss.
    for q=cfg.qBounds
        [x,e]=refine1D(rs/rate.scale,@(r) objective([q/1e-9,r])); accept([q/1e-9,x],e);
    end
    for r=rate.bounds
        [x,e]=refine1D(qs/1e-9,@(q) objective([q,r/rate.scale])); accept([x,r/rate.scale],e);
    end
elseif cfg.fitDarcy
    [x,e]=refine1D(qs/1e-9,@(q) objective([q,rate.value/rate.scale])); accept([x,rate.value/rate.scale],e);
elseif rate.fit
    [x,e]=refine1D(rs/rate.scale,@(r) objective([cfg.qMagnitude/1e-9,r])); accept([cfg.qMagnitude/1e-9,x],e);
end
candidate=cache(key(bestX));
result= candidate;
result.config=cfg; result.data=data; result.z=model.z;
result.sourceFm=model.sourceFm; result.sourceAge=model.sourceAge;
result.diagnostics=computeDiagnostics(model,candidate);
result.history=array2table(history,'VariableNames',{'qMagnitude_m_s',rate.historyName,'E_C','E_F','innerExitflag'});
result.outerFlags=outerFlags; result.outerScore=bestScore;
result.qAtUpperBound=cfg.fitDarcy && abs(result.qMagnitude-cfg.qBounds(2))<1e-13;
result.rateAtUpperBound=rate.fit && abs(result.reactionAmplitude-rate.bounds(2))<rate.scale*2e-4;
result.kLossAtUpperBound=result.rateAtUpperBound;
result.elapsedSeconds=toc(timer); result.matlabVersion=version;
if result.productionFit.exitflag<=0 || any(result.outerFlags<=0)
    warning('methane:OptimizationNotConverged', ...
        ['The selected production fit or an outer optimization search stopped ' ...
         'without convergence. Treat the results as provisional; review ' ...
         'result.productionFit.exitflag and result.outerFlags, and increase ' ...
         'the relevant optimization limits before interpreting the fit.']);
end
if cfg.verbose
    fprintf('%s / %s / %s: q=%.5g m/s %s=%.5g E_C=%.6g E_F=%.6g dissolved export=%.5g mmol/m2/day (%.1fs)\n', ...
        cfg.ageScenario,cfg.productionFamily,cfg.fmTarget,result.qMagnitude,rate.historyName,result.reactionAmplitude, ...
        result.E_C,result.E_F,result.diagnostics.upwardFlux*86400*1000,result.elapsedSeconds);
end

    function value=objective(x)
        cacheKey=key(x);
        if isKey(cache,cacheKey), c=cache(cacheKey);
        else
            c=evaluateMethaneCandidate(model,x(1)*1e-9,x(2)*rate.scale);
            cache(cacheKey)=c;
            history(end+1,:)=[c.qMagnitude,c.reactionAmplitude,c.E_C,c.E_F,c.productionFit.exitflag];
        end
        switch cfg.outerObjective
            case 'Fm', value=c.E_F;
            case 'concentration', value=c.E_C;
            otherwise, error('methane:OuterObjective','Unknown objective.');
        end
    end
    function accept(x,e)
        if e<bestScore, bestScore=e; bestX=x; end
    end
    function k=key(x)
        k=sprintf('%.14g,%.14g',x(1),x(2));
    end
    function [xbest,ebest]=refine1D(grid,fun)
        e=arrayfun(fun,grid); [ebest,bestGridIndex]=min(e); xbest=grid(bestGridIndex);
        opts=optimset('Display','off','TolX',cfg.outerTolerance);
        for n=2:numel(grid)-1
            if e(n)<=e(n-1) && e(n)<=e(n+1)
                [x,v,flag]=fminbnd(fun,grid(n-1),grid(n+1),opts);
                outerFlags(end+1)=flag;
                if v<ebest, xbest=x; ebest=v; end
            end
        end
        % Refine adjacent to the best endpoint as well.
        if numel(grid)>1 && (bestGridIndex==1 || bestGridIndex==numel(grid))
            if bestGridIndex==1, bracket=grid(1:2); else, bracket=grid(end-1:end); end
            [x,v,flag]=fminbnd(fun,bracket(1),bracket(2),opts);
            outerFlags(end+1)=flag;
            if v<ebest, xbest=x; ebest=v; end
        end
    end
end
