% version that includes optimization of the outer parameters

function [q_opt, E_tot, E_C, E_F, sM, VV, WW1, WW2] = ...
    optimizeProductionProfile(p, TRp, refUh, ...
                              sedDepth, sedFm, methConcDepth, methConcData, ...
                              metDepth, metFm, ...
                              betaF)
%OPTIMIZEPRODUCTIONPROFILE  Solve inverse problem for production profile.
%
% Fits the reduced production coefficients q >= 0 by nonlinear least squares:
%   E(q) = ||C - VV*q||^2 + betaF * ||F - (WW2*q)./(WW1*q)||^2
%
% Inputs are production-shape parameters and measurement data.
%
% Outputs
%   q_opt : optimized reduced production coefficients (length = size(sM,1))
%   E_tot : total objective = E_C + betaF*E_F
%   E_C   : concentration misfit = ||C - VV*q||^2
%   E_F   : Fm misfit           = ||F - (WW2*q)./(WW1*q)||^2


    % ---------------------------------------------------------------------
    % Outer optimization variables: x = [z_max, width]
    % ---------------------------------------------------------------------
    lb_outer = [0.0, 0.1];
    ub_outer = [0.6, 0.6];
    x0_outer = [0.2, 0.3];

    % Outer objective (calls inner solve)
    outerObj = @(x) outer_objective(x, p, TRp, refUh, ...
                                     sedDepth, sedFm, methConcDepth, methConcData, ...
                                     metDepth, metFm, betaF);

    opts_outer = optimoptions('patternsearch', ...
        'Display','final', ...
        'UseCompletePoll',true, ...
        'UseCompleteSearch',false, ...
        'MaxIterations', 50);

    [x_best, ~] = patternsearch(outerObj, x0_outer,[],[],[],[], ...
                                lb_outer, ub_outer, [], opts_outer);

    z_max = x_best(1);
    width = x_best(2);

    % ---------------------------------------------------------------------
    % Final inner solve at optimal outer parameters
    % ---------------------------------------------------------------------
    [q_opt, E_tot, E_C, E_F, sM, VV, WW1, WW2] = ...
        inner_solve(z_max, width, p, TRp, refUh, ...
                    sedDepth, sedFm, methConcDepth, methConcData, ...
                    metDepth, metFm, betaF);
end

function E = outer_objective(x, p, TRp, refUh, ...
                             sedDepth, sedFm, methConcDepth, methConcData, ...
                             metDepth, metFm, betaF)

    z_max = x(1);
    width = x(2);

    [~, E_tot] = inner_solve(z_max, width, p, TRp, refUh, ...
                             sedDepth, sedFm, methConcDepth, methConcData, ...
                             metDepth, metFm, betaF);

    E = E_tot;
end


function [q_opt, E_tot, E_C, E_F, sM, VV, WW1, WW2] = ...
    inner_solve(z_max, width, p, TRp, refUh, ...
                sedDepth, sedFm, methConcDepth, methConcData, ...
                metDepth, metFm, betaF)

    % Build production shape mapping
    sM = productionShape(p, z_max, width);

    guess0 = 0.2;
    nDof = size(refUh, 1);
    q_full0 = guess0 * ones(nDof, 1);
    q0 = sM * q_full0;

    V = sampleMatrixP1(methConcDepth, TRp);
    W = sampleMatrixP1(metDepth, TRp);

    VV  = V * refUh * sM';
    WW1 = W * refUh * sM';

    sedFmFcn = @(z) interp1(sedDepth, sedFm, z, 'linear', 'extrap');
    Fm_nodes = sedFmFcn(TRp);
    
    WW2 = W * refUh * diag(Fm_nodes) * sM';

    C = methConcData(:);
    F = metFm(:);

    lb = zeros(size(q0));
    ub = [];

    fun = @(q) residual_prod(q, VV, WW1, WW2, C, F, betaF);

    options = optimoptions('lsqnonlin', ...
        'Display','off', ...
        'FunctionTolerance',1e-20, ...
        'StepTolerance',1e-20, ...
        'OptimalityTolerance',1e-20, ...
        'MaxIterations',2000, ...
        'MaxFunctionEvaluations',2e5);

    q_opt = lsqnonlin(fun, q0, lb, ub, options);

    modelC = VV * q_opt;
    rC = C - modelC;
    E_C = rC' * rC;

    denom = WW1 * q_opt;
    num   = WW2 * q_opt;

    eps_denom = 1e-12;
    small = abs(denom) < eps_denom;
    denom(small) = eps_denom * sign(denom(small) + (denom(small) == 0));

    modelF = num ./ denom;
    rF = F - modelF;
    E_F = rF' * rF;

    E_tot = E_C + betaF * E_F;
end

function r = residual_prod(q, VV, WW1, WW2, C, F, betaF)
%RESIDUAL_PROD Residual vector so that lsqnonlin minimizes ||r||^2.

    modelC = VV * q;

    denom = WW1 * q;
    num   = WW2 * q;

    % Guard against division by near-zero
    eps_denom = 1e-12;
    small = abs(denom) < eps_denom;
    denom(small) = eps_denom * sign(denom(small) + (denom(small) == 0));

    modelF = num ./ denom;

    rC = C - modelC;
    rF = F - modelF;

    r = [rC; sqrt(betaF) * rF];
end