function V = sampleMatrixP1(zq, TRp)
%SAMPLEMATRIXP1  Build a sampling matrix for P1 (piecewise linear) 1D FE.
%
%   V*u gives the values of the nodal FE function u at query depths zq by
%   linear interpolation on the mesh nodes TRp. No assumption on zq spacing.
%
% Inputs
%   zq  : query depths (scalar or vector), in meters
%   TRp : nodal coordinates (N+1 x 1), in meters, increasing
%
% Output
%   V   : sparse matrix of size (length(zq) x length(TRp))

    nodes = TRp(:);
    nNodes = numel(nodes);

    zq = zq(:);
    nQ = numel(zq);

    % Preallocate (each row has at most 2 nonzeros)
    I = zeros(2*nQ, 1);
    J = zeros(2*nQ, 1);
    S = zeros(2*nQ, 1);
    nnzCount = 0;

    for r = 1:nQ
        z = zq(r);

        % Clamp outside domain to nearest endpoint node
        if z <= nodes(1)
            nnzCount = nnzCount + 1;
            I(nnzCount) = r; J(nnzCount) = 1; S(nnzCount) = 1;
            continue;
        elseif z >= nodes(end)
            nnzCount = nnzCount + 1;
            I(nnzCount) = r; J(nnzCount) = nNodes; S(nnzCount) = 1;
            continue;
        end

        % Find element [i, i+1] such that nodes(i) <= z < nodes(i+1)
        i = find(nodes <= z, 1, 'last');
        if i >= nNodes
            i = nNodes - 1;
        end

        x1 = nodes(i);
        x2 = nodes(i+1);
        h  = x2 - x1;

        % Linear interpolation weights
        w2 = (z - x1) / h;
        w1 = 1 - w2;

        nnzCount = nnzCount + 1;
        I(nnzCount) = r; J(nnzCount) = i;   S(nnzCount) = w1;

        nnzCount = nnzCount + 1;
        I(nnzCount) = r; J(nnzCount) = i+1; S(nnzCount) = w2;
    end

    % Build sparse matrix
    V = sparse(I(1:nnzCount), J(1:nnzCount), S(1:nnzCount), nQ, nNodes);
end