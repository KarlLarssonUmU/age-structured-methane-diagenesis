% version which is smooth at z=0

% inputs:
% z     : vector of z-values to evaluate
% z_max : z-value for peak
% width : width of peak support
function y = productionShape2(z, z_max, width)
width = z_max + width;

% fz: piecewise polynomial with f^(6) = 0 on each interval.
% zm = z_max, zm+width in comments = width in code
%
% Interval 1: [0, zm], f^(6) = 0  -> quintic polynomial
%   left BC at z = 0:
%       f(0)      = 0
%       f^(4)(0)  = 0
%       f^(5)(0)  = 0
%   right BC at z = zm:
%       f(zm)     = 1
%       f'(zm)    = 0
%       f'''(zm)  = 0
%
% Interval 2: (zm, zm+width], f^(6) = 0  -> quintic polynomial
%   left BC at z = zm:
%       f(zm)           = 1
%       f'(zm)          = 0
%       f'''(zm)        = 0
%   right BC at z = zm+width:
%       f(zm+width)     = 0
%       f'(zm+width)    = 0
%       f''(zm+width)   = 0
%
% For z > zm+width (and z < 0): f(z) = 0.
%
% Inputs:
%   z      : evaluation points (scalar, vector, or array)
%   zm     : position of the maximum (must be > 0)
%   width  : distance from zm to the right endpoint (must be > 0)

z = reshape(z,1,[]);

    if z_max <= 0
        error('zm must be > 0.');
    end
    if width <= 0
        error('width must be > 0.');
    end

    z2 = width;

    % Initialize output (zero outside [0, z2])
    y = zeros(size(z));

    %% --- Interval 1: [0, zm] ---
    % Polynomial: p1(z) = a0 + a1 z + a2 z^2 + a3 z^3 + a4 z^4 + a5 z^5

    A1 = zeros(6,6);
    b1 = zeros(6,1);

    % 1) f(0) = 0
    A1(1,:) = drow(0, 0);
    b1(1)   = 0;

    % 2) f'(0) = 0
    A1(2,:) = drow(0, 1);
    b1(2)   = 0;

    % 3) f''(0) = 0
    A1(3,:) = drow(0, 2);
    b1(3)   = 0;

    % 4) f(zm) = 1
    A1(4,:) = drow(z_max, 0);
    b1(4)   = 1;

    % 5) f'(zm) = 0
    A1(5,:) = drow(z_max, 1);
    b1(5)   = 0;

    % 6) f'''(zm) = 0
    A1(6,:) = drow(z_max, 3);
    b1(6)   = 0;

    a1 = A1 \ b1;   % coefficients [a0; a1; ...; a5] for interval 1

    %% --- Interval 2: (zm, zm+width] ---
    % Polynomial: p2(z) = b0 + b1 z + ... + b5 z^5

    A2 = zeros(6,6);
    b2 = zeros(6,1);

    % Left BC at z = zm
    % 1) f(zm) = 1
    A2(1,:) = drow(z_max, 0);
    b2(1)   = 1;

    % 2) f'(zm) = 0
    A2(2,:) = drow(z_max, 1);
    b2(2)   = 0;

    % 3) f'''(zm) = 0
    A2(3,:) = drow(z_max, 3);
    b2(3)   = 0;

    % Right BC at z = z2 = zm + width
    % 4) f(z2) = 0
    A2(4,:) = drow(z2, 0);
    b2(4)   = 0;

    % 5) f'(z2) = 0
    A2(5,:) = drow(z2, 1);
    b2(5)   = 0;

    % 6) f''(z2) = 0
    A2(6,:) = drow(z2, 2);
    b2(6)   = 0;

    a2 = A2 \ b2;   % coefficients [b0; b1; ...; b5] for interval 2

    %% --- Evaluate piecewise ---

    y2 = y; % new basis

    % Interval 1: [0, zm]
    idx1 = (z >= 0) & (z <= z_max);
    if any(idx1)
        zz = z(idx1);
        y(idx1) = polyval_flip(a1, zz);
        y2(idx1) = ones(size(polyval_flip(a1, zz)));
    end

    % Interval 2: (zm, z2]
    idx2 = (z > z_max) & (z <= z2);
    if any(idx2)
        zz = z(idx2);
        y(idx2) = polyval_flip(a2, zz);
        y2(idx2) = polyval_flip(a2, zz);
    end

    % For z < 0 or z > z2, y remains zero from initialization.

    % Add shape for background production
    y = [y; ones(size(y))];

end

function row = drow(x, k)
% drow: row vector such that row * a = d^k/dx^k P(x),
% where P(x) = a0 + a1 x + ... + a5 x^5 and a = [a0; ...; a5].
% k = derivative order (0..5).

    row = zeros(1,6);
    for j = 0:5
        if j < k
            row(j+1) = 0;
        else
            % factor j*(j-1)*...*(j-k+1)
            coeff = 1;
            for m = 0:k-1
                coeff = coeff * (j - m);
            end
            row(j+1) = coeff * x^(j-k);
        end
    end
end

function y = polyval_flip(a, x)
% Evaluate P(x) = a0 + a1 x + ... + a5 x^5 for vector a (ascending powers).
    y = zeros(size(x));
    % Horner's method for ascending coefficients:
    % P(x) = (...((a5*x + a4)*x + a3)*x + a2)*x + a1)*x + a0
    for k = 5:-1:0
        y = y .* x + a(k+1);
    end
end