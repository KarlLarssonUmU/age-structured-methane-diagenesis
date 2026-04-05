% defining list of elements:
% basis polynomials and coordinates in reference domain [0,1]

function [coeff1d,coords,order] = polyBasisTypes1d(type)

switch type
    case {'p0'}
        order = 0; coords = 1/2;
        coeff1d = 1;
    case {'p1','h1'}
        order = 1; coords = linspace(0,1,order+1);
        coeff1d = [-1,1;1,0];
    case {'p2'}
        order = 2; coords = linspace(0,1,order+1);
        coeff1d = [2,-3,1;-4,4,0;2,-1,0];
    otherwise
        error(strcat(type,' polynomial not implemented.'))
end
coords=coords';
end