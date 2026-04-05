% s - local coordinate in [0,1]
% h - element length
% type - element type
% n - order of derivative
function [psin] = basis1d(s,h,type,n)
[coeff1d,~,order] = polyBasisTypes1d(type);
s = s(:);
np = size(s,1);
pows = order-n:-1:0; cows = factorial(pows+n)./factorial(pows);
basiss = h^(-n)*[cows.*s.^pows, zeros(np,n)];
psin = basiss*coeff1d';
end