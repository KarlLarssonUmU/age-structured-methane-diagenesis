function y = productionShape(z,peakDepth,rightWidth)
%PRODUCTIONSHAPE Original quintic bump in well-conditioned local coordinates.
% Returns [bump; background]; rightWidth is peak-to-end distance.
assert(isfinite(peakDepth) && peakDepth>=0 && isfinite(rightWidth) && rightWidth>0, ...
    'Peak depth must be nonnegative and width positive.');
z = z(:)';
if peakDepth==0
    % Exact surface endpoint: the smooth descending limb, without division by zero.
    s=min(max(z/rightWidth,0),1);
    bump=(1-s).^3.*(1+3*s+(8/3)*s.^2);
    bump(z<0)=0;
    y=[bump;ones(size(z))];
    return
end
bump = zeros(size(z));
left = z>=0 & z<=peakDepth;
s = z(left)/peakDepth;
bump(left) = (20*s.^3-25*s.^4+8*s.^5)/3;
right = z>peakDepth & z<=peakDepth+rightWidth;
s = (z(right)-peakDepth)/rightWidth;
bump(right) = 1-(10/3)*s.^2+5*s.^4-(8/3)*s.^5;
y = [max(bump,0);ones(size(z))];
end
