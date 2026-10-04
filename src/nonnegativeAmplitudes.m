function [amplitudes,errorC] = nonnegativeAmplitudes(G,y)
%NONNEGATIVEAMPLITUDES Exact enumeration of active sets for TWO columns.
assert(size(G,2)==2);
norms = vecnorm(G);
assert(all(norms>0) && all(isfinite(G),'all'));
H = G./norms;
candidates = zeros(2,3);
candidates(1,2) = max(0,H(:,1)'*y);
candidates(2,3) = max(0,H(:,2)'*y);
if rcond(H'*H)>1e-12
    interior = H\y;
    if all(interior>=0), candidates(:,end+1) = interior; end
end
errors = sum((H*candidates-y).^2,1);
[errorC,i] = min(errors);
amplitudes = candidates(:,i)./norms';
end
