% Conservative weak form of ( -coeff*u' + v*u )' + kcoeff*u = f.
% Depth increases to the right; v is a signed downward transport velocity.
% For upward Darcy magnitude q, pass v=-q. Boundary flux is handled
% separately, leaving zero total flux as the natural bottom condition.
function [A,M,B,TRp,TRt,p,t] = assemble1dAdvectionDiffusion(a,b,N,type,coeff,kcoeff,f,v)
if ~exist('f','var')
    f = @(x) 0;
end
if ~exist('v','var')
    v = @(x) 0;
end

TRp = linspace(a,b,N+1)'; TRt = [1:N;2:N+1]';
[p,t] = meshAddPolyBasis1d(TRp,TRt,type);

% assemble stiffness matrix
nel = size(t,1); nbf = size(t,2);
ndof = max(t,[],'all'); A = sparse(ndof,ndof); M = sparse(ndof,ndof);
B = zeros(ndof,1);
[wvec,cmat]=quadratureGaussPoints(2*nbf,1);
for k=1:nel
    p1 = TRp(TRt(k,1),:); p2 = TRp(TRt(k,2),:); h = p2 - p1;
    Ak = zeros(nbf,nbf); Mk = zeros(nbf,nbf); Bk = zeros(nbf,1);
    for q=1:length(wvec)
        s = cmat(q,1)'; % local quadrature point
        x = (1-s)*p1 + s*p2; % global quadrature point
        dx = h*wvec(q); % global quadrature weight
        psi = basis1d(s,h,type,0); psin = basis1d(s,h,type,1);
        Ak = Ak + dx*(coeff(x)*(psin'*psin) - v(x)*(psin'*psi) + kcoeff(x)*(psi'*psi));
        Mk = Mk + dx*(psi'*psi);
        Bk = Bk + dx*psi'*f(x); % simple load vector for test
    end    
    dofs = t(k,:);
    A(dofs,dofs) = A(dofs,dofs) + Ak; M(dofs,dofs) = M(dofs,dofs) + Mk;
    B(dofs) = B(dofs) + Bk;
end

end
