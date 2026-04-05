% generates connectivity list for polyBasis element types
% map really only needed for periodic (and perfectly fitted) domains in 2D
function [p,t] = meshAddPolyBasis1d(pTR,tTR,type)

% local coordinates for nodal points
[~,B] = polyBasisTypes1d(type);

nel=size(tTR,1); nbf=size(B,1);
tnew=zeros(nel,nbf); pnew=zeros(nel*nbf,1);
cnt=1;
for k=1:nel
    globCoords=pTR(tTR(k,1))*(1-B) + pTR(tTR(k,2))*B;
    pnew(cnt:cnt+nbf-1,:)=globCoords;
    tnew(k,1:end)=(cnt:cnt+nbf-1);
    cnt=cnt+nbf;
end

% removing duplicate-nodes
[p,t]=meshRemoveDGv2(pnew,tnew);

end

function [p,t]=meshRemoveDGv2(p,tDG,tol)
if ~exist('tol', 'var')
    tol = 1e-12; % relative tolerance
end
[p,~,DG2CG]=uniquetol(p,tol,'ByRows',true);
t=zeros(size(tDG));
for i=1:size(t,1)
    for j=1:size(t,2)
        t(i,j)=DG2CG(tDG(i,j));
    end
end
end