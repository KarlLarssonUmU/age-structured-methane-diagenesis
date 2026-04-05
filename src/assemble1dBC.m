function [A,B] = assemblePoissonBC1d(TRp,TRt,p,t,type,k_val,c_air)

% assemble stiffness matrix
ndof = max(max(t)); A = sparse(ndof,ndof); B = zeros(ndof,1);

% left side of first element
k = 1; p1 = TRp(TRt(k,1),:); p2 = TRp(TRt(k,2),:); h = p2 - p1; s = 0;
psi = basis1d(s,h,type,0); psin = basis1d(s,h,type,1); dofs = t(k,:);
A(dofs,dofs) = A(dofs,dofs) + k_val*(psi'*psi);
B(dofs) = B(dofs) + k_val*(psi'*c_air);

% % right side of last element
% k = nel; p1 = TRp(TRt(k,1),:); p2 = TRp(TRt(k,2),:); h = p2 - p1; s = 1;
% psi = basis1d(s,h,type,0); psin = basis1d(s,h,type,1); dofs = t(k,:);
% A(dofs,dofs) = A(dofs,dofs) - psin'*psi - psi'*psin + beta/h*(psi'*psi);
% B(dofs) = B(dofs) - psin'*gD(2) + beta/h*(psi'*gD(2));

end