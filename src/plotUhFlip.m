function plotUhFlip(TRp,TRt,t,Uh,type)

X = zeros(0,1); U = zeros(0,1);
for k=1:size(TRt,1)
    p1 = TRp(TRt(k,1),:); p2 = TRp(TRt(k,2),:); h = p2 - p1;
    s = linspace(0,1,20)'; % local evaluation points
    x = (1-s)*p1 + s*p2; % global evaluation points
    psi = basis1d(s,h,type,0);
    dofs = t(k,:);
    u = psi*Uh(dofs);
    X=[X;x]; U=[U;u];
end

plot(U,X,'LineWidth',2)

end