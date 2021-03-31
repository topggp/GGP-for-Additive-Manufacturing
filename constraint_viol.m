% This function tests for violation of constraints of volume fraction and
% 3DP constraints. Redundant output has been given to satisfy
% requirement of FMINCON function.
% c = value of constraint
% ceq = Value of equality constraints = NULL
% dcdx = Value of supplied constraint gradient
% gradceq = Value of equality constraints gradient = NULL
function [c,ceq,dc,gradceq] = constraint_viol(X,nY,nelx,nely,np,Yc,xPhys,dE_dX,dE_dL,dE_dm,edofMat,U,KE,E,emptyelts,fullelts,drho_dX,drho_dL,drho_dm)
% [d,dd] = density1(X,ypos,nelx,nely,nX,nY);
% c=sum(d(:))-mMax;
% dcdx=sum(dd); %add ' for fmincon optimiser (without ALM constraints only !!)
ceq = [];
gradceq = [];
 ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx);
    c = sum(sum((E).*ce));
    v=mean(xPhys(:));
    dc_dE = -ce;
    dc_dE(emptyelts) = 0;
    dc_dE(fullelts) = 0;
    dc_dX=dE_dX*dc_dE(:);
    dc_dL=dE_dL*dc_dE(:);
    dc_dm=dE_dm*dc_dE(:);
    dc=zeros(size(X));
    dc(1:2:end)=dc_dX;
    dc(2:2:end)=dc_dL;
%     dc(3:3:end)=dc_dm;  %%
    dv_dxPhys = ones(nely,nelx)/nelx/nely;
    dv_dxPhys(emptyelts) = 0;
    dv_dxPhys(fullelts) = 0;
    dv_dX=drho_dX*dv_dxPhys(:);
    dv_dL=drho_dL*dv_dxPhys(:);
    dv_dm=drho_dm*dv_dxPhys(:);
    dv=zeros(size(X));
    dv(1:2:end)=dv_dX;
    dv(2:2:end)=dv_dL;
%% ALM constraints
% e = 0.05;       % Measure of confidence. Choose appropriately
% p = log(2*nX*(nY-1))/e;
p = 4;
almc = zeros(2*np*(nY-1),1);
a = zeros(8*np*(nY-1),3);
for i = 1:2:2*np*(nY-1)
    almc(i) = X(i)-X(i+1)/2-X(i+2*np)+X(i+2*np+1)/2-3;
    almc(i+1) = X(i+2*np)+X(i+2*np+1)/2-X(i)-X(i+1)/2-3;
    a(4*i-3,:) = [i i 1];
    a(4*i-2,:) = [i i+1 -0.5];
    a(4*i-1,:) = [i i+2*np -1];
    a(4*i,:) = [i i+2*np+1 0.5];
    a(4*i+1,:) = [i+1 i -1];
    a(4*i+2,:) = [i+1 i+1 -0.5];
    a(4*i+3,:) = [i+1 i+2*np 1];
    a(4*i+4,:) = [i+1 i+2*np+1 0.5];
end
dalmcdx = sparse(a(:,1),a(:,2),a(:,3),2*np*(nY-1),2*np*nY);
c = [c; log(sum(exp(p*almc))/(2*np*(nY-1)))/p]; % Add ' for fmincon optimiser
dc = [dc; exp(p*almc)'*dalmcdx/sum(exp(p*almc))]; % Add ' for fmincon optimiser
end