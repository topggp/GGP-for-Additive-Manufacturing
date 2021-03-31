function [W,dW_dX,dW_dL,dW_dh,dW_dtheta0,dW_dy0]=Wgp_modified2(x,y,Xc,p,np,nY,Yk,nely,nelx)
Xk=Xc(1:2:2*np*nY);        % Extraire tous les Xk
Lk=Xc(2:2:2*np*nY);        % Extraire tous les Lk
h=Xc((2*np*nY)+1:end-np-2);
theta0 = Xc(end-1);
y0 = Xc(end);
Xk=reshape(Xk,nY,np);
Lk=reshape(Lk,nY,np);
% Yk1=Yk;
Yk=reshape(Yk,nY,np);

% Gauss points sqrt(nelx^2+nely^2)/2+
N_GP=length(x);
xg = nelx/2+cos(theta0)*repmat(x-nelx/2,1,np*(nY-1))+sin(theta0)*repmat(y-nely/2,1,np*(nY-1));
yg = nely/2-y0-sin(theta0)*repmat(x-nelx/2,1,np*(nY-1))+cos(theta0)*repmat(y-nely/2,1,np*(nY-1));
dxg_theta0 = -sin(theta0)*repmat(x-nelx/2,1,np*(nY-1))'+cos(theta0)*repmat(y-nely/2,1,np*(nY-1))';
dyg_theta0 = -cos(theta0)*repmat(x-nelx/2,1,np*(nY-1))'-sin(theta0)*repmat(y-nely/2,1,np*(nY-1))';
dyg_y0 = -ones(size(yg))';

% xk & dxk_dXk
xk=Xk(1:end-1,:);
dxk_dXk=zeros(size(xk,1),size(Xk,1));
dxk_dXk(1:end,1:(end-1))=eye(size(xk,1));
dxk_dXk=kron(eye(np),sparse(dxk_dXk));
% xk+1 & dxk+1_dXk
xk1=Xk(2:end,:);
dxk1_dXk=zeros(size(xk1,1),size(Xk,1));
dxk1_dXk(1:end,2:(end))=eye(size(xk1,1));
dxk1_dXk=kron(eye(np),sparse(dxk1_dXk));% xk+1 & dxk+1_dXk
% yk & yk+1
yk=Yk(1:end-1,:);
yk1=Yk(2:end,:);
% lk & dlk_dLk
lk=Lk(1:end-1,:);
dlk_dLk=zeros(size(lk,1),size(Lk,1));
dlk_dLk(1:end,1:(end-1))=eye(size(lk,1));
dlk_dLk=kron(eye(np),sparse(dlk_dLk));
% lk+1 & dlk+1_dLk
lk1=Lk(2:end,:);
dlk1_dLk=zeros(size(lk,1),size(Lk,1));
dlk1_dLk(1:end,2:(end))=eye(size(lk1,1));
dlk1_dLk=kron(eye(np),sparse(dlk1_dLk));% ly & uy
lx=repmat(xk(:)-lk(:)/2,1,N_GP)'+(yg-repmat(yk(:),1,N_GP)')./repmat(yk1(:)-yk(:),1,N_GP)'.*repmat((xk1(:)-lk1(:)/2)-(xk(:)-lk(:)/2),1,N_GP)';
ux=repmat(xk(:)+lk(:)/2,1,N_GP)'+(yg-repmat(yk(:),1,N_GP)')./repmat(yk1(:)-yk(:),1,N_GP)'.*repmat((xk1(:)+lk1(:)/2)-(xk(:)+lk(:)/2),1,N_GP)';

% component angle lt and rt
slt = repmat((yk1(:)-yk(:))./(((xk1(:)-lk1(:)/2-xk(:)+lk(:)/2).^2+(yk1(:)-yk(:)).^2).^0.5),1,N_GP);
srt = repmat((yk1(:)-yk(:))./(((xk1(:)+lk1(:)/2-xk(:)-lk(:)/2).^2+(yk1(:)-yk(:)).^2).^0.5),1,N_GP);

cd1 = (2*((yk(:)-yk1(:)).^2+(lk(:)/2-lk1(:)/2-xk(:)+xk1(:)).^2).^(3/2));
dslt_dxk = -((yk(:)-yk1(:)).*(lk(:)-lk1(:)-2*xk(:)+2*xk1(:)))./cd1;
dslt_dxk = sparse(diag(dslt_dxk));
% dslt_dxk1 = ((yk(:)-yk1(:)).*(lk(:)-lk1(:)-2*xk(:)+2*xk1(:)))./(2*((yk(:)-yk1(:)).^2+(lk(:)/2-lk1(:)/2-xk(:)+xk1(:)).^2).^(3/2));
dslt_dxk1 = -dslt_dxk;
dslt_dlk = ((yk(:)-yk1(:)).*(lk(:)/2-lk1(:)/2-xk(:)+xk1(:)))./cd1;
dslt_dlk = sparse(diag(dslt_dlk));
dslt_dlk1 = -dslt_dlk;
cd2 = (2*((yk(:) - yk1(:)).^2 + (lk(:)/2 - lk1(:)/2 - xk(:) + xk1(:)).^2).^(3/2));
dsrt_dxk = -((yk(:) - yk1(:)).*(lk(:) - lk1(:) - 2*xk(:) + 2*xk1(:)))./cd2;
dsrt_dxk = sparse(diag(dsrt_dxk));
dsrt_dxk1 = - dsrt_dxk;
dsrt_dlk = ((yk(:) - yk1(:)).*(lk(:)/2 - lk1(:)/2 - xk(:) + xk1(:)))./cd2;
dsrt_dlk = sparse(diag(dsrt_dlk));
dsrt_dlk1 = - dsrt_dlk;

% derivé ly et uy 2704x90
dlx_xk= 1-(yg-repmat(yk(:),1,N_GP)')./repmat(yk1(:)-yk(:),1,N_GP)';
dlx_xk1= (yg-repmat(yk(:),1,N_GP)')./repmat(yk1(:)-yk(:),1,N_GP)';
dlx_Xk= dlx_xk*dxk_dXk+dlx_xk1*dxk1_dXk;
dlx_lk= (-1+(yg-repmat(yk(:),1,N_GP)')./repmat(yk1(:)-yk(:),1,N_GP)')/2;
dlx_lk1= -(yg-repmat(yk(:),1,N_GP)')./repmat(yk1(:)-yk(:),1,N_GP)'/2;
dlx_Lk= dlx_lk*dlk_dLk+dlx_lk1*dlk1_dLk;
dlx_theta0= repmat(((xk1(:)-lk1(:)/2)-(xk(:)-lk(:)/2))./(yk1(:)-yk(:)),1,N_GP).*dyg_theta0;
dlx_y0= repmat(((xk1(:)-lk1(:)/2)-(xk(:)-lk(:)/2))./(yk1(:)-yk(:)),1,N_GP).*dyg_y0;
dux_xk= 1-(yg-repmat(yk(:),1,N_GP)')./repmat(yk1(:)-yk(:),1,N_GP)';
dux_xk1= (yg-repmat(yk(:),1,N_GP)')./repmat(yk1(:)-yk(:),1,N_GP)';
dux_Xk= dux_xk*dxk_dXk+dux_xk1*dxk1_dXk;
dux_lk= (1-(yg-repmat(yk(:),1,N_GP)')./repmat(yk1(:)-yk(:),1,N_GP)')/2;
dux_lk1= (yg-repmat(yk(:),1,N_GP)')./repmat(yk1(:)-yk(:),1,N_GP)'/2;
dux_Lk= dux_lk*dlk_dLk+dux_lk1*dlk1_dLk;
dux_theta0= repmat(((xk1(:)-lk1(:)/2)-(xk(:)-lk(:)/2))./(yk1(:)-yk(:)),1,N_GP).*dyg_theta0;
dux_y0= repmat(((xk1(:)-lk1(:)/2)-(xk(:)-lk(:)/2))./(yk1(:)-yk(:)),1,N_GP).*dyg_y0;
% duy_Xk= repmat(dxk_dXk(:),1,np)+reshape(repmat((yg-yk(:))./(yk1(:)-yk(:)),nY,1),length(idy),[]).*repmat((dxk1_dXk(:)-(dxk_dXk(:))),1,np);
% duy_Lk= repmat(dlk_dLk(:)/2,1,np)+reshape(repmat((yg-yk(:))./(yk1(:)-yk(:)),nY,1),length(idy),[]).*repmat((dlk1_dLk(:)/2)-dlk_dLk(:)/2,1,np);
% duy_theta0= ((xk1(:)+lk1(:)/2)-(xk(:)+lk(:)/2))'./(yk1(:)-yk(:))'.*dyg_theta0;
% duy_y0= ((xk1(:)+lk1(:)/2)-(xk(:)+lk(:)/2))'./(yk1(:)-yk(:))'.*dyg_y0;
% b=nely*h*cos(theta0)-y0;
ly= repmat(yk(:),1,N_GP)';
uy= repmat(yk1(:),1,N_GP)';
b=sqrt(nelx^2+nely^2)*repmat(kron(h(:),ones(nY-1,1)),1,size(y,1))';
db_h =sqrt(nelx^2+nely^2)*ones(size(h,1)*(nY-1),size(y,1)); % 5x1
db_theta0=zeros(size(xg))';
db_y0=zeros(size(xg))';
switch p.method
    case 'MMC'
        alp = p.alp;
        alp2 = 2*alp;
        epsi=p.epsi;
        bet=p.bet;
        c = and((yg <= repmat(yk1(:)',N_GP,1)),(yg <= repmat(yk1(:)',N_GP,1)));
%         midx = (ux-lx)/2;
%         dmidx_dux = 0.5;
%         dmidx_dlx = -0.5;
%         
%         midy = repmat(Yk(:),1,N_GP)';
%         
%         rho = sqrt((xg-midx).^2+(yg-midy).^2);
%         drho_dmidx = (midx-xg)./(rho+(rho==0));
%         drho_dxg = (xg-midx)./(rho+(rho==0));
%         drho_dyg = (yg-midy)./(rho+(rho==0));
%         
%         phi=atan2(yg-midy,xg-midx)-pi/2;
%         dphi_yg = (xg-midx)./(rho.^2+(rho==0));
%         dphi_midx =((-midy+yg)./(rho.^2+(rho==0)));
%         dphi_xg = -dphi_midx;
%         
%         upsi=sqrt(rho.^2+1/4-rho.*abs(cos(phi))).*(((rho.*cos(phi)).^2)>=(1/4))+~(((rho.*cos(phi)).^2)>=(1/4)).*abs(rho.*sin(phi));
%         dupsi_drho=(2*rho-abs(cos(phi)))/2./(upsi+(upsi==0)).*((((rho.*cos(phi)).^2)>=(1/4)))+~(((rho.*cos(phi)).^2)>=(1/4)).*abs(sin(phi));
%         dupsi_dphi=(rho.*sign(cos(phi)).*sin(phi))/2./(upsi+(upsi==0)).*((((rho.*cos(phi)).^2)>=(1/4)))+~(((rho.*cos(phi)).^2)>=(1/4)).*rho.*sign(sin(phi)).*cos(phi);
%         
        % local chi
        chi0_fun = @(x,y) c.*(1-(4*(x-(ux+lx)/2).^2./(ux-lx).^2).^alp-((y).^2./(b.^2)).^(alp2));%CHANGE
        chi0 = chi0_fun(xg,yg)';
        chi0(chi0==0) = -1;
        dchi0_ux = (-(4*alp*(lx - xg).*((lx + ux - 2*xg).^2./(lx - ux).^2).^(alp - 1).*(lx + ux - 2*xg))./(lx - ux).^3)';%CHANGE
        dchi0_lx = ((4*alp*(ux - xg).*((lx + ux - 2*xg).^2./(lx - ux).^2).^(alp - 1).*(lx + ux - 2*xg))./(lx - ux).^3)';%CHANGE
        dchi0_b = (2*alp*(y.^2./b.^2).^(alp)./b)';%CHANGE
%         dchi0_xg = ; %CHANGE
%         dchi0_yg = ; %CHANGE
        [chi,dchi]=Aggregation_Pi(chi0,p);
        chi(chi<=-1e6)=-1e6; %1x2704
        dchi_ux = dchi0_ux.*dchi;
        dchi_lx = dchi0_lx.*dchi;
        dchi_b = dchi0_b.*dchi;
%         dchi_xg = dchi0_xg.*dchi; %VERIFY
%         dchi_yg = dchi0_yg.*dchi; %VERIFY
        W=(chi>epsi)+(chi<=epsi&chi>=-epsi).*(3/4*(1-bet)*(chi/epsi-chi.^3/3/epsi^3)+(1+bet)/2)+(chi<-epsi)*bet;
        dW_dchi=-3/4*(1/epsi-chi.^2/epsi^3).*(bet-1).*(abs(chi)<epsi);
        dW_chi3 = repmat(dW_dchi,np*(nY-1),1);
        dW_ux=dchi_ux.*dW_chi3;
        dW_lx=dchi_lx.*dW_chi3;
        dW_b=dchi_b.*dW_chi3;
        dW_xg=dchi_xg.*dW_chi3;
        dW_yg=dchi_yg.*dW_chi3;
        dW_dX=dxk_dXk'*(dW_ux.*dux_xk'+dW_lx.*dlx_xk')+dxk1_dXk'*(dW_ux.*dux_xk1'+dW_lx.*dlx_xk1');
        dW_dL=dlk_dLk'*(dW_ux.*dux_lk'+dW_lx.*dlx_lk')+dlk1_dLk'*(dW_ux.*dux_lk1'+dW_lx.*dlx_lk1');
        dW_dh=reshape(sum(reshape(db_h.*dW_b,nY-1,[])),np,[]);
        dW_dy0=sum(dW_ux.*dux_y0+dW_lx.*dlx_y0+dW_yg.*dyg_y0);
        dW_dtheta0=sum(dW_ux.*dux_theta0+dW_lx.*dlx_theta0+dW_xg.*dxg_theta0+dW_yg.*dyg_theta0);
    case 'GP'
        deltamin=p.deltamin;
        r=p.r;
        tv1=(-xg+lx)'.*slt;
        tv2=(xg-ux)'.*srt;
        tv3=(yg-b)';
        tv4=(-yg+ly)';
        tv5=(yg-uy)';
        %derivé tv1 et tv2
        dtv1_xk=+dlx_xk.*slt' + (-xg+lx)*dslt_dxk;
        dtv1_xk1=+dlx_xk1.*slt' + (-xg+lx)*dslt_dxk1;
        dtv2_xk=-dux_xk.*srt' + (xg-ux)*dsrt_dxk;
        dtv2_xk1=-dux_xk1.*srt' + (xg-ux)*dsrt_dxk1;
        dtv1_lk=+dlx_lk.*slt' + (-xg+lx)*dslt_dlk;
        dtv1_lk1=+dlx_lk1.*slt' + (-xg+lx)*dslt_dlk1;
        dtv2_lk=-dux_lk.*srt' + (xg-ux)*dsrt_dlk;
        dtv2_lk1=-dux_lk1.*srt' + (xg-ux)*dsrt_dlk1;
        dtv1_theta0=-dxg_theta0+dlx_theta0;
        dtv2_theta0=dxg_theta0-dux_theta0;
        dtv1_y0=+dlx_y0;
        dtv2_y0=-dux_y0;
        dtv3_h=-db_h;
        dtv3_theta0=dyg_theta0-db_theta0;
        dtv3_y0=dyg_y0-db_y0;
        dtv4_theta0=-dyg_theta0;
        dtv4_y0=-dyg_y0;        
        dtv5_theta0=dyg_theta0;
        dtv5_y0=dyg_y0;
        %
        
        [zetavar_micro,dzetavar_micro]=Aggregation_Pi([tv1(:) tv2(:) tv3(:) tv4(:) tv5(:)]',p);
        zetavar = reshape(zetavar_micro,np*(nY-1),size(y,1));
        dzetavar1 = reshape(dzetavar_micro(1,:),np*(nY-1),size(y,1));
        dzetavar2 = reshape(dzetavar_micro(2,:),np*(nY-1),size(y,1));
        dzetavar3 = reshape(dzetavar_micro(3,:),np*(nY-1),size(y,1));        
        dzetavar4 = reshape(dzetavar_micro(4,:),np*(nY-1),size(y,1));
        dzetavar5 = reshape(dzetavar_micro(5,:),np*(nY-1),size(y,1));
%         zetavar = zeros(np,size(y,1));
%         dzetavar = zeros(np*(nY-1),size(y,1));

%         zetavar = reshape(zetavar,np,size(y,1));
%         dzetavar1 = reshape(dzetavar,np*(nY-1),size(y,1)).*dzetavar1_micro;
%         dzetavar2 = reshape(dzetavar,np*(nY-1),size(y,1)).*dzetavar2_micro;
%         dzetavar3 = reshape(dzetavar,np*(nY-1),size(y,1)).*dzetavar3_micro;        
%         dzetavar4 = reshape(dzetavar,np*(nY-1),size(y,1)).*dzetavar4_micro;
%         dzetavar5 = reshape(dzetavar,np*(nY-1),size(y,1)).*dzetavar5_micro;
        dzetavar_xk = dzetavar1.*dtv1_xk' + dzetavar2.*dtv2_xk';
        dzetavar_xk1 = dzetavar1.*dtv1_xk1' + dzetavar2.*dtv2_xk1';
        dzetavar_lk = dzetavar1.*dtv1_lk' + dzetavar2.*dtv2_lk';
        dzetavar_lk1 = dzetavar1.*dtv1_lk1' + dzetavar2.*dtv2_lk1';
        dzetavar_h = dzetavar3.*dtv3_h;
        dzetavar_theta0 = dzetavar1.*dtv1_theta0 + dzetavar2.*dtv2_theta0 + dzetavar3.*dtv3_theta0 + dzetavar4.*dtv4_theta0 + dzetavar5.*dtv5_theta0;
        dzetavar_y0 = dzetavar1.*dtv1_y0 + dzetavar2.*dtv2_y0 + dzetavar3.*dtv3_y0 + dzetavar4.*dtv4_y0 + dzetavar5.*dtv5_y0;
        deltaiel_micro=(1/pi/r^2*(r^2*acos(zetavar/r)-zetavar.*sqrt(r^2-zetavar.^2))).*(abs(zetavar)<=r)+((zetavar<-r));
        ddeltaiel_micro_dzetavar=(-2*sqrt(r^2-zetavar.^2)/pi/r^2).*(abs(zetavar)<=r);
        p.aggregation='asymptotic';
        deltaiel=zeros(np,N_GP);
        ddeltaiel = zeros(np*(nY-1),N_GP);
        for i=1:np
            [deltaiel_i,ddeltaiel_i]=Aggregation_Pi(deltaiel_micro((nY-1)*(i-1)+1:(nY-1)*i,:),p);
            [deltaiel_i,ddelta_smooth] = smooth_sat(deltaiel_i,p);
            deltaiel(i,:) = deltaiel_i;
            ddeltaiel((nY-1)*(i-1)+1:(nY-1)*i,:) = ddeltaiel_i.*ddelta_smooth;
        end
        ddetlaiel_dzetavar = ddeltaiel_micro_dzetavar.*ddeltaiel;
        W=deltamin+(1-deltamin)*deltaiel;
        dW_ddeltaiel=(1-deltamin);
        dW_dX=dW_ddeltaiel*(dxk_dXk'*(ddetlaiel_dzetavar.*dzetavar_xk) + dxk1_dXk'*(ddetlaiel_dzetavar.*dzetavar_xk1));
        dW_dL=dW_ddeltaiel*(dlk_dLk'*(ddetlaiel_dzetavar.*dzetavar_lk) + dlk1_dLk'*(ddetlaiel_dzetavar.*dzetavar_lk1));
        dW_dh=dW_ddeltaiel*reshape(sum(reshape(ddetlaiel_dzetavar.*dzetavar_h,nY-1,[])),np,[]);
        dW_dtheta0=dW_ddeltaiel*reshape(sum(reshape(ddetlaiel_dzetavar.*dzetavar_theta0,nY-1,[])),np,[]);
        dW_dy0=dW_ddeltaiel*reshape(sum(reshape(ddetlaiel_dzetavar.*dzetavar_y0,nY-1,[])),np,[]);
    case 'MNA'
        gt=p.sigma;
        
        w =@(x) (0.5-(15/(16*gt))*x+(5/(8*gt^3))*x.^3-(3/(16*gt^5))*x.^5).*(x>=-gt&x<=gt)+(x<-gt);
        dw =@(x) (-15/(16*gt) +3*(5/(8*gt^3))*x.^2 -5*(3/(16*gt^5))*x.^4).*(x>=-gt&x<=gt);
        tv1=(-xg+lx)';
        tv2=(xg-ux)';
        tv3=(yg-b)';
        tv4=(-yg+ly)';
        tv5=(yg-uy)';
        Wel=w(tv1).*w(tv2).*w(tv3).*w(tv4).*w(tv5);
        p.aggregation='asymptotic';
%         W=zeros(np,N_GP);
%         dW = zeros(np*(nY-1),N_GP);
%         for i=1:np
%             [W_i,dW_i]=Aggregation_Pi(Wel((nY-1)*(i-1)+1:(nY-1)*i,:),p);
%             [W_i,dW_smooth] = smooth_sat(W_i,p);
%             W(i,:) = W_i;
%             dW((nY-1)*(i-1)+1:(nY-1)*i,:) = dW_i.*dW_smooth;
%         end
        
        [W,dW]=Aggregation_Pi(reshape(Wel,nY-1,[]),p);
        [W,dW_smooth] = smooth_sat(W,p);
        dW = dW.*repmat(dW_smooth,nY-1,1);
        W = reshape(W,np,[]);
        dW = reshape(dW,np*(nY-1),[]);
% %         W=Wel;
% %         dW=ones(size(W));
        %derivé tv1 et tv2
        %derivé tv1 et tv2
        
        dtv1_xk=+dlx_xk';
        dtv1_xk1=+dlx_xk1';
        dtv2_xk=-dux_xk';
        dtv2_xk1=-dux_xk1';
        dtv1_lk=+dlx_lk';
        dtv1_lk1=+dlx_lk1';
        dtv2_lk=-dux_lk';
        dtv2_lk1=-dux_lk1';
        dtv1_theta0=-dxg_theta0+dlx_theta0;
        dtv2_theta0=dxg_theta0-dux_theta0;
        dtv1_y0=+dlx_y0;
        dtv2_y0=-dux_y0;
        dtv3_h=-db_h;
        dtv3_theta0=dyg_theta0-db_theta0;
        dtv3_y0=dyg_y0-db_y0;
        dtv4_theta0=-dyg_theta0;
        dtv4_y0=-dyg_y0;        
        dtv5_theta0=dyg_theta0;
        dtv5_y0=dyg_y0;
        %
        dW345 = dW.*w(tv3).*w(tv4).*w(tv5);
        dW_xk =dW345.*((dw(tv1).*w(tv2).*dtv1_xk+ dw(tv2).*w(tv1).*dtv2_xk));
        dW_lk =dW345.*((dw(tv1).*w(tv2).*dtv1_lk+ dw(tv2).*w(tv1).*dtv2_lk));
        dW_xk1 =dW345.*((dw(tv1).*w(tv2).*dtv1_xk1+dw(tv2).*w(tv1).*dtv2_xk1)) ;
        dW_lk1 =dW345.*((dw(tv1).*w(tv2).*dtv1_lk1+ dw(tv2).*w(tv1).*dtv2_lk1));
        dW_dX=dxk_dXk'*dW_xk+dxk1_dXk'*dW_xk1;
        dW_dL=dlk_dLk'*dW_lk+dlk1_dLk'*dW_lk1;
        W2345 = w(tv2).*w(tv3).*w(tv4).*w(tv5);
        W1345 = w(tv1).*w(tv3).*w(tv4).*w(tv5);
        W1245 = w(tv1).*w(tv2).*w(tv4).*w(tv5);
        W1235 = w(tv1).*w(tv2).*w(tv3).*w(tv5);
        W1234 = w(tv1).*w(tv2).*w(tv3).*w(tv4);
        dW_dh = reshape(sum(reshape(dW.*dw(tv3).*dtv3_h.*W1245,nY-1,[])),np,[]);
        dW_dtheta0 =reshape(sum(reshape(dW.*(dw(tv1).*W2345.*dtv1_theta0+dw(tv2).*W1345.*dtv2_theta0+W1245.*dw(tv3).*dtv3_theta0+W1235.*dw(tv4).*dtv4_theta0+W1234.*dw(tv5).*dtv5_theta0),nY-1,[])),np,[]);
        dW_dy0 =reshape(sum(reshape(dW.*(dw(tv1).*W2345.*dtv1_y0+dw(tv2).*W1345.*dtv2_y0+W1245.*dw(tv3).*dtv3_y0+W1235.*dw(tv4).*dtv4_y0+W1234.*dw(tv5).*dtv5_y0),nY-1,[])),np,[]);
end