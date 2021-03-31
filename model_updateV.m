function [rho,drho_ddelta,drho_dm]=model_updateV(delta,p,X,np,nY)

m=X((2*np*nY)+np+1:end-2);
nc=size(delta,1);
m=repmat(m(:),1,size(delta,2));
switch p.method
    case 'MMC'
%         %update the Young Modulus on the base of delta
        rho=delta(1,1:end);
        drho_ddelta=ones(size(delta));
%         drho_ddelta=repmat(drho_ddelta,size(m,1),1);
        drho_dm=0*m;
    case 'GP'
        hatdelta=delta.*m.^p.gammav;
        [rho,drho_dhatdelta]=Aggregation_Pi(hatdelta,p);
        if p.saturation
        [rho,ds]=smooth_sat(rho,p,nc);
        drho_dhatdelta=ds.*drho_dhatdelta;
        end
        dhatdelta_ddelta=m.^p.gammav;
        dhatdelta_dm=p.gammav*delta.*m.^(p.gammav-1);
        drho_ddelta=dhatdelta_ddelta.*drho_dhatdelta;
        drho_dm=drho_dhatdelta.*dhatdelta_dm;
    case 'MNA'
        hatdelta=delta.*m.^p.gammav;
        [rho,drho_dhatdelta]=Aggregation_Pi(hatdelta,p);
        if p.saturation
        [rho,ds]=smooth_sat(rho,p,nc);
        drho_dhatdelta=repmat(ds,size(drho_dhatdelta,1),1).*drho_dhatdelta;
        end
        dhatdelta_ddelta=m.^p.gammav;
        dhatdelta_dm=p.gammav*delta.*m.^(p.gammav-1);
        drho_ddelta=dhatdelta_ddelta.*drho_dhatdelta;
        drho_dm=drho_dhatdelta.*dhatdelta_dm;
end