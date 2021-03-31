function [E,dE_ddelta,dE_dm]=model_updateM(delta,p,X,np,nY)

m=X((2*np*nY)+np+1:end-2);
nc=size(delta,1);
m=repmat(m(:),1,size(delta,2));
switch p.method
    case 'MMC'
        %update the Young Modulus on the base of delta
        E=p.E0*delta(1,1:end);
        dE_ddelta=p.E0*ones(size(delta));
%         dE_ddelta=repmat(dE_ddelta,size(m,1),1);
        dE_dm=0*m;
    case 'GP'
        hatdelta=delta.*m.^p.gammac;
        [E,dE_dhatdelta]=Aggregation_Pi(hatdelta,p);
        if p.saturation
        [E,ds]=smooth_sat(E,p,nc);
        dE_dhatdelta=ds.*dE_dhatdelta;
        end
        E=E.*p.E0;
        dhatdelta_ddelta=m.^p.gammac;
        dhatdelta_dm=p.gammac*delta.*m.^(p.gammac-1);
        dE_ddelta=p.E0*dhatdelta_ddelta.*dE_dhatdelta;
        dE_dm=p.E0*dE_dhatdelta.*dhatdelta_dm;
    case 'MNA'
        hatdelta=delta.*m.^p.gammac;
        [rho,drho_dhatdelta]=Aggregation_Pi(hatdelta,p);
        if p.saturation
        [rho,ds]=smooth_sat(rho,p,nc);
        drho_dhatdelta=repmat(ds,size(drho_dhatdelta,1),1).*drho_dhatdelta;
        end
        E=rho.^p.penalty.*(p.E0-p.Emin)+p.Emin;
        dhatdelta_ddelta=m.^p.gammac;
        dhatdelta_dm=p.gammac*delta.*m.^(p.gammac-1);
        dE_ddelta=p.penalty*(p.E0-p.Emin)*dhatdelta_ddelta.*drho_dhatdelta.*repmat(rho,size(drho_dhatdelta,1),1).^(p.penalty-1);
        dE_dm=p.penalty*(p.E0-p.Emin)*dhatdelta_dm.*drho_dhatdelta.*repmat(rho,size(drho_dhatdelta,1),1).^(p.penalty-1);
end