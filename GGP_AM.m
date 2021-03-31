%% Main m file for Generalized geometric projection Dec,2018
% moa = Maximum overhang angle permitted in degrees
% mm = number of constraints
% BL = Bridge length specifier
%step 1) Problem set-up
function GGP_AM(nelx,nely,volfrac,BC,method,np,mm,BL,theta0,y0)
close all
tic
p.method = method;
settings = method;
% volfrac=0.4;
% settings='MNA';
stopping_criteria='kktnorm';
switch settings
    case 'GGP'
%         nelx = 52;
%         nely = 52;
%         np = 5;    % Number of deformable elements along x
        nY = 18;   % Number of deformable elements along y
%         BC='L-shape';%L-shape %Short_Cantiliever%MBB
%         p.method='MNA';%MMC%MNA 0%GP
        q=1;%q=1
        p.zp=1 ;% parameter for p-norm/mean regularization
        p.alp=1; %parameter for MMC
        p.epsi=0.7;% parameter for MMC
        p.bet=1e-3; %parameter for MMC
        p.deltamin=1e-6; %parameter for GP
        p.r=.5;%parameter for GP
        minh=3;
        p.sigma=1;%parameter for MNA
        p.gammav=1;%parameter for GP
        p.gammac=3;%parameter for GP
        p.penalty=3;%parameter for MNA
        p.aggregation='KSl'; %parameter for the aggregation function to be used
        % IE= Induced Exponential % KS= KS function %KSl= lowerbound KS function
        % p-norm %p-mean
        p.ka=10; % parameter for the aggregation constant
        p.saturation=true; % switch for saturation
%         ncx=1; % number of components in the x direction
%         ncy=1; % number of components in the y direction
        Ngp=1; % number of Gauss point per sampling window
        R=sqrt(1)/2; % radius of the sampling window (infty norm)
        initial_d=0.5;
    case 'MMC'
%         nelx = 100;
%         nely = 100;
%         np = 1;                 % Number of deformable elements along x
        nY = 18;   % Number of deformable elements along y
%         BC='Short_Cantilever';%L-shape %Short_Cantiliever
%         p.method='MMC';%MMC%MNA %GP
        q=3;%q=1
        p.zp=1 ;% parameter for p-norm/mean regularization
        p.alp=2; %parameter for MMC
        p.epsi=0.7;% parameter for MMC
        p.bet=1e-3; %parameter for MMC
        p.deltamin=1e-6; %parameter for GP
        p.r=1.5;%parameter for GP
        minh=2;
        p.sigma=1.5;%parameter for MNA
        p.gammav=1;%parameter for GP
        p.gammac=3;%parameter for GP
        p.penalty=3;%parameter for MNA
        p.aggregation='KSl'; %parameter for the aggregation function to be used
        % IE= Induced Exponential % KS= KS function %KSl= lowerbound KS function
        % p-norm %p-mean
        p.ka=4; % parameter for the aggregation constant
        p.saturation=false; % switch for saturation
%         ncx=1; % number of components in the x direction
%         ncy=1; % number of components in the y direction
        Ngp=1; % number of Gauss point per sampling window
        R=sqrt(3)/2; % radius of the sampling window (infty norm)
        initial_d=1; % component initial mass
    case 'MNA'
%         nelx = 100;
%         nely = 100;
%         np = 5;                 % Number of deformable elements along x
        nY = 18;   % Number of deformable elements along y
%         BC='MBB';%L-shape %Short_Cantiliever
%         p.method='MNA';%MMC%MNA %GP
        q=1;%q=1
        p.zp=1 ;% parameter for p-norm/mean regularization
        p.alp=1; %parameter for MMC
        p.epsi=0.7;% parameter for MMC
        p.bet=1e-3; %parameter for MMC
        p.deltamin=1e-6; %parameter for GP
        p.r=3;%parameter for GP
        minh=3;
        p.sigma=3;%parameter for MNA
        p.gammav=1;%parameter for GP
        p.gammac=1;%parameter for GP 3
        p.penalty=3;%parameter for MNA 3
        p.aggregation='KSl'; %parameter for the aggregation function to be used
        % IE= Induced Exponential % KS= KS function %KSl= lowerbound KS function
        % p-norm %p-mean
        p.ka=10; % parameter for the aggregation constant
        p.saturation=true; % switch for saturation
%         ncx=1; % number of components in the x direction
%         ncy=1; % number of components in the y direction
        Ngp=1; % number of Gauss point per sampling window
        R=sqrt(1)/2; % radius of the sampling window (infty norm)
        initial_d=0.5;
    case 'GP'
%         nelx = 100;
%         nely = 100;
%         np = 5;                 % Number of deformable elements along x
        nY = 18;   % Number of deformable elements along y
%         BC='Short_Cantilever';%L-shape %Short_Cantiliever
%         p.method='GP';%MMC%MNA %GP
        q=1;%q=1
        p.zp=1 ;% parameter for p-norm/mean regularization
        p.alp=1; %parameter for MMC
        p.epsi=0.7;% parameter for MMC
        p.bet=1e-3; %parameter for MMC
        p.deltamin=1e-6; %parameter for GP
        p.r=1.5;%parameter for GP
        minh=3;
        p.sigma=1.5;%parameter for MNA
        p.gammav=1;%parameter for GP
        p.gammac=3;%parameter for GP
        p.penalty=3;%parameter for MNA
        p.aggregation='KS'; %parameter for the aggregation function to be used
        % IE= Induced Exponential % KS= KS function %KSl= lowerbound KS function
        % p-norm %p-mean
        p.ka=40; % parameter for the aggregation constant
        p.saturation=false; % switch for saturation
%         ncx=1; % number of components in the x direction
%         ncy=1; % number of components in the y direction
        Ngp=1; % number of Gauss point per sampling window
        R=sqrt(1)/2; % radius of the sampling window (infty norm)
        initial_d=0.5;
end
cross_starting_guess=true;
%% Generate a *folder* and a prefix to save images *optimization history*:

rs=replace(num2str(R,'%3.2f'),'.','_');
vf = replace(num2str(volfrac,'%3.2f'),'.','_');
folder_name=[settings,'_OptimHis_',BC,'_',p.method,'_Volfrac_',vf,'_nelx_',num2str(nelx),...
    '_nely_',num2str(nely),'_ncx&y_',num2str(np),'_',num2str(nY),'_BL_',num2str(BL),'_theta0_',...
    num2str(theta0/pi*180),'_NoC_',num2str(mm),'_y0_',num2str(y0),'_Ngp_',num2str(Ngp),'_SC_',stopping_criteria];
image_prefix=[settings,'_',BC,'_',p.method,'_Volfrac_',vf,'_nelx_',num2str(nelx),'_nely_',num2str(nely),'_R_',rs,'_Ngp_',num2str(Ngp)];
mkdir(folder_name)
Path=[folder_name,'/'];
%% MATERIAL PROPERTIES
p.E0 = 1;
p.Emin = 1e-6;
nu = 0.3;
%% PREPARE FINITE ELEMENT ANALYSIS
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
KE = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);
nele = nelx*nely;
nnodes = (nelx+1)*(nely+1);
ndof = 2*nnodes;
nodenrs = reshape(1:nnodes,1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);
U = zeros(ndof,1);
%define the nodal coordinates
[Yy,Xx]=find(nodenrs);
Yy=nely+1-Yy;
Xx=Xx-1;
% Element connectivity
enodeMat=edofMat(:,[2,4,6,8])/2;
% compute the centroid coordinates
xc=mean(Xx(enodeMat'));
yc=mean(Yy(enodeMat'));
centroid_coordinate=[xc(:),yc(:)];
a=-R;
b=R;
[gpc,wc]=lgwt(Ngp,a,b);
[gpcx,gpcy]=meshgrid(gpc,gpc);
gauss_weight=wc*wc';
gpcx=reshape((repmat(gpcx(:),1,size(centroid_coordinate,1)))',[],1);
gpcy=reshape((repmat(gpcy(:),1,size(centroid_coordinate,1)))',[],1);
gauss_weight=reshape((repmat(gauss_weight(:),1,size(centroid_coordinate,1)))',[],1);
cc=repmat(centroid_coordinate,Ngp^2,1);
gauss_point=cc+[gpcx,gpcy];
[ugp,~,idgp]=unique(gauss_point,'rows');
%% DEFINE LOADS AND SUPPORTS
switch BC
    case 'MBB'
        excitation_node=1;excitation_direction=2;
        amplitude=-1;
        F = sparse(2*(excitation_node-1)+excitation_direction,1,amplitude,2*(nely+1)*(nelx+1),1);
        fixednodes=[find(Xx==min(Xx));(nelx+1)*(nely+1)];fixed_dir=[ones(nely+1,1);2];
        fixeddofs=2*(fixednodes-1)+fixed_dir;
        emptyelts=[]; fullelts = [];
    case 'Short_Cantilever'
        excitation_node=find((Xx==max(Xx))&(Yy==fix(0.5*min(Yy)+0.5*max(Yy))));excitation_direction=2;
        amplitude=-1;
        F = sparse(2*(excitation_node-1)+excitation_direction,1,amplitude,2*(nely+1)*(nelx+1),1);
        fixednodes=repmat(find(Xx==min(Xx)),2,1);fixed_dir=[ones(nely+1,1);2*ones(nely+1,1)];
        fixeddofs=2*(fixednodes-1)+fixed_dir(:);
        emptyelts=[]; fullelts = [];
    case 'L_Shape'
        excitation_node=find((Xx==max(Xx))&(Yy==fix(0.5*min(Yy)+0.5*max(Yy))));excitation_direction=2;
        amplitude=-1;
        F = sparse(2*(excitation_node-1)+excitation_direction,1,amplitude,2*(nely+1)*(nelx+1),1);
        fixednodes=repmat(find(Yy==max(Yy)),2,1);fixed_dir=[ones(nelx+1,1),2*ones(nelx+1,1)];
        fixeddofs=2*(fixednodes-1)+fixed_dir(:);
        emptyelts=find(xc>=(((max(Xx)+min(Xx))/2))&(yc>=((max(Yy)+min(Yy))/2)));
        fullelts = [];
    otherwise
        error('BC string should be a valid entry: ''MBB'',''L-Shape'',''Short_Cantiliever''')
end
mult = 1-numel(emptyelts)/nele;
volfrac = volfrac*mult; 
alldofs = [1:2*(nely+1)*(nelx+1)];
freedofs = setdiff(alldofs,fixeddofs);
%% INITIALIZE ITERATION
% define the initial guess components
% xp=nelx/(np+1):nelx/(np+1):np/(np+1)*nelx;
xp = linspace(nelx/4,3*nelx/4,np);
% yp=linspace(0,sqrt(nely^2+nelx^2)+1,nY);
yp = linspace(0,sqrt(nelx^2+nely^2), nY);
%
[xx,yy]=meshgrid(xp,yp);
if cross_starting_guess
    Xc=xx(:);
%     Xc(2:2:end)=Xc(2:2:end)+10;
    Yc=yy(:);
    Lc=5*ones(size(Xc));
%     Lc=(nelx/4/np)*ones(size(Xc));
end
% Mc=initial_d*ones(size([Xc;Yc]));
m=(linspace(1,1,np))';
h=((linspace(1,1,np))/sqrt(1))';
% h= [ 0 1 1 1 1 ]';
Xg=reshape([Xc,Lc]',[],1);
% theta0 = -pi/2;
% y0 = -25;%-25;
Xg=[Xg;h;m;theta0;y0];
[W,~,~,~,~,~] = Wgp_modified2(ugp(:,1),ugp(:,2),Xg,p,np,nY,Yc,nely,nelx);
delta=sum(reshape(W(:,idgp).*repmat(gauss_weight(:)',size(W,1),1),size(W,1),...
    [],Ngp^2),3)./sum(reshape(repmat(gauss_weight(:)',size(W,1),1),size(W,1),[],Ngp^2),3);
[rho,~,~]=model_updateV(delta,p,Xg,np,nY);
vfadj = mean(rho(:))/volfrac;
%% To obtain full component height with adjusted length
% Xg(2:2:2*length(Lc)) = Xg(2:2:2*length(Lc))/vfadj;
%% To obtain partial component height with initial length
if vfadj > 1
    Xg(2*length(Lc)+1:1:2*length(Lc)+np) = Xg(2*length(Lc)+1:1:2*length(Lc)+np)/vfadj;
end
%% Adimensional variable so that all variables are between 0 and 1
Xl=-sqrt((max(Xx)-min(Xx)).^2+(max(Yy)-min(Yy)).^2)*ones(size(Xc));
Xu=sqrt((max(Xx)-min(Xx)).^2+(max(Yy)-min(Yy)).^2)*ones(size(Xc))+1;
Ll=minh*ones(size(Xc));
Lu=(nelx/np)*ones(size(Xc));
hl=0.2*ones(size(h));
hu=ones(size(h));
ml=zeros(size(m));
mu=ones(size(m));
thetal=-pi/2;
thetau=pi/2;
yl=-nely/2;
yu=nely/2;
lower_bound=reshape([Xl,Ll]',[],1);lower_bound= [lower_bound ; hl; ml; thetal; yl];
upper_bound=reshape([Xu,Lu]',[],1);upper_bound= [upper_bound ; hu ; mu; thetau; yu];
X=(Xg-lower_bound)./(upper_bound-lower_bound);
%% Iteration initialisation

loop = 0;
change = 1;
% mm = 3;
n = length(X(:));
epsimin = 0.0000001;
eeen    = ones(n,1);
eeem    = ones(mm,1);
zeron   = zeros(n,1);
zerom   = zeros(mm,1);
xval    = X(:);
xold1   = xval;
xold2   = xval;
xmin    = zeron;
xmax    = eeen;
low     = xmin;
upp     = xmax;
C       = 1000*eeem;
d       = 0*eeem;
a0      = 1;
a       = zerom;
outeriter = 0;
maxoutit  = 2000;
kkttol  =0.001;
changetol=0.001;
kktnorm = kkttol+10;
outit = 0;
almc = 0; blc = 0;
mtv = 0; mbl = 0;
%% START ITERATION
cvec=zeros(maxoutit,1);
vvec=cvec;ovvec=cvec;gvec=cvec;pvec=cvec;
plot_rate=5;
% transition=50;
% change_of_formultion=200;
% change_of_formultion2=20000;
% active_elements=setdiff(1:nelx*nely,[emptyelts(:);fullelts(:)]);

switch stopping_criteria
    case 'kktnorm'
        stop_cond=outit < maxoutit && kktnorm>kkttol;
    case 'change'
        stop_cond=outit < maxoutit &&change>changetol;
end
%% Write logs in file for tracking command window outputs
fname = fullfile(folder_name, image_prefix);
if exist([fname,'.txt']) == 2
    f1 = fopen([fname,'.txt'],'at+');
else
    f1 = fopen([fname,'.txt'],'wt+');
end
if exist([image_prefix,'.mat']) == 2
        load([image_prefix,'.mat']);
end
storeXg = zeros(length(Xg),2000);
% storexPhys = zeros(nelx*nely,2000);
%% Start loop
while  stop_cond
    %     change>0.001&&
    outit   = outit+1;
    outeriter = outeriter+1;
    %     if outit == transition
    %          p.sigma= 1;p.saturation=true;
    %     end
    %% Project component on DZ
    % rotation
    % ugp_rot = ugp*[cos(alpha) sin(alpha);sin(alpha) cos(alpha)];
    % translation
    % ugp_rot(:,2) = ugp_rot(:,2)+DELTAY;
    if outeriter==80
        p.sigma=3;%parameter for MNA
        p.gammac=2;%parameter for GP 3
        p.penalty=2;%parameter for MNA 3
    elseif outeriter==160
        p.gammac=2;%parameter for GP 3
        p.penalty=2;%parameter for MNA 3
    end
    [W,dW_dX,dW_dL,dW_dh,dW_dtheta0,dW_dy0]=Wgp_modified2(ugp(:,1),ugp(:,2),Xg,p,np,nY,Yc,nely,nelx);
    %generalized projection
    delta=sum(reshape(W(:,idgp).*repmat(gauss_weight(:)',size(W,1),1),size(W,1),[],Ngp^2),3)./sum(reshape(repmat(gauss_weight(:)',size(W,1),1),size(W,1),[],Ngp^2),3);
    ddelta_dX=sum(reshape(full(dW_dX(:,idgp).*repmat(gauss_weight(:)',size(dW_dX,1),1)),size(dW_dX,1),[],Ngp^2),3)./sum(reshape(repmat(gauss_weight(:)',size(dW_dX,1),1),size(dW_dX,1),[],Ngp^2),3);
    ddelta_dL=sum(reshape(full(dW_dL(:,idgp).*repmat(gauss_weight(:)',size(dW_dX,1),1)),size(dW_dX,1),[],Ngp^2),3)./sum(reshape(repmat(gauss_weight(:)',size(dW_dX,1),1),size(dW_dX,1),[],Ngp^2),3);
    ddelta_dh=sum(reshape(dW_dh(:,idgp).*repmat(gauss_weight(:)',size(h,1),1),size(h,1),[],Ngp^2),3)./sum(reshape(repmat(gauss_weight(:)',size(h,1),1),size(h,1),[],Ngp^2),3);
    delta_c=sum(reshape(W(:,idgp).^q.*repmat(gauss_weight(:)',size(W,1),1),size(W,1),[],Ngp^2),3)./sum(reshape(repmat(gauss_weight(:)',size(W,1),1),size(W,1),[],Ngp^2),3);
    switch settings
        case 'MMC'
            ddelta_dtheta0=sum(reshape(full(dW_dtheta0(:,idgp).*repmat(gauss_weight(:)',size(dW_dtheta0,1),1)),size(dW_dtheta0,1),[],Ngp^2),3)./sum(reshape(repmat(gauss_weight(:)',size(dW_dtheta0,1),1),size(dW_dtheta0,1),[],Ngp^2),3);
            ddelta_dy0=sum(reshape(full(dW_dy0(:,idgp).*repmat(gauss_weight(:)',size(dW_dy0,1),1)),size(dW_dy0,1),[],Ngp^2),3)./sum(reshape(repmat(gauss_weight(:)',size(dW_dy0,1),1),size(dW_dy0,1),[],Ngp^2),3);
            ddelta_c_dX=sum(reshape(full(q*dW_dX(:,idgp).*repmat(W(:,idgp)'.^(q-1),1,nY*size(h,1))'.*repmat(gauss_weight(:)',size(dW_dX,1),1)),size(dW_dX,1),[],Ngp^2),3)./sum(reshape(repmat(gauss_weight(:)',size(dW_dX,1),1),size(dW_dX,1),[],Ngp^2),3);
            ddelta_c_dL=sum(reshape(full(q*dW_dL(:,idgp).*repmat(W(:,idgp)'.^(q-1),1,nY*size(h,1))'.*repmat(gauss_weight(:)',size(dW_dX,1),1)),size(dW_dX,1),[],Ngp^2),3)./sum(reshape(repmat(gauss_weight(:)',size(dW_dX,1),1),size(dW_dX,1),[],Ngp^2),3);
            ddelta_c_dh=sum(reshape(q*dW_dh(:,idgp).*repmat(W(:,idgp),size(h,1),1).^(q-1).*repmat(gauss_weight(:)',size(h,1),1),size(h,1),[],Ngp^2),3)./sum(reshape(repmat(gauss_weight(:)',size(h,1),1),size(h,1),[],Ngp^2),3);
            ddelta_c_dtheta0=sum(reshape(q*dW_dtheta0(:,idgp).*repmat(W(:,idgp),size(dW_dtheta0,1),1).^(q-1).*repmat(gauss_weight(:)',size(dW_dtheta0,1),1),size(dW_dtheta0,1),[],Ngp^2),3)./sum(reshape(repmat(gauss_weight(:)',size(dW_dtheta0,1),1),size(dW_dtheta0,1),[],Ngp^2),3);
            ddelta_c_dy0=sum(reshape(q*dW_dy0(:,idgp).*repmat(W(:,idgp),size(dW_dy0,1),1).^(q-1).*repmat(gauss_weight(:)',size(dW_dy0,1),1),size(dW_dy0,1),[],Ngp^2),3)./sum(reshape(repmat(gauss_weight(:)',size(dW_dy0,1),1),size(dW_dy0,1),[],Ngp^2),3);
        otherwise
            ddelta_dtheta0=sum(reshape(dW_dtheta0(:,idgp).*repmat(gauss_weight(:)',size(h,1),1),size(h,1),[],Ngp^2),3)./sum(reshape(repmat(gauss_weight(:)',size(h,1),1),size(h,1),[],Ngp^2),3);
            ddelta_dy0=sum(reshape(dW_dy0(:,idgp).*repmat(gauss_weight(:)',size(h,1),1),size(h,1),[],Ngp^2),3)./sum(reshape(repmat(gauss_weight(:)',size(h,1),1),size(h,1),[],Ngp^2),3);
            ddelta_c_dX=sum(reshape(full(q*dW_dX(:,idgp).*reshape(repmat(W(:,idgp)'.^(q-1),nY,1),size(W,2),[])'.*repmat(gauss_weight(:)',size(dW_dX,1),1)),size(dW_dX,1),[],Ngp^2),3)./sum(reshape(repmat(gauss_weight(:)',size(dW_dX,1),1),size(dW_dX,1),[],Ngp^2),3);
            ddelta_c_dL=sum(reshape(full(q*dW_dL(:,idgp).*reshape(repmat(W(:,idgp)'.^(q-1),nY,1),size(W,2),[])'.*repmat(gauss_weight(:)',size(dW_dX,1),1)),size(dW_dX,1),[],Ngp^2),3)./sum(reshape(repmat(gauss_weight(:)',size(dW_dX,1),1),size(dW_dX,1),[],Ngp^2),3);
            ddelta_c_dh=sum(reshape(q*dW_dh(:,idgp).*W(:,idgp).^(q-1).*repmat(gauss_weight(:)',size(W,1),1),size(W,1),[],Ngp^2),3)./sum(reshape(repmat(gauss_weight(:)',size(W,1),1),size(W,1),[],Ngp^2),3);
            ddelta_c_dtheta0=sum(reshape(q*dW_dtheta0(:,idgp).*W(:,idgp).^(q-1).*repmat(gauss_weight(:)',size(W,1),1),size(W,1),[],Ngp^2),3)./sum(reshape(repmat(gauss_weight(:)',size(W,1),1),size(W,1),[],Ngp^2),3);
            ddelta_c_dy0=sum(reshape(q*dW_dy0(:,idgp).*W(:,idgp).^(q-1).*repmat(gauss_weight(:)',size(W,1),1),size(W,1),[],Ngp^2),3)./sum(reshape(repmat(gauss_weight(:)',size(W,1),1),size(W,1),[],Ngp^2),3);
    end
            %model update
    [E,dE,dE_dm]=model_updateM(delta_c,p,X,np,nY);
    [rho,drho_ddelta,drho_dm]=model_updateV(delta,p,X,np,nY);
%     dE=reshape(repmat(dE',nY,1),size(ddelta_c_dX,2),[])';
%     drho_ddelta=reshape(repmat(drho_ddelta',nY,1),size(ddelta_c_dX,2),[])';
%     drho_dm=reshape(repmat(drho_dm',nY,1),size(ddelta_c_dX,2),[])';
%     dE_dm=reshape(repmat(dE_dm',nY,1),size(ddelta_c_dX,2),[])';
    
    dE_dh=dE.*ddelta_c_dh;
    dE_dtheta0=sum(dE.*ddelta_c_dtheta0,1);
    dE_dy0=sum(dE.*ddelta_c_dy0,1);
       
    drho_dh=drho_ddelta.*ddelta_dh;
    drho_dtheta0=sum(drho_ddelta.*ddelta_dtheta0,1);
    drho_dy0=sum(drho_ddelta.*ddelta_dy0,1);
    if ~strcmp(method,'MMC')
        dE=reshape(repmat(dE',nY,1),size(ddelta_c_dX,2),[])';
        drho_ddelta=reshape(repmat(drho_ddelta',nY,1),size(ddelta_c_dX,2),[])';
    end
    dE_dX=dE.*ddelta_c_dX;
    dE_dL=dE.*ddelta_c_dL; 
    drho_dX=drho_ddelta.*ddelta_dX;
    drho_dL=drho_ddelta.*ddelta_dL;
    xPhys=full(reshape(rho(:),nely,nelx));%
    E=full(reshape(E(:),nely,nelx));
    %passive elements
    xPhys(emptyelts) = 0;
    xPhys(fullelts) = 1;
    E(emptyelts) = p.Emin;
    E(fullelts) = p.E0;
    %% FE-ANALYSIS
    sK = reshape(KE(:)*(E(:)'),64*nelx*nely,1);
    K = sparse(iK,jK,sK); K = (K+K')/2;
    U(freedofs) = K(freedofs,freedofs)\F(freedofs);
    %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx);
    c = sum(sum((E).*ce));
    if mm==2
        [almc,dalmc,mtv,blc,~,~] = ALM_constraint1(Xg,nely,nY,np,p,BL);
    elseif mm==3
        [almc,dalmc,mtv, blc, dblc, mbl] = ALM_constraint(Xg,nely,nY,np,p,BL);
        dblc(:) = dblc(:).*(upper_bound(:)-lower_bound(:));
    end
    if mm > 1
        dalmc(:)=dalmc(:).*(upper_bound(:)-lower_bound(:));
    end
    v=mean(xPhys(:));
    dc_dE = -ce;
    dc_dE(emptyelts) = 0;
    dc_dE(fullelts) = 0;
    dc_dX=dE_dX*dc_dE(:);
    dc_dL=dE_dL*dc_dE(:);
    dc_dh=dE_dh*dc_dE(:);
    dc_dm=dE_dm*dc_dE(:);
    dc_dtheta0=dE_dtheta0*dc_dE(:);
    dc_dy0=dE_dy0*dc_dE(:);
    dc=zeros(size(X));
    dc(1:2:2*np*nY)=dc_dX;
    dc(2:2:2*np*nY)=dc_dL;
    dc((2*np*nY)+1:end-np-2)=dc_dh ;
    dc((2*np*nY)+np+1:end-2)=dc_dm;
    dc(end-1)=dc_dtheta0;
    dc(end)=dc_dy0;
    dv_dxPhys = ones(nely,nelx)/nelx/nely;
    dv_dxPhys(emptyelts) = 0;
    dv_dxPhys(fullelts) = 0;
    dv_dX=drho_dX*dv_dxPhys(:);
    dv_dL=drho_dL*dv_dxPhys(:);
    dv_dh=drho_dh*dv_dxPhys(:);
    dv_dm=drho_dm*dv_dxPhys(:);
    dv_dtheta0=drho_dtheta0*dv_dxPhys(:);
    dv_dy0=drho_dy0*dv_dxPhys(:);
    dv=zeros(size(X));
    dv(1:2:2*np*nY)=dv_dX;
    dv(2:2:2*np*nY)=dv_dL;
    dv((2*np*nY)+1:end-np-2)=dv_dh ;
    dv((2*np*nY)+np+1:end-2)=dv_dm;
    dv(end-1)=dv_dtheta0;
    dv(end)=dv_dy0;
    cvec(outit)=c;vvec(outit)=v;
    f0val=log(c+1);
    df0dx=(dc(:)/(c+1).*(upper_bound(:)-lower_bound(:)));
    if mm == 1
        fval =(v-volfrac)/volfrac*100;
        dfdx = dv(:)'/volfrac*100.*((upper_bound(:)-lower_bound(:))');
    elseif mm == 2
        fval=[(v-volfrac)/volfrac*100; almc];%
        dfdx=[dv(:)'/volfrac*100.*((upper_bound(:)-lower_bound(:))'); dalmc(:)'];
    else
        fval=[(v-volfrac)/volfrac*100; almc; blc];%
        dfdx=[dv(:)'/volfrac*100.*((upper_bound(:)-lower_bound(:))'); dalmc(:)'; dblc(:)'];
    end
%     if mm > 2
%         fval = [fval;blc];
%         dfdx = [dfdx; dblc(:)'];
%     end
    innerit=0;
    outvector1 = [outeriter innerit f0val fval'];
    outvector2 = xval;
    c0=c;
    x0=xPhys;
    %% Greyness Level
    gl = 4/nele*sum(xPhys(:).*(1-xPhys(:)));
    %% PRINT RESULTS
%     fprintf(f1, ' It.:%5i Obj.:%4.3e Vol.:%7.3f kktnorm.:%7.3f almc.:%7.3f  mtv.:%7.3f Theta: %5.2f Y0: %5.2f mbl: %3.2f ch.:%7.3f GL: %4.3f T: %7.2f\n',outit,c, ...
%         mean(xPhys(:)),kktnorm,almc,mtv,theta0/pi*180,y0,mbl,change, gl, toc);
    fprintf(' It.:%5i Obj.:%4.3e Vol.:%7.3f kktnorm.:%7.3f almc.:%7.3f  mtv.:%7.3f Theta: %5.2f Y0: %5.2f mbl: %3.2f ch.:%7.3f GL: %4.3f T: %7.2f\n',outit,c, ...
        mean(xPhys(:)),kktnorm,almc,mtv,theta0/pi*180,y0,mbl,change, gl, toc);
    %% Store variables
    storeXg(:,outit) = Xg;
%     storexPhys(:,outit) = xPhys(:);
    %% Plot Section
    if outit > 101
        plot_rate = 20;
    end
    if or(mod(outit,plot_rate)==0,outit == 1)
        figure(1)
%         hold on
        subplot(2,1,1)
        plot(1:outit,cvec(1:outit),'bo','MarkerFaceColor','b')
        title(['Convergence  Compliance =',num2str(c,'%4.3e'),', iter = ', num2str(outit)])
        grid on
        legend(['C =',num2str(c,'%4.3e'),', iter = ', num2str(outit)],'Location','Best')
        xlabel('iter')
        axis([0 outit 0 1e3])
        
        subplot(2,1,2)
        plot(1:outit,vvec(1:outit)*100,'ro','MarkerFaceColor','r')
        grid on
        legend(['V = ',num2str(mean(xPhys(:))*100),', iter = ', num2str(outit)],'Location','Best')
        xlabel('iter')
        title(['Convergence volfrac = ',num2str(mean(xPhys(:))*100),', iter = ', num2str(outit)])
        print([Path,'convergence_',num2str(outit,'%04d')],'-dpng')
        axis([0 outit 0 Inf])
%         hold off
        print([Path,image_prefix,'convergence'],'-dpng')
        %% PLOT DENSITIES
        figure(2)
        map=colormap(gray);
        map=map(end:-1:1,:);
        caxis([0 1])
        patchplot2 = patch('Vertices',[Xx,Yy],'Faces',edofMat(:,[2,4,6,8])/2,'FaceVertexCData',(1-xPhys(:))*[1 1 1],'FaceColor','flat','EdgeColor','none'); axis equal; axis off; ;hold on
        hold on
        fill([min(Xx),max(Xx),max(Xx),min(Xx)],[min(Yy),min(Yy),max(Yy),max(Yy)],'w','FaceAlpha',0.)
        scatter(Xx(fixednodes(fixed_dir==1)),Yy(fixednodes(fixed_dir==1)),'>b','filled')
        scatter(Xx(fixednodes(fixed_dir==2)),Yy(fixednodes(fixed_dir==2)),'^b','filled')
        scal=10;
        quiver(Xx(excitation_node),Yy(excitation_node)+scal*(excitation_direction==2),excitation_direction==1,-(excitation_direction==2),scal,'r','Linewidth',2)
        %     title('density plot')
        colormap(map)
        colorbar
        axis([min(Xx),max(Xx),min(Yy),max(Yy)])
        drawnow
        hold off
        print([Path,'density_',num2str(outit-1,'%04d')],'-dpng')
%%       Component Plot
        figure(3)
        Xk=Xg(1:2:2*np*nY);
        Lk=Xg(2:2:2*np*nY);
        m=Xg((2*np*nY)+np+1:end-2);
        h=Xg((2*np*nY)+1:end-np-2);
        theta0 = Xg(end-1);
        y0 = Xg(end);
        y1=linspace(0,sqrt(nelx^2 + nely^2),nY);
        y2=linspace(sqrt(nelx^2 + nely^2),0,nY);
        yn= repmat(reshape([y1 ;y2]',nY*2,1),1,np);
%         yn= reshape(yn,(nY)*2,1);
%         yn=repmat(yn,1,np);
%         theta0=0;
        lev_h1=fix(h(:)'*nY);
        lev_h2=ceil(h(:)'*nY);
        Xh=Xk(lev_h1+((1:np)-1)*nY).*(1-(h*nY-fix(nY*h)))+Xk(lev_h2+((1:np)-1)*nY).*(1-(ceil(nY*h)-h*nY));
        Lh=Lk(lev_h1+((1:np)-1)*nY).*(1-(h*nY-fix(nY*h)))+Lk(lev_h2+((1:np)-1)*nY).*(1-(ceil(nY*h)-h*nY));
%         Xh=Xk(lev_h1+((1:np)-1)*nY).*(1-(h*nY-fix(nY*h)))+Xk(lev_h2+((1:np)-1)*nY).*(1-(ceil(nY*h)-h*nY));
%         Lh=Lk(lev_h1+((1:np)-1)*nY).*(1-(h*nY-fix(nY*h)))+Lk(lev_h2+((1:np)-1)*nY).*(1-(ceil(nY*h)-h*nY));
        subsy=repmat(h(:)'*sqrt(nelx^2 + nely^2),2*nY,1);
        yt1=reshape(repmat(h(:)'*sqrt(nelx^2 + nely^2),nY,1),[],1);
        Y1=reshape(repmat(linspace(0,sqrt(nelx^2 + nely^2),nY)',1,np),[],1);
        subsX=reshape(repmat(Xh(:)',nY,1),[],1);
        subsL=reshape(repmat(Lh(:)',nY,1),[],1);
        Lk(Y1>=yt1)=subsL(Y1>=yt1);
        Xk(Y1>=yt1)=subsX(Y1>=yt1);
        yn(yn>=subsy)=subsy(yn>=subsy);
        xn1=reshape((Xk - Lk./2),nY,np);
        xn_inv=reshape((Xk + Lk./2),nY,np);
        xn_inv=flip(xn_inv);
        xn=cat(1,xn1,xn_inv);
        x_comp = nelx/2+cos(theta0)*(xn-nelx/2)-sin(theta0)*(yn-nely/2+y0); 
        y_comp = nely/2+sin(theta0)*(xn-nelx/2)+cos(theta0)*(yn-nely/2+y0); 
%         xn=xn(:);
%        [dd]=norato_bar(xi,Eta,repmat(Lc(:),1,size(cc,2)),repmat(Yc(:),1,size(cc,2)));
%         xn=repmat(Xc,1,size(cc,2))+dd.*cc;
%         yn=repmat(Yc,1,size(cc,2))+dd.*ss;
        tolshow=0.1;
        Shown_compo=find(m>tolshow);
        fill([min(Xx),max(Xx),max(Xx),min(Xx)],[min(Yy),min(Yy),max(Yy),max(Yy)],'w','FaceAlpha',0.)
        hold on
        fill(x_comp(:,Shown_compo),y_comp(:,Shown_compo),m(Shown_compo),'FaceAlpha',0.5)
        if strcmp(BC,'L-shape')
            fill([fix((min(Xx)+max(Xx))/2),max(Xx),max(Xx),fix((min(Xx)+max(Xx))/2)],[fix((min(Yy)+max(Yy))/2),fix((min(Yy)+max(Yy))/2),max(Yy),max(Yy)],'w')
        end
        %     fill([X1(Shown_compo),X2(Shown_compo),X3(Shown_compo),X4(Shown_compo)]',[Y1(Shown_compo),Y2(Shown_compo),Y3(Shown_compo),Y4(Shown_compo)]',[Mc(Shown_compo)],'FaceAlpha',0.5)
        caxis([0,1])
        colormap 'jet'
        axis equal; axis off;
        hold on
        scatter(Xx(fixednodes(fixed_dir==1)),Yy(fixednodes(fixed_dir==1)),'>b','filled')
        scatter(Xx(fixednodes(fixed_dir==2)),Yy(fixednodes(fixed_dir==2)),'^b','filled')
        scal=10;
        quiver(Xx(excitation_node),Yy(excitation_node)+scal*(excitation_direction==2),excitation_direction==1,-(excitation_direction==2),scal,'r','Linewidth',2)
        colorbar
        axis([min(Xx),max(Xx),min(Yy),max(Yy)])
        title('component plot')
        drawnow
        print([Path,'component_',num2str(outit-1,'%04d')],'-dpng')
        hold off
    end
    %% MMA code optimization
    [X,ymma,zmma,lam,xsi,eta,mu,zet,S,low,upp] = ...
        mmasub(mm,n,outeriter,xval,xmin,xmax,xold1,xold2, ...
        f0val,df0dx,fval,dfdx,low,upp,a0,a,C,d);
    xold2 = xold1;
    xold1 = xval;
    xval  = X;
    Xg=lower_bound+(upper_bound-lower_bound).*X;
    change=norm(xval-xold1);
    %% %% The residual vector of the KKT conditions is calculated:
    [residu,kktnorm,residumax] = ...
        kktcheck(mm,n,X,ymma,zmma,lam,xsi,eta,mu,zet,S, ...
        xmin,xmax,df0dx,fval,dfdx,a0,a,C,d);
    outvector1 = [outeriter innerit f0val fval(:)'];
    outvector2 = xval;
    switch stopping_criteria
        case 'kktnorm'
            stop_cond=outit < maxoutit && kktnorm>kkttol;
        case 'change'
            stop_cond=outit < maxoutit &&change>changetol;
    end
     %
% figure(2)   %%% essai plot
%     fill([min(Xx),max(Xx),max(Xx),min(Xx)],[min(Yy),min(Yy),max(Yy),max(Yy)],'w','FaceAlpha',0.)
%         hold on
%     plotConfig([Xg],Yc,np,nY,nely); axis equal; axis off; drawnow; axis([0 nelx 0 nely]);
%     scatter(Xx(fixednodes(fixed_dir==1)),Yy(fixednodes(fixed_dir==1)),'>b','filled')
%     scatter(Xx(fixednodes(fixed_dir==2)),Yy(fixednodes(fixed_dir==2)),'^b','filled')
%     if strcmp(BC,'L-shape')
%             fill([fix((min(Xx)+max(Xx))/2),max(Xx),max(Xx),fix((min(Xx)+max(Xx))/2)],[fix((min(Yy)+max(Yy))/2),fix((min(Yy)+max(Yy))/2),max(Yy),max(Yy)],'w')
%     end
%     scal=10;
%     quiver(Xx(excitation_node),Yy(excitation_node)+scal*(excitation_direction==2),excitation_direction==1,-(excitation_direction==2),scal,'r','Linewidth',2)
%     print([Path,'component_',num2str(outit,'%03d')],'-dpng')
%     drawnow
%         axis([min(Xx),max(Xx),min(Yy),max(Yy)])
% hold off
%     %
    if toc/3600 > 23
        save([image_prefix,'.mat']);
    end
    
end
%% PLOT DENSITIES
figure(1)
map=colormap(gray);
map=map(end:-1:1,:);
caxis([0 1])
patchplot2 = patch('Vertices',[Xx,Yy],'Faces',edofMat(:,[2,4,6,8])/2,'FaceVertexCData',(1-xPhys(:))*[1 1 1],'FaceColor','flat','EdgeColor','none'); axis equal; axis off; ;hold on
hold on
fill([min(Xx),max(Xx),max(Xx),min(Xx)],[min(Yy),min(Yy),max(Yy),max(Yy)],'w','FaceAlpha',0.)
scatter(Xx(fixednodes(fixed_dir==1)),Yy(fixednodes(fixed_dir==1)),'>b','filled')
scatter(Xx(fixednodes(fixed_dir==2)),Yy(fixednodes(fixed_dir==2)),'^b','filled')
scal=10;
quiver(Xx(excitation_node),Yy(excitation_node)+scal*(excitation_direction==2),excitation_direction==1,-(excitation_direction==2),scal,'r','Linewidth',2)
%     title('density plot')
colormap(map)
colorbar
drawnow
hold off
axis([min(Xx),max(Xx),min(Yy),max(Yy)])
FEM_enhanceDisplay
print([Path,'density_',num2str(outit,'%04d')],'-dpng')
fprintf(f1,'Elapsed Time: %10.3f',toc);
fclose(f1);
matname = fullfile(folder_name, 'desvar.mat');
save(matname,'storeXg');
end
%% Component Plot
%         figure(2)
%         Xk=Xg(1:2:2*np*nY);
%         Lk=Xg(2:2:2*np*nY);
%         m=Xg((2*np*nY)+np+1:end);
%         h=Xg((2*np*nY)+1:end-np);
%         y1=linspace(0,nely,nY);
%         y2=linspace(nely,0,nY);
%         yn= [y1 ;y2]';
%         yn= reshape(yn,(nY)*2,1);
%         yn=repmat(yn,1,np);
%         
%         lev_h1=fix(h(:)'*nY);
%         lev_h2=ceil(h(:)'*nY);
%         Xh=Xk(lev_h1+((1:np)-1)*nY).*(1-(h*nY-fix(nY*h)))+Xk(lev_h2+((1:np)-1)*nY).*(1-(ceil(nY*h)-h*nY));
%         Lh=Lk(lev_h1+((1:np)-1)*nY).*(1-(h*nY-fix(nY*h)))+Lk(lev_h2+((1:np)-1)*nY).*(1-(ceil(nY*h)-h*nY));
%         subsy=repmat(h(:)'*nely,2*nY,1);
%         yt1=reshape(repmat(h(:)'*nely,nY,1),[],1);
%         Y1=reshape(repmat(linspace(0,nely,nY)',1,np),[],1);
%         subsX=reshape(repmat(Xh(:)',nY,1),[],1);
%         subsL=reshape(repmat(Lh(:)',nY,1),[],1);
%         Lk(Y1>=yt1)=subsL(Y1>=yt1);
%         Xk(Y1>=yt1)=subsX(Y1>=yt1);
%         yn(yn>=subsy)=subsy(yn>=subsy);
%         xn1=reshape((Xk - Lk./2),nY,np);
%         xn_inv=reshape((Xk + Lk./2),nY,np);
%         xn_inv=flip(xn_inv);
%         xn=cat(1,xn1,xn_inv);
%          
% %         xn=xn(:);
% %        [dd]=norato_bar(xi,Eta,repmat(Lc(:),1,size(cc,2)),repmat(Yc(:),1,size(cc,2)));
% %         xn=repmat(Xc,1,size(cc,2))+dd.*cc;
% %         yn=repmat(Yc,1,size(cc,2))+dd.*ss;
%         tolshow=0.1;
%         Shown_compo=find(m>tolshow);
%         fill([min(Xx),max(Xx),max(Xx),min(Xx)],[min(Yy),min(Yy),max(Yy),max(Yy)],'w','FaceAlpha',0.)
%         hold on
%         fill(xn(:,Shown_compo),yn(:,Shown_compo),m(Shown_compo),'FaceAlpha',0.5)
%         if strcmp(BC,'L-shape')
%             fill([fix((min(Xx)+max(Xx))/2),max(Xx),max(Xx),fix((min(Xx)+max(Xx))/2)],[fix((min(Yy)+max(Yy))/2),fix((min(Yy)+max(Yy))/2),max(Yy),max(Yy)],'w')
%         end
%         %     fill([X1(Shown_compo),X2(Shown_compo),X3(Shown_compo),X4(Shown_compo)]',[Y1(Shown_compo),Y2(Shown_compo),Y3(Shown_compo),Y4(Shown_compo)]',[Mc(Shown_compo)],'FaceAlpha',0.5)
%         caxis([0,1])
%         colormap 'jet'
%         axis equal; axis off;
%         hold on
%         scatter(Xx(fixednodes(fixed_dir==1)),Yy(fixednodes(fixed_dir==1)),'>b','filled')
%         scatter(Xx(fixednodes(fixed_dir==2)),Yy(fixednodes(fixed_dir==2)),'^b','filled')
%         scal=10;
%         quiver(Xx(excitation_node),Yy(excitation_node)+scal*(excitation_direction==2),excitation_direction==1,-(excitation_direction==2),scal,'r','Linewidth',2)
%         colorbar
%         %     title('component plot')
%         drawnow
%         axis([min(Xx),max(Xx),min(Yy),max(Yy)])
%         %     print([Path,'component_',num2str(outit-1,'%03d')],'-dpng')
%         hold off
% %     title('component plot')
% drawnow
% axis([min(Xx),max(Xx),min(Yy),max(Yy)])
% FEM_enhanceDisplay
% print([Path,'component_',num2str(outit,'%03d')],'-dpng')
% hold off
% figure(3)
% subplot(2,1,1)
% plot(1:outit,cvec(1:outit),'bo','MarkerFaceColor','b')
% % title(['Convergence  Compliance =',num2str(c,'%4.3e'),', iter = ', num2str(outit)])
% grid on
% hold on
% scatter(outit,c,'k','fill')
% hold off
% text(outit,c,['C =',num2str(c,'%4.2f'),' at iteration ', num2str(outit)],'VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10,'FontWeight','bold')
% % legend(['C =',num2str(c,'%4.3e'),', iter = ', num2str(outit)],'Location','Best')
% xlabel('iter')
% ylabel('C')
% FEM_enhanceDisplay
% % axis([0 outit 0 1e3])
% subplot(2,1,2)
% plot(1:outit,vvec(1:outit)*100,'ro','MarkerFaceColor','r')
% grid on
% hold on
% scatter(outit,mean(xPhys(:))*100,'k','fill')
% hold off
% text(outit,mean(xPhys(:))*100,['V = ',num2str(mean(xPhys(:))*100,'%4.2f'),'% at iteration ', num2str(outit)],'VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10,'FontWeight','bold')
% % legend(['V = ',num2str(mean(xPhys(:))*100),', iter = ', num2str(outit)],'Location','Best')
% xlabel('iter')
% ylabel('V [%]')
% FEM_enhanceDisplay
% % axis([0 outit 0 Inf])
% % title(['Convergence volfrac = ',num2str(mean(xPhys(:))*100),', iter = ', num2str(outit)])
% print([Path,'convergence_',num2str(outit,'%03d')],'-dpng')

