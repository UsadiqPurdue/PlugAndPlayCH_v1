function [u,Err,TV,Fwd,Cost]=ATV_ROF_3D_2(f,mu,lambda,Niter)

%=================================================================
% function u=ATV_ROF_3D(f,mu,lambda,Niter)
%
% 3D Anisotropic ROF denoising
% Version:
% -v1.2 - 06/20/3013
% -v1.0 - 02/27/2011
%
% This function performs the minimization of the Isotropic TV - L2 
% Rudin-Osher-Fatemi model for 3D data by Split Bregman Iterations
%
% f = noisy 3D datacube
% mu = regularization parameter
% lambda = regularization parameter for the Bregman variable
% Niter = maximum number of iteration
%
% Author: Jerome Gilles
% Institution: UCLA - Math Department
% email: jegilles@math.ucla.edu
%
%=================================================================

[M,N,P]=size(f);

f=double(f);
dx=zeros(M,N,P);
dy=zeros(M,N,P);
dz=zeros(M,N,P);
bx=zeros(M,N,P);
by=zeros(M,N,P);
bz=zeros(M,N,P);
u=f;
Z=zeros(M,N,P);

Mask=zeros(M,N,P);
Mask(1,1,1) = 1;
FMask=fftn(Mask);

%Fourier Laplacian mask initialization
D = zeros(M,N,P);
Lap=zeros(3,3,3);

Lap(:,:,1)=[0,0,0;0, 1,0;0,0,0];
Lap(:,:,2)=[0,1,0;1,-6,1;0,1,0];
Lap(:,:,3)=[0,0,0;0, 1,0;0,0,0];

D([end,1,2],[end,1,2],[end,1,2]) = Lap;
clear Lap;
FD=fftn(D);

%Fourier constant initialization
FW=((mu/lambda)*abs(FMask).^2-real(FD)).^-1;
FF=(mu/lambda)*conj(FMask).*fftn(f);

err=sum(sum(sum(f.^2)));
tol=1e-3*err;
tol=1e-8*err;
K=0;
while ((err>tol) && (K<Niter)),
    K=K+1;
    tx=dx-bx;
    ty=dy-by;
    tz=dz-bz;
    
    up=u;
    %Update u
    u=real(ifftn(FW.*(FF-fftn(tx-tx(:,[1,1:N-1],:)+ty-ty([1,1:M-1],:,:)+tz-tz(:,:,[1,1:P-1])))));
    ux=u-u(:,[1,1:N-1],:);
    uy=u-u([1,1:M-1],:,:);
    uz=u-u(:,:,[1,1:P-1]);
        
    %Update dx
    tmpx=ux+bx;
    dx=sign(tmpx).*max(Z,abs(tmpx)-1/lambda);
    
    %Update dy
    tmpy=uy+by;
    dy=sign(tmpy).*max(Z,abs(tmpy)-1/lambda);

    %Update dz
    tmpz=uz+bz;
    dz=sign(tmpz).*max(Z,abs(tmpz)-1/lambda);

    %Update bx and by and bz
    bx=tmpx-dx;
    by=tmpy-dy;
    bz=tmpz-dz;
   
    err=sum(sum(sum((up-u).^2)));
    
    Err(K)=err;
    
    TV(K)=TV3D(u);
    
    Fwd(K)=sum(sum(sum((u-f).^2)));
    
    Cost(K)=mu/lambda*Fwd(K)+TV(K);
    
end