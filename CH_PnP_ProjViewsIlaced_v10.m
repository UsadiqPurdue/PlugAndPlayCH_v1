
function []=CH_PnP_ProjViewsIlaced_v10(Ang,r,K,FinSF,MINTH,MAXTH,ExpCat,JFunc,gamma,sig_H,mode,varargin)

% v7 runs with random initial conditions
% _v6 does worse FBP by increasing r and cleans up some code

% Output in: 'Denoise_AL_Poly_For_PnP_120_angles_full_den';
% Recon works partially
% Includes Poly, Does PnP after denoising for 5 iterations(no crit on reg)
% Log for cost
% Some data also saved in 'Denoise_AL_Poly_For_PnP_120_angles_full_den'

switch mode
    case 8
        sig_H=varargin{1};
    case 2
        tvmu=varargin{1};
end

%% Generate Cahn-Hilliard Data
clc;
close all;

randn('seed',5);
% Variables: CHANGE NAME AT EACH RUN

DenAlgos={'BM4D','TV3D','','','','','','Cahn-Hilliard'};

Nt=r*K*FinSF;

s=128;

s=256;

ss_rate=r;

rr=ceil(sqrt(2*s^2));

n=Nt+1;

iter=30;
name='AL_correct_theta';
a=200e4;b=25e4;

a=800/50;b=100/50;
sig=6e-2;


lam0=1e5;
lam0=1e8;

a=800/3;b=100/3;
%% Input and Output Folders


RDir=strcat('/Users/USadiq/Results/CHPhantomDataResults/',ExpCat,'/');


mkdir(RDir);

DestFName=strcat('IViewsBRev_Ang',num2str(Ang),'r',num2str(r),'K',num2str(K),'Nt',num2str(n),'MaxTH',num2str(MAXTH),'gamma',num2str(gamma),'sig_H',num2str(sig_H));


%% Initialize Data

switch mode
    case 8
        ResultsFolder=strcat(RDir,'/',date,'_',char(DenAlgos(mode)),'/',DestFName,'/SigH_',num2str(sig_H),'/');
    case 2
 ResultsFolder=strcat(RDir,'/',date,'_',char(DenAlgos(mode)),'/',DestFName,'/Mu_',num2str(tvmu),'/');
    case 1
        ResultsFolder=strcat(RDir,'/',date,'_',char(DenAlgos(mode)),'/',DestFName);
end

mkdir(ResultsFolder);

SimRunOn=date;

save(strcat(ResultsFolder,'InputParams.mat'));


%% Generate A Matrix


theta=Ang;tt=Ang;beta=0.002;


m=rr*theta;

Ax=zeros(rr,theta);



%% Generate Phantom

%CH_gen_x_fft_recon_01_circ256_subsample(s,ss_rate*n,sig,a/ss_rate,b/ss_rate,15,5,ss_rate);

[img_in,img_orig]=CH_gen_x_fft_recon_01_circ256_subsample_v9(s,n,sig,a,b,0,4,ss_rate);
%n=n-1;



daten=date;
diary(strcat('PnP_180_est_theta_',daten));

ups_rate=1;
%% Initialize sino3, geom3, Amatrix 

             

[SubFrame,UpSampledMap]=CHPlugAndPlayInitProgViewsILaced(img_orig,ResultsFolder,rr,theta,lam0,ups_rate,K,FinSF,MINTH,MAXTH);

%[sino3,map_image3,params3,geom,Amatrix,AmatrixZ,d3]=PlugAndPlayInitCustomGeomSubsProgViewsMultiK(img_orig,rr,theta,lam0,ups_rate,K,FinSF);


[SubFrame]=SetKParamsIL(SubFrame,mode);


%% Do reconstruction using Plug-n-Play
fprintf(strcat('Saving results in ',ResultsFolder));


% For doing without theta estimation Call ADMM_Core3_CH_01_Poly_Java_SubSample_varss_r3

close all;

%[map_image3,params3,FwdPrj,denObj,diffxv,muAll]=ADMM_Core3_CH_v10(img_orig,SubFrame,ResultsFolder,FinSF,JFunc,gamma,sig_H,a*r,b*r,mode);

switch mode
    case 2
       [~,~,x_hat,v_hat,Admm_d,Admm_h]=ADMM_Core3_CH_v11(img_orig,SubFrame,ResultsFolder,FinSF,JFunc,sig,sig_H,a,b,mode,15,tvmu);    
    case 1
       [~,~,x_hat,v_hat,Admm_d,Admm_h]=ADMM_Core3_CH_v11(img_orig,SubFrame,ResultsFolder,FinSF,JFunc,sig,sig_H,a,b,mode,15);
end
1;


