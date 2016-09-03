function [map_image3,params3,x_hat,v_hat,Admm_costd,tv,diffxv,muAll,cost_all_ICD,...
    cost_all_H_ICD,cost_all_f_ICD,ICDiters]=ADMM_Core3_CH_v11(img_orig3,...
    SubFrame,foldername,FinSF,JFunc,sig,sig_H,a,b,mode,NIterPnP,varargin)

% Compares with QGGMRF and other denoising techniques

% V5: Code management: Introduced new structures.


%Performs ADMM based minimization of min_{x,v} (1/2)||y-Ax||^2 + s(v)
% Inputs
% map_image is an inital estimate for x
% sino is a structure which contains all the sinogram information
% geom contains the geometry information like number of pixels etc.
% Amatrix = A stored as a sparse matrix
% params contains parameters associated with the optimization
% kparams contains prior model parameters
% mode represents what prior is used in s(v)
%    0 - kSVD
%    1 - BM3D
%    2 - TV
%    3 - PLOW
%    4 - qGGMRF
%    5 - Otsu's segmentation
%    6 - SMAP segmentation %DOES NOT WORK 5/4/2013
%    7 - MAP segmentation


iterdn=3;


drct=cd;

mkdir(foldername);


%% Naive initialization

PnPIters=25; XVDiffThresh=5e-5;

if(length(varargin)>1)
    TV_mu=varargin{1};
end

TV_lambda=100;TV_it=20;

params=SubFrame(1).Frames(1).Settings.params;

kparams=SubFrame(1).Frames(1).Settings.kparams;

tempf=SubFrame.Frames;

l=length(SubFrame)*length(tempf);

for ll=1:l
    
    SFNo=ceil(ll/FinSF);
    FNo=rem(ll-1,FinSF)+1;
    
    map_image3(:,:,ll)=SubFrame(SFNo).Frames(FNo).Settings.params.v;
    sino3(ll,1)=SubFrame(SFNo).Frames(FNo).Projs;
    
    Angs(ll,:)=SubFrame(SFNo).Frames(FNo).Angles;
    kparams3(ll,1)=SubFrame(SFNo).Frames(FNo).Settings.kparams;
    params3(ll,1)=SubFrame(SFNo).Frames(FNo).Settings.params;
    geom(ll,1)=SubFrame(SFNo).Frames(FNo).Settings.geom;
    
    d3(ll,:,:)=SubFrame(SFNo).Frames(FNo).D;
    
end



[m n l]=size(map_image3);

iter=1;

%% END OF DEBUG STATEMENTS

params.max_iter=PnPIters;

XVDiff=Inf;

RMSE=4500;

ICDIters=0;HUpds=0;


for PnPit=1:NIterPnP% && stop_crit > params.threshold)
    
    fprintf('PnP Iteration %d',iter);
    
    
    if(iter > 1) %After the first outer iteration just do 1 iteration
        %    params.num_iter=1;
        
        if(isfield(kparams,'initdict') && ischar(kparams.initdict))
            kparams.initdict=dict; %for kSVD use the previous dictionary to start
        end
    end
    
    %% Forward Projection block
    
    
    for ll=1:l
        
        if(iter==1)
            FBPImg(:,:,ll)=map_image3(:,:,ll);
        end
        
        SFNo=ceil(ll/FinSF);
        FNo=rem(ll-1,FinSF)+1;
        
        map_image=map_image3(:,:,ll);
        
        sino=sino3(ll,1);
        kparams=kparams3(ll,1);
        params=params3(ll,1);
        
        params.lambda=1/(2*sig^2);
        
        x_til(iter,:,:,ll)=params.v-params.u;
        
        %         if(rem(ll,ups_rate)==1)
        %            map_image3(:,:,ll) = DataTermOpt(map_image,sino,d3(ll,:,:),geom,params,Amatrix);
        %            1;
        %         else
        %            map_image3(:,:,ll) = DataTermOpt(map_image,sino,d3(ll,:,:),geom,params,AmatrixZ);
        %         end
        
        AMat=strcat('A_SF',num2str(SFNo),'_F',num2str(FNo));
        
        load(strcat(foldername,AMat,'.mat'));
        
        %         fprintf('\nAngles:%d\t',Angs(ll));
        %         fprintf('\nAMatrix:%s\t',AMat);
        
        map_image3(:,:,ll) = DataTermOpt(map_image,sino,(1/(2*RMSE^2))*d3(ll,:,:),geom,params,Amatrix);
        
        %[imout,err,tv,lambda] = perform_tv_denoising(kparams.x,kparams);
        
        
        
        clear(strcat(foldername,AMat,'.mat'));
        
        x_hat(iter,:,:,ll)=map_image3(:,:,ll);
        
        
        [~,~,costd(ll)]= ADMM_data_cost(map_image,sino,d3(ll,:,:),params,Amatrix);
        
        options.null=0;
        
        
        
    end
    
    
    save(strcat(foldername,'x_til.mat'),'x_til');
    
    save(strcat(foldername,'x_hat.mat'),'x_hat');
    
    %% De-noising
    
    for ll=1:l
        
        kparams3(ll,1).x = map_image3(:,:,ll) + params3(ll,1).u;
        
        img_in(:,:,ll)=map_image3(:,:,ll) + params3(ll,1).u;
        
        kparams=kparams3(ll,1);
        
        kparams.maxval=1.0;
        
        v_tilda(iter,:,:,ll)=img_in(:,:,ll);
        
        switch mode
            
            case 4
                imout3(:,:,ll) = qGGMRFdenoise(img_in(:,:,ll),kparams);
                
        end
        
        
    end
    
    switch mode
        
        case 2
            
            tv(iter)=TV3D(map_image3);
            Admm_costd(iter)=sum(sum(costd));
            
            [imout3,TVTmp]=ATV_ROF_3D_2(img_in,TV_mu,TV_lambda,TV_it);
            TVDen(iter,1:length(TVTmp))=TVTmp;
            
            imout3(imout3>1)=1.0;
            imout3(imout3<0)=0.0;
            
        case 1
            
            imout3 = bm4d(img_in,'Gauss',0,'lc',0,0);
            tv=[];
            Admm_costd(iter)=sum(sum(costd));
            
        case 8
            
            H00_poly=denoise_cost_Hf_impl_1_0528_01(map_image3,s,n,a,b);
            Hc=sum(sum(sum(H00_poly.^2)));
            costHx(iter)=(1/100)*sqrt(Hc/(m*n*l));
            
            Admm_costd(iter)=sum(sum(costd));
            tv(iter)=costHx(iter);
            
            [imout3,cost,costH,costF]=...
                Denoise_AL_Poly_For_PnP_Java128_v3_HUpdate_Costs_v10...
                (JFunc,1,img_in,img_orig3,m,l,iterdn,sig,sig_H,foldername,...
                a,b,300,ICDIters,HUpds);
            
            
    end
    
    save(strcat(foldername,'v_til.mat'),'v_tilda');
    
    save(strcat(foldername,'PnP_InitData.mat'));
    
    fprintf('Denoising iteration %d...\n',iter);
    v_hat(iter,:,:,:)=imout3;
    save(strcat(foldername,'v_hat.mat'),'v_hat');
    
    %% PnP Continues
    
    for ll=1:l
        prev_v3(:,:,ll)=params3(ll,1).v;
        params3(ll,1).v = imout3(:,:,ll);
        
        vv=params3(ll,1).v;
        params3(ll,1).u = params3(ll,1).u+(map_image3(:,:,ll)-vv);
        
        diffxv(iter,ll,:,:)=(map_image3(:,:,ll)-vv);
        
        mu(iter,:,:,ll)=params3(ll,1).u;
        
    end
    
    save(strcat(foldername,'mu.mat'),'mu');
    
    eps_primal=0;
    eps_dual=0;
    
    for ll=1:l
        prev_v=prev_v3(:,:,ll);
        params=params3(ll,1);
        map_image=map_image3(:,:,ll);
        eps_primal = eps_primal+sum(abs(params.v(:)-map_image(:)))./sum(abs(map_image(:)));
        eps_dual = eps_dual+sum(abs(params.v(:)-prev_v(:)))./sum((abs(prev_v(:))));
        
    end
    
    eps_dual = sum(abs(params.v(:)-prev_v(:)))./sum((abs(prev_v(:))));
    
    stop_crit = (eps_primal+eps_dual)/2;
    
    XVDiff=sum(sum(sum((x_hat(iter,:,:,:)-v_hat(iter,:,:,:)).^2)));
    
    XVRes(iter)=XVDiff;
    
    iter=iter+1;
    
    params.max_iter=PnPIters;
    
    
    
end

save(strcat(foldername,'PnPCosts.mat'),'tv','Admm_costd','TVDen');

%PlotPnPCost(foldername);

%GenDenoisingCostsPlot(foldername);

%MakeVideo(img_orig3,map_image3,foldername);

diary off;

1;
