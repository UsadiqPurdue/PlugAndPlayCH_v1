function [SubFrame]=SetKParamsIL(SubFrame,mode)

load('AL_Data.mat');%% Denoising parameter for ALL algorithms

params=SubFrame(1).Frames(1).Settings.params;

for sf=1:length(SubFrame)
    for f=1:length(SubFrame(sf).Frames)
        
    kparams.sigma=sqrt(params.beta/params.lambda);
    %10.6066;%10.6066;
    
    %% Denoising parameters - KSVD
    if(mode == 0 || mode == 1)
        %kparams.x=(map_image./max(map_image(:))).*255;
        kparams.blocksize = 4;
        kparams.dictsize =256;
        kparams.maxval=255;
        kparams.trainnum = 3600;%40000;
        kparams.iternum =10;
        kparams.memusage='high';
    end
    %% These are only for TV denoising
    if(mode == 2)
        kparams.verb = 0;
        kparams.display = 0;
        kparams.niter = 100;    % number of iterations
        kparams.c_TV = .0998;
        %kparams.etgt = (kparams.sigma/255)*(max(m,n));
        kparams.lambda = kparams.c_TV*kparams.sigma^2;%2./kparams.sigma^2; % initial regularization
        kparams.lambda = 0.25;
    end
    %% qGGMRF
    if(mode == 4)
        kparams.p=2;
        kparams.q=1.2;
        kparams.c=1e-2;
        kparams.sigmax= 0.29;%0.5940*(params.lambda^(1/kparams.p));
        kparams.niter = 10; % CHANGE THIS TO AVOID CONFLICT
        kparams.verbose = 0;
        display('Bad variable names for num iter - denoising!');
    end
    %% General Segmentation parameter
    if(mode ==7)
        kparams.num_class=5;
        
        %% MAP Segmentation
        kparams.max_iter=10;
        kparams.filter = params.filter;
        kparams.beta = 4;
        kparams.debug = 0;
        kparams.sigma_sq = 1/(params.lambda);
        kparams.rand_ord = 0;%Regular ICM or random order ICM
    end
    
    SubFrame(sf).Frames(f).Settings.kparams=kparams;
end
end