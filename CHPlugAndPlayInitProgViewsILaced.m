function [SubFrame,map_image3,MapImag2UpSample2,map_image_orig]=CHPlugAndPlayInitProgViewsILaced(map_Object,foldername,rr,n_theta,lam0,ups_rate,K,FinSF,MINTH,MAXTH)

[m n ll]=size(map_Object);

stepd=(MAXTH -  MINTH)/((FinSF-1)/FinSF+(n_theta-1));

dt=stepd/FinSF;


maxtheta=(n_theta-1)/2;

theta_u=MAXTH;

geom.n_x=m;
geom.n_y=n;
geom.x_0=-m/2;
geom.y_0=-n/2;
geom.delta=1;


% [sino,ANG_D,ANG_R]=SetSinoGeom(MAXTH,MINTH,OFFSET,n_theta,rr);
%     
% Amatrix=ComputeAMatrix(sino,geom);
% AmatrixZ=AMatrixZero(sino,geom);

params.max_iter = 100; %max outer iterations for ADMM
params.threshold = 5e-6;% Stopping criteria for algorithm
%% Other paramters - ICD

params.num_iter=5;
params.filter=[1/12 1/6 1/12;1/6 0 1/6;1/12 1/6 1/12];
params.lambda = 1/20;%1/75;%Lagrange multiplier
params.u = zeros(m,n); %Augmented Lagrange vector
params.v   = zeros(m,n); %Auxilarry variable
params.verbose=1;
params.beta = 1.50;
params.xray=0;
if(params.xray) %If this is true, we should have a dose value stored
    params.dose = Dose;
end

% No of subframes=K

sf=0;

f=0;
rad=floor(m/2);

THETA.MAXTH=MAXTH;
THETA.MINTH=MINTH;
THETA.n_theta=n_theta;

SFOrder='bitrev';

bitorder=bitrevorder([0:FinSF-1]);
RRecon=1;
count=1;

for i=1:ups_rate:ll

 
        image_i(:,:)=map_Object(:,:,i);
        
        if(strcmp(SFOrder,'bitrev'))
            offset=-bitorder(mod(f,FinSF)+1)*dt;
        else
            
            offset=-rem(f,FinSF)*dt;
        end
        %% Forward Project and add noise
        
        [sino,ANG_D,ANG_R]=SetSinoGeom(MAXTH,MINTH,offset,rr,stepd);
        
        Amatrix=ComputeAMatrix(sino,geom);
        
        Ax=forward_project_v2(image_i,sino,Amatrix);
        
        lam=lam0*exp(-Ax);
    
        l=abs(lam+sqrt(lam).*randn(n_theta,rr));
        
        
        sinocounts=-log(l./lam0);
        
        sino.counts=sinocounts;
        
        SubFrame(sf+1).Frames(f+1).Angles=ANG_D;
        SubFrame(sf+1).Frames(f+1).Projs=sino;
        %SubFrame(sf+1).Frames(f+1).AMatrix=Amatrix;
        SubFrame(sf+1).Frames(f+1).D=l;
        
        AMat=strcat('A_SF',num2str(sf+1),'_F',num2str(f+1));
        
        save(strcat(foldername,AMat,'.mat'),'Amatrix');
        
%         fprintf('\tAngles:');
%         fprintf('%d',ANG_D);
%         fprintf('\n');
%         fprintf('\tAMatrix:%s\n',AMat);
        
        % Write Parameters for Frame
        
        SubFrame(sf+1).Frames(f+1).Settings.params=params;
        
        SubFrame(sf+1).Frames(f+1).Settings.params.vorig=iradon(sinocounts',ANG_D,'cubic','Hamming',m);
        
        SubFrame(sf+1).Frames(f+1).Settings.params.v=circ_section(SubFrame(sf+1).Frames(f+1).Settings.params.vorig,m,n,rad-6,0.08); 
        
        %SubFrame(sf+1).Frames(f+1).Settings.params.v=iradon(sinocounts',ANG_D,'cubic','Hamming',m);
        
        SubFrame(sf+1).Frames(f+1).Settings.geom=geom;
        
        f=f+1;
        
        map_image_orig(:,:,count)=iradon(sinocounts',ANG_D,'cubic','Hamming',m);
        
        count=count+1;
        %% Interleave sinograms to construct a subframe
        
        if(f==FinSF)       
            [SubFrame(sf+1).SFRecon.map_image3(:,:,:),SubFrame(sf+1).SFRecon.sino3(1,1),SubFrame(sf+1).SFRecon.d3(:,:,:),SubFrame(sf+1).SFRecon.ANG_D(:,:,:)]=InterleaveRecon(SubFrame(sf+1),FinSF,THETA,rr,K,m);

            ImgNewCode(:,:,:)=InterleaveReconCustom(SubFrame(sf+1),FinSF,RRecon,THETA,rr,K,m);
            
            ReconsInSF=FinSF/RRecon;
            
            ImgNew((sf*ReconsInSF)+1:(sf+1)*ReconsInSF,:,:)=ImgNewCode(:,:,:);
            
            f=0;
            sf=sf+1;

        end
        
        
  
    
end

% Data must have:
% map_image3,sino3,geom,Amatrix

%map_image3=(map_image3-min(min(min(map_image3))))./(max(max(max(map_image3)))-min(min(min(map_image3))));


count=1;

for j=1:K
        MapImag2UpSample(:,:,j)=SubFrame(j).SFRecon.map_image3;
end

kk=size(ImgNew);

for j=1:kk
        MapImag2UpSample2(:,:,j)=ImgNew(j,:,:);
end

if(K>1) % Make sure there are at least two data points
    
    %ImgoutUpsample=ImageInterp3D(MapImag2UpSample,ll);
    map_image3(:,:,1:ll)=ImageInterp3D(MapImag2UpSample,ll/K);
    %map_image3(:,:,1:sf*ups_rate)=ImgoutUpsample(:,:,1:sf*ups_rate);
else
    display('Upsampling not possible');
    map_image3=MapImag2UpSample;
end


figure;imagesc(SubFrame(1).SFRecon.map_image3);colormap(gray);colorbar;title('Actual FBP');
figure;imagesc(SubFrame(1).Frames(1).Settings.params.v);colormap(gray);colorbar;title('My FBP');

end