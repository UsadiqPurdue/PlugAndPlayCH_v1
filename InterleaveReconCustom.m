function [image_i,sino,lam,ANG_D]=InterleaveReconCustom(SFData,FinSF,r,THETA,rr,lam0,m)

ReconsInSF=FinSF/r;

MAXTH=THETA.MAXTH;
MINTH=THETA.MINTH;

n_theta=THETA.n_theta;

fno=1;

for i=1:ReconsInSF
    
    for j=1:r
        
        % Club the projections together
        AllAngs(j,1:n_theta)=SFData.Frames(fno).Angles;
        AllProjs(j,1:n_theta,:)=SFData.Frames(fno).Projs.counts;
        %counts(j:r:r*n_theta,:)=SFData.Frames(f).Projs.counts;
        fno=fno+1;
        
    end
    
    Angs=reshape(AllAngs,1,r*n_theta);
    
    [~,~,np]=size(AllProjs);
    
    Projs=reshape(AllProjs,r*n_theta,np);
    
    image_i(i,:,:)=iradon(Projs',Angs,'cubic','Hamming',m);
    
    stepd=Angs(1,2)-Angs(1,1);

    [sino,ANG_D,ANG_R]=SetSinoGeom(MAXTH,MINTH,0,rr,stepd);

    sino.counts=Projs;
    
    sino(i,:)=sino;

    lam(i,:,:)=lam0*exp(-Projs);
    
end


end