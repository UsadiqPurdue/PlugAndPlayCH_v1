function [image_i,sino,lam,ANG_D]=InterleaveRecon(SFData,FinSF,THETA,rr,lam0,m)

MAXTH=THETA.MAXTH;
MINTH=THETA.MINTH;

n_theta=THETA.n_theta;

stepd=(MAXTH-MINTH)/((n_theta*FinSF)-1);

[sino,ANG_D,ANG_R]=SetSinoGeom(MAXTH,MINTH,0,rr,stepd);

ssize=length(SFData.Frames(1).Projs.counts);

for f=1:FinSF
    
counts(f:FinSF:FinSF*n_theta,:)=SFData.Frames(f).Projs.counts;

end

image_i(:,:)=iradon(counts',ANG_D,'cubic','Hamming',m);

sino.counts=counts;

lam(:,:)=lam0*exp(-counts);

end