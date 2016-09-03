
function [DenImg,Cost,CostH,CostF]= Denoise_AL_Poly_For_PnP_Java128_v3_HUpdate_Costs_v10(Jfunc,Filtering,img_in,img_orig,s,n,iterdn,sig,sig_H,foldname,a,b,inner1,ICDiters,HUpds)

% v6: The Denoising Java method does L_a<L_this ONLY if min<this<max

javaaddpath('/Users/USadiq/Work/Code-META/JavaCode/JProjectTest/build/classes');
% AL with

% Circorrection
% Fast
% Handles inner iterations to minimize L(x)
% Varies eps to adjust for precision
% Correct theta

randn('seed',1);
% Variables: CHANGE NAME AT EACH RUN


iter=iterdn;

lambda=1/(2*sig_H^2);


%%

filt2=fspecial('Gaussian',3,3);

if(Filtering)
    
    imf=imfilter(img_in,filt2);
    
else
    
    imf=img_in;
    
end

sig_H=1;

lambda=1/(2*sig_H^2);

mmin_n=min(min(min(img_orig)));
mmax_n=max(max(max(img_orig)));


mu_poly=zeros(s,s,n-1);

a_p=a;b_p=b;

img_x_poly=imf;

clear img_lp;

f1=Debug;

H00_poly=denoise_cost_Hf_impl_1_0528_01(img_x_poly,s,n,a_p,b_p);

tic

%inl=javaMethod('DenoiseLoopWithHUpdateCostsAHICD', f1,img_x_poly,img_in,n,s,a_p,b_p,sig,lambda, mmin_n,mmax_n,H00_poly,mu_poly,inner1,iter,100,3,5.0);

%inl=javaMethod(Jfunc, f1,img_x_poly,img_in,n,s,a_p,b_p,sig,lambda, mmin_n,mmax_n,H00_poly,mu_poly,inner1,iter);
HUpds=25;
HThresh=4e-5;



switch Jfunc
    
    case 'ObjectDenoiseTemp'
        
        ICDiters=100;
        HUpds=5;
        jDenObjThresh30=javaMethod('ObjectDenoiseTemp',f1, img_x_poly,img_x_poly,n,s,a,b,sig,lambda, mmin_n,mmax_n,mu_poly,ICDiters,HUpds);     
        
    case 'ObjectDenoiseCustNoH'
        
        ICDiters=250;
        HUpds=1;
        jDenObjThresh30=javaMethod('ObjectDenoiseTemp',f1, img_x_poly,img_x_poly,n,s,a,b,sig,lambda, mmin_n,mmax_n,mu_poly,ICDiters,HUpds);
    
    case 'ObjectDenoiseMore'
        
        ICDiters=500;
        HUpds=1;
        jDenObjThresh30=javaMethod('ObjectDenoiseTemp',f1, img_x_poly,img_x_poly,n,s,a,b,sig,lambda, mmin_n,mmax_n,mu_poly,ICDiters,HUpds);
   
    case 'ObjectDenoiseMore'
        
        ICDiters=400;
        HUpds=1;
        jDenObjThresh30=javaMethod('ObjectDenoiseTemp',f1, img_x_poly,img_x_poly,n,s,a,b,sig,lambda, mmin_n,mmax_n,mu_poly,ICDiters,HUpds);

    case 'ObjectDenoiseEverMore'
        
        ICDiters=200;
        HUpds=4;
        jDenObjThresh30=javaMethod('ObjectDenoiseTemp',f1, img_x_poly,img_x_poly,n,s,a,b,sig,lambda, mmin_n,mmax_n,mu_poly,ICDiters,HUpds);

    case 'CircObjectDenoiseMore'
        
        [img_x_poly]=FitInVessel(img_x_poly,s,n);
        
        ICDiters=400;
        HUpds=1;
        jDenObjThresh30=javaMethod('ObjectDenoiseTemp',f1, img_x_poly,img_x_poly,n,s,a,b,sig,lambda, mmin_n,mmax_n,mu_poly,ICDiters,HUpds);

    otherwise
        jDenObjThresh30=javaMethod(Jfunc, f1,img_x_poly,img_x_poly,n,s,a,b,sig,lambda, mmin_n,mmax_n,mu_poly,HThresh,HUpds);
        
        
end
% Bug.img=inl.img_x;Bug.H=inl.H00;Bug.mu=inl.mu;
% Bug.Index.i=inl.bug_i+1;Bug.Index.j=inl.bug_j+1;
% Bug.Index.k=inl.bug_k+1;
% Bug.y=inl.y;

DenImg=jDenObjThresh30.img_x;

CostH=jDenObjThresh30.jcostH;
CostF=jDenObjThresh30.jcostF;
Cost=jDenObjThresh30.jcost;

toc

end