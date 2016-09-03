function [Hx,AA,BB,CC]= denoise_cost_Hf_impl_1_0528_01(x,siz,n,a,b)
pix=siz*siz;

Hx=zeros(siz,siz,n-1);

% a=49;
% b=10;

for i=1:n-1
    dt = 10^(-1.475); %unit in second
    dt=1*10^(-3);
    dt=1;
    
    xp=x(:,:,i);
    xn=x(:,:,i+1);
    %dt=1;
    dis = 1.46e-3; % 1.46e-3 unit in millimeter    1.46 unit in micro
    
    dis=1e-3;%modified for testing
    dis=1;
    
    
    step=1;
    
    l = [0 1 0;1 -4 1;0 1 0]/(step*dis)^2;
    ll = [0  0  1  0  0;
        0  2 -8  2  0;
        1 -8 20 -8  1;
        0  2 -8  2  0;
        0  0  1  0  0]/(step*dis)^4;
    
    un=reshape(xn,pix,1);
    
    %th1=[eps c b a] with the functional: au3+bu2+cu
    %th1=1*10^(-4.5)*[0.6*10^(-6) 2 -6 4]';
    th1=[0.5 2 -6 4]'*1e-2;
    
    fa=xn;
    
    
    d2_1 = circfilter2(ll,fa,5);     
    
    drhs=-(2*fa)-(4*xp.^3)+(6*xp.^2);
    
    B=circfilter2(l,drhs,3);
    
%     d_1 = circfilter2(l,fa,3);
%     
%     d_03=circfilter2(l,xp.^3,3);
%     d_0=circfilter2(l,xp,3);
    
    rhs=(a*d2_1)+b*B;
    Hx(:,:,i)=(xn-xp)./dt+rhs;
    
    %B=-(2*d_1)-(d_03)+(3*d_0);
    
    AA(:,:,i)=d2_1;
    BB(:,:,i)=B;
    CC(:,:,i)=(xn-xp)./dt;
    
    % lu2_a = circfilter2(l,fa.^2,3); lu2_b = circfilter2(l,fa.^2,3);
    % lu3_a = circfilter2(l,fa.^3,3); lu3_b = circfilter2(l,fa.^3,3);
    %
    % lu_t=circfilter2(ll,fa,5);
    % llu_a = circfilter2(ll,fa,5);   llu_b = circfilter2(l,lu_b,3); % Changed here by Usman because this term needs to be filtered twice
    %
    % Beta(:,1) = reshape(-lu_t,pix,1);%-0.5*(llu_a+llu_b)
    % Beta(:,2) = reshape(lu_b,pix,1);%0.5*(lu_a+lu_b)
    % Beta(:,3) = reshape(lu2_b,pix,1);%0.5*(lu2_a+lu2_b)
    % Beta(:,4) = reshape(lu3_b,pix,1);%0.5*(lu3_a+lu3_b)
    % %reshape(lu3_b,pix,1);%0.5*(lu3_a+lu3_b)
    %
    % cons=1;
    %
    % Hx((i-1)*pix+1:i*pix,:)=(un-reshape(xp,pix,1))./dt-(Beta*th1);
    
end
end