%
% Cahn-Hilliard integrator.  Solves the CH equation
%   with natural and no-flux boundary conditions.
%   Nonlinear term: f(u) = u - u^3
%
% Numerical scheme;
%
%    Eyre's linearly stabilized CH integration scheme
%
% Notes;
%
%    1. Uses a fixed time step.
%
%    2. Uses discrete cosine transform to invert
%       the linear update matrix by performing
%       a spectral decomposition of the matrix.
%
%    3. Calculates spinodal decomposition from
%       random initial conditions with mean value
%       close to zero.
%
%    4. With the given parameters, there are roughly
%       10 grid points in a transition layer.
%
%    5. The timestep probably should be shorted
%       initially and lengthed later for a more
%       accurate solution.
%
%    6. The intermediate output is time and max(|U|)
%
%clearvars -except img_b img_orig_rate img_better img_stack

% spatial dimensions
function [img_n,img_s]=CH_gen_x_fft_recon_01_circ256_subsample_v9(s,n,sigma,a1,b1,limit,seed,ss_rate)

%prev version: _v3

format long;
N = s;  M = s;
delx = 1/(M-1);
delx=1;
%delx=1.27055;
delx2 = delx^2;
x = (0:delx:delx*M)';

img_s=zeros(N,M,n/ss_rate-1);
img_n=zeros(N,M,n/ss_rate-1);

randn('seed',seed);
rand('seed',seed);

% graphics parameters

visual_update = ss_rate;
type_update = 10;

% time parameters

offs1=0;

t = 0;
delt = 0.00005;
delt=1;
%delt=4189;
ntmax = limit+n+1;

img_no=1;

% CH parameters

% epsilon = 4;
% eps2 = epsilon^2;

% time-step parameter used in Eyre's scheme

a = 2;

% update parameters

lam1 = delt/delx2;
lam2 = a1*lam1/delx2;

% unscaled eigenvalues of the laplacian (nuemann bc)
h=[0,1,0;1,-4,1;0,1,0];
Leig=ifftshift(freqz2(h,M,N));


%Leig  = (((2*cos(pi*(0:N-1)'/(N-1)))-2)*ones(1,M)) + ...
%         (ones(N,1)*((2*cos(pi*(0:M-1)/(M-1)))-2));


% scaled eigenvalues of stabilized CH update matrix

CHeig = ones(N,M) - (b1*a*lam1*Leig) + (lam2*Leig.*Leig);

% scaled eigenvalues of the laplacian

Seig = lam1*Leig;

rad=floor(s/2);

% random initial conditions

offs=200;
U = rand(N,M);
%U=(img_stack(1:M,1:N,1)*2.0)-1;
hat_U = fft2(U);

% main loop

it = 0; j=0;
t = 0.0;
f1=figure;


fname=strcat('CH_256_fast_test');

AviName=strcat(fname,'.avi');
vidObj=VideoWriter(AviName);
vidObj.FrameRate=4;


 
open(vidObj);

fig=figure;

while it < ntmax
    
    it;
    if rem(it,type_update) == 0
        [it t max(max(abs(U)))];
    end
    
    if(it==limit)
         U=circ_section(U,s,s,rad-6,0.08);  
         %imwrite(U,'Phantom_g8.tiff');
        
        %U=imread('Phantom_g8.tiff');
        %U=double(U)./255.0;
        
    end
    
    if rem(it,visual_update) == 0 && it>=limit+offs1 && it<=limit+n+offs1
        

        x=[0:delx:delx*M];
        y=[0:delx:delx*N];
        
        imagesc(x,y,U);
        colormap(gray);
        colorbar;
        title(strcat('CH Solution at t=',num2str(it*delt),'s'));
        
        
        
        xlabel('\mu m');
        ylabel('\mu m');
        
        img_s(:,:,img_no)=U;
        img_n(:,:,img_no)=U+sigma.*randn(s,s);
        
        
        F=getframe(fig);
        writeVideo(vidObj,F);

        if it == 0
            mov = moviein(ntmax/visual_update,gcf);
        end
        j=j+1; mov(:,j) = getframe(gcf);
        
        img_no=img_no+1;
    end
    
    
    hat_U = fft2(U);
    % Update the solution
    
    it = it+1;
    t = t+delt;
    
    % compute the shifted nonlinear term
    
    fU = (4*U.*U.*U) - (6*U.^2);
    
    % compute the right hand side in tranform space
    
    hat_rhs = hat_U + b1*(Seig.*fft2(fU));
    
    % compute the updated solution in tranform space
    
    hat_U = hat_rhs./CHeig;
    
    % invert the cosine transform
    
    U = ifft2(hat_U);
    
    
end
close(vidObj);

img_s(:,:,size(img_s,3))=[];
end
