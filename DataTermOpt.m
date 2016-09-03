function [new_image,cost_d,cost_p]=DataTermOpt(map_image,sino,d,geom,params,Amatrix)
% Optimizes the function : data mismatch term + augmented vector
% min_x ((1/2)||y-Ax||_D^{2} + lambda ||x-(v-rho)||^2 )
% Inputs: map_image : Initial value of x 
%        sino : sinogram structure containing data
%        geom : geometry structure contianing geom info of object 
%        params : prior model parameters + algorithm parameters structure
%        Amatrix : A stored as a userdefined sparse vector 

[m n]=size(map_image);
%d=reshape(sino.counts',1,sino.n_t*sino.n_theta);%Diagonal covariance matrix entries

d=reshape(d,1,sino.n_t*sino.n_theta);
Ax=forward_project_v2(map_image,sino,Amatrix); %returns a n_theta X n_t matrix

if(isfield(params,'xray') && params.xray == 1)
    sino.counts=-log(sino.counts./params.dose);
end

e=(sino.counts-Ax)'; %n_t X n_theta
e=reshape(e,1,sino.n_t*sino.n_theta);

params.num_iter;
if(params.verbose)
cost = zeros(1,params.num_iter);
end

% if(~isfield(params,'xray') || params.xray == 0) %if this is not transmission data use the corret noise model
%     for i=1:length(d)
%         if(d(i)~=0)
%             d(i)=1/d(i);
%         else
%             d(i)=0;
%         end
%     end
% end

%TempK = params.v-params.u; %A temporary variable used to store v-u

tic;

perCostChange=0;cost_pd=1e5;

for k=1:params.num_iter
    
    %computing cost after each iteration
    if(params.verbose)
        cost(k) = Cost_DataTermAL(e,d,map_image,params);
    end
    prev_img = map_image;
    
    [map_image,e,cost_d(k,:),cost_p(k,:),cost_T(k,:)]=HomogenousUpdate(map_image,Amatrix,e,d,params);
    
    
    %perCostChange=(cost_p-cost_pd)/cost_pd;
    
%     if(mod(k,2) ~= 0) %alternate between homogenous and non-homogenous updates
%         [map_image,e]=HomogenousUpdate(map_image,Amatrix,e,d,params);
%     else
%         K=20;
%         [map_image,e]=NonHomogenousUpdate(map_image,Amatrix,e,d,params,K,UpdateMap);
%     end
    
    UpdateMap=abs(map_image-prev_img);
    
    cost_pd=cost_d;

       
end

timeopt=toc;
fprintf('Time taken by Data Term Opt. %d\n',timeopt);


new_image = map_image;





