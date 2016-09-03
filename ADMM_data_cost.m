function [cost,costp,costd]=ADMM_data_cost(map_image,sino,d,params,Amatrix)


[m n]=size(map_image);
%d=reshape(sino.counts',1,sino.n_t*sino.n_theta);%Diagonal covariance matrix entries

d=reshape(d,1,sino.n_t*sino.n_theta);
Ax=forward_project_v2(map_image,sino,Amatrix); %returns a n_theta X n_t matrix

if(isfield(params,'xray') && params.xray == 1)
    sino.counts=-log(sino.counts./params.dose);
end

e=(sino.counts-Ax)'; %n_t X n_theta
e=reshape(e,1,sino.n_t*sino.n_theta);

[cost,costp,costd]=ADMM_costs(e,d,map_image,params);

end