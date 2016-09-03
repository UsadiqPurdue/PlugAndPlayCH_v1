function [cost,costp,costd]=ADMM_costs(error,weight,x,params)
% Cost assosiated with the data term optimization in Augmented Lagrangian
% error holds the current error as a vector
% weight is the weight assosiated with the error vector
% x is the current 2-D image
% params is a struct which holds the AL parameters: lambda, v and rho

costd = ((error.*weight)*error')/2;
temp = x-params.v+params.u;
costp = (params.lambda/2)*sum(sum((temp.*temp)));
cost = costd+costp;

end