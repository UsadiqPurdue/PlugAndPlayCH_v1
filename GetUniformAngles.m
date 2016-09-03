function [ANG_D,ANG_R,step]=GetUniformAngles(max,min,offset,step)

if(max-min>=180)
    display('Some angles will be repeated');
end
    

%step=(max -  min)/((FinSF-1)/FinSF+(n_t-1));

ANG_D=min-offset:step:max-offset;

ANG_R=ANG_D.*pi/180;

    
end