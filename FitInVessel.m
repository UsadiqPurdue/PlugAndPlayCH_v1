function [img_out]=FitInVessel(img_in,s,n)

rad=floor(s/2);

for i=1:n
    
    U=img_in(:,:,i);
    
    img_out(:,:,i)=circ_section(U,s,s,rad,0.08); 
end

end