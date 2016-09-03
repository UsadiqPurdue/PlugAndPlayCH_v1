function [img_out]=circ_section(img_in,m,n,r,value)

    img_out=ones(m,n)*value;
    
    cnt_x=floor(n/2);
    cnt_y=floor(m/2);
    
    for i=1:m
        
        for j=1:n
            
            i1=abs(i-cnt_x);
            j1=abs(j-cnt_y);
            
            if(i1^2+j1^2<=r^2)
                img_out(i,j)=img_in(i,j);
            end
        end
    end

    
end