
function imgout=ImageInterp3D(img_x_poly,rate)

ss_rate=rate;
[s,~,n]=size(img_x_poly);

imgout=zeros(s,s,n*ss_rate);
for ii=1:s
    for jj=1:s
        
        imxy3(1:n)=img_x_poly(ii,jj,:);
        inc=(n*ss_rate-1)/(n-1);
        imgout(ii,jj,:)=interp1([1:inc:ss_rate*n],imxy3,[1:ss_rate*n]);
    end
end

end
