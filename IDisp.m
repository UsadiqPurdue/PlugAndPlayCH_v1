function []=IDisp(x,varargin)

figure;imagesc(x,[0,1]);colormap(gray);colorbar;
format_image_for_publication(gcf);

if(length(varargin)>0)
    title(varargin{1});
end

end

