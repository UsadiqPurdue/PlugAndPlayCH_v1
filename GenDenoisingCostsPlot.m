function []=GenDenoisingCostsPlot(ResultsFolder)

close all;

load(strcat(ResultsFolder,'PnP_InitData.mat'));
XLabel='PnP Iterations';



load(strcat(ResultsFolder,'v_hat.mat'));
[iter,s,~,n] =size(v_hat);

pix=s*s;


vHP=strcat(ResultsFolder,'vH.tif');
xHP=strcat(ResultsFolder,'xH.tif');
xtP=strcat(ResultsFolder,'xT.tif');
vtP=strcat(ResultsFolder,'vT.tif');
muP=strcat(ResultsFolder,'mu.tif');

[costHvh]=HPlot(v_hat,GetTitle('\hat{v}',sig,sig_H),XLabel,s,n,a,b,iter);
clear 'v_hat';
saveas(gcf,vHP);


load(strcat(ResultsFolder,'x_hat.mat'));
[costHxh]=HPlot(x_hat,GetTitle('\hat{x}',sig,sig_H),XLabel,s,n,a,b,iter);
saveas(gcf,xHP);
clear 'x_hat';


load(strcat(ResultsFolder,'x_til.mat'));
[costHxtil]=HPlot(x_til,GetTitle('\tilde{x}',sig,sig_H),XLabel,s,n,a,b,iter);
saveas(gcf,xtP);
clear 'x_til';

load(strcat(ResultsFolder,'v_til.mat'));
[costHvtil]=HPlot(v_tilda,GetTitle('\tilde{v}',sig,sig_H),XLabel,s,n,a,b,iter);
saveas(gcf,vtP);
clear 'v_tilda';

load(strcat(ResultsFolder,'v_hat.mat'));
load(strcat(ResultsFolder,'x_hat.mat'));

for i=1:iter
    xvdiff(i)=sum(sum(sum((x_hat(i,:,:,:)-v_hat(i,:,:,:)).^2)));
end

clear 'v_hat';
clear 'x_hat';

plot([1:iter],(xvdiff/(pix*n)),'-b','LineWidth',1.5);
format_image_for_publication(gcf);
title(strcat('$$(\hat{x} - \hat{v})^2 ,\gamma=',num2str(sig),', \sigma_H=',num2str(sig_H),'$$'),'Interpreter','latex');
xlabel('iterations','Interpreter','latex');
saveas(gcf,muP);

% for i=1:iter
%     img_x_poly(:,:,:)=x_hat(i,:,:,:);
%     H00_poly=denoise_cost_Hf_impl_1_0528_01(img_x_poly,s,n,a,b);
%     Hc=sum(sum(sum(H00_poly.^2)));
%     costHx(i)=(1/100)*sqrt(Hc/(pix*n));
% end
% 
% inner1=iter;
% plot([1:inner1],(1/100)*sqrt(costHx(1:inner1)/(pix*n)),'-b','LineWidth',1.5);
% format_image_for_publication(gcf);
% title('$$RMSH(\hat{x})$$','Interpreter','Latex');
% xlabel('PnP iterations','Interpreter','latex');

load(strcat(ResultsFolder,'1Costs.mat'));

%load('/Users/USadiq/Results/CHPhantomDataResults/NoHUpdates/InterlacedViews_Ang45r1K8Nt33MaxTH179_30-Sep-2015/1Costs.mat');
inner1=100;
pix=s*s;
figure;
plot([1:inner1],(1/100)*sqrt(cost_H(1:inner1)/(pix*n)),'-b','LineWidth',1.5);
format_image_for_publication(gcf);
title('RMSH/s','Interpreter','latex');
xlabel('iterations','Interpreter','latex');







end
