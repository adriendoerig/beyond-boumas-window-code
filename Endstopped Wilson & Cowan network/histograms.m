% Creates histograms for the data. May need some fiddling
% depending on the data set.

P = SigmoidML(FitParams,CrossCorr);

figure;
subplot(3,1,1);
BarH = bar(1:size(Data,2),Data(1,:)); set(BarH,'FaceColor','k');
title('Octagon psych fct');
hold on; errorbar(1:size(Data,2),Data(1,:),Data(2,:),'k.'); hold off;

subplot(3,1,2);
BarH = bar(1:size(Data,2)-1,P); set(BarH,'FaceColor','k');

for ind = 1:size(Img,3)-1
    subplot(3,size(Img,3)-1,ind+2*(size(Img,3)-1));
    imagesc(Img(:,:,ind+1));
    axis image;
    axis off;
end
colormap gray;


Residuals = sqrt((P-Data(1,2:length(Data(1,:)))).^2);